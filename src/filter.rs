/// Pipeline: extract child k-mers, scan parents, annotate proband-unique read counts.
///
/// For each variant:
/// 1. Extract k-mers from child reads that overlap the variant position.
/// 2. Scan parent BAM/CRAM files for those k-mers (whole-file, aligner-agnostic),
///    or load pre-built k-mer databases from disk.
/// 3. Count child reads with at least one k-mer not found in either parent.
/// 4. Annotate as KMER_PROBAND_UNIQUE.
use ahash::AHashSet;
use log::{debug, info, warn};
use std::time::Instant;

use crate::counter;
use crate::kmer_db::KmerDb;
use crate::metrics::FilterSummary;
use crate::vcf_io::{self, VariantAnnotation};

/// Configuration for the pipeline.
#[derive(Debug, Clone)]
pub struct FilterConfig {
    pub kmer_size: usize,
    pub min_baseq: u8,
    pub min_mapq: u8,
    pub max_reads_per_locus: usize,
    pub debug_kmers: bool,
    pub threads: usize,
}

/// Per-variant child read k-mers (u128-encoded).
struct VariantReadKmers {
    record_idx: usize,
    per_read_kmers: Vec<AHashSet<u128>>,
}

/// Main pipeline entry point.
pub fn run_filter(
    child_bam: &str,
    mother_bam: &str,
    father_bam: &str,
    ref_path: &str,
    vcf_path: &str,
    output_path: &str,
    config: &FilterConfig,
    mother_db_path: Option<&str>,
    father_db_path: Option<&str>,
) -> Result<FilterSummary, Box<dyn std::error::Error>> {
    let start_time = Instant::now();

    // 1. Load variants
    info!("Reading variants from {}", vcf_path);
    let variants = vcf_io::read_variants(vcf_path)?;
    info!("Loaded {} candidate variants", variants.len());

    if variants.is_empty() {
        vcf_io::write_annotated_vcf(vcf_path, output_path, &[])?;
        return Ok(FilterSummary {
            total_variants: 0,
            kmer_size: config.kmer_size,
            runtime_seconds: start_time.elapsed().as_secs_f64(),
        });
    }

    // 2. Extract per-read k-mers from child at each variant
    info!("Extracting child read k-mers per variant...");
    let mut variant_read_kmers: Vec<VariantReadKmers> = Vec::with_capacity(variants.len());

    for var in &variants {
        let per_read = counter::extract_reads_kmers_at_variant(
            child_bam,
            ref_path,
            &var.chrom,
            var.pos,
            var.ref_allele.len(),
            config.kmer_size,
            config.min_mapq,
            config.min_baseq,
            config.max_reads_per_locus,
            config.threads,
        )
        .unwrap_or_else(|e| {
            warn!(
                "Child k-mer extraction error for {}:{}: {}",
                var.chrom,
                var.pos + 1,
                e
            );
            Vec::new()
        });

        variant_read_kmers.push(VariantReadKmers {
            record_idx: var.record_idx,
            per_read_kmers: per_read,
        });
    }

    // 3. Collect all child k-mers for parent scanning
    let all_child_kmers: AHashSet<u128> = variant_read_kmers
        .iter()
        .flat_map(|vrk| vrk.per_read_kmers.iter().flat_map(|s| s.iter().copied()))
        .collect();
    info!(
        "Collected {} distinct child k-mers for parent scanning",
        all_child_kmers.len()
    );

    // 4. Load or scan parent k-mer databases
    let mother_db = load_or_scan_parent(
        "mother",
        mother_bam,
        ref_path,
        mother_db_path,
        &all_child_kmers,
        config,
    )?;
    let father_db = load_or_scan_parent(
        "father",
        father_bam,
        ref_path,
        father_db_path,
        &all_child_kmers,
        config,
    )?;

    // 5. Count proband-unique reads per variant
    info!("Computing proband-unique read counts...");
    let mut annotations: Vec<VariantAnnotation> = Vec::with_capacity(variant_read_kmers.len());

    for vrk in &variant_read_kmers {
        let proband_unique = vrk
            .per_read_kmers
            .iter()
            .filter(|read_kmers| {
                read_kmers
                    .iter()
                    .any(|km| !mother_db.contains(*km) && !father_db.contains(*km))
            })
            .count() as i32;

        if config.debug_kmers {
            debug!(
                "Variant record_idx={} - reads={}, proband_unique_reads={}",
                vrk.record_idx,
                vrk.per_read_kmers.len(),
                proband_unique,
            );
        }

        annotations.push(VariantAnnotation {
            record_idx: vrk.record_idx,
            kmer_k: config.kmer_size as i32,
            proband_unique,
        });
    }

    annotations.sort_by_key(|a| a.record_idx);

    let summary = FilterSummary {
        total_variants: annotations.len(),
        kmer_size: config.kmer_size,
        runtime_seconds: start_time.elapsed().as_secs_f64(),
    };

    // 6. Write output VCF
    info!("Writing annotated VCF to {}", output_path);
    vcf_io::write_annotated_vcf(vcf_path, output_path, &annotations)?;

    Ok(summary)
}

/// Load a pre-built parent k-mer database from disk, or scan the BAM file.
///
/// If `db_path` is provided:
///   - If the file exists, load the database from disk (avoiding a BAM scan).
///   - If the file does not exist, scan the BAM, build the database, and save it.
/// If `db_path` is None, perform a targeted scan for only the child k-mers.
fn load_or_scan_parent(
    label: &str,
    bam_path: &str,
    ref_path: &str,
    db_path: Option<&str>,
    target_kmers: &AHashSet<u128>,
    config: &FilterConfig,
) -> Result<KmerDb, Box<dyn std::error::Error>> {
    match db_path {
        Some(path) if std::path::Path::new(path).exists() => {
            info!("Loading {} k-mer database from {}", label, path);
            let db = KmerDb::load(path)?;
            info!(
                "{} database: loaded {} k-mers (k={})",
                label,
                db.len(),
                db.kmer_size()
            );
            if db.kmer_size() != config.kmer_size {
                return Err(format!(
                    "{} database k-mer size ({}) does not match config ({})",
                    label,
                    db.kmer_size(),
                    config.kmer_size
                )
                .into());
            }
            Ok(db)
        }
        Some(path) => {
            info!(
                "Building {} k-mer database (whole-file scan) and saving to {}",
                label, path
            );
            let counts = counter::count_all_kmers_whole_file(
                bam_path,
                ref_path,
                config.kmer_size,
                config.min_baseq,
                config.threads,
            )?;
            info!(
                "{} scan complete: {} distinct k-mers",
                label,
                counts.len()
            );
            let db = KmerDb::from_counts(config.kmer_size, counts);
            db.save(path)?;
            info!("{} database saved to {}", label, path);
            Ok(db)
        }
        None => {
            info!("Scanning {} BAM/CRAM (whole-file, targeted)...", label);
            let counts = counter::count_kmers_whole_file(
                bam_path,
                ref_path,
                target_kmers,
                config.kmer_size,
                config.min_baseq,
                config.threads,
            )?;
            info!(
                "{} scan: found {} distinct child k-mers",
                label,
                counts.len()
            );
            Ok(KmerDb::from_counts(config.kmer_size, counts))
        }
    }
}
