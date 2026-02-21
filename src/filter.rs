/// Pipeline: extract child k-mers, scan parents, annotate proband-unique read counts.
///
/// For each variant:
/// 1. Extract k-mers from child reads that overlap the variant position.
/// 2. Scan parent BAM/CRAM files for those k-mers (whole-file, aligner-agnostic).
/// 3. Count child reads with at least one k-mer not found in either parent.
/// 4. Annotate as KMER_PROBAND_UNIQUE.
use ahash::AHashSet;
use log::{debug, info, warn};
use std::time::Instant;

use crate::counter;
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
}

/// Per-variant child read k-mers.
struct VariantReadKmers {
    record_idx: usize,
    per_read_kmers: Vec<AHashSet<Vec<u8>>>,
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
    let all_child_kmers: AHashSet<Vec<u8>> = variant_read_kmers
        .iter()
        .flat_map(|vrk| vrk.per_read_kmers.iter().flat_map(|s| s.iter().cloned()))
        .collect();
    info!(
        "Collected {} distinct child k-mers for parent scanning",
        all_child_kmers.len()
    );

    // 4. Scan each parent's entire BAM/CRAM
    info!("Scanning mother BAM/CRAM (whole-file)...");
    let mother_kmers = counter::count_kmers_whole_file(
        mother_bam,
        ref_path,
        &all_child_kmers,
        config.kmer_size,
        config.min_baseq,
    )?;
    info!(
        "Mother scan: found {} distinct child k-mers",
        mother_kmers.len()
    );

    info!("Scanning father BAM/CRAM (whole-file)...");
    let father_kmers = counter::count_kmers_whole_file(
        father_bam,
        ref_path,
        &all_child_kmers,
        config.kmer_size,
        config.min_baseq,
    )?;
    info!(
        "Father scan: found {} distinct child k-mers",
        father_kmers.len()
    );

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
                    .any(|km| !mother_kmers.contains_key(km) && !father_kmers.contains_key(km))
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
