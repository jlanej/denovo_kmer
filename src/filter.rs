/// Scoring, filtering, and main pipeline logic.
///
/// The pipeline uses a two-phase approach for parent vs child counting:
/// - **Child**: Region-based counting around the variant locus (trusts alignment).
/// - **Parents**: Whole-file aligner-agnostic scanning of ALL reads, with all
///   ALT k-mers batched together for a single pass per parent file. This catches
///   inherited variants even when parent reads are mismapped to different loci.
use ahash::AHashSet;
use log::{debug, info, warn};
use rayon::prelude::*;
use std::sync::Mutex;
use std::time::Instant;

use crate::counter;
use crate::haplotype::{self, Variant};
use crate::metrics::FilterSummary;
use crate::vcf_io::{self, VariantAnnotation};

/// Configuration for the filtering pipeline.
#[derive(Debug, Clone)]
pub struct FilterConfig {
    pub kmer_size: usize,
    pub min_baseq: u8,
    pub min_mapq: u8,
    pub max_reads_per_locus: usize,
    pub min_child_alt: u32,
    pub max_parent_alt: u32,
    pub min_child_alt_ratio: f64,
    pub max_variant_size: usize,
    pub window: usize,
    pub debug_kmers: bool,
}

/// Precomputed k-mer information for a single variant.
struct VariantKmerInfo {
    /// VCF record index for output tracking
    record_idx: usize,
    /// Chromosome
    chrom: String,
    /// 0-based position
    pos: i64,
    /// Length of reference allele
    ref_allele_len: usize,
    /// Canonical ALT-unique k-mers for this variant
    alt_kmer_set: AHashSet<Vec<u8>>,
    /// Canonical REF-unique k-mers for this variant
    ref_kmer_set: AHashSet<Vec<u8>>,
    /// All distinct canonical k-mers from child reads at this variant (reference-free)
    child_read_kmers: AHashSet<Vec<u8>>,
    /// Whether all ALT k-mers are low complexity
    low_complexity: bool,
    /// Whether variant type/size is unsupported
    not_implemented: bool,
}

/// Result of scoring a single variant.
#[derive(Debug)]
struct ScoreResult {
    child_alt_count: u32,
    child_ref_count: u32,
    mother_alt_count: u32,
    father_alt_count: u32,
    proband_unique: u32,
    low_complexity: bool,
    not_implemented: bool,
}

/// Build haplotype k-mers for all variants, loading reference per contig.
fn build_all_variant_kmers(
    variants: &[Variant],
    ref_path: &str,
    config: &FilterConfig,
) -> Result<Vec<VariantKmerInfo>, Box<dyn std::error::Error>> {
    let mut by_contig: ahash::AHashMap<String, Vec<&Variant>> = ahash::AHashMap::new();
    for var in variants {
        by_contig.entry(var.chrom.clone()).or_default().push(var);
    }

    let mut infos: Vec<VariantKmerInfo> = Vec::with_capacity(variants.len());

    for (contig, contig_vars) in &by_contig {
        let ref_seq = load_reference_contig(ref_path, contig)
            .map_err(|e| -> Box<dyn std::error::Error> { e })?;

        for var in contig_vars {
            if var.size() > config.max_variant_size {
                infos.push(VariantKmerInfo {
                    record_idx: var.record_idx,
                    chrom: var.chrom.clone(),
                    pos: var.pos,
                    ref_allele_len: var.ref_allele.len(),
                    alt_kmer_set: AHashSet::new(),
                    ref_kmer_set: AHashSet::new(),
                    child_read_kmers: AHashSet::new(),
                    low_complexity: false,
                    not_implemented: true,
                });
                continue;
            }

            let hap_kmers = haplotype::build_haplotype_kmers(&ref_seq, var, config.kmer_size);
            let alt_set: AHashSet<Vec<u8>> = hap_kmers.alt_kmers.into_iter().collect();
            let ref_set: AHashSet<Vec<u8>> = hap_kmers.ref_kmers.into_iter().collect();
            let low_complexity = hap_kmers.low_complexity || alt_set.is_empty();

            infos.push(VariantKmerInfo {
                record_idx: var.record_idx,
                chrom: var.chrom.clone(),
                pos: var.pos,
                ref_allele_len: var.ref_allele.len(),
                alt_kmer_set: alt_set,
                ref_kmer_set: ref_set,
                child_read_kmers: AHashSet::new(),
                low_complexity,
                not_implemented: false,
            });
        }
    }

    Ok(infos)
}

/// Determine FILTER status from scores.
fn determine_filter(score: &ScoreResult, config: &FilterConfig) -> String {
    if score.not_implemented {
        return vcf_io::FILTER_NOT_IMPLEMENTED.to_string();
    }
    if score.low_complexity {
        return vcf_io::FILTER_LOW_COMPLEXITY.to_string();
    }

    // Check parent contamination first
    if score.mother_alt_count > config.max_parent_alt
        || score.father_alt_count > config.max_parent_alt
    {
        return vcf_io::FILTER_PARENT_ALT.to_string();
    }

    // Check child ALT support
    if score.child_alt_count < config.min_child_alt {
        return vcf_io::FILTER_NO_CHILD_ALT.to_string();
    }

    // Check ALT ratio
    let alt_ratio =
        score.child_alt_count as f64 / (score.child_ref_count as f64 + 1.0);
    if alt_ratio < config.min_child_alt_ratio {
        return vcf_io::FILTER_NO_CHILD_ALT.to_string();
    }

    vcf_io::FILTER_OK.to_string()
}

/// Load reference sequence for a contig from an indexed FASTA.
fn load_reference_contig(
    ref_path: &str,
    contig: &str,
) -> Result<Vec<u8>, Box<dyn std::error::Error + Send + Sync>> {
    let reader = rust_htslib::faidx::Reader::from_path(ref_path)
        .map_err(|e| format!("Cannot open reference {}: {}", ref_path, e))?;
    // Fetch the whole contig by requesting up to a large end position.
    // htslib clamps the request to the actual contig length.
    let seq = reader
        .fetch_seq_string(contig, 0, 1 << 30)
        .map_err(|e| format!("Cannot fetch contig {} from {}: {}", contig, ref_path, e))?;
    Ok(seq.into_bytes())
}

/// Main pipeline entry point.
///
/// Pipeline phases:
/// 1. Load variants and build haplotype k-mers for all variants.
/// 2. Collect all ALT k-mers into a unified set and scan each parent's
///    entire BAM/CRAM in ONE pass (aligner-agnostic, no mapQ filter).
/// 3. For each variant, count child k-mers via region-based queries and
///    look up parent counts from the whole-file scan results.
/// 4. Score, annotate, and write output VCF.
pub fn run_filter(
    child_bam: &str,
    mother_bam: &str,
    father_bam: &str,
    ref_path: &str,
    vcf_path: &str,
    output_path: &str,
    config: &FilterConfig,
    threads: usize,
    _regions: Option<&str>,
) -> Result<FilterSummary, Box<dyn std::error::Error>> {
    let start_time = Instant::now();

    // Configure thread pool
    if threads > 1 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .unwrap_or_else(|_| {
                warn!("Could not set thread pool size, using default");
            });
    }

    // Phase 1: Load variants and build haplotype k-mers
    info!("Reading variants from {}", vcf_path);
    let variants = vcf_io::read_variants(vcf_path)?;
    info!("Loaded {} candidate variants", variants.len());

    if variants.is_empty() {
        vcf_io::write_annotated_vcf(vcf_path, output_path, &[])?;
        return Ok(FilterSummary {
            total_variants: 0,
            kmer_size: config.kmer_size,
            runtime_seconds: start_time.elapsed().as_secs_f64(),
            ..Default::default()
        });
    }

    info!("Building haplotype k-mers for {} variants...", variants.len());
    let mut var_infos = build_all_variant_kmers(&variants, ref_path, config)?;

    // Phase 1.5: Extract all k-mers from child reads per variant (reference-free).
    // This enables RUFUS-like proband-unique k-mer annotation: k-mers present
    // in the child's reads but absent from both parents.
    info!("Extracting child read k-mers per variant (reference-free)...");
    for vi in var_infos.iter_mut() {
        if vi.not_implemented || vi.low_complexity {
            continue;
        }
        let start = (vi.pos - config.window as i64).max(0);
        let end = vi.pos + vi.ref_allele_len as i64 + config.window as i64;

        match counter::extract_all_kmers_in_region(
            child_bam,
            ref_path,
            &vi.chrom,
            start,
            end,
            config.kmer_size,
            config.min_mapq,
            config.min_baseq,
            config.max_reads_per_locus,
        ) {
            Ok(kmers) => {
                vi.child_read_kmers = kmers;
            }
            Err(e) => {
                warn!(
                    "Child k-mer extraction error for {}:{}: {}",
                    vi.chrom,
                    vi.pos + 1,
                    e
                );
            }
        }
    }

    // Phase 2: Whole-file parent scanning (aligner-agnostic)
    // Collect ALL ALT k-mers plus ALL child read k-mers into one set for
    // efficient batched scanning. This supports both existing ALT-based filtering
    // and the new proband-unique k-mer annotation.
    let all_target_kmers: AHashSet<Vec<u8>> = var_infos
        .iter()
        .filter(|vi| !vi.low_complexity && !vi.not_implemented)
        .flat_map(|vi| {
            vi.alt_kmer_set
                .iter()
                .chain(vi.child_read_kmers.iter())
                .cloned()
        })
        .collect();
    info!(
        "Collected {} unique target k-mers for parent scanning (ALT + child read k-mers)",
        all_target_kmers.len()
    );

    // Scan each parent's entire BAM/CRAM in ONE pass. No mapQ filter is applied
    // so that reads mismapped to other loci are still checked for ALT k-mers.
    info!("Scanning mother BAM/CRAM (whole-file, aligner-agnostic)...");
    let mother_kmer_counts = counter::count_kmers_whole_file(
        mother_bam,
        ref_path,
        &all_target_kmers,
        config.kmer_size,
        config.min_baseq,
    )?;
    info!(
        "Mother scan complete: found {} distinct target k-mers in reads",
        mother_kmer_counts.len()
    );

    info!("Scanning father BAM/CRAM (whole-file, aligner-agnostic)...");
    let father_kmer_counts = counter::count_kmers_whole_file(
        father_bam,
        ref_path,
        &all_target_kmers,
        config.kmer_size,
        config.min_baseq,
    )?;
    info!(
        "Father scan complete: found {} distinct target k-mers in reads",
        father_kmer_counts.len()
    );

    // Phase 3: Per-variant child counting (region-based) and scoring.
    // Child reads are queried by region since we trust the child alignment
    // for the candidate variant loci.
    info!("Counting child k-mers (region-based) and scoring variants...");
    let annotations: Mutex<Vec<VariantAnnotation>> = Mutex::new(Vec::new());

    var_infos.par_iter().for_each(|vi| {
        let score = if vi.not_implemented {
            ScoreResult {
                child_alt_count: 0,
                child_ref_count: 0,
                mother_alt_count: 0,
                father_alt_count: 0,
                proband_unique: 0,
                low_complexity: false,
                not_implemented: true,
            }
        } else if vi.low_complexity {
            ScoreResult {
                child_alt_count: 0,
                child_ref_count: 0,
                mother_alt_count: 0,
                father_alt_count: 0,
                proband_unique: 0,
                low_complexity: true,
                not_implemented: false,
            }
        } else {
            // Child: region-based counting around the variant locus
            let start = (vi.pos - config.window as i64).max(0);
            let end = vi.pos + vi.ref_allele_len as i64 + config.window as i64;

            let child_alt_count = counter::count_target_kmers_in_region(
                child_bam,
                ref_path,
                &vi.chrom,
                start,
                end,
                &vi.alt_kmer_set,
                config.kmer_size,
                config.min_mapq,
                config.min_baseq,
                config.max_reads_per_locus,
            )
            .unwrap_or_else(|e| {
                warn!(
                    "Child ALT counting error for {}:{}: {}",
                    vi.chrom,
                    vi.pos + 1,
                    e
                );
                0
            });

            let child_ref_count = counter::count_target_kmers_in_region(
                child_bam,
                ref_path,
                &vi.chrom,
                start,
                end,
                &vi.ref_kmer_set,
                config.kmer_size,
                config.min_mapq,
                config.min_baseq,
                config.max_reads_per_locus,
            )
            .unwrap_or_else(|e| {
                warn!(
                    "Child REF counting error for {}:{}: {}",
                    vi.chrom,
                    vi.pos + 1,
                    e
                );
                0
            });

            // Parents: look up counts from the whole-file scan results
            let mother_alt_count: u32 = vi
                .alt_kmer_set
                .iter()
                .map(|km| mother_kmer_counts.get(km).copied().unwrap_or(0))
                .sum();
            let father_alt_count: u32 = vi
                .alt_kmer_set
                .iter()
                .map(|km| father_kmer_counts.get(km).copied().unwrap_or(0))
                .sum();

            // Proband-unique k-mers: child read k-mers not found in either parent
            let proband_unique: u32 = vi
                .child_read_kmers
                .iter()
                .filter(|km| {
                    !mother_kmer_counts.contains_key(*km)
                        && !father_kmer_counts.contains_key(*km)
                })
                .count() as u32;

            if config.debug_kmers {
                debug!(
                    "Variant {}:{} - child_alt={}, child_ref={}, mother_alt={}, father_alt={}, alt_kmers={}, ref_kmers={}, child_read_kmers={}, proband_unique={}",
                    vi.chrom,
                    vi.pos + 1,
                    child_alt_count,
                    child_ref_count,
                    mother_alt_count,
                    father_alt_count,
                    vi.alt_kmer_set.len(),
                    vi.ref_kmer_set.len(),
                    vi.child_read_kmers.len(),
                    proband_unique,
                );
            }

            ScoreResult {
                child_alt_count,
                child_ref_count,
                mother_alt_count,
                father_alt_count,
                proband_unique,
                low_complexity: false,
                not_implemented: false,
            }
        };

        let filter = determine_filter(&score, config);
        let ann = VariantAnnotation {
            record_idx: vi.record_idx,
            filter: filter.clone(),
            kmer_k: config.kmer_size as i32,
            child_alt: score.child_alt_count as i32,
            child_ref: score.child_ref_count as i32,
            mother_alt: score.mother_alt_count as i32,
            father_alt: score.father_alt_count as i32,
            proband_unique: score.proband_unique as i32,
            status: filter,
        };

        annotations.lock().unwrap().push(ann);
    });

    let mut anns = annotations.into_inner().unwrap();
    // Sort by record index to maintain VCF order
    anns.sort_by_key(|a| a.record_idx);

    // Compute summary
    let mut summary = FilterSummary {
        total_variants: anns.len(),
        kmer_size: config.kmer_size,
        ..Default::default()
    };
    for ann in &anns {
        match ann.filter.as_str() {
            vcf_io::FILTER_OK => summary.dn_kmer_ok += 1,
            vcf_io::FILTER_NO_CHILD_ALT => summary.dn_kmer_no_child_alt += 1,
            vcf_io::FILTER_PARENT_ALT => summary.dn_kmer_parent_alt += 1,
            vcf_io::FILTER_LOW_COMPLEXITY => summary.dn_kmer_low_complexity += 1,
            vcf_io::FILTER_NOT_IMPLEMENTED => summary.dn_kmer_not_implemented += 1,
            _ => {}
        }
    }
    summary.runtime_seconds = start_time.elapsed().as_secs_f64();

    // Phase 4: Write output VCF
    info!("Writing annotated VCF to {}", output_path);
    vcf_io::write_annotated_vcf(vcf_path, output_path, &anns)?;

    Ok(summary)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_determine_filter_ok() {
        let config = FilterConfig {
            kmer_size: 31,
            min_baseq: 20,
            min_mapq: 20,
            max_reads_per_locus: 200,
            min_child_alt: 3,
            max_parent_alt: 1,
            min_child_alt_ratio: 0.1,
            max_variant_size: 50,
            window: 500,
            debug_kmers: false,
        };

        let score = ScoreResult {
            child_alt_count: 10,
            child_ref_count: 15,
            mother_alt_count: 0,
            father_alt_count: 0,
            proband_unique: 5,
            low_complexity: false,
            not_implemented: false,
        };
        assert_eq!(determine_filter(&score, &config), vcf_io::FILTER_OK);
    }

    #[test]
    fn test_determine_filter_no_child_alt() {
        let config = FilterConfig {
            kmer_size: 31,
            min_baseq: 20,
            min_mapq: 20,
            max_reads_per_locus: 200,
            min_child_alt: 3,
            max_parent_alt: 1,
            min_child_alt_ratio: 0.1,
            max_variant_size: 50,
            window: 500,
            debug_kmers: false,
        };

        let score = ScoreResult {
            child_alt_count: 1,
            child_ref_count: 20,
            mother_alt_count: 0,
            father_alt_count: 0,
            proband_unique: 0,
            low_complexity: false,
            not_implemented: false,
        };
        assert_eq!(
            determine_filter(&score, &config),
            vcf_io::FILTER_NO_CHILD_ALT
        );
    }

    #[test]
    fn test_determine_filter_parent_alt() {
        let config = FilterConfig {
            kmer_size: 31,
            min_baseq: 20,
            min_mapq: 20,
            max_reads_per_locus: 200,
            min_child_alt: 3,
            max_parent_alt: 1,
            min_child_alt_ratio: 0.1,
            max_variant_size: 50,
            window: 500,
            debug_kmers: false,
        };

        let score = ScoreResult {
            child_alt_count: 10,
            child_ref_count: 15,
            mother_alt_count: 5,
            father_alt_count: 0,
            proband_unique: 0,
            low_complexity: false,
            not_implemented: false,
        };
        assert_eq!(
            determine_filter(&score, &config),
            vcf_io::FILTER_PARENT_ALT
        );
    }

    #[test]
    fn test_determine_filter_low_complexity() {
        let config = FilterConfig {
            kmer_size: 31,
            min_baseq: 20,
            min_mapq: 20,
            max_reads_per_locus: 200,
            min_child_alt: 3,
            max_parent_alt: 1,
            min_child_alt_ratio: 0.1,
            max_variant_size: 50,
            window: 500,
            debug_kmers: false,
        };

        let score = ScoreResult {
            child_alt_count: 0,
            child_ref_count: 0,
            mother_alt_count: 0,
            father_alt_count: 0,
            proband_unique: 0,
            low_complexity: true,
            not_implemented: false,
        };
        assert_eq!(
            determine_filter(&score, &config),
            vcf_io::FILTER_LOW_COMPLEXITY
        );
    }

    #[test]
    fn test_determine_filter_not_implemented() {
        let config = FilterConfig {
            kmer_size: 31,
            min_baseq: 20,
            min_mapq: 20,
            max_reads_per_locus: 200,
            min_child_alt: 3,
            max_parent_alt: 1,
            min_child_alt_ratio: 0.1,
            max_variant_size: 50,
            window: 500,
            debug_kmers: false,
        };

        let score = ScoreResult {
            child_alt_count: 0,
            child_ref_count: 0,
            mother_alt_count: 0,
            father_alt_count: 0,
            proband_unique: 0,
            low_complexity: false,
            not_implemented: true,
        };
        assert_eq!(
            determine_filter(&score, &config),
            vcf_io::FILTER_NOT_IMPLEMENTED
        );
    }

    #[test]
    fn test_determine_filter_alt_ratio() {
        let config = FilterConfig {
            kmer_size: 31,
            min_baseq: 20,
            min_mapq: 20,
            max_reads_per_locus: 200,
            min_child_alt: 3,
            max_parent_alt: 1,
            min_child_alt_ratio: 0.1,
            max_variant_size: 50,
            window: 500,
            debug_kmers: false,
        };

        // ALT ratio too low: 3 / (100 + 1) = 0.0297
        let score = ScoreResult {
            child_alt_count: 3,
            child_ref_count: 100,
            mother_alt_count: 0,
            father_alt_count: 0,
            proband_unique: 0,
            low_complexity: false,
            not_implemented: false,
        };
        assert_eq!(
            determine_filter(&score, &config),
            vcf_io::FILTER_NO_CHILD_ALT
        );
    }
}
