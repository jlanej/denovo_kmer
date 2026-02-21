/// Scoring, filtering, and main pipeline logic.
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

/// Result of scoring a single variant.
#[derive(Debug)]
struct ScoreResult {
    child_alt_count: u32,
    child_ref_count: u32,
    mother_alt_count: u32,
    father_alt_count: u32,
    low_complexity: bool,
    not_implemented: bool,
}

/// Score a single variant against all three trio members.
fn score_variant(
    var: &Variant,
    ref_seq: &[u8],
    child_bam: &str,
    mother_bam: &str,
    father_bam: &str,
    ref_path: &str,
    config: &FilterConfig,
) -> Result<ScoreResult, Box<dyn std::error::Error + Send + Sync>> {
    // Check variant size
    if var.size() > config.max_variant_size {
        return Ok(ScoreResult {
            child_alt_count: 0,
            child_ref_count: 0,
            mother_alt_count: 0,
            father_alt_count: 0,
            low_complexity: false,
            not_implemented: true,
        });
    }

    // Build haplotype k-mers
    let hap_kmers = haplotype::build_haplotype_kmers(ref_seq, var, config.kmer_size);

    if hap_kmers.low_complexity {
        return Ok(ScoreResult {
            child_alt_count: 0,
            child_ref_count: 0,
            mother_alt_count: 0,
            father_alt_count: 0,
            low_complexity: true,
            not_implemented: false,
        });
    }

    if hap_kmers.alt_kmers.is_empty() {
        debug!(
            "No ALT-unique k-mers for {}:{}",
            var.chrom,
            var.pos + 1
        );
        return Ok(ScoreResult {
            child_alt_count: 0,
            child_ref_count: 0,
            mother_alt_count: 0,
            father_alt_count: 0,
            low_complexity: true,
            not_implemented: false,
        });
    }

    let alt_kmer_set: AHashSet<Vec<u8>> = hap_kmers.alt_kmers.into_iter().collect();
    let ref_kmer_set: AHashSet<Vec<u8>> = hap_kmers.ref_kmers.into_iter().collect();

    // Define search region
    let start = (var.pos - config.window as i64).max(0);
    let end = var.pos + var.ref_allele.len() as i64 + config.window as i64;

    // Count child ALT k-mers
    let child_alt_count = counter::count_target_kmers_in_region(
        child_bam,
        ref_path,
        &var.chrom,
        start,
        end,
        &alt_kmer_set,
        config.kmer_size,
        config.min_mapq,
        config.min_baseq,
        config.max_reads_per_locus,
    )
    .map_err(|e| format!("Child ALT counting error: {}", e))?;

    // Count child REF k-mers
    let child_ref_count = counter::count_target_kmers_in_region(
        child_bam,
        ref_path,
        &var.chrom,
        start,
        end,
        &ref_kmer_set,
        config.kmer_size,
        config.min_mapq,
        config.min_baseq,
        config.max_reads_per_locus,
    )
    .map_err(|e| format!("Child REF counting error: {}", e))?;

    // Count parent ALT k-mers - search the region around the variant
    // For a more thorough aligner-free approach, we search the parents' reads
    // in the same region. The window parameter controls how wide we look.
    let mother_alt_count = counter::count_target_kmers_in_region(
        mother_bam,
        ref_path,
        &var.chrom,
        start,
        end,
        &alt_kmer_set,
        config.kmer_size,
        config.min_mapq,
        config.min_baseq,
        config.max_reads_per_locus,
    )
    .map_err(|e| format!("Mother ALT counting error: {}", e))?;

    let father_alt_count = counter::count_target_kmers_in_region(
        father_bam,
        ref_path,
        &var.chrom,
        start,
        end,
        &alt_kmer_set,
        config.kmer_size,
        config.min_mapq,
        config.min_baseq,
        config.max_reads_per_locus,
    )
    .map_err(|e| format!("Father ALT counting error: {}", e))?;

    if config.debug_kmers {
        debug!(
            "Variant {}:{} - child_alt={}, child_ref={}, mother_alt={}, father_alt={}, alt_kmers={}, ref_kmers={}",
            var.chrom,
            var.pos + 1,
            child_alt_count,
            child_ref_count,
            mother_alt_count,
            father_alt_count,
            alt_kmer_set.len(),
            ref_kmer_set.len(),
        );
    }

    Ok(ScoreResult {
        child_alt_count,
        child_ref_count,
        mother_alt_count,
        father_alt_count,
        low_complexity: false,
        not_implemented: false,
    })
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

    // Read variants
    info!("Reading variants from {}", vcf_path);
    let variants = vcf_io::read_variants(vcf_path)?;
    info!("Loaded {} candidate variants", variants.len());

    if variants.is_empty() {
        // Write empty annotated VCF
        vcf_io::write_annotated_vcf(vcf_path, output_path, &[])?;
        return Ok(FilterSummary {
            total_variants: 0,
            kmer_size: config.kmer_size,
            runtime_seconds: start_time.elapsed().as_secs_f64(),
            ..Default::default()
        });
    }

    // Group variants by contig for efficient reference loading
    let mut by_contig: ahash::AHashMap<String, Vec<&Variant>> = ahash::AHashMap::new();
    for var in &variants {
        by_contig
            .entry(var.chrom.clone())
            .or_default()
            .push(var);
    }

    // Process variants
    let annotations: Mutex<Vec<VariantAnnotation>> = Mutex::new(Vec::new());
    let errors: Mutex<Vec<String>> = Mutex::new(Vec::new());

    let contigs: Vec<(String, Vec<&Variant>)> = by_contig.into_iter().collect();

    contigs.into_par_iter().for_each(|(contig, contig_vars)| {
        info!(
            "Processing {} variants on {}",
            contig_vars.len(),
            contig
        );

        // Load reference for this contig
        let ref_seq = match load_reference_contig(ref_path, &contig) {
            Ok(seq) => seq,
            Err(e) => {
                errors
                    .lock()
                    .unwrap()
                    .push(format!("Failed to load contig {}: {}", contig, e));
                return;
            }
        };

        for var in &contig_vars {
            let score = match score_variant(
                var, &ref_seq, child_bam, mother_bam, father_bam, ref_path, config,
            ) {
                Ok(s) => s,
                Err(e) => {
                    warn!(
                        "Error scoring variant {}:{}: {}",
                        var.chrom,
                        var.pos + 1,
                        e
                    );
                    // Default to not implemented on error
                    ScoreResult {
                        child_alt_count: 0,
                        child_ref_count: 0,
                        mother_alt_count: 0,
                        father_alt_count: 0,
                        low_complexity: false,
                        not_implemented: true,
                    }
                }
            };

            let filter = determine_filter(&score, config);
            let ann = VariantAnnotation {
                record_idx: var.record_idx,
                filter: filter.clone(),
                kmer_k: config.kmer_size as i32,
                child_alt: score.child_alt_count as i32,
                child_ref: score.child_ref_count as i32,
                mother_alt: score.mother_alt_count as i32,
                father_alt: score.father_alt_count as i32,
                status: filter,
            };

            annotations.lock().unwrap().push(ann);
        }
    });

    let errs = errors.into_inner().unwrap();
    if !errs.is_empty() {
        let msg = errs.join("; ");
        return Err(format!("{} contig(s) failed to process: {}", errs.len(), msg).into());
    }

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

    // Write output VCF
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
            low_complexity: false,
            not_implemented: false,
        };
        assert_eq!(
            determine_filter(&score, &config),
            vcf_io::FILTER_NO_CHILD_ALT
        );
    }
}
