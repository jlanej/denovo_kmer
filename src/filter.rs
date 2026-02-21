use anyhow::{Context, Result};
use log::info;
use rustc_hash::FxHashSet;

use crate::counter::{extract_reads_kmers_at_variant, scan_parent_bam};
use crate::metrics::{write_metrics, Metrics};
use crate::vcf_io::{read_vcf, write_vcf};

/// Configuration for the de novo k-mer filter pipeline.
#[derive(Debug, Clone)]
pub struct FilterConfig {
    pub child_bam: String,
    pub mother_bam: String,
    pub father_bam: String,
    pub ref_fasta: String,
    pub vcf_path: String,
    pub output_path: String,
    pub metrics_path: Option<String>,
    pub kmer_size: usize,
    pub min_baseq: u8,
    pub min_mapq: u8,
    pub debug_kmers: bool,
    pub threads: usize,
}

/// Run the de novo k-mer filter pipeline.
pub fn run_filter(config: &FilterConfig) -> Result<()> {
    info!("Reading VCF: {}", config.vcf_path);
    let mut variants =
        read_vcf(&config.vcf_path).with_context(|| "Failed to read input VCF")?;

    info!("Loaded {} candidate variants", variants.len());

    let mut metrics = Metrics::new(config.kmer_size);
    metrics.total_variants = variants.len();

    // Step 1: Extract child k-mers at each variant position
    info!("Extracting child k-mers at variant positions...");
    let mut per_variant_read_kmers: Vec<Vec<FxHashSet<u128>>> = Vec::with_capacity(variants.len());

    for variant in &variants {
        let read_kmers = extract_reads_kmers_at_variant(
            &config.child_bam,
            &variant.chrom,
            variant.pos,
            config.kmer_size,
            config.min_baseq,
            config.min_mapq,
            config.threads,
        )
        .with_context(|| {
            format!(
                "Failed to extract child k-mers at {}:{}",
                variant.chrom, variant.pos
            )
        })?;

        per_variant_read_kmers.push(read_kmers);
    }

    // Step 2: Collect all child k-mers into a single set for parent scanning
    let mut all_child_kmers = FxHashSet::default();
    for read_kmers in &per_variant_read_kmers {
        for read_set in read_kmers {
            for kmer in read_set {
                all_child_kmers.insert(*kmer);
            }
        }
    }

    metrics.total_child_kmers = all_child_kmers.len();
    info!(
        "Collected {} unique child k-mers across all variants",
        all_child_kmers.len()
    );

    // Step 3: Scan parent BAMs for child k-mers (whole-file scan)
    info!("Scanning mother BAM: {}", config.mother_bam);
    let mother_found = scan_parent_bam(
        &config.mother_bam,
        &all_child_kmers,
        config.kmer_size,
        config.threads,
    )
    .with_context(|| "Failed to scan mother BAM")?;

    info!("Scanning father BAM: {}", config.father_bam);
    let father_found = scan_parent_bam(
        &config.father_bam,
        &all_child_kmers,
        config.kmer_size,
        config.threads,
    )
    .with_context(|| "Failed to scan father BAM")?;

    // Union of parent k-mers
    let parent_kmers: FxHashSet<u128> = mother_found.union(&father_found).copied().collect();
    metrics.parent_found_kmers = parent_kmers.len();

    info!(
        "Found {} child k-mers in parents ({} mother, {} father)",
        parent_kmers.len(),
        mother_found.len(),
        father_found.len()
    );

    // Step 4: Count proband-unique reads for each variant
    info!("Computing proband-unique read counts...");
    for (i, variant) in variants.iter_mut().enumerate() {
        let read_kmers = &per_variant_read_kmers[i];

        let mut unique_read_count = 0i32;
        for read_set in read_kmers {
            // A read is proband-unique if at least one of its k-mers
            // overlapping the variant position is absent from both parents
            let has_unique = read_set.iter().any(|kmer| !parent_kmers.contains(kmer));
            if has_unique {
                unique_read_count += 1;
            }
        }

        variant.proband_unique = Some(unique_read_count);

        if unique_read_count > 0 {
            metrics.variants_with_unique_kmers += 1;
        }

        if config.debug_kmers {
            let total_kmers: usize = read_set_total(read_kmers);
            let unique_kmers: usize = read_kmers
                .iter()
                .flat_map(|s| s.iter())
                .filter(|k| !parent_kmers.contains(k))
                .count();
            info!(
                "  {}:{} reads={} unique_reads={} kmers={} unique_kmers={}",
                variant.chrom,
                variant.pos + 1,
                read_kmers.len(),
                unique_read_count,
                total_kmers,
                unique_kmers,
            );
        }
    }

    // Step 5: Write output VCF
    info!("Writing output VCF: {}", config.output_path);
    write_vcf(&config.output_path, &config.vcf_path, &variants)
        .with_context(|| "Failed to write output VCF")?;

    // Step 6: Write metrics if requested
    if let Some(ref metrics_path) = config.metrics_path {
        info!("Writing metrics: {}", metrics_path);
        write_metrics(metrics_path, &metrics)
            .with_context(|| "Failed to write metrics")?;
    }

    info!(
        "Done. {}/{} variants have proband-unique k-mers",
        metrics.variants_with_unique_kmers, metrics.total_variants
    );

    Ok(())
}

fn read_set_total(read_kmers: &[FxHashSet<u128>]) -> usize {
    read_kmers.iter().map(|s| s.len()).sum()
}
