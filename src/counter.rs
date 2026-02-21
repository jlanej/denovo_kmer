use anyhow::{Context, Result};
use log::debug;
use rust_htslib::bam::{self, Read as BamRead};
use rustc_hash::FxHashSet;

use crate::kmer::{extract_kmers_u128, extract_kmers_with_qual_u128};

/// Extract per-read k-mer sets from child BAM at a variant position.
/// Only k-mers whose genomic span includes the variant position are returned.
/// Returns a Vec of per-read k-mer sets.
pub fn extract_reads_kmers_at_variant(
    bam_path: &str,
    chrom: &str,
    var_pos: i64,
    k: usize,
    min_baseq: u8,
    min_mapq: u8,
    threads: usize,
) -> Result<Vec<FxHashSet<u128>>> {
    let mut bam = bam::IndexedReader::from_path(bam_path)
        .with_context(|| format!("Failed to open BAM: {}", bam_path))?;
    bam.set_threads(threads).ok();

    // Determine the target region: we need reads that could produce k-mers overlapping var_pos
    let tid = {
        let header = bam.header().clone();
        let tid = header
            .tid(chrom.as_bytes())
            .with_context(|| format!("Chromosome {} not found in BAM header", chrom))?;
        tid
    };

    let fetch_start = (var_pos - k as i64 + 1).max(0);
    let fetch_end = var_pos + k as i64;
    bam.fetch((tid, fetch_start as u64, fetch_end as u64))
        .with_context(|| format!("Failed to fetch region {}:{}-{}", chrom, fetch_start, fetch_end))?;

    let mut per_read_kmers: Vec<FxHashSet<u128>> = Vec::new();

    for result in bam.records() {
        let record = result?;

        // Skip unmapped, secondary, supplementary, duplicate, QC-fail
        if record.is_unmapped()
            || record.is_secondary()
            || record.is_supplementary()
            || record.is_duplicate()
            || record.is_quality_check_failed()
        {
            continue;
        }

        if record.mapq() < min_mapq {
            continue;
        }

        let read_pos = record.pos(); // 0-based leftmost
        let seq = record.seq().as_bytes();
        let qual = record.qual().to_vec();

        let all_kmers = extract_kmers_with_qual_u128(&seq, &qual, k, min_baseq);

        // Filter: only k-mers whose genomic span overlaps the variant position
        let mut read_variant_kmers = FxHashSet::default();
        for (offset, canonical) in &all_kmers {
            let kmer_genome_start = read_pos + *offset as i64;
            let kmer_genome_end = kmer_genome_start + k as i64 - 1;
            if kmer_genome_start <= var_pos && kmer_genome_end >= var_pos {
                read_variant_kmers.insert(*canonical);
            }
        }

        if !read_variant_kmers.is_empty() {
            per_read_kmers.push(read_variant_kmers);
        }
    }

    debug!(
        "Extracted k-mers from {} reads at {}:{}",
        per_read_kmers.len(),
        chrom,
        var_pos
    );

    Ok(per_read_kmers)
}

/// Collect all unique k-mers from the per-read k-mer sets.
pub fn collect_all_child_kmers(per_read_kmers: &[FxHashSet<u128>]) -> FxHashSet<u128> {
    let mut all = FxHashSet::default();
    for read_kmers in per_read_kmers {
        for kmer in read_kmers {
            all.insert(*kmer);
        }
    }
    all
}

/// Scan an entire parent BAM/CRAM file, returning the set of child k-mers found
/// in any parent read. No mapping quality filter is applied.
pub fn scan_parent_bam(
    bam_path: &str,
    child_kmers: &FxHashSet<u128>,
    k: usize,
    threads: usize,
) -> Result<FxHashSet<u128>> {
    if child_kmers.is_empty() {
        return Ok(FxHashSet::default());
    }

    let mut bam = bam::Reader::from_path(bam_path)
        .with_context(|| format!("Failed to open parent BAM: {}", bam_path))?;
    bam.set_threads(threads).ok();

    let mut found = FxHashSet::default();

    for result in bam.records() {
        let record = result?;

        if record.is_unmapped()
            || record.is_secondary()
            || record.is_supplementary()
            || record.is_duplicate()
            || record.is_quality_check_failed()
        {
            continue;
        }

        let seq = record.seq().as_bytes();
        let read_kmers = extract_kmers_u128(&seq, k);

        for (_, canonical) in &read_kmers {
            if child_kmers.contains(canonical) {
                found.insert(*canonical);
            }
        }

        // Early exit if all child k-mers have been found
        if found.len() == child_kmers.len() {
            break;
        }
    }

    debug!(
        "Parent scan of {}: found {}/{} child k-mers",
        bam_path,
        found.len(),
        child_kmers.len()
    );

    Ok(found)
}
