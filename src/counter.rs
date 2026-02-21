/// Read k-mer extraction and counting from BAM/CRAM files.
use ahash::{AHashMap, AHashSet};
use log::debug;
use rust_htslib::bam::{self, Read};

use crate::kmer;

/// Extract per-read k-mer sets from child reads at a variant position.
///
/// For each read overlapping the variant, extracts only the canonical k-mers
/// whose genomic span includes the variant position. Returns one set per read.
/// K-mers are encoded as u128 values using 2-bit encoding (supports k up to 64).
pub fn extract_reads_kmers_at_variant(
    bam_path: &str,
    ref_path: &str,
    chrom: &str,
    var_pos: i64,
    ref_allele_len: usize,
    k: usize,
    min_mapq: u8,
    min_baseq: u8,
    max_reads: usize,
    threads: usize,
) -> Result<Vec<AHashSet<u128>>, Box<dyn std::error::Error>> {
    let mut per_read_kmers: Vec<AHashSet<u128>> = Vec::new();

    let mut reader = bam::IndexedReader::from_path(bam_path)?;
    reader.set_reference(ref_path)?;
    reader.set_threads(threads)?;

    let tid = reader
        .header()
        .tid(chrom.as_bytes())
        .ok_or_else(|| format!("Contig '{}' not found in {}", chrom, bam_path))?;

    // Fetch reads overlapping the variant span
    let fetch_end = var_pos + ref_allele_len.max(1) as i64;
    reader.fetch(bam::FetchDefinition::Region(tid as i32, var_pos, fetch_end))?;

    let mut read_count = 0u64;
    let mut record = bam::Record::new();

    while let Some(result) = reader.read(&mut record) {
        result?;

        if record.mapq() < min_mapq {
            continue;
        }
        let flags = record.flags();
        if flags & 0x904 != 0 {
            continue;
        }
        if flags & 0x400 != 0 {
            continue;
        }

        read_count += 1;
        if read_count as usize > max_reads {
            debug!(
                "Reached max reads ({}) for {}:{}",
                max_reads, chrom, var_pos + 1
            );
            break;
        }

        let read_pos = record.pos();
        let seq = record.seq().as_bytes();
        let qual = record.qual();

        if seq.len() < k {
            continue;
        }

        // K-mer at read offset i spans genomic [read_pos + i, read_pos + i + k - 1].
        // Variant spans [var_pos, var_pos + ref_allele_len - 1].
        // Overlap when: read_pos + i + k - 1 >= var_pos AND read_pos + i <= var_pos + ref_allele_len - 1
        let min_offset = ((var_pos - read_pos - k as i64 + 1).max(0)) as usize;
        let var_end = var_pos + ref_allele_len.max(1) as i64 - 1;
        if var_end < read_pos {
            continue;
        }
        let max_offset = (var_end - read_pos) as usize;
        let max_kmer_offset = seq.len() - k;

        let mut read_kmer_set = AHashSet::new();
        for i in min_offset..=max_offset.min(max_kmer_offset) {
            let kmer_seq = &seq[i..i + k];
            let kmer_qual = &qual[i..i + k];
            if kmer_qual.iter().all(|&q| q >= min_baseq) {
                if let Some(encoded) = kmer::encode_kmer_canonical(kmer_seq, k) {
                    read_kmer_set.insert(encoded);
                }
            }
        }

        if !read_kmer_set.is_empty() {
            per_read_kmers.push(read_kmer_set);
        }
    }

    Ok(per_read_kmers)
}

/// Scan the entire BAM/CRAM for target k-mers (whole-file, aligner-agnostic).
///
/// Used for parent samples to check all reads for presence of child k-mers,
/// regardless of mapping position. No mapQ filter is applied.
///
/// Returns a map from canonical k-mer (u128) -> read count.
pub fn count_kmers_whole_file(
    bam_path: &str,
    ref_path: &str,
    target_kmers: &AHashSet<u128>,
    k: usize,
    min_baseq: u8,
    threads: usize,
) -> Result<AHashMap<u128, u32>, Box<dyn std::error::Error>> {
    let mut counts: AHashMap<u128, u32> = AHashMap::new();

    if target_kmers.is_empty() {
        return Ok(counts);
    }

    let mut reader = bam::Reader::from_path(bam_path)?;
    reader.set_reference(ref_path)?;
    reader.set_threads(threads)?;

    let mut record = bam::Record::new();
    let mut total_reads = 0u64;
    let mut seen_in_read: AHashSet<u128> = AHashSet::new();

    while let Some(result) = reader.read(&mut record) {
        result?;

        let flags = record.flags();
        // Skip secondary, QC-fail, duplicate, supplementary
        if flags & 0xF00 != 0 {
            continue;
        }

        if record.seq_len() == 0 {
            continue;
        }

        total_reads += 1;

        let seq = record.seq().as_bytes();
        let qual = record.qual();

        let read_kmers = kmer::extract_kmers_with_qual_u128(&seq, &qual, k, min_baseq);

        seen_in_read.clear();
        for km in &read_kmers {
            if target_kmers.contains(km) && seen_in_read.insert(*km) {
                *counts.entry(*km).or_insert(0) += 1;
            }
        }
    }

    debug!(
        "Whole-file scan of {}: processed {} reads, found {} distinct target k-mers",
        bam_path,
        total_reads,
        counts.len()
    );

    Ok(counts)
}

/// Scan the entire BAM/CRAM and count ALL k-mers (for building a complete database).
///
/// This is used to pre-build a parent k-mer database that can be saved to disk
/// and reused across multiple analyses, avoiding repeated whole-file scans.
pub fn count_all_kmers_whole_file(
    bam_path: &str,
    ref_path: &str,
    k: usize,
    min_baseq: u8,
    threads: usize,
) -> Result<AHashMap<u128, u32>, Box<dyn std::error::Error>> {
    let mut counts: AHashMap<u128, u32> = AHashMap::new();

    let mut reader = bam::Reader::from_path(bam_path)?;
    reader.set_reference(ref_path)?;
    reader.set_threads(threads)?;

    let mut record = bam::Record::new();
    let mut total_reads = 0u64;
    let mut seen_in_read: AHashSet<u128> = AHashSet::new();

    while let Some(result) = reader.read(&mut record) {
        result?;

        let flags = record.flags();
        if flags & 0xF00 != 0 {
            continue;
        }

        if record.seq_len() == 0 {
            continue;
        }

        total_reads += 1;

        let seq = record.seq().as_bytes();
        let qual = record.qual();

        let read_kmers = kmer::extract_kmers_with_qual_u128(&seq, &qual, k, min_baseq);

        seen_in_read.clear();
        for km in &read_kmers {
            if seen_in_read.insert(*km) {
                *counts.entry(*km).or_insert(0) += 1;
            }
        }
    }

    debug!(
        "Whole-file scan of {}: processed {} reads, counted {} distinct k-mers",
        bam_path,
        total_reads,
        counts.len()
    );

    Ok(counts)
}
