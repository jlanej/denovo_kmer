/// Read k-mer extraction and counting from BAM/CRAM files.
use ahash::{AHashMap, AHashSet};
use log::debug;
use rust_htslib::bam::{self, Read};

use crate::kmer;

/// Extract per-read k-mer sets from child reads at a variant position.
///
/// For each read overlapping the variant, extracts only the canonical k-mers
/// whose genomic span includes the variant position. Returns one set per read.
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
) -> Result<Vec<AHashSet<Vec<u8>>>, Box<dyn std::error::Error>> {
    let mut per_read_kmers: Vec<AHashSet<Vec<u8>>> = Vec::new();

    let mut reader = bam::IndexedReader::from_path(bam_path)?;
    reader.set_reference(ref_path)?;

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
        let qual = record.qual().to_vec();

        // Mask low-quality bases
        let masked: Vec<u8> = seq
            .iter()
            .zip(qual.iter())
            .map(|(&base, &q)| if q >= min_baseq { base } else { b'N' })
            .collect();

        if masked.len() < k {
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
        let max_kmer_offset = masked.len() - k;

        let mut read_kmer_set = AHashSet::new();
        for i in min_offset..=max_offset.min(max_kmer_offset) {
            let kmer_seq = &masked[i..i + k];
            if kmer::is_valid_dna(kmer_seq) {
                read_kmer_set.insert(kmer::canonical_kmer(kmer_seq));
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
/// Returns a map from canonical k-mer -> read count.
pub fn count_kmers_whole_file(
    bam_path: &str,
    ref_path: &str,
    target_kmers: &AHashSet<Vec<u8>>,
    k: usize,
    min_baseq: u8,
) -> Result<AHashMap<Vec<u8>, u32>, Box<dyn std::error::Error>> {
    let mut counts: AHashMap<Vec<u8>, u32> = AHashMap::new();

    if target_kmers.is_empty() {
        return Ok(counts);
    }

    let mut reader = bam::Reader::from_path(bam_path)?;
    reader.set_reference(ref_path)?;

    let mut record = bam::Record::new();
    let mut total_reads = 0u64;

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
        let qual = record.qual().to_vec();

        let read_kmers = kmer::extract_kmers_with_qual(&seq, &qual, k, min_baseq);

        let mut seen_in_read: AHashSet<Vec<u8>> = AHashSet::new();
        for km in &read_kmers {
            if target_kmers.contains(km) && seen_in_read.insert(km.clone()) {
                *counts.entry(km.clone()).or_insert(0) += 1;
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
