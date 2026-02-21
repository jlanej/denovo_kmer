/// Read k-mer counting: fetch reads from BAM/CRAM around variant loci and count k-mer matches.
use ahash::AHashMap;
use log::debug;
use rust_htslib::bam::{self, Read};

use crate::kmer;

/// Count how many reads contain each of the target k-mers in a region of a BAM/CRAM file.
///
/// Returns a map from canonical k-mer -> count of reads containing that k-mer.
///
/// Parameters:
/// - `bam_path`: path to BAM/CRAM file
/// - `ref_path`: path to reference FASTA (needed for CRAM)
/// - `chrom`: chromosome/contig name
/// - `start`: 0-based start of region
/// - `end`: 0-based exclusive end of region
/// - `target_kmers`: set of canonical k-mers to look for
/// - `k`: k-mer size
/// - `min_mapq`: minimum mapping quality
/// - `min_baseq`: minimum base quality
/// - `max_reads`: maximum reads to process (downsampling guard)
pub fn count_kmers_in_reads(
    bam_path: &str,
    ref_path: &str,
    chrom: &str,
    start: i64,
    end: i64,
    target_kmers: &ahash::AHashSet<Vec<u8>>,
    k: usize,
    min_mapq: u8,
    min_baseq: u8,
    max_reads: usize,
) -> Result<AHashMap<Vec<u8>, u32>, Box<dyn std::error::Error>> {
    let mut counts: AHashMap<Vec<u8>, u32> = AHashMap::new();

    if target_kmers.is_empty() {
        return Ok(counts);
    }

    let mut reader = bam::IndexedReader::from_path(bam_path)?;
    reader.set_reference(ref_path)?;

    let tid = reader
        .header()
        .tid(chrom.as_bytes())
        .ok_or_else(|| format!("Contig '{}' not found in {}", chrom, bam_path))?;

    reader.fetch(bam::FetchDefinition::Region(
        tid as i32,
        start,
        end,
    ))?;

    let mut read_count = 0u64;
    let mut record = bam::Record::new();

    while let Some(result) = reader.read(&mut record) {
        result?;

        // Skip reads with low mapping quality
        if record.mapq() < min_mapq {
            continue;
        }
        // Skip unmapped, duplicate, secondary, supplementary
        let flags = record.flags();
        if flags & 0x904 != 0 {
            // unmapped(4) | secondary(256) | supplementary(2048)
            continue;
        }
        if flags & 0x400 != 0 {
            // duplicate
            continue;
        }

        read_count += 1;
        if read_count as usize > max_reads {
            debug!(
                "Reached max reads ({}) for {}:{}-{}",
                max_reads, chrom, start, end
            );
            break;
        }

        let seq = record.seq().as_bytes();
        let qual = record.qual().to_vec();

        // Extract k-mers from this read with quality filtering
        let read_kmers = kmer::extract_kmers_with_qual(&seq, &qual, k, min_baseq);

        // Count matches against target k-mers (count each target at most once per read)
        let mut seen_in_read: ahash::AHashSet<Vec<u8>> = ahash::AHashSet::new();
        for km in &read_kmers {
            if target_kmers.contains(km) && seen_in_read.insert(km.clone()) {
                *counts.entry(km.clone()).or_insert(0) += 1;
            }
        }
    }

    Ok(counts)
}

/// Scan the entire BAM/CRAM for target k-mers (aligner-free whole-file search).
/// This is used for parent samples to do a thorough search for alt k-mers across the whole file.
///
/// Returns a map from canonical k-mer -> total read count.
///
/// This function iterates through ALL reads in the file, which is expensive but
/// provides maximum sensitivity for detecting inherited variants.
pub fn count_kmers_whole_file(
    bam_path: &str,
    ref_path: &str,
    target_kmers: &ahash::AHashSet<Vec<u8>>,
    k: usize,
    min_mapq: u8,
    min_baseq: u8,
) -> Result<AHashMap<Vec<u8>, u32>, Box<dyn std::error::Error>> {
    let mut counts: AHashMap<Vec<u8>, u32> = AHashMap::new();

    if target_kmers.is_empty() {
        return Ok(counts);
    }

    let mut reader = bam::IndexedReader::from_path(bam_path)?;
    reader.set_reference(ref_path)?;

    let mut record = bam::Record::new();

    while let Some(result) = reader.read(&mut record) {
        result?;

        if record.mapq() < min_mapq {
            continue;
        }
        let flags = record.flags();
        if flags & 0xD04 != 0 {
            // unmapped(0x4) | secondary(0x100) | duplicate(0x400) | supplementary(0x800)
            continue;
        }

        let seq = record.seq().as_bytes();
        let qual = record.qual().to_vec();

        let read_kmers = kmer::extract_kmers_with_qual(&seq, &qual, k, min_baseq);

        let mut seen_in_read: ahash::AHashSet<Vec<u8>> = ahash::AHashSet::new();
        for km in &read_kmers {
            if target_kmers.contains(km) && seen_in_read.insert(km.clone()) {
                *counts.entry(km.clone()).or_insert(0) += 1;
            }
        }
    }

    Ok(counts)
}

/// Count k-mers in reads from a specific region, returning total count for the target set.
/// This is a simpler API that returns just the sum of counts for all target k-mers.
pub fn count_target_kmers_in_region(
    bam_path: &str,
    ref_path: &str,
    chrom: &str,
    start: i64,
    end: i64,
    target_kmers: &ahash::AHashSet<Vec<u8>>,
    k: usize,
    min_mapq: u8,
    min_baseq: u8,
    max_reads: usize,
) -> Result<u32, Box<dyn std::error::Error>> {
    let counts =
        count_kmers_in_reads(bam_path, ref_path, chrom, start, end, target_kmers, k, min_mapq, min_baseq, max_reads)?;
    Ok(counts.values().sum())
}
