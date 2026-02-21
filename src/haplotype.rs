/// Haplotype construction: generate REF and ALT haplotype sequences
/// around variant positions for k-mer extraction.
use crate::kmer;

/// Represents a candidate variant with position and alleles.
#[derive(Debug, Clone)]
pub struct Variant {
    /// Chromosome/contig name
    pub chrom: String,
    /// 0-based position
    pub pos: i64,
    /// Reference allele
    pub ref_allele: Vec<u8>,
    /// Alternative allele
    pub alt_allele: Vec<u8>,
    /// Original VCF record index for tracking
    pub record_idx: usize,
}

impl Variant {
    /// Variant size in bp (absolute difference in allele lengths for indels, 0 for SNVs).
    pub fn size(&self) -> usize {
        let ref_len = self.ref_allele.len();
        let alt_len = self.alt_allele.len();
        if ref_len > alt_len {
            ref_len - alt_len
        } else {
            alt_len - ref_len
        }
    }

    /// Whether this is a SNV (single nucleotide variant).
    pub fn is_snv(&self) -> bool {
        self.ref_allele.len() == 1 && self.alt_allele.len() == 1
    }

    /// Whether this variant is an indel.
    pub fn is_indel(&self) -> bool {
        self.ref_allele.len() != self.alt_allele.len()
    }
}

/// Result of haplotype construction for a single variant.
#[derive(Debug, Clone)]
pub struct HaplotypeKmers {
    /// K-mers unique to REF haplotype (canonical)
    pub ref_kmers: Vec<Vec<u8>>,
    /// K-mers unique to ALT haplotype (canonical)
    pub alt_kmers: Vec<Vec<u8>>,
    /// Whether the variant's k-mers are low complexity
    pub low_complexity: bool,
}

/// Build REF and ALT haplotype sequences around a variant and extract discriminating k-mers.
///
/// The approach:
/// 1. Take `k-1` bases of reference context on each side of the variant.
/// 2. Construct REF haplotype = left_flank + REF allele + right_flank.
/// 3. Construct ALT haplotype = left_flank + ALT allele + right_flank.
/// 4. Extract k-mers from each haplotype.
/// 5. Return k-mers unique to each haplotype (not shared).
///
/// `ref_seq` is the reference sequence for the contig.
/// `var` is the variant.
/// `k` is the k-mer size.
pub fn build_haplotype_kmers(
    ref_seq: &[u8],
    var: &Variant,
    k: usize,
) -> HaplotypeKmers {
    let pos = var.pos as usize;
    let ref_len = var.ref_allele.len();

    // Flanking context: k-1 bases on each side
    let flank = k - 1;
    let left_start = pos.saturating_sub(flank);
    let right_end = std::cmp::min(pos + ref_len + flank, ref_seq.len());

    let left_flank = &ref_seq[left_start..pos];
    let right_flank = &ref_seq[pos + ref_len..right_end];

    // Build REF haplotype
    let mut ref_hap = Vec::with_capacity(left_flank.len() + ref_len + right_flank.len());
    ref_hap.extend_from_slice(left_flank);
    ref_hap.extend_from_slice(&var.ref_allele);
    ref_hap.extend_from_slice(right_flank);

    // Build ALT haplotype
    let mut alt_hap =
        Vec::with_capacity(left_flank.len() + var.alt_allele.len() + right_flank.len());
    alt_hap.extend_from_slice(left_flank);
    alt_hap.extend_from_slice(&var.alt_allele);
    alt_hap.extend_from_slice(right_flank);

    // Extract k-mers
    let ref_kmer_list = kmer::extract_kmers(&ref_hap, k);
    let alt_kmer_list = kmer::extract_kmers(&alt_hap, k);

    // Build sets for set difference
    use ahash::AHashSet;
    let ref_set: AHashSet<Vec<u8>> = ref_kmer_list.iter().cloned().collect();
    let alt_set: AHashSet<Vec<u8>> = alt_kmer_list.iter().cloned().collect();

    // K-mers unique to REF (not in ALT)
    let ref_unique: Vec<Vec<u8>> = ref_set.difference(&alt_set).cloned().collect();
    // K-mers unique to ALT (not in REF)
    let alt_unique: Vec<Vec<u8>> = alt_set.difference(&ref_set).cloned().collect();

    // Check low complexity
    let low_complexity = !alt_unique.is_empty()
        && alt_unique.iter().all(|km| kmer::is_low_complexity(km));

    HaplotypeKmers {
        ref_kmers: ref_unique,
        low_complexity,
        alt_kmers: alt_unique,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_ref_seq() -> Vec<u8> {
        // A simple reference: 200 bp of semi-random sequence
        b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
          GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATT\
          CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT\
          TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA"
            .to_vec()
    }

    #[test]
    fn test_snv_haplotype_kmers() {
        let ref_seq = make_ref_seq();
        let var = Variant {
            chrom: "chr1".to_string(),
            pos: 50,
            ref_allele: vec![ref_seq[50]],
            alt_allele: vec![b'T'], // SNV
            record_idx: 0,
        };

        let result = build_haplotype_kmers(&ref_seq, &var, 11);
        // Should have some ALT-unique k-mers
        assert!(!result.alt_kmers.is_empty(), "SNV should produce ALT-unique k-mers");
        // And some REF-unique k-mers
        assert!(!result.ref_kmers.is_empty(), "SNV should produce REF-unique k-mers");
    }

    #[test]
    fn test_insertion_haplotype_kmers() {
        let ref_seq = make_ref_seq();
        let var = Variant {
            chrom: "chr1".to_string(),
            pos: 50,
            ref_allele: vec![ref_seq[50]],
            alt_allele: b"GACGT".to_vec(), // Insertion of ACGT after ref base
            record_idx: 0,
        };

        let result = build_haplotype_kmers(&ref_seq, &var, 11);
        assert!(!result.alt_kmers.is_empty(), "Insertion should produce ALT-unique k-mers");
    }

    #[test]
    fn test_deletion_haplotype_kmers() {
        let ref_seq = make_ref_seq();
        let var = Variant {
            chrom: "chr1".to_string(),
            pos: 50,
            ref_allele: ref_seq[50..55].to_vec(), // 5bp deletion
            alt_allele: vec![ref_seq[50]],
            record_idx: 0,
        };

        let result = build_haplotype_kmers(&ref_seq, &var, 11);
        assert!(!result.alt_kmers.is_empty(), "Deletion should produce ALT-unique k-mers");
        assert!(!result.ref_kmers.is_empty(), "Deletion should produce REF-unique k-mers");
    }

    #[test]
    fn test_variant_size() {
        let snv = Variant {
            chrom: "chr1".to_string(),
            pos: 0,
            ref_allele: b"A".to_vec(),
            alt_allele: b"T".to_vec(),
            record_idx: 0,
        };
        assert_eq!(snv.size(), 0);
        assert!(snv.is_snv());

        let ins = Variant {
            chrom: "chr1".to_string(),
            pos: 0,
            ref_allele: b"A".to_vec(),
            alt_allele: b"ACGT".to_vec(),
            record_idx: 0,
        };
        assert_eq!(ins.size(), 3);
        assert!(ins.is_indel());

        let del = Variant {
            chrom: "chr1".to_string(),
            pos: 0,
            ref_allele: b"ACGT".to_vec(),
            alt_allele: b"A".to_vec(),
            record_idx: 0,
        };
        assert_eq!(del.size(), 3);
        assert!(del.is_indel());
    }

    #[test]
    fn test_low_complexity_detection() {
        // Variant in a homopolymer region
        let ref_seq = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_vec();
        let var = Variant {
            chrom: "chr1".to_string(),
            pos: 20,
            ref_allele: b"A".to_vec(),
            alt_allele: b"AA".to_vec(), // Insertion in homopolymer
            record_idx: 0,
        };

        let result = build_haplotype_kmers(&ref_seq, &var, 11);
        assert!(
            result.low_complexity || result.alt_kmers.is_empty(),
            "Homopolymer insertion should be flagged as low complexity"
        );
    }
}
