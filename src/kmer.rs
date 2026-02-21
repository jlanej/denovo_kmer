/// K-mer operations: encoding, canonicalization, reverse complement, complexity checks.
///
/// K-mers are encoded as u128 values using 2-bit encoding (A=0, C=1, G=2, T=3),
/// supporting k-mer sizes up to 64. Canonical form is the lexicographically smaller
/// of the k-mer and its reverse complement.

/// Complement a single DNA base.
#[inline]
pub fn complement(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        _ => b'N',
    }
}

/// Reverse complement of a DNA sequence.
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

/// Canonical k-mer: the lexicographically smaller of the k-mer and its reverse complement.
/// Both inputs and output are uppercase ASCII DNA bytes.
pub fn canonical_kmer(kmer: &[u8]) -> Vec<u8> {
    let rc = reverse_complement(kmer);
    let upper: Vec<u8> = kmer.iter().map(|b| b.to_ascii_uppercase()).collect();
    if upper <= rc {
        upper
    } else {
        rc
    }
}

/// Check if a sequence contains only valid DNA bases (ACGT) and is non-empty.
pub fn is_valid_dna(seq: &[u8]) -> bool {
    !seq.is_empty()
        && seq
            .iter()
            .all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'))
}

/// Extract all k-mers of size k from a sequence, returning canonical forms.
/// Skips k-mers containing non-ACGT bases.
pub fn extract_kmers(seq: &[u8], k: usize) -> Vec<Vec<u8>> {
    if seq.len() < k {
        return Vec::new();
    }
    let mut result = Vec::with_capacity(seq.len() - k + 1);
    for i in 0..=seq.len() - k {
        let kmer = &seq[i..i + k];
        if is_valid_dna(kmer) {
            result.push(canonical_kmer(kmer));
        }
    }
    result
}

/// Check if a k-mer is low complexity (e.g., homopolymer or dinucleotide repeat).
/// Returns true if the k-mer is deemed non-informative.
pub fn is_low_complexity(kmer: &[u8]) -> bool {
    if kmer.is_empty() {
        return true;
    }

    let upper: Vec<u8> = kmer.iter().map(|b| b.to_ascii_uppercase()).collect();

    // Check homopolymer: all same base
    if upper.iter().all(|&b| b == upper[0]) {
        return true;
    }

    // Check dinucleotide repeat: pattern of length 2 repeating
    if kmer.len() >= 4 {
        let is_dinuc = upper
            .iter()
            .enumerate()
            .all(|(i, &b)| b == upper[i % 2]);
        if is_dinuc {
            return true;
        }
    }

    // Check if dominated by a single base (>80%)
    let mut counts = [0u32; 4];
    for &b in &upper {
        match b {
            b'A' => counts[0] += 1,
            b'C' => counts[1] += 1,
            b'G' => counts[2] += 1,
            b'T' => counts[3] += 1,
            _ => {}
        }
    }
    let max_count = *counts.iter().max().unwrap();
    let total = upper.len() as u32;
    if max_count * 100 > total * 80 {
        return true;
    }

    false
}

/// Extract k-mers from a read, applying base quality filtering.
/// Bases below min_baseq are masked as 'N', which causes k-mers spanning them to be skipped.
pub fn extract_kmers_with_qual(seq: &[u8], qual: &[u8], k: usize, min_baseq: u8) -> Vec<Vec<u8>> {
    let masked: Vec<u8> = seq
        .iter()
        .zip(qual.iter())
        .map(|(&base, &q)| if q >= min_baseq { base } else { b'N' })
        .collect();
    extract_kmers(&masked, k)
}

// --- Compact u128-encoded k-mer operations (k <= 64) ---
//
// Inspired by Jellyfish's 2-bit encoding: A=0, C=1, G=2, T=3.
// Using u128 allows k-mers up to 64 bases (128 bits / 2 bits per base).

/// Encode a single DNA base as a 2-bit value. Returns None for non-ACGT bases.
#[inline]
pub fn encode_base(base: u8) -> Option<u128> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Compute the reverse complement of a u128-encoded k-mer of length k.
///
/// Uses Jellyfish-style bit-parallel reversal: complement all 2-bit pairs,
/// then reverse their order within the 128-bit word.
#[inline]
pub fn revcomp_u128(val: u128, k: usize) -> u128 {
    // Complement: XOR all 2-bit pairs with 0b11
    let mut rc = !val;
    // Reverse 2-bit pairs using parallel swap at each power-of-2 stride
    rc = ((rc >> 2) & 0x33333333_33333333_33333333_33333333)
        | ((rc & 0x33333333_33333333_33333333_33333333) << 2);
    rc = ((rc >> 4) & 0x0F0F0F0F_0F0F0F0F_0F0F0F0F_0F0F0F0F)
        | ((rc & 0x0F0F0F0F_0F0F0F0F_0F0F0F0F_0F0F0F0F) << 4);
    rc = ((rc >> 8) & 0x00FF00FF_00FF00FF_00FF00FF_00FF00FF)
        | ((rc & 0x00FF00FF_00FF00FF_00FF00FF_00FF00FF) << 8);
    rc = ((rc >> 16) & 0x0000FFFF_0000FFFF_0000FFFF_0000FFFF)
        | ((rc & 0x0000FFFF_0000FFFF_0000FFFF_0000FFFF) << 16);
    rc = ((rc >> 32) & 0x00000000_FFFFFFFF_00000000_FFFFFFFF)
        | ((rc & 0x00000000_FFFFFFFF_00000000_FFFFFFFF) << 32);
    rc = (rc >> 64) | (rc << 64);
    // Shift right to align the k bases in the lower 2k bits
    rc >> (128 - 2 * k)
}

/// Compute the canonical (min of forward and revcomp) u128-encoded k-mer.
#[inline]
pub fn canonical_kmer_u128(val: u128, k: usize) -> u128 {
    val.min(revcomp_u128(val, k))
}

/// Encode a k-mer slice as a canonical u128 value. Returns None if any base is non-ACGT.
#[inline]
pub fn encode_kmer_canonical(seq: &[u8], k: usize) -> Option<u128> {
    let mut val = 0u128;
    for &b in seq {
        val = (val << 2) | encode_base(b)?;
    }
    Some(canonical_kmer_u128(val, k))
}

/// Extract canonical k-mers as u128 from a sequence using a rolling window.
/// Requires k <= 64. Skips windows containing non-ACGT bases.
pub fn extract_kmers_u128(seq: &[u8], k: usize) -> Vec<u128> {
    if seq.len() < k || k == 0 || k > 64 {
        return Vec::new();
    }
    let mask: u128 = if k == 64 { u128::MAX } else { (1u128 << (2 * k)) - 1 };
    let mut result = Vec::with_capacity(seq.len() - k + 1);
    let mut fwd: u128 = 0;
    let mut valid_run: usize = 0;

    for &base in seq {
        if let Some(enc) = encode_base(base) {
            fwd = ((fwd << 2) | enc) & mask;
            valid_run += 1;
        } else {
            valid_run = 0;
            fwd = 0;
        }
        if valid_run >= k {
            result.push(canonical_kmer_u128(fwd, k));
        }
    }
    result
}

/// Extract canonical k-mers as u128 with base quality filtering.
/// Bases below min_baseq are treated as invalid (like N), breaking the k-mer window.
pub fn extract_kmers_with_qual_u128(
    seq: &[u8],
    qual: &[u8],
    k: usize,
    min_baseq: u8,
) -> Vec<u128> {
    if seq.len() < k || k == 0 || k > 64 || seq.len() != qual.len() {
        return Vec::new();
    }
    let mask: u128 = if k == 64 { u128::MAX } else { (1u128 << (2 * k)) - 1 };
    let mut result = Vec::with_capacity(seq.len() - k + 1);
    let mut fwd: u128 = 0;
    let mut valid_run: usize = 0;

    for (i, &base) in seq.iter().enumerate() {
        if qual[i] >= min_baseq {
            if let Some(enc) = encode_base(base) {
                fwd = ((fwd << 2) | enc) & mask;
                valid_run += 1;
            } else {
                valid_run = 0;
                fwd = 0;
            }
        } else {
            valid_run = 0;
            fwd = 0;
        }
        if valid_run >= k {
            result.push(canonical_kmer_u128(fwd, k));
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_complement() {
        assert_eq!(complement(b'A'), b'T');
        assert_eq!(complement(b'T'), b'A');
        assert_eq!(complement(b'C'), b'G');
        assert_eq!(complement(b'G'), b'C');
        assert_eq!(complement(b'N'), b'N');
        assert_eq!(complement(b'a'), b'T');
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement(b"ATCG"), b"CGAT");
        assert_eq!(reverse_complement(b"GATTACA"), b"TGTAATC");
    }

    #[test]
    fn test_canonical_kmer() {
        // ACGT is a palindrome
        assert_eq!(canonical_kmer(b"ACGT"), b"ACGT");
        // Between AAAA and TTTT, AAAA is smaller
        assert_eq!(canonical_kmer(b"TTTT"), b"AAAA");
        assert_eq!(canonical_kmer(b"AAAA"), b"AAAA");
        // Between ATCG and CGAT, ATCG is smaller
        assert_eq!(canonical_kmer(b"ATCG"), b"ATCG");
        assert_eq!(canonical_kmer(b"CGAT"), b"ATCG");
    }

    #[test]
    fn test_extract_kmers() {
        let seq = b"ACGTACGT";
        let kmers = extract_kmers(seq, 4);
        assert_eq!(kmers.len(), 5); // 8 - 4 + 1

        // All should be canonical
        for k in &kmers {
            let rc = reverse_complement(k);
            assert!(k <= &rc);
        }
    }

    #[test]
    fn test_extract_kmers_with_n() {
        let seq = b"ACNGTAC";
        let kmers = extract_kmers(seq, 3);
        // Positions: ACN (skip), CNG (skip), NGT (skip), GTA, TAC
        assert_eq!(kmers.len(), 2);
    }

    #[test]
    fn test_is_low_complexity() {
        assert!(is_low_complexity(b"AAAAAAA"));
        assert!(is_low_complexity(b"TTTTTTT"));
        assert!(is_low_complexity(b"ATATATATAT"));
        assert!(is_low_complexity(b"AAAAAAAAC")); // >80% A
        assert!(!is_low_complexity(b"ACGTACGTACG"));
        assert!(!is_low_complexity(b"GATTACAGATTACA"));
    }

    #[test]
    fn test_extract_kmers_with_qual() {
        let seq = b"ACGTACGT";
        let qual = vec![30, 30, 30, 5, 30, 30, 30, 30]; // position 3 is low quality
        let kmers = extract_kmers_with_qual(seq, &qual, 4, 20);
        // Position 3 is masked => k-mers at positions 0,1,2,3 contain N
        // Only positions 4..8 = "ACGT" => 1 k-mer at position 4
        assert_eq!(kmers.len(), 1);
    }

    #[test]
    fn test_is_valid_dna() {
        assert!(is_valid_dna(b"ACGT"));
        assert!(is_valid_dna(b"acgt"));
        assert!(!is_valid_dna(b"ACGN"));
        assert!(!is_valid_dna(b""));
    }

    #[test]
    fn test_empty_sequence() {
        assert_eq!(extract_kmers(b"", 5).len(), 0);
        assert_eq!(extract_kmers(b"ACG", 5).len(), 0);
    }

    // --- u128-encoded k-mer tests ---

    #[test]
    fn test_revcomp_u128() {
        // ACGT (A=0,C=1,G=2,T=3) => 0b00_01_10_11 = 0x1B
        // revcomp(ACGT) = ACGT (palindrome)
        let fwd = 0b00_01_10_11u128; // ACGT
        assert_eq!(revcomp_u128(fwd, 4), fwd);

        // AAAA = 0b00_00_00_00 = 0
        // revcomp(AAAA) = TTTT = 0b11_11_11_11 = 0xFF
        assert_eq!(revcomp_u128(0, 4), 0xFF);
    }

    #[test]
    fn test_canonical_kmer_u128_palindrome() {
        // ACGT is a palindrome, canonical should equal forward
        let fwd = 0b00_01_10_11u128;
        assert_eq!(canonical_kmer_u128(fwd, 4), fwd);
    }

    #[test]
    fn test_extract_kmers_u128_matches_vec() {
        // Verify u128 and Vec<u8> implementations produce the same number of k-mers
        let seq = b"ACGTACGT";
        let k = 4;
        let vec_kmers = extract_kmers(seq, k);
        let u128_kmers = extract_kmers_u128(seq, k);
        assert_eq!(vec_kmers.len(), u128_kmers.len());
    }

    #[test]
    fn test_extract_kmers_u128_with_n() {
        let seq = b"ACNGTAC";
        let kmers = extract_kmers_u128(seq, 3);
        // Same as Vec version: positions ACN, CNG, NGT skipped; GTA and TAC valid
        assert_eq!(kmers.len(), 2);
    }

    #[test]
    fn test_extract_kmers_u128_empty() {
        assert_eq!(extract_kmers_u128(b"", 5).len(), 0);
        assert_eq!(extract_kmers_u128(b"ACG", 5).len(), 0);
    }

    #[test]
    fn test_extract_kmers_with_qual_u128() {
        let seq = b"ACGTACGT";
        let qual = vec![30, 30, 30, 5, 30, 30, 30, 30];
        let kmers = extract_kmers_with_qual_u128(seq, &qual, 4, 20);
        // Same as Vec version: position 3 masked, only 1 k-mer at position 4
        assert_eq!(kmers.len(), 1);
    }

    #[test]
    fn test_u128_canonical_consistency() {
        // Verify that the u128 canonical form of a sequence and its reverse complement
        // produce the same value
        let seq = b"GATTACA";
        let rc = b"TGTAATC";
        let k = 7;
        let fwd_kmers = extract_kmers_u128(seq, k);
        let rc_kmers = extract_kmers_u128(rc, k);
        assert_eq!(fwd_kmers.len(), 1);
        assert_eq!(rc_kmers.len(), 1);
        assert_eq!(fwd_kmers[0], rc_kmers[0], "Canonical form should be identical for seq and its revcomp");
    }

    #[test]
    fn test_u128_large_k() {
        // Test with k=33 which exceeds u64 range but fits in u128
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 40 bases
        let kmers = extract_kmers_u128(seq, 33);
        assert_eq!(kmers.len(), 40 - 33 + 1);
    }
}
