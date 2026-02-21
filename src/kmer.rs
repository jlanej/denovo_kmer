/// K-mer operations: encoding, canonicalization, reverse complement, complexity checks.

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
}
