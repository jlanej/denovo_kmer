/// 2-bit encoding: A=0, C=1, G=2, T=3
#[inline]
pub fn encode_base(b: u8) -> Option<u128> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Reverse complement of a k-mer encoded as u128.
#[inline]
pub fn revcomp_u128(kmer: u128, k: usize) -> u128 {
    let mut rc: u128 = 0;
    let mut fwd = kmer;
    for _ in 0..k {
        let base = fwd & 3;
        rc = (rc << 2) | (3 - base);
        fwd >>= 2;
    }
    rc
}

/// Canonical k-mer: lexicographically smaller of k-mer and its reverse complement.
#[inline]
pub fn canonical_kmer_u128(kmer: u128, k: usize) -> u128 {
    let rc = revcomp_u128(kmer, k);
    kmer.min(rc)
}

/// Extract canonical k-mers from a sequence, filtering by base quality.
/// Returns a Vec of (offset, canonical_kmer) pairs.
pub fn extract_kmers_with_qual_u128(
    seq: &[u8],
    qual: &[u8],
    k: usize,
    min_baseq: u8,
) -> Vec<(usize, u128)> {
    if seq.len() < k || qual.len() < k {
        return Vec::new();
    }

    let mask: u128 = if k < 64 { (1u128 << (2 * k)) - 1 } else { u128::MAX };
    let mut result = Vec::new();
    let mut kmer: u128 = 0;
    let mut valid = 0usize;

    for (i, (&base, &q)) in seq.iter().zip(qual.iter()).enumerate() {
        if q < min_baseq {
            valid = 0;
            kmer = 0;
            continue;
        }
        match encode_base(base) {
            Some(bits) => {
                kmer = ((kmer << 2) | bits) & mask;
                valid += 1;
            }
            None => {
                valid = 0;
                kmer = 0;
                continue;
            }
        }
        if valid >= k {
            let offset = i + 1 - k;
            result.push((offset, canonical_kmer_u128(kmer, k)));
        }
    }
    result
}

/// Extract canonical k-mers from a sequence without quality filtering.
pub fn extract_kmers_u128(seq: &[u8], k: usize) -> Vec<(usize, u128)> {
    if seq.len() < k {
        return Vec::new();
    }
    let mask: u128 = if k < 64 { (1u128 << (2 * k)) - 1 } else { u128::MAX };
    let mut result = Vec::new();
    let mut kmer: u128 = 0;
    let mut valid = 0usize;

    for (i, &base) in seq.iter().enumerate() {
        match encode_base(base) {
            Some(bits) => {
                kmer = ((kmer << 2) | bits) & mask;
                valid += 1;
            }
            None => {
                valid = 0;
                kmer = 0;
                continue;
            }
        }
        if valid >= k {
            let offset = i + 1 - k;
            result.push((offset, canonical_kmer_u128(kmer, k)));
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_base() {
        assert_eq!(encode_base(b'A'), Some(0));
        assert_eq!(encode_base(b'C'), Some(1));
        assert_eq!(encode_base(b'G'), Some(2));
        assert_eq!(encode_base(b'T'), Some(3));
        assert_eq!(encode_base(b'a'), Some(0));
        assert_eq!(encode_base(b'N'), None);
    }

    #[test]
    fn test_revcomp() {
        // ACG -> CGT
        // A=0, C=1, G=2 → encoded as 0b_00_01_10 = 6
        let kmer = 0b_00_01_10u128; // ACG
        let rc = revcomp_u128(kmer, 3);
        // CGT: C=1, G=2, T=3 → 0b_01_10_11 = 27
        assert_eq!(rc, 0b_01_10_11u128);
    }

    #[test]
    fn test_canonical_is_min() {
        let kmer = 0b_00_01_10u128; // ACG
        let rc = revcomp_u128(kmer, 3); // CGT
        let canon = canonical_kmer_u128(kmer, 3);
        assert_eq!(canon, kmer.min(rc));
    }

    #[test]
    fn test_extract_kmers_simple() {
        let seq = b"ACGT";
        let k = 3;
        let kmers = extract_kmers_u128(seq, k);
        assert_eq!(kmers.len(), 2); // ACG at 0, CGT at 1
    }

    #[test]
    fn test_extract_kmers_with_n() {
        let seq = b"ACNGT";
        let k = 3;
        let kmers = extract_kmers_u128(seq, k);
        // N resets, so only NGT wouldn't work, and after N we need k valid bases
        // After N: G, T → only 2 valid, not enough for k=3
        assert_eq!(kmers.len(), 0);
    }

    #[test]
    fn test_extract_kmers_qual_filter() {
        let seq = b"ACGT";
        let qual = &[30u8, 30, 10, 30]; // third base below threshold
        let k = 3;
        let kmers = extract_kmers_with_qual_u128(seq, qual, k, 20);
        // Position 2 (G, q=10) resets. After reset: T alone, not enough.
        // ACG: bases 0,1,2 → qual 30,30,10 → base at index 2 has q=10 < 20, resets
        // So only first k-mer attempt fails at index 2
        assert_eq!(kmers.len(), 0);
    }

    #[test]
    fn test_extract_kmers_all_pass_qual() {
        let seq = b"ACGT";
        let qual = &[30u8, 30, 30, 30];
        let k = 3;
        let kmers = extract_kmers_with_qual_u128(seq, qual, k, 20);
        assert_eq!(kmers.len(), 2);
    }

    #[test]
    fn test_palindrome_canonical() {
        // AATT is its own reverse complement
        // A=0,A=0,T=3,T=3 → 0b00_00_11_11 = 15
        let kmer = 0b_00_00_11_11u128;
        let canon = canonical_kmer_u128(kmer, 4);
        let rc = revcomp_u128(kmer, 4);
        assert_eq!(kmer, rc); // palindrome
        assert_eq!(canon, kmer);
    }
}
