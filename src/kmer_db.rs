/// Disk-based k-mer database inspired by Jellyfish.
///
/// Provides a compact, sorted binary file format for k-mer count databases that
/// supports memory-mapped lookups via binary search. This enables:
///
/// - **Persistence**: Pre-built parent k-mer databases can be saved to disk and
///   reused across multiple analyses, avoiding expensive repeated BAM/CRAM scans.
/// - **Large k-mer support**: u128 encoding supports k-mers up to 64 bases.
/// - **Efficient lookups**: Binary search on sorted entries gives O(log n) lookups
///   with excellent cache behavior.
///
/// # File format (little-endian)
///
/// ```text
/// [8 bytes]  magic: b"KMDB0001"
/// [4 bytes]  k-mer size (u32)
/// [4 bytes]  reserved/padding (u32, currently 0)
/// [8 bytes]  entry count (u64)
/// [N * 20 bytes] entries: sorted by kmer
///   [16 bytes] canonical k-mer (u128 LE)
///   [4 bytes]  read count (u32 LE)
/// ```
use ahash::AHashMap;
use log::info;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};

const MAGIC: &[u8; 8] = b"KMDB0001";
const HEADER_SIZE: usize = 8 + 4 + 4 + 8; // magic + k + reserved + count
const ENTRY_SIZE: usize = 16 + 4; // u128 + u32

/// A k-mer count database that supports disk persistence and efficient lookups.
///
/// Can be constructed from an in-memory AHashMap (e.g., after scanning a BAM),
/// saved to disk, and loaded back for efficient binary-search lookups.
pub struct KmerDb {
    kmer_size: usize,
    /// Sorted array of (canonical_kmer, count) pairs
    entries: Vec<(u128, u32)>,
}

impl KmerDb {
    /// Create a KmerDb from an in-memory k-mer count map.
    pub fn from_counts(kmer_size: usize, counts: AHashMap<u128, u32>) -> Self {
        let mut entries: Vec<(u128, u32)> = counts.into_iter().collect();
        entries.sort_unstable_by_key(|(k, _)| *k);
        KmerDb { kmer_size, entries }
    }

    /// Save the database to a binary file.
    pub fn save(&self, path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        // Header
        writer.write_all(MAGIC)?;
        writer.write_all(&(self.kmer_size as u32).to_le_bytes())?;
        writer.write_all(&0u32.to_le_bytes())?; // reserved
        writer.write_all(&(self.entries.len() as u64).to_le_bytes())?;

        // Entries
        for &(kmer, count) in &self.entries {
            writer.write_all(&kmer.to_le_bytes())?;
            writer.write_all(&count.to_le_bytes())?;
        }

        writer.flush()?;
        info!(
            "Saved k-mer database: {} entries, k={}, file={}",
            self.entries.len(),
            self.kmer_size,
            path
        );
        Ok(())
    }

    /// Load a database from a binary file into memory.
    pub fn load(path: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let file = File::open(path)?;
        let file_len = file.metadata()?.len() as usize;
        let mut reader = BufReader::new(file);

        // Read and validate header
        let mut magic = [0u8; 8];
        reader.read_exact(&mut magic)?;
        if &magic != MAGIC {
            return Err(format!(
                "Invalid k-mer database file: bad magic in {}",
                path
            )
            .into());
        }

        let mut buf4 = [0u8; 4];
        reader.read_exact(&mut buf4)?;
        let kmer_size = u32::from_le_bytes(buf4) as usize;

        reader.read_exact(&mut buf4)?; // reserved

        let mut buf8 = [0u8; 8];
        reader.read_exact(&mut buf8)?;
        let n_entries = u64::from_le_bytes(buf8) as usize;

        // Validate file size
        let expected_size = HEADER_SIZE + n_entries * ENTRY_SIZE;
        if file_len != expected_size {
            return Err(format!(
                "K-mer database file size mismatch: expected {} bytes, got {} bytes",
                expected_size, file_len
            )
            .into());
        }

        // Read entries
        let mut entries = Vec::with_capacity(n_entries);
        let mut kmer_buf = [0u8; 16];
        let mut count_buf = [0u8; 4];

        for _ in 0..n_entries {
            reader.read_exact(&mut kmer_buf)?;
            reader.read_exact(&mut count_buf)?;
            let kmer = u128::from_le_bytes(kmer_buf);
            let count = u32::from_le_bytes(count_buf);
            entries.push((kmer, count));
        }

        info!(
            "Loaded k-mer database: {} entries, k={}, file={}",
            entries.len(),
            kmer_size,
            path
        );

        Ok(KmerDb { kmer_size, entries })
    }

    /// Check if a k-mer is present in the database.
    #[inline]
    pub fn contains(&self, kmer: u128) -> bool {
        self.entries
            .binary_search_by_key(&kmer, |(k, _)| *k)
            .is_ok()
    }

    /// Get the count for a k-mer, or None if not present.
    #[inline]
    pub fn get(&self, kmer: u128) -> Option<u32> {
        self.entries
            .binary_search_by_key(&kmer, |(k, _)| *k)
            .ok()
            .map(|idx| self.entries[idx].1)
    }

    /// Number of entries in the database.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether the database is empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// The k-mer size used when building this database.
    pub fn kmer_size(&self) -> usize {
        self.kmer_size
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kmerdb_roundtrip() {
        let mut counts = AHashMap::new();
        counts.insert(42u128, 5u32);
        counts.insert(100u128, 10u32);
        counts.insert(1u128, 1u32);

        let db = KmerDb::from_counts(21, counts);
        assert_eq!(db.len(), 3);
        assert_eq!(db.kmer_size(), 21);

        // Lookups
        assert!(db.contains(42));
        assert!(db.contains(100));
        assert!(db.contains(1));
        assert!(!db.contains(99));
        assert_eq!(db.get(42), Some(5));
        assert_eq!(db.get(100), Some(10));
        assert_eq!(db.get(99), None);

        // Save and reload
        let dir = tempfile::TempDir::new().unwrap();
        let path = dir.path().join("test.kmdb");
        let path_str = path.to_str().unwrap();

        db.save(path_str).unwrap();
        let loaded = KmerDb::load(path_str).unwrap();

        assert_eq!(loaded.len(), 3);
        assert_eq!(loaded.kmer_size(), 21);
        assert!(loaded.contains(42));
        assert!(loaded.contains(100));
        assert!(loaded.contains(1));
        assert!(!loaded.contains(99));
        assert_eq!(loaded.get(42), Some(5));
        assert_eq!(loaded.get(1), Some(1));
    }

    #[test]
    fn test_kmerdb_empty() {
        let db = KmerDb::from_counts(31, AHashMap::new());
        assert!(db.is_empty());
        assert_eq!(db.len(), 0);
        assert!(!db.contains(42));
        assert_eq!(db.get(42), None);
    }

    #[test]
    fn test_kmerdb_large_kmer_values() {
        // Test with large u128 values (simulating k > 32)
        let mut counts = AHashMap::new();
        let large_kmer = (1u128 << 100) | 0xDEAD_BEEF;
        counts.insert(large_kmer, 3);
        counts.insert(u128::MAX - 1, 7);

        let db = KmerDb::from_counts(51, counts);
        assert!(db.contains(large_kmer));
        assert_eq!(db.get(large_kmer), Some(3));
        assert!(db.contains(u128::MAX - 1));
        assert_eq!(db.get(u128::MAX - 1), Some(7));

        // Roundtrip through disk
        let dir = tempfile::TempDir::new().unwrap();
        let path = dir.path().join("large.kmdb");
        db.save(path.to_str().unwrap()).unwrap();
        let loaded = KmerDb::load(path.to_str().unwrap()).unwrap();
        assert_eq!(loaded.kmer_size(), 51);
        assert!(loaded.contains(large_kmer));
        assert_eq!(loaded.get(u128::MAX - 1), Some(7));
    }

    #[test]
    fn test_kmerdb_invalid_magic() {
        let dir = tempfile::TempDir::new().unwrap();
        let path = dir.path().join("bad.kmdb");
        std::fs::write(&path, b"BADMAGIC").unwrap();
        assert!(KmerDb::load(path.to_str().unwrap()).is_err());
    }
}
