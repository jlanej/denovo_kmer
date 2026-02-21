use anyhow::Result;
use serde::Serialize;
use std::fs;

#[derive(Debug, Clone, Serialize)]
pub struct Metrics {
    pub total_variants: usize,
    pub variants_with_unique_kmers: usize,
    pub total_child_kmers: usize,
    pub parent_found_kmers: usize,
    pub kmer_size: usize,
}

impl Metrics {
    pub fn new(kmer_size: usize) -> Self {
        Metrics {
            total_variants: 0,
            variants_with_unique_kmers: 0,
            total_child_kmers: 0,
            parent_found_kmers: 0,
            kmer_size,
        }
    }
}

pub fn write_metrics(path: &str, metrics: &Metrics) -> Result<()> {
    let json = serde_json::to_string_pretty(metrics)?;
    fs::write(path, json)?;
    Ok(())
}
