/// Summary metrics for the filtering run.
use serde::Serialize;
use std::fs;

#[derive(Debug, Clone, Serialize, Default)]
pub struct FilterSummary {
    pub total_variants: usize,
    pub dn_kmer_ok: usize,
    pub dn_kmer_no_child_alt: usize,
    pub dn_kmer_parent_alt: usize,
    pub dn_kmer_low_complexity: usize,
    pub dn_kmer_not_implemented: usize,
    pub kmer_size: usize,
    pub runtime_seconds: f64,
}

impl FilterSummary {
    pub fn write_json(&self, path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let json = serde_json::to_string_pretty(self)?;
        fs::write(path, json)?;
        Ok(())
    }
}
