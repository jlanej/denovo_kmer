/// Summary metrics for the annotation run.
use serde::Serialize;
use std::fs;

#[derive(Debug, Clone, Serialize, Default)]
pub struct FilterSummary {
    pub total_variants: usize,
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
