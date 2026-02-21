use clap::Parser;
use log::{error, info};
use std::process;

use kmer_denovo::filter;

/// Proband-unique k-mer annotation for de novo variant analysis in trio WGS data
#[derive(Parser, Debug)]
#[command(name = "kmer-denovo", version, about)]
struct Cli {
    /// Child BAM/CRAM file
    #[arg(long)]
    child: String,

    /// Mother BAM/CRAM file
    #[arg(long)]
    mother: String,

    /// Father BAM/CRAM file
    #[arg(long)]
    father: String,

    /// Reference FASTA file (with .fai index)
    #[arg(long, short = 'r')]
    ref_fasta: String,

    /// Input VCF with candidate variants
    #[arg(long)]
    vcf: String,

    /// Output annotated VCF
    #[arg(long, short = 'o')]
    output: String,

    /// Output summary metrics JSON
    #[arg(long)]
    metrics: Option<String>,

    /// K-mer size (must be odd)
    #[arg(long, short = 'k', default_value = "31")]
    kmer_size: usize,

    /// Minimum base quality for read k-mers
    #[arg(long, default_value = "20")]
    min_baseq: u8,

    /// Minimum mapping quality for child reads
    #[arg(long, default_value = "20")]
    min_mapq: u8,

    /// Maximum reads per locus
    #[arg(long, default_value = "200")]
    max_reads_per_locus: usize,

    /// Enable per-variant debug output
    #[arg(long)]
    debug_kmers: bool,
}

fn main() {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
    let cli = Cli::parse();

    if cli.kmer_size % 2 == 0 {
        error!("K-mer size must be odd, got {}", cli.kmer_size);
        process::exit(1);
    }
    if cli.kmer_size < 11 || cli.kmer_size > 63 {
        error!(
            "K-mer size must be between 11 and 63, got {}",
            cli.kmer_size
        );
        process::exit(1);
    }

    let config = filter::FilterConfig {
        kmer_size: cli.kmer_size,
        min_baseq: cli.min_baseq,
        min_mapq: cli.min_mapq,
        max_reads_per_locus: cli.max_reads_per_locus,
        debug_kmers: cli.debug_kmers,
    };

    info!("kmer-denovo v{}", env!("CARGO_PKG_VERSION"));
    info!("K={}", cli.kmer_size);
    info!("Child: {}", cli.child);
    info!("Mother: {}", cli.mother);
    info!("Father: {}", cli.father);
    info!("VCF: {}", cli.vcf);
    info!("Reference: {}", cli.ref_fasta);

    match filter::run_filter(
        &cli.child,
        &cli.mother,
        &cli.father,
        &cli.ref_fasta,
        &cli.vcf,
        &cli.output,
        &config,
    ) {
        Ok(summary) => {
            info!(
                "Complete: {} variants annotated",
                summary.total_variants
            );

            if let Some(metrics_path) = &cli.metrics {
                if let Err(e) = summary.write_json(metrics_path) {
                    error!("Failed to write metrics: {}", e);
                    process::exit(1);
                }
                info!("Metrics written to {}", metrics_path);
            }
        }
        Err(e) => {
            error!("Pipeline failed: {}", e);
            process::exit(1);
        }
    }
}
