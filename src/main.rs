use clap::Parser;
use log::{error, info};
use std::process;

mod counter;
mod filter;
mod haplotype;
mod kmer;
mod metrics;
mod vcf_io;

/// De novo variant filtering using k-mer evidence from trio WGS data
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

    /// Input VCF with candidate de novo variants
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

    /// Minimum mapping quality for reads
    #[arg(long, default_value = "20")]
    min_mapq: u8,

    /// Maximum reads to sample per locus per sample
    #[arg(long, default_value = "200")]
    max_reads_per_locus: usize,

    /// Number of threads
    #[arg(long, short = 't', default_value = "1")]
    threads: usize,

    /// Minimum child ALT k-mer count
    #[arg(long, default_value = "3")]
    min_child_alt: u32,

    /// Maximum parent ALT k-mer count
    #[arg(long, default_value = "1")]
    max_parent_alt: u32,

    /// Minimum child ALT ratio: ALT / (REF + 1)
    #[arg(long, default_value = "0.1")]
    min_child_alt_ratio: f64,

    /// Maximum variant size (bp) to process
    #[arg(long, default_value = "50")]
    max_variant_size: usize,

    /// BED file to restrict processing regions
    #[arg(long)]
    regions: Option<String>,

    /// Enable per-variant debug output
    #[arg(long)]
    debug_kmers: bool,

    /// Window size around variant for read fetching (bp)
    #[arg(long, default_value = "500")]
    window: usize,
}

fn main() {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
    let cli = Cli::parse();

    if cli.kmer_size % 2 == 0 {
        error!("K-mer size must be odd, got {}", cli.kmer_size);
        process::exit(1);
    }
    if cli.kmer_size < 11 || cli.kmer_size > 63 {
        error!("K-mer size must be between 11 and 63, got {}", cli.kmer_size);
        process::exit(1);
    }

    let config = filter::FilterConfig {
        kmer_size: cli.kmer_size,
        min_baseq: cli.min_baseq,
        min_mapq: cli.min_mapq,
        max_reads_per_locus: cli.max_reads_per_locus,
        min_child_alt: cli.min_child_alt,
        max_parent_alt: cli.max_parent_alt,
        min_child_alt_ratio: cli.min_child_alt_ratio,
        max_variant_size: cli.max_variant_size,
        window: cli.window,
        debug_kmers: cli.debug_kmers,
    };

    info!("kmer-denovo filter v{}", env!("CARGO_PKG_VERSION"));
    info!("K={}, threads={}", cli.kmer_size, cli.threads);
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
        cli.threads,
        cli.regions.as_deref(),
    ) {
        Ok(summary) => {
            info!("Filtering complete: {} variants processed", summary.total_variants);
            info!(
                "  OK={}, NO_CHILD_ALT={}, PARENT_ALT={}, LOW_COMPLEXITY={}, NOT_IMPLEMENTED={}",
                summary.dn_kmer_ok,
                summary.dn_kmer_no_child_alt,
                summary.dn_kmer_parent_alt,
                summary.dn_kmer_low_complexity,
                summary.dn_kmer_not_implemented,
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
