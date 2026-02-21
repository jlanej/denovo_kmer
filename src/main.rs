use anyhow::Result;
use clap::Parser;
use log::info;

use denovo_kmer::filter::{run_filter, FilterConfig};

/// De novo variant curation via k-mer analysis.
///
/// Annotates candidate de novo variants with proband-unique k-mer counts
/// by comparing child reads against whole-file parent BAM scans.
#[derive(Parser, Debug)]
#[command(name = "kmer-denovo", version, about)]
struct Cli {
    /// Child BAM/CRAM file (indexed)
    #[arg(long)]
    child: String,

    /// Mother BAM/CRAM file (index not required; scanned in full)
    #[arg(long)]
    mother: String,

    /// Father BAM/CRAM file (index not required; scanned in full)
    #[arg(long)]
    father: String,

    /// Reference FASTA with .fai index
    #[arg(long = "ref-fasta", short = 'r')]
    ref_fasta: String,

    /// Input VCF with candidate variants
    #[arg(long)]
    vcf: String,

    /// Output annotated VCF
    #[arg(long, short = 'o')]
    output: String,

    /// Output summary metrics JSON file
    #[arg(long)]
    metrics: Option<String>,

    /// K-mer size
    #[arg(long = "kmer-size", short = 'k', default_value = "31")]
    kmer_size: usize,

    /// Minimum base quality for read k-mers
    #[arg(long = "min-baseq", default_value = "20")]
    min_baseq: u8,

    /// Minimum mapping quality for child reads
    #[arg(long = "min-mapq", default_value = "20")]
    min_mapq: u8,

    /// Enable per-variant debug output
    #[arg(long = "debug-kmers", default_value = "false")]
    debug_kmers: bool,

    /// Number of threads for BAM decompression
    #[arg(long, default_value = "4")]
    threads: usize,
}

fn main() -> Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let cli = Cli::parse();

    info!("kmer-denovo v{}", env!("CARGO_PKG_VERSION"));
    info!("K-mer size: {}", cli.kmer_size);

    let config = FilterConfig {
        child_bam: cli.child,
        mother_bam: cli.mother,
        father_bam: cli.father,
        ref_fasta: cli.ref_fasta,
        vcf_path: cli.vcf,
        output_path: cli.output,
        metrics_path: cli.metrics,
        kmer_size: cli.kmer_size,
        min_baseq: cli.min_baseq,
        min_mapq: cli.min_mapq,
        debug_kmers: cli.debug_kmers,
        threads: cli.threads,
    };

    run_filter(&config)?;

    Ok(())
}
