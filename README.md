# kmer-denovo

De novo variant filtering using k-mer evidence from trio WGS data.

## Overview

`kmer-denovo` takes a VCF of candidate de novo variants and BAM/CRAM files from
a child–mother–father trio. For each variant, it constructs k-mers representing
the REF and ALT haplotypes, then counts those k-mers in reads from each family
member. Variants are annotated with k-mer evidence counts and assigned FILTER
tags indicating whether they are likely true de novos, inherited, or
uninformative.

**Key design**: Parent samples are scanned in an **aligner-agnostic** manner
across their entire BAM/CRAM (one pass per parent), so inherited variants are
detected even when parent reads are mismapped to a different locus. The child
is queried region-by-region since we trust its alignment at the candidate loci.

## Installation

### Requirements

- Rust 1.70+ (with cargo)
- C compiler (gcc/clang) for htslib compilation
- zlib, libbz2, liblzma, libcurl development headers

### Build

```bash
# Install system dependencies (Ubuntu/Debian)
sudo apt-get install -y build-essential zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev

# Build
cargo build --release

# The binary will be at target/release/kmer-denovo
```

## Usage

```bash
kmer-denovo \
  --child child.bam \
  --mother mother.bam \
  --father father.bam \
  --ref-fasta GRCh38.fa \
  --vcf candidates.vcf \
  --output annotated.vcf \
  --metrics summary.json \
  --threads 4 \
  --kmer-size 31
```

### Required Arguments

| Argument | Description |
|----------|-------------|
| `--child` | Child BAM/CRAM file (indexed) |
| `--mother` | Mother BAM/CRAM file (indexed) |
| `--father` | Father BAM/CRAM file (indexed) |
| `--ref-fasta` / `-r` | Reference FASTA with .fai index |
| `--vcf` | Input VCF with candidate de novo variants |
| `--output` / `-o` | Output annotated VCF |

### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--metrics` | - | Output summary metrics JSON file |
| `--kmer-size` / `-k` | 31 | K-mer size (must be odd, 11–63) |
| `--threads` / `-t` | 1 | Number of threads |
| `--min-baseq` | 20 | Minimum base quality for read k-mers |
| `--min-mapq` | 20 | Minimum mapping quality for reads |
| `--max-reads-per-locus` | 200 | Max reads per locus per sample |
| `--min-child-alt` | 3 | Minimum child ALT k-mer count |
| `--max-parent-alt` | 1 | Maximum parent ALT k-mer count |
| `--min-child-alt-ratio` | 0.1 | Minimum ALT/(REF+1) ratio |
| `--max-variant-size` | 50 | Maximum variant size (bp) |
| `--window` | 500 | Read fetch window around variant |
| `--debug-kmers` | false | Enable per-variant debug output |

## Output

### Annotated VCF

The output VCF contains all input variants with added:

- **FILTER tags**: `DN_KMER_OK`, `DN_KMER_NO_CHILD_ALT`, `DN_KMER_PARENT_ALT`,
  `DN_KMER_LOW_COMPLEXITY`, `DN_KMER_NOT_IMPLEMENTED`
- **INFO fields**: `KMER_K`, `KMER_CHILD_ALT`, `KMER_CHILD_REF`,
  `KMER_MOTHER_ALT`, `KMER_FATHER_ALT`, `KMER_STATUS`

### Summary Metrics

JSON file with counts by filter status and runtime:

```json
{
  "total_variants": 1000,
  "dn_kmer_ok": 42,
  "dn_kmer_no_child_alt": 800,
  "dn_kmer_parent_alt": 100,
  "dn_kmer_low_complexity": 50,
  "dn_kmer_not_implemented": 8,
  "kmer_size": 31,
  "runtime_seconds": 120.5
}
```

## Algorithm

See [docs/algorithm.md](docs/algorithm.md) for detailed algorithm description,
threshold definitions, and parameter tuning guidance.

## Testing

```bash
# Run all tests
cargo test

# Run only unit tests
cargo test --lib

# Run integration tests
cargo test --test integration_test
```

## License

MIT
