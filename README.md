# kmer-denovo

Proband-unique k-mer annotation for de novo variant analysis in trio WGS data.

## Overview

`kmer-denovo` takes a VCF of candidate variants and BAM/CRAM files from a
child–mother–father trio. For each variant, it counts the number of child reads
that contain at least one k-mer overlapping the variant position that is
**not found in either parent's reads**. This count is annotated as
`KMER_PROBAND_UNIQUE`.

**Key design**: Parent samples are scanned **whole-file** in an aligner-agnostic
manner (one pass per parent), so inherited variants are detected even when
parent reads are mismapped to a different locus.

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
  --kmer-size 31
```

### Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--child` | required | Child BAM/CRAM file (indexed) |
| `--mother` | required | Mother BAM/CRAM file (indexed) |
| `--father` | required | Father BAM/CRAM file (indexed) |
| `--ref-fasta` / `-r` | required | Reference FASTA with .fai index |
| `--vcf` | required | Input VCF with candidate variants |
| `--output` / `-o` | required | Output annotated VCF |
| `--metrics` | - | Output summary metrics JSON file |
| `--kmer-size` / `-k` | 31 | K-mer size (must be odd, 11–63) |
| `--min-baseq` | 20 | Minimum base quality for read k-mers |
| `--min-mapq` | 20 | Minimum mapping quality for child reads |
| `--max-reads-per-locus` | 200 | Max reads per locus |
| `--debug-kmers` | false | Enable per-variant debug output |

## Output

### Annotated VCF

Each variant is annotated with:

- **`KMER_K`**: K-mer size used
- **`KMER_PROBAND_UNIQUE`**: Number of child reads with at least one unique
  k-mer overlapping this position (not found in either parent)

### Summary Metrics

```json
{
  "total_variants": 1000,
  "kmer_size": 31,
  "runtime_seconds": 120.5
}
```

## Algorithm

For each variant position:

1. **Extract child k-mers**: Fetch child reads overlapping the variant.
   Extract only k-mers whose genomic span includes the variant position.
   K-mers are canonicalized (lexicographically smaller of k-mer and its
   reverse complement).

2. **Scan parents (whole-file)**: All child k-mers are collected into a
   single set. Each parent's entire BAM/CRAM is scanned in one pass,
   checking every read for those k-mers. No mapping quality filter is
   applied so mismapped parent reads are still detected.

3. **Count proband-unique reads**: For each variant, count the number of
   child reads where at least one k-mer overlapping the variant position
   is absent from both parents.

## Testing

```bash
cargo test
```

## License

MIT
