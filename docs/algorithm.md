# K-mer De Novo Variant Filtering: Algorithm Description

## Overview

`kmer-denovo` uses k-mer evidence from whole-genome sequencing (WGS) reads to
filter candidate de novo variants in a trio (child + mother + father). The core
idea is that a true de novo mutation in the child will produce ALT-supporting
k-mers in the child's reads but **not** in either parent's reads.

## K-mer Construction

For each candidate variant in the input VCF:

1. **Haplotype construction**: Build two haplotype sequences centered on the variant:
   - **REF haplotype** = `left_flank` + `REF allele` + `right_flank`
   - **ALT haplotype** = `left_flank` + `ALT allele` + `right_flank`
   - Flanking context is `k - 1` bases from the reference on each side.

2. **K-mer extraction**: Extract all k-mers from each haplotype.

3. **Discriminating k-mers**: Compute set differences:
   - **REF-unique k-mers** = k-mers in REF haplotype but not in ALT haplotype
   - **ALT-unique k-mers** = k-mers in ALT haplotype but not in REF haplotype
   These are the "informative" k-mers that distinguish the two alleles.

4. **Canonicalization**: All k-mers are stored as their canonical form
   (lexicographically smaller of the k-mer and its reverse complement).
   This ensures strand-independent matching.

## Read Counting

For each trio member (child, mother, father):

1. **Fetch reads** overlapping a window around the variant locus.
2. **Quality filter**: Skip reads with mapping quality below `--min-mapq`.
   Mask individual bases below `--min-baseq` as N (which breaks k-mers
   spanning those positions).
3. **Extract k-mers** from each read's sequence.
4. **Count matches** against the ALT-unique and REF-unique k-mer sets.
   Each target k-mer is counted at most once per read.

## Filtering Decision

The following thresholds are applied (all configurable via CLI):

| Filter | Condition | Default |
|--------|-----------|---------|
| `DN_KMER_OK` | Child ALT ≥ threshold AND parents ALT ≤ threshold AND ALT ratio ≥ threshold | child_alt ≥ 3, parent_alt ≤ 1, ratio ≥ 0.1 |
| `DN_KMER_NO_CHILD_ALT` | Child ALT < threshold OR ALT ratio too low | child_alt < 3 |
| `DN_KMER_PARENT_ALT` | Either parent ALT > threshold | parent_alt > 1 |
| `DN_KMER_LOW_COMPLEXITY` | All ALT-unique k-mers are low complexity | - |
| `DN_KMER_NOT_IMPLEMENTED` | Variant size exceeds max or unsupported type | size > 50 bp |

The ALT ratio is defined as:
```
ALT_ratio = child_ALT_count / (child_REF_count + 1)
```

## Low Complexity Detection

A k-mer is considered low complexity if any of:
- It is a homopolymer (all same base, e.g., `AAAAAAA`)
- It is a perfect dinucleotide repeat (e.g., `ATATATATAT`)
- More than 80% of bases are a single nucleotide

Variants where **all** ALT-unique k-mers are low complexity are flagged
`DN_KMER_LOW_COMPLEXITY`, as the evidence is unreliable.

## VCF Annotations

The output VCF includes the following INFO fields:

| Field | Type | Description |
|-------|------|-------------|
| `KMER_K` | Integer | K-mer size used |
| `KMER_CHILD_ALT` | Integer | Count of reads in child with ALT k-mers |
| `KMER_CHILD_REF` | Integer | Count of reads in child with REF k-mers |
| `KMER_MOTHER_ALT` | Integer | Count of reads in mother with ALT k-mers |
| `KMER_FATHER_ALT` | Integer | Count of reads in father with ALT k-mers |
| `KMER_STATUS` | String | Filter status label |

## Parameter Tuning Guide

- **`--kmer-size` (default: 31)**: Odd values only. Larger k gives more
  specificity but requires longer reads. Recommended range: 21–51.
- **`--min-child-alt` (default: 3)**: Minimum ALT k-mer read support in child.
  Lower values increase sensitivity but also false positives.
- **`--max-parent-alt` (default: 1)**: Maximum ALT k-mer count allowed in
  either parent. Set to 0 for strict filtering.
- **`--min-child-alt-ratio` (default: 0.1)**: Prevents calling de novos with
  very low variant allele fraction.
- **`--window` (default: 500)**: Read fetch window around each variant.
  Increase for structural variants or repetitive regions.

## Limitations

- **Variant size**: Only SNVs and small indels (≤ 50 bp by default) are
  supported. Larger structural variants get `DN_KMER_NOT_IMPLEMENTED`.
- **Repetitive regions**: K-mers in highly repetitive regions may match
  multiple genomic locations, leading to false counts.
- **Coverage dependence**: Low-coverage samples may not have enough reads
  to meet the `--min-child-alt` threshold.
- **Multiallelic sites**: Each ALT allele is processed independently.
- **Symbolic alleles**: Breakend and symbolic alleles (e.g., `<DEL>`) are skipped.
