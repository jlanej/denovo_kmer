# Proband-Unique K-mer Annotation Algorithm

## Overview

For each candidate variant, `kmer-denovo` counts the number of child reads
containing at least one k-mer overlapping the variant position that is absent
from both parents. This reference-free approach is similar to the strategy
used by RUFUS for identifying de novo mutations through k-mer set difference.

## Pipeline Steps

### 1. Extract Child K-mers

For each variant at position `var_pos` with reference allele length `ref_len`:

- Fetch child reads whose alignment overlaps `[var_pos, var_pos + ref_len)`.
- For each read at alignment position `read_pos`, extract k-mers at offsets
  `i` where the k-mer spans the variant:
  `read_pos + i + k - 1 >= var_pos` AND `read_pos + i <= var_pos + ref_len - 1`
- Apply base quality filtering: bases below `--min-baseq` are masked as N,
  breaking any k-mer that spans them.
- Each k-mer is stored in canonical form (lexicographically smaller of the
  k-mer and its reverse complement).
- Result: one set of k-mers per read, per variant.

### 2. Scan Parents (Whole-File)

- All child k-mers from all variants are collected into a single hash set.
- Each parent BAM/CRAM is scanned in **one pass** using a sequential reader.
- No mapping quality filter is applied (aligner-agnostic).
- For each primary, non-duplicate read, k-mers are extracted and checked
  against the target set.
- This catches inherited variants even when parent reads are mismapped.

### 3. Count Proband-Unique Reads

For each variant, count child reads where at least one k-mer overlapping
the variant position is **not found in either parent**.

This count is annotated as `KMER_PROBAND_UNIQUE` in the output VCF.

## Read Filtering

| Sample | Filter | Notes |
|--------|--------|-------|
| Child | mapQ ≥ `--min-mapq` | Trust child alignment at candidate loci |
| Child | Skip unmapped, secondary, supplementary, duplicate | Standard flags |
| Parent | No mapQ filter | Aligner-agnostic whole-file scan |
| Parent | Skip secondary, QC-fail, duplicate, supplementary | Standard flags |
| All | Base quality ≥ `--min-baseq` | Low-quality bases masked as N |

## K-mer Canonicalization

All k-mers use canonical form: the lexicographically smaller of the k-mer
and its reverse complement. This ensures strand-independent matching between
child and parent reads.
