/// VCF I/O: reading candidate variants and writing annotated output.
use log::{debug, warn};
use rust_htslib::bcf::{self, Read};

use crate::haplotype::Variant;

/// Read candidate variants from a VCF file.
pub fn read_variants(vcf_path: &str) -> Result<Vec<Variant>, Box<dyn std::error::Error>> {
    let mut reader = bcf::Reader::from_path(vcf_path)?;
    let header = reader.header().clone();
    let mut variants = Vec::new();
    let mut idx = 0usize;

    for result in reader.records() {
        let record = result?;
        let chrom = std::str::from_utf8(header.rid2name(record.rid().unwrap())?)?.to_string();
        let pos = record.pos(); // 0-based

        let alleles: Vec<Vec<u8>> = record.alleles().iter().map(|a| a.to_vec()).collect();
        if alleles.len() < 2 {
            warn!("Skipping record at {}:{} with < 2 alleles", chrom, pos + 1);
            continue;
        }

        let ref_allele = alleles[0].clone();

        // Process each ALT allele
        for alt in &alleles[1..] {
            // Skip symbolic alleles
            if alt.starts_with(b"<") || alt.contains(&b'[') || alt.contains(&b']') {
                debug!(
                    "Skipping symbolic allele at {}:{}: {}",
                    chrom,
                    pos + 1,
                    String::from_utf8_lossy(alt)
                );
                continue;
            }

            variants.push(Variant {
                chrom: chrom.clone(),
                pos,
                ref_allele: ref_allele.clone(),
                alt_allele: alt.clone(),
                record_idx: idx,
            });
        }
        idx += 1;
    }

    Ok(variants)
}

/// FILTER and INFO field names used by this tool.
pub const FILTER_OK: &str = "DN_KMER_OK";
pub const FILTER_NO_CHILD_ALT: &str = "DN_KMER_NO_CHILD_ALT";
pub const FILTER_PARENT_ALT: &str = "DN_KMER_PARENT_ALT";
pub const FILTER_LOW_COMPLEXITY: &str = "DN_KMER_LOW_COMPLEXITY";
pub const FILTER_NOT_IMPLEMENTED: &str = "DN_KMER_NOT_IMPLEMENTED";

pub const INFO_KMER_K: &str = "KMER_K";
pub const INFO_CHILD_ALT: &str = "KMER_CHILD_ALT";
pub const INFO_CHILD_REF: &str = "KMER_CHILD_REF";
pub const INFO_MOTHER_ALT: &str = "KMER_MOTHER_ALT";
pub const INFO_FATHER_ALT: &str = "KMER_FATHER_ALT";
pub const INFO_STATUS: &str = "KMER_STATUS";
pub const INFO_PROBAND_UNIQUE: &str = "KMER_PROBAND_UNIQUE";

/// Per-variant annotation result.
#[derive(Debug, Clone)]
pub struct VariantAnnotation {
    pub record_idx: usize,
    pub filter: String,
    pub kmer_k: i32,
    pub child_alt: i32,
    pub child_ref: i32,
    pub mother_alt: i32,
    pub father_alt: i32,
    pub proband_unique: i32,
    pub status: String,
}

/// Write an annotated VCF from input VCF with added FILTER/INFO fields.
pub fn write_annotated_vcf(
    input_vcf: &str,
    output_vcf: &str,
    annotations: &[VariantAnnotation],
) -> Result<(), Box<dyn std::error::Error>> {
    let reader = bcf::Reader::from_path(input_vcf)?;
    let input_header = reader.header().clone();

    // Build output header with new fields
    let mut header = bcf::Header::from_template(&input_header);

    // Add FILTER lines
    header.push_record(
        format!(
            "##FILTER=<ID={},Description=\"Sufficient child ALT k-mer support, no parent ALT\">",
            FILTER_OK
        )
        .as_bytes(),
    );
    header.push_record(
        format!(
            "##FILTER=<ID={},Description=\"Insufficient child ALT k-mer evidence\">",
            FILTER_NO_CHILD_ALT
        )
        .as_bytes(),
    );
    header.push_record(
        format!(
            "##FILTER=<ID={},Description=\"Parent ALT k-mers exceed threshold\">",
            FILTER_PARENT_ALT
        )
        .as_bytes(),
    );
    header.push_record(
        format!(
            "##FILTER=<ID={},Description=\"K-mers deemed non-informative due to low complexity\">",
            FILTER_LOW_COMPLEXITY
        )
        .as_bytes(),
    );
    header.push_record(
        format!(
            "##FILTER=<ID={},Description=\"Unsupported variant type or size\">",
            FILTER_NOT_IMPLEMENTED
        )
        .as_bytes(),
    );

    // Add INFO lines
    header.push_record(
        format!("##INFO=<ID={},Number=1,Type=Integer,Description=\"K-mer size used\">", INFO_KMER_K)
            .as_bytes(),
    );
    header.push_record(
        format!(
            "##INFO=<ID={},Number=1,Type=Integer,Description=\"Child ALT k-mer count\">",
            INFO_CHILD_ALT
        )
        .as_bytes(),
    );
    header.push_record(
        format!(
            "##INFO=<ID={},Number=1,Type=Integer,Description=\"Child REF k-mer count\">",
            INFO_CHILD_REF
        )
        .as_bytes(),
    );
    header.push_record(
        format!(
            "##INFO=<ID={},Number=1,Type=Integer,Description=\"Mother ALT k-mer count\">",
            INFO_MOTHER_ALT
        )
        .as_bytes(),
    );
    header.push_record(
        format!(
            "##INFO=<ID={},Number=1,Type=Integer,Description=\"Father ALT k-mer count\">",
            INFO_FATHER_ALT
        )
        .as_bytes(),
    );
    header.push_record(
        format!(
            "##INFO=<ID={},Number=1,Type=String,Description=\"K-mer filter status\">",
            INFO_STATUS
        )
        .as_bytes(),
    );
    header.push_record(
        format!(
            "##INFO=<ID={},Number=1,Type=Integer,Description=\"Number of proband-unique k-mers overlapping this position (not found in parents)\">",
            INFO_PROBAND_UNIQUE
        )
        .as_bytes(),
    );

    let mut writer = bcf::Writer::from_path(output_vcf, &header, false, bcf::Format::Vcf)?;

    // Re-open reader for records
    let mut reader2 = bcf::Reader::from_path(input_vcf)?;

    // Build annotation lookup by record index
    use ahash::AHashMap;
    let ann_map: AHashMap<usize, &VariantAnnotation> =
        annotations.iter().map(|a| (a.record_idx, a)).collect();

    let mut idx = 0usize;
    for result in reader2.records() {
        let record = result?;

        // Translate record to new header
        let mut new_record = writer.empty_record();
        // Copy core fields
        new_record.set_rid(record.rid());
        new_record.set_pos(record.pos());

        // Copy ID
        let id = record.id();
        new_record.set_id(&id)?;

        // Copy alleles
        let alleles: Vec<Vec<u8>> = record.alleles().iter().map(|a| a.to_vec()).collect();
        let allele_refs: Vec<&[u8]> = alleles.iter().map(|a| a.as_slice()).collect();
        new_record.set_alleles(&allele_refs)?;

        // Copy QUAL
        new_record.set_qual(record.qual());

        if let Some(ann) = ann_map.get(&idx) {
            // Push INFO fields
            new_record.push_info_integer(INFO_KMER_K.as_bytes(), &[ann.kmer_k])?;
            new_record.push_info_integer(INFO_CHILD_ALT.as_bytes(), &[ann.child_alt])?;
            new_record.push_info_integer(INFO_CHILD_REF.as_bytes(), &[ann.child_ref])?;
            new_record.push_info_integer(INFO_MOTHER_ALT.as_bytes(), &[ann.mother_alt])?;
            new_record.push_info_integer(INFO_FATHER_ALT.as_bytes(), &[ann.father_alt])?;
            new_record.push_info_string(INFO_STATUS.as_bytes(), &[ann.status.as_bytes()])?;
            new_record.push_info_integer(INFO_PROBAND_UNIQUE.as_bytes(), &[ann.proband_unique])?;

            // Set FILTER
            let filter_id = writer.header().name_to_id(ann.filter.as_bytes())?;
            new_record.set_filters(&[&filter_id])?;
        } else {
            // No annotation for this record (should not happen in normal operation).
            // Retain PASS as default since we cannot copy filters across headers easily.
            let pass_id = writer.header().name_to_id(b"PASS")?;
            new_record.set_filters(&[&pass_id])?;
        }

        writer.write(&new_record)?;
        idx += 1;
    }

    Ok(())
}
