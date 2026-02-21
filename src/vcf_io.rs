/// VCF I/O: reading candidate variants and writing annotated output.
use log::{debug, warn};
use rust_htslib::bcf::{self, Read};

/// Represents a candidate variant.
#[derive(Debug, Clone)]
pub struct Variant {
    pub chrom: String,
    /// 0-based position
    pub pos: i64,
    pub ref_allele: Vec<u8>,
    pub alt_allele: Vec<u8>,
    /// Original VCF record index for tracking
    pub record_idx: usize,
}

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

pub const INFO_KMER_K: &str = "KMER_K";
pub const INFO_PROBAND_UNIQUE: &str = "KMER_PROBAND_UNIQUE";

/// Per-variant annotation.
#[derive(Debug, Clone)]
pub struct VariantAnnotation {
    pub record_idx: usize,
    pub kmer_k: i32,
    pub proband_unique: i32,
}

/// Write an annotated VCF with KMER_PROBAND_UNIQUE INFO field.
pub fn write_annotated_vcf(
    input_vcf: &str,
    output_vcf: &str,
    annotations: &[VariantAnnotation],
) -> Result<(), Box<dyn std::error::Error>> {
    let reader = bcf::Reader::from_path(input_vcf)?;
    let input_header = reader.header().clone();

    let mut header = bcf::Header::from_template(&input_header);

    header.push_record(
        format!(
            "##INFO=<ID={},Number=1,Type=Integer,Description=\"K-mer size used\">",
            INFO_KMER_K
        )
        .as_bytes(),
    );
    header.push_record(
        format!(
            "##INFO=<ID={},Number=1,Type=Integer,Description=\"Number of proband reads with at least one unique k-mer overlapping this position (not found in parents)\">",
            INFO_PROBAND_UNIQUE
        )
        .as_bytes(),
    );

    let mut writer = bcf::Writer::from_path(output_vcf, &header, false, bcf::Format::Vcf)?;

    let mut reader2 = bcf::Reader::from_path(input_vcf)?;

    use ahash::AHashMap;
    let ann_map: AHashMap<usize, &VariantAnnotation> =
        annotations.iter().map(|a| (a.record_idx, a)).collect();

    let mut idx = 0usize;
    for result in reader2.records() {
        let record = result?;

        let mut new_record = writer.empty_record();
        new_record.set_rid(record.rid());
        new_record.set_pos(record.pos());

        let id = record.id();
        new_record.set_id(&id)?;

        let alleles: Vec<Vec<u8>> = record.alleles().iter().map(|a| a.to_vec()).collect();
        let allele_refs: Vec<&[u8]> = alleles.iter().map(|a| a.as_slice()).collect();
        new_record.set_alleles(&allele_refs)?;

        new_record.set_qual(record.qual());

        if let Some(ann) = ann_map.get(&idx) {
            new_record.push_info_integer(INFO_KMER_K.as_bytes(), &[ann.kmer_k])?;
            new_record.push_info_integer(INFO_PROBAND_UNIQUE.as_bytes(), &[ann.proband_unique])?;
        }

        let pass_id = writer.header().name_to_id(b"PASS")?;
        new_record.set_filters(&[&pass_id])?;

        writer.write(&new_record)?;
        idx += 1;
    }

    Ok(())
}
