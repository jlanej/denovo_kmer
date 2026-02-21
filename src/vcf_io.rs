use anyhow::{Context, Result};
use rust_htslib::bcf::{self, Read as BcfRead};

/// Represents a variant from the input VCF.
#[derive(Debug, Clone)]
pub struct Variant {
    pub chrom: String,
    pub pos: i64, // 0-based
    pub id: String,
    pub ref_allele: String,
    pub alt_alleles: Vec<String>,
    pub qual: Option<f32>,
    /// Number of proband-unique reads (annotated by the pipeline)
    pub proband_unique: Option<i32>,
}

/// Read variants from a VCF file.
pub fn read_vcf(path: &str) -> Result<Vec<Variant>> {
    let mut reader = bcf::Reader::from_path(path)
        .with_context(|| format!("Failed to open VCF: {}", path))?;

    let header_view = reader.header().clone();
    let mut variants = Vec::new();

    for result in reader.records() {
        let record = result?;
        let chrom = String::from_utf8_lossy(
            header_view.rid2name(record.rid().unwrap()).unwrap(),
        )
        .to_string();

        let pos = record.pos(); // 0-based

        let id = {
            let id_bytes = record.id();
            let id_str = String::from_utf8_lossy(&id_bytes).to_string();
            if id_str == "." { String::new() } else { id_str }
        };

        let alleles: Vec<String> = record
            .alleles()
            .iter()
            .map(|a| String::from_utf8_lossy(a).to_string())
            .collect();

        let ref_allele = alleles.first().cloned().unwrap_or_default();
        let alt_alleles: Vec<String> = alleles.into_iter().skip(1).collect();

        let qual = {
            let q = record.qual();
            if q.is_nan() { None } else { Some(q) }
        };

        variants.push(Variant {
            chrom,
            pos,
            id,
            ref_allele,
            alt_alleles,
            qual,
            proband_unique: None,
        });
    }

    Ok(variants)
}

/// Write annotated variants to a VCF file.
/// Re-reads the input VCF header to build the output header with added annotations.
pub fn write_vcf(
    output_path: &str,
    input_vcf_path: &str,
    variants: &[Variant],
) -> Result<()> {
    // Re-open input VCF to get header
    let reader = bcf::Reader::from_path(input_vcf_path)
        .with_context(|| format!("Failed to re-open input VCF: {}", input_vcf_path))?;

    let mut header = bcf::Header::from_template(reader.header());
    header.push_record(
        b"##INFO=<ID=KMER_PROBAND_UNIQUE,Number=1,Type=Integer,Description=\"Number of child reads with at least one k-mer absent from both parents\">"
    );

    let mut writer = bcf::Writer::from_path(output_path, &header, true, bcf::Format::Vcf)
        .with_context(|| format!("Failed to create output VCF: {}", output_path))?;

    for variant in variants {
        let mut record = writer.empty_record();

        let rid = writer
            .header()
            .name2rid(variant.chrom.as_bytes())
            .with_context(|| format!("Chromosome {} not in output header", variant.chrom))?;
        record.set_rid(Some(rid));
        record.set_pos(variant.pos);

        if !variant.id.is_empty() {
            record.set_id(variant.id.as_bytes())?;
        }

        // Set alleles
        let mut allele_bytes: Vec<&[u8]> = Vec::new();
        allele_bytes.push(variant.ref_allele.as_bytes());
        for alt in &variant.alt_alleles {
            allele_bytes.push(alt.as_bytes());
        }
        record.set_alleles(&allele_bytes)?;

        if let Some(q) = variant.qual {
            record.set_qual(q);
        }

        // Set proband unique count
        if let Some(count) = variant.proband_unique {
            record
                .push_info_integer(b"KMER_PROBAND_UNIQUE", &[count])
                .with_context(|| "Failed to write KMER_PROBAND_UNIQUE")?;
        }

        writer.write(&record)?;
    }

    Ok(())
}
