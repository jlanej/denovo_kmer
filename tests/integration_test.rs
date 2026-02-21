use std::fs;
use std::path::PathBuf;

use denovo_kmer::filter::{run_filter, FilterConfig};
use rust_htslib::bam;
use rust_htslib::bam::header::{Header as BamHeader, HeaderRecord};
use rust_htslib::bam::record::{Cigar, CigarString, Record as BamRecord};
use rust_htslib::bcf;

/// Create a unique temporary directory for this test run.
fn test_dir(name: &str) -> PathBuf {
    let dir = PathBuf::from(format!(
        "/tmp/denovo_kmer_test_{}_{name}",
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    dir
}

/// Build a BAM header with a single reference sequence (chr1, length 1000).
fn bam_header() -> BamHeader {
    let mut header = BamHeader::new();
    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.6");
    hd.push_tag(b"SO", "coordinate");
    header.push_record(&hd);

    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", "chr1");
    sq.push_tag(b"LN", 1000);
    header.push_record(&sq);
    header
}

/// Write a BAM file containing the given reads.
/// Each element is (sequence, 0-based start position).
fn write_bam(path: &str, reads: &[(&[u8], i64)]) {
    let header = bam_header();
    let mut writer = bam::Writer::from_path(path, &header, bam::Format::Bam).unwrap();

    for (i, (seq, pos)) in reads.iter().enumerate() {
        let mut rec = BamRecord::new();
        let cigar = CigarString(vec![Cigar::Match(seq.len() as u32)]);
        let qual = vec![30u8; seq.len()];
        rec.set(format!("read{i}").as_bytes(), Some(&cigar), seq, &qual);
        rec.set_flags(0); // mapped, primary, forward strand
        rec.set_pos(*pos);
        rec.set_tid(0); // chr1
        rec.set_mapq(60);
        writer.write(&rec).unwrap();
    }
}

/// Write a minimal VCF with one variant record.
fn write_vcf(path: &str, chrom: &str, pos: i64, ref_allele: &[u8], alt_allele: &[u8]) {
    let mut header = bcf::Header::new();
    header.push_record(format!("##contig=<ID={chrom},length=1000>").as_bytes());

    let mut writer =
        bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
    let mut record = writer.empty_record();
    let rid = writer.header().name2rid(chrom.as_bytes()).unwrap();
    record.set_rid(Some(rid));
    record.set_pos(pos);
    record.set_alleles(&[ref_allele, alt_allele]).unwrap();
    writer.write(&record).unwrap();
}

/// Parse `KMER_PROBAND_UNIQUE=<int>` from a plain-text VCF INFO column.
fn parse_proband_unique(output_path: &str) -> Vec<i32> {
    let content = fs::read_to_string(output_path).unwrap();
    content
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            assert!(fields.len() >= 8, "VCF record has too few columns");
            let info = fields[7];
            let val_str = info
                .split("KMER_PROBAND_UNIQUE=")
                .nth(1)
                .unwrap_or_else(|| panic!("KMER_PROBAND_UNIQUE missing in INFO: {info}"));
            val_str
                .split(';')
                .next()
                .unwrap()
                .parse::<i32>()
                .unwrap()
        })
        .collect()
}

#[test]
fn test_pipeline_detects_denovo_unique_kmers() {
    // ---- set up paths ----
    let dir = test_dir("e2e");
    let child_bam = dir.join("child.bam");
    let mother_bam = dir.join("mother.bam");
    let father_bam = dir.join("father.bam");
    let vcf_in = dir.join("input.vcf");
    let vcf_out = dir.join("output.vcf");
    let ref_fa = dir.join("ref.fa");

    // Dummy reference FASTA (the pipeline requires the path but does not read it).
    fs::write(&ref_fa, ">chr1\nACGTACGT\n").unwrap();

    // ---- design synthetic reads ----
    // Variant: chr1 pos 100 (0-based internally; VCF uses 1-based = 101).
    //
    // Read layout (start pos 95, length 15):
    //   offset:  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14
    //   genome: 95 96 97 98 99 100 101 102 103 104 …
    //   child:   A  C  G  T  A  T   T   C   G   T  A  A  C  G  T
    //   parent:  A  C  G  T  A  A   T   C   G   T  A  A  C  G  T
    //
    // With k = 5 the five k-mers overlapping position 100 (offsets 1-5)
    // all contain the variant base and thus differ between child and
    // parents.  The pipeline should therefore flag the child read as
    // proband-unique (KMER_PROBAND_UNIQUE ≥ 1).

    let child_seq: &[u8] = b"ACGTATTCGTAACGT"; // T at pos 100
    let parent_seq: &[u8] = b"ACGTAATCGTAACGT"; // A at pos 100

    // ---- write BAM files ----
    write_bam(child_bam.to_str().unwrap(), &[(child_seq, 95)]);
    // Child BAM needs an index for region-based fetching.
    bam::index::build(
        child_bam.to_str().unwrap(),
        None::<&str>,
        bam::index::Type::Bai,
        1,
    )
    .unwrap();

    write_bam(mother_bam.to_str().unwrap(), &[(parent_seq, 95)]);
    write_bam(father_bam.to_str().unwrap(), &[(parent_seq, 95)]);

    // ---- write input VCF ----
    write_vcf(vcf_in.to_str().unwrap(), "chr1", 100, b"A", b"T");

    // ---- run pipeline ----
    let config = FilterConfig {
        child_bam: child_bam.to_string_lossy().into_owned(),
        mother_bam: mother_bam.to_string_lossy().into_owned(),
        father_bam: father_bam.to_string_lossy().into_owned(),
        ref_fasta: ref_fa.to_string_lossy().into_owned(),
        vcf_path: vcf_in.to_string_lossy().into_owned(),
        output_path: vcf_out.to_string_lossy().into_owned(),
        metrics_path: None,
        kmer_size: 5,
        min_baseq: 20,
        min_mapq: 20,
        debug_kmers: false,
        threads: 1,
    };

    run_filter(&config).unwrap();

    // ---- verify ----
    let counts = parse_proband_unique(vcf_out.to_str().unwrap());
    assert_eq!(counts.len(), 1, "expected exactly one variant in output");
    assert!(
        counts[0] > 0,
        "expected KMER_PROBAND_UNIQUE > 0, got {}",
        counts[0]
    );

    // Cleanup
    let _ = fs::remove_dir_all(&dir);
}
