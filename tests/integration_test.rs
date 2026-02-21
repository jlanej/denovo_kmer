/// Integration tests for proband-unique k-mer annotation pipeline.
///
/// Tests create synthetic BAM/VCF data and verify that the pipeline correctly
/// annotates variants with the number of child reads containing at least one
/// k-mer overlapping the variant position that is absent from both parents.
use rust_htslib::bam::{self, header, Header, Writer};
use rust_htslib::bcf::{self, Read as BcfRead};
use std::io::Write;
use tempfile::TempDir;

/// Create a synthetic reference FASTA with a single contig.
fn create_reference(dir: &std::path::Path) -> String {
    let ref_path = dir.join("ref.fa");
    let fai_path = dir.join("ref.fa.fai");

    let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
               GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATT\
               CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT\
               TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA\
               ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
               GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATT\
               CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT\
               TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA\
               ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
               GATTACAGATTACAGATTACAGATTACAGATTACAGATT";

    let mut f = std::fs::File::create(&ref_path).unwrap();
    writeln!(f, ">chr1").unwrap();
    for chunk in seq.as_bytes().chunks(80) {
        f.write_all(chunk).unwrap();
        writeln!(f).unwrap();
    }
    drop(f);

    let seq_len = seq.len();
    let mut fai = std::fs::File::create(&fai_path).unwrap();
    writeln!(fai, "chr1\t{}\t6\t80\t81", seq_len).unwrap();
    drop(fai);

    ref_path.to_str().unwrap().to_string()
}

/// Create a synthetic BAM file with reads at specified positions.
fn create_bam(
    dir: &std::path::Path,
    name: &str,
    _ref_path: &str,
    reads: &[(i64, &[u8])],
) -> String {
    let bam_path = dir.join(format!("{}.bam", name));

    let mut hdr = Header::new();
    let mut sq = header::HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", "chr1");
    sq.push_tag(b"LN", 509);
    hdr.push_record(&sq);

    let mut writer = Writer::from_path(&bam_path, &hdr, bam::Format::Bam).unwrap();

    for (i, (pos, seq)) in reads.iter().enumerate() {
        let mut rec = bam::Record::new();
        let read_name = format!("read_{}", i);

        let cigar = vec![bam::record::CigarString(vec![bam::record::Cigar::Match(
            seq.len() as u32,
        )])];

        let qual: Vec<u8> = vec![40; seq.len()];

        rec.set(read_name.as_bytes(), Some(&cigar[0]), seq, &qual);
        rec.set_pos(*pos);
        rec.set_tid(0);
        rec.set_mapq(60);
        rec.set_flags(0);

        writer.write(&rec).unwrap();
    }

    drop(writer);

    bam::index::build(&bam_path, None, bam::index::Type::Bai, 1).unwrap();

    bam_path.to_str().unwrap().to_string()
}

/// Create a VCF with a candidate SNV at position 100 (0-based).
fn create_vcf(dir: &std::path::Path, ref_path: &str) -> String {
    let vcf_path = dir.join("candidates.vcf");

    let mut content = String::new();
    content.push_str("##fileformat=VCFv4.2\n");
    content.push_str(&format!("##reference={}\n", ref_path));
    content.push_str("##contig=<ID=chr1,length=509>\n");
    content.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    // SNV at position 101 (1-based) = position 100 (0-based)
    content.push_str("chr1\t101\t.\tA\tT\t100\tPASS\t.\n");

    std::fs::write(&vcf_path, &content).unwrap();
    vcf_path.to_str().unwrap().to_string()
}

/// Create reads spanning position 100, with some carrying an ALT base.
fn make_reads_with_alt(count: usize, alt_count: usize) -> Vec<(i64, Vec<u8>)> {
    let ref_seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
               GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATT\
               CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT\
               TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA\
               ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
               GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATT\
               CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT\
               TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA\
               ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
               GATTACAGATTACAGATTACAGATTACAGATTACAGATT";

    let read_len = 100;
    let var_pos = 100usize;

    let mut reads = Vec::new();

    for i in 0..count {
        let read_start = if var_pos >= 50 { var_pos - 50 } else { 0 };
        let read_end = std::cmp::min(read_start + read_len, ref_seq.len());
        let mut read_seq = ref_seq[read_start..read_end].to_vec();

        if i < alt_count {
            let offset = var_pos - read_start;
            read_seq[offset] = b'T'; // ALT allele
        }

        reads.push((read_start as i64, read_seq));
    }

    reads
}

#[test]
fn test_denovo_proband_unique_reads() {
    // Child has ALT reads, parents have only REF → KMER_PROBAND_UNIQUE > 0
    let dir = TempDir::new().unwrap();
    let ref_path = create_reference(dir.path());

    let child_reads = make_reads_with_alt(20, 10);
    let child_reads_ref: Vec<(i64, &[u8])> =
        child_reads.iter().map(|(p, s)| (*p, s.as_slice())).collect();
    let child_bam = create_bam(dir.path(), "child", &ref_path, &child_reads_ref);

    let mother_reads = make_reads_with_alt(20, 0);
    let mother_reads_ref: Vec<(i64, &[u8])> =
        mother_reads.iter().map(|(p, s)| (*p, s.as_slice())).collect();
    let mother_bam = create_bam(dir.path(), "mother", &ref_path, &mother_reads_ref);

    let father_reads = make_reads_with_alt(20, 0);
    let father_reads_ref: Vec<(i64, &[u8])> =
        father_reads.iter().map(|(p, s)| (*p, s.as_slice())).collect();
    let father_bam = create_bam(dir.path(), "father", &ref_path, &father_reads_ref);

    let vcf_path = create_vcf(dir.path(), &ref_path);
    let output_path = dir.path().join("output.vcf").to_str().unwrap().to_string();

    let config = kmer_denovo::filter::FilterConfig {
        kmer_size: 21,
        min_baseq: 20,
        min_mapq: 20,
        max_reads_per_locus: 200,
        debug_kmers: true,
        threads: 1,
    };

    let summary = kmer_denovo::filter::run_filter(
        &child_bam,
        &mother_bam,
        &father_bam,
        &ref_path,
        &vcf_path,
        &output_path,
        &config,
        None,
        None,
    )
    .unwrap();

    assert_eq!(summary.total_variants, 1);

    // Read back and verify KMER_PROBAND_UNIQUE
    let mut reader = bcf::Reader::from_path(&output_path).unwrap();
    let mut record = reader.empty_record();
    let _ = reader.read(&mut record).unwrap();
    let proband_unique = record
        .info(b"KMER_PROBAND_UNIQUE")
        .integer()
        .unwrap()
        .unwrap()[0];

    // Should equal the number of ALT reads (10), since only those reads
    // have k-mers overlapping the variant that parents don't have
    assert_eq!(
        proband_unique, 10,
        "Expected KMER_PROBAND_UNIQUE == 10 (one per ALT read), got {}",
        proband_unique
    );

    // Verify metrics
    let metrics_path = dir.path().join("metrics.json").to_str().unwrap().to_string();
    summary.write_json(&metrics_path).unwrap();
    let metrics_content = std::fs::read_to_string(&metrics_path).unwrap();
    assert!(metrics_content.contains("\"total_variants\": 1"));
}

#[test]
fn test_inherited_zero_proband_unique() {
    // Child and mother both have ALT → KMER_PROBAND_UNIQUE == 0
    let dir = TempDir::new().unwrap();
    let ref_path = create_reference(dir.path());

    let child_reads = make_reads_with_alt(20, 10);
    let child_reads_ref: Vec<(i64, &[u8])> =
        child_reads.iter().map(|(p, s)| (*p, s.as_slice())).collect();
    let child_bam = create_bam(dir.path(), "child", &ref_path, &child_reads_ref);

    let mother_reads = make_reads_with_alt(20, 10);
    let mother_reads_ref: Vec<(i64, &[u8])> =
        mother_reads.iter().map(|(p, s)| (*p, s.as_slice())).collect();
    let mother_bam = create_bam(dir.path(), "mother", &ref_path, &mother_reads_ref);

    let father_reads = make_reads_with_alt(20, 0);
    let father_reads_ref: Vec<(i64, &[u8])> =
        father_reads.iter().map(|(p, s)| (*p, s.as_slice())).collect();
    let father_bam = create_bam(dir.path(), "father", &ref_path, &father_reads_ref);

    let vcf_path = create_vcf(dir.path(), &ref_path);
    let output_path = dir.path().join("output.vcf").to_str().unwrap().to_string();

    let config = kmer_denovo::filter::FilterConfig {
        kmer_size: 21,
        min_baseq: 20,
        min_mapq: 20,
        max_reads_per_locus: 200,
        debug_kmers: false,
        threads: 1,
    };

    let summary = kmer_denovo::filter::run_filter(
        &child_bam,
        &mother_bam,
        &father_bam,
        &ref_path,
        &vcf_path,
        &output_path,
        &config,
        None,
        None,
    )
    .unwrap();

    assert_eq!(summary.total_variants, 1);

    let mut reader = bcf::Reader::from_path(&output_path).unwrap();
    let mut record = reader.empty_record();
    let _ = reader.read(&mut record).unwrap();
    let proband_unique = record
        .info(b"KMER_PROBAND_UNIQUE")
        .integer()
        .unwrap()
        .unwrap()[0];
    assert_eq!(
        proband_unique, 0,
        "Expected KMER_PROBAND_UNIQUE == 0 for inherited variant, got {}",
        proband_unique
    );
}

#[test]
fn test_mismapped_parent_detected() {
    // Child has ALT at pos 100. Mother has same ALT reads mapped to pos 300.
    // Whole-file parent scan should find the k-mers → KMER_PROBAND_UNIQUE == 0.
    let dir = TempDir::new().unwrap();
    let ref_path = create_reference(dir.path());

    let child_reads = make_reads_with_alt(20, 10);
    let child_reads_ref: Vec<(i64, &[u8])> =
        child_reads.iter().map(|(p, s)| (*p, s.as_slice())).collect();
    let child_bam = create_bam(dir.path(), "child", &ref_path, &child_reads_ref);

    // Mother: same ALT-containing reads, but mapped to pos 300
    let mother_alt_reads = make_reads_with_alt(20, 10);
    let mother_mismapped: Vec<(i64, Vec<u8>)> = mother_alt_reads
        .into_iter()
        .map(|(_pos, seq)| (300i64, seq))
        .collect();
    let mother_reads_ref: Vec<(i64, &[u8])> = mother_mismapped
        .iter()
        .map(|(p, s)| (*p, s.as_slice()))
        .collect();
    let mother_bam = create_bam(dir.path(), "mother", &ref_path, &mother_reads_ref);

    let father_reads = make_reads_with_alt(20, 0);
    let father_reads_ref: Vec<(i64, &[u8])> =
        father_reads.iter().map(|(p, s)| (*p, s.as_slice())).collect();
    let father_bam = create_bam(dir.path(), "father", &ref_path, &father_reads_ref);

    let vcf_path = create_vcf(dir.path(), &ref_path);
    let output_path = dir.path().join("output.vcf").to_str().unwrap().to_string();

    let config = kmer_denovo::filter::FilterConfig {
        kmer_size: 21,
        min_baseq: 20,
        min_mapq: 20,
        max_reads_per_locus: 200,
        debug_kmers: true,
        threads: 1,
    };

    let summary = kmer_denovo::filter::run_filter(
        &child_bam,
        &mother_bam,
        &father_bam,
        &ref_path,
        &vcf_path,
        &output_path,
        &config,
        None,
        None,
    )
    .unwrap();

    assert_eq!(summary.total_variants, 1);

    let mut reader = bcf::Reader::from_path(&output_path).unwrap();
    let mut record = reader.empty_record();
    let _ = reader.read(&mut record).unwrap();
    let proband_unique = record
        .info(b"KMER_PROBAND_UNIQUE")
        .integer()
        .unwrap()
        .unwrap()[0];
    assert_eq!(
        proband_unique, 0,
        "Expected KMER_PROBAND_UNIQUE == 0: whole-file scan should find ALT k-mers \
         in mother even though mapped to different position, got {}",
        proband_unique
    );
}

#[test]
fn test_kmer_db_persistence() {
    // Test that building and loading parent k-mer databases produces the same
    // results as the direct BAM scan approach.
    let dir = TempDir::new().unwrap();
    let ref_path = create_reference(dir.path());

    let child_reads = make_reads_with_alt(20, 10);
    let child_reads_ref: Vec<(i64, &[u8])> =
        child_reads.iter().map(|(p, s)| (*p, s.as_slice())).collect();
    let child_bam = create_bam(dir.path(), "child", &ref_path, &child_reads_ref);

    let mother_reads = make_reads_with_alt(20, 0);
    let mother_reads_ref: Vec<(i64, &[u8])> =
        mother_reads.iter().map(|(p, s)| (*p, s.as_slice())).collect();
    let mother_bam = create_bam(dir.path(), "mother", &ref_path, &mother_reads_ref);

    let father_reads = make_reads_with_alt(20, 0);
    let father_reads_ref: Vec<(i64, &[u8])> =
        father_reads.iter().map(|(p, s)| (*p, s.as_slice())).collect();
    let father_bam = create_bam(dir.path(), "father", &ref_path, &father_reads_ref);

    let vcf_path = create_vcf(dir.path(), &ref_path);

    let config = kmer_denovo::filter::FilterConfig {
        kmer_size: 21,
        min_baseq: 20,
        min_mapq: 20,
        max_reads_per_locus: 200,
        debug_kmers: false,
        threads: 1,
    };

    // First run: build databases (files don't exist yet → scan BAM and save)
    let mother_db_path = dir.path().join("mother.kmdb").to_str().unwrap().to_string();
    let father_db_path = dir.path().join("father.kmdb").to_str().unwrap().to_string();
    let output1 = dir.path().join("output1.vcf").to_str().unwrap().to_string();

    let summary1 = kmer_denovo::filter::run_filter(
        &child_bam,
        &mother_bam,
        &father_bam,
        &ref_path,
        &vcf_path,
        &output1,
        &config,
        Some(&mother_db_path),
        Some(&father_db_path),
    )
    .unwrap();

    // Verify databases were created
    assert!(std::path::Path::new(&mother_db_path).exists());
    assert!(std::path::Path::new(&father_db_path).exists());

    // Second run: load databases from disk (files exist → skip BAM scan)
    let output2 = dir.path().join("output2.vcf").to_str().unwrap().to_string();

    let summary2 = kmer_denovo::filter::run_filter(
        &child_bam,
        &mother_bam,
        &father_bam,
        &ref_path,
        &vcf_path,
        &output2,
        &config,
        Some(&mother_db_path),
        Some(&father_db_path),
    )
    .unwrap();

    assert_eq!(summary1.total_variants, summary2.total_variants);

    // Both runs should produce KMER_PROBAND_UNIQUE == 10
    let mut reader = bcf::Reader::from_path(&output2).unwrap();
    let mut record = reader.empty_record();
    let _ = reader.read(&mut record).unwrap();
    let proband_unique = record
        .info(b"KMER_PROBAND_UNIQUE")
        .integer()
        .unwrap()
        .unwrap()[0];
    assert_eq!(
        proband_unique, 10,
        "Expected KMER_PROBAND_UNIQUE == 10 from loaded database, got {}",
        proband_unique
    );
}
