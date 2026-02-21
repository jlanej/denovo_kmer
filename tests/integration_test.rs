/// Integration tests using synthetic BAM/VCF data to verify end-to-end pipeline behavior.
///
/// These tests create small synthetic BAM files with reads containing known variants,
/// and verify that the filtering pipeline correctly classifies them.
use rust_htslib::bam::{self, header, Header, Writer};
use std::io::Write;
use tempfile::TempDir;

/// Create a synthetic reference FASTA with a single contig.
fn create_reference(dir: &std::path::Path) -> String {
    let ref_path = dir.join("ref.fa");
    let fai_path = dir.join("ref.fa.fai");

    // 500 bp reference sequence with mixed bases
    let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
               GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATT\
               CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT\
               TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA\
               ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
               GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATT\
               CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT\
               TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA\
               ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
               GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATT";

    let mut f = std::fs::File::create(&ref_path).unwrap();
    writeln!(f, ">chr1").unwrap();
    // Write sequence in 80-char lines
    for chunk in seq.as_bytes().chunks(80) {
        f.write_all(chunk).unwrap();
        writeln!(f).unwrap();
    }
    drop(f);

    // Create .fai index
    let seq_len = seq.len();
    let mut fai = std::fs::File::create(&fai_path).unwrap();
    // chr1\tseq_len\toffset\tbases_per_line\tbytes_per_line
    // Offset = 6 (">chr1\n" = 6 bytes)
    writeln!(fai, "chr1\t{}\t6\t80\t81", seq_len).unwrap();
    drop(fai);

    ref_path.to_str().unwrap().to_string()
}

/// Create a synthetic BAM file with reads around a specific position.
/// `alt_base` controls whether reads carry the ALT allele at position 100.
/// Returns the path to the BAM file.
fn create_bam(
    dir: &std::path::Path,
    name: &str,
    _ref_path: &str,
    reads: &[(i64, &[u8])], // (pos, sequence) pairs
) -> String {
    let bam_path = dir.join(format!("{}.bam", name));

    // Build BAM header
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

        let qual: Vec<u8> = vec![40; seq.len()]; // High quality

        rec.set(
            read_name.as_bytes(),
            Some(&cigar[0]),
            seq,
            &qual,
        );
        rec.set_pos(*pos);
        rec.set_tid(0);
        rec.set_mapq(60);
        rec.set_flags(0); // Mapped, primary

        writer.write(&rec).unwrap();
    }

    drop(writer);

    // Index the BAM (reads are written in sorted order already)
    bam::index::build(&bam_path, None, bam::index::Type::Bai, 1).unwrap();

    bam_path.to_str().unwrap().to_string()
}

/// Create a VCF with a candidate SNV at position 100 (0-based).
fn create_vcf(dir: &std::path::Path, ref_path: &str) -> String {
    let vcf_path = dir.join("candidates.vcf");

    let mut content = String::new();
    content.push_str("##fileformat=VCFv4.2\n");
    content.push_str(&format!(
        "##reference={}\n",
        ref_path
    ));
    content.push_str("##contig=<ID=chr1,length=509>\n");
    content.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    // SNV at position 101 (1-based) = position 100 (0-based)
    // Reference base at pos 100 is 'A' (from the sequence)
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
               GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATT";

    let read_len = 100;
    let var_pos = 100usize; // 0-based position of the variant in the reference

    let mut reads = Vec::new();

    for i in 0..count {
        let read_start = if var_pos >= 50 { var_pos - 50 } else { 0 };
        let read_end = std::cmp::min(read_start + read_len, ref_seq.len());
        let mut read_seq = ref_seq[read_start..read_end].to_vec();

        // For the first `alt_count` reads, introduce the ALT base
        if i < alt_count {
            let offset = var_pos - read_start;
            read_seq[offset] = b'T'; // ALT allele
        }

        reads.push((read_start as i64, read_seq));
    }

    reads
}

#[test]
fn test_denovo_snv_detected() {
    // Scenario: Child has ALT support, parents have REF only -> DN_KMER_OK
    let dir = TempDir::new().unwrap();
    let ref_path = create_reference(dir.path());

    // Child: 10 reads with ALT, 10 with REF
    let child_reads = make_reads_with_alt(20, 10);
    let child_reads_ref: Vec<(i64, &[u8])> =
        child_reads.iter().map(|(p, s)| (*p, s.as_slice())).collect();
    let child_bam = create_bam(dir.path(), "child", &ref_path, &child_reads_ref);

    // Mother: 20 reads all REF
    let mother_reads = make_reads_with_alt(20, 0);
    let mother_reads_ref: Vec<(i64, &[u8])> =
        mother_reads.iter().map(|(p, s)| (*p, s.as_slice())).collect();
    let mother_bam = create_bam(dir.path(), "mother", &ref_path, &mother_reads_ref);

    // Father: 20 reads all REF
    let father_reads = make_reads_with_alt(20, 0);
    let father_reads_ref: Vec<(i64, &[u8])> =
        father_reads.iter().map(|(p, s)| (*p, s.as_slice())).collect();
    let father_bam = create_bam(dir.path(), "father", &ref_path, &father_reads_ref);

    let vcf_path = create_vcf(dir.path(), &ref_path);
    let output_path = dir.path().join("output.vcf").to_str().unwrap().to_string();
    let metrics_path = dir.path().join("metrics.json").to_str().unwrap().to_string();

    let config = kmer_denovo::filter::FilterConfig {
        kmer_size: 21,
        min_baseq: 20,
        min_mapq: 20,
        max_reads_per_locus: 200,
        min_child_alt: 3,
        max_parent_alt: 1,
        min_child_alt_ratio: 0.1,
        max_variant_size: 50,
        window: 500,
        debug_kmers: true,
    };

    let summary = kmer_denovo::filter::run_filter(
        &child_bam,
        &mother_bam,
        &father_bam,
        &ref_path,
        &vcf_path,
        &output_path,
        &config,
        1,
        None,
    )
    .unwrap();

    assert_eq!(summary.total_variants, 1);
    assert_eq!(summary.dn_kmer_ok, 1, "Expected DN_KMER_OK for true de novo");

    // Also write and check metrics
    summary.write_json(&metrics_path).unwrap();
    let metrics_content = std::fs::read_to_string(&metrics_path).unwrap();
    assert!(metrics_content.contains("\"dn_kmer_ok\": 1"));
}

#[test]
fn test_inherited_variant_detected() {
    // Scenario: Child has ALT, mother also has ALT -> DN_KMER_PARENT_ALT
    let dir = TempDir::new().unwrap();
    let ref_path = create_reference(dir.path());

    // Child: 10 reads with ALT
    let child_reads = make_reads_with_alt(20, 10);
    let child_reads_ref: Vec<(i64, &[u8])> =
        child_reads.iter().map(|(p, s)| (*p, s.as_slice())).collect();
    let child_bam = create_bam(dir.path(), "child", &ref_path, &child_reads_ref);

    // Mother: 10 reads with ALT (inherited!)
    let mother_reads = make_reads_with_alt(20, 10);
    let mother_reads_ref: Vec<(i64, &[u8])> =
        mother_reads.iter().map(|(p, s)| (*p, s.as_slice())).collect();
    let mother_bam = create_bam(dir.path(), "mother", &ref_path, &mother_reads_ref);

    // Father: all REF
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
        min_child_alt: 3,
        max_parent_alt: 1,
        min_child_alt_ratio: 0.1,
        max_variant_size: 50,
        window: 500,
        debug_kmers: false,
    };

    let summary = kmer_denovo::filter::run_filter(
        &child_bam,
        &mother_bam,
        &father_bam,
        &ref_path,
        &vcf_path,
        &output_path,
        &config,
        1,
        None,
    )
    .unwrap();

    assert_eq!(summary.total_variants, 1);
    assert_eq!(
        summary.dn_kmer_parent_alt, 1,
        "Expected DN_KMER_PARENT_ALT for inherited variant"
    );
}

#[test]
fn test_no_child_alt_support() {
    // Scenario: Child has insufficient ALT reads -> DN_KMER_NO_CHILD_ALT
    let dir = TempDir::new().unwrap();
    let ref_path = create_reference(dir.path());

    // Child: only 1 read with ALT (below threshold of 3)
    let child_reads = make_reads_with_alt(20, 1);
    let child_reads_ref: Vec<(i64, &[u8])> =
        child_reads.iter().map(|(p, s)| (*p, s.as_slice())).collect();
    let child_bam = create_bam(dir.path(), "child", &ref_path, &child_reads_ref);

    // Parents: all REF
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
        min_child_alt: 3,
        max_parent_alt: 1,
        min_child_alt_ratio: 0.1,
        max_variant_size: 50,
        window: 500,
        debug_kmers: false,
    };

    let summary = kmer_denovo::filter::run_filter(
        &child_bam,
        &mother_bam,
        &father_bam,
        &ref_path,
        &vcf_path,
        &output_path,
        &config,
        1,
        None,
    )
    .unwrap();

    assert_eq!(summary.total_variants, 1);
    assert_eq!(
        summary.dn_kmer_no_child_alt, 1,
        "Expected DN_KMER_NO_CHILD_ALT for insufficient child support"
    );
}

#[test]
fn test_mismapped_parent_detected() {
    // Scenario: Child has ALT reads at pos 100. Mother has reads carrying the
    // same ALT k-mers, but they are *mapped to a different position* (pos 300).
    // A region-based parent search would miss these, but the whole-file
    // aligner-agnostic scan must detect them -> DN_KMER_PARENT_ALT.
    let dir = TempDir::new().unwrap();
    let ref_path = create_reference(dir.path());

    // Child: 10 reads with ALT at pos 100
    let child_reads = make_reads_with_alt(20, 10);
    let child_reads_ref: Vec<(i64, &[u8])> =
        child_reads.iter().map(|(p, s)| (*p, s.as_slice())).collect();
    let child_bam = create_bam(dir.path(), "child", &ref_path, &child_reads_ref);

    // Mother: reads containing the same ALT sequence but mapped to pos 300.
    // Generated separately from child reads, then remapped to a distant position
    // to simulate a mismapping scenario.
    let mother_alt_reads = make_reads_with_alt(20, 10);
    let mother_mismapped: Vec<(i64, Vec<u8>)> = mother_alt_reads
        .into_iter()
        .map(|(_pos, seq)| (300i64, seq)) // map to a distant position
        .collect();
    let mother_reads_ref: Vec<(i64, &[u8])> =
        mother_mismapped.iter().map(|(p, s)| (*p, s.as_slice())).collect();
    let mother_bam = create_bam(dir.path(), "mother", &ref_path, &mother_reads_ref);

    // Father: all REF at the correct position
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
        min_child_alt: 3,
        max_parent_alt: 1,
        min_child_alt_ratio: 0.1,
        max_variant_size: 50,
        window: 500,
        debug_kmers: true,
    };

    let summary = kmer_denovo::filter::run_filter(
        &child_bam,
        &mother_bam,
        &father_bam,
        &ref_path,
        &vcf_path,
        &output_path,
        &config,
        1,
        None,
    )
    .unwrap();

    assert_eq!(summary.total_variants, 1);
    assert_eq!(
        summary.dn_kmer_parent_alt, 1,
        "Expected DN_KMER_PARENT_ALT: whole-file scan should find ALT k-mers \
         in mother even though her reads are mapped to a different position"
    );
}
