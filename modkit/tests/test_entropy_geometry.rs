use rust_htslib::bam::{
    self,
    header::HeaderRecord,
    record::{Aux, Cigar, CigarString},
};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn write_reference(root: &Path, sequence: &str) -> PathBuf {
    let reference = root.join("reference.fa");
    fs::write(&reference, format!(">chr1\n{sequence}\n")).unwrap();
    fs::write(
        root.join("reference.fa.fai"),
        format!(
            "chr1\t{}\t6\t{}\t{}\n",
            sequence.len(),
            sequence.len(),
            sequence.len() + 1
        ),
    )
    .unwrap();
    reference
}

fn write_bam(
    root: &Path,
    sequence: &str,
    mm_tag: &str,
    ml_count: usize,
    strands: &[bool],
) -> PathBuf {
    let bam_path = root.join("reads.bam");
    let mut header = bam::Header::new();
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", "chr1").push_tag(b"LN", sequence.len());
    header.push_record(&sq);

    let cigar = CigarString(vec![Cigar::Match(sequence.len() as u32)]);
    let mut writer =
        bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
    for (idx, reverse) in strands.iter().copied().enumerate() {
        let mut record = bam::Record::new();
        record.set(
            format!("read-{idx}").as_bytes(),
            Some(&cigar),
            sequence.as_bytes(),
            &vec![30; sequence.len()],
        );
        record.set_tid(0);
        record.set_pos(0);
        record.set_mapq(60);
        if reverse {
            record.set_flags(16);
        }
        let ml = vec![255u8; ml_count];
        record.push_aux(b"MM", Aux::String(mm_tag)).unwrap();
        record.push_aux(b"ML", Aux::ArrayU8((&ml[..]).into())).unwrap();
        record.push_aux(b"MN", Aux::U32(sequence.len() as u32)).unwrap();
        record.push_aux(b"NM", Aux::U32(0)).unwrap();
        writer.write(&record).unwrap();
    }
    drop(writer);
    bam::index::build(&bam_path, None, bam::index::Type::Bai, 1).unwrap();
    bam_path
}

fn write_two_contig_conflict_fixture(root: &Path) -> (PathBuf, PathBuf) {
    let reference = root.join("two-contigs.fa");
    fs::write(&reference, ">chr1\nCG\n>chr2\nCGCG\n").unwrap();
    fs::write(
        root.join("two-contigs.fa.fai"),
        "chr1\t2\t6\t2\t3\nchr2\t4\t15\t4\t5\n",
    )
    .unwrap();

    let bam_path = root.join("two-contigs.bam");
    let mut header = bam::Header::new();
    for (name, length) in [("chr1", 2), ("chr2", 4)] {
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", name).push_tag(b"LN", length);
        header.push_record(&sq);
    }
    let mut writer =
        bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
    for (tid, sequence, mm_tag, ml_count) in [
        (0, "CG", "C+m?,0;", 1),
        (1, "CGCG", "C+m?,0,0;", 2),
    ] {
        let cigar = CigarString(vec![Cigar::Match(sequence.len() as u32)]);
        let mut record = bam::Record::new();
        record.set(
            format!("read-{tid}").as_bytes(),
            Some(&cigar),
            sequence.as_bytes(),
            &vec![30; sequence.len()],
        );
        record.set_tid(tid);
        record.set_pos(0);
        record.set_mapq(60);
        let ml = vec![255u8; ml_count];
        record.push_aux(b"MM", Aux::String(mm_tag)).unwrap();
        record.push_aux(b"ML", Aux::ArrayU8((&ml[..]).into())).unwrap();
        record.push_aux(b"MN", Aux::U32(sequence.len() as u32)).unwrap();
        record.push_aux(b"NM", Aux::U32(0)).unwrap();
        writer.write(&record).unwrap();
    }
    drop(writer);
    bam::index::build(&bam_path, None, bam::index::Type::Bai, 1).unwrap();
    (reference, bam_path)
}

fn run_entropy(
    bam: &Path,
    reference: &Path,
    output: &Path,
    window_size: usize,
    regions: Option<&Path>,
) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_modkit"));
    command.args([
        "entropy",
        "--in-bam",
        bam.to_str().unwrap(),
        "--out-bed",
        output.to_str().unwrap(),
        "--ref",
        reference.to_str().unwrap(),
        "--motif",
        "CG",
        "0",
        "--num-positions",
        "1",
        "--window-size",
        &window_size.to_string(),
        "--min-coverage",
        "1",
        "--max-filtered-positions",
        "0",
        "--no-filtering",
        "--threads",
        "1",
        "--io-threads",
        "1",
        "--suppress-progress",
        "--force",
    ]);
    if let Some(regions) = regions {
        command.args([
            "--regions",
            regions.to_str().unwrap(),
            "--prefix",
            "anchor",
        ]);
    }
    command.output().unwrap()
}

#[test]
fn cg_rows_are_invariant_to_search_window_size() {
    let temp_dir = tempfile::tempdir().unwrap();
    let reference = write_reference(temp_dir.path(), "CGTACG");
    let bam =
        write_bam(temp_dir.path(), "CGTACG", "C+m?,0,0;", 2, &[false, true]);
    let expected = concat!(
        "chr1\t0\t1\t0\t+\t1\n",
        "chr1\t1\t2\t0\t-\t1\n",
        "chr1\t4\t5\t0\t+\t1\n",
        "chr1\t5\t6\t0\t-\t1\n",
    );

    for window_size in [1, 2, 3, 6, 100] {
        let output = temp_dir.path().join(format!("entropy-{window_size}.bed"));
        let result = run_entropy(&bam, &reference, &output, window_size, None);
        assert!(
            result.status.success(),
            "window_size={window_size}:\n{}",
            String::from_utf8_lossy(&result.stderr)
        );
        assert_eq!(
            fs::read(output).unwrap(),
            expected.as_bytes(),
            "window_size={window_size}"
        );
    }
}

#[test]
fn region_owning_only_the_anchor_uses_reference_context() {
    let temp_dir = tempfile::tempdir().unwrap();
    let reference = write_reference(temp_dir.path(), "AACGAA");
    let bam = write_bam(temp_dir.path(), "AACGAA", "C+m?,0;", 1, &[false]);
    let regions = temp_dir.path().join("regions.bed");
    fs::write(&regions, "chr1\t2\t3\tanchor-only\n").unwrap();
    let output_dir = temp_dir.path().join("entropy-regions");

    let result = run_entropy(&bam, &reference, &output_dir, 1, Some(&regions));
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(
        fs::read(output_dir.join("anchor_windows.bedgraph")).unwrap(),
        b"chr1\t2\t3\t0\t+\t1\n"
    );
    let region_rows =
        fs::read_to_string(output_dir.join("anchor_regions.bed")).unwrap();
    let fields = region_rows.trim_end().split('\t').collect::<Vec<_>>();
    assert_eq!(&fields[..4], &["chr1", "2", "3", "anchor-only"]);
    assert_eq!(fields[5], "+");
    assert_eq!(&fields[12..], &["1", "0"]);
}

#[test]
fn conflicting_combined_motif_partners_fail_before_output_creation() {
    let temp_dir = tempfile::tempdir().unwrap();
    let reference = write_reference(temp_dir.path(), "CGCG");
    let bam = write_bam(
        temp_dir.path(),
        "CGCG",
        "C+m?,0,0;",
        2,
        &[false],
    );
    let output = temp_dir.path().join("conflict.bed");
    let result = Command::new(env!("CARGO_BIN_EXE_modkit"))
        .args([
            "entropy",
            "--in-bam",
            bam.to_str().unwrap(),
            "--out-bed",
            output.to_str().unwrap(),
            "--ref",
            reference.to_str().unwrap(),
            "--motif",
            "CG",
            "0",
            "--motif",
            "CGCG",
            "0",
            "--combine-strands",
            "--num-positions",
            "1",
            "--window-size",
            "4",
            "--min-coverage",
            "1",
            "--no-filtering",
            "--threads",
            "1",
            "--io-threads",
            "1",
            "--suppress-progress",
        ])
        .output()
        .unwrap();

    assert!(!result.status.success());
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(
        stderr.contains("conflicting combined-strand motif partners"),
        "{stderr}"
    );
    assert!(stderr.contains("chr1:0"), "{stderr}");
    assert!(!output.exists());
}

#[test]
fn conflict_on_later_contig_fails_before_output_creation() {
    let temp_dir = tempfile::tempdir().unwrap();
    let (reference, bam) =
        write_two_contig_conflict_fixture(temp_dir.path());
    let output = temp_dir.path().join("later-conflict.bed");
    let result = Command::new(env!("CARGO_BIN_EXE_modkit"))
        .args([
            "entropy",
            "--in-bam",
            bam.to_str().unwrap(),
            "--out-bed",
            output.to_str().unwrap(),
            "--ref",
            reference.to_str().unwrap(),
            "--motif",
            "CG",
            "0",
            "--motif",
            "CGCG",
            "0",
            "--combine-strands",
            "--num-positions",
            "1",
            "--window-size",
            "4",
            "--min-coverage",
            "1",
            "--no-filtering",
            "--threads",
            "1",
            "--io-threads",
            "1",
            "--suppress-progress",
        ])
        .output()
        .unwrap();

    assert!(!result.status.success());
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(
        stderr.contains("conflicting combined-strand motif partners"),
        "{stderr}"
    );
    assert!(stderr.contains("chr2:0"), "{stderr}");
    assert!(!output.exists());
}
