use rust_htslib::bam::{
    self,
    header::HeaderRecord,
    record::{Aux, Cigar, CigarString},
};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

fn write_reference(root: &Path) -> PathBuf {
    let reference = root.join("reference.fa");
    fs::write(&reference, ">chr1\nAACGAA\n").unwrap();
    fs::write(root.join("reference.fa.fai"), "chr1\t6\t6\t6\t7\n").unwrap();
    reference
}

fn write_bam(root: &Path) -> PathBuf {
    let bam_path = root.join("reads.bam");
    let mut header = bam::Header::new();
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", "chr1").push_tag(b"LN", 6);
    header.push_record(&sq);

    let cigar = CigarString(vec![Cigar::Match(6)]);
    let mut record = bam::Record::new();
    record.set(b"read-0", Some(&cigar), b"AACGAA", &[30; 6]);
    record.set_tid(0);
    record.set_pos(0);
    record.set_mapq(60);
    record.push_aux(b"MM", Aux::String("C+m?,0;")).unwrap();
    record.push_aux(b"ML", Aux::ArrayU8((&[255][..]).into())).unwrap();
    record.push_aux(b"MN", Aux::U32(6)).unwrap();
    record.push_aux(b"NM", Aux::U32(0)).unwrap();

    let mut writer =
        bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
    writer.write(&record).unwrap();
    drop(writer);
    bam::index::build(&bam_path, None, bam::index::Type::Bai, 1).unwrap();
    bam_path
}

#[test]
fn singleton_region_emits_exact_summary_statistics() {
    let temp_dir = tempfile::tempdir().unwrap();
    let reference = write_reference(temp_dir.path());
    let bam = write_bam(temp_dir.path());
    let regions = temp_dir.path().join("regions.bed");
    fs::write(&regions, "chr1\t2\t3\tsingleton\n").unwrap();
    let output_dir = temp_dir.path().join("output");

    let result = Command::new(env!("CARGO_BIN_EXE_modkit"))
        .args([
            "entropy",
            "--in-bam",
            bam.to_str().unwrap(),
            "--out-bed",
            output_dir.to_str().unwrap(),
            "--ref",
            reference.to_str().unwrap(),
            "--base",
            "C",
            "--num-positions",
            "1",
            "--window-size",
            "1",
            "--min-coverage",
            "1",
            "--max-filtered-positions",
            "0",
            "--filter-threshold",
            "0",
            "--threads",
            "1",
            "--io-threads",
            "1",
            "--regions",
            regions.to_str().unwrap(),
            "--prefix",
            "singleton",
            "--suppress-progress",
        ])
        .output()
        .unwrap();

    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let windows = fs::read_to_string(
        output_dir.join("singleton_windows.bedgraph"),
    )
    .unwrap();
    let window_fields = windows.trim_end().split('\t').collect::<Vec<_>>();
    // Motif/window interval normalization belongs to issue #681. This test
    // isolates the regional summary and only requires one successful window
    // with the expected biological value and coverage.
    assert_eq!(window_fields[0], "chr1");
    assert_eq!(window_fields[1], "2");
    assert_eq!(&window_fields[3..], &["0", "+", "1"]);
    assert_eq!(
        fs::read(output_dir.join("singleton_regions.bed")).unwrap(),
        b"chr1\t2\t3\tsingleton\t0\t+\t0\t0\t0\t1\t1\t1\t1\t0\n"
    );
}
