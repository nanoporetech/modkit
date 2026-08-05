use rust_htslib::bam::{
    self,
    header::HeaderRecord,
    record::{Aux, Cigar, CigarString},
};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

fn write_reference(temp_dir: &Path) -> PathBuf {
    let reference = temp_dir.join("reference.fa");
    fs::write(&reference, ">chr1\nC\n").unwrap();
    let fai = PathBuf::from(format!("{}.fai", reference.display()));
    fs::write(fai, "chr1\t1\t6\t1\t2\n").unwrap();
    reference
}

fn write_ten_code_bam(temp_dir: &Path, name: &str, codes: &[char]) -> PathBuf {
    let bam_path = temp_dir.join(format!("{name}.bam"));
    let mut header = bam::Header::new();
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", "chr1").push_tag(b"LN", 1);
    header.push_record(&sq);

    let cigar = CigarString(vec![Cigar::Match(1)]);
    let mut writer =
        bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
    for (idx, code) in codes.iter().enumerate() {
        let mut record = bam::Record::new();
        record.set(
            format!("read-{idx:02}").as_bytes(),
            Some(&cigar),
            b"C",
            &[30],
        );
        record.set_tid(0);
        record.set_pos(0);
        record.set_mapq(60);
        record.set_flags(0);
        let mm = format!("C+{code}?,0;");
        record.push_aux(b"MM", Aux::String(&mm)).unwrap();
        record.push_aux(b"ML", Aux::ArrayU8((&[255][..]).into())).unwrap();
        record.push_aux(b"MN", Aux::U32(1)).unwrap();
        writer.write(&record).unwrap();
    }
    drop(writer);
    bam::index::build(&bam_path, None, bam::index::Type::Bai, 1).unwrap();
    bam_path
}

fn run_entropy(input: &Path, reference: &Path, output: &Path, threads: usize) {
    let result = Command::new(env!("CARGO_BIN_EXE_modkit"))
        .args([
            "entropy",
            "--in-bam",
            input.to_str().unwrap(),
            "--out-bed",
            output.to_str().unwrap(),
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
            "--no-filtering",
            "--threads",
            &threads.to_string(),
            "--io-threads",
            "1",
            "--suppress-progress",
        ])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "entropy failed:\n{}",
        String::from_utf8_lossy(&result.stderr)
    );
}

#[test]
fn ten_code_cli_is_stable_across_encounter_order_and_threads() {
    let temp_dir = tempfile::tempdir().unwrap();
    let reference = write_reference(temp_dir.path());
    let codes = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j'];
    let forward = write_ten_code_bam(temp_dir.path(), "forward", &codes);
    let mut reversed_codes = codes;
    reversed_codes.reverse();
    let reversed =
        write_ten_code_bam(temp_dir.path(), "reversed", &reversed_codes);

    let mut expected = None;
    for (order, input) in [("forward", &forward), ("reversed", &reversed)] {
        for threads in [1, 4] {
            let output = temp_dir.path().join(format!("{order}-{threads}.bed"));
            run_entropy(input, &reference, &output, threads);
            let observed = fs::read(&output).unwrap();
            if let Some(expected) = expected.as_ref() {
                assert_eq!(&observed, expected);
            } else {
                expected = Some(observed);
            }
        }
    }

    let output = String::from_utf8(expected.unwrap()).unwrap();
    let rows = output.lines().collect::<Vec<_>>();
    assert_eq!(rows.len(), 1);
    let fields = rows[0].split('\t').collect::<Vec<_>>();
    assert_eq!(fields.len(), 6);
    assert_eq!(fields[5], "10");
    let entropy = fields[3].parse::<f32>().unwrap();
    assert!((entropy - 10f32.log2()).abs() < 0.000_01);
}
