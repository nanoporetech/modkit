use rust_htslib::bam::{self, Read};
use std::collections::HashSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const UNMAPPED_BAM: &str =
    "../tests/resources/bc_anchored_10_reads.unmapped.bam";
const UNMAPPED_CRAM: &str =
    "../tests/resources/bc_anchored_10_reads_unmapped.cram";
const UNMAPPED_CRAM_INDEX: &str =
    "../tests/resources/bc_anchored_10_reads_unmapped.cram.crai";
const EMPTY_THRESHOLD_ERROR: &str =
    "cannot calculate automatic thresholds because no modification \
     probabilities were sampled";

struct InputCopies {
    serial_bam: PathBuf,
    indexed_bam: PathBuf,
    serial_cram: PathBuf,
    indexed_cram: PathBuf,
    protected_bytes: Vec<(PathBuf, Vec<u8>)>,
}

impl InputCopies {
    fn new(temp_dir: &Path) -> Self {
        let serial_bam = temp_dir.join("serial.bam");
        let indexed_bam = temp_dir.join("indexed.bam");
        let serial_cram = temp_dir.join("serial.cram");
        let indexed_cram = temp_dir.join("indexed.cram");
        let indexed_cram_index = temp_dir.join("indexed.cram.crai");

        fs::copy(UNMAPPED_BAM, &serial_bam).unwrap();
        fs::copy(UNMAPPED_BAM, &indexed_bam).unwrap();
        bam::index::build(&indexed_bam, None, bam::index::Type::Bai, 1)
            .unwrap();
        let indexed_bam_index = temp_dir.join("indexed.bam.bai");
        assert!(indexed_bam_index.is_file());

        fs::copy(UNMAPPED_CRAM, &serial_cram).unwrap();
        fs::copy(UNMAPPED_CRAM, &indexed_cram).unwrap();
        fs::copy(UNMAPPED_CRAM_INDEX, &indexed_cram_index).unwrap();

        assert!(!temp_dir.join("serial.bam.bai").exists());
        assert!(!temp_dir.join("serial.cram.crai").exists());

        let protected_bytes = [
            &serial_bam,
            &indexed_bam,
            &indexed_bam_index,
            &serial_cram,
            &indexed_cram,
            &indexed_cram_index,
        ]
        .into_iter()
        .map(|path| (path.clone(), fs::read(path).unwrap()))
        .collect();

        Self {
            serial_bam,
            indexed_bam,
            serial_cram,
            indexed_cram,
            protected_bytes,
        }
    }

    fn assert_unchanged(&self) {
        for (path, expected) in &self.protected_bytes {
            assert_eq!(
                fs::read(path).unwrap(),
                *expected,
                "input or index changed: {}",
                path.display()
            );
        }
    }
}

fn run_extract(
    command: &str,
    input: &Path,
    output: &Path,
    serial: bool,
    mode_args: &[&str],
) -> Output {
    let mut args = vec![
        "extract".to_string(),
        command.to_string(),
        input.to_str().unwrap().to_string(),
        output.to_str().unwrap().to_string(),
        "--threads".to_string(),
        "1".to_string(),
        "--io-threads".to_string(),
        "1".to_string(),
        "--no-headers".to_string(),
        "--suppress-progress".to_string(),
        "--force".to_string(),
    ];
    if command == "calls" {
        args.extend(["--sample-num-reads".to_string(), "10".to_string()]);
    }
    if serial {
        args.push("--ignore-index".to_string());
    }
    args.extend(mode_args.iter().map(|arg| arg.to_string()));
    Command::new(env!("CARGO_BIN_EXE_modkit")).args(args).output().unwrap()
}

fn assert_success(output: &Output) {
    assert!(
        output.status.success(),
        "command failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn assert_exact_output(output: &Output, path: &Path, expected_rows: usize) {
    assert_success(output);
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("processed 10 reads"),
        "unmapped reads were not processed exactly once:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let bytes = fs::read(path).unwrap();
    assert_eq!(
        std::str::from_utf8(&bytes).unwrap().lines().count(),
        expected_rows
    );
    let read_ids = std::str::from_utf8(&bytes)
        .unwrap()
        .lines()
        .map(|line| line.split('\t').next().unwrap())
        .collect::<HashSet<_>>();
    assert_eq!(read_ids.len(), 10);
}

#[test]
fn indexed_zero_target_bam_and_cram_match_serial_extract_modes() {
    let temp_dir = tempfile::tempdir().unwrap();
    let inputs = InputCopies::new(temp_dir.path());

    for (command, mode, mode_args, expected_rows) in [
        ("full", "full", &[][..], 218),
        ("calls", "automatic", &[][..], 109),
        ("calls", "explicit", &["--filter-threshold", "0"][..], 109),
        ("calls", "no-filter", &["--no-filtering"][..], 109),
    ] {
        let mut expected_bytes = None;
        for (format, serial_input, indexed_input) in [
            ("bam", &inputs.serial_bam, &inputs.indexed_bam),
            ("cram", &inputs.serial_cram, &inputs.indexed_cram),
        ] {
            for (route, input, serial) in [
                ("serial", serial_input, true),
                ("indexed", indexed_input, false),
            ] {
                let output_path = temp_dir
                    .path()
                    .join(format!("{format}-{mode}-{route}.tsv"));
                let output = run_extract(
                    command,
                    input,
                    &output_path,
                    serial,
                    mode_args,
                );
                assert_exact_output(&output, &output_path, expected_rows);
                let observed = fs::read(&output_path).unwrap();
                if let Some(expected) = expected_bytes.as_ref() {
                    assert_eq!(
                        &observed, expected,
                        "{format} {mode} {route} output differed"
                    );
                } else {
                    expected_bytes = Some(observed);
                }
            }
        }
    }

    inputs.assert_unchanged();
}

fn write_without_modification_tags(input: &Path, output: &Path) {
    let mut reader = bam::Reader::from_path(input).unwrap();
    let header = bam::Header::from_template(reader.header());
    let mut writer =
        bam::Writer::from_path(output, &header, bam::Format::Bam).unwrap();
    for record in reader.records() {
        let mut record = record.unwrap();
        for tag in [b"MM", b"ML", b"Mm", b"Ml"] {
            if record.aux(tag).is_ok() {
                record.remove_aux(tag).unwrap();
            }
        }
        writer.write(&record).unwrap();
    }
}

#[test]
fn indexed_zero_target_automatic_threshold_preserves_empty_sample_contract() {
    let temp_dir = tempfile::tempdir().unwrap();
    let serial_bam = temp_dir.path().join("no-observations-serial.bam");
    let indexed_bam = temp_dir.path().join("no-observations-indexed.bam");
    write_without_modification_tags(Path::new(UNMAPPED_BAM), &serial_bam);
    fs::copy(&serial_bam, &indexed_bam).unwrap();
    bam::index::build(&indexed_bam, None, bam::index::Type::Bai, 1).unwrap();
    let indexed_bam_index =
        temp_dir.path().join("no-observations-indexed.bam.bai");
    let protected_bytes = [
        (&serial_bam, fs::read(&serial_bam).unwrap()),
        (&indexed_bam, fs::read(&indexed_bam).unwrap()),
        (&indexed_bam_index, fs::read(&indexed_bam_index).unwrap()),
    ];

    for (route, input, serial) in [
        ("serial", serial_bam.as_path(), true),
        ("indexed", indexed_bam.as_path(), false),
    ] {
        let output_path =
            temp_dir.path().join(format!("automatic-{route}.tsv"));
        let output = run_extract("calls", input, &output_path, serial, &[]);
        assert!(!output.status.success());
        assert!(
            String::from_utf8_lossy(&output.stderr)
                .contains(EMPTY_THRESHOLD_ERROR),
            "unexpected {route} automatic-threshold error:\n{}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(!output_path.exists());
    }

    for mode_args in [&["--filter-threshold", "0"][..], &["--no-filtering"][..]]
    {
        let serial_output = temp_dir.path().join(format!(
            "control-serial-{}.tsv",
            mode_args[0].trim_start_matches('-')
        ));
        let indexed_output = temp_dir.path().join(format!(
            "control-indexed-{}.tsv",
            mode_args[0].trim_start_matches('-')
        ));
        assert_success(&run_extract(
            "calls",
            &serial_bam,
            &serial_output,
            true,
            mode_args,
        ));
        assert_success(&run_extract(
            "calls",
            &indexed_bam,
            &indexed_output,
            false,
            mode_args,
        ));
        assert_eq!(fs::read(&serial_output).unwrap(), Vec::<u8>::new());
        assert_eq!(fs::read(&indexed_output).unwrap(), Vec::<u8>::new());
    }

    for (path, expected) in protected_bytes {
        assert_eq!(fs::read(path).unwrap(), expected);
    }
}
