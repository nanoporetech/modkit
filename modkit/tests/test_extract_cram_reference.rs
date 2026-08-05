use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const BAM: &str = "../tests/resources/bc_anchored_10_reads.sorted.bam";
const CRAM: &str = "../tests/resources/bc_anchored_10_reads.sorted.cram";
const REFERENCE: &str = "../tests/resources/CGI_ladder_3.6kb_ref.fa";
const WRONG_REFERENCE: &str = "../tests/resources/genome_for_phased_test.fasta";
const UNMAPPED_BAM: &str =
    "../tests/resources/bc_anchored_10_reads.unmapped.bam";
const UNMAPPED_CRAM: &str =
    "../tests/resources/bc_anchored_10_reads_unmapped.cram";

fn run_modkit(args: &[String], empty_reference_cache: Option<&Path>) -> Output {
    let exe = Path::new(env!("CARGO_BIN_EXE_modkit"));
    assert!(exe.exists());
    let mut command = Command::new(exe);
    command.args(args);
    if let Some(empty_reference_cache) = empty_reference_cache {
        let reference_pattern = empty_reference_cache.join("%2s/%2s/%s");
        command
            .env("REF_PATH", &reference_pattern)
            .env("REF_CACHE", &reference_pattern);
    }
    command.output().unwrap()
}

fn run_extract(
    command: &str,
    input: &Path,
    output: &Path,
    reference: Option<&Path>,
    extra_args: &[&str],
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
    if let Some(reference) = reference {
        args.extend([
            "--reference".to_string(),
            reference.to_str().unwrap().to_string(),
        ]);
    }
    args.extend(extra_args.iter().map(|arg| arg.to_string()));
    run_modkit(&args, None)
}

fn assert_success(output: &Output) {
    assert!(
        output.status.success(),
        "command failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn assert_parity(expected: &Path, observed: &Path, expected_rows: usize) {
    let expected = fs::read(expected).unwrap();
    let observed = fs::read(observed).unwrap();
    assert_eq!(
        row_count(&expected),
        expected_rows,
        "unexpected control row count"
    );
    assert_eq!(
        row_count(&observed),
        expected_rows,
        "unexpected CRAM row count"
    );
    assert_eq!(expected, observed);
}

fn row_count(output: &[u8]) -> usize {
    std::str::from_utf8(output).unwrap().lines().count()
}

fn unindexed_copies(temp_dir: &Path) -> (PathBuf, PathBuf) {
    let bam = temp_dir.join("unindexed.bam");
    let cram = temp_dir.join("unindexed.cram");
    fs::copy(BAM, &bam).unwrap();
    fs::copy(CRAM, &cram).unwrap();
    assert!(!bam.with_extension("bam.bai").exists());
    assert!(!cram.with_extension("cram.crai").exists());
    (bam, cram)
}

#[test]
fn mapped_cram_matches_bam_for_full_and_calls() {
    let temp_dir = tempfile::tempdir().unwrap();
    let (unindexed_bam, unindexed_cram) = unindexed_copies(temp_dir.path());
    let reference = Path::new(REFERENCE);

    let full_bam_indexed = temp_dir.path().join("full.bam.indexed.tsv");
    let full_cram_indexed = temp_dir.path().join("full.cram.indexed.tsv");
    let full_bam_serial = temp_dir.path().join("full.bam.serial.tsv");
    let full_cram_serial = temp_dir.path().join("full.cram.serial.tsv");
    for (input, output, extra_args) in [
        (Path::new(BAM), full_bam_indexed.as_path(), &[][..]),
        (Path::new(CRAM), full_cram_indexed.as_path(), &[][..]),
        (
            unindexed_bam.as_path(),
            full_bam_serial.as_path(),
            &["--ignore-index"][..],
        ),
        (
            unindexed_cram.as_path(),
            full_cram_serial.as_path(),
            &["--ignore-index"][..],
        ),
    ] {
        assert_success(&run_extract(
            "full",
            input,
            output,
            Some(reference),
            extra_args,
        ));
    }
    assert_parity(&full_bam_indexed, &full_cram_indexed, 218);
    assert_parity(&full_bam_serial, &full_cram_serial, 218);
    assert_parity(&full_bam_indexed, &full_bam_serial, 218);

    for (mode, mode_args) in [
        ("automatic", &[][..]),
        ("explicit", &["--filter-threshold", "0"][..]),
        ("no-filter", &["--no-filtering"][..]),
    ] {
        let bam_indexed =
            temp_dir.path().join(format!("calls.{mode}.bam.indexed.tsv"));
        let cram_indexed =
            temp_dir.path().join(format!("calls.{mode}.cram.indexed.tsv"));
        let bam_serial =
            temp_dir.path().join(format!("calls.{mode}.bam.serial.tsv"));
        let cram_serial =
            temp_dir.path().join(format!("calls.{mode}.cram.serial.tsv"));

        for (input, output, extra_args) in [
            (Path::new(BAM), bam_indexed.as_path(), mode_args),
            (Path::new(CRAM), cram_indexed.as_path(), mode_args),
            (unindexed_bam.as_path(), bam_serial.as_path(), mode_args),
            (unindexed_cram.as_path(), cram_serial.as_path(), mode_args),
        ] {
            assert_success(&run_extract(
                "calls",
                input,
                output,
                Some(reference),
                extra_args,
            ));
        }
        assert_parity(&bam_indexed, &cram_indexed, 109);
        assert_parity(&bam_serial, &cram_serial, 109);
        assert_parity(&bam_indexed, &bam_serial, 109);
    }
}

#[test]
fn indexed_cram_num_reads_schedule_uses_reference() {
    let temp_dir = tempfile::tempdir().unwrap();
    for (command, mode_args, expected_rows) in [
        ("full", &[][..], 218),
        ("calls", &["--filter-threshold", "0"][..], 109),
    ] {
        let output_path = temp_dir.path().join(format!("{command}.tsv"));
        let mut args = mode_args.to_vec();
        args.extend(["--num-reads", "10042"]);
        assert_success(&run_extract(
            command,
            Path::new(CRAM),
            &output_path,
            Some(Path::new(REFERENCE)),
            &args,
        ));
        let output = fs::read(&output_path).unwrap();
        assert_eq!(row_count(&output), expected_rows);
    }
}

#[test]
fn unmapped_cram_without_reference_matches_bam() {
    let temp_dir = tempfile::tempdir().unwrap();
    for (command, mode_args, expected_rows) in
        [("full", &[][..], 218), ("calls", &["--no-filtering"][..], 109)]
    {
        let bam_output = temp_dir.path().join(format!("{command}.bam.tsv"));
        let cram_output = temp_dir.path().join(format!("{command}.cram.tsv"));
        let mut args = mode_args.to_vec();
        args.push("--ignore-index");
        assert_success(&run_extract(
            command,
            Path::new(UNMAPPED_BAM),
            &bam_output,
            None,
            &args,
        ));
        assert_success(&run_extract(
            command,
            Path::new(UNMAPPED_CRAM),
            &cram_output,
            None,
            &args,
        ));
        assert_parity(&bam_output, &cram_output, expected_rows);
    }
}

#[test]
fn reference_errors_do_not_mutate_extract_outputs() {
    let temp_dir = tempfile::tempdir().unwrap();
    let empty_reference_cache = temp_dir.path().join("empty-reference-cache");
    fs::create_dir(&empty_reference_cache).unwrap();
    let disguised_cram = temp_dir.path().join("reference-dependent.bam");
    fs::copy(CRAM, &disguised_cram).unwrap();
    let missing_reference = temp_dir.path().join("missing-reference.fa");

    for (reference_name, reference) in [
        ("none", None),
        ("missing", Some(missing_reference.as_path())),
        ("wrong", Some(Path::new(WRONG_REFERENCE))),
    ] {
        for command in ["full", "calls"] {
            for preexisting in [false, true] {
                let output_path = temp_dir.path().join(format!(
                    "{command}-{reference_name}-{}.tsv",
                    if preexisting { "existing" } else { "absent" }
                ));
                if preexisting {
                    fs::write(&output_path, b"sentinel\n").unwrap();
                }
                let mut args = vec![
                    "extract".to_string(),
                    command.to_string(),
                    disguised_cram.to_str().unwrap().to_string(),
                    output_path.to_str().unwrap().to_string(),
                    "--threads".to_string(),
                    "1".to_string(),
                    "--io-threads".to_string(),
                    "1".to_string(),
                    "--suppress-progress".to_string(),
                    "--force".to_string(),
                ];
                if command == "calls" {
                    args.push("--no-filtering".to_string());
                }
                if let Some(reference) = reference {
                    args.extend([
                        "--reference".to_string(),
                        reference.to_str().unwrap().to_string(),
                    ]);
                }
                let output =
                    run_modkit(&args, Some(empty_reference_cache.as_path()));
                assert!(
                    !output.status.success(),
                    "{command} unexpectedly accepted {reference_name} reference"
                );
                if preexisting {
                    assert_eq!(fs::read(&output_path).unwrap(), b"sentinel\n");
                } else {
                    assert!(!output_path.exists());
                }
            }
        }
    }
}
