use std::process::{Command, Output};

fn run_modkit(args: &[String]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_modkit"))
        .args(args)
        .output()
        .expect("failed to run modkit")
}

fn assert_sampling_fraction_parse_error(args: Vec<String>) {
    let output = run_modkit(&args);
    let stderr = String::from_utf8(output.stderr).unwrap();

    assert_eq!(
        output.status.code(),
        Some(2),
        "expected argument parsing to fail for {args:?}, stderr: {stderr}"
    );
    assert!(
        stderr.contains(
            "sampling fraction must be a finite number in the inclusive \
             range [0, 1]"
        ),
        "unexpected stderr for {args:?}: {stderr}"
    );
}

#[test]
fn invalid_sampling_fraction_fails_during_argument_parsing() {
    let temp_dir = tempfile::tempdir().unwrap();
    let missing_bam = temp_dir.path().join("input-that-does-not-exist.bam");
    let missing_reference =
        temp_dir.path().join("reference-that-does-not-exist.fa");
    let out_bed = temp_dir.path().join("pileup.bed");
    let out_calls = temp_dir.path().join("calls.tsv");
    let out_bam = temp_dir.path().join("calls.bam");
    let missing_bam = missing_bam.to_str().unwrap().to_string();
    let missing_reference = missing_reference.to_str().unwrap().to_string();

    let cases = [
        vec![
            "pileup".to_string(),
            missing_bam.clone(),
            out_bed.to_str().unwrap().to_string(),
            "--sampling-frac".to_string(),
            "1.1".to_string(),
        ],
        vec![
            "pileup-hemi".to_string(),
            missing_bam.clone(),
            "--reference".to_string(),
            missing_reference,
            "--cpg".to_string(),
            "--sampling-frac".to_string(),
            "1.1".to_string(),
        ],
        vec![
            "extract".to_string(),
            "calls".to_string(),
            missing_bam.clone(),
            out_calls.to_str().unwrap().to_string(),
            "--sampling-frac".to_string(),
            "1.1".to_string(),
        ],
        vec![
            "sample-probs".to_string(),
            missing_bam.clone(),
            "--sample-frac".to_string(),
            "1.1".to_string(),
        ],
        vec![
            "modbam".to_string(),
            "sample-probs".to_string(),
            missing_bam.clone(),
            "--sample-frac".to_string(),
            "1.1".to_string(),
        ],
        vec![
            "summary".to_string(),
            missing_bam.clone(),
            "--sampling-frac".to_string(),
            "1.1".to_string(),
        ],
        vec![
            "modbam".to_string(),
            "summary".to_string(),
            missing_bam.clone(),
            "--sampling-frac".to_string(),
            "1.1".to_string(),
        ],
        vec![
            "call-mods".to_string(),
            missing_bam.clone(),
            out_bam.to_str().unwrap().to_string(),
            "--sampling-frac".to_string(),
            "1.1".to_string(),
        ],
        vec![
            "modbam".to_string(),
            "call-mods".to_string(),
            missing_bam,
            out_bam.to_str().unwrap().to_string(),
            "--sampling-frac".to_string(),
            "1.1".to_string(),
        ],
    ];

    for args in cases {
        assert_sampling_fraction_parse_error(args);
    }

    for output_path in [&out_bed, &out_calls, &out_bam] {
        assert!(
            !output_path.exists(),
            "command created {output_path:?} before argument validation"
        );
    }
}

#[test]
fn boundary_sampling_fractions_reach_runtime_validation() {
    let temp_dir = tempfile::tempdir().unwrap();
    let missing_bam = temp_dir.path().join("missing.bam");
    let missing_bam = missing_bam.to_str().unwrap().to_string();

    for fraction in ["0", "1"] {
        let args = vec![
            "summary".to_string(),
            missing_bam.clone(),
            "--sampling-frac".to_string(),
            fraction.to_string(),
        ];
        let output = run_modkit(&args);
        let stderr = String::from_utf8(output.stderr).unwrap();

        assert_eq!(
            output.status.code(),
            Some(1),
            "expected {fraction} to parse and reach missing-input handling, \
             stderr: {stderr}"
        );
        assert!(
            !stderr.contains("sampling fraction must"),
            "boundary value {fraction} failed parsing: {stderr}"
        );
    }
}

#[test]
fn sampling_fraction_help_states_allowed_range() {
    let cases = [
        vec!["pileup", "--help"],
        vec!["pileup-hemi", "--help"],
        vec!["extract", "calls", "--help"],
        vec!["sample-probs", "--help"],
        vec!["summary", "--help"],
        vec!["call-mods", "--help"],
    ];

    for args in cases {
        let args = args.into_iter().map(String::from).collect::<Vec<_>>();
        let output = run_modkit(&args);
        let stdout = String::from_utf8(output.stdout).unwrap();
        let normalized =
            stdout.split_whitespace().collect::<Vec<_>>().join(" ");

        assert!(output.status.success(), "help failed for {args:?}");
        assert!(
            normalized.contains(
                "Must be a finite value in the inclusive range [0, 1]"
            ),
            "sampling range missing from help for {args:?}: {stdout}"
        );
    }
}
