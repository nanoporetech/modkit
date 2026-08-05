use std::fs;
use std::path::Path;
use std::process::{Command, Output};

const BAM: &str = "../tests/resources/bc_anchored_10_reads.sorted.bam";
const CRAM: &str = "../tests/resources/bc_anchored_10_reads.sorted.cram";
const REFERENCE: &str = "../tests/resources/CGI_ladder_3.6kb_ref.fa";
const REFERENCE_INDEPENDENT_CRAM: &str =
    "../tests/resources/bc_anchored_10_reads_unmapped.cram";

fn run_modkit(args: &[&str]) -> Output {
    let exe = Path::new(env!("CARGO_BIN_EXE_modkit"));
    assert!(exe.exists());
    Command::new(exe).args(args).output().unwrap()
}

fn run_modkit_without_reference_resolution(
    args: &[&str],
    empty_reference_cache: &Path,
) -> Output {
    let exe = Path::new(env!("CARGO_BIN_EXE_modkit"));
    assert!(exe.exists());
    let reference_pattern = empty_reference_cache.join("%2s/%2s/%s");
    Command::new(exe)
        .args(args)
        .env("REF_PATH", &reference_pattern)
        .env("REF_CACHE", &reference_pattern)
        .output()
        .unwrap()
}

fn assert_success(output: &Output) {
    assert!(
        output.status.success(),
        "command failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn normalized_sam_records(path: &Path) -> Vec<String> {
    fs::read_to_string(path)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| {
            let mut fields = line.split('\t').collect::<Vec<_>>();
            let mut aux = fields.split_off(11);
            aux.sort_unstable();
            fields.extend(aux);
            fields.join("\t")
        })
        .collect()
}

fn threshold_messages(output: &Output) -> Vec<&str> {
    std::str::from_utf8(&output.stderr)
        .unwrap()
        .lines()
        .filter(|line| {
            line.contains("Threshold of")
                || line.contains("calculated thresholds:")
        })
        .collect()
}

fn assert_entropy_parity(bam_output: &Path, cram_output: &Path) {
    let bam = fs::read_to_string(bam_output).unwrap();
    let cram = fs::read_to_string(cram_output).unwrap();
    let bam_lines = bam.lines().collect::<Vec<_>>();
    let cram_lines = cram.lines().collect::<Vec<_>>();
    assert_eq!(bam_lines.len(), cram_lines.len());
    for (bam_line, cram_line) in bam_lines.iter().zip(cram_lines) {
        let bam_fields = bam_line.split('\t').collect::<Vec<_>>();
        let cram_fields = cram_line.split('\t').collect::<Vec<_>>();
        assert_eq!(bam_fields.len(), 6);
        assert_eq!(cram_fields.len(), 6);
        assert_eq!(&bam_fields[..3], &cram_fields[..3]);
        assert_eq!(&bam_fields[4..], &cram_fields[4..]);
        let bam_entropy = bam_fields[3].parse::<f64>().unwrap();
        let cram_entropy = cram_fields[3].parse::<f64>().unwrap();
        assert!(
            (bam_entropy - cram_entropy).abs() <= 1e-6,
            "entropy differed: {bam_entropy} versus {cram_entropy}"
        );
    }
}

#[test]
fn call_and_adjust_mods_match_bam_and_cram() {
    let temp_dir = tempfile::tempdir().unwrap();
    for command in ["call-mods", "adjust-mods"] {
        let help = run_modkit(&[command, "--help"]);
        assert_success(&help);
        assert!(String::from_utf8_lossy(&help.stdout)
            .contains("-r, --reference <REFERENCE_FASTA>"));
    }

    let call_bam = temp_dir.path().join("call.bam.sam");
    let call_cram = temp_dir.path().join("call.cram.sam");

    let bam_output = run_modkit(&[
        "call-mods",
        BAM,
        call_bam.to_str().unwrap(),
        "--output-sam",
        "--sampling-frac",
        "1",
        "--seed",
        "7",
        "--threads",
        "1",
        "--sampling-interval-size",
        "20",
        "--reference",
        REFERENCE,
        "--suppress-progress",
    ]);
    let cram_output = run_modkit(&[
        "call-mods",
        CRAM,
        call_cram.to_str().unwrap(),
        "--output-sam",
        "--sampling-frac",
        "1",
        "--seed",
        "7",
        "--threads",
        "1",
        "--sampling-interval-size",
        "20",
        "--ref",
        REFERENCE,
        "--suppress-progress",
    ]);
    assert_success(&bam_output);
    assert_success(&cram_output);
    let bam_records = normalized_sam_records(&call_bam);
    assert_eq!(bam_records.len(), 10);
    assert_eq!(bam_records, normalized_sam_records(&call_cram));

    let adjust_bam = temp_dir.path().join("adjust.bam.sam");
    let adjust_cram = temp_dir.path().join("adjust.cram.sam");
    let bam_output = run_modkit(&[
        "adjust-mods",
        BAM,
        adjust_bam.to_str().unwrap(),
        "--output-sam",
        "--filter-probs",
        "--filter-threshold",
        "0.55",
        "--threads",
        "1",
        "-r",
        REFERENCE,
        "--suppress-progress",
    ]);
    let cram_output = run_modkit(&[
        "adjust-mods",
        CRAM,
        adjust_cram.to_str().unwrap(),
        "--output-sam",
        "--filter-probs",
        "--filter-threshold",
        "0.55",
        "--threads",
        "1",
        "--reference",
        REFERENCE,
        "--suppress-progress",
    ]);
    assert_success(&bam_output);
    assert_success(&cram_output);
    let bam_records = normalized_sam_records(&adjust_bam);
    assert_eq!(bam_records.len(), 10);
    assert_eq!(bam_records, normalized_sam_records(&adjust_cram));

    let sampled_adjust_cram = temp_dir.path().join("adjust.sampled.cram.sam");
    let output = run_modkit(&[
        "adjust-mods",
        CRAM,
        sampled_adjust_cram.to_str().unwrap(),
        "--output-sam",
        "--filter-probs",
        "--num-reads",
        "10042",
        "--threads",
        "1",
        "--sampling-interval-size",
        "20",
        "--ref",
        REFERENCE,
        "--suppress-progress",
    ]);
    assert_success(&output);
    assert_eq!(normalized_sam_records(&sampled_adjust_cram).len(), 10);
}

#[test]
fn pileup_hemi_and_entropy_match_bam_and_cram() {
    let temp_dir = tempfile::tempdir().unwrap();
    let hemi_bam = temp_dir.path().join("hemi.bam.bed");
    let hemi_cram = temp_dir.path().join("hemi.cram.bed");
    let include_bed = "../tests/resources/include-pos-1-site.bed";

    let run_hemi = |input: &str, output: &Path| {
        run_modkit(&[
            "pileup-hemi",
            input,
            "--out-bed",
            output.to_str().unwrap(),
            "--ref",
            REFERENCE,
            "--cpg",
            "--region",
            "oligo_1512_adapters",
            "--include-bed",
            include_bed,
            "--sampling-frac",
            "1",
            "--seed",
            "7",
            "--threads",
            "1",
            "--interval-size",
            "20",
            "--sampling-interval-size",
            "20",
            "--suppress-progress",
        ])
    };
    let bam_output = run_hemi(BAM, &hemi_bam);
    let cram_output = run_hemi(CRAM, &hemi_cram);
    assert_success(&bam_output);
    assert_success(&cram_output);
    let bam_thresholds = threshold_messages(&bam_output);
    let cram_thresholds = threshold_messages(&cram_output);
    assert!(!bam_thresholds.is_empty());
    assert_eq!(bam_thresholds, cram_thresholds);
    assert!(String::from_utf8_lossy(&bam_output.stderr)
        .contains("Processed ~10 reads"));
    assert!(String::from_utf8_lossy(&cram_output.stderr)
        .contains("Processed ~10 reads"));
    assert_eq!(fs::read(&hemi_bam).unwrap(), fs::read(&hemi_cram).unwrap());

    let entropy_bam = temp_dir.path().join("entropy.bam.bed");
    let entropy_cram = temp_dir.path().join("entropy.cram.bed");
    let run_entropy = |input: &str, output: &Path| {
        run_modkit(&[
            "entropy",
            "--in-bam",
            input,
            "--out-bed",
            output.to_str().unwrap(),
            "--ref",
            REFERENCE,
            "--cpg",
            "--min-coverage",
            "1",
            "--num-reads",
            "10042",
            "--threads",
            "1",
            "--io-threads",
            "1",
            "--suppress-progress",
            "--force",
        ])
    };
    let bam_output = run_entropy(BAM, &entropy_bam);
    let cram_output = run_entropy(CRAM, &entropy_cram);
    assert_success(&bam_output);
    assert_success(&cram_output);
    let bam_thresholds = threshold_messages(&bam_output);
    let cram_thresholds = threshold_messages(&cram_output);
    assert!(!bam_thresholds.is_empty());
    assert_eq!(bam_thresholds, cram_thresholds);
    assert!(String::from_utf8_lossy(&bam_output.stderr)
        .contains("6 windows processed successfully"));
    assert!(String::from_utf8_lossy(&cram_output.stderr)
        .contains("6 windows processed successfully"));
    assert_entropy_parity(&entropy_bam, &entropy_cram);
}

#[test]
fn missing_external_reference_fails_before_output_without_rejecting_all_cram() {
    let temp_dir = tempfile::tempdir().unwrap();
    let cache_dir = temp_dir.path().join("empty-reference-cache");
    fs::create_dir(&cache_dir).unwrap();
    let disguised_cram = temp_dir.path().join("reference-dependent.bam");
    fs::copy(CRAM, &disguised_cram).unwrap();
    let disguised_cram = disguised_cram.to_str().unwrap();

    for (command, transform_args) in [
        ("call-mods", vec!["--filter-threshold", "0.5"]),
        ("adjust-mods", vec!["--ignore", "m"]),
    ] {
        for preexisting_output in [false, true] {
            let output_path = temp_dir.path().join(format!(
                "{command}-{}.sam",
                if preexisting_output { "existing" } else { "absent" }
            ));
            if preexisting_output {
                fs::write(&output_path, b"sentinel\n").unwrap();
            }
            let mut args = vec![
                command,
                disguised_cram,
                output_path.to_str().unwrap(),
                "--output-sam",
                "--threads",
                "1",
                "--suppress-progress",
            ];
            args.extend(transform_args.iter().copied());
            let output =
                run_modkit_without_reference_resolution(&args, &cache_dir);
            assert!(!output.status.success());
            assert!(String::from_utf8_lossy(&output.stderr)
                .contains("failed to decode CRAM input"));
            if preexisting_output {
                assert_eq!(fs::read(&output_path).unwrap(), b"sentinel\n");
            } else {
                assert!(!output_path.exists());
            }
        }
    }

    for (command, transform_args) in [
        ("call-mods", vec!["--filter-threshold", "0.5"]),
        ("adjust-mods", vec!["--ignore", "m"]),
    ] {
        let output_path = temp_dir
            .path()
            .join(format!("reference-independent-{command}.sam"));
        let mut args = vec![
            command,
            REFERENCE_INDEPENDENT_CRAM,
            output_path.to_str().unwrap(),
            "--output-sam",
            "--threads",
            "1",
            "--suppress-progress",
        ];
        args.extend(transform_args);
        let output = run_modkit_without_reference_resolution(&args, &cache_dir);
        assert_success(&output);
        assert!(output_path.exists());
        assert_eq!(normalized_sam_records(&output_path).len(), 10);
    }

    let hemi_output = temp_dir.path().join("missing-ref-hemi.bed");
    let output = run_modkit(&[
        "pileup-hemi",
        CRAM,
        "--out-bed",
        hemi_output.to_str().unwrap(),
        "--cpg",
    ]);
    assert!(!output.status.success());
    assert!(!hemi_output.exists());

    let entropy_output = temp_dir.path().join("missing-ref-entropy.bed");
    let output = run_modkit(&[
        "entropy",
        "--in-bam",
        CRAM,
        "--out-bed",
        entropy_output.to_str().unwrap(),
        "--cpg",
    ]);
    assert!(!output.status.success());
    assert!(!entropy_output.exists());
}
