use std::path::Path;
use std::process::Command;

#[test]
fn pileup_hemi_rejects_an_empty_automatic_threshold_sample_before_output() {
    let temp_dir = tempfile::tempdir().unwrap();
    let output_path = temp_dir.path().join("pileup.bed");
    let executable = Path::new(env!("CARGO_BIN_EXE_modkit"));

    let output = Command::new(executable)
        .args([
            "pileup-hemi",
            "../tests/resources/duplex_modcalls_sort.bam",
            "--out-bed",
            output_path.to_str().unwrap(),
            "--ref",
            "../tests/resources/GRCh38_chr20.fa",
            "--motif",
            "CG",
            "0",
            "--region",
            "chr20:22,613,835-22,640,468",
            "--sampling-frac",
            "0",
            "--suppress-progress",
        ])
        .output()
        .expect("failed to run modkit pileup-hemi");

    assert!(!output.status.success(), "an empty sample must fail");
    assert!(
        String::from_utf8_lossy(&output.stderr).contains(
            "cannot calculate automatic thresholds because no modification \
             probabilities were sampled"
        ),
        "unexpected stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        !output_path.exists(),
        "failure must occur before creating the scientific output"
    );
}
