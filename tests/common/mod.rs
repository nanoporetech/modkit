use anyhow::{anyhow, Result as AnyhowResult};
use mod_kit::mod_bam::CollapseMethod;
use mod_kit::summarize::{summarize_modbam, ModSummary};
use mod_kit::threshold_mod_caller::MultipleThresholdModCaller;
use std::fs::File;
use std::io::Read;
use std::path::Path;
use std::process::Output;

pub fn run_modkit(args: &[&str]) -> AnyhowResult<Output> {
    let exe = Path::new(env!("CARGO_BIN_EXE_modkit"));
    assert!(exe.exists());

    let output = std::process::Command::new(exe)
        .args(args)
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .spawn()?
        .wait_with_output()?;
    if output.status.success() {
        Ok(output)
    } else {
        Err(anyhow!("failed to run {:?}", args))
    }
}

fn run_summary<'a>(
    bam_fp: &str,
    interval_size: u32,
    collapse_method: Option<&CollapseMethod>,
) -> AnyhowResult<ModSummary<'a>> {
    let threads = 1usize;
    let pool = rayon::ThreadPoolBuilder::new().num_threads(1).build()?;
    pool.install(|| {
        summarize_modbam(
            &Path::new(bam_fp).to_path_buf(),
            threads,
            interval_size,
            None,
            None,
            None,
            None,
            0.1, // doesn't matter
            Some(MultipleThresholdModCaller::new_passthrough()),
            None,
            collapse_method,
            true,
        )
    })
}

pub fn run_simple_summary(
    bam_fp: &str,
    interval_size: u32,
) -> AnyhowResult<ModSummary> {
    run_summary(bam_fp, interval_size, None)
}

pub fn run_simple_summary_with_collapse_method<'a>(
    bam_fp: &str,
    interval_size: u32,
    collapse_method: &CollapseMethod,
) -> AnyhowResult<ModSummary<'a>> {
    run_summary(bam_fp, interval_size, Some(collapse_method))
}

pub fn check_against_expected_text_file(output_fp: &str, expected_fp: &str) {
    let test = {
        let mut fh = File::open(output_fp).unwrap();
        let mut buff = String::new();
        fh.read_to_string(&mut buff).unwrap();
        buff
    };
    let expected = {
        // this file was hand-checked for correctness.
        let mut fh = File::open(expected_fp).unwrap();
        let mut buff = String::new();
        fh.read_to_string(&mut buff).unwrap();
        buff
    };

    similar_asserts::assert_eq!(test, expected);
}
