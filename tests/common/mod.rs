use anyhow::{anyhow, Result as AnyhowResult};
use std::path::Path;
use std::process::Output;

pub fn run_modkit(exe: &Path, args: &[&str]) -> AnyhowResult<Output> {
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
