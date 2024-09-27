use crate::common::run_modkit;

mod common;

#[test]
fn test_localise_helps() {
    let _ = run_modkit(&["localize", "--help"])
        .expect("failed to run modkit localize help");
    let _ = run_modkit(&["localise", "--help"])
        .expect("failed to run modkit localise help");
}
