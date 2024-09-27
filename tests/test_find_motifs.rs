use crate::common::run_modkit;

mod common;

#[test]
fn test_motif_helps() {
    let _ = run_modkit(&["motif", "search", "--help"])
        .expect("failed to run modkit motif search help");
    let _ = run_modkit(&["motif", "evaluate", "--help"])
        .expect("failed to run modkit motif evaluate help");
    let _ = run_modkit(&["motif", "refine", "--help"])
        .expect("failed to run modkit motif refine help");
}
