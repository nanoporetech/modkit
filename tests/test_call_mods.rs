use crate::common::run_modkit;

mod common;

#[test]
fn test_call_mods_same_positions() {
    let mod_call_out_bam =
        std::env::temp_dir().join("test_call_mods_same_positions_mod_call.bam");
    // run_modkit([
    //     "call-mods",
    //     "data/thresholding/ecoli_reg.sorted.bam"
    // ])
}

#[test]
fn test_call_mods_same_pileup() {}
