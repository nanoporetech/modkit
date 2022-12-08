use rust_htslib::bam;
use rust_htslib::bam::Read;

#[test]
fn test_help() {
    let workdir = std::env::var("CARGO_MANIFEST_DIR").unwrap();
    let exe = std::path::Path::new(&workdir).join("target/debug/mod_flatten");
    assert!(exe.exists());

    let help = std::process::Command::new(exe)
        .arg("--help")
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .spawn()
        .unwrap()
        .wait_with_output()
        .unwrap();
    assert!(help.status.success());
}

fn run_mod_flatten(args: &[&str]) {
    let exe = std::path::Path::new(env!("CARGO_BIN_EXE_mod_flatten"));
    assert!(exe.exists());

    let output = std::process::Command::new(exe)
        .args(args)
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .spawn()
        .unwrap()
        .wait_with_output()
        .unwrap();
    assert!(output.status.success());
}

fn test_output(input_path: &str, output_path: &str) {
    let temp_file = std::env::temp_dir().join(output_path);
    let args = [input_path, temp_file.to_str().unwrap()];
    run_mod_flatten(&args);
    assert!(temp_file.exists());

    let mut test_bam = bam::Reader::from_path(temp_file).unwrap();
    let mut ref_bam = bam::Reader::from_path(input_path).unwrap();
    for (test_res, ref_res) in test_bam.records().zip(ref_bam.records()) {
        let test_record = test_res.unwrap();
        let ref_record = ref_res.unwrap();
        assert_eq!(ref_record, test_record);
    }
}

#[test]
fn test_canonical() {
    test_output("tests/resources/ref_out_C_auto.bam", "test_C.bam");
}

#[test]
fn test_methyl() {
    test_output("tests/resources/ref_out_5mC_auto.bam", "test_5mC.bam");
}
