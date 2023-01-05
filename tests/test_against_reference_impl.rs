use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::process::Output;

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

fn run_mod_flatten(args: &[&str]) -> Output {
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
    output
}

fn test_output(input_path: &str, output_path: &str, check_file_path: &str) {
    let temp_file = std::env::temp_dir().join(output_path);
    let args = [input_path, temp_file.to_str().unwrap()];
    run_mod_flatten(&args);
    assert!(temp_file.exists());

    let mut test_bam = bam::Reader::from_path(temp_file).unwrap();
    let mut ref_bam = bam::Reader::from_path(check_file_path).unwrap();
    for (test_res, ref_res) in test_bam.records().zip(ref_bam.records()) {
        let test_record = test_res.unwrap();
        let ref_record = ref_res.unwrap();
        assert_eq!(ref_record, test_record);
    }
}

#[test]
fn test_canonical() {
    test_output(
        "tests/resources/input_C.bam",
        "test_C.bam",
        "tests/resources/ref_out_C_auto.bam",
    );
}

#[test]
fn test_methyl() {
    test_output(
        "tests/resources/input_5mC.bam",
        "test_5mC.bam",
        "tests/resources/ref_out_5mC_auto.bam",
    );
}

#[test]
fn test_no_tags() {
    let temp_file = std::env::temp_dir().join("test_out_no_tags.bam");
    run_mod_flatten(&[
        "tests/resources/input_C_no_tags.bam",
        temp_file.to_str().unwrap(),
    ]);
}

// #[test]
// fn test_mod_pileup_processor() {
//     let bam_fp = "tests/resources/fwd_rev_modbase_records.sorted.bam";
//     let fasta_fp = "tests/resources/CGI_ladder_3.6kb_ref.fa";
//     let header = bam::Reader::from_path(bam_fp)
//         .map_err(|e| e.to_string())
//         .map(|reader| reader.header().to_owned())
//         .unwrap();
//     let tids = (0..header.target_count())
//         .map(|tid| {
//             let size = header.target_len(tid).unwrap() as u32;
//             let ref_name = String::from_utf8(header.tid2name(tid).to_vec()).unwrap();
//             (tid, size, ref_name)
//         })
//         .collect::<Vec<(u32, u32, String)>>();
//
//     let fasta_fp = "tests/resources/CGI_ladder_3.6kb_ref.fa";
//     let mut fasta_reader = faidx::Reader::from_path(fasta_fp).unwrap();
//     let reference_sequences = tids
//         .iter()
//         .map(|(tid, size, name)| {
//             let size = *size as usize;
//             let seq = fasta_reader.fetch_seq_string(name, 0, size).unwrap();
//             (name.to_string(), seq)
//         })
//         .collect::<HashMap<String, String>>();
//
//     let contig = "oligo_1512_adapters";
//     let dna = reference_sequences.get(contig).unwrap();
//     let motif = CpGModificationMotif;
//     let code = &HydroxyMethylCytosineCode;
//
//     let intervals = IntervalChunks::new(dna.len() as u32, 50, motif.required_overlap());
//     for (start, end) in intervals {
//         dbg!(start, end);
//         let mut processor = ModBasePileupProcessor::new(
//             "tests/resources/fwd_rev_modbase_records.sorted.bam",
//             0,
//             start,
//             end,
//         )
//             .unwrap();
//         let pileup = processor.process_region(code, &motif, contig, dna);
//         let counts = pileup
//             .features
//             .iter_counts(code)
//             .map(|pileup| {
//                 let key = (pileup.position, pileup.strand);
//                 (key, pileup)
//             })
//             .collect::<HashMap<_, _>>();
//     }
// }
