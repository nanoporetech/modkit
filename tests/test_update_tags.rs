use crate::common::run_modkit;
use rust_htslib::bam::{self, record::Aux, Read};

mod common;

#[test]
fn test_update_tags_implicit_no_probabilities() {
    let temp_file = std::env::temp_dir()
        .join("test_update_tags_implicit_no_probabilities_updated_tags.bam");
    let mut reader = bam::Reader::from_path(
        "tests/resources/single_read_old_tags_no_probs.bam",
    )
    .unwrap();
    let record = reader.records().next().map(|r| r.unwrap()).unwrap();
    let old_tags = match record.aux(&['M' as u8, 'M' as u8]).unwrap() {
        Aux::String(raw_mm) => raw_mm.to_string(),
        _ => panic!("wrong type"),
    };
    assert_eq!(old_tags, "C+h;C+m;");

    run_modkit(&[
        "update-tags",
        "tests/resources/single_read_old_tags_no_probs.bam",
        temp_file.to_str().unwrap(),
    ])
    .expect("should run update-tags");

    let mut reader =
        bam::Reader::from_path(temp_file.to_str().unwrap()).unwrap();
    let record = reader.records().next().map(|r| r.unwrap()).unwrap();
    let old_tags = match record.aux(&['M' as u8, 'M' as u8]).unwrap() {
        Aux::String(raw_mm) => raw_mm.to_string(),
        _ => panic!("wrong type"),
    };
    assert_eq!(old_tags, "C+h.;C+m.;");
}
