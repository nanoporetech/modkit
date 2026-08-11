use std::ffi::CString;
use std::io::Write;
use std::path::{Path, PathBuf};

const VALID_ROWS: [&str; 4] = [
    "chr1\t100\t101\tm\t4\t+\t100\t101\t255,0,0\t4 25.00 1 3 0 0 0 0 0\n",
    "chr1\t200\t201\tm\t6\t+\t200\t201\t255,0,0\t6 50.00 3 3 0 0 0 0 0\n",
    "chr1\t300\t301\tm\t8\t-\t300\t301\t255,0,0\t8 50.00 4 4 0 0 0 0 0\n",
    "chr1\t400\t401\tm\t10\t-\t400\t401\t255,0,0\t10 20.00 2 8 0 0 0 0 0\n",
];

const INVALID_ROW: &str =
    "chr1\t200\t201\tm\tnot-a-number\t+\t200\t201\t255,0,0\t6 50.00 3 3 0 0 0 0 0\n";

pub const REGIONS: &str = "chr1\t100\t101\n\
                           chr1\t200\t201\n\
                           chr1\t300\t301\n\
                           chr1\t400\t401\n\
                           absent\t0\t1\n";

pub const GENOME_SIZES: &str = "chr1\t1000\nabsent\t10\n";

pub const OUTPUT_SENTINEL: &[u8] = b"existing output must survive\n";

pub fn make_indexed_bedmethyl(directory: &Path, malformed: bool) -> PathBuf {
    let path = directory.join(if malformed {
        "malformed.bed.gz"
    } else {
        "valid.bed.gz"
    });
    let mut writer = rust_htslib::bgzf::Writer::from_path(&path).unwrap();
    for (index, row) in VALID_ROWS.iter().enumerate() {
        let row = if malformed && index == 1 { INVALID_ROW } else { row };
        writer.write_all(row.as_bytes()).unwrap();
    }
    writer.flush().unwrap();
    drop(writer);

    let path_cstr = CString::new(path.to_str().unwrap()).unwrap();
    let result = unsafe {
        rust_htslib::htslib::tbx_index_build(
            path_cstr.as_ptr(),
            0,
            &rust_htslib::htslib::tbx_conf_bed,
        )
    };
    assert_eq!(result, 0, "failed to index {}", path.display());
    assert!(PathBuf::from(format!("{}.tbi", path.display())).is_file());
    path
}
