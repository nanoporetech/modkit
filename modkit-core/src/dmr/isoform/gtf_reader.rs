use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum InputKind {
    Plain,
    Bgzip,
    GzipButNotBgzip,
}

/// Detect whether a file is BGZF/bgzip-compressed.
///
/// BGZF is gzip with:
/// - gzip magic bytes: 1f 8b
/// - compression method: 8
/// - FEXTRA flag set
/// - an extra subfield named "BC" with length 2
fn detect_input_kind<P: AsRef<Path>>(path: P) -> io::Result<InputKind> {
    let mut file = File::open(path)?;

    let mut header = [0u8; 12];
    let n = file.read(&mut header)?;

    if n < 2 || header[0] != 0x1f || header[1] != 0x8b {
        return Ok(InputKind::Plain);
    }

    if n < 12 {
        return Ok(InputKind::GzipButNotBgzip);
    }

    let compression_method = header[2];
    let flags = header[3];

    if compression_method != 8 {
        return Ok(InputKind::GzipButNotBgzip);
    }

    // GZIP FEXTRA flag.
    if flags & 0x04 == 0 {
        return Ok(InputKind::GzipButNotBgzip);
    }

    let xlen = u16::from_le_bytes([header[10], header[11]]) as usize;
    let mut extra = vec![0u8; xlen];
    file.read_exact(&mut extra)?;

    let mut i = 0;
    while i + 4 <= extra.len() {
        let si1 = extra[i];
        let si2 = extra[i + 1];
        let slen = u16::from_le_bytes([extra[i + 2], extra[i + 3]]) as usize;
        i += 4;

        if i + slen > extra.len() {
            break;
        }

        if si1 == b'B' && si2 == b'C' && slen == 2 {
            return Ok(InputKind::Bgzip);
        }

        i += slen;
    }

    Ok(InputKind::GzipButNotBgzip)
}

/// Open either a plain-text GTF or a bgzip-compressed GTF as a line reader.
///
/// This intentionally rejects ordinary gzip files that are not BGZF,
/// so `.gz` files made with regular `gzip` do not get silently treated as text.
pub(super) fn open_gtf_reader<P: AsRef<Path>>(
    path: P,
) -> anyhow::Result<Box<dyn BufRead>> {
    let path = path.as_ref();
    let kind = detect_input_kind(path)?;
    let file = File::open(path)?;

    match kind {
        InputKind::Plain => Ok(Box::new(BufReader::new(file))),
        InputKind::Bgzip => {
            // BGZF is a series of gzip-compatible blocks.
            // MultiGzDecoder handles concatenated gzip members.
            let decoder = gzp::BgzfSyncReader::new(file);
            Ok(Box::new(BufReader::new(decoder)))
        }
        InputKind::GzipButNotBgzip => Ok(Box::new(BufReader::new(
            rust_htslib::bgzf::Reader::from_path(path)?,
        ))),
    }
}
