use regex::Regex;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

fn iupac_to_regex(pattern: &str) -> String {
    let mut regex = String::new();
    for c in pattern.chars() {
        regex.push_str(match c {
            'A' => "A",
            'C' => "C",
            'G' => "G",
            'T' => "T",
            'U' => "U",
            'M' => "[AC]",
            'R' => "[AG]",
            'W' => "[AT]",
            'S' => "[CG]",
            'Y' => "[CT]",
            'K' => "[GT]",
            'V' => "[ACG]",
            'H' => "[ACT]",
            'D' => "[AGT]",
            'B' => "[CGT]",
            'X' => "[ACGT]",
            'N' => "[ACGT]",
            _ => panic!("Invalid IUPAC code: {}", c),
        });
    }
    regex
}

fn motif_rev_comp(motif: &str) -> String {
    let mut reverse_complement = motif.chars().rev().collect::<String>();
    reverse_complement = reverse_complement
        .chars()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            'U' => 'A',
            '[' => ']',
            ']' => '[',
            _ => c,
        })
        .collect();
    reverse_complement
}

fn process_record(
    header: &str,
    seq: &str,
    re: &Regex,
    offset: usize,
    rc_re: &Regex,
    rc_offset: usize,
) {
    let mut motif_hits = vec![];
    // if reverse complement pattern is the same, only search forward pattern
    // and avoid sort
    if re.as_str() == rc_re.as_str() {
        for m in re.find_iter(seq) {
            if offset <= rc_offset {
                motif_hits.push((m.start() + offset, "+"));
                motif_hits.push((m.start() + rc_offset, "-"));
            } else {
                motif_hits.push((m.start() + rc_offset, "-"));
                motif_hits.push((m.start() + offset, "+"));
            }
        }
    } else {
        for m in re.find_iter(seq) {
            motif_hits.push((m.start() + offset, "+"));
        }
        for m in rc_re.find_iter(seq) {
            motif_hits.push((m.start() + rc_offset, "-"));
        }
        motif_hits.sort();
    }

    // get contig name
    let (_, rest) = header.split_at(1);
    let ctg = rest.split_whitespace().next().unwrap_or("");
    for (pos, strand) in motif_hits {
        println!("{}\t{}\t{}\t.\t.\t{}", ctg, pos, pos + 1, strand);
    }
}

pub fn motif_bed(path: &PathBuf, motif_raw: &str, offset: usize) {
    let motif = iupac_to_regex(&motif_raw);
    let re = Regex::new(&motif).unwrap();
    let rc_motif = motif_rev_comp(&motif);
    let rc_re = Regex::new(&rc_motif).unwrap();
    let rc_offset = motif_raw.len().checked_sub(offset + 1).unwrap();

    let file = std::fs::File::open(path).expect("Could not open file");
    let reader = BufReader::new(file);

    let mut seq = String::new();
    let mut header = String::new();
    for line in reader.lines() {
        let line = line.expect("Could not read line");
        if line.starts_with(">") {
            if !seq.is_empty() {
                process_record(&header, &seq, &re, offset, &rc_re, rc_offset);
            }
            header = line;
            seq.clear();
        } else {
            seq.push_str(&line);
        }
    }
    if !seq.is_empty() {
        process_record(&header, &seq, &re, offset, &rc_re, rc_offset);
    }
}
