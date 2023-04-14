use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use anyhow::{anyhow, Result as AnyhowResult};
use bio::io::fasta::Reader as FastaReader;
use derive_new::new;
use indicatif::{
    ParallelProgressIterator, ProgressBar, ProgressIterator, ProgressStyle,
};
use rayon::prelude::*;
use regex::Regex;

use crate::util::{get_spinner, ReferenceRecord, Strand};

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

#[derive(Debug, new)]
pub struct RegexMotif {
    forward_pattern: Regex,
    reverse_pattern: Regex,
    pub forward_offset: usize,
    pub reverse_offset: usize,
    pub length: usize,
}

impl RegexMotif {
    pub fn parse_string(raw_motif: &str, offset: usize) -> AnyhowResult<Self> {
        let length = raw_motif.len();
        let motif = iupac_to_regex(raw_motif);
        let re = Regex::new(&motif)?;
        let rc_motif = motif_rev_comp(&motif);
        let rc_re = Regex::new(&rc_motif)?;
        let rc_offset = raw_motif
            .len()
            .checked_sub(offset + 1)
            .ok_or(anyhow!("motif not long enough for offset {}", offset))?;
        Ok(Self::new(re, rc_re, offset, rc_offset, length))
    }

    fn is_palendrome(&self) -> bool {
        self.forward_pattern.as_str() == self.reverse_pattern.as_str()
    }
}

fn find_motif_hits(
    seq: &str,
    regex_motif: &RegexMotif,
) -> Vec<(usize, Strand)> {
    let mut motif_hits = vec![];
    // if reverse complement pattern is the same, only search forward pattern
    // and avoid sort
    if regex_motif.is_palendrome() {
        for m in regex_motif.forward_pattern.find_iter(seq) {
            if regex_motif.forward_offset <= regex_motif.reverse_offset {
                motif_hits.push((
                    m.start() + regex_motif.forward_offset,
                    Strand::Positive,
                ));
                motif_hits.push((
                    m.start() + regex_motif.reverse_offset,
                    Strand::Negative,
                ));
            } else {
                motif_hits.push((
                    m.start() + regex_motif.reverse_offset,
                    Strand::Negative,
                ));
                motif_hits.push((
                    m.start() + regex_motif.forward_offset,
                    Strand::Positive,
                ));
            }
        }
    } else {
        for m in regex_motif.forward_pattern.find_iter(seq) {
            motif_hits.push((
                m.start() + regex_motif.forward_offset,
                Strand::Positive,
            ));
        }
        for m in regex_motif.reverse_pattern.find_iter(seq) {
            motif_hits.push((
                m.start() + regex_motif.reverse_offset,
                Strand::Negative,
            ));
        }
        motif_hits.sort_by(|(x_pos, _), (y_pos, _)| x_pos.cmp(&y_pos));
    }

    motif_hits
}

fn process_record(header: &str, seq: &str, regex_motif: &RegexMotif) {
    let motif_hits = find_motif_hits(seq, regex_motif);
    // get contig name
    let (_, rest) = header.split_at(1);
    let ctg = rest.split_whitespace().next().unwrap_or("");
    for (pos, strand) in motif_hits {
        println!("{}\t{}\t{}\t.\t.\t{}", ctg, pos, pos + 1, strand.to_char());
    }
}

pub fn motif_bed(path: &PathBuf, motif_raw: &str, offset: usize, mask: bool) {
    let motif = iupac_to_regex(&motif_raw);
    let re = Regex::new(&motif).unwrap();
    let rc_motif = motif_rev_comp(&motif);
    let rc_re = Regex::new(&rc_motif).unwrap();
    let rc_offset = motif_raw.len().checked_sub(offset + 1).unwrap();

    // todo make a test, then refactor to use parse_str
    let regex_motif =
        RegexMotif::new(re, rc_re, offset, rc_offset, motif_raw.len());

    let file = std::fs::File::open(path).expect("Could not open file");
    let reader = BufReader::new(file);

    let mut seq = String::new();
    let mut header = String::new();
    for line in reader.lines() {
        let line = line.expect("Could not read line");
        if line.starts_with(">") {
            if !seq.is_empty() {
                process_record(&header, &seq, &regex_motif);
            }
            header = line;
            seq.clear();
        } else {
            if mask {
                seq.push_str(&line);
            } else {
                seq.push_str(&line.to_ascii_uppercase())
            }
        }
    }
    if !seq.is_empty() {
        process_record(&header, &seq, &regex_motif);
    }
}

pub struct MotifLocations {
    tid_to_motif_positions: HashMap<u32, HashMap<u32, Strand>>,
    motif: RegexMotif,
}

impl MotifLocations {
    pub fn from_fasta(
        fasta_fp: &PathBuf,
        regex_motif: RegexMotif,
        name_to_tid: &HashMap<&str, u32>,
        mask: bool,
    ) -> AnyhowResult<Self> {
        let reader = FastaReader::from_file(fasta_fp)?;
        let records_progress = get_spinner();
        records_progress.set_message("Reading reference sequences");
        let seqs_and_target_ids = reader
            .records()
            .progress_with(records_progress)
            .filter_map(|res| res.ok())
            .filter_map(|record| {
                name_to_tid.get(record.id()).map(|tid| (record, *tid))
            })
            .filter_map(|(record, tid)| {
                String::from_utf8(record.seq().to_vec())
                    .map(|s| if mask { s } else { s.to_ascii_uppercase() })
                    .ok()
                    .map(|s| (s, tid))
            })
            .collect::<Vec<(String, u32)>>();

        let motif_progress = ProgressBar::new(seqs_and_target_ids.len() as u64);
        let sty = ProgressStyle::with_template(
            "[{elapsed_precise}] {bar:40.red/yellow} {pos:>7}/{len:7} {msg}",
        )
        .unwrap()
        .progress_chars("##-");
        motif_progress.set_style(sty);
        motif_progress.set_message(format!(
            "finding {} motifs",
            regex_motif.forward_pattern.as_str()
        ));
        let tid_to_motif_positions = seqs_and_target_ids
            .into_par_iter()
            .progress_with(motif_progress)
            .map(|(seq, tid)| {
                let positions = find_motif_hits(&seq, &regex_motif)
                    .into_iter()
                    .map(|(pos, strand)| (pos as u32, strand))
                    .collect::<HashMap<u32, Strand>>();
                (tid, positions)
            })
            .collect();

        Ok(Self {
            tid_to_motif_positions,
            motif: regex_motif,
        })
    }

    pub fn filter_reference_records(
        &self,
        reference_records: Vec<ReferenceRecord>,
    ) -> Vec<ReferenceRecord> {
        reference_records
            .into_iter()
            .filter(|target| {
                self.tid_to_motif_positions.contains_key(&target.tid)
            })
            .collect()
    }

    /// will panic if not available, caller needs to check ahead of time
    pub fn get_locations_unchecked(
        &self,
        target_id: u32,
    ) -> &HashMap<u32, Strand> {
        self.tid_to_motif_positions.get(&target_id).unwrap()
    }

    pub fn motif_length(&self) -> usize {
        self.motif.length
    }

    pub fn motif(&self) -> &RegexMotif {
        &self.motif
    }
}

#[cfg(test)]
mod motif_bed_tests {
    use crate::motif_bed::RegexMotif;

    #[test]
    fn test_regex_motif() {
        let regex_motif = RegexMotif::parse_string("CCWGG", 1).unwrap();
        assert_eq!(regex_motif.reverse_offset, 3);
        let regex_motif = RegexMotif::parse_string("CG", 0).unwrap();
        assert_eq!(regex_motif.reverse_offset, 1);
    }
}
