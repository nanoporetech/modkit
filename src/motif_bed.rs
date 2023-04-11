use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use anyhow::{anyhow, Context, Result as AnyhowResult};
use bio::io::fasta::Reader as FastaReader;
use derive_new::new;
use indicatif::{ParallelProgressIterator, ProgressIterator};
use log::debug;
use rayon::prelude::*;
use regex::{Match, Regex};

use crate::util::{
    get_master_progress_bar, get_spinner, ReferenceRecord, Strand,
};

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

struct OverlappingPatternIterator<'a> {
    text: &'a str,
    re: &'a Regex,
    start: usize,
}

impl<'a> Iterator for OverlappingPatternIterator<'a> {
    type Item = Match<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.start >= self.text.len() {
            return None;
        }
        match self.re.find_at(self.text, self.start) {
            Some(m) => {
                self.start = m.start() + 1;
                Some(m)
            }
            None => None,
        }
    }
}

#[derive(Debug)]
pub struct OverlappingRegex {
    inner: Regex,
}

impl OverlappingRegex {
    fn new(pattern: &str) -> Result<Self, regex::Error> {
        Regex::new(pattern).map(|re| Self { inner: re })
    }

    fn find_iter<'a>(
        &'a self,
        text: &'a str,
    ) -> OverlappingPatternIterator<'a> {
        OverlappingPatternIterator {
            text: &text,
            re: &self.inner,
            start: 0,
        }
    }

    fn as_str(&self) -> &str {
        self.inner.as_str()
    }
}

#[derive(Debug, new)]
pub struct RegexMotif {
    forward_pattern: OverlappingRegex,
    reverse_pattern: OverlappingRegex,
    pub forward_offset: usize,
    pub reverse_offset: usize,
    pub length: usize,
}

impl RegexMotif {
    pub fn parse_string(raw_motif: &str, offset: usize) -> AnyhowResult<Self> {
        let length = raw_motif.len();
        let motif = iupac_to_regex(raw_motif);
        let re = OverlappingRegex::new(&motif)?;
        let rc_motif = motif_rev_comp(&motif);
        let rc_re = OverlappingRegex::new(&rc_motif)?;
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

pub fn motif_bed(path: &PathBuf, motif_raw: &str, offset: usize, mask: bool) -> AnyhowResult<()> {
    let motif = iupac_to_regex(&motif_raw);
    let re = OverlappingRegex::new(&motif)?;
    let rc_motif = motif_rev_comp(&motif);
    let rc_re = OverlappingRegex::new(&rc_motif)?;
    let rc_offset = motif_raw
        .len()
        .checked_sub(offset + 1)
        .ok_or(anyhow!("invalid offset for motif"))?;

    let regex_motif =
        RegexMotif::new(re, rc_re, offset, rc_offset, motif_raw.len());

    let file =
        std::fs::File::open(path).context(format!("could not open file"))?;
    let reader = BufReader::new(file);

    let mut seq = String::new();
    let mut header = String::new();
    for line in reader.lines().filter_map(|res| match res {
        Ok(l) => Some(l),
        Err(e) => {
            debug!("could not read line {}", e.to_string());
            None
        }
    }) {
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
    Ok(())
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

        let motif_progress = get_master_progress_bar(seqs_and_target_ids.len());
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
    use crate::motif_bed::{find_motif_hits, RegexMotif};
    use crate::util::Strand;

    #[test]
    fn test_regex_motif() {
        let regex_motif = RegexMotif::parse_string("CCWGG", 1).unwrap();
        assert_eq!(regex_motif.reverse_offset, 3);
        let regex_motif = RegexMotif::parse_string("CG", 0).unwrap();
        assert_eq!(regex_motif.reverse_offset, 1);
    }

    #[test]
    fn test_overlapping_motifs() {
        let regex_motif = RegexMotif::parse_string("CHH", 0).unwrap();
        //               01234567
        let dna = "AACCCCTG";
        //                 CCC
        //                  CCC
        //                   CCT
        let hits = find_motif_hits(&dna, &regex_motif);
        assert_eq!(
            hits,
            vec![
                (2, Strand::Positive),
                (3, Strand::Positive),
                (4, Strand::Positive),
            ]
        );
        let dna = "ACCTAG";
        let hits = find_motif_hits(&dna, &regex_motif);
        assert_eq!(
            hits,
            vec![
                (1, Strand::Positive),
                (2, Strand::Positive),
                (5, Strand::Negative),
            ]
        );
    }
}
