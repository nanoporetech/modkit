use std::collections::HashMap;
use std::fmt::{Display, Formatter};
use std::path::PathBuf;

use anyhow::{anyhow, Context, Result as AnyhowResult};
use bio::io::fasta::Reader as FastaReader;
use derive_new::new;
use indicatif::{MultiProgress, ParallelProgressIterator, ProgressIterator};
use log::debug;
use rayon::prelude::*;
use regex::{Match, Regex};
use rustc_hash::{FxHashMap, FxHashSet};

use crate::position_filter::StrandedPositionFilter;
use crate::util::{
    get_master_progress_bar, get_spinner, get_ticker, ReferenceRecord, Strand,
    StrandRule,
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
    pub raw_motif: String,
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
        Ok(Self::new(
            re,
            rc_re,
            offset,
            rc_offset,
            length,
            raw_motif.to_owned(),
        ))
    }

    pub(crate) fn is_palendrome(&self) -> bool {
        self.forward_pattern.as_str() == self.reverse_pattern.as_str()
    }

    pub(crate) fn offset(&self) -> i32 {
        self.reverse_offset as i32 - self.forward_offset as i32
    }

    pub(crate) fn negative_strand_position(
        &self,
        positive_position: u32,
    ) -> Option<u32> {
        if !self.is_palendrome() {
            None
        } else {
            let pos = positive_position as i64;
            let offset = self.offset() as i64;
            let adj = pos + offset;
            if adj < 0 {
                None
            } else {
                Some(adj as u32)
            }
        }
    }
}

impl Display for RegexMotif {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{},{}", self.raw_motif, self.forward_offset)
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

fn process_record(header: &str, seq: &str, regex_motif: &RegexMotif) -> usize {
    let motif_hits = find_motif_hits(seq, regex_motif);
    let n_hits = motif_hits.len();
    for (pos, strand) in motif_hits {
        println!(
            "{}\t{}\t{}\t.\t.\t{}",
            header,
            pos,
            pos + 1,
            strand.to_char()
        );
    }
    n_hits
}

pub fn motif_bed(
    path: &PathBuf,
    motif_raw: &str,
    offset: usize,
    mask: bool,
) -> AnyhowResult<()> {
    let motif = iupac_to_regex(&motif_raw);
    let re = OverlappingRegex::new(&motif)
        .context("failed to make forward regex pattern")?;
    let rc_motif = motif_rev_comp(&motif);
    let rc_re = OverlappingRegex::new(&rc_motif)
        .context("failed to make reverse complement regex pattern")?;
    let rc_offset = motif_raw
        .len()
        .checked_sub(offset + 1)
        .ok_or(anyhow!("invalid offset for motif"))?;

    let regex_motif = RegexMotif::new(
        re,
        rc_re,
        offset,
        rc_offset,
        motif_raw.len(),
        motif_raw.to_owned(),
    );

    let reader =
        FastaReader::from_file(path).context("failed to open FASTA")?;

    // prog bar stuff
    let master_pb = MultiProgress::new();
    let records_progress = master_pb.add(get_spinner());
    records_progress.set_message("Reading reference sequences");
    let motifs_progress = master_pb.add(get_ticker());
    motifs_progress.set_message(format!(
        "{} motifs found",
        regex_motif.forward_pattern.as_str()
    ));
    // end prog bar stuff

    reader
        .records()
        .progress_with(records_progress)
        .filter_map(|r| match r {
            Ok(r) => Some(r),
            Err(e) => {
                debug!("failed to read record, {}", e.to_string());
                None
            }
        })
        .filter_map(|record| {
            let seq = String::from_utf8(record.seq().to_vec());
            match seq {
                Ok(s) => Some((s, record)),
                Err(e) => {
                    let header = record.id();
                    debug!(
                        "sequence for {header} failed UTF-8 conversion, {}",
                        e.to_string()
                    );
                    None
                }
            }
        })
        .for_each(|(seq, record)| {
            let seq = if mask { seq } else { seq.to_ascii_uppercase() };
            let n_hits = process_record(record.id(), &seq, &regex_motif);
            motifs_progress.inc(n_hits as u64);
        });

    motifs_progress.finish_and_clear();
    Ok(())
}

pub struct MultipleMotifLocations {
    pub(crate) motif_locations: Vec<MotifLocations>,
    /// mapping of target to mapping of position to vector of which motif locations
    /// (in `self.motif_locations`) have a hit at that position
    position_lookup: FxHashMap<u32, FxHashMap<(u32, Strand), Vec<usize>>>,
}

impl MultipleMotifLocations {
    pub fn new(motif_locations: Vec<MotifLocations>) -> Self {
        let position_lookup = motif_locations.iter().enumerate().fold(
            FxHashMap::<u32, FxHashMap<(u32, Strand), Vec<usize>>>::default(),
            |mut acc, (idx, mls)| {
                mls.tid_to_motif_positions.iter().for_each(
                    |(target_id, positions)| {
                        let positions_for_target = acc
                            .entry(*target_id)
                            .or_insert(FxHashMap::default());
                        positions.iter().for_each(|(position, strand_rule)| {
                            match strand_rule {
                                StrandRule::Positive => {
                                    let k = (*position, Strand::Positive);
                                    positions_for_target
                                        .entry(k)
                                        .or_insert(Vec::new())
                                        .push(idx)
                                }
                                StrandRule::Negative => {
                                    let k = (*position, Strand::Negative);
                                    positions_for_target
                                        .entry(k)
                                        .or_insert(Vec::new())
                                        .push(idx)
                                }
                                StrandRule::Both => {
                                    for k in [
                                        (*position, Strand::Positive),
                                        (*position, Strand::Negative),
                                    ] {
                                        positions_for_target
                                            .entry(k)
                                            .or_insert(Vec::new())
                                            .push(idx)
                                    }
                                }
                            }
                        });
                    },
                );
                acc
            },
        );

        Self {
            motif_locations,
            position_lookup,
        }
    }

    pub fn motifs_for_position(
        &self,
        target_id: u32,
        position: u32,
        strand: Strand,
    ) -> Option<Vec<(usize, &MotifLocations)>> {
        self.position_lookup
            .get(&target_id)
            .and_then(|positions| positions.get(&(position, strand)))
            .map(|idxs| {
                idxs.iter()
                    .filter_map(|&i| {
                        self.motif_locations.get(i).map(|ml| (i, ml))
                    })
                    .collect()
            })
    }

    pub fn motif_idxs_for_position(
        &self,
        target_id: u32,
        position: u32,
        strand: Strand,
    ) -> Option<&Vec<usize>> {
        self.position_lookup
            .get(&target_id)
            .and_then(|positions| positions.get(&(position, strand)))
    }
}

#[derive(Debug)]
pub struct MotifLocations {
    tid_to_motif_positions: FxHashMap<u32, FxHashMap<u32, StrandRule>>,
    motif: RegexMotif,
}

impl MotifLocations {
    pub fn from_fasta(
        fasta_fp: &PathBuf,
        regex_motif: RegexMotif,
        name_to_tid: &HashMap<&str, u32>,
        mask: bool,
        position_filter: Option<&StrandedPositionFilter>,
        master_progress_bar: &MultiProgress,
    ) -> AnyhowResult<Self> {
        let reader = FastaReader::from_file(fasta_fp)?;

        let records_progress = master_progress_bar.add(get_spinner());
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

        let motif_progress = master_progress_bar
            .add(get_master_progress_bar(seqs_and_target_ids.len()));
        motif_progress.set_message(format!("finding {} motifs", regex_motif));
        let tid_to_motif_positions = seqs_and_target_ids
            .into_par_iter()
            .progress_with(motif_progress)
            .map(|(seq, tid)| {
                let positions = find_motif_hits(&seq, &regex_motif)
                    .into_iter() // todo into_par_iter?
                    .filter_map(|(pos, strand)| {
                        if let Some(position_filter) = position_filter {
                            if position_filter
                                .contains(tid as i32, pos as u64, strand)
                            {
                                Some((pos as u32, strand))
                            } else {
                                None
                            }
                        } else {
                            Some((pos as u32, strand))
                        }
                    })
                    .fold(
                        FxHashMap::<u32, StrandRule>::default(),
                        |mut acc, (pos, strand)| {
                            if let Some(strand_rule) = acc.get_mut(&pos) {
                                *strand_rule = strand_rule.absorb(strand);
                            } else {
                                acc.insert(pos, strand.into());
                            }
                            acc
                        },
                    );
                (tid, positions)
            })
            .collect();

        Ok(Self {
            tid_to_motif_positions,
            motif: regex_motif,
        })
    }

    pub(crate) fn references_with_hits(&self) -> FxHashSet<u32> {
        self.tid_to_motif_positions.keys().copied().collect()
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
    ) -> &FxHashMap<u32, StrandRule> {
        self.tid_to_motif_positions.get(&target_id).unwrap()
    }

    pub fn motif_length(&self) -> usize {
        self.motif.length
    }

    pub fn motif(&self) -> &RegexMotif {
        &self.motif
    }

    pub fn targets_to_positions(
        &self,
    ) -> &FxHashMap<u32, FxHashMap<u32, StrandRule>> {
        &self.tid_to_motif_positions
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

        let motif = RegexMotif::parse_string("CGCG", 2).unwrap();
        assert_eq!(motif.offset(), -1);
    }

    #[test]
    fn test_motif_hits() {
        let seq = "AACGCGAACGCGA";
        let motif = RegexMotif::parse_string("CGCG", 2).unwrap();
        assert_eq!(motif.offset(), -1);
        let hits = find_motif_hits(seq, &motif);
        let expected = vec![
            (3, Strand::Negative),
            (4, Strand::Positive),
            (9, Strand::Negative),
            (10, Strand::Positive),
        ];
        assert_eq!(&hits, &expected);
        for pos in hits
            .iter()
            .filter(|(_, strand)| *strand == Strand::Positive)
            .map(|(p, _)| *p as u32)
        {
            let negative_strand_pos = motif
                .negative_strand_position(pos)
                .expect("should find position");
            let _found = hits
                .iter()
                .find(|(p, strand)| {
                    *p as u32 == negative_strand_pos
                        && *strand == Strand::Negative
                })
                .expect("should find negative strand position");
        }
        assert!(motif.negative_strand_position(0).is_none())
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

    #[test]
    fn test_motif_palendrome() {
        let chh = RegexMotif::parse_string("CHH", 0).unwrap();
        assert!(!chh.is_palendrome());
        let cg = RegexMotif::parse_string("CG", 0).unwrap();
        assert!(cg.is_palendrome());
        let c = RegexMotif::parse_string("C", 0).unwrap();
        assert!(!c.is_palendrome());
        let gatc = RegexMotif::parse_string("GATC", 1).unwrap();
        assert!(gatc.is_palendrome());
    }
}
