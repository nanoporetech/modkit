use std::collections::{BTreeMap, HashMap};
use std::fmt::{Display, Formatter};
use std::path::PathBuf;

use anyhow::{anyhow, bail, Context, Result as AnyhowResult};
use derive_new::new;
use indicatif::{MultiProgress, ProgressIterator};
use itertools::Itertools;
use log::{debug, info};
use rayon::prelude::*;
use regex::{Match, Regex};
use rustc_hash::FxHashMap;

use crate::fasta::HtsLibFastaRecords;
use crate::monoid::Moniod;
use crate::position_filter::StrandedPositionFilter;
use crate::util::{
    get_master_progress_bar, get_spinner, get_ticker, ReferenceRecord, Strand,
    StrandRule,
};

fn iupac_to_regex(pattern: &str) -> anyhow::Result<String> {
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
            _ => bail!("Invalid IUPAC code: {}", c),
        });
    }
    Ok(regex)
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

pub(crate) struct OverlappingPatternIterator<'a> {
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

#[derive(Debug, Clone)]
pub struct OverlappingRegex {
    inner: Regex,
}

impl OverlappingRegex {
    pub(crate) fn as_str(&self) -> &str {
        self.inner.as_str()
    }

    fn new(pattern: &str) -> Result<Self, regex::Error> {
        Regex::new(pattern).map(|re| Self { inner: re })
    }

    pub(crate) fn find_iter<'a>(
        &'a self,
        text: &'a str,
    ) -> OverlappingPatternIterator<'a> {
        OverlappingPatternIterator { text: &text, re: &self.inner, start: 0 }
    }
}

#[derive(Debug, new, Copy, Clone)]
pub struct MotifInfo {
    pub forward_offset: usize,
    pub reverse_offset: usize,
    pub length: usize,
    pub is_palendrome: bool,
}

impl MotifInfo {
    pub(crate) fn offset(&self) -> i32 {
        self.reverse_offset as i32 - self.forward_offset as i32
    }

    pub(crate) fn negative_strand_position(
        &self,
        positive_position: u32,
    ) -> Option<u32> {
        if !self.is_palendrome {
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

#[derive(Debug, Clone, new)]
pub struct RegexMotif {
    pub(crate) forward_pattern: OverlappingRegex,
    pub(crate) reverse_pattern: OverlappingRegex,
    pub motif_info: MotifInfo,
    pub raw_motif: String,
}

impl RegexMotif {
    pub fn from_raw_parts(
        raw_motif_parts: &[String],
        cpg: bool,
    ) -> AnyhowResult<Vec<Self>> {
        if raw_motif_parts
            .chunks(2)
            .map(|chunk| (&chunk[0], &chunk[1]))
            .counts()
            .values()
            .max()
            .unwrap_or(&1)
            > &1
        {
            bail!("cannot have the same motif more than once")
        }
        let mut raw_motif_parts = raw_motif_parts.to_owned();
        if cpg {
            if raw_motif_parts.chunks(2).any(|motif| motif == ["CG", "0"]) {
                info!(
                    "CG 0 motif already, don't need --cpg and --motif CG 0, \
                     ignoring --cpg"
                );
            } else {
                info!("--cpg flag received, adding CG, 0 to motifs");
                raw_motif_parts
                    .extend_from_slice(&["CG".to_string(), "0".to_string()]);
            }
        }
        raw_motif_parts
            .chunks(2)
            .map(|c| {
                let motif = &c[0];
                let focus_base = &c[1];
                focus_base
                    .parse::<usize>()
                    .map_err(|e| {
                        anyhow!("couldn't parse focus base, {}", e.to_string())
                    })
                    .and_then(|focus_base| {
                        RegexMotif::parse_string(motif.as_str(), focus_base)
                    })
            })
            .collect::<Result<Vec<RegexMotif>, anyhow::Error>>()
    }

    pub fn parse_string(raw_motif: &str, offset: usize) -> AnyhowResult<Self> {
        let length = raw_motif.len();
        if length == 1 {
            match raw_motif {
                "A" | "C" | "G" | "T" => {}
                _ => bail!(
                    "degenerate bases are not supported as single base \
                     motifs, must be 'A', 'C', 'G', or 'T'."
                ),
            };
        };
        let motif = iupac_to_regex(raw_motif)?;
        let re = OverlappingRegex::new(&motif)?;
        let rc_motif = motif_rev_comp(&motif);
        let rc_re = OverlappingRegex::new(&rc_motif)?;
        let rc_offset = raw_motif
            .len()
            .checked_sub(offset + 1)
            .ok_or(anyhow!("motif not long enough for offset {}", offset))?;
        let motif_info = MotifInfo::new(
            offset,
            rc_offset,
            length,
            re.as_str() == rc_re.as_str(),
        );
        Ok(Self::new(re, rc_re, motif_info, raw_motif.to_owned()))
    }

    pub(crate) fn is_palendrome(&self) -> bool {
        self.forward_pattern.as_str() == self.reverse_pattern.as_str()
    }

    #[inline(always)]
    pub(crate) fn length(&self) -> usize {
        self.motif_info.length
    }

    #[inline(always)]
    pub(crate) fn forward_offset(&self) -> usize {
        self.motif_info.forward_offset
    }

    #[inline(always)]
    pub(crate) fn reverse_offset(&self) -> usize {
        self.motif_info.reverse_offset
    }

    #[cfg(test)]
    fn offset(&self) -> i32 {
        self.motif_info.offset()
    }

    pub(crate) fn find_hits(&self, seq: &str) -> Vec<(usize, Strand)> {
        find_motif_hits(seq, &self)
    }
}

impl Display for RegexMotif {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{},{}", self.raw_motif, self.forward_offset())
    }
}

pub(crate) fn find_single_bases(
    seq: &str,
    regex_motif: &RegexMotif,
) -> Vec<(usize, Strand)> {
    let haystack = seq.as_bytes();
    let (fw, rv) = match regex_motif.forward_pattern.as_str() {
        "A" => ('A', 'T'),
        "C" => ('C', 'G'),
        "G" => ('G', 'C'),
        "T" => ('T', 'A'),
        // todo refactor this into a compile-time check
        _ => unreachable!(
            "RegexMotif cannot be constructed from non DNA-base single letter \
             motifs"
        ),
    };
    let (fw, rv) = (fw as u8, rv as u8);
    memchr::memchr2_iter(fw, rv, haystack)
        .map(|pos| {
            if haystack[pos] == fw {
                (pos, Strand::Positive)
            } else {
                (pos, Strand::Negative)
            }
        })
        .collect()
}

pub(crate) fn find_motif_hits(
    seq: &str,
    regex_motif: &RegexMotif,
) -> Vec<(usize, Strand)> {
    let mut motif_hits = vec![];
    // if reverse complement pattern is the same, only search forward pattern
    // and avoid sort
    if regex_motif.is_palendrome() {
        for m in regex_motif.forward_pattern.find_iter(seq) {
            if regex_motif.forward_offset() <= regex_motif.reverse_offset() {
                motif_hits.push((
                    m.start() + regex_motif.forward_offset(),
                    Strand::Positive,
                ));
                motif_hits.push((
                    m.start() + regex_motif.reverse_offset(),
                    Strand::Negative,
                ));
            } else {
                motif_hits.push((
                    m.start() + regex_motif.reverse_offset(),
                    Strand::Negative,
                ));
                motif_hits.push((
                    m.start() + regex_motif.forward_offset(),
                    Strand::Positive,
                ));
            }
        }
    } else if regex_motif.length() == 1 {
        let mut single_base_sites = find_single_bases(seq, regex_motif);
        motif_hits.append(&mut single_base_sites);
    } else {
        for m in regex_motif.forward_pattern.find_iter(seq) {
            motif_hits.push((
                m.start() + regex_motif.forward_offset(),
                Strand::Positive,
            ));
        }
        for m in regex_motif.reverse_pattern.find_iter(seq) {
            motif_hits.push((
                m.start() + regex_motif.reverse_offset(),
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
    let motif = iupac_to_regex(&motif_raw)?;
    let re = OverlappingRegex::new(&motif)
        .context("failed to make forward regex pattern")?;
    let rc_motif = motif_rev_comp(&motif);
    let rc_re = OverlappingRegex::new(&rc_motif)
        .context("failed to make reverse complement regex pattern")?;
    let rc_offset = motif_raw
        .len()
        .checked_sub(offset + 1)
        .ok_or(anyhow!("invalid offset for motif"))?;

    let motif_info = MotifInfo::new(
        offset,
        rc_offset,
        motif_raw.len(),
        re.as_str() == rc_re.as_str(),
    );
    let regex_motif =
        RegexMotif::new(re, rc_re, motif_info, motif_raw.to_owned());

    let reader =
        HtsLibFastaRecords::from_file(path).context("failed to open FASTA")?;

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
        .into_iter()
        .progress_with(records_progress)
        .filter_map(|r| match r {
            Ok(r) => Some(r),
            Err(e) => {
                debug!("failed to read record, {}", e.to_string());
                None
            }
        })
        .for_each(|(read_id, seq)| {
            let seq = if mask { seq } else { seq.to_ascii_uppercase() };
            let n_hits = process_record(&read_id, &seq, &regex_motif);
            motifs_progress.inc(n_hits as u64);
        });

    motifs_progress.finish_and_clear();
    Ok(())
}

/// A wrapper for a collection of MotifLocations
pub struct MultipleMotifLocations {
    pub(crate) motif_locations: Vec<MotifLocations>,
}

impl MultipleMotifLocations {
    pub fn new(motif_locations: Vec<MotifLocations>) -> Self {
        Self { motif_locations }
    }

    pub fn motifs_at_position(
        &self,
        target_id: u32,
        position: u32,
    ) -> Vec<&RegexMotif> {
        self.motif_locations
            .iter()
            .filter(|ml| {
                ml.tid_to_motif_positions
                    .get(&target_id)
                    .and_then(|pos| pos.get(&position))
                    .is_some()
            })
            .map(|ml| &ml.motif)
            .collect()
    }

    fn indices_and_motifs_at_position(
        &self,
        target_id: u32,
        position: u32,
        strand: Strand,
    ) -> Vec<(usize, &MotifLocations)> {
        self.motif_locations
            .iter()
            .enumerate()
            .filter(|(_idx, mls)| {
                mls.tid_to_motif_positions
                    .get(&target_id)
                    .and_then(|positions| positions.get(&position))
                    .map(|strand_rule| strand_rule.covers(strand))
                    .unwrap_or(false)
            })
            .collect::<Vec<(usize, &MotifLocations)>>()
    }

    // todo refactor this so that the Option is unnecessary
    pub fn motifs_at_position_nonempty(
        &self,
        target_id: u32,
        position: u32,
        strand: Strand,
    ) -> Option<Vec<(usize, &MotifLocations)>> {
        let x =
            self.indices_and_motifs_at_position(target_id, position, strand);
        if x.is_empty() {
            None
        } else {
            Some(x)
        }
    }

    // todo refactor this so that the Option is unnecessary
    pub fn motif_idx_at_position_nonempty(
        &self,
        target_id: u32,
        position: u32,
        strand: Strand,
    ) -> Option<Vec<usize>> {
        let idxs = self
            .indices_and_motifs_at_position(target_id, position, strand)
            .into_iter()
            .map(|(idx, _)| idx)
            .collect::<Vec<usize>>();
        if idxs.is_empty() {
            None
        } else {
            Some(idxs)
        }
    }
}

pub fn get_masked_sequences(
    fasta_fp: &PathBuf,
    name_to_tid: &HashMap<&str, u32>,
    mask: bool,
    master_progress_bar: &MultiProgress,
) -> anyhow::Result<Vec<(String, u32)>> {
    let reader = HtsLibFastaRecords::from_file(fasta_fp)?;

    let records_progress = master_progress_bar.add(get_ticker());
    records_progress.set_message("Reading reference sequences");

    Ok(reader
        .into_iter()
        .progress_with(records_progress)
        .filter_map(|res| res.ok())
        .map(|(read_id, read_sequence)| {
            let seq = if mask {
                read_sequence
            } else {
                read_sequence.to_ascii_uppercase()
            };
            (read_id, seq)
        })
        .filter_map(|(read_id, record)| {
            name_to_tid.get(read_id.as_str()).map(|tid| (record, *tid))
        })
        .collect::<Vec<(String, u32)>>())
}

#[derive(Debug, new)]
pub struct MotifLocations {
    pub tid_to_motif_positions: FxHashMap<u32, BTreeMap<u32, StrandRule>>,
    pub motif: RegexMotif,
}

impl MotifLocations {
    pub fn from_sequences(
        regex_motif: RegexMotif,
        position_filter: Option<&StrandedPositionFilter<()>>,
        sequences_and_ids: &[(String, u32)],
        master_progress_bar: &MultiProgress,
    ) -> AnyhowResult<Self> {
        let motif_progress = master_progress_bar
            .add(get_master_progress_bar(sequences_and_ids.len()));
        motif_progress.set_message(format!("finding {} motifs", regex_motif));
        let sequences_and_ids = sequences_and_ids
            .iter()
            .sorted_by(|(s, _), (p, _)| s.len().cmp(&p.len()))
            .collect::<Vec<_>>();

        let tid_to_motif_positions = sequences_and_ids
            .into_par_iter()
            .filter(|(_seq, tid)| {
                position_filter
                    .map(|x| x.contains_chrom_id(&(*tid as i64)))
                    .unwrap_or(true)
            })
            .map(|(seq, tid)| {
                let now = std::time::Instant::now();
                let positions = find_motif_hits(&seq, &regex_motif)
                    .into_par_iter()
                    .filter_map(|(pos, strand)| {
                        if let Some(position_filter) = position_filter {
                            if position_filter.contains(
                                *tid as i32,
                                pos as u64,
                                strand,
                            ) {
                                Some((pos as u32, strand))
                            } else {
                                None
                            }
                        } else {
                            Some((pos as u32, strand))
                        }
                    })
                    .fold(
                        || BTreeMap::<u32, StrandRule>::new(),
                        |mut acc, (pos, strand)| {
                            if let Some(strand_rule) = acc.get_mut(&pos) {
                                *strand_rule = strand_rule.absorb(strand);
                            } else {
                                acc.insert(pos, strand.into());
                            }
                            acc
                        },
                    )
                    .reduce(
                        || BTreeMap::<u32, StrandRule>::new(),
                        |a, b| a.into_iter().chain(b).collect(),
                    );
                let duration = now.elapsed().as_millis() as f64 / 1000f64;
                let rate = seq.len() as f64 / (duration * 1000_000f64);
                debug!(
                    "motif {} has {} positions in tid {}, took {:.4}s ({rate} \
                     kb/msec)",
                    &regex_motif,
                    positions.len(),
                    tid,
                    duration
                );

                motif_progress.inc(1);
                (*tid, positions)
            })
            .collect();

        Ok(Self { tid_to_motif_positions, motif: regex_motif })
    }

    pub fn from_fasta(
        fasta_fp: &PathBuf,
        regex_motif: RegexMotif,
        name_to_tid: &HashMap<&str, u32>,
        mask: bool,
        position_filter: Option<&StrandedPositionFilter<()>>,
        master_progress_bar: &MultiProgress,
    ) -> AnyhowResult<Self> {
        let seqs_and_target_ids = get_masked_sequences(
            &fasta_fp,
            &name_to_tid,
            mask,
            master_progress_bar,
        )?;
        Self::from_sequences(
            regex_motif,
            position_filter,
            &seqs_and_target_ids,
            master_progress_bar,
        )
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
    ) -> &BTreeMap<u32, StrandRule> {
        self.tid_to_motif_positions.get(&target_id).unwrap()
    }

    pub fn motif_length(&self) -> usize {
        self.motif.length()
    }

    pub fn motif(&self) -> &RegexMotif {
        &self.motif
    }
}

pub(crate) struct MotifPositionLookup {
    inner: FxHashMap<(u32, usize, Strand), Vec<usize>>,
    motifs: Vec<RegexMotif>,
}

impl MotifPositionLookup {
    pub(crate) fn new(
        tid_motif_positions: FxHashMap<(u32, usize), Vec<(usize, Strand)>>,
        motifs: Vec<RegexMotif>,
    ) -> Self {
        let inner = tid_motif_positions
            .into_par_iter()
            .fold(
                || FxHashMap::default(),
                |mut agg, ((tid, motif_idx), positions)| {
                    for (pos, strand) in positions {
                        agg.entry((tid, pos, strand))
                            .or_insert_with(Vec::new)
                            .push(motif_idx);
                    }
                    agg
                },
            )
            .reduce(|| FxHashMap::default(), |a, b| a.op(b));

        Self { inner, motifs }
    }

    pub(crate) fn get_motif_hits(
        &self,
        tid: u32,
        position: usize,
        strand: Strand,
    ) -> Option<String> {
        let k = (tid, position, strand);
        self.inner.get(&k).map(|ixs| {
            ixs.iter().map(|i| self.motifs[*i].to_string()).join(";")
        })
    }
}

#[cfg(test)]
mod motif_bed_tests {
    use crate::motifs::motif_bed::{find_motif_hits, RegexMotif};
    use crate::util::Strand;

    #[test]
    fn test_regex_motif() {
        let regex_motif = RegexMotif::parse_string("CCWGG", 1).unwrap();
        assert_eq!(regex_motif.forward_offset(), 1);
        assert_eq!(regex_motif.reverse_offset(), 3);
        assert_eq!(regex_motif.length(), 5);
        let regex_motif = RegexMotif::parse_string("CG", 0).unwrap();
        assert_eq!(regex_motif.reverse_offset(), 1);

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
                .motif_info
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
        assert!(motif.motif_info.negative_strand_position(0).is_none())
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
    fn test_motif_palindrome() {
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
