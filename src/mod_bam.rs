use crate::errs::{InputError, RunError};
use crate::mod_base_code::{DnaBase, ModCode};
use crate::util;
use crate::util::{get_tag, record_is_secondary, Strand};

use std::cmp::Ordering;

use crate::position_filter::StrandedPositionFilter;
use derive_new::new;
use itertools::Itertools;
use log::debug;
use nom::bytes::complete::tag;
use nom::character::complete::{digit1, multispace0};
use nom::multi::separated_list1;
use nom::IResult;
use rust_htslib::bam;
use rust_htslib::bam::record::Aux;
use rustc_hash::FxHashMap;
use std::collections::{HashMap, HashSet};

pub(crate) struct TrackingModRecordIter<'a, T: bam::Read> {
    records: bam::Records<'a, T>,
    skip_unmapped: bool,
    pub(crate) num_used: usize,
    pub(crate) num_skipped: usize,
    pub(crate) num_failed: usize,
}

impl<'a, T: bam::Read> TrackingModRecordIter<'a, T> {
    pub(crate) fn new(
        records: bam::Records<'a, T>,
        skip_unmapped: bool,
    ) -> Self {
        Self {
            records,
            skip_unmapped,
            num_used: 0,
            num_skipped: 0,
            num_failed: 0,
        }
    }
}

impl<'a, T: bam::Read> Iterator for &mut TrackingModRecordIter<'a, T> {
    type Item = (bam::Record, String, ModBaseInfo);

    fn next(&mut self) -> Option<Self::Item> {
        let mut ret: Option<Self::Item> = None;
        while let Some(result) = self.records.next() {
            match result {
                Ok(record) => {
                    let record_name =
                        String::from_utf8(record.qname().to_vec())
                            .unwrap_or("utf-decode-failed".to_string());
                    if record_is_secondary(&record)
                        || (record.is_unmapped() && self.skip_unmapped)
                    {
                        self.num_skipped += 1;
                        continue;
                    } else {
                        if record.seq_len() == 0 {
                            debug!(
                                "record {record_name} has zero length sequence"
                            );
                            self.num_failed += 1;
                            continue;
                        } else {
                            match ModBaseInfo::new_from_record(&record) {
                                Ok(modbase_info) => {
                                    if modbase_info.is_empty() {
                                        self.num_skipped += 1;
                                        debug!(
                                            "record {record_name} has no base \
                                        modification information, skipping"
                                        );
                                        continue;
                                    } else {
                                        self.num_used += 1;
                                        ret = Some((
                                            record,
                                            record_name,
                                            modbase_info,
                                        ));
                                        break;
                                    }
                                }
                                Err(e) => match e {
                                    RunError::BadInput(e) => {
                                        debug!(
                                    "record {record_name} has improper data, {}",
                                    e.to_string()
                                );
                                        self.num_failed += 1;
                                        continue;
                                    }
                                    RunError::Failed(e) => {
                                        debug!(
                                    "record {record_name} failed to extract \
                                    mod base info, {}",
                                    e.to_string()
                                );
                                        self.num_failed += 1;
                                        continue;
                                    }
                                    RunError::Skipped(reason) => {
                                        debug!(
                                            "record {record_name} skipped, {}",
                                            reason.to_string()
                                        );
                                        self.num_skipped += 1;
                                        continue;
                                    }
                                },
                            }
                        }
                    }
                }
                Err(e) => {
                    debug!(
                        "failed to read record from bam information, {}",
                        e.to_string()
                    );
                }
            }
        }
        ret
    }
}

pub(crate) fn filter_records_iter<T: bam::Read>(
    records: bam::Records<T>,
) -> impl Iterator<Item = (bam::Record, ModBaseInfo)> + '_ {
    records
        // skip records that fail to parse htslib (todo this could be cleaned up)
        .filter_map(|res| res.ok())
        // skip non-primary
        .filter(|record| !record_is_secondary(&record))
        // skip records with empty sequences
        .filter(|record| record.seq_len() > 0)
        .filter_map(|record| {
            ModBaseInfo::new_from_record(&record).ok().and_then(
                |mod_base_info| {
                    if mod_base_info.is_empty() {
                        None
                    } else {
                        Some((record, mod_base_info))
                    }
                },
            )
        })
}

pub const MM_TAGS: [&str; 2] = ["MM", "Mm"];
pub const ML_TAGS: [&str; 2] = ["ML", "Ml"];

pub type RawModCode = char;

pub struct RawModTags {
    raw_mm: String,
    raw_ml: Vec<u16>,
    mm_style: &'static str,
    ml_style: &'static str,
}

impl RawModTags {
    #[cfg(test)]
    fn new(raw_mm: &str, raw_ml: &[u16], new_style: bool) -> Self {
        let mm_style = if new_style { MM_TAGS[0] } else { MM_TAGS[1] };
        let ml_style = if new_style { ML_TAGS[0] } else { ML_TAGS[1] };
        Self {
            raw_mm: raw_mm.to_owned(),
            raw_ml: raw_ml.to_vec(),
            mm_style,
            ml_style,
        }
    }

    pub fn get_raw_mm(&self) -> &str {
        &self.raw_mm
    }

    pub fn get_raw_ml(&self) -> &[u16] {
        &self.raw_ml
    }

    pub fn mm_is_new_style(&self) -> bool {
        self.mm_style == MM_TAGS[0]
    }

    pub fn ml_is_new_style(&self) -> bool {
        self.ml_style == ML_TAGS[0]
    }
}

#[derive(Debug, Clone)]
pub enum CollapseMethod {
    /// ModCode is the modified base to remove
    ReNormalize(RawModCode),
    /// ModCode is the modified base to remove
    ReDistribute(RawModCode),
    /// Convert one mod base to another
    Convert {
        from: HashSet<RawModCode>,
        to: RawModCode,
    },
}

impl CollapseMethod {
    pub fn parse_str(
        raw: &str,
        mod_code: RawModCode,
    ) -> Result<Self, InputError> {
        match raw {
            "norm" => Ok(Self::ReNormalize(mod_code)),
            "dist" => Ok(Self::ReDistribute(mod_code)),
            _ => Err(InputError::new(&format!("bad collapse method: {}", raw))),
        }
    }
}

#[derive(Debug, Copy, Clone, Eq, PartialEq, Hash)]
pub enum SkipMode {
    /// '?' mode, no probability for a position means we have no information
    /// about base modifications at that position
    Ambiguous,
    /// '.' mode, no probability means the base is canonical (or predicted
    /// canonical).
    ProbModified,
    /// Same as `ProbModified` except the BAM record does not specify the
    /// actual mode.
    ImplicitProbModified,
}

impl SkipMode {
    fn parse(raw_mode: char) -> Result<Self, InputError> {
        match raw_mode {
            '?' => Ok(Self::Ambiguous),
            '.' => Ok(Self::ProbModified),
            _ => Err(InputError::new(&format!("unknown mode {}", raw_mode))),
        }
    }

    fn char(&self) -> Option<char> {
        match self {
            Self::Ambiguous => Some('?'),
            Self::ProbModified => Some('.'),
            Self::ImplicitProbModified => None,
        }
    }
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub enum BaseModCall {
    Canonical(f32),
    Modified(f32, ModCode),
    Filtered,
}

impl Eq for BaseModCall {}

impl Ord for BaseModCall {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other)
            .expect("should not have NaN probability")
    }
}

impl PartialOrd for BaseModCall {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        let p_a = match self {
            BaseModCall::Canonical(p) | BaseModCall::Modified(p, _) => Some(*p),
            BaseModCall::Filtered => None,
        };
        let p_b = match other {
            BaseModCall::Canonical(p) | BaseModCall::Modified(p, _) => Some(*p),
            BaseModCall::Filtered => None,
        };
        p_a.partial_cmp(&p_b)
    }
}

#[derive(Debug, PartialEq, Clone)]
pub struct BaseModProbs {
    probs: FxHashMap<char, f32>,
    // skip_mode: SkipMode,
    // strand: Strand,
}

impl BaseModProbs {
    pub fn new(mod_code: char, prob: f32) -> Self {
        Self {
            probs: [(mod_code, prob)].into_iter().collect(),
        }
    }

    pub fn insert_base_mod_prob(&mut self, mod_code: char, prob: f32) {
        (*self.probs.entry(mod_code).or_insert(0f32)) += prob;
    }

    pub fn argmax_base_mod_call(&self) -> anyhow::Result<BaseModCall> {
        let canonical_prob = self.canonical_prob();
        let max_mod_prob = self
            .iter_probs()
            .max_by(|(_, p), (_, q)| p.partial_cmp(q).unwrap());
        let base_mod_call = if let Some((mod_code, mod_prob)) = max_mod_prob {
            let mod_code = ModCode::parse_raw_mod_code(*mod_code)?;
            if *mod_prob > canonical_prob {
                BaseModCall::Modified(*mod_prob, mod_code)
            } else {
                BaseModCall::Canonical(canonical_prob)
            }
        } else {
            BaseModCall::Canonical(canonical_prob)
        };
        Ok(base_mod_call)
    }

    pub fn canonical_prob(&self) -> f32 {
        1f32 - self.probs.values().sum::<f32>()
    }

    // todo(arand): these methods should be removed/renamed to be more useful
    pub fn iter_probs(&self) -> impl Iterator<Item = (&char, &f32)> {
        self.probs.iter()
    }

    pub fn iter_mut_probs(&mut self) -> impl Iterator<Item = &mut f32> {
        self.probs.iter_mut().map(|(_, p)| p)
    }

    pub fn iter_mut(&mut self) -> impl Iterator<Item = (&char, &mut f32)> {
        self.probs.iter_mut()
    }

    pub(crate) fn into_collapsed(
        self,
        method: &CollapseMethod,
    ) -> BaseModProbs {
        let canonical_prob = self.canonical_prob();
        match method {
            CollapseMethod::ReNormalize(mod_to_collapse) => {
                let marginal_collapsed_prob = self
                    .iter_probs()
                    .filter(|(mod_code, _prob)| *mod_code != mod_to_collapse)
                    .collect::<Vec<(&char, &f32)>>();
                let total_marginal_collapsed_prob = marginal_collapsed_prob
                    .iter()
                    .map(|(_, p)| *p)
                    .sum::<f32>()
                    + canonical_prob;

                let probs = marginal_collapsed_prob
                    .into_iter()
                    .map(|(&mod_code, &mod_prob)| {
                        let collapsed_prob =
                            mod_prob / total_marginal_collapsed_prob;
                        (mod_code, collapsed_prob)
                    })
                    .collect();
                Self { probs }
            }
            CollapseMethod::ReDistribute(mod_to_collapse) => {
                let marginal_prob = self
                    .iter_probs()
                    .filter_map(|(mod_code, prob)| {
                        if mod_code == mod_to_collapse {
                            Some(*prob)
                        } else {
                            None
                        }
                    })
                    .sum::<f32>();
                let other_mods = self
                    .iter_probs()
                    .filter(|(mod_code, _prob)| *mod_code != mod_to_collapse)
                    .collect::<Vec<(&char, &f32)>>();

                let n_other_mods = other_mods.len() as f32 + 1f32; // plus 1 for the canonical base
                let prob_to_redistribute = marginal_prob / n_other_mods;

                let mut check_total = 0f32;
                let (probs, check_total) = other_mods
                    .into_iter()
                    .map(|(&mod_code, prob)| {
                        let new_prob = prob + prob_to_redistribute;
                        check_total += new_prob;
                        (mod_code, new_prob)
                    })
                    .fold(
                        (FxHashMap::default(), 0f32),
                        |(mut probs, total), (mod_code, new_prob)| {
                            probs.insert(mod_code, new_prob);
                            (probs, total + new_prob)
                        },
                    );
                if check_total - 100f32 > 0.00001 {
                    debug!(
                        "total probability {check_total} did not re-normalize"
                    )
                }

                Self { probs }
            }
            CollapseMethod::Convert { from, to } => {
                let (probs, converted_prob) = self.iter_probs().fold(
                    (FxHashMap::default(), 0f32),
                    |(mut probs, converted_prob), (raw_mod_code, prob)| {
                        if from.contains(&raw_mod_code) {
                            (probs, converted_prob + prob)
                        } else {
                            probs.insert(*raw_mod_code, *prob);
                            (probs, converted_prob)
                        }
                    },
                );
                let mut new_base_mod_probs = Self { probs };

                if converted_prob > 0f32 {
                    new_base_mod_probs
                        .insert_base_mod_prob(*to, converted_prob);
                }

                new_base_mod_probs
            }
        }
    }

    fn combine(&mut self, other: Self) {
        for (mod_code, prob) in other.iter_probs() {
            (*self.probs.entry(*mod_code).or_insert(0f32)) += *prob;
        }
    }
}

pub struct DeltaListConverter {
    cumulative_counts: Vec<u32>,
    pub(crate) canonical_base: char,
}

impl DeltaListConverter {
    pub fn new_from_record(
        record: &bam::Record,
        canonical_base: char,
    ) -> Result<Self, RunError> {
        let seq = util::get_forward_sequence(&record)?;

        Ok(Self::new(&seq, canonical_base))
    }

    pub fn new(read_sequence: &str, base: char) -> Self {
        let cumulative_counts = read_sequence
            .chars()
            .scan(0, |count, nt| {
                if nt == base {
                    *count = *count + 1;
                }
                Some(*count)
            })
            .collect::<Vec<u32>>();

        debug_assert_eq!(cumulative_counts.len(), read_sequence.len());
        Self {
            cumulative_counts,
            canonical_base: base,
        }
    }

    pub fn to_positions(
        &self,
        delta_list: &[u32],
    ) -> Result<Vec<usize>, InputError> {
        let mut finger = 0usize;
        let mut n_skips = 0u32;
        let mut positions = Vec::with_capacity(delta_list.len());
        for d in delta_list {
            if finger >= self.cumulative_counts.len() {
                return Err(InputError::new("malformed MM delta list"));
            }
            assert!(finger < self.cumulative_counts.len());
            while self.cumulative_counts[finger] <= (*d + n_skips) {
                finger += 1;
                debug_assert!(
                    finger < self.cumulative_counts.len(),
                    "{:?} >= {:?},\ndelta_list: {:?}\ncumulative counts: {:?}",
                    finger,
                    self.cumulative_counts.len(),
                    delta_list,
                    self.cumulative_counts
                );
                if finger >= self.cumulative_counts.len() {
                    return Err(InputError::new("malformed MM delta list"));
                }
            }
            positions.push(finger);
            n_skips += d + 1;
        }
        Ok(positions)
    }

    pub fn to_delta_list(&self, positions: &[usize]) -> Vec<u32> {
        let mut last = 0;
        let mut delta_list = Vec::new();
        for pos in positions {
            let cumulative_count = self.cumulative_counts[*pos];
            let d = cumulative_count - last - 1;
            delta_list.push(d);
            last = cumulative_count
        }
        delta_list
    }
}

#[inline]
fn prob_to_qual(prob: f32) -> u8 {
    if prob == 1.0f32 {
        255u8
    } else {
        let p = prob * 256f32;
        let q = p.floor() as u8;
        q
    }
}

fn quals_to_probs(quals: &mut [f32]) {
    let arch = pulp::Arch::new();
    arch.dispatch(|| {
        for q in quals {
            let qual = *q;
            *q = (qual + 0.5f32) / 256f32;
        }
    });
}

/// Container for the information in the MM and ML tags
#[derive(Debug, Eq, PartialEq)]
pub struct BaseModPositions {
    pub(crate) canonical_base: char,
    mode: SkipMode,
    strand: Strand,
    mod_base_codes: Vec<char>,
    delta_list: Vec<u32>,
}

fn parse_int_list<'a>(input: &'a str) -> IResult<&str, Vec<u32>> {
    separated_list1(tag(","), |input: &'a str| {
        let (input, _) = multispace0(input)?;
        let (input, num) = digit1(input)?;
        let num = num.parse::<u32>().unwrap();
        let (input, _) = multispace0(input)?;
        Ok((input, num))
    })(input)
}

impl BaseModPositions {
    pub fn parse(mod_positions: &str) -> Result<Self, InputError> {
        let mut parts = mod_positions.split(',');
        let mut header = parts
            .nth(0)
            .ok_or(InputError::new(
                "failed to get leader for base mod position line",
            ))?
            .chars();

        let canonical_base = header
            .nth(0)
            .ok_or(InputError::new("failed to get canonical base"))?;

        let raw_stand = header
            .nth(0)
            .ok_or(InputError::new("failed to get strand"))?;

        let strand = Strand::parse_char(raw_stand)?;

        let mut mod_base_codes = Vec::new();
        let mut mode: Option<SkipMode> = None;

        let mut offset = 2usize;
        while let Some(c) = header.next() {
            match c {
                '?' | '.' => {
                    mode = Some(SkipMode::parse(c).unwrap());
                    offset += 1;
                }
                _ => {
                    mod_base_codes.push(c);
                    offset += 1;
                }
            }
        }
        // default to the "old version"
        let mode = mode.unwrap_or(SkipMode::ImplicitProbModified);

        let delta_list = if offset + 1 <= mod_positions.len() {
            let (_, raw_delta_list) = mod_positions.split_at(offset + 1);
            let (_, delta_list) =
                parse_int_list(raw_delta_list).map_err(|e| {
                    InputError::from(format!(
                        "invalid MM delta list, {}",
                        e.to_string()
                    ))
                })?;
            delta_list
        } else {
            vec![]
        };

        Ok(Self {
            canonical_base,
            mod_base_codes,
            mode,
            strand,
            delta_list,
        })
    }

    fn stride(&self) -> usize {
        self.mod_base_codes.len()
    }

    fn size(&self) -> usize {
        self.delta_list.len() * self.mod_base_codes.len()
    }

    fn is_positive_strand(&self) -> bool {
        self.strand == Strand::Positive
    }
}

fn combine_positions_to_probs(
    record: &bam::Record,
    agg: &mut SeqPosBaseModProbs,
    to_add: SeqPosBaseModProbs,
) -> Result<(), InputError> {
    if agg.skip_mode != to_add.skip_mode {
        let record_name =
            util::get_query_name_string(record).unwrap_or("???".to_string());
        Err(InputError::new(&format!(
            "record: {record_name}, two skip modes ({} and {}) do not match",
            agg.skip_mode.char().unwrap_or('.'),
            to_add.skip_mode.char().unwrap_or('.')
        )))
    } else {
        for (position, base_mod_probs) in
            to_add.pos_to_base_mod_probs.into_iter()
        {
            if let Some(probs) = agg.pos_to_base_mod_probs.get_mut(&position) {
                probs.combine(base_mod_probs);
            } else {
                agg.pos_to_base_mod_probs.insert(position, base_mod_probs);
            }
        }

        Ok(())
    }
}

// pub type SeqPosBaseModProbs = HashMap<usize, BaseModProbs>;
/// Mapping of _forward sequence_ position to `BaseModProbs`.
#[derive(PartialEq, Debug, Clone)]
pub struct SeqPosBaseModProbs {
    /// The `.` or `?` or implied mode, see `SkipMode`.
    pub skip_mode: SkipMode,
    /// Mapping of _forward_ sequence position to the predicted base
    /// modification probabilities for that position.
    pub pos_to_base_mod_probs: FxHashMap<usize, BaseModProbs>,
}

impl SeqPosBaseModProbs {
    // todo(arand) derive new?
    pub(crate) fn new(
        pos_to_base_mod_probs: FxHashMap<usize, BaseModProbs>,
        skip_mode: SkipMode,
    ) -> Self {
        Self {
            skip_mode,
            pos_to_base_mod_probs,
        }
    }
    fn new_empty(skip_mode: SkipMode) -> Self {
        Self {
            skip_mode,
            pos_to_base_mod_probs: FxHashMap::default(),
        }
    }

    /// Remove positions that are outside the bounds of the `EdgeFilter`.
    /// Returning None means that either the read was too short or all
    /// of the positions were filtered out.
    pub(crate) fn edge_filter_positions(
        self,
        edge_filter: &EdgeFilter,
        read_length: usize,
    ) -> Option<Self> {
        if read_length <= edge_filter.edge_filter_start {
            None
        } else {
            read_length
                .checked_sub(edge_filter.edge_filter_end)
                .and_then(|edge_filter_end| {
                    let pos_to_base_mod_probs = self
                        .pos_to_base_mod_probs
                        .into_iter()
                        .filter(|(pos, _)| {
                            *pos >= edge_filter.edge_filter_start
                                && *pos < edge_filter_end
                        })
                        .collect::<FxHashMap<usize, BaseModProbs>>();
                    if pos_to_base_mod_probs.is_empty() {
                        // all positions filtered out
                        None
                    } else {
                        Some(Self::new(pos_to_base_mod_probs, self.skip_mode))
                    }
                })
        }
    }

    pub(crate) fn filter_positions(
        self,
        edge_filter: Option<&EdgeFilter>,
        position_filter: Option<&StrandedPositionFilter>,
        only_mapped: bool,
        aligned_pairs: &FxHashMap<usize, u64>,
        mod_strand: Strand,
        record: &bam::Record,
    ) -> Option<Self> {
        let read_length = record.seq_len();
        let (edge_filter_start, edge_filter_end) =
            if let Some(edge_filter) = edge_filter {
                if read_length <= edge_filter.edge_filter_start {
                    return None;
                }
                match read_length.checked_sub(edge_filter.edge_filter_end) {
                    None => return None,
                    Some(l) => (edge_filter.edge_filter_start, l),
                }
            } else {
                (0, record.seq_len())
            };

        let probs = self
            .pos_to_base_mod_probs
            .into_iter()
            .filter(|(q_pos, _)| {
                let edge_keep =
                    *q_pos >= edge_filter_start && *q_pos < edge_filter_end;
                let only_mapped_keep = if only_mapped {
                    aligned_pairs.contains_key(q_pos)
                } else {
                    true
                };
                let position_keep = match position_filter {
                    Some(position_filter) => aligned_pairs
                        .get(q_pos)
                        .map(|ref_pos| {
                            let reference_strand =
                                match (mod_strand, record.is_reverse()) {
                                    (Strand::Positive, false) => {
                                        Strand::Positive
                                    }
                                    (Strand::Positive, true) => {
                                        Strand::Negative
                                    }
                                    (Strand::Negative, false) => {
                                        Strand::Negative
                                    }
                                    (Strand::Negative, true) => {
                                        Strand::Positive
                                    }
                                };

                            position_filter.contains(
                                record.tid(),
                                *ref_pos,
                                reference_strand,
                            )
                        })
                        .unwrap_or(false),
                    None => true,
                };

                edge_keep && only_mapped_keep && position_keep
            })
            .collect::<FxHashMap<usize, BaseModProbs>>();
        if probs.is_empty() {
            None
        } else {
            Some(Self::new(probs, self.skip_mode))
        }
    }
}

// todo(arand) remove, or put behind cfg(test)
pub fn extract_mod_probs(
    record: &bam::Record,
    raw_mm: &str,
    mod_quals: &[u16],
    converter: &DeltaListConverter,
) -> Result<SeqPosBaseModProbs, InputError> {
    // warn!("[deprecation warning] this method should not be called in production code");
    let splited = raw_mm.split(";");
    let mut positions_to_probs =
        SeqPosBaseModProbs::new_empty(SkipMode::Ambiguous);
    let mut pointer = 0usize;
    for mod_positions in splited {
        if mod_positions.len() == 0 {
            continue;
        }
        let base_mod_positions = BaseModPositions::parse(mod_positions)?;
        if base_mod_positions.canonical_base == converter.canonical_base {
            let base_mod_probs = get_base_mod_probs(
                &base_mod_positions,
                &mod_quals,
                pointer,
                converter,
            )
            .unwrap(); // todo(arand) remove this unwrap
            combine_positions_to_probs(
                record,
                &mut positions_to_probs,
                base_mod_probs,
            )?;
        }
        pointer +=
            base_mod_positions.delta_list.len() * base_mod_positions.stride();
    }

    Ok(positions_to_probs)
}

fn get_base_mod_probs(
    base_mod_positions: &BaseModPositions,
    mod_quals: &[u16],
    pointer: usize,
    converter: &DeltaListConverter,
) -> Result<SeqPosBaseModProbs, InputError> {
    let positions = converter.to_positions(&base_mod_positions.delta_list)?;
    let probs = {
        let mut probs = mod_quals[pointer..pointer + base_mod_positions.size()]
            .iter()
            .map(|qual| *qual as f32)
            .collect::<Vec<f32>>();
        quals_to_probs(&mut probs);
        probs
    };

    let mut positions_to_probs = FxHashMap::<usize, BaseModProbs>::default();
    let stride = base_mod_positions.stride();
    debug_assert_eq!(probs.len() / stride, positions.len());
    for (chunk, position) in probs.chunks(stride).zip(positions) {
        assert_eq!(chunk.len(), stride);
        for (i, mod_base_code) in
            base_mod_positions.mod_base_codes.iter().enumerate()
        {
            let prob = chunk[i];
            if let Some(base_mod_probs) = positions_to_probs.get_mut(&position)
            {
                base_mod_probs.insert_base_mod_prob(*mod_base_code, prob);
            } else {
                positions_to_probs
                    .insert(position, BaseModProbs::new(*mod_base_code, prob));
            }
        }
    }

    Ok(SeqPosBaseModProbs::new(
        positions_to_probs,
        base_mod_positions.mode,
    ))
}

pub fn collapse_mod_probs(
    positions_to_probs: SeqPosBaseModProbs,
    method: &CollapseMethod,
) -> SeqPosBaseModProbs {
    let collapsed_positions_to_probs = positions_to_probs
        .pos_to_base_mod_probs
        .into_iter()
        .map(|(pos, mod_base_probs)| {
            (pos, mod_base_probs.into_collapsed(method))
        })
        .collect();
    SeqPosBaseModProbs {
        pos_to_base_mod_probs: collapsed_positions_to_probs,
        skip_mode: positions_to_probs.skip_mode,
    }
}

pub fn format_mm_ml_tag(
    positions_to_probs: SeqPosBaseModProbs,
    strand: Strand,
    converter: &DeltaListConverter,
) -> (String, Vec<u8>) {
    let canonical_base = converter.canonical_base;
    let mut mod_code_to_position =
        HashMap::<(char, Strand), Vec<(usize, f32)>>::new();

    for (position, mod_base_probs) in positions_to_probs.pos_to_base_mod_probs {
        for (mod_base_code, mod_base_prob) in mod_base_probs.iter_probs() {
            let entry = mod_code_to_position
                .entry((*mod_base_code, strand))
                .or_insert(Vec::new());
            entry.push((position, *mod_base_prob));
        }
    }

    let mut mm_tag = String::new();
    let mut ml_tag = Vec::new();
    let skip_mode_label = positions_to_probs
        .skip_mode
        .char()
        .map(|s| s.to_string())
        .unwrap_or("".to_string());
    if mod_code_to_position.is_empty() {
        let raw_mod_code = ModCode::get_ambig_modcode(canonical_base)
            .map(|mod_code| mod_code.char())
            .unwrap_or(canonical_base);

        mm_tag.push_str(&format!(
            "{}{}{}{};",
            canonical_base,
            strand.to_char(),
            raw_mod_code,
            skip_mode_label
        ));
    } else {
        // todo(arand) this should emit C+hm style tags when possible
        for ((mod_code, strand), mut positions_and_probs) in
            mod_code_to_position.into_iter().sorted_by(
                |((mc_a, s_a), _), ((mc_b, s_b), _)| match mc_a.cmp(mc_b) {
                    Ordering::Equal => s_a.cmp(&s_b),
                    ordering @ _ => ordering,
                },
            )
        {
            positions_and_probs
                .sort_by(|(x_pos, _), (y_pos, _)| x_pos.cmp(&y_pos));
            let header = format!(
                "{}{}{}{},", // C+m?,
                canonical_base,
                strand.to_char(),
                mod_code,
                skip_mode_label
            );
            let positions = positions_and_probs
                .iter()
                .map(|(pos, _prob)| *pos)
                .collect::<Vec<usize>>();
            let delta_list = converter.to_delta_list(&positions);
            let delta_list = delta_list
                .into_iter()
                .map(|d| d.to_string())
                .collect::<Vec<String>>()
                .join(",");
            mm_tag.push_str(&header);
            mm_tag.push_str(&delta_list);
            mm_tag.push(';');
            let quals = positions_and_probs
                .iter()
                .map(|(_pos, prob)| prob_to_qual(*prob))
                .collect::<Vec<u8>>();
            ml_tag.extend(quals.into_iter());
        }
    }

    (mm_tag, ml_tag)
}

/// tag keys should be the new then old tags, for example ["MM", "Mm"].
fn parse_mm_tag(mm_aux: &Aux, tag_key: &str) -> Result<String, RunError> {
    match mm_aux {
        Aux::String(s) => Ok(s.to_string()),
        _ => Err(RunError::new_input_error(format!(
            "incorrect {} tag, should be string",
            tag_key
        ))),
    }
}

/// tag keys should be the new then old tags, for example ["ML", "Ml"].
fn parse_ml_tag(ml_aux: &Aux, tag_key: &str) -> Result<Vec<u16>, RunError> {
    match ml_aux {
        Aux::ArrayU8(arr) => Ok(arr.iter().map(|x| x as u16).collect()),
        _ => Err(RunError::new_input_error(format!(
            "invalid {} tag, expected array",
            tag_key
        ))),
    }
}

pub fn get_mm_tag_from_record(
    record: &bam::Record,
) -> Option<Result<(String, &'static str), RunError>> {
    get_tag::<String>(&record, &MM_TAGS, &parse_mm_tag)
}

pub fn get_ml_tag_from_record(
    record: &bam::Record,
) -> Option<Result<(Vec<u16>, &'static str), RunError>> {
    get_tag::<Vec<u16>>(&record, &ML_TAGS, &parse_ml_tag)
}

pub fn parse_raw_mod_tags(
    record: &bam::Record,
) -> Option<Result<RawModTags, RunError>> {
    let mm = get_mm_tag_from_record(record);
    let ml = get_ml_tag_from_record(record);
    match (mm, ml) {
        (None, _) | (_, None) => None,
        (Some(Ok((raw_mm, mm_style))), Some(Ok((raw_ml, ml_style)))) => {
            Some(Ok(RawModTags {
                raw_mm,
                raw_ml,
                mm_style,
                ml_style,
            }))
        }
        (Some(Err(err)), _) => Some(Err(RunError::new_input_error(format!(
            "MM tag malformed {}",
            err.to_string()
        )))),
        (_, Some(Err(err))) => Some(Err(RunError::new_input_error(format!(
            "ML tag malformed {}",
            err.to_string()
        )))),
    }
}

pub struct ModBaseInfo {
    pos_seq_base_mod_probs: HashMap<char, SeqPosBaseModProbs>,
    neg_seq_base_mod_probs: HashMap<char, SeqPosBaseModProbs>,
    converters: HashMap<char, DeltaListConverter>,
    pub mm_style: &'static str,
    pub ml_style: &'static str,
}

impl ModBaseInfo {
    pub fn new_from_record(record: &bam::Record) -> Result<Self, RunError> {
        let raw_mod_tags = match parse_raw_mod_tags(record) {
            Some(Ok(raw_mod_tags)) => raw_mod_tags,
            Some(Err(run_error)) => {
                return Err(run_error);
            }
            None => {
                return Err(RunError::new_skipped("no mod tags"));
            }
        };

        let forward_sequence = util::get_forward_sequence(record)?;
        Self::new(&raw_mod_tags, &forward_sequence, record)
    }

    pub fn new(
        raw_mod_tags: &RawModTags,
        forward_seq: &str,
        record: &bam::Record,
    ) -> Result<Self, RunError> {
        let mm = &raw_mod_tags.raw_mm;
        let raw_ml = &raw_mod_tags.raw_ml;

        let mut pos_seq_base_mod_probs = HashMap::new();
        let mut converters = HashMap::new();
        let mut neg_seq_base_mod_probs = HashMap::new();
        let mut pointer = 0usize;
        for raw_mm in mm.split(';').filter(|raw_mm| !raw_mm.is_empty()) {
            let base_mod_positions = BaseModPositions::parse(raw_mm)?;
            let converter = converters
                .entry(base_mod_positions.canonical_base)
                .or_insert(DeltaListConverter::new(
                    forward_seq,
                    base_mod_positions.canonical_base,
                ));
            let base_mod_probs = get_base_mod_probs(
                &base_mod_positions,
                &raw_ml,
                pointer,
                &converter,
            )?;

            let seq_base_mod_probs = if base_mod_positions.is_positive_strand()
            {
                &mut pos_seq_base_mod_probs
            } else {
                &mut neg_seq_base_mod_probs
            };

            if let Some(positions_to_probs) =
                seq_base_mod_probs.get_mut(&base_mod_positions.canonical_base)
            {
                combine_positions_to_probs(
                    record,
                    positions_to_probs,
                    base_mod_probs,
                )?;
            } else {
                seq_base_mod_probs
                    .insert(base_mod_positions.canonical_base, base_mod_probs);
            }

            pointer += base_mod_positions.delta_list.len()
                * base_mod_positions.stride();
        }

        Ok(Self {
            pos_seq_base_mod_probs,
            neg_seq_base_mod_probs,
            converters,
            mm_style: raw_mod_tags.mm_style,
            ml_style: raw_mod_tags.ml_style,
        })
    }

    pub fn into_iter_base_mod_probs(
        self,
    ) -> (
        HashMap<char, DeltaListConverter>,
        impl Iterator<Item = (char, Strand, SeqPosBaseModProbs)>,
    ) {
        // todo(arand) change the Item here to include the converter
        let pos_iter = self.pos_seq_base_mod_probs.into_iter().map(
            |(canonical_base, seq_pos_base_mod_probs)| {
                (canonical_base, Strand::Positive, seq_pos_base_mod_probs)
            },
        );
        let neg_iter = self.neg_seq_base_mod_probs.into_iter().map(
            |(canonical_base, seq_pos_base_mod_probs)| {
                (canonical_base, Strand::Negative, seq_pos_base_mod_probs)
            },
        );
        (self.converters, pos_iter.chain(neg_iter))
    }

    pub(crate) fn is_empty(&self) -> bool {
        let n_probs = self
            .pos_seq_base_mod_probs
            .values()
            .chain(self.neg_seq_base_mod_probs.values())
            .map(|p| p.pos_to_base_mod_probs.len())
            .sum::<usize>();
        n_probs == 0
    }

    pub fn iter_seq_base_mod_probs(
        &self,
    ) -> impl Iterator<Item = (&char, Strand, &SeqPosBaseModProbs)> {
        let pos_iter = self.pos_seq_base_mod_probs.iter().map(
            |(canonical_base, seq_pos_base_mod_probs)| {
                (canonical_base, Strand::Positive, seq_pos_base_mod_probs)
            },
        );
        let neg_iter = self.neg_seq_base_mod_probs.iter().map(
            |(canonical_base, seq_pos_base_mod_probs)| {
                (canonical_base, Strand::Negative, seq_pos_base_mod_probs)
            },
        );
        pos_iter.chain(neg_iter)
    }
}

pub fn get_canonical_bases_with_mod_calls(
    record: &bam::Record,
) -> Result<Vec<DnaBase>, RunError> {
    match parse_raw_mod_tags(record) {
        Some(Ok(raw_mod_tags)) => raw_mod_tags
            .raw_mm
            .split(';')
            .filter_map(|raw_mm| {
                if raw_mm.is_empty() {
                    None
                } else {
                    Some(BaseModPositions::parse(raw_mm).and_then(
                        |base_mod_positions| {
                            DnaBase::parse(base_mod_positions.canonical_base)
                                .map_err(|e| InputError::new(&e.to_string()))
                        },
                    ))
                }
            })
            .collect::<Result<HashSet<DnaBase>, InputError>>()
            .map(|canonical_bases| {
                canonical_bases.into_iter().collect::<Vec<DnaBase>>()
            })
            .map_err(|input_err| input_err.into()),
        Some(Err(e)) => Err(e),
        None => Ok(Vec::new()),
    }
}

pub fn base_mod_probs_from_record(
    record: &bam::Record,
    converter: &DeltaListConverter,
) -> Result<SeqPosBaseModProbs, RunError> {
    let (mm, ml) = match parse_raw_mod_tags(record) {
        Some(Ok(raw_mod_tags)) => (raw_mod_tags.raw_mm, raw_mod_tags.raw_ml),
        Some(Err(run_error)) => {
            return Err(run_error);
        }
        None => {
            return Err(RunError::new_skipped("no mod tags"));
        }
    };

    extract_mod_probs(record, &mm, &ml, &converter)
        .map_err(|input_err| RunError::BadInput(input_err))
}

#[derive(new, Debug)]
pub struct EdgeFilter {
    // TODO(arand) in the CLIs only one option is allowed, i.e.
    // both `edge_filter_start` and `edge_filter_end` will be the
    // same. If this changes make sure there are suitable tests.
    pub(crate) edge_filter_start: usize,
    pub(crate) edge_filter_end: usize,
}

#[cfg(test)]
mod mod_bam_tests {
    use super::*;
    use crate::util::get_aligned_pairs_forward;
    use rust_htslib::bam::Read;
    use rustc_hash::FxHashSet;
    use std::collections::BTreeSet;
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    use std::path::Path;

    fn qual_to_prob(qual: u16) -> f32 {
        let q = qual as f32;
        (q + 0.5f32) / 256f32
    }

    // first implementation that does not account for multiple mods in the MM tag (i.e. C and A))
    pub fn get_mod_probs_for_query_positions(
        mm: &str,
        canonical_base: char,
        mod_quals: &[u16],
        converter: &DeltaListConverter,
    ) -> Result<HashMap<usize, BaseModProbs>, InputError> {
        // todo move this outside this function. should handle the case where mods for another base
        //  come first and the offset of mod_quals is already handled
        let filtered_mod_positions = mm
            .split(';')
            .filter(|positions| positions.starts_with(canonical_base))
            .collect::<Vec<&str>>();

        let mut probs_for_positions = HashMap::<usize, BaseModProbs>::new();
        let mut prob_array_idx = 0usize;
        for mod_positions in filtered_mod_positions {
            let mut parts = mod_positions.split(',');
            let mut header = parts
                .nth(0)
                .ok_or(InputError::new(
                    "failed to get leader for base mod position line",
                ))?
                .chars();

            let raw_stand = header
                .nth(1)
                .ok_or(InputError::new("failed to get strand"))?;

            // TODO handle duplex
            let _strand = Strand::parse_char(raw_stand).unwrap();

            let mut mod_base_codes = Vec::new();
            let mut _mode: Option<char> = None;

            while let Some(c) = header.next() {
                match c {
                    '?' | '.' => {
                        _mode = Some(c);
                    }
                    _ => mod_base_codes.push(c),
                }
            }

            // taking the liberty to think that a read wouldn't be larger
            // than 2**32 - 1 bases long
            let delta_list = parts
                .into_iter()
                .map(|raw_pos| raw_pos.parse::<u32>())
                .collect::<Result<Vec<u32>, _>>()
                .map_err(|e| {
                    InputError::new(&format!(
                        "failed to parse position list, {}",
                        e.to_string()
                    ))
                })?;

            let positions = converter.to_positions(&delta_list)?;
            for mod_base in mod_base_codes {
                for pos in &positions {
                    let qual = mod_quals[prob_array_idx];
                    let prob = qual_to_prob(qual);
                    if let Some(base_mod_probs) =
                        probs_for_positions.get_mut(pos)
                    {
                        base_mod_probs.insert_base_mod_prob(mod_base, prob);
                    } else {
                        probs_for_positions
                            .insert(*pos, BaseModProbs::new(mod_base, prob));
                    }
                    // consume from the ML array
                    prob_array_idx += 1;
                }
            }
        }

        Ok(probs_for_positions)
    }

    #[test]
    fn test_delta_list_to_positions() {
        let canonical_base = 'C';
        let read_sequence = "ACCGCCGTCGTCG";
        let converter = DeltaListConverter::new(read_sequence, canonical_base);

        let ds = [1, 1, 0];
        let expected = [2, 5, 8];
        let obs = converter.to_positions(&ds).unwrap();
        assert_eq!(&obs, &expected);
        let obs = converter.to_delta_list(&expected);
        assert_eq!(obs, &ds);

        let ds = [3, 0, 0];
        let expected = [5, 8, 11];
        let obs = converter.to_positions(&ds).unwrap();
        assert_eq!(&obs, &expected);
        let obs = converter.to_delta_list(&expected);
        assert_eq!(obs, &ds);

        let ds = [3, 1];
        let expected = [5, 11];
        let obs = converter.to_positions(&ds).unwrap();
        assert_eq!(&obs, &expected);
        let obs = converter.to_delta_list(&expected);
        assert_eq!(obs, &ds);
    }

    #[test]
    fn test_mod_prob_collapse() {
        let probs = vec![('h', 0.85), ('m', 0.10)].into_iter().collect();

        let mod_base_probs = BaseModProbs { probs };
        let collapsed = mod_base_probs
            .clone()
            .into_collapsed(&CollapseMethod::ReDistribute('h'));
        assert_eq!(
            collapsed.probs,
            vec![('m', 0.52500004)]
                .into_iter()
                .collect::<FxHashMap<char, f32>>()
        );
        let collapsed = mod_base_probs
            .clone()
            .into_collapsed(&CollapseMethod::ReNormalize('h'));
        assert_eq!(
            collapsed.probs,
            vec![('m', 0.6666669)]
                .into_iter()
                .collect::<FxHashMap<char, f32>>()
        );

        let collapsed = mod_base_probs
            .clone()
            .into_collapsed(&CollapseMethod::ReNormalize('a'));
        assert_eq!(&collapsed, &mod_base_probs);
        let collapsed = mod_base_probs
            .clone()
            .into_collapsed(&CollapseMethod::ReDistribute('a'));
        assert_eq!(&collapsed, &mod_base_probs);
    }

    #[test]
    fn test_mod_prob_collapse_norm_examples() {
        let probs = vec![('h', 0.05273438), ('m', 0.03320312)]
            .into_iter()
            .collect();

        let mod_base_probs = BaseModProbs { probs };
        let collapsed =
            mod_base_probs.into_collapsed(&CollapseMethod::ReNormalize('h'));
        assert_eq!(
            collapsed.probs,
            vec![('m', 0.035051543)]
                .into_iter()
                .collect::<FxHashMap<char, f32>>()
        );
    }

    #[test]
    fn test_mod_prob_collapse_dist_examples() {
        let probs = vec![('h', 0.05273438), ('m', 0.03320312)]
            .into_iter()
            .collect();
        let mod_base_probs = BaseModProbs { probs };
        let collapsed =
            mod_base_probs.into_collapsed(&CollapseMethod::ReDistribute('h'));
        assert_eq!(
            collapsed.probs,
            vec![('m', 0.059570313)]
                .into_iter()
                .collect::<FxHashMap<char, f32>>()
        );
    }

    #[test]
    fn test_mod_prob_convert() {
        let probs = vec![('h', 0.10), ('m', 0.75)]
            .into_iter()
            .collect::<FxHashMap<char, f32>>();
        let mod_base_probs = BaseModProbs {
            probs: probs.clone(),
        };

        let collapsed =
            mod_base_probs.into_collapsed(&CollapseMethod::Convert {
                from: HashSet::from(['h']),
                to: 'C',
            });
        assert_eq!(
            collapsed.probs,
            vec![('m', 0.75), ('C', 0.10)]
                .into_iter()
                .collect::<FxHashMap<char, f32>>()
        );
        let mod_base_probs = BaseModProbs {
            probs: probs.clone(),
        };
        let collapsed =
            mod_base_probs.into_collapsed(&CollapseMethod::Convert {
                from: HashSet::from(['h', 'm']),
                to: 'C',
            });
        assert_eq!(
            collapsed.probs,
            vec![('C', 0.85)]
                .into_iter()
                .collect::<FxHashMap<char, f32>>()
        );
    }

    #[test]
    fn test_mod_prob_convert_sums_prob() {
        let probs = vec![('h', 0.10), ('m', 0.75)].into_iter().collect();
        let mod_base_probs = BaseModProbs { probs };
        let collapsed =
            mod_base_probs.into_collapsed(&CollapseMethod::Convert {
                from: HashSet::from(['h']),
                to: 'm',
            });
        assert_eq!(
            collapsed.probs,
            vec![('m', 0.85)]
                .into_iter()
                .collect::<FxHashMap<char, f32>>()
        );
    }

    #[test]
    fn test_mod_prob_convert_noop() {
        let probs = vec![('h', 0.10), ('m', 0.75)]
            .into_iter()
            .collect::<FxHashMap<char, f32>>();
        let mod_base_probs = BaseModProbs {
            probs: probs.clone(),
        };
        let collapsed =
            mod_base_probs.into_collapsed(&CollapseMethod::Convert {
                from: HashSet::from(['a']),
                to: 'A',
            });
        assert_eq!(collapsed.probs, probs);
    }

    #[test]
    fn test_mod_prob_combine() {
        let a_probs = vec![('h', 0.05273438), ('m', 0.03320312)]
            .into_iter()
            .collect();
        let mut a = BaseModProbs { probs: a_probs };
        let b_probs = vec![('m', 0.03320312)].into_iter().collect();
        let b = BaseModProbs { probs: b_probs };
        a.combine(b);
        assert_eq!(
            &a.probs,
            &vec![('h', 0.05273438), ('m', 0.06640624)]
                .into_iter()
                .collect::<FxHashMap<char, f32>>()
        );

        let a_probs = vec![('m', 0.03320312)].into_iter().collect();
        let b_probs = vec![('h', 0.05273438)].into_iter().collect();

        let mut a = BaseModProbs { probs: a_probs };

        let b = BaseModProbs { probs: b_probs };
        a.combine(b);
        assert_eq!(
            &a.probs,
            &[('m', 0.03320312), ('h', 0.05273438)]
                .into_iter()
                .collect::<FxHashMap<char, f32>>()
        );
    }

    #[test]
    fn test_parse_mm_tag() {
        let tag =
            "C+h?,5,2,1,3,1,2,3,1,2,1,11,5;C+m?,5,2,1,3,1,2,3,1,2,1,11,5;";
        let dna = "ATGTGCCTGCTGGACATGTTTATGCTCGTCTACTTCGTTCAGTTACGTATTGCTCCAG\
            CGCTCGAACTGTAGCCGCTGCTGCTGGGTGAAGTTGTGGCGGTACACGAGCTCCGCCGGCTGCAGCAGCTTC\
            TCCCCATCCTGGCGCTTCTCCCCGAGCAATTGGTG";
        let mod_quals = vec![
            197, 13, 156, 1, 3, 5, 9, 26, 8, 1, 0, 13, 10, 67, 1, 0, 1, 0, 5,
            5, 5, 0, 0, 8,
        ];

        let converter = DeltaListConverter::new(dna, 'C');
        let positions_to_probs =
            get_mod_probs_for_query_positions(tag, 'C', &mod_quals, &converter)
                .unwrap();
        assert_eq!(positions_to_probs.len(), 12);
    }

    #[test]
    fn test_format_mm_ml_tags() {
        let canonical_base = 'C';
        let read_sequence = "ACCGCCGTCGTCG";
        let converter = DeltaListConverter::new(read_sequence, canonical_base);

        let positions_and_probs = vec![
            (5, BaseModProbs::new('m', 0.9)),
            (2, BaseModProbs::new('m', 0.1)),
            (8, BaseModProbs::new('m', 0.2)),
        ]
        .into_iter()
        .collect::<FxHashMap<usize, BaseModProbs>>();

        let seq_pos_base_mod_probs =
            SeqPosBaseModProbs::new(positions_and_probs, SkipMode::Ambiguous);
        let (mm, ml) = format_mm_ml_tag(
            seq_pos_base_mod_probs,
            Strand::Positive,
            &converter,
        );
        assert_eq!(mm, "C+m?,1,1,0;");
        assert_eq!(ml, vec![25, 230, 51]);

        let skip_mode = SkipMode::ProbModified;
        let positions_and_probs = vec![
            (5, BaseModProbs::new('m', 0.9)),
            (2, BaseModProbs::new('m', 0.1)),
            (8, BaseModProbs::new('m', 0.2)),
        ]
        .into_iter()
        .collect::<FxHashMap<usize, BaseModProbs>>();

        let seq_pos_base_mod_probs =
            SeqPosBaseModProbs::new(positions_and_probs, skip_mode);
        let (mm, ml) = format_mm_ml_tag(
            seq_pos_base_mod_probs,
            Strand::Positive,
            &converter,
        );
        assert_eq!(mm, "C+m.,1,1,0;");
        assert_eq!(ml, vec![25, 230, 51]);
    }

    #[test]
    fn test_mod_parse_base_positions() {
        let raw_positions = "C+h?,5,2,1,3,1,2,3,1,2,1,11,5;";
        let base_mod_positions =
            BaseModPositions::parse(raw_positions).unwrap();
        let expected = BaseModPositions {
            canonical_base: 'C',
            mode: SkipMode::Ambiguous,
            strand: Strand::Positive,
            mod_base_codes: vec!['h'],
            delta_list: vec![5, 2, 1, 3, 1, 2, 3, 1, 2, 1, 11, 5],
        };

        assert_eq!(base_mod_positions, expected);

        let raw_positions = "C+m,5,2,1,3,1,2,3,1,2,1,11,5;";
        let base_mod_positions =
            BaseModPositions::parse(raw_positions).unwrap();
        let expected = BaseModPositions {
            canonical_base: 'C',
            mode: SkipMode::ImplicitProbModified,
            strand: Strand::Positive,
            mod_base_codes: vec!['m'],
            delta_list: vec![5, 2, 1, 3, 1, 2, 3, 1, 2, 1, 11, 5],
        };

        assert_eq!(base_mod_positions, expected);
        let raw_positions = "C+m.,5,2,1,3,1,2,3,1,2,1,11,5;";
        let base_mod_positions =
            BaseModPositions::parse(raw_positions).unwrap();
        let expected = BaseModPositions {
            canonical_base: 'C',
            mode: SkipMode::ProbModified,
            strand: Strand::Positive,
            mod_base_codes: vec!['m'],
            delta_list: vec![5, 2, 1, 3, 1, 2, 3, 1, 2, 1, 11, 5],
        };
        assert_eq!(base_mod_positions, expected);
    }

    #[test]
    fn test_get_base_mod_probs() {
        let dna = "GATCGACTACGTCGA";
        let tag = "C+hm?,0,1,0;";
        let quals = vec![1, 200, 1, 200, 1, 200];
        let canonical_base = 'C';
        let converter = DeltaListConverter::new(dna, canonical_base);

        let positions_to_probs =
            extract_mod_probs(&bam::Record::new(), tag, &quals, &converter)
                .unwrap();

        assert_eq!(positions_to_probs.pos_to_base_mod_probs.len(), 3);
        let mut found_positions = Vec::new();
        for (position, base_mod_probs) in
            positions_to_probs.pos_to_base_mod_probs.iter()
        {
            found_positions.push(position);
            let quals = base_mod_probs
                .probs
                .iter()
                .map(|(_, p)| prob_to_qual(*p))
                .collect::<Vec<_>>();
            assert_eq!(&quals, &[1, 200]);
        }

        let tag = "C+h?,0,1,0;C+m?,0,1,0;";
        let quals = vec![1, 1, 1, 200, 200, 200];

        let positions_to_probs_1 =
            extract_mod_probs(&bam::Record::new(), tag, &quals, &converter)
                .unwrap();

        assert_eq!(positions_to_probs, positions_to_probs_1);
    }

    #[test]
    fn test_extract_positions_to_probs() {
        let dna = "GATCGACTACGTCGA";
        let tag = "C+h?,0,1,0;A+a?,0,1,0;C+m?,0,1,0;";
        let quals = vec![1, 1, 1, 200, 200, 200, 1, 1, 1];
        let canonical_base = 'C';
        let converter = DeltaListConverter::new(dna, canonical_base);

        let positions_to_probs =
            extract_mod_probs(&bam::Record::new(), tag, &quals, &converter)
                .unwrap();
        assert_eq!(positions_to_probs.pos_to_base_mod_probs.len(), 3);
        for (_pos, base_mod_probs) in
            positions_to_probs.pos_to_base_mod_probs.iter()
        {
            assert_eq!(
                &base_mod_probs.probs,
                &[('h', 0.005859375), ('m', 0.005859375)]
                    .into_iter()
                    .collect::<FxHashMap<char, f32>>()
            );
        }

        let tag = "C+hm?,0,1,0;A+a?,0,1,0;";
        let quals = vec![1, 1, 1, 1, 1, 1, 200, 200, 200];
        let positions_to_probs_comb =
            extract_mod_probs(&bam::Record::new(), tag, &quals, &converter)
                .unwrap();
        assert_eq!(
            positions_to_probs_comb.pos_to_base_mod_probs.len(),
            positions_to_probs.pos_to_base_mod_probs.len()
        );
        for (position, base_mod_probs) in
            positions_to_probs_comb.pos_to_base_mod_probs.iter()
        {
            let other = positions_to_probs
                .pos_to_base_mod_probs
                .get(position)
                .unwrap();
            assert_eq!(base_mod_probs.probs, other.probs);
        }
    }

    #[test]
    fn test_mod_base_info() {
        // Preamble, make a short DNA and the converters
        let dna = "GATCGACTACGTCGA";
        let c_converter = DeltaListConverter::new(dna, 'C');
        let a_converter = DeltaListConverter::new(dna, 'A');

        // these tags only have 1 canonical base, parse these and these are the
        // expected values for the rest of the test
        let c_tag = "C+hm?,0,1,0;";
        let a_tag = "A+a?,0,1,0;";
        let c_quals = vec![1, 100, 1, 100, 1, 100]; // interleaved!
        let a_quals = vec![200, 200, 200];

        let c_expected_seq_pos_base_mod_probs = extract_mod_probs(
            &bam::Record::new(),
            c_tag,
            &c_quals,
            &c_converter,
        )
        .unwrap();
        let a_expected_seq_pos_base_mod_probs = extract_mod_probs(
            &bam::Record::new(),
            a_tag,
            &a_quals,
            &a_converter,
        )
        .unwrap();

        let c_expected_probs = vec![('h', 0.005859375), ('m', 0.39257813)]
            .into_iter()
            .collect();
        let c_expected = BaseModProbs {
            probs: c_expected_probs,
        };
        let a_expected_probs = vec![('a', 0.7832031)].into_iter().collect();
        let a_expected = BaseModProbs {
            probs: a_expected_probs,
        };
        assert_eq!(
            c_expected_seq_pos_base_mod_probs
                .pos_to_base_mod_probs
                .keys()
                .map(|c| *c)
                .collect::<HashSet<usize>>(),
            HashSet::from([9, 12, 3])
        );
        assert_eq!(
            a_expected_seq_pos_base_mod_probs
                .pos_to_base_mod_probs
                .keys()
                .map(|c| *c)
                .collect::<HashSet<usize>>(),
            HashSet::from([1, 8, 14])
        );

        for base_mod_probs in c_expected_seq_pos_base_mod_probs
            .pos_to_base_mod_probs
            .values()
        {
            assert_eq!(base_mod_probs, &c_expected);
        }

        for base_mod_probs in a_expected_seq_pos_base_mod_probs
            .pos_to_base_mod_probs
            .values()
        {
            assert_eq!(base_mod_probs, &a_expected);
        }

        // test with the tag/quals separated
        let tag = "C+h?,0,1,0;A+a?,0,1,0;C+m?,0,1,0;";
        let quals = vec![1, 1, 1, 200, 200, 200, 100, 100, 100];
        let raw_mod_tags = RawModTags::new(tag, &quals, true);
        let obs_mod_base_info =
            ModBaseInfo::new(&raw_mod_tags, dna, &bam::Record::new()).unwrap();
        assert_eq!(
            obs_mod_base_info.pos_seq_base_mod_probs.get(&'C').unwrap(),
            &c_expected_seq_pos_base_mod_probs
        );
        assert_eq!(
            obs_mod_base_info.pos_seq_base_mod_probs.get(&'A').unwrap(),
            &a_expected_seq_pos_base_mod_probs
        );

        let obs_base_mod_probs =
            extract_mod_probs(&bam::Record::new(), tag, &quals, &c_converter)
                .unwrap();
        assert_eq!(&obs_base_mod_probs, &c_expected_seq_pos_base_mod_probs);
        let obs_base_mod_probs =
            extract_mod_probs(&bam::Record::new(), tag, &quals, &a_converter)
                .unwrap();
        assert_eq!(&obs_base_mod_probs, &a_expected_seq_pos_base_mod_probs);

        let tag = "C+h?,0,1,0;C+m?,0,1,0;A+a?,0,1,0;";
        let quals = vec![1, 1, 1, 100, 100, 100, 200, 200, 200];
        let raw_mod_tags = RawModTags::new(tag, &quals, true);
        let obs_mod_base_info =
            ModBaseInfo::new(&raw_mod_tags, dna, &bam::Record::new()).unwrap();
        assert_eq!(
            obs_mod_base_info.pos_seq_base_mod_probs.get(&'C').unwrap(),
            &c_expected_seq_pos_base_mod_probs
        );
        assert_eq!(
            obs_mod_base_info.pos_seq_base_mod_probs.get(&'A').unwrap(),
            &a_expected_seq_pos_base_mod_probs
        );

        let obs_base_mod_probs =
            extract_mod_probs(&bam::Record::new(), tag, &quals, &c_converter)
                .unwrap();
        assert_eq!(&obs_base_mod_probs, &c_expected_seq_pos_base_mod_probs);
        let obs_base_mod_probs =
            extract_mod_probs(&bam::Record::new(), tag, &quals, &a_converter)
                .unwrap();
        assert_eq!(&obs_base_mod_probs, &a_expected_seq_pos_base_mod_probs);

        // test with the mods "combined"
        let tag = "C+hm?,0,1,0;A+a?,0,1,0;";
        let quals = vec![1, 100, 1, 100, 1, 100, 200, 200, 200];
        let raw_mod_tags = RawModTags::new(tag, &quals, true);
        let obs_mod_base_info =
            ModBaseInfo::new(&raw_mod_tags, dna, &bam::Record::new()).unwrap();
        assert_eq!(
            obs_mod_base_info.pos_seq_base_mod_probs.get(&'C').unwrap(),
            &c_expected_seq_pos_base_mod_probs
        );
        assert_eq!(
            obs_mod_base_info.pos_seq_base_mod_probs.get(&'A').unwrap(),
            &a_expected_seq_pos_base_mod_probs
        );
        let obs_base_mod_probs =
            extract_mod_probs(&bam::Record::new(), tag, &quals, &c_converter)
                .unwrap();
        assert_eq!(&obs_base_mod_probs, &c_expected_seq_pos_base_mod_probs);
        let obs_base_mod_probs =
            extract_mod_probs(&bam::Record::new(), tag, &quals, &a_converter)
                .unwrap();
        assert_eq!(&obs_base_mod_probs, &a_expected_seq_pos_base_mod_probs);
    }

    #[test]
    fn test_duplex_modbase_info() {
        //               g c CG c  gg CG CG
        let dna = "GACTCGACTGGACGTCGA";
        let tag = "C+h?,1,1,0;C+m?,1,1,0;G-h?,1,2,0;G-m?,1,2,0";
        let quals = vec![100, 100, 100, 1, 1, 1, 150, 150, 150, 2, 2, 2];
        let tags = RawModTags::new(tag, &quals, true);
        let info = ModBaseInfo::new(&tags, dna, &bam::Record::new()).unwrap();
        let (_converters, iterator) = info.into_iter_base_mod_probs();
        for (c, strand, probs) in iterator {
            if c == 'C' {
                assert_eq!(strand, Strand::Positive);
                assert_eq!(
                    probs
                        .pos_to_base_mod_probs
                        .keys()
                        .map(|i| *i)
                        .collect::<HashSet<usize>>(),
                    HashSet::from([12, 15, 4])
                );
                let expected_probs =
                    vec![('h', 0.39257813), ('m', 0.005859375)]
                        .into_iter()
                        .collect();
                let expected_modbase_probs = BaseModProbs {
                    probs: expected_probs,
                };
                for mod_probs in probs.pos_to_base_mod_probs.values() {
                    assert_eq!(mod_probs, &expected_modbase_probs);
                }
            } else {
                assert_eq!(c, 'G');
                assert_eq!(strand, Strand::Negative);
                assert_eq!(
                    probs
                        .pos_to_base_mod_probs
                        .keys()
                        .map(|i| *i)
                        .collect::<HashSet<usize>>(),
                    HashSet::from([13, 16, 5])
                );
                let expected_probs = vec![('h', 0.5878906), ('m', 0.009765625)]
                    .into_iter()
                    .collect();
                let expected_modbase_probs = BaseModProbs {
                    probs: expected_probs,
                };
                for mod_probs in probs.pos_to_base_mod_probs.values() {
                    assert_eq!(mod_probs, &expected_modbase_probs);
                }
            }
        }
    }

    #[test]
    fn test_base_modcall_equality() {
        let a = BaseModCall::Canonical(1.0);
        let b = BaseModCall::Canonical(1.0);
        let c = BaseModCall::Modified(0.8, ModCode::a);
        let d = BaseModCall::Modified(0.7, ModCode::a);
        let e = BaseModCall::Filtered;
        assert_eq!(a, b);
        assert!(a > c);
        assert!(c > d);
        assert!(d > e);
    }

    #[test]
    fn test_seq_pos_base_mod_probs_edge_filter() {
        //               012345678901234
        let dna = "GATCGACTACGTCGA";
        let tag = "C+h?,0,1,0;A+a?,0,1,0;C+m?,0,1,0;";
        let quals = vec![1, 1, 1, 200, 200, 200, 100, 100, 100];
        let raw_mod_tags = RawModTags::new(tag, &quals, true);
        let mut obs_mod_base_info =
            ModBaseInfo::new(&raw_mod_tags, dna, &bam::Record::new()).unwrap();
        let c_seq_base_mod_probs = obs_mod_base_info
            .pos_seq_base_mod_probs
            .remove(&'C')
            .unwrap();
        let expected_pos = vec![3, 9, 12];
        let obs_pos = c_seq_base_mod_probs
            .pos_to_base_mod_probs
            .keys()
            .map(|p| *p)
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        assert_eq!(expected_pos, obs_pos);
        let edge_filter = EdgeFilter::new(4, 4);
        let c_seq_base_mod_probs = c_seq_base_mod_probs
            .edge_filter_positions(&edge_filter, dna.len())
            .unwrap();
        let obs_pos = c_seq_base_mod_probs
            .pos_to_base_mod_probs
            .keys()
            .map(|p| *p)
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        assert_eq!(vec![9], obs_pos);

        // trim larger than read
        let edge_filter = EdgeFilter::new(50, 50);
        let mut obs_mod_base_info =
            ModBaseInfo::new(&raw_mod_tags, dna, &bam::Record::new()).unwrap();
        let c_seq_base_mod_probs = obs_mod_base_info
            .pos_seq_base_mod_probs
            .remove(&'C')
            .unwrap();
        let c_seq_base_mod_probs =
            c_seq_base_mod_probs.edge_filter_positions(&edge_filter, dna.len());
        assert!(c_seq_base_mod_probs.is_none());

        // trim with mod call _at_ the position to be trimmed
        let mut obs_mod_base_info =
            ModBaseInfo::new(&raw_mod_tags, dna, &bam::Record::new()).unwrap();
        let c_seq_base_mod_probs = obs_mod_base_info
            .pos_seq_base_mod_probs
            .remove(&'C')
            .unwrap();
        let expected_pos = vec![3, 9, 12];
        let obs_pos = c_seq_base_mod_probs
            .pos_to_base_mod_probs
            .keys()
            .map(|p| *p)
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        assert_eq!(expected_pos, obs_pos);
        let edge_filter = EdgeFilter::new(3, 3);
        let c_seq_base_mod_probs = c_seq_base_mod_probs
            .edge_filter_positions(&edge_filter, dna.len())
            .unwrap();
        let obs_pos = c_seq_base_mod_probs
            .pos_to_base_mod_probs
            .keys()
            .map(|p| *p)
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        assert_eq!(vec![3, 9], obs_pos);

        // todo edge_filter_start larger than read length
        //  edge_filter_end larger than read length
        //  asymmetric edge filter with reverse record
    }

    #[test]
    fn test_seq_pos_base_mod_probs_filter_positions() {
        let mut reader = bam::Reader::from_path(
            "tests/resources/bc_anchored_10_reads.sorted.bam",
        )
        .unwrap();
        let header = reader.header().to_owned();
        let records = reader.records();
        let chrom_to_tid = (0..header.target_count())
            .map(|tid| {
                (
                    String::from_utf8(header.tid2name(tid).to_vec()).unwrap(),
                    tid,
                )
            })
            .collect::<HashMap<String, u32>>();

        let position_bed_fp = "tests/resources/CGI_ladder_3.6kb_ref_CG.bed";
        let position_filter = StrandedPositionFilter::from_bed_file(
            &Path::new(position_bed_fp).to_path_buf(),
            &chrom_to_tid.iter().map(|(k, v)| (k.as_str(), *v)).collect(),
            true,
        )
        .unwrap();

        let mut pos_positions = FxHashSet::default();
        let mut neg_positions = FxHashSet::default();
        for line in BufReader::new(File::open(position_bed_fp).unwrap())
            .lines()
            .map(|r| r.unwrap())
        {
            let parts = line.split_whitespace().collect::<Vec<&str>>();
            assert_eq!(parts.len(), 6);
            assert_eq!(parts[0], "oligo_1512_adapters");
            let pos = parts[1].parse::<u64>().unwrap();
            match parts[5] {
                "+" => assert!(pos_positions.insert(pos)),
                "-" => assert!(neg_positions.insert(pos)),
                _ => panic!("illegal strand in BED"),
            }
        }

        let mod_base_info_iter = filter_records_iter(records);
        for (record, mod_base_info) in mod_base_info_iter {
            let aligned_pairs = get_aligned_pairs_forward(&record)
                .filter_map(|pair| pair.ok())
                .collect::<FxHashMap<usize, u64>>();

            let (_converters, base_mod_probs_iter) =
                mod_base_info.into_iter_base_mod_probs();
            for (_primary_base, mod_strand, seq_pos_mod_base_probs) in
                base_mod_probs_iter
            {
                let seq_pos_mod_base_probs = seq_pos_mod_base_probs
                    .filter_positions(
                        None,
                        Some(&position_filter),
                        true,
                        &aligned_pairs,
                        mod_strand,
                        &record,
                    )
                    .unwrap();
                let positions_to_check = if record.is_reverse() {
                    &neg_positions
                } else {
                    &pos_positions
                };
                for (q_pos, _) in seq_pos_mod_base_probs.pos_to_base_mod_probs {
                    let r_pos = aligned_pairs.get(&q_pos).unwrap();
                    assert!(positions_to_check.contains(r_pos));
                }
            }
        }
    }

    #[test]
    fn test_mod_bam_modbase_info_empty() {
        let dna = "GATCGACTACGTCGA";
        let tag = "C+h?;C+m?;";
        let quals = vec![];
        let raw_mod_tags = RawModTags::new(tag, &quals, true);
        let obs_mod_base_info =
            ModBaseInfo::new(&raw_mod_tags, dna, &bam::Record::new()).unwrap();
        assert!(obs_mod_base_info.is_empty());
    }
}
