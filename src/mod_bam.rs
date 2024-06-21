use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::fmt::{Debug, Display, Formatter};
use std::hash::Hash;

use anyhow::bail;
use derive_new::new;
use itertools::{Itertools, PeekingNext};
use log::{debug, error};
use nom::bytes::complete::tag;
use nom::character::complete::{digit1, multispace0};
use nom::multi::separated_list1;
use nom::IResult;
use rust_htslib::bam;
use rust_htslib::bam::record::Aux;
use rustc_hash::FxHashMap;

use crate::errs::{InputError, RunError};
use crate::mod_base_code::{DnaBase, ModCodeRepr};
use crate::util;
use crate::util::{
    get_query_name_string, get_tag, record_is_not_primary, Strand,
};

pub(crate) struct TrackingModRecordIter<'a, T: bam::Read> {
    records: bam::Records<'a, T>,
    skip_unmapped: bool,
    allow_non_primary: bool,
    pub(crate) num_used: usize,
    pub(crate) num_skipped: usize,
    pub(crate) num_failed: usize,
}

impl<'a, T: bam::Read> TrackingModRecordIter<'a, T> {
    pub(crate) fn new(
        records: bam::Records<'a, T>,
        skip_unmapped: bool,
        allow_non_primary: bool,
    ) -> Self {
        Self {
            records,
            skip_unmapped,
            allow_non_primary,
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
                    let should_skip = {
                        let based_on_primary = record_is_not_primary(&record)
                            && !self.allow_non_primary;
                        let based_on_unmapped =
                            record.is_unmapped() && self.skip_unmapped;
                        based_on_primary || based_on_unmapped
                    };
                    if should_skip {
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
                                             modification information, \
                                             skipping"
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
                                            "record {record_name} has \
                                             improper data, {}",
                                            e.to_string()
                                        );
                                        self.num_failed += 1;
                                        continue;
                                    }
                                    RunError::Failed(e) => {
                                        debug!(
                                            "record {record_name} failed to \
                                             extract mod base info, {}",
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

// todo deprecate this function or move it into the tracking iterator above
pub(crate) fn filter_records_iter<T: bam::Read>(
    records: bam::Records<T>,
) -> impl Iterator<Item = (bam::Record, ModBaseInfo)> + '_ {
    records
        .filter_map(|res| match res {
            Ok(rec) => Some(rec),
            Err(e) => {
                debug!("failed to read record from BAM, {}", e.to_string());
                None
            }
        })
        // skip non-primary
        .filter(|record| !record_is_not_primary(&record))
        // skip records with empty sequences
        .filter(|record| {
            if record.seq_len() > 0 {
                true
            } else {
                let query_name = get_query_name_string(&record)
                    .unwrap_or("'UTF-8 decode failure'".to_string());
                debug!("record {query_name} has empty seq");
                false
            }
        })
        .filter_map(|record| match ModBaseInfo::new_from_record(&record) {
            Ok(modbase_info) => {
                if modbase_info.is_empty() {
                    let query_name = get_query_name_string(&record)
                        .unwrap_or("'UTF-8 decode failure'".to_string());
                    debug!("{query_name} modbase info empty");
                    None
                } else {
                    Some((record, modbase_info))
                }
            }
            Err(e) => {
                let query_name = get_query_name_string(&record)
                    .unwrap_or("'UTF-8 decode failure'".to_string());
                error!(
                    "failed to get modbase info for record {query_name}, {}",
                    e.to_string()
                );
                None
            }
        })
}

pub const MM_TAGS: [&str; 2] = ["MM", "Mm"];
pub const ML_TAGS: [&str; 2] = ["ML", "Ml"];
pub const MN_TAG: &str = "MN";

pub struct RawModTags {
    raw_mm: String,
    raw_ml: Vec<u16>,
    mn_length: Option<usize>,
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
            mn_length: None,
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

    pub fn get_mn_length(&self) -> Option<usize> {
        self.mn_length
    }
}

#[derive(Debug, Clone)]
pub enum CollapseMethod {
    /// ModCode is the modified base to remove
    ReNormalize(ModCodeRepr),
    /// ModCode is the modified base to remove
    ReDistribute(ModCodeRepr),
    /// Convert one mod base to another
    Convert { from: HashSet<ModCodeRepr>, to: ModCodeRepr },
}

impl CollapseMethod {
    pub fn parse_str(
        raw: &str,
        mod_code: ModCodeRepr,
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
    Explicit,
    /// '.' mode, a.k.a 'implicit' no probability means the base is canonical
    /// (or predicted canonical).
    ProbModified,
    /// Same as `ProbModified` except the BAM record does not specify the
    /// actual mode.
    ImplicitProbModified,
}

impl SkipMode {
    fn parse(raw_mode: char) -> Result<Self, InputError> {
        match raw_mode {
            '?' => Ok(Self::Explicit),
            '.' => Ok(Self::ProbModified),
            _ => Err(InputError::new(&format!("unknown mode {}", raw_mode))),
        }
    }

    fn char(&self) -> Option<char> {
        match self {
            Self::Explicit => Some('?'),
            Self::ProbModified => Some('.'),
            Self::ImplicitProbModified => None,
        }
    }
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub enum BaseModCall {
    Canonical(f32),
    Modified(f32, ModCodeRepr),
    Filtered,
}

impl Eq for BaseModCall {}

impl Ord for BaseModCall {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).expect("should not have NaN probability")
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

impl BaseModCall {
    pub fn is_canonical(&self) -> bool {
        match self {
            BaseModCall::Canonical(_) => true,
            _ => false,
        }
    }

    pub fn is_match_modcall(&self, mod_code: &ModCodeRepr) -> bool {
        match self {
            BaseModCall::Modified(_, code) => code == mod_code,
            _ => false,
        }
    }
}

#[derive(new, Debug, PartialEq, Clone)]
pub struct BaseModProbs {
    probs: FxHashMap<ModCodeRepr, f32>,
    pub inferred: bool,
    // skip_mode: SkipMode,
    // strand: Strand,
}

impl BaseModProbs {
    pub fn new_init<T: Into<ModCodeRepr>>(mod_code: T, prob: f32) -> Self {
        Self {
            probs: FxHashMap::from_iter([(mod_code.into(), prob)]),
            inferred: false,
        }
    }

    pub fn new_inferred_canonical<
        'a,
        T: Into<ModCodeRepr> + Copy + Hash + 'a,
        IT: Iterator<Item = &'a T>,
    >(
        mod_codes: IT,
    ) -> Self {
        let probs = mod_codes.map(|code| ((*code).into(), 0f32)).collect();
        Self { probs, inferred: true }
    }

    pub fn insert_base_mod_prob(&mut self, mod_code: ModCodeRepr, prob: f32) {
        (*self.probs.entry(mod_code).or_insert(0f32)) += prob;
    }

    pub fn argmax_base_mod_call(&self) -> BaseModCall {
        let canonical_prob = self.canonical_prob();
        let max_mod_prob = self
            .iter_probs()
            .max_by(|(_, p), (_, q)| p.partial_cmp(q).unwrap());
        let base_mod_call = if let Some((mod_code, mod_prob)) = max_mod_prob {
            // let mod_code = ModCode::parse_raw_mod_code(*mod_code)?;
            if *mod_prob > canonical_prob {
                BaseModCall::Modified(*mod_prob, *mod_code)
            } else {
                BaseModCall::Canonical(canonical_prob)
            }
        } else {
            BaseModCall::Canonical(canonical_prob)
        };
        base_mod_call
    }

    pub fn canonical_prob(&self) -> f32 {
        1f32 - self.probs.values().sum::<f32>()
    }

    fn iter_codes(&self) -> impl Iterator<Item = &ModCodeRepr> {
        self.probs.keys()
    }

    // todo(arand): these methods should be removed/renamed to be more useful
    pub fn iter_probs(&self) -> impl Iterator<Item = (&ModCodeRepr, &f32)> {
        self.probs.iter()
    }

    pub fn iter_mut_probs(&mut self) -> impl Iterator<Item = &mut f32> {
        self.probs.iter_mut().map(|(_, p)| p)
    }

    pub fn iter_mut(
        &mut self,
    ) -> impl Iterator<Item = (&ModCodeRepr, &mut f32)> {
        self.probs.iter_mut()
    }

    pub(crate) fn into_collapsed(
        self,
        method: &CollapseMethod,
    ) -> BaseModProbs {
        let canonical_prob = self.canonical_prob();
        let inferred = self.inferred;
        match method {
            CollapseMethod::ReNormalize(mod_to_collapse) => {
                let marginal_collapsed_prob = self
                    .iter_probs()
                    .filter(|(mod_code, _prob)| *mod_code != mod_to_collapse)
                    .collect::<Vec<(&ModCodeRepr, &f32)>>();
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
                Self { probs, inferred }
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
                    .collect::<Vec<(&ModCodeRepr, &f32)>>();

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

                Self { probs, inferred }
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
                let mut new_base_mod_probs = Self { probs, inferred };

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
        Self { cumulative_counts, canonical_base: base }
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
    mod_base_codes: Vec<ModCodeRepr>,
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
            .chars()
            .peekable();

        let canonical_base = header
            .nth(0)
            .ok_or(InputError::new("failed to get canonical base"))?;

        let raw_stand =
            header.nth(0).ok_or(InputError::new("failed to get strand"))?;

        let strand = Strand::parse_char(raw_stand)?;

        let mut mod_base_codes = Vec::new();
        let mut mode: Option<SkipMode> = None;

        let mut offset = 2usize;
        let mut seen_chebi = false;

        let is_chebi =
            header.peek().map(|c| c.is_ascii_digit()).unwrap_or(false);
        if is_chebi {
            let mut agg = Vec::new();
            while let Some(d) = header.peeking_next(|d| d.is_ascii_digit()) {
                agg.push(d)
            }
            offset += agg.len();
            let raw_chebi = agg.into_iter().collect::<String>();
            let chebi_code = raw_chebi.parse::<u32>().map_err(|e| {
                InputError(format!(
                    "illegal chEBI code {raw_chebi}, {}",
                    e.to_string()
                ))
            })?;
            mod_base_codes.push(ModCodeRepr::ChEbi(chebi_code));
            seen_chebi = true;
        }

        while let Some(c) = header.next() {
            match c {
                '?' | '.' => {
                    mode = Some(SkipMode::parse(c).unwrap());
                    offset += 1;
                }
                _ => {
                    if c.is_ascii_digit() {
                        return Err(InputError::new(
                            "cannot have digit mod code, illegal MM tag \
                             {mod_positions}",
                        ));
                    } else {
                        if seen_chebi {
                            return Err(InputError(
                                "cannot combine chEBI codes and regular \
                                 codes, {header}"
                                    .to_string(),
                            ));
                        }
                        mod_base_codes.push(ModCodeRepr::Code(c));
                        offset += 1;
                    }
                }
            }
        }
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

        Ok(Self { canonical_base, mod_base_codes, mode, strand, delta_list })
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
            get_query_name_string(record).unwrap_or("???".to_string());
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
#[derive(PartialEq, Debug, Clone, new)]
pub struct SeqPosBaseModProbs {
    /// The `.` or `?` or implied mode, see `SkipMode`.
    pub skip_mode: SkipMode,
    /// Mapping of _forward_ sequence position to the predicted base
    /// modification probabilities for that position.
    pub pos_to_base_mod_probs: FxHashMap<usize, BaseModProbs>,
}

impl SeqPosBaseModProbs {
    fn new_empty(skip_mode: SkipMode) -> Self {
        Self { skip_mode, pos_to_base_mod_probs: FxHashMap::default() }
    }

    /// Remove positions that are outside the bounds of the `EdgeFilter`.
    /// Returning None means that either the read was too short or all
    /// of the positions were filtered out.
    pub(crate) fn edge_filter_positions(
        self,
        edge_filter: &EdgeFilter,
        read_length: usize,
    ) -> Option<Self> {
        if edge_filter.read_can_be_trimmed(read_length) {
            let pos_to_base_mod_probs = self
                .pos_to_base_mod_probs
                .into_iter()
                .filter(|(pos, _)| {
                    match edge_filter.keep_position(*pos, read_length) {
                        Ok(b) => b,
                        Err(_) => {
                            // shouldn't really happen,
                            false
                        }
                    }
                })
                .collect::<FxHashMap<usize, BaseModProbs>>();
            if pos_to_base_mod_probs.is_empty() {
                // all positions filtered out
                None
            } else {
                Some(Self::new(SkipMode::Explicit, pos_to_base_mod_probs))
            }
        } else {
            None
        }
    }

    // adds the implicit canonical calls when the mode is appropriate.
    fn add_implicit_mod_calls(
        self,
        delta_list: &[u32],
        all_mod_codes: &HashSet<ModCodeRepr>,
    ) -> Self {
        if self.skip_mode == SkipMode::ProbModified
            || self.skip_mode == SkipMode::ImplicitProbModified
        {
            let (probs, _) = delta_list.iter().enumerate().fold(
                (self.pos_to_base_mod_probs, 0u32),
                |(mut acc, cum_sum), (pos, d)| {
                    if *d > cum_sum {
                        acc.entry(pos).or_insert_with(|| {
                            BaseModProbs::new_inferred_canonical(
                                all_mod_codes.iter(),
                            )
                        });
                    }
                    (acc, *d)
                },
            );
            Self::new(self.skip_mode, probs)
        } else {
            self
        }
    }

    pub(crate) fn into_collapsed(
        self,
        collapse_method: &CollapseMethod,
    ) -> Self {
        let skip_mode = self.skip_mode;
        let pos_to_base_mod_probs = self
            .pos_to_base_mod_probs
            .into_iter()
            .map(|(pos, probs)| (pos, probs.into_collapsed(collapse_method)))
            .collect::<FxHashMap<usize, BaseModProbs>>();
        Self { skip_mode, pos_to_base_mod_probs }
    }

    pub(crate) fn set_skip_mode(&mut self, skip_mode: SkipMode) {
        self.skip_mode = skip_mode;
    }

    pub(crate) fn get_skip_mode(&self) -> SkipMode {
        self.skip_mode
    }

    /// removes implicit canonical probs and sets mode to Ambiguous, this is
    /// helpful when the initial mode is not provided
    /// `[SkipMode::ImplicitProbModified]`
    pub(crate) fn remove_implicit_probs(self) -> Self {
        let probs = self
            .pos_to_base_mod_probs
            .into_iter()
            .filter(|(_pos, base_mod_probs)| !base_mod_probs.inferred)
            .collect();
        Self::new(SkipMode::Explicit, probs)
    }
}

// todo(arand) remove, or put behind cfg(test)
pub fn extract_mod_probs(
    record: &bam::Record,
    raw_mm: &str,
    mod_quals: &[u16],
    converter: &DeltaListConverter,
) -> Result<SeqPosBaseModProbs, InputError> {
    // warn!("[deprecation warning] this method should not be called in
    // production code");
    let splited = raw_mm.split(";");
    let mut positions_to_probs =
        SeqPosBaseModProbs::new_empty(SkipMode::Explicit);
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
                positions_to_probs.insert(
                    position,
                    BaseModProbs::new_init(*mod_base_code, prob),
                );
            }
        }
    }

    Ok(SeqPosBaseModProbs::new(base_mod_positions.mode, positions_to_probs))
}

pub fn format_mm_ml_tag(
    positions_to_probs: SeqPosBaseModProbs,
    strand: Strand,
    converter: &DeltaListConverter,
) -> (String, Vec<u8>) {
    let canonical_base = converter.canonical_base;
    let skip_mode = positions_to_probs.skip_mode;
    let mut mod_code_to_position =
        HashMap::<(ModCodeRepr, Strand), Vec<(usize, f32)>>::new();

    for (position, mod_base_probs) in positions_to_probs.pos_to_base_mod_probs {
        // don't write down inferred base mod probs.
        if mod_base_probs.inferred
            && (skip_mode == SkipMode::ProbModified
                || skip_mode == SkipMode::ImplicitProbModified)
        {
            // add mod codes so that they are added later
            for mod_base_code in mod_base_probs.iter_codes() {
                mod_code_to_position
                    .entry((*mod_base_code, strand))
                    .or_insert_with(|| Vec::new());
            }
        } else {
            for (mod_base_code, mod_base_prob) in mod_base_probs.iter_probs() {
                let entry = mod_code_to_position
                    .entry((*mod_base_code, strand))
                    .or_insert(Vec::new());
                entry.push((position, *mod_base_prob));
            }
        }
    }

    let mut mm_tag = String::new();
    let mut ml_tag = Vec::new();
    let skip_mode_label =
        skip_mode.char().map(|s| s.to_string()).unwrap_or("".to_string());
    if mod_code_to_position.is_empty() {
        mm_tag.push_str(&format!(
            "{}{}{}{};",
            canonical_base,
            strand.to_char(),
            canonical_base, // "any mod" for a base is the same char as itself
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
            let header = {
                let mut header = format!(
                    "{}{}{}{}", // C+m?,
                    canonical_base,
                    strand.to_char(),
                    mod_code,
                    skip_mode_label
                );

                if !positions_and_probs.is_empty() {
                    // don't want to add this comma if there aren't any probs..
                    header.push(',');
                }
                header
            };
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
) -> Result<(String, &'static str), RunError> {
    get_tag::<String>(&record, &MM_TAGS, &parse_mm_tag)
}

pub fn get_ml_tag_from_record(
    record: &bam::Record,
) -> Result<(Vec<u16>, &'static str), RunError> {
    get_tag::<Vec<u16>>(&record, &ML_TAGS, &parse_ml_tag)
}

#[inline]
fn get_mn_tag_from_record(
    record: &bam::Record,
) -> Result<Option<usize>, RunError> {
    match record.aux(MN_TAG.as_bytes()) {
        Ok(Aux::U8(x)) => Ok(Some(x as usize)),
        Ok(Aux::U16(x)) => Ok(Some(x as usize)),
        Ok(Aux::U32(x)) => Ok(Some(x as usize)),
        Ok(Aux::I8(x)) => Ok(Some(x as usize)),
        Ok(Aux::I16(x)) => Ok(Some(x as usize)),
        Ok(Aux::I32(x)) => Ok(Some(x as usize)),
        Ok(_) => Err(RunError::new_input_error("MN invalid type")),
        Err(rust_htslib::errors::Error::BamAuxTagNotFound) => Ok(None),
        Err(e) => Err(RunError::new_failed(format!(
            "failed to parse MN tag, {}",
            e.to_string()
        ))),
    }
}

#[inline]
fn check_mn_tag_correct(
    record: &bam::Record,
    mn_tag: Option<usize>,
) -> Result<(), RunError> {
    match mn_tag {
        Some(l) if l != record.seq_len() => {
            return Err(RunError::new_input_error(format!(
                "MN tag length {} and seq length {} don't match",
                l,
                record.seq_len()
            )));
        }
        _ => {}
    }
    if record_is_not_primary(&record) && mn_tag.is_none() {
        return Err(RunError::new_skipped(
            "non-primary alignments must have MN tag",
        ));
    }
    Ok(())
}

fn validate_mn_tag_on_record(
    record: &bam::Record,
) -> Result<Option<usize>, RunError> {
    let mn_tag_value = get_mn_tag_from_record(record)?;
    check_mn_tag_correct(record, mn_tag_value).map(|_| mn_tag_value)
}

pub fn parse_raw_mod_tags(
    record: &bam::Record,
) -> Result<RawModTags, RunError> {
    let (raw_mm, mm_style) =
        get_mm_tag_from_record(record).map_err(|e| match e {
            RunError::Skipped(_) => e,
            _ => RunError::new_input_error("MM tag malformed"),
        })?;
    let (raw_ml, ml_style) =
        get_ml_tag_from_record(record).map_err(|e| match e {
            RunError::Skipped(_) => e,
            _ => RunError::new_input_error("ML tag malformed"),
        })?;
    let mn = validate_mn_tag_on_record(record)?;
    Ok(RawModTags { raw_mm, raw_ml, mn_length: mn, mm_style, ml_style })
}

pub struct ModBaseInfo {
    pub pos_seq_base_mod_probs: HashMap<char, SeqPosBaseModProbs>,
    pub neg_seq_base_mod_probs: HashMap<char, SeqPosBaseModProbs>,
    converters: HashMap<char, DeltaListConverter>,
    pub mm_style: &'static str,
    pub ml_style: &'static str,
}

impl ModBaseInfo {
    pub fn new_from_record(record: &bam::Record) -> Result<Self, RunError> {
        let raw_mod_tags = parse_raw_mod_tags(record)?;
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
        let mut strand_observed_mod_codes = HashMap::new();
        let mut pointer = 0usize;
        for raw_mm in mm.split(';').filter(|raw_mm| !raw_mm.is_empty()) {
            let base_mod_positions = BaseModPositions::parse(raw_mm)?;
            let converter = converters
                .entry(base_mod_positions.canonical_base)
                .or_insert_with(|| {
                    DeltaListConverter::new(
                        forward_seq,
                        base_mod_positions.canonical_base,
                    )
                });
            let strand = if base_mod_positions.is_positive_strand() {
                Strand::Positive
            } else {
                Strand::Negative
            };
            for mod_code in base_mod_positions.mod_base_codes.iter() {
                strand_observed_mod_codes
                    .entry(strand)
                    .or_insert(HashMap::new())
                    .entry(base_mod_positions.canonical_base)
                    .or_insert(HashSet::new())
                    .insert(*mod_code);
            }

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
        let pos_seq_base_mod_probs = Self::add_implicit_calls(
            pos_seq_base_mod_probs,
            &converters,
            strand_observed_mod_codes
                .get(&Strand::Positive)
                .unwrap_or(&HashMap::new()),
        );
        let neg_seq_base_mod_probs = Self::add_implicit_calls(
            neg_seq_base_mod_probs,
            &converters,
            strand_observed_mod_codes
                .get(&Strand::Negative)
                .unwrap_or(&HashMap::new()),
        );

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
        HashMap<char, DeltaListConverter>, // todo make this DnaBase
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
        self.pos_seq_base_mod_probs
            .values()
            .chain(self.neg_seq_base_mod_probs.values())
            .all(|p| p.pos_to_base_mod_probs.is_empty())
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

    #[inline(always)]
    fn add_implicit_calls(
        base_to_seq_base_mod_probs: HashMap<char, SeqPosBaseModProbs>,
        converters: &HashMap<char, DeltaListConverter>,
        all_mod_codes: &HashMap<char, HashSet<ModCodeRepr>>,
    ) -> HashMap<char, SeqPosBaseModProbs> {
        base_to_seq_base_mod_probs
            .into_iter()
            .map(|(primary_base, seq_base_mod_probs)| {
                let corrected = seq_base_mod_probs.add_implicit_mod_calls(
                    &converters
                        .get(&primary_base)
                        .expect("somehow missing delta list for {primary_base}")
                        .cumulative_counts,
                    all_mod_codes.get(&primary_base).unwrap_or(&HashSet::new()),
                );
                (primary_base, corrected)
            })
            .collect()
    }
}

#[cfg(test)]
pub fn base_mod_probs_from_record(
    record: &bam::Record,
    converter: &DeltaListConverter,
) -> Result<SeqPosBaseModProbs, RunError> {
    let (mm, ml) =
        parse_raw_mod_tags(record).map(|tags| (tags.raw_mm, tags.raw_ml))?;
    extract_mod_probs(record, &mm, &ml, &converter)
        .map_err(|input_err| RunError::BadInput(input_err))
}

#[derive(new, Debug)]
pub struct EdgeFilter {
    edge_filter_start: usize,
    edge_filter_end: usize,
    inverted: bool,
}

impl EdgeFilter {
    pub(crate) fn keep_position(
        &self,
        position: usize,
        read_length: usize,
    ) -> anyhow::Result<bool> {
        if !self.read_can_be_trimmed(read_length) {
            bail!(
                "read length not suitable for edge filter with start trim {}, \
                 end trim {} and read length {}, there should be a check \
                 before this call",
                self.edge_filter_start,
                self.edge_filter_end,
                read_length
            );
        } else if self.inverted {
            let before_start = position < self.edge_filter_start;
            let after_end = position >= (read_length - self.edge_filter_end);
            Ok(before_start || after_end)
        } else {
            let after_start = position >= self.edge_filter_start;
            let before_end = position < (read_length - self.edge_filter_end);
            Ok(after_start && before_end)
        }
    }

    #[inline]
    pub(crate) fn read_can_be_trimmed(&self, read_length: usize) -> bool {
        !(read_length <= self.edge_filter_start
            || read_length <= self.edge_filter_end)
    }
}

// DuplexPattern, None means it's a canonical call ('-')
#[derive(Copy, Clone, Debug, Hash, Eq, PartialEq, Ord, PartialOrd)]
pub enum DuplexModCodeRepr {
    Canonical,
    Code(char),
    ChEbi(u32),
}

impl From<ModCodeRepr> for DuplexModCodeRepr {
    fn from(value: ModCodeRepr) -> Self {
        match value {
            ModCodeRepr::Code(c) => Self::Code(c),
            ModCodeRepr::ChEbi(n) => Self::ChEbi(n),
        }
    }
}

impl Display for DuplexModCodeRepr {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ChEbi(x) => {
                write!(f, "{}", x)
            }
            Self::Code(x) => write!(f, "{}", x),
            Self::Canonical => {
                write!(f, "-")
            }
        }
    }
}

const CANONICAL_DUPLEX_PATTERN: [DuplexModCodeRepr; 2] =
    [DuplexModCodeRepr::Canonical, DuplexModCodeRepr::Canonical];

pub type DuplexPattern = [DuplexModCodeRepr; 2];

#[derive(Hash, Eq, PartialEq, Debug)]
pub(crate) enum DuplexModCall {
    ModCall { pattern: DuplexPattern, primary_base: char },
    Filtered { primary_base: char },
    NoCall { primary_base: char },
}

impl DuplexModCall {
    pub(crate) fn from_base_mod_calls(
        pos_base_mod_call: BaseModCall,
        neg_base_mod_call: BaseModCall,
        primary_base: char,
    ) -> Self {
        match (pos_base_mod_call, neg_base_mod_call) {
            (BaseModCall::Canonical(_), BaseModCall::Canonical(_)) => {
                Self::ModCall {
                    pattern: CANONICAL_DUPLEX_PATTERN,
                    primary_base,
                }
            }
            (BaseModCall::Canonical(_), BaseModCall::Modified(_, mod_code)) => {
                Self::ModCall {
                    pattern: [DuplexModCodeRepr::Canonical, mod_code.into()],
                    primary_base,
                }
            }
            (BaseModCall::Modified(_, mod_code), BaseModCall::Canonical(_)) => {
                Self::ModCall {
                    pattern: [mod_code.into(), DuplexModCodeRepr::Canonical],
                    primary_base,
                }
            }
            (
                BaseModCall::Modified(_, mod_code_pos),
                BaseModCall::Modified(_, mod_code_neg),
            ) => Self::ModCall {
                pattern: [mod_code_pos.into(), mod_code_neg.into()],
                primary_base,
            },
            (_, BaseModCall::Filtered) | (BaseModCall::Filtered, _) => {
                Self::Filtered { primary_base }
            }
        }
    }

    pub(crate) fn primary_base(&self) -> char {
        match self {
            Self::ModCall { pattern: _, primary_base } => *primary_base,
            Self::NoCall { primary_base } => *primary_base,
            Self::Filtered { primary_base } => *primary_base,
        }
    }

    pub(crate) fn is_canonical(&self) -> bool {
        match self {
            Self::ModCall { pattern, primary_base: _ } => {
                pattern == &CANONICAL_DUPLEX_PATTERN
            }
            _ => false,
        }
    }

    pub(crate) fn is_mod_call(&self) -> bool {
        match self {
            Self::ModCall { pattern, primary_base: _ } => {
                pattern != &CANONICAL_DUPLEX_PATTERN
            }
            _ => false,
        }
    }

    pub(crate) fn pattern(&self) -> Option<DuplexPattern> {
        match self {
            Self::ModCall { pattern, primary_base: _ } => Some(*pattern),
            _ => None,
        }
    }

    pub(crate) fn is_nocall(&self) -> bool {
        match self {
            Self::NoCall { primary_base: _ } => true,
            _ => false,
        }
    }

    pub(crate) fn is_filtered(&self) -> bool {
        match self {
            Self::Filtered { primary_base: _ } => true,
            _ => false,
        }
    }

    pub(crate) fn into_combined(self) -> Self {
        if let Some(pattern) = self.pattern() {
            if self.is_canonical() {
                self
            } else {
                let x = if pattern[0] == DuplexModCodeRepr::Canonical {
                    DuplexModCodeRepr::Canonical
                } else {
                    let any_mod_code = ModCodeRepr::any_mod_code(
                        &DnaBase::parse(self.primary_base()).unwrap(),
                    );
                    any_mod_code.into()
                };
                let y = if pattern[1] == DuplexModCodeRepr::Canonical {
                    DuplexModCodeRepr::Canonical
                } else {
                    let any_mod_code = ModCodeRepr::any_mod_code(
                        &DnaBase::parse(self.primary_base()).unwrap(),
                    );
                    any_mod_code.into()
                };
                let pattern = [x, y];
                Self::ModCall { pattern, primary_base: self.primary_base() }
            }
        } else {
            self
        }
    }
}

#[cfg(test)]
mod mod_bam_tests {
    use super::*;
    use std::collections::BTreeSet;

    fn qual_to_prob(qual: u16) -> f32 {
        let q = qual as f32;
        (q + 0.5f32) / 256f32
    }

    // first implementation that does not account for multiple mods in the MM
    // tag (i.e. C and A)) pub fn get_mod_probs_for_query_positions(
    //     mm: &str,
    //     canonical_base: char,
    //     mod_quals: &[u16],
    //     converter: &DeltaListConverter,
    // ) -> Result<HashMap<usize, BaseModProbs>, InputError> {
    //     // todo move this outside this function. should handle the case where
    // mods for another base     //  come first and the offset of mod_quals
    // is already handled     let filtered_mod_positions = mm
    //         .split(';')
    //         .filter(|positions| positions.starts_with(canonical_base))
    //         .collect::<Vec<&str>>();
    //
    //     let mut probs_for_positions = HashMap::<usize, BaseModProbs>::new();
    //     let mut prob_array_idx = 0usize;
    //     for mod_positions in filtered_mod_positions {
    //         let mut parts = mod_positions.split(',');
    //         let mut header = parts
    //             .nth(0)
    //             .ok_or(InputError::new(
    //                 "failed to get leader for base mod position line",
    //             ))?
    //             .chars();
    //
    //         let raw_stand = header
    //             .nth(1)
    //             .ok_or(InputError::new("failed to get strand"))?;
    //
    //         // TODO handle duplex
    //         let _strand = Strand::parse_char(raw_stand).unwrap();
    //
    //         let mut mod_base_codes = Vec::new();
    //         let mut _mode: Option<char> = None;
    //
    //         while let Some(c) = header.next() {
    //             match c {
    //                 '?' | '.' => {
    //                     _mode = Some(c);
    //                 }
    //                 _ => mod_base_codes.push(c),
    //             }
    //         }
    //
    //         // taking the liberty to think that a read wouldn't be larger
    //         // than 2**32 - 1 bases long
    //         let delta_list = parts
    //             .into_iter()
    //             .map(|raw_pos| raw_pos.parse::<u32>())
    //             .collect::<Result<Vec<u32>, _>>()
    //             .map_err(|e| {
    //                 InputError::new(&format!(
    //                     "failed to parse position list, {}",
    //                     e.to_string()
    //                 ))
    //             })?;
    //
    //         let positions = converter.to_positions(&delta_list)?;
    //         for mod_base in mod_base_codes {
    //             for pos in &positions {
    //                 let qual = mod_quals[prob_array_idx];
    //                 let prob = qual_to_prob(qual);
    //                 if let Some(base_mod_probs) =
    //                     probs_for_positions.get_mut(pos)
    //                 {
    //                     base_mod_probs.insert_base_mod_prob(mod_base, prob);
    //                 } else {
    //                     probs_for_positions.insert(
    //                         *pos,
    //                         BaseModProbs::new_init(mod_base, prob),
    //                     );
    //                 }
    //                 // consume from the ML array
    //                 prob_array_idx += 1;
    //             }
    //         }
    //     }
    //
    //     Ok(probs_for_positions)
    // }

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
        let probs = HashMap::from([('h'.into(), 0.85), ('m'.into(), 0.10)])
            .into_iter()
            .collect();

        let mod_base_probs = BaseModProbs { probs, inferred: false };
        let collapsed = mod_base_probs
            .clone()
            .into_collapsed(&CollapseMethod::ReDistribute('h'.into()));
        assert_eq!(
            collapsed.probs,
            vec![('m'.into(), 0.52500004)]
                .into_iter()
                .collect::<FxHashMap<ModCodeRepr, f32>>()
        );
        let collapsed = mod_base_probs
            .clone()
            .into_collapsed(&CollapseMethod::ReNormalize('h'.into()));
        assert_eq!(
            collapsed.probs,
            vec![('m'.into(), 0.6666669)]
                .into_iter()
                .collect::<FxHashMap<ModCodeRepr, f32>>()
        );

        let collapsed = mod_base_probs
            .clone()
            .into_collapsed(&CollapseMethod::ReNormalize('a'.into()));
        assert_eq!(&collapsed, &mod_base_probs);
        let collapsed = mod_base_probs
            .clone()
            .into_collapsed(&CollapseMethod::ReDistribute('a'.into()));
        assert_eq!(&collapsed, &mod_base_probs);
    }

    #[test]
    fn test_mod_prob_collapse_norm_examples() {
        let probs = vec![('h'.into(), 0.05273438), ('m'.into(), 0.03320312)]
            .into_iter()
            .collect();

        let mod_base_probs = BaseModProbs { probs, inferred: false };
        let collapsed = mod_base_probs
            .into_collapsed(&CollapseMethod::ReNormalize('h'.into()));
        assert_eq!(
            collapsed.probs,
            vec![('m'.into(), 0.035051543)]
                .into_iter()
                .collect::<FxHashMap<ModCodeRepr, f32>>()
        );
    }

    #[test]
    fn test_mod_prob_collapse_dist_examples() {
        let probs = vec![('h'.into(), 0.05273438), ('m'.into(), 0.03320312)]
            .into_iter()
            .collect();
        let mod_base_probs = BaseModProbs { probs, inferred: false };
        let collapsed = mod_base_probs
            .into_collapsed(&CollapseMethod::ReDistribute('h'.into()));
        assert_eq!(
            collapsed.probs,
            vec![('m'.into(), 0.059570313)]
                .into_iter()
                .collect::<FxHashMap<ModCodeRepr, f32>>()
        );
    }

    #[test]
    fn test_mod_prob_convert() {
        let inferred = false;
        let probs = vec![('h'.into(), 0.10), ('m'.into(), 0.75)]
            .into_iter()
            .collect::<FxHashMap<ModCodeRepr, f32>>();
        let mod_base_probs = BaseModProbs { probs: probs.clone(), inferred };

        let collapsed =
            mod_base_probs.into_collapsed(&CollapseMethod::Convert {
                from: HashSet::from(['h'.into()]),
                to: 'C'.into(),
            });
        assert_eq!(
            collapsed.probs,
            vec![('m'.into(), 0.75), ('C'.into(), 0.10)]
                .into_iter()
                .collect::<FxHashMap<ModCodeRepr, f32>>()
        );
        let mod_base_probs = BaseModProbs { probs: probs.clone(), inferred };
        let collapsed =
            mod_base_probs.into_collapsed(&CollapseMethod::Convert {
                from: HashSet::from(['h'.into(), 'm'.into()]),
                to: 'C'.into(),
            });
        assert_eq!(
            collapsed.probs,
            vec![('C'.into(), 0.85)]
                .into_iter()
                .collect::<FxHashMap<ModCodeRepr, f32>>()
        );
    }

    #[test]
    fn test_mod_prob_convert_sums_prob() {
        let probs =
            vec![('h'.into(), 0.10), ('m'.into(), 0.75)].into_iter().collect();
        let mod_base_probs = BaseModProbs { probs, inferred: false };
        let collapsed =
            mod_base_probs.into_collapsed(&CollapseMethod::Convert {
                from: HashSet::from(['h'.into()]),
                to: 'm'.into(),
            });
        assert_eq!(
            collapsed.probs,
            vec![('m'.into(), 0.85)]
                .into_iter()
                .collect::<FxHashMap<ModCodeRepr, f32>>()
        );
    }

    #[test]
    fn test_mod_prob_convert_noop() {
        let probs = vec![('h'.into(), 0.10), ('m'.into(), 0.75)]
            .into_iter()
            .collect::<FxHashMap<ModCodeRepr, f32>>();
        let mod_base_probs =
            BaseModProbs { probs: probs.clone(), inferred: false };
        let collapsed =
            mod_base_probs.into_collapsed(&CollapseMethod::Convert {
                from: HashSet::from(['a'.into()]),
                to: 'A'.into(),
            });
        assert_eq!(collapsed.probs, probs);
    }

    #[test]
    fn test_mod_prob_combine() {
        let inferred = false;
        let a_probs = vec![('h'.into(), 0.05273438), ('m'.into(), 0.03320312)]
            .into_iter()
            .collect();
        let mut a = BaseModProbs { probs: a_probs, inferred };
        let b_probs = vec![('m'.into(), 0.03320312)].into_iter().collect();
        let b = BaseModProbs { probs: b_probs, inferred };
        a.combine(b);
        assert_eq!(
            &a.probs,
            &vec![('h'.into(), 0.05273438), ('m'.into(), 0.06640624)]
                .into_iter()
                .collect::<FxHashMap<ModCodeRepr, f32>>()
        );

        let a_probs = vec![('m'.into(), 0.03320312)].into_iter().collect();
        let b_probs = vec![('h'.into(), 0.05273438)].into_iter().collect();

        let mut a = BaseModProbs { probs: a_probs, inferred };

        let b = BaseModProbs { probs: b_probs, inferred };
        a.combine(b);
        assert_eq!(
            &a.probs,
            &[('m'.into(), 0.03320312), ('h'.into(), 0.05273438)]
                .into_iter()
                .collect::<FxHashMap<ModCodeRepr, f32>>()
        );
    }

    #[test]
    #[ignore = "old implementation"]
    fn test_parse_mm_tag() {
        let _tag =
            "C+h?,5,2,1,3,1,2,3,1,2,1,11,5;C+m?,5,2,1,3,1,2,3,1,2,1,11,5;";
        // let dna = "
        // ATGTGCCTGCTGGACATGTTTATGCTCGTCTACTTCGTTCAGTTACGTATTGCTCCAG\
        //     CGCTCGAACTGTAGCCGCTGCTGCTGGGTGAAGTTGTGGCGGTACACGAGCTCCGCCGGCTGCAGCAGCTTC\
        //     TCCCCATCCTGGCGCTTCTCCCCGAGCAATTGGTG";
        // let mod_quals = vec![
        //     197, 13, 156, 1, 3, 5, 9, 26, 8, 1, 0, 13, 10, 67, 1, 0, 1, 0, 5,
        //     5, 5, 0, 0, 8,
        // ];
        //
        // let converter = DeltaListConverter::new(dna, 'C');
        // let positions_to_probs =
        //     get_mod_probs_for_query_positions(tag, 'C', &mod_quals,
        // &converter)         .unwrap();
        // assert_eq!(positions_to_probs.len(), 12);
        unimplemented!()
    }

    #[test]
    fn test_format_mm_ml_tags() {
        let canonical_base = 'C';
        //_________________________-12-34--5--6
        let read_sequence = "ACCGCCGTCGTCG";
        //_________________________0123456789012
        let converter = DeltaListConverter::new(read_sequence, canonical_base);

        let positions_and_probs = vec![
            (5, BaseModProbs::new_init('m', 0.9)),
            (2, BaseModProbs::new_init('m', 0.1)),
            (8, BaseModProbs::new_init('m', 0.2)),
        ]
        .into_iter()
        .collect::<FxHashMap<usize, BaseModProbs>>();

        let seq_pos_base_mod_probs =
            SeqPosBaseModProbs::new(SkipMode::Explicit, positions_and_probs);
        let (mm, ml) = format_mm_ml_tag(
            seq_pos_base_mod_probs,
            Strand::Positive,
            &converter,
        );
        assert_eq!(mm, "C+m?,1,1,0;");
        assert_eq!(ml, vec![25, 230, 51]);

        let skip_mode = SkipMode::ProbModified;
        let positions_and_probs = vec![
            (5, BaseModProbs::new_init('m', 0.9)),
            (2, BaseModProbs::new_init('m', 0.1)),
            (8, BaseModProbs::new_init('m', 0.2)),
        ]
        .into_iter()
        .collect::<FxHashMap<usize, BaseModProbs>>();

        let seq_pos_base_mod_probs =
            SeqPosBaseModProbs::new(skip_mode, positions_and_probs);
        let (mm, ml) = format_mm_ml_tag(
            seq_pos_base_mod_probs.clone(),
            Strand::Positive,
            &converter,
        );
        assert_eq!(mm, "C+m.,1,1,0;");
        assert_eq!(ml, vec![25, 230, 51]);
        let mod_codes =
            vec![ModCodeRepr::Code('m')].into_iter().collect::<HashSet<_>>();

        let seq_pos_base_mod_probs = seq_pos_base_mod_probs
            .add_implicit_mod_calls(&converter.cumulative_counts, &mod_codes);
        let (mm, ml) = format_mm_ml_tag(
            seq_pos_base_mod_probs.clone(),
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
            mode: SkipMode::Explicit,
            strand: Strand::Positive,
            mod_base_codes: vec!['h'.into()],
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
            mod_base_codes: vec!['m'.into()],
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
            mod_base_codes: vec!['m'.into()],
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
                &[('h'.into(), 0.005859375), ('m'.into(), 0.005859375)]
                    .into_iter()
                    .collect::<FxHashMap<ModCodeRepr, f32>>()
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
            let other =
                positions_to_probs.pos_to_base_mod_probs.get(position).unwrap();
            assert_eq!(base_mod_probs.probs, other.probs);
        }
    }

    #[test]
    fn test_mod_base_info() {
        let inferred = false;
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

        let c_expected_probs =
            vec![('h'.into(), 0.005859375), ('m'.into(), 0.39257813)]
                .into_iter()
                .collect();
        let c_expected = BaseModProbs { probs: c_expected_probs, inferred };
        let a_expected_probs =
            vec![('a'.into(), 0.7832031)].into_iter().collect();
        let a_expected = BaseModProbs { probs: a_expected_probs, inferred };
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

        for base_mod_probs in
            c_expected_seq_pos_base_mod_probs.pos_to_base_mod_probs.values()
        {
            assert_eq!(base_mod_probs, &c_expected);
        }

        for base_mod_probs in
            a_expected_seq_pos_base_mod_probs.pos_to_base_mod_probs.values()
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
        let inferred = false;
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
                    vec![('h'.into(), 0.39257813), ('m'.into(), 0.005859375)]
                        .into_iter()
                        .collect();
                let expected_modbase_probs =
                    BaseModProbs { probs: expected_probs, inferred };
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
                let expected_probs =
                    vec![('h'.into(), 0.5878906), ('m'.into(), 0.009765625)]
                        .into_iter()
                        .collect();
                let expected_modbase_probs =
                    BaseModProbs { probs: expected_probs, inferred };
                for mod_probs in probs.pos_to_base_mod_probs.values() {
                    assert_eq!(mod_probs, &expected_modbase_probs);
                }
            }
        }
    }

    #[test]
    fn test_duplex_modbase_info_implicit() {
        //               g c CG c gg CG CG
        let dna = "GACTCGACTGGACGTCGA";
        //               012345678901234567
        let tag = "C+h.,1,1,0;C+m.,1,1,0;G-h.,1,2,0;G-m.,1,2,0";
        let quals = vec![100, 100, 100, 1, 1, 1, 150, 150, 150, 2, 2, 2];
        let raw_mod_tags = RawModTags::new(tag, &quals, true);
        let obs_mod_base_info =
            ModBaseInfo::new(&raw_mod_tags, dna, &bam::Record::new()).unwrap();

        let top_strand_mods =
            obs_mod_base_info.pos_seq_base_mod_probs.get(&'C').unwrap();
        assert_eq!(top_strand_mods.pos_to_base_mod_probs.len(), 5);
        for pos in [2, 7] {
            let base_mod_probs =
                top_strand_mods.pos_to_base_mod_probs.get(&pos).unwrap();
            assert!(base_mod_probs.inferred);
            let probs = base_mod_probs.canonical_prob();
            assert_eq!(probs, 1.0f32);
        }

        let bottom_strand_mods =
            obs_mod_base_info.neg_seq_base_mod_probs.get(&'G').unwrap();
        assert_eq!(bottom_strand_mods.pos_to_base_mod_probs.len(), 6);
        for pos in [0, 9, 10] {
            let base_mod_probs =
                bottom_strand_mods.pos_to_base_mod_probs.get(&pos).unwrap();
            assert!(base_mod_probs.inferred);
            let probs = base_mod_probs.canonical_prob();
            assert_eq!(probs, 1.0f32);
        }
    }

    #[test]
    fn test_base_modcall_equality() {
        let a = BaseModCall::Canonical(1.0);
        let b = BaseModCall::Canonical(1.0);
        let c = BaseModCall::Modified(0.8, 'a'.into());
        let d = BaseModCall::Modified(0.7, 'a'.into());
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
        let c_seq_base_mod_probs =
            obs_mod_base_info.pos_seq_base_mod_probs.remove(&'C').unwrap();
        let expected_pos = vec![3, 9, 12];
        let obs_pos = c_seq_base_mod_probs
            .pos_to_base_mod_probs
            .keys()
            .map(|p| *p)
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        assert_eq!(expected_pos, obs_pos);
        let edge_filter = EdgeFilter::new(4, 4, false);
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
        let edge_filter = EdgeFilter::new(50, 50, false);
        let mut obs_mod_base_info =
            ModBaseInfo::new(&raw_mod_tags, dna, &bam::Record::new()).unwrap();
        let c_seq_base_mod_probs =
            obs_mod_base_info.pos_seq_base_mod_probs.remove(&'C').unwrap();
        let c_seq_base_mod_probs =
            c_seq_base_mod_probs.edge_filter_positions(&edge_filter, dna.len());
        assert!(c_seq_base_mod_probs.is_none());

        // trim with mod call _at_ the position to be trimmed
        let mut obs_mod_base_info =
            ModBaseInfo::new(&raw_mod_tags, dna, &bam::Record::new()).unwrap();
        let c_seq_base_mod_probs =
            obs_mod_base_info.pos_seq_base_mod_probs.remove(&'C').unwrap();
        let expected_pos = vec![3, 9, 12];
        let obs_pos = c_seq_base_mod_probs
            .pos_to_base_mod_probs
            .keys()
            .map(|p| *p)
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        assert_eq!(expected_pos, obs_pos);
        let edge_filter = EdgeFilter::new(3, 3, false);
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
    fn test_mod_bam_modbase_info_empty() {
        let dna = "GATCGACTACGTCGA";
        let tag = "C+h?;C+m?;";
        let quals = vec![];
        let raw_mod_tags = RawModTags::new(tag, &quals, true);
        let obs_mod_base_info =
            ModBaseInfo::new(&raw_mod_tags, dna, &bam::Record::new()).unwrap();
        assert!(obs_mod_base_info.is_empty());
        let tag = "C+h.;C+m.;";
        let raw_mod_tags = RawModTags::new(tag, &quals, true);
        let obs_mod_base_info =
            ModBaseInfo::new(&raw_mod_tags, dna, &bam::Record::new()).unwrap();
        assert!(!obs_mod_base_info.is_empty());
        //               g c CG c gg CG CG
        let dna = "GACTCGACTGGACGTCGA";
        //               012345678901234567
        let tag = "C+h?;C+m?;G-h?;G-m?;";
        // let quals = vec![100, 100, 100, 1, 1, 1, 150, 150, 150, 2, 2, 2];
        let raw_mod_tags = RawModTags::new(tag, &quals, true);
        let obs_mod_base_info =
            ModBaseInfo::new(&raw_mod_tags, dna, &bam::Record::new()).unwrap();
        assert!(obs_mod_base_info.is_empty());
        let tag = "C+h.;C+m.;G-h.;G-m.;";
        let raw_mod_tags = RawModTags::new(tag, &quals, true);
        let obs_mod_base_info =
            ModBaseInfo::new(&raw_mod_tags, dna, &bam::Record::new()).unwrap();
        assert!(!obs_mod_base_info.is_empty());
        assert_eq!(obs_mod_base_info.pos_seq_base_mod_probs.len(), 1);

        let top_strand_mods =
            obs_mod_base_info.pos_seq_base_mod_probs.get(&'C').unwrap();
        assert_eq!(top_strand_mods.pos_to_base_mod_probs.len(), 5);
        for pos in [2, 4, 7, 12, 15] {
            let base_mod_probs =
                top_strand_mods.pos_to_base_mod_probs.get(&pos).unwrap();
            assert!(base_mod_probs.inferred);
            let probs = base_mod_probs.canonical_prob();
            assert_eq!(probs, 1.0f32);
        }

        let bottom_strand_mods =
            obs_mod_base_info.neg_seq_base_mod_probs.get(&'G').unwrap();
        assert_eq!(bottom_strand_mods.pos_to_base_mod_probs.len(), 6);
        for pos in [0, 5, 9, 10, 13, 16] {
            let base_mod_probs =
                bottom_strand_mods.pos_to_base_mod_probs.get(&pos).unwrap();
            assert!(base_mod_probs.inferred);
            let probs = base_mod_probs.canonical_prob();
            assert_eq!(probs, 1.0f32);
        }
    }
}
