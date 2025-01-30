use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::fmt::{Debug, Display, Formatter};
use std::hash::Hash;

use anyhow::bail;
use derive_new::new;
use itertools::{Itertools, PeekingNext};
use log::debug;
use nom::bytes::complete::tag;
use nom::character::complete::{digit1, multispace0};
use nom::multi::separated_list1;
use nom::IResult;
use rayon::prelude::*;
use rust_htslib::bam;
use rust_htslib::bam::record::Aux;
use rustc_hash::FxHashMap;

use crate::errs::{ConflictError, MkError, MkResult};
use crate::find_motifs::iupac::nt_bytes;
use crate::mod_base_code::{DnaBase, ModCodeRepr, ParseChar};
use crate::util::{
    get_forward_sequence, get_tag, record_is_not_primary, Strand,
};

const MAX_PROB: f32 = 1.01f32;
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
                                "{record_name}: {}",
                                MkError::EmptyReadSequence
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
                                Err(e) => {
                                    debug!("{record_name}: {e}");
                                    self.num_failed += 1;
                                }
                            }
                        }
                    }
                }
                Err(e) => {
                    debug!(
                        "failed to read record from bam information, {}",
                        e.to_string()
                    );
                    self.num_failed += 1;
                }
            }
        }
        ret
    }
}

pub(crate) struct ModBaseInfoRecordTracker<
    I: Iterator<Item = rust_htslib::errors::Result<bam::record::Record>>,
> {
    // total: usize,
    num_errors: usize,
    records: I,
}

pub(crate) trait WithModBaseInfos<
    I: Iterator<Item = rust_htslib::errors::Result<bam::record::Record>>,
>
{
    fn with_mod_base_info(self) -> ModBaseInfoRecordTracker<Self>
    where
        Self: Iterator<Item = rust_htslib::errors::Result<bam::record::Record>>
            + Sized,
    {
        ModBaseInfoRecordTracker {
            // total: 0,
            num_errors: 0,
            records: self,
        }
    }
}

impl<I: Iterator<Item = rust_htslib::errors::Result<bam::record::Record>>>
    WithModBaseInfos<I> for I
{
}

impl<I: Iterator<Item = rust_htslib::errors::Result<bam::record::Record>>>
    Iterator for ModBaseInfoRecordTracker<I>
{
    type Item = (bam::Record, ModBaseInfo);

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(r) = self.records.next() {
            match r {
                Ok(record) => {
                    if record_is_not_primary(&record) || record.seq_len() == 0 {
                        continue;
                    }
                    match ModBaseInfo::new_from_record(&record) {
                        Ok(modbase_info) => {
                            if modbase_info.is_empty() {
                                continue;
                            } else {
                                return Some((record, modbase_info));
                            }
                        }
                        Err(e) => match e {
                            MkError::AuxMissing
                            | MkError::NoModifiedBaseInformation
                            | MkError::EmptyReadSequence => continue,
                            _ => {
                                self.num_errors += 1;
                                continue;
                            }
                        },
                    }
                }
                Err(_e) => {
                    self.num_errors += 1;
                    continue;
                }
            }
        }
        if self.num_errors > 0 {
            debug!(
                "{} records failed, consider checking mod base tags",
                self.num_errors
            );
        }
        None
    }
}

// todo deprecate this function or move it into the tracking iterator above
#[cfg(test)]
pub(crate) fn filter_records_iter<T: bam::Read>(
    records: bam::Records<T>,
) -> impl Iterator<Item = (bam::Record, ModBaseInfo)> + '_ {
    use crate::util::get_query_name_string;
    use log::error;
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

// TODO go back and make a ModTag struct that contains this info and the MMtag
// info
pub struct RawModTags {
    pub raw_mm: String,
    pub raw_ml: Vec<u16>,
    pub mn_length: Option<usize>,
    pub mm_style: &'static str,
    pub ml_style: &'static str,
}

impl RawModTags {
    pub fn new_from_record(record: &bam::Record) -> MkResult<Self> {
        parse_raw_mod_tags(record)
    }

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
    pub fn parse_str(raw: &str, mod_code: ModCodeRepr) -> MkResult<Self> {
        match raw {
            "norm" => Ok(Self::ReNormalize(mod_code)),
            "dist" => Ok(Self::ReDistribute(mod_code)),
            _ => Err(MkError::InvalidCollapseMethod),
        }
    }
}

#[derive(Debug, Copy, Clone, Eq, PartialEq, Hash, Ord, PartialOrd)]
pub enum SkipMode {
    /// '?' mode, no probability for a position means we have no information
    /// about base modifications at that position
    Explicit,
    /// '.' mode, a.k.a. 'implicit' no probability means the base is canonical
    /// (or predicted canonical).
    ImplicitUnmodified,
    /// Same as `ImplicitUnmodified` except the BAM record does not specify the
    /// actual mode.
    DefaultImplicitUnmodified,
}

impl SkipMode {
    fn parse(raw_mode: char) -> MkResult<Self> {
        match raw_mode {
            '?' => Ok(Self::Explicit),
            '.' => Ok(Self::ImplicitUnmodified),
            _ => Err(MkError::InvalidSkipMode),
        }
    }

    fn char(&self) -> Option<char> {
        match self {
            Self::Explicit => Some('?'),
            Self::ImplicitUnmodified => Some('.'),
            Self::DefaultImplicitUnmodified => None,
        }
    }

    fn is_implicit(&self) -> bool {
        self != &Self::Explicit
    }
}

impl Display for SkipMode {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        if let Some(c) = self.char() {
            write!(f, "{c}")
        } else {
            write!(f, "default-implicit-unmodified")
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
    /// Probability of each modification code
    probs: FxHashMap<ModCodeRepr, f32>,
    /// This position is implicitly inferred to be unmodified
    pub inferred_unmodified: bool,
}

impl BaseModProbs {
    pub fn new_init<T: Into<ModCodeRepr>>(mod_code: T, prob: f32) -> Self {
        Self {
            probs: FxHashMap::from_iter([(mod_code.into(), prob)]),
            inferred_unmodified: false,
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
        Self { probs, inferred_unmodified: true }
    }

    #[inline]
    pub(crate) fn add_base_mod_prob(
        &mut self,
        mod_code: ModCodeRepr,
        prob: f32,
    ) -> MkResult<()> {
        if self.inferred_unmodified && prob > 0f32 {
            return Err(MkError::Conflict(
                ConflictError::InferredSumGreaterThanOne,
            ));
        } else {
            let q = self.probs.entry(mod_code).or_insert(0f32);
            if *q + prob > MAX_PROB {
                // warn!("{q} {prob}");
                return Err(MkError::Conflict(
                    ConflictError::ProbaGreaterThanOne,
                ));
            } else {
                *q += prob;
                Ok(())
            }
        }
    }

    pub(crate) fn add_inferred_canonical<
        'a,
        T: Into<ModCodeRepr> + Copy + Hash + Debug + 'a,
        IT: Iterator<Item = &'a T>,
    >(
        &mut self,
        mod_codes: IT,
    ) -> MkResult<()> {
        if self.inferred_unmodified {
            for code in mod_codes {
                let prev = self.probs.insert((*code).into(), 0f32);
                if prev.map(|x| x > 0f32).unwrap_or(false) {
                    return Err(MkError::Conflict(
                        ConflictError::InferredSumGreaterThanOne,
                    ));
                }
            }
            Ok(())
        } else {
            Ok(())
        }
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
        let inferred = self.inferred_unmodified;
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
                Self { probs, inferred_unmodified: inferred }
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

                Self { probs, inferred_unmodified: inferred }
            }
            CollapseMethod::Convert { from, to } => {
                let (probs, converted_prob) = self.iter_probs().fold(
                    (FxHashMap::default(), 0f32),
                    |(mut probs, converted_prob), (raw_mod_code, prob)| {
                        // todo clean up?
                        if from.contains(&raw_mod_code) {
                            (probs, converted_prob + prob)
                        } else {
                            probs.insert(*raw_mod_code, *prob);
                            (probs, converted_prob)
                        }
                    },
                );
                let mut new_base_mod_probs =
                    Self { probs, inferred_unmodified: inferred };

                if converted_prob > 0f32 {
                    // safe because inferred is always canonical
                    new_base_mod_probs
                        .add_base_mod_prob(*to, converted_prob)
                        .unwrap();
                }

                new_base_mod_probs
            }
        }
    }

    fn combine_checked(&mut self, other: Self) -> MkResult<()> {
        if self.inferred_unmodified != other.inferred_unmodified {
            return Err(MkError::Conflict(
                ConflictError::ExplicitConflictInferred,
            ));
        }
        for (mod_code, &prob) in other.iter_probs() {
            let q = self.probs.entry(*mod_code).or_insert(0f32);
            *q += prob;
        }

        self.check()
    }

    fn check(&self) -> MkResult<()> {
        let x = self.probs.values().sum::<f32>();
        if x > MAX_PROB {
            // due to "centering" you can get a tiny but more than 1
            // warn!("{:?}, {x}", &self.probs);
            // for p in self.probs.values() {
            //     let q = prob_to_qual(*p);
            //     warn!("{p} {q}");
            // }
            Err(MkError::Conflict(ConflictError::ProbaGreaterThanOne))
        } else {
            Ok(())
        }
    }
}

pub struct DeltaListConverter {
    pub(crate) cumulative_counts: Vec<u32>,
    fundamental_base: FundamentalBase,
}

impl DeltaListConverter {
    /// "forward sequence" is the read sequence in the orientation that the
    /// instrument acquired it
    fn new(forward_sequence: &[u8], fundamental_base: FundamentalBase) -> Self {
        if fundamental_base == FundamentalBase::N {
            Self { cumulative_counts: Vec::new(), fundamental_base }
        } else {
            let cumulative_counts = forward_sequence
                .iter()
                .scan(0, |count, nt| {
                    if fundamental_base.matches(*nt) {
                        *count = *count + 1;
                    }
                    Some(*count)
                })
                .collect::<Vec<u32>>();

            debug_assert_eq!(cumulative_counts.len(), forward_sequence.len());
            Self { cumulative_counts, fundamental_base }
        }
    }

    pub fn new_base(forward_sequence: &[u8], base: DnaBase) -> Self {
        let fb = match base {
            DnaBase::A => FundamentalBase::A,
            DnaBase::C => FundamentalBase::C,
            DnaBase::G => FundamentalBase::G,
            DnaBase::T => FundamentalBase::T,
        };
        Self::new(forward_sequence, fb)
    }

    #[inline]
    fn to_positions_specific(
        &self,
        delta_list: &[u32],
    ) -> MkResult<Vec<usize>> {
        let mut finger = 0usize;
        let mut n_skips = 0u32;
        let mut positions = Vec::with_capacity(delta_list.len());
        for d in delta_list {
            if finger >= self.cumulative_counts.len() {
                return Err(MkError::InvalidMm(
                    "delta list refers to positions beyond end of seq"
                        .to_string(),
                ));
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
                    return Err(MkError::InvalidMm(
                        "delta list refers to positions beyond end of seq"
                            .to_string(),
                    ));
                }
            }
            positions.push(finger);
            n_skips += d + 1;
        }
        Ok(positions)
    }

    pub fn to_positions(
        &self,
        delta_list: &[u32],
        forward_seq: &[u8],
    ) -> MkResult<Vec<usize>> {
        if self.fundamental_base == FundamentalBase::N {
            if delta_list.is_empty() {
                Ok(Vec::new())
            } else {
                let first = delta_list[0] as usize;
                let lim = forward_seq.len();
                let (pos, _) =
                    delta_list.iter().skip(1).map(|x| *x as usize).try_fold(
                        (vec![first], first),
                        |(mut agg, last), next| {
                            let next_pos = last + next + 1usize;
                            if next_pos >= lim {
                                Err(MkError::InvalidMm(
                                    "refers to positions beyond end of seq"
                                        .to_string(),
                                ))
                            } else {
                                agg.push(next_pos);
                                Ok((agg, next_pos))
                            }
                        },
                    )?;
                Ok(pos)
            }
        } else {
            self.to_positions_specific(delta_list)
        }
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

fn positions_to_delta_list<'a, IT: Iterator<Item = &'a usize>>(
    positions: IT,
    cumulative_counts: &[u32],
) -> Vec<u32> {
    let mut last = 0;
    let mut delta_list = Vec::new();
    for pos in positions {
        let cumulative_count = cumulative_counts[*pos];
        let d = cumulative_count - last - 1;
        delta_list.push(d);
        last = cumulative_count
    }
    delta_list
}

#[inline]
pub fn prob_to_qual(prob: f32) -> u8 {
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

#[derive(Hash, Eq, PartialEq, Copy, Clone, Debug)]
enum FundamentalBase {
    A,
    C,
    G,
    T,
    U,
    N,
}

impl ParseChar for FundamentalBase {
    fn parse_char(c: char) -> MkResult<Self>
    where
        Self: Sized,
    {
        match c {
            'A' => Ok(Self::A),
            'C' => Ok(Self::C),
            'G' => Ok(Self::G),
            'T' => Ok(Self::T),
            'U' => Ok(Self::U),
            'N' => Ok(Self::N),
            _ => Err(MkError::InvalidMm(format!(
                "invalid symbol for fundamental base, {c}"
            ))),
        }
    }

    fn char(&self) -> char {
        match self {
            FundamentalBase::A => 'A',
            FundamentalBase::C => 'C',
            FundamentalBase::G => 'G',
            FundamentalBase::T => 'T',
            FundamentalBase::U => 'U',
            FundamentalBase::N => 'N',
        }
    }
}

impl FundamentalBase {
    fn matches(&self, b: u8) -> bool {
        match self {
            FundamentalBase::A => b == nt_bytes::A,
            FundamentalBase::C => b == nt_bytes::C,
            FundamentalBase::G => b == nt_bytes::G,
            FundamentalBase::T => b == nt_bytes::T,
            FundamentalBase::U => b == nt_bytes::T,
            FundamentalBase::N => true,
        }
    }
}

/// Container for the information in the MM and ML tags
#[derive(Debug, Eq, PartialEq)]
pub struct MmTagInfo {
    fundamental_base: FundamentalBase,
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

impl MmTagInfo {
    pub(crate) fn has_positions(&self) -> bool {
        if self.mode.is_implicit() {
            true
        } else {
            !self.delta_list.is_empty()
        }
    }

    pub(crate) fn parse_mm_tag(mm_tag: &str) -> MkResult<Vec<Self>> {
        // todo go back and see if making this parallel is faster
        mm_tag
            .split(';')
            .filter(|raw| !raw.is_empty())
            .map(|raw_mm| Self::parse(raw_mm))
            .collect::<MkResult<Vec<MmTagInfo>>>()
    }

    pub(crate) fn parse(mod_positions: &str) -> MkResult<Self> {
        let mut parts = mod_positions.split(',');
        let mut header = parts
            .nth(0)
            .ok_or(MkError::InvalidMm(
                "failed to get leader for base mod position line".to_string(),
            ))?
            .chars()
            .peekable();

        let fundamental_base = header
            .nth(0)
            .ok_or(MkError::InvalidMm(
                "failed to get canonical base".to_string(),
            ))
            .and_then(|c| FundamentalBase::parse_char(c))?;

        let raw_stand = header
            .nth(0)
            .ok_or(MkError::InvalidMm("failed to get strand".to_string()))?;

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
                MkError::InvalidMm(format!(
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
                        return Err(MkError::InvalidMm(format!(
                            "cannot have digit mod code, illegal MM tag \
                             {mod_positions}",
                        )));
                    } else {
                        if seen_chebi {
                            return Err(MkError::InvalidMm(
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
        let mode = mode.unwrap_or(SkipMode::DefaultImplicitUnmodified);

        let delta_list = if offset + 1 <= mod_positions.len() {
            let (_, raw_delta_list) = mod_positions.split_at(offset + 1);
            let (_, delta_list) =
                parse_int_list(raw_delta_list).map_err(|e| {
                    MkError::InvalidMm(format!(
                        "invalid MM delta list, {}",
                        e.to_string()
                    ))
                })?;
            delta_list
        } else {
            vec![]
        };

        Ok(Self { fundamental_base, mod_base_codes, mode, strand, delta_list })
    }

    fn stride(&self) -> usize {
        self.mod_base_codes.len()
    }

    fn size(&self) -> usize {
        self.delta_list.len() * self.mod_base_codes.len()
    }

    fn is_implicit(&self) -> bool {
        self.mode.is_implicit()
    }

    pub(crate) fn header(&self) -> String {
        if let Some(mode_repr) = self.mode.char() {
            format!(
                "{}{}{mode_repr}",
                self.fundamental_base.char(),
                self.strand.to_char()
            )
        } else {
            format!("{}{}", self.fundamental_base.char(), self.strand.to_char())
        }
    }
}

fn combine_positions_to_probs(
    agg: &mut SeqPosBaseModProbs,
    to_add: SeqPosBaseModProbs,
) -> MkResult<()> {
    if agg.skip_mode != to_add.skip_mode {
        agg.skip_mode = SkipMode::ImplicitUnmodified;
    }

    for (position, base_mod_probs) in to_add.pos_to_base_mod_probs.into_iter() {
        if let Some(probs) = agg.pos_to_base_mod_probs.get_mut(&position) {
            probs.combine_checked(base_mod_probs)?;
        } else {
            agg.pos_to_base_mod_probs.insert(position, base_mod_probs);
        }
    }

    Ok(())
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
                    edge_filter.keep_position(*pos, read_length).unwrap_or_else(
                        |_| {
                            // shouldn't really happen,
                            false
                        },
                    )
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
    #[cfg(test)]
    fn add_implicit_mod_calls(
        self,
        delta_list: &[u32],
        all_mod_codes: &HashSet<ModCodeRepr>,
    ) -> Self {
        if self.skip_mode == SkipMode::ImplicitUnmodified
            || self.skip_mode == SkipMode::DefaultImplicitUnmodified
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
            .filter(|(_pos, base_mod_probs)| {
                !base_mod_probs.inferred_unmodified
            })
            .collect();
        Self::new(SkipMode::Explicit, probs)
    }
}

// todo(arand) remove, or put behind cfg(test)
#[cfg(test)]
pub fn extract_mod_probs(
    _record: &bam::Record,
    forward_seq: &[u8],
    raw_mm: &str,
    mod_quals: &[u16],
    converter: &DeltaListConverter,
) -> MkResult<SeqPosBaseModProbs> {
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
        let mm_tag_info = MmTagInfo::parse(mod_positions)?;
        if mm_tag_info.fundamental_base == converter.fundamental_base {
            let mut base_mod_probs = get_base_mod_probs(
                &mm_tag_info,
                &mod_quals,
                pointer,
                &forward_seq,
                converter,
            )
            .unwrap();
            let dna_base =
                DnaBase::parse(mm_tag_info.fundamental_base.char()).unwrap();
            let seq_pos_base_mod_probs =
                base_mod_probs.remove(&dna_base).unwrap();
            combine_positions_to_probs(
                &mut positions_to_probs,
                seq_pos_base_mod_probs,
            )?;
        }
        pointer += mm_tag_info.delta_list.len() * mm_tag_info.stride();
    }

    Ok(positions_to_probs)
}

fn get_base_mod_probs(
    mm_tag_info: &MmTagInfo,
    mod_quals: &[u16],
    pointer: usize,
    forward_sequence: &[u8],
    converter: &DeltaListConverter,
) -> MkResult<FxHashMap<DnaBase, SeqPosBaseModProbs>> {
    let positions =
        converter.to_positions(&mm_tag_info.delta_list, forward_sequence)?;
    let end = pointer + mm_tag_info.size();
    if end > mod_quals.len() {
        return Err(MkError::InvalidMl(format!(
            "ML array too short, need {end} have {}",
            mod_quals.len()
        )));
    }
    let probs = {
        let mut probs = mod_quals[pointer..end]
            .iter()
            .map(|qual| *qual as f32)
            .collect::<Vec<f32>>();
        quals_to_probs(&mut probs);
        probs
    };

    let mut base_to_mod_probs =
        FxHashMap::<DnaBase, SeqPosBaseModProbs>::default();
    let stride = mm_tag_info.stride();
    debug_assert_eq!(probs.len() / stride, positions.len());
    for (chunk, position) in probs.chunks(stride).zip(positions) {
        assert_eq!(chunk.len(), stride);
        // safe due to check in .to_positions above
        let dna_base = DnaBase::try_from(forward_sequence[position])?;
        let seq_pos_base_mod_probs = base_to_mod_probs
            .entry(dna_base)
            .or_insert_with(|| SeqPosBaseModProbs::new_empty(mm_tag_info.mode));
        for (i, mod_base_code) in mm_tag_info.mod_base_codes.iter().enumerate()
        {
            let prob = chunk[i];
            if let Some(base_mod_probs) =
                seq_pos_base_mod_probs.pos_to_base_mod_probs.get_mut(&position)
            {
                base_mod_probs.add_base_mod_prob(*mod_base_code, prob)?;
            } else {
                seq_pos_base_mod_probs.pos_to_base_mod_probs.insert(
                    position,
                    BaseModProbs::new_init(*mod_base_code, prob),
                );
            }
        }
    }

    if mm_tag_info.is_implicit() {
        let mut cum_sum = 0u32;
        for (pos, x) in converter.cumulative_counts.iter().enumerate() {
            if *x > cum_sum {
                // a panic (OOB) here is an algorithm error
                let base = DnaBase::try_from(forward_sequence[pos])?;
                let seq_pos_base_mod_probs =
                    base_to_mod_probs.entry(base).or_insert_with(|| {
                        SeqPosBaseModProbs::new_empty(mm_tag_info.mode)
                    });
                if let Some(base_mod_probs) =
                    seq_pos_base_mod_probs.pos_to_base_mod_probs.get_mut(&pos)
                {
                    base_mod_probs.add_inferred_canonical(
                        mm_tag_info.mod_base_codes.iter(),
                    )?;
                } else {
                    seq_pos_base_mod_probs.pos_to_base_mod_probs.insert(
                        pos,
                        BaseModProbs::new_inferred_canonical(
                            mm_tag_info.mod_base_codes.iter(),
                        ),
                    );
                }
            }
            cum_sum = *x;
        }
    }

    Ok(base_to_mod_probs)
}

// todo put this in the DeltaListCoverter, it will be a better API to have
//  converter.format_mm_ml_tags(seq_pos_base_mod_probs)
pub fn format_mm_ml_tag(
    positions_to_probs: SeqPosBaseModProbs,
    primary_base: DnaBase,
    cumulative_counts: &[u32],
    strand: Strand,
) -> (String, Vec<u8>) {
    let skip_mode = positions_to_probs.skip_mode;
    let mut mod_code_to_position =
        HashMap::<(ModCodeRepr, Strand), Vec<(usize, f32)>>::new();

    for (position, mod_base_probs) in positions_to_probs.pos_to_base_mod_probs {
        // don't write down inferred base mod probs.
        if mod_base_probs.inferred_unmodified && skip_mode.is_implicit() {
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
            primary_base,
            strand.to_char(),
            primary_base, // "any mod" for a base is the same char as itself
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
                    primary_base,
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
            let positions = positions_and_probs.iter().map(|(pos, _prob)| pos);
            // todo maybe make this part of DeltaListConverter? (like it was)
            let delta_list =
                positions_to_delta_list(positions, cumulative_counts)
                    .into_iter()
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
fn parse_mm_tag(mm_aux: &Aux) -> MkResult<String> {
    match mm_aux {
        Aux::String(s) => Ok(s.to_string()),
        _ => Err(MkError::InvalidMm("wrong type".to_string())),
    }
}

/// tag keys should be the new then old tags, for example ["ML", "Ml"].
fn parse_ml_tag(ml_aux: &Aux) -> MkResult<Vec<u16>> {
    match ml_aux {
        Aux::ArrayU8(arr) => Ok(arr.iter().map(|x| x as u16).collect()),
        _ => Err(MkError::InvalidMl("wrong type".to_string())),
    }
}

pub fn get_mm_tag_from_record(
    record: &bam::Record,
) -> MkResult<(String, &'static str)> {
    get_tag::<String>(&record, &MM_TAGS, &parse_mm_tag)
}

pub fn get_ml_tag_from_record(
    record: &bam::Record,
) -> MkResult<(Vec<u16>, &'static str)> {
    get_tag::<Vec<u16>>(&record, &ML_TAGS, &parse_ml_tag)
}

#[inline]
fn get_mn_tag_from_record(record: &bam::Record) -> MkResult<Option<usize>> {
    match record.aux(MN_TAG.as_bytes()) {
        Ok(Aux::U8(x)) => Ok(Some(x as usize)),
        Ok(Aux::U16(x)) => Ok(Some(x as usize)),
        Ok(Aux::U32(x)) => Ok(Some(x as usize)),
        Ok(Aux::I8(x)) => Ok(Some(x as usize)),
        Ok(Aux::I16(x)) => Ok(Some(x as usize)),
        Ok(Aux::I32(x)) => Ok(Some(x as usize)),
        Ok(_) => Err(MkError::InvalidMn("wrong type".to_string())),
        Err(rust_htslib::errors::Error::BamAuxTagNotFound) => Ok(None),
        Err(e) => Err(MkError::HtsLibError(e)),
    }
}

#[inline]
fn check_mn_tag_correct(
    record: &bam::Record,
    mn_tag: Option<usize>,
) -> MkResult<()> {
    match mn_tag {
        Some(l) if l != record.seq_len() => {
            return Err(MkError::InvalidMn(format!(
                "MN tag length {} and seq length {} don't match",
                l,
                record.seq_len()
            )));
        }
        _ => {}
    }
    if record_is_not_primary(&record) && mn_tag.is_none() {
        return Err(MkError::NonPrimaryMissingMn);
    }
    Ok(())
}

fn validate_mn_tag_on_record(record: &bam::Record) -> MkResult<Option<usize>> {
    let mn_tag_value = get_mn_tag_from_record(record)?;
    check_mn_tag_correct(record, mn_tag_value).map(|_| mn_tag_value)
}

#[inline]
fn parse_raw_mod_tags(record: &bam::Record) -> MkResult<RawModTags> {
    let (raw_mm, mm_style) =
        get_mm_tag_from_record(record).map_err(|e| match e {
            MkError::AuxMissing => MkError::MmMissing,
            e @ _ => e,
        })?;
    let (raw_ml, ml_style) =
        get_ml_tag_from_record(record).map_err(|e| match e {
            MkError::AuxMissing => MkError::MlMissing,
            e @ _ => e,
        })?;
    let mn = validate_mn_tag_on_record(record)?;
    Ok(RawModTags { raw_mm, raw_ml, mn_length: mn, mm_style, ml_style })
}

pub struct ModBaseInfo {
    pub pos_seq_base_mod_probs: HashMap<DnaBase, SeqPosBaseModProbs>,
    pub neg_seq_base_mod_probs: HashMap<DnaBase, SeqPosBaseModProbs>,
    converters: HashMap<DnaBase, DeltaListConverter>,
    pub mm_style: &'static str,
    pub ml_style: &'static str,
}

impl ModBaseInfo {
    pub fn new_from_record(record: &bam::Record) -> MkResult<Self> {
        let raw_mod_tags = parse_raw_mod_tags(record)?;
        let forward_sequence = get_forward_sequence(record);
        let mm_tag_infos = MmTagInfo::parse_mm_tag(&raw_mod_tags.raw_mm)?;
        Self::new(&mm_tag_infos, &raw_mod_tags, &forward_sequence)
    }

    pub fn new(
        tag_infos: &[MmTagInfo],
        raw_mod_tags: &RawModTags,
        forward_seq: &[u8],
    ) -> MkResult<Self> {
        let raw_ml = &raw_mod_tags.raw_ml;

        // todo make these DnaBase keys..
        let mut pos_seq_base_mod_probs =
            HashMap::<DnaBase, SeqPosBaseModProbs>::new();
        let mut neg_seq_base_mod_probs =
            HashMap::<DnaBase, SeqPosBaseModProbs>::new();

        let mut converters = HashMap::new();
        let mut pointer = 0usize;
        for mm_tag_info in tag_infos {
            let converter = converters
                .entry(mm_tag_info.fundamental_base)
                .or_insert_with(|| {
                    DeltaListConverter::new(
                        forward_seq,
                        mm_tag_info.fundamental_base,
                    )
                });

            // implicit probs are added here!
            let base_mod_probs = get_base_mod_probs(
                &mm_tag_info,
                &raw_ml,
                pointer,
                forward_seq,
                &converter,
            )?;

            let seq_base_mod_probs = match mm_tag_info.strand {
                Strand::Positive => &mut pos_seq_base_mod_probs,
                Strand::Negative => &mut neg_seq_base_mod_probs,
            };

            for (base, to_add) in base_mod_probs {
                let agg = seq_base_mod_probs.entry(base).or_insert_with(|| {
                    SeqPosBaseModProbs::new_empty(mm_tag_info.mode)
                });
                combine_positions_to_probs(agg, to_add)?
            }

            pointer += mm_tag_info.delta_list.len() * mm_tag_info.stride();
        }

        let mut converters = converters
            .into_iter()
            .filter_map(
                |(fundamental_base, converter)| match fundamental_base {
                    FundamentalBase::A => Some((DnaBase::A, converter)),
                    FundamentalBase::C => Some((DnaBase::C, converter)),
                    FundamentalBase::G => Some((DnaBase::G, converter)),
                    FundamentalBase::T => Some((DnaBase::T, converter)),
                    FundamentalBase::U => Some((DnaBase::T, converter)),
                    FundamentalBase::N => None,
                },
            )
            .collect::<HashMap<DnaBase, DeltaListConverter>>();

        let extra_bases = pos_seq_base_mod_probs
            .keys()
            .chain(neg_seq_base_mod_probs.keys())
            .filter(|b| !converters.contains_key(b))
            .unique()
            .map(|b| FundamentalBase::parse_char(b.char()))
            .collect::<MkResult<Vec<FundamentalBase>>>()?;
        let extra_converters = extra_bases
            .into_par_iter()
            .map(|fb| DeltaListConverter::new(forward_seq, fb))
            .map(|converter| {
                (
                    DnaBase::parse(converter.fundamental_base.char()).unwrap(),
                    converter,
                )
            })
            .collect::<HashMap<DnaBase, DeltaListConverter>>();
        converters.extend(extra_converters.into_iter());

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
        HashMap<DnaBase, DeltaListConverter>, // todo make this DnaBase
        impl Iterator<Item = (DnaBase, Strand, SeqPosBaseModProbs)>,
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
    ) -> impl Iterator<Item = (DnaBase, Strand, &SeqPosBaseModProbs)> {
        let pos_iter = self.pos_seq_base_mod_probs.iter().map(
            |(canonical_base, seq_pos_base_mod_probs)| {
                (*canonical_base, Strand::Positive, seq_pos_base_mod_probs)
            },
        );
        let neg_iter = self.neg_seq_base_mod_probs.iter().map(
            |(canonical_base, seq_pos_base_mod_probs)| {
                (*canonical_base, Strand::Negative, seq_pos_base_mod_probs)
            },
        );
        pos_iter.chain(neg_iter)
    }
}

#[cfg(test)]
pub fn base_mod_probs_from_record(
    record: &bam::Record,
    converter: &DeltaListConverter,
) -> MkResult<SeqPosBaseModProbs> {
    let (mm, ml) =
        parse_raw_mod_tags(record).map(|tags| (tags.raw_mm, tags.raw_ml))?;
    let forward_seq = get_forward_sequence(record);
    extract_mod_probs(record, &forward_seq, &mm, &ml, &converter)
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
        let canonical_base = FundamentalBase::C;
        let read_sequence = "ACCGCCGTCGTCG";
        let converter =
            DeltaListConverter::new(read_sequence.as_bytes(), canonical_base);

        let ds = [1, 1, 0];
        let expected = [2, 5, 8];
        let obs =
            converter.to_positions(&ds, read_sequence.as_bytes()).unwrap();
        assert_eq!(&obs, &expected);
        let obs = converter.to_delta_list(&expected);
        assert_eq!(obs, &ds);

        let ds = [3, 0, 0];
        let expected = [5, 8, 11];
        let obs =
            converter.to_positions(&ds, read_sequence.as_bytes()).unwrap();
        assert_eq!(&obs, &expected);
        let obs = converter.to_delta_list(&expected);
        assert_eq!(obs, &ds);

        let ds = [3, 1];
        let expected = [5, 11];
        let obs =
            converter.to_positions(&ds, read_sequence.as_bytes()).unwrap();
        assert_eq!(&obs, &expected);
        let obs = converter.to_delta_list(&expected);
        assert_eq!(obs, &ds);
    }

    #[test]
    fn test_mod_prob_collapse() {
        let probs = HashMap::from([('h'.into(), 0.85), ('m'.into(), 0.10)])
            .into_iter()
            .collect();

        let mod_base_probs = BaseModProbs { probs, inferred_unmodified: false };
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

        let mod_base_probs = BaseModProbs { probs, inferred_unmodified: false };
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
        let mod_base_probs = BaseModProbs { probs, inferred_unmodified: false };
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
        let mod_base_probs = BaseModProbs {
            probs: probs.clone(),
            inferred_unmodified: inferred,
        };

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
        let mod_base_probs = BaseModProbs {
            probs: probs.clone(),
            inferred_unmodified: inferred,
        };
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
        let mod_base_probs = BaseModProbs { probs, inferred_unmodified: false };
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
            BaseModProbs { probs: probs.clone(), inferred_unmodified: false };
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
        let mut a =
            BaseModProbs { probs: a_probs, inferred_unmodified: inferred };
        let b_probs = vec![('m'.into(), 0.03320312)].into_iter().collect();
        let b = BaseModProbs { probs: b_probs, inferred_unmodified: inferred };
        a.combine_checked(b).unwrap();
        assert_eq!(
            &a.probs,
            &vec![('h'.into(), 0.05273438), ('m'.into(), 0.06640624)]
                .into_iter()
                .collect::<FxHashMap<ModCodeRepr, f32>>()
        );

        let a_probs = vec![('m'.into(), 0.03320312)].into_iter().collect();
        let b_probs = vec![('h'.into(), 0.05273438)].into_iter().collect();

        let mut a =
            BaseModProbs { probs: a_probs, inferred_unmodified: inferred };

        let b = BaseModProbs { probs: b_probs, inferred_unmodified: inferred };
        a.combine_checked(b).unwrap();
        assert_eq!(
            &a.probs,
            &[('m'.into(), 0.03320312), ('h'.into(), 0.05273438)]
                .into_iter()
                .collect::<FxHashMap<ModCodeRepr, f32>>()
        );
    }

    #[test]
    fn test_format_mm_ml_tags() {
        let canonical_base = FundamentalBase::C;
        //_________________________-12-34--5--6
        let read_sequence = "ACCGCCGTCGTCG";
        //_________________________0123456789012
        let converter =
            DeltaListConverter::new(read_sequence.as_bytes(), canonical_base);

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
            DnaBase::C,
            &converter.cumulative_counts,
            Strand::Positive,
        );
        assert_eq!(mm, "C+m?,1,1,0;");
        assert_eq!(ml, vec![25, 230, 51]);

        let skip_mode = SkipMode::ImplicitUnmodified;
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
            DnaBase::C,
            &converter.cumulative_counts,
            Strand::Positive,
        );
        assert_eq!(mm, "C+m.,1,1,0;");
        assert_eq!(ml, vec![25, 230, 51]);
        let mod_codes =
            vec![ModCodeRepr::Code('m')].into_iter().collect::<HashSet<_>>();

        let seq_pos_base_mod_probs = seq_pos_base_mod_probs
            .add_implicit_mod_calls(&converter.cumulative_counts, &mod_codes);
        let (mm, ml) = format_mm_ml_tag(
            seq_pos_base_mod_probs.clone(),
            DnaBase::C,
            &converter.cumulative_counts,
            Strand::Positive,
        );
        assert_eq!(mm, "C+m.,1,1,0;");
        assert_eq!(ml, vec![25, 230, 51]);
    }

    #[test]
    fn test_mod_parse_base_positions() {
        let raw_positions = "C+h?,5,2,1,3,1,2,3,1,2,1,11,5;";
        let base_mod_positions = MmTagInfo::parse(raw_positions).unwrap();
        let expected = MmTagInfo {
            fundamental_base: FundamentalBase::C,
            mode: SkipMode::Explicit,
            strand: Strand::Positive,
            mod_base_codes: vec!['h'.into()],
            delta_list: vec![5, 2, 1, 3, 1, 2, 3, 1, 2, 1, 11, 5],
        };

        assert_eq!(base_mod_positions, expected);

        let raw_positions = "C+m,5,2,1,3,1,2,3,1,2,1,11,5;";
        let base_mod_positions = MmTagInfo::parse(raw_positions).unwrap();
        let expected = MmTagInfo {
            fundamental_base: FundamentalBase::C,
            mode: SkipMode::DefaultImplicitUnmodified,
            strand: Strand::Positive,
            mod_base_codes: vec!['m'.into()],
            delta_list: vec![5, 2, 1, 3, 1, 2, 3, 1, 2, 1, 11, 5],
        };

        assert_eq!(base_mod_positions, expected);
        let raw_positions = "C+m.,5,2,1,3,1,2,3,1,2,1,11,5;";
        let base_mod_positions = MmTagInfo::parse(raw_positions).unwrap();
        let expected = MmTagInfo {
            fundamental_base: FundamentalBase::C,
            mode: SkipMode::ImplicitUnmodified,
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
        let canonical_base = FundamentalBase::C;
        let converter = DeltaListConverter::new(dna.as_bytes(), canonical_base);
        // let mut record = bam::Record::new();
        // record.set(b"read", None, dna.as_bytes(), &vec![255u8; dna.len()]);
        let positions_to_probs = extract_mod_probs(
            &bam::Record::new(),
            dna.as_bytes(),
            tag,
            &quals,
            &converter,
        )
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

        let positions_to_probs_1 = extract_mod_probs(
            &bam::Record::new(),
            dna.as_bytes(),
            tag,
            &quals,
            &converter,
        )
        .unwrap();

        assert_eq!(positions_to_probs, positions_to_probs_1);
    }

    #[test]
    fn test_extract_positions_to_probs() {
        let dna = "GATCGACTACGTCGA";
        let tag = "C+h?,0,1,0;A+a?,0,1,0;C+m?,0,1,0;";
        let quals = vec![1, 1, 1, 200, 200, 200, 1, 1, 1];
        let canonical_base = FundamentalBase::C;
        let converter = DeltaListConverter::new(dna.as_bytes(), canonical_base);

        let positions_to_probs = extract_mod_probs(
            &bam::Record::new(),
            dna.as_bytes(),
            tag,
            &quals,
            &converter,
        )
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
        let positions_to_probs_comb = extract_mod_probs(
            &bam::Record::new(),
            dna.as_bytes(),
            tag,
            &quals,
            &converter,
        )
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
        let c_converter =
            DeltaListConverter::new(dna.as_bytes(), FundamentalBase::C);
        let a_converter =
            DeltaListConverter::new(dna.as_bytes(), FundamentalBase::A);

        // these tags only have 1 canonical base, parse these and these are the
        // expected values for the rest of the test
        let c_tag = "C+hm?,0,1,0;";
        let a_tag = "A+a?,0,1,0;";
        let c_quals = vec![1, 100, 1, 100, 1, 100]; // interleaved!
        let a_quals = vec![200, 200, 200];

        let c_expected_seq_pos_base_mod_probs = extract_mod_probs(
            &bam::Record::new(),
            dna.as_bytes(),
            c_tag,
            &c_quals,
            &c_converter,
        )
        .unwrap();
        let a_expected_seq_pos_base_mod_probs = extract_mod_probs(
            &bam::Record::new(),
            dna.as_bytes(),
            a_tag,
            &a_quals,
            &a_converter,
        )
        .unwrap();

        let c_expected_probs =
            vec![('h'.into(), 0.005859375), ('m'.into(), 0.39257813)]
                .into_iter()
                .collect();
        let c_expected = BaseModProbs {
            probs: c_expected_probs,
            inferred_unmodified: inferred,
        };
        let a_expected_probs =
            vec![('a'.into(), 0.7832031)].into_iter().collect();
        let a_expected = BaseModProbs {
            probs: a_expected_probs,
            inferred_unmodified: inferred,
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
        let mm_tag_infos = MmTagInfo::parse_mm_tag(tag).unwrap();
        let obs_mod_base_info =
            ModBaseInfo::new(&mm_tag_infos, &raw_mod_tags, dna.as_bytes())
                .unwrap();
        assert_eq!(
            obs_mod_base_info.pos_seq_base_mod_probs.get(&DnaBase::C).unwrap(),
            &c_expected_seq_pos_base_mod_probs
        );
        assert_eq!(
            obs_mod_base_info.pos_seq_base_mod_probs.get(&DnaBase::A).unwrap(),
            &a_expected_seq_pos_base_mod_probs
        );

        let obs_base_mod_probs = extract_mod_probs(
            &bam::Record::new(),
            dna.as_bytes(),
            tag,
            &quals,
            &c_converter,
        )
        .unwrap();
        assert_eq!(&obs_base_mod_probs, &c_expected_seq_pos_base_mod_probs);
        let obs_base_mod_probs = extract_mod_probs(
            &bam::Record::new(),
            dna.as_bytes(),
            tag,
            &quals,
            &a_converter,
        )
        .unwrap();
        assert_eq!(&obs_base_mod_probs, &a_expected_seq_pos_base_mod_probs);

        let tag = "C+h?,0,1,0;C+m?,0,1,0;A+a?,0,1,0;";
        let mm_tag_infos = MmTagInfo::parse_mm_tag(tag).unwrap();
        let quals = vec![1, 1, 1, 100, 100, 100, 200, 200, 200];
        let raw_mod_tags = RawModTags::new(tag, &quals, true);
        let obs_mod_base_info =
            ModBaseInfo::new(&mm_tag_infos, &raw_mod_tags, dna.as_bytes())
                .unwrap();
        assert_eq!(
            obs_mod_base_info.pos_seq_base_mod_probs.get(&DnaBase::C).unwrap(),
            &c_expected_seq_pos_base_mod_probs
        );
        assert_eq!(
            obs_mod_base_info.pos_seq_base_mod_probs.get(&DnaBase::A).unwrap(),
            &a_expected_seq_pos_base_mod_probs
        );

        let obs_base_mod_probs = extract_mod_probs(
            &bam::Record::new(),
            dna.as_bytes(),
            tag,
            &quals,
            &c_converter,
        )
        .unwrap();
        assert_eq!(&obs_base_mod_probs, &c_expected_seq_pos_base_mod_probs);
        let obs_base_mod_probs = extract_mod_probs(
            &bam::Record::new(),
            dna.as_bytes(),
            tag,
            &quals,
            &a_converter,
        )
        .unwrap();
        assert_eq!(&obs_base_mod_probs, &a_expected_seq_pos_base_mod_probs);

        // test with the mods "combined"
        let tag = "C+hm?,0,1,0;A+a?,0,1,0;";
        let mm_tag_infos = MmTagInfo::parse_mm_tag(tag).unwrap();
        let quals = vec![1, 100, 1, 100, 1, 100, 200, 200, 200];
        let raw_mod_tags = RawModTags::new(tag, &quals, true);
        let obs_mod_base_info =
            ModBaseInfo::new(&mm_tag_infos, &raw_mod_tags, dna.as_bytes())
                .unwrap();
        assert_eq!(
            obs_mod_base_info.pos_seq_base_mod_probs.get(&DnaBase::C).unwrap(),
            &c_expected_seq_pos_base_mod_probs
        );
        assert_eq!(
            obs_mod_base_info.pos_seq_base_mod_probs.get(&DnaBase::A).unwrap(),
            &a_expected_seq_pos_base_mod_probs
        );
        let obs_base_mod_probs = extract_mod_probs(
            &bam::Record::new(),
            dna.as_bytes(),
            tag,
            &quals,
            &c_converter,
        )
        .unwrap();
        assert_eq!(&obs_base_mod_probs, &c_expected_seq_pos_base_mod_probs);
        let obs_base_mod_probs = extract_mod_probs(
            &bam::Record::new(),
            dna.as_bytes(),
            tag,
            &quals,
            &a_converter,
        )
        .unwrap();
        assert_eq!(&obs_base_mod_probs, &a_expected_seq_pos_base_mod_probs);
    }

    #[test]
    fn test_duplex_modbase_info() {
        //               g c CG c  gg CG CG
        let dna = "GACTCGACTGGACGTCGA";
        let tag = "C+h?,1,1,0;C+m?,1,1,0;G-h?,1,2,0;G-m?,1,2,0";
        let mm_tag_infos = MmTagInfo::parse_mm_tag(tag).unwrap();
        let quals = vec![100, 100, 100, 1, 1, 1, 150, 150, 150, 2, 2, 2];
        let tags = RawModTags::new(tag, &quals, true);
        let info =
            ModBaseInfo::new(&mm_tag_infos, &tags, dna.as_bytes()).unwrap();
        let inferred = false;
        let (_converters, iterator) = info.into_iter_base_mod_probs();
        for (c, strand, probs) in iterator {
            if c == DnaBase::C {
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
                let expected_modbase_probs = BaseModProbs {
                    probs: expected_probs,
                    inferred_unmodified: inferred,
                };
                for mod_probs in probs.pos_to_base_mod_probs.values() {
                    assert_eq!(mod_probs, &expected_modbase_probs);
                }
            } else {
                assert_eq!(c, DnaBase::G);
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
                let expected_modbase_probs = BaseModProbs {
                    probs: expected_probs,
                    inferred_unmodified: inferred,
                };
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
        let mm_tag_infos = MmTagInfo::parse_mm_tag(tag).unwrap();
        let quals = vec![100, 100, 100, 1, 1, 1, 150, 150, 150, 2, 2, 2];
        let raw_mod_tags = RawModTags::new(tag, &quals, true);
        let obs_mod_base_info =
            ModBaseInfo::new(&mm_tag_infos, &raw_mod_tags, dna.as_bytes())
                .unwrap();

        let top_strand_mods =
            obs_mod_base_info.pos_seq_base_mod_probs.get(&DnaBase::C).unwrap();
        assert_eq!(top_strand_mods.pos_to_base_mod_probs.len(), 5);
        for pos in [2, 7] {
            let base_mod_probs =
                top_strand_mods.pos_to_base_mod_probs.get(&pos).unwrap();
            assert!(base_mod_probs.inferred_unmodified);
            let probs = base_mod_probs.canonical_prob();
            assert_eq!(probs, 1.0f32);
        }

        let bottom_strand_mods =
            obs_mod_base_info.neg_seq_base_mod_probs.get(&DnaBase::G).unwrap();
        assert_eq!(bottom_strand_mods.pos_to_base_mod_probs.len(), 6);
        for pos in [0, 9, 10] {
            let base_mod_probs =
                bottom_strand_mods.pos_to_base_mod_probs.get(&pos).unwrap();
            assert!(base_mod_probs.inferred_unmodified);
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
        let mm_tag_infos = MmTagInfo::parse_mm_tag(tag).unwrap();
        let quals = vec![1, 1, 1, 200, 200, 200, 100, 100, 100];
        let raw_mod_tags = RawModTags::new(tag, &quals, true);
        let mut obs_mod_base_info =
            ModBaseInfo::new(&mm_tag_infos, &raw_mod_tags, dna.as_bytes())
                .unwrap();
        let c_seq_base_mod_probs = obs_mod_base_info
            .pos_seq_base_mod_probs
            .remove(&DnaBase::C)
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
            ModBaseInfo::new(&mm_tag_infos, &raw_mod_tags, dna.as_bytes())
                .unwrap();
        let c_seq_base_mod_probs = obs_mod_base_info
            .pos_seq_base_mod_probs
            .remove(&DnaBase::C)
            .unwrap();
        let c_seq_base_mod_probs =
            c_seq_base_mod_probs.edge_filter_positions(&edge_filter, dna.len());
        assert!(c_seq_base_mod_probs.is_none());

        // trim with mod call _at_ the position to be trimmed
        let mut obs_mod_base_info =
            ModBaseInfo::new(&mm_tag_infos, &raw_mod_tags, dna.as_bytes())
                .unwrap();
        let c_seq_base_mod_probs = obs_mod_base_info
            .pos_seq_base_mod_probs
            .remove(&DnaBase::C)
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
        let quals = vec![];
        let tag = "C+h?;C+m?;";
        let mm_tag_infos = MmTagInfo::parse_mm_tag(tag).unwrap();
        let raw_mod_tags = RawModTags::new(tag, &quals, true);
        let obs_mod_base_info =
            ModBaseInfo::new(&mm_tag_infos, &raw_mod_tags, dna.as_bytes())
                .unwrap();
        assert!(obs_mod_base_info.is_empty());
        let tag = "C+h.;C+m.;";
        let mm_tag_infos = MmTagInfo::parse_mm_tag(tag).unwrap();
        let raw_mod_tags = RawModTags::new(tag, &quals, true);
        let obs_mod_base_info =
            ModBaseInfo::new(&mm_tag_infos, &raw_mod_tags, dna.as_bytes())
                .unwrap();
        assert!(!obs_mod_base_info.is_empty());
        //               g c CG c gg CG CG
        let dna = "GACTCGACTGGACGTCGA";
        //               012345678901234567
        let tag = "C+h?;C+m?;G-h?;G-m?;";
        // let quals = vec![100, 100, 100, 1, 1, 1, 150, 150, 150, 2, 2, 2];
        let mm_tag_infos = MmTagInfo::parse_mm_tag(tag).unwrap();
        let raw_mod_tags = RawModTags::new(tag, &quals, true);
        let obs_mod_base_info =
            ModBaseInfo::new(&mm_tag_infos, &raw_mod_tags, dna.as_bytes())
                .unwrap();
        assert!(obs_mod_base_info.is_empty());
        let tag = "C+h.;C+m.;G-h.;G-m.;";
        let mm_tag_infos = MmTagInfo::parse_mm_tag(tag).unwrap();
        let raw_mod_tags = RawModTags::new(tag, &quals, true);
        let obs_mod_base_info =
            ModBaseInfo::new(&mm_tag_infos, &raw_mod_tags, dna.as_bytes())
                .unwrap();
        assert!(!obs_mod_base_info.is_empty());
        assert_eq!(obs_mod_base_info.pos_seq_base_mod_probs.len(), 1);

        let top_strand_mods =
            obs_mod_base_info.pos_seq_base_mod_probs.get(&DnaBase::C).unwrap();
        assert_eq!(top_strand_mods.pos_to_base_mod_probs.len(), 5);
        for pos in [2, 4, 7, 12, 15] {
            let base_mod_probs =
                top_strand_mods.pos_to_base_mod_probs.get(&pos).unwrap();
            assert!(base_mod_probs.inferred_unmodified);
            let probs = base_mod_probs.canonical_prob();
            assert_eq!(probs, 1.0f32);
        }

        let bottom_strand_mods =
            obs_mod_base_info.neg_seq_base_mod_probs.get(&DnaBase::G).unwrap();
        assert_eq!(bottom_strand_mods.pos_to_base_mod_probs.len(), 6);
        for pos in [0, 5, 9, 10, 13, 16] {
            let base_mod_probs =
                bottom_strand_mods.pos_to_base_mod_probs.get(&pos).unwrap();
            assert!(base_mod_probs.inferred_unmodified);
            let probs = base_mod_probs.canonical_prob();
            assert_eq!(probs, 1.0f32);
        }
    }

    #[test]
    fn test_quals_and_probs() {
        let qs = (0u8..=255u8).collect::<Vec<u8>>();
        let mut ps = qs.iter().map(|q| *q as f32).collect::<Vec<f32>>();
        quals_to_probs(&mut ps);
        let qs2 = ps.into_iter().map(|p| prob_to_qual(p)).collect::<Vec<u8>>();
        assert_eq!(qs, qs2);
    }

    #[test]
    fn test_delta_list_converter_n_base() {
        let dna = "GCGGATTTCTGAGTTTG";
        let delta_list = vec![5u32, 0, 0, 1, 3, 0, 0];
        let converter =
            DeltaListConverter::new(dna.as_bytes(), FundamentalBase::N);
        assert!(converter.cumulative_counts.is_empty());
        let positions =
            converter.to_positions(&delta_list, dna.as_bytes()).unwrap();
        assert_eq!(positions, vec![5, 6, 7, 9, 13, 14, 15]);
    }

    #[test]
    fn test_generic_mm_tags() {
        //               00000123344445677;
        let dna = "GCGGATTTCTGAGTTTG";
        let mm = "N+b?,5,0,0,1,3,0,0;";
        let mm_tag_infos = MmTagInfo::parse_mm_tag(mm).unwrap();
        let mut record = bam::Record::new();
        record.set(b"test", None, dna.as_bytes(), &vec![255; dna.len()]);
        let raw_mm_tags = RawModTags::new(mm, &vec![255u16; 7], true);
        let modbase_info =
            ModBaseInfo::new(&mm_tag_infos, &raw_mm_tags, dna.as_bytes())
                .unwrap();
        let t_converter = modbase_info.converters.get(&DnaBase::T).unwrap();
        let expected = vec![0, 0, 0, 0, 0, 1, 2, 3, 3, 4, 4, 4, 4, 5, 6, 7, 7];
        assert_eq!(&expected, &t_converter.cumulative_counts);
        assert_eq!(modbase_info.converters.len(), 1);
    }

    #[test]
    fn test_generic_mm_tags_multibase() {
        //         T     00000123344445677;
        //               01111111222222222
        let dna = "GCGGATTTCTGAGTTTG";
        let mm = "N+b?,1,3,0,0,1,3,0,0;";
        let mm_tag_infos = MmTagInfo::parse_mm_tag(mm).unwrap();
        let mut record = bam::Record::new();
        record.set(b"test", None, dna.as_bytes(), &vec![255; dna.len()]);
        let raw_mm_tags = RawModTags::new(mm, &vec![255u16; 8], true);
        let modbase_info =
            ModBaseInfo::new(&mm_tag_infos, &raw_mm_tags, dna.as_bytes())
                .unwrap();
        let t_converter = modbase_info.converters.get(&DnaBase::T).unwrap();
        let expected = vec![0, 0, 0, 0, 0, 1, 2, 3, 3, 4, 4, 4, 4, 5, 6, 7, 7];
        assert_eq!(&expected, &t_converter.cumulative_counts);
        let c_converter = modbase_info.converters.get(&DnaBase::C).unwrap();
        let expected = vec![0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2];
        assert_eq!(&expected, &c_converter.cumulative_counts);

        assert_eq!(modbase_info.converters.len(), 2);
    }

    #[test]
    fn test_generic_mm_tags_multibase_conflict() {
        let dna = "GCGGATTTCTGAGTTTG";
        let mm = "N+b?,1,3,0,0,1,3,0,0;C+m?,0;";
        let mm_tag_infos = MmTagInfo::parse_mm_tag(mm).unwrap();
        let mut record = bam::Record::new();
        record.set(b"test", None, dna.as_bytes(), &vec![255; dna.len()]);
        let raw_mm_tags = RawModTags::new(mm, &vec![255u16; 9], true);
        let modbase_info =
            ModBaseInfo::new(&mm_tag_infos, &raw_mm_tags, dna.as_bytes());
        assert!(modbase_info.is_err());
        let mm = "C+m.;N+b?,1,3,0,0,1,3,0,0;";
        let mm_tag_infos = MmTagInfo::parse_mm_tag(mm).unwrap();
        let raw_mm_tags = RawModTags::new(mm, &vec![255u16; 8], true);
        let modbase_info =
            ModBaseInfo::new(&mm_tag_infos, &raw_mm_tags, dna.as_bytes());
        assert!(modbase_info.is_err());
        // todo add test for specific error
        // if let Err(e) = modbase_info {
        //     dbg!(e);
        // }
    }

    #[test]
    fn test_generic_mm_tags_mixed_modes() {
        let dna = "CATCACA";
        let mm = "N+b?,0,1;C+m.,0,1;";
        let mm_tag_infos = MmTagInfo::parse_mm_tag(mm).unwrap();
        let raw_mm_tags = RawModTags::new(mm, &vec![200u16, 255, 50, 0], true);
        let mut record = bam::Record::new();
        record.set(b"test", None, dna.as_bytes(), &vec![255; dna.len()]);
        let _modbase_info =
            ModBaseInfo::new(&mm_tag_infos, &raw_mm_tags, dna.as_bytes())
                .unwrap();
        // let c_probs = modbase_info.pos_seq_base_mod_probs.get(&'C').unwrap();
        // dbg!(&c_probs.pos_to_base_mod_probs);
    }

    #[test]
    fn test_generic_mm_tags_inferred_conflict() {
        let dna = "CATCACA";
        let mm = "C+mh.,0,1;C+h.,0";
        let mm_tag_infos = MmTagInfo::parse_mm_tag(mm).unwrap();
        let raw_mm_tags = RawModTags::new(mm, &vec![200, 0, 0, 200, 25], true);
        let mut record = bam::Record::new();
        record.set(b"test", None, dna.as_bytes(), &vec![255; dna.len()]);
        let parse_result =
            ModBaseInfo::new(&mm_tag_infos, &raw_mm_tags, dna.as_bytes());

        match parse_result {
            Err(MkError::Conflict(ConflictError::ExplicitConflictInferred)) => {
                Ok(())
            }
            _ => Err("should fail"),
        }
        .unwrap();
    }
}
