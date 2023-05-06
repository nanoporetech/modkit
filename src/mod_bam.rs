use crate::errs::{InputError, RunError};
use crate::mod_base_code::{DnaBase, ModCode};
use crate::util;
use crate::util::{get_tag, Strand};
use indexmap::{indexset, IndexSet};
use std::cmp::Ordering;

use log::debug;
use rust_htslib::bam;
use rust_htslib::bam::record::Aux;
use std::collections::{HashMap, HashSet};

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
    // todo(arand) simplify to just use a hashmap
    pub(crate) mod_codes: IndexSet<char>,
    probs: Vec<f32>,
    // skip_mode: SkipMode,
    // strand: Strand,
}

impl BaseModProbs {
    pub(crate) fn new(mod_code: char, prob: f32) -> Self {
        let mod_codes = indexset! {mod_code};
        let probs = vec![prob];
        Self { mod_codes, probs }
    }

    pub(crate) fn insert_base_mod_prob(&mut self, mod_code: char, prob: f32) {
        if let Some(idx) = self.mod_codes.get_index_of(&mod_code) {
            self.probs[idx] += prob;
        } else {
            self.mod_codes.insert(mod_code);
            self.probs.push(prob);
        }
    }

    pub fn base_mod_call(&self) -> anyhow::Result<BaseModCall> {
        let canonical_prob = self.canonical_prob();
        // todo(arand) use iterprobs here
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
        1f32 - self.probs.iter().sum::<f32>()
    }

    pub fn iter_probs(&self) -> impl Iterator<Item = (&char, &f32)> {
        self.mod_codes.iter().zip(self.probs.iter())
    }

    pub fn iter_mut_probs(&mut self) -> impl Iterator<Item = &mut f32> {
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

                let mut mod_codes = IndexSet::new();
                let mut probs = Vec::new();
                for (mod_code, mod_prob) in marginal_collapsed_prob {
                    let collapsed_prob =
                        mod_prob / total_marginal_collapsed_prob;
                    mod_codes.insert(*mod_code);
                    probs.push(collapsed_prob)
                }

                Self { mod_codes, probs }
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
                    .collect::<Vec<_>>();
                let n_other_mods = other_mods.len() as f32 + 1f32; // plus 1 for the canonical base
                let prob_to_redistribute = marginal_prob / n_other_mods;

                let mut check_total = 0f32;
                let mut mod_codes = IndexSet::new();
                let mut probs = Vec::new();
                for (mod_code, prob) in other_mods {
                    let new_prob = prob + prob_to_redistribute;
                    check_total += new_prob;
                    mod_codes.insert(*mod_code);
                    probs.push(new_prob);
                }
                if check_total - 100f32 > 0.00001 {
                    debug!(
                        "total probability {check_total} did not re-normalize"
                    )
                }

                Self { mod_codes, probs }
            }
            CollapseMethod::Convert { from, to } => {
                let mut probs = Vec::new();
                let mut mod_codes = IndexSet::new();

                let mut converted_prob = 0f32;
                for (raw_mod_code, prob) in self.iter_probs() {
                    // if we're converting from, add to the accumulator
                    if from.contains(&raw_mod_code) {
                        converted_prob += prob;
                    } else {
                        // keep track as before
                        probs.push(*prob);
                        mod_codes.insert(*raw_mod_code);
                    }
                }

                let mut new_base_mod_probs = Self { probs, mod_codes };

                if converted_prob > 0f32 {
                    new_base_mod_probs
                        .insert_base_mod_prob(*to, converted_prob);
                }

                new_base_mod_probs
            }
        }
    }

    fn combine(&mut self, other: Self) {
        for (mod_code, prob) in
            other.mod_codes.into_iter().zip(other.probs.into_iter())
        {
            if let Some(idx) = self.mod_codes.get_index_of(&mod_code) {
                self.probs[idx] += prob
            } else {
                self.mod_codes.insert(mod_code);
                self.probs.push(prob);
            }
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
            let qual = *q as f32;
            *q = (qual + 0.5f32) / 256f32;
        }
    });
}

#[derive(Debug, Eq, PartialEq)]
pub(crate) struct BaseModPositions {
    pub(crate) canonical_base: char,
    mode: SkipMode,
    strand: Strand,
    mod_base_codes: Vec<char>,
    delta_list: Vec<u32>,
}

impl BaseModPositions {
    pub(crate) fn parse(mod_positions: &str) -> Result<Self, InputError> {
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

        while let Some(c) = header.next() {
            match c {
                '?' | '.' => {
                    mode = Some(SkipMode::parse(c).unwrap());
                }
                _ => mod_base_codes.push(c),
            }
        }
        // default to the "old version"
        let mode = mode.unwrap_or(SkipMode::ImplicitProbModified);

        // taking the liberty to think that a read wouldn't be larger
        // than 2**32 - 1 bases long
        let delta_list = parts
            .into_iter()
            .map(|raw_pos| raw_pos.replace(";", "").parse::<u32>())
            .collect::<Result<Vec<u32>, _>>()
            .map_err(|e| {
                InputError::new(&format!(
                    "failed to parse position list, {}",
                    e.to_string()
                ))
            })?;

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
#[derive(PartialEq, Debug)]
pub struct SeqPosBaseModProbs {
    /// The `.` or `?` or implied mode, see `SkipMode`.
    pub skip_mode: SkipMode,
    /// Mapping of _forward_ sequence position to the predicted base
    /// modification probabilities for that position.
    pub pos_to_base_mod_probs: HashMap<usize, BaseModProbs>,
}

impl SeqPosBaseModProbs {
    fn new(
        pos_to_base_mod_probs: HashMap<usize, BaseModProbs>,
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
            pos_to_base_mod_probs: HashMap::new(),
        }
    }
}

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

    let mut positions_to_probs = HashMap::<usize, BaseModProbs>::new();
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
            mod_code_to_position.into_iter()
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
        if self.pos_seq_base_mod_probs.is_empty()
            && self.neg_seq_base_mod_probs.is_empty()
        {
            debug_assert!(self.converters.is_empty());
            true
        } else {
            false
        }
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

#[cfg(test)]
mod mod_bam_tests {
    use super::*;

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
        let mod_base_probs = BaseModProbs {
            mod_codes: indexset! {'h', 'm'},
            probs: vec![0.85, 0.10],
        };
        let collapsed = mod_base_probs
            .clone()
            .into_collapsed(&CollapseMethod::ReDistribute('h'));
        assert_eq!(collapsed.probs, vec![0.52500004]);
        assert_eq!(collapsed.mod_codes, indexset! {'m'});
        let collapsed = mod_base_probs
            .clone()
            .into_collapsed(&CollapseMethod::ReNormalize('h'));
        assert_eq!(collapsed.probs, vec![0.6666669]);
        assert_eq!(collapsed.mod_codes, indexset! {'m'});

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
        let mod_base_probs = BaseModProbs {
            mod_codes: indexset! {'h', 'm'},
            probs: vec![0.05273438, 0.03320312],
        };
        let collapsed =
            mod_base_probs.into_collapsed(&CollapseMethod::ReNormalize('h'));
        assert_eq!(collapsed.probs, vec![0.035051543]);
        assert_eq!(collapsed.mod_codes, indexset! {'m'});
    }

    #[test]
    fn test_mod_prob_collapse_dist_examples() {
        let mod_base_probs = BaseModProbs {
            mod_codes: indexset! {'h', 'm'},
            probs: vec![0.05273438, 0.03320312],
        };
        let collapsed =
            mod_base_probs.into_collapsed(&CollapseMethod::ReDistribute('h'));
        assert_eq!(collapsed.probs, vec![0.059570313]);
        assert_eq!(collapsed.mod_codes, indexset! {'m'});
    }

    #[test]
    fn test_mod_prob_convert() {
        let mod_base_probs = BaseModProbs {
            mod_codes: indexset! {'h', 'm'},
            probs: vec![0.10, 0.75],
        };
        let collapsed =
            mod_base_probs.into_collapsed(&CollapseMethod::Convert {
                from: HashSet::from(['h']),
                to: 'C',
            });
        assert_eq!(collapsed.probs, vec![0.75, 0.10]);
        assert_eq!(collapsed.mod_codes, indexset! {'m', 'C'});

        let mod_base_probs = BaseModProbs {
            mod_codes: indexset! {'h', 'm'},
            probs: vec![0.10, 0.75],
        };
        let collapsed =
            mod_base_probs.into_collapsed(&CollapseMethod::Convert {
                from: HashSet::from(['h', 'm']),
                to: 'C',
            });
        assert_eq!(collapsed.probs, vec![0.85]);
        assert_eq!(collapsed.mod_codes, indexset! {'C'});
    }

    #[test]
    fn test_mod_prob_convert_sums_prob() {
        let mod_base_probs = BaseModProbs {
            mod_codes: indexset! {'h', 'm'},
            probs: vec![0.10, 0.75],
        };
        let collapsed =
            mod_base_probs.into_collapsed(&CollapseMethod::Convert {
                from: HashSet::from(['h']),
                to: 'm',
            });
        assert_eq!(collapsed.probs, vec![0.85]);
        assert_eq!(collapsed.mod_codes, indexset! {'m'});
    }

    #[test]
    fn test_mod_prob_convert_noop() {
        let mod_base_probs = BaseModProbs {
            mod_codes: indexset! {'h', 'm'},
            probs: vec![0.10, 0.75],
        };
        let collapsed =
            mod_base_probs.into_collapsed(&CollapseMethod::Convert {
                from: HashSet::from(['a']),
                to: 'A',
            });
        assert_eq!(collapsed.probs, vec![0.10, 0.75]);
        assert_eq!(collapsed.mod_codes, indexset! {'h', 'm'});
    }

    #[test]
    fn test_mod_prob_combine() {
        let mut a = BaseModProbs {
            mod_codes: indexset! {'h', 'm'},
            probs: vec![0.05273438, 0.03320312],
        };
        let b = BaseModProbs {
            mod_codes: indexset! {'m'},
            probs: vec![0.03320312],
        };
        a.combine(b);
        assert_eq!(&a.probs, &[0.05273438, 0.06640624]);

        let mut a = BaseModProbs {
            mod_codes: indexset! {'m'},
            probs: vec![0.03320312],
        };

        let b = BaseModProbs {
            mod_codes: indexset! {'h'},
            probs: vec![0.05273438],
        };
        a.combine(b);
        assert_eq!(&a.probs, &[0.03320312, 0.05273438]);
        assert_eq!(&a.mod_codes, &indexset! {'m', 'h'});
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
        .collect::<HashMap<usize, BaseModProbs>>();

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
        // let strand = Strand::Positive;
        let positions_and_probs = vec![
            (5, BaseModProbs::new('m', 0.9)),
            (2, BaseModProbs::new('m', 0.1)),
            (8, BaseModProbs::new('m', 0.2)),
        ]
        .into_iter()
        .collect::<HashMap<usize, BaseModProbs>>();

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
                .map(|p| prob_to_qual(*p))
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
            assert_eq!(&base_mod_probs.probs, &[0.005859375, 0.005859375]);
            assert_eq!(&base_mod_probs.mod_codes, &indexset! { 'h', 'm' });
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
            assert_eq!(base_mod_probs.mod_codes, other.mod_codes);
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

        let c_expected = BaseModProbs {
            mod_codes: indexset! { 'h', 'm' },
            probs: vec![0.005859375, 0.39257813],
        };
        let a_expected = BaseModProbs {
            mod_codes: indexset! {'a'},
            probs: vec![0.7832031],
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
                let expected_modbase_probs = BaseModProbs {
                    probs: vec![0.39257813, 0.005859375],
                    mod_codes: indexset! { 'h', 'm'},
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
                let expected_modbase_probs = BaseModProbs {
                    probs: vec![0.5878906, 0.009765625],
                    mod_codes: indexset! { 'h', 'm'},
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
}
