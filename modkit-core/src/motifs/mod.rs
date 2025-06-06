use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap};
use std::fmt::{Display, Formatter};
use std::io::{BufRead, BufReader};
use std::ops::Range;
use std::path::PathBuf;
use std::sync::{Arc, RwLock};

use anyhow::{anyhow, bail, Context};
use bitvec::prelude::*;
use derive_new::new;
use indicatif::{MultiProgress, ParallelProgressIterator, ProgressIterator};
use itertools::{iproduct, Either, Itertools};
use log::{error, info, warn};
use nom::character::complete::multispace1;
use nom::IResult;
use prettytable as pt;
use prettytable::row;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use tracing::debug;

use crate::dmr::bedmethyl::BedMethylLine;
use crate::errs::MkError;
use crate::fasta::HtsLibFastaRecords;
use crate::mod_base_code::{DnaBase, ModCodeRepr, MOD_CODE_TO_DNA_BASE};
use crate::motifs::args::KnownMotifsArgs;
use crate::motifs::iupac::nt_bytes::BASES;
use crate::motifs::iupac::IupacBase;
use crate::parsing_utils::consume_string;
use crate::util::{get_subroutine_progress_bar, get_ticker, StrandRule};

mod args;
pub(crate) mod iupac;
pub mod motif_bed;
pub mod subcommand;
pub(super) mod util;

type KmerRef<'a> = &'a [u8];
type RawBase = u8;

#[derive(Debug)]
enum Action {
    Found,
    Refined,
    Discard,
}

impl Display for Action {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let l = match self {
            Action::Found => "found",
            Action::Discard => "discard",
            Action::Refined => "refined",
        };
        write!(f, "{l}")
    }
}

#[derive(Debug)]
enum Stage {
    Seeded,
    Seedless,
    Search,
}

impl Display for Stage {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let l = match self {
            Stage::Seeded => "seeded",
            Stage::Seedless => "seedless",
            Stage::Search => "search",
        };
        write!(f, "{l}")
    }
}

#[derive(Copy, Clone, Debug)]
enum SearchConfig {
    FullSearch,
    TopFrac {
        frac: f32,
        min_seeds: usize,
        max_seeds: usize,
    },
    BatchAndNarrow {
        frac: f32,
        min_seeds: usize,
        max_seeds: usize,
        max_iters: Option<usize>,
    },
    TimeoutAndNarrow {
        frac: f32,
        min_seeds: usize,
        max_seeds: usize,
        total_time: std::time::Duration,
        max_iters: Option<usize>,
    },
    Timeout {
        batch_size: usize,
        total_time: std::time::Duration,
    },
}

fn get_interval(focus_position: usize, context: &[usize; 2]) -> Range<usize> {
    // todo should be a result and/or need wrap this into a ContextBase struct
    // that can be validated on construction time
    (focus_position - context[0])..(focus_position + context[1] + 1)
}

#[derive(Debug, Clone, Eq, PartialEq, Hash)]
struct MultiSequence {
    mod_code: ModCodeRepr,
    seq: BTreeMap<i8, IupacBase>, // might need to make this a faster thing
}

impl MultiSequence {
    fn new(mod_code: ModCodeRepr, seq: BTreeMap<i8, IupacBase>) -> Self {
        Self { mod_code, seq }
    }

    fn from_iter<'a, IT: Iterator<Item = KmerRef<'a>>>(
        iter: IT,
        context_bases: [usize; 2],
        mod_code: ModCodeRepr,
    ) -> Self {
        // let mut raw_kmers = Vec::new();
        let kmer_len = context_bases.iter().sum::<usize>() + 1usize;
        let mut seq = BTreeMap::<i8, IupacBase>::default();
        for kmer in iter {
            assert_eq!(kmer.len(), kmer_len);
            // raw_kmers.push(kmer);
            for left_pos in 0..context_bases[0] {
                let motif_pos = left_pos as i8 - context_bases[0] as i8;
                if let Some(&nt) = kmer.get(left_pos) {
                    seq.entry(motif_pos)
                        .or_insert_with(|| IupacBase::from_base_unchecked(nt))
                        .add_mut(nt);
                }
            }
            for right_pos in
                (0..context_bases[1]).map(|i| i + context_bases[0] + 1)
            {
                let motif_pos = (right_pos - context_bases[0]) as i8;
                if let Some(&nt) = kmer.get(right_pos) {
                    seq.entry(motif_pos)
                        .or_insert_with(|| IupacBase::from_base_unchecked(nt))
                        .add_mut(nt);
                }
            }
        }
        let seq = seq
            .into_iter()
            .filter(|(_pos, base)| *base != IupacBase::N)
            .collect::<BTreeMap<i8, IupacBase>>();
        let this = Self { seq, mod_code };
        this
    }

    fn matches(&self, seq: &[u8], focus_position: usize) -> bool {
        let fp = focus_position as i8;
        self.seq.iter().all(|(&relative_pos, base)| {
            let offset = fp + relative_pos;
            assert!(
                offset >= 0,
                "offset should still be > 0 {offset}, {focus_position} \
                 {relative_pos}"
            );
            let offset = offset as usize;
            base.matches(seq[offset])
        })
    }

    fn clean(&mut self) {
        self.seq.retain(|_pos, base| match *base {
            IupacBase::Hole | IupacBase::N => false,
            _ => true,
        });
    }

    fn is_superset(&self, other: &Self) -> bool {
        if self.mod_code != other.mod_code {
            return false;
        }
        if self.seq == other.seq {
            return true;
        }
        if self.seq.is_empty() {
            return false;
        } else if other.seq.is_empty() {
            return true;
        }
        if self.seq.len() < other.seq.len() {
            return false;
        }
        assert!(self.seq.len() >= other.seq.len());

        let xs = self.seq.keys().collect::<FxHashSet<&i8>>();
        let ys = other.seq.keys().collect::<FxHashSet<&i8>>();
        if xs.is_superset(&ys) {
            let matches = self
                .seq
                .iter()
                .map(|(x, a)| {
                    other
                        .seq
                        .get(x)
                        .copied()
                        .map(|b| (*a, b))
                        .unwrap_or_else(|| (*a, IupacBase::N))
                })
                .all(|(a, b)| a.is_superset(&b));
            // original implementation left here for reference.. In this case a
            // motif such as GS[a]TC would be a superset of C[a]TG,
            // but that's not exactly true so I changed to the above
            // implementaton with no regression on H. pylori
            // let matches = self
            //     .seq
            //     .iter()
            //     .filter_map(|(x, a)| other.seq.get(x).map(|b| (a, b)))
            //     .all(|(a, b)| a.is_superset(b));
            matches
        } else {
            false
        }
    }

    #[inline]
    fn get_string_bookends(&self) -> (String, String) {
        let mut before_slots = {
            let size = self
                .seq
                .keys()
                .filter(|x| **x < 0i8)
                .map(|x| x.abs() as usize)
                .max();
            if let Some(s) = size {
                vec![IupacBase::N; s]
            } else {
                Vec::new()
            }
        };
        let mut after_slots = {
            let size = self
                .seq
                .keys()
                .filter(|x| **x > 0i8)
                .map(|x| x.abs() as usize)
                .max();
            if let Some(s) = size {
                vec![IupacBase::N; s]
            } else {
                Vec::new()
            }
        };
        for (pos, base) in self.seq.iter() {
            if *pos < 0i8 {
                assert!(!before_slots.is_empty());
                let offset = before_slots
                    .len()
                    .checked_sub(pos.abs() as usize)
                    .expect("should be within range");
                before_slots[offset] = *base
            } else {
                assert!(*pos > 0i8);
                let offset =
                    (*pos as usize).checked_sub(1).expect("should be >= 1");
                after_slots[offset] = *base;
            }
        }
        let before =
            before_slots.into_iter().map(|b| b.to_char()).collect::<String>();
        let after =
            after_slots.into_iter().map(|b| b.to_char()).collect::<String>();
        (before, after)
    }

    fn format_seq(&self, canonical_base: DnaBase) -> String {
        let (before, after) = self.get_string_bookends();
        let middle = canonical_base.char();
        format!("{before}{middle}{after}")
    }

    fn get_offset(&self) -> usize {
        match self.seq.keys().min() {
            Some(leftmost) if *leftmost < 0i8 => leftmost.abs() as usize,
            _ => 0,
        }
    }

    fn to_bits(&self, context_size: [usize; 2]) -> BitVec {
        let focus_position = context_size[0] as i8;
        let total_size = (context_size[0] + context_size[1])
            * iupac::iupac_offsets::ALPHABET_SIZE;
        let mut bv = bitvec![0; total_size];
        let iter = self
            .seq
            .iter()
            .map(|(motif_pos, base)| {
                if *motif_pos > 0i8 {
                    (*motif_pos - 1i8, base)
                } else {
                    (*motif_pos, base)
                }
            })
            .map(|(motif_pos, base)| {
                let idx = motif_pos + focus_position;
                assert!(idx >= 0);
                (idx as usize, base.to_offsets())
            });

        for (idx, offsets) in iter {
            let i = idx * iupac::iupac_offsets::ALPHABET_SIZE;
            for &offset in offsets {
                bv.set(i + offset, true);
            }
        }

        bv
    }

    fn edit_distance(&self, other: &Self, context_size: [usize; 2]) -> usize {
        (self.to_bits(context_size) ^ other.to_bits(context_size)).count_ones()
    }
}

impl Display for MultiSequence {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let (before, after) = self.get_string_bookends();
        let middle = self.mod_code.to_string();
        write!(f, "{before}[{middle}]{after}")
    }
}

struct KmerTable {
    focus_position: usize,
    counts: FxHashMap<Vec<u8>, FxHashMap<ModCodeRepr, u32>>,
}

impl KmerTable {
    fn new_empty(context_bases: &[u64; 2]) -> Self {
        let focus_position = context_bases[0] as usize;
        Self { focus_position, counts: FxHashMap::default() }
    }

    fn add(&mut self, kmer: &[u8], mod_code: ModCodeRepr) {
        // decided to use this api instead of .entry so that I can have a slice
        // in the method signature
        if let Some(counts_per_mod) = self.counts.get_mut(kmer) {
            *counts_per_mod.entry(mod_code).or_insert(0u32) += 1;
        } else {
            let kmer = kmer.to_vec();
            self.counts
                .entry(kmer)
                .or_insert(FxHashMap::default())
                .insert(mod_code, 1);
        }
    }

    fn sliced_kmers(
        &self,
        canonical_base: DnaBase,
        mod_code: &ModCodeRepr,
        mask: &FxHashSet<KmerRef>,
        context: &[usize; 2],
    ) -> FxHashMap<&[u8], u32> {
        let can_base = canonical_base.as_byte();
        self.counts
            .iter()
            .filter(|(kmer, _)| !mask.contains(kmer.as_slice()))
            .filter(|(kmer, _)| {
                let focus_base = kmer[self.focus_position];
                focus_base == can_base
            })
            .filter_map(|(kmer, counts)| {
                if let Some(count) = counts.get(mod_code) {
                    let interval = get_interval(self.focus_position, context);
                    let sliced = &kmer[interval];
                    Some((sliced, *count))
                } else {
                    None
                }
            })
            .fold(FxHashMap::default(), |mut agg, (kmer, count)| {
                *agg.entry(kmer).or_insert(0) += count;
                agg
            })
    }

    fn count_total_matches(
        &self,
        motif: &EnrichedMotif,
        focus_position: usize,
        canonical_base: RawBase,
    ) -> u64 {
        self.counts
            .par_iter()
            // todo a performance improvement would be to organize the contexts
            //  by their canonical base  ahead of time
            .filter(|(kmer, _)| kmer[focus_position] == canonical_base)
            .filter_map(|(kmer, counts)| {
                counts.get(&motif.multi_sequence.mod_code).and_then(|count| {
                    if motif
                        .multi_sequence
                        .matches(kmer.as_slice(), focus_position)
                    {
                        Some(*count as u64)
                    } else {
                        None
                    }
                })
            })
            .sum::<u64>()
    }

    fn count_total_not_matching(
        &self,
        motif: &EnrichedMotif,
        focus_position: usize,
        total_matching: u64,
    ) -> u64 {
        let total_potential_count = self
            .counts
            .par_iter()
            .filter(|(kmer, _)| {
                kmer[focus_position] == motif.canonical_base.as_byte()
            })
            .filter_map(|(_kmer, counts)| {
                counts.get(&motif.multi_sequence.mod_code).map(|x| *x as u64)
            })
            .sum::<u64>();
        assert!(
            total_potential_count >= total_matching,
            "total potential should be >= total matching"
        );
        total_potential_count.checked_sub(total_matching).unwrap()
    }

    fn count_matches(
        &self,
        motif: &EnrichedMotif,
        mask: &FxHashSet<KmerRef>,
        focus_position: usize,
        canonical_base: RawBase,
    ) -> u64 {
        self.counts
            .par_iter()
            .filter(|(kmer, _)| !mask.contains(kmer.as_slice()))
            // todo a performance improvement would be to organize the contexts
            // by their canonical base  ahead of time
            .filter(|(kmer, _)| kmer[focus_position] == canonical_base)
            .filter_map(|(kmer, counts)| {
                counts.get(&motif.multi_sequence.mod_code).and_then(|count| {
                    if motif
                        .multi_sequence
                        .matches(kmer.as_slice(), focus_position)
                    {
                        Some(*count as u64)
                    } else {
                        None
                    }
                })
            })
            .sum::<u64>()
    }

    fn get_matching_contexts(
        &self,
        motif: &EnrichedMotif,
        mask: &FxHashSet<KmerRef>,
        focus_position: usize,
        canonical_base: u8,
    ) -> Vec<(KmerRef, u32)> {
        self.counts
            .par_iter()
            .filter(|(kmer, _)| !mask.contains(kmer.as_slice()))
            // todo a performance improvement would be to organize the contexts
            // by their canonical base  ahead of time
            .filter(|(kmer, _)| kmer[focus_position] == canonical_base)
            .filter_map(|(kmer, counts)| {
                counts
                    .get(&motif.multi_sequence.mod_code)
                    .map(|count| (kmer.as_slice(), *count))
            })
            .filter(|(kmer, _)| {
                motif.multi_sequence.matches(kmer, focus_position)
            })
            .collect()
    }
}

#[derive(new, Default)]
struct KmerMask<'a> {
    high_mod_mask: FxHashSet<KmerRef<'a>>,
    low_mod_mask: FxHashSet<KmerRef<'a>>,
}

impl<'a> KmerMask<'a> {
    fn update_with_check(
        self,
        motifs: &[EnrichedMotif],
        mod_db: &'a KmerModificationDb,
        mod_code: ModCodeRepr,
        stage: Stage,
    ) -> Either<Self, Self> {
        let prev_high = self.high_mod_mask;
        let prev_low = self.low_mod_mask;
        let new_high = motifs
            .par_iter()
            .flat_map(|motif| {
                mod_db
                    .get_high_mod_contexts(motif, &prev_high)
                    .into_iter()
                    .map(|(kmer, _)| kmer)
                    .collect::<FxHashSet<KmerRef>>()
            })
            .collect::<FxHashSet<KmerRef>>();
        let new_low = motifs
            .par_iter()
            .flat_map(|motif| {
                mod_db
                    .get_low_mod_contexts(motif, &prev_low)
                    .into_iter()
                    .map(|(kmer, _)| kmer)
                    .collect::<FxHashSet<KmerRef>>()
            })
            .collect::<FxHashSet<KmerRef>>();

        let count_high_removed = new_high.len();
        let count_low_removed = new_low.len();
        // todo could add mod_code and stage here..
        debug!(
            mod_code = mod_code.to_string(),
            stage = stage.to_string(),
            "removing {} contexts from the high list",
            count_high_removed
        );
        debug!(
            mod_code = mod_code.to_string(),
            stage = stage.to_string(),
            "removing {} contexts from the low list",
            count_low_removed
        );
        let high_mod_mask = prev_high
            .into_iter()
            .chain(new_high)
            .collect::<FxHashSet<KmerRef>>();
        let low_mod_mask =
            prev_low.into_iter().chain(new_low).collect::<FxHashSet<KmerRef>>();
        let new_mask = Self::new(high_mod_mask, low_mod_mask);
        match count_high_removed.checked_add(count_low_removed) {
            Some(x) if x == 0 => Either::Left(new_mask),
            _ => Either::Right(new_mask),
        }
    }
}

#[derive(new)]
struct PositionBools {
    high_bools: FxHashMap<(usize, RawBase), BitVec<usize, Lsb0>>,
    low_bools: FxHashMap<(usize, RawBase), BitVec<usize, Lsb0>>,
    n_high: usize,
    n_low: usize,
    high_total: u32,
    low_total: u32,
}

impl PositionBools {
    #[inline]
    fn get_high_empty(&self) -> BitVec {
        bitvec![usize, Lsb0; 0; self.n_high]
    }

    #[inline]
    fn get_low_empty(&self) -> BitVec {
        bitvec![usize, Lsb0; 0; self.n_low]
    }

    fn get_counts(
        &self,
        kmer: &[RawBase],
        positions: &[usize],
    ) -> ContextsCounts {
        let high_matches = kmer
            .iter()
            .zip(positions.iter())
            .map(|(b, p)| {
                let k = (*p, *b);
                self.high_bools
                    .get(&k)
                    .map(|b| b.clone())
                    .unwrap_or(self.get_high_empty())
            })
            .reduce(|acc, next| acc & next)
            .unwrap()
            .count_ones() as u32;

        let low_matches = kmer
            .iter()
            .zip(positions.iter())
            .map(|(b, p)| {
                let k = (*p, *b);
                self.low_bools
                    .get(&k)
                    .map(|b| b.clone())
                    .unwrap_or(self.get_low_empty())
            })
            .reduce(|acc, next| acc & next)
            .unwrap()
            .count_ones() as u32;
        let high_not_matches = self
            .high_total
            .checked_sub(high_matches)
            .expect("high matches should not be more than total");
        let low_not_matches = self
            .low_total
            .checked_sub(low_matches)
            .expect("low matches should not be more than total");
        ContextsCounts::new(
            high_matches,
            high_not_matches,
            low_matches,
            low_not_matches,
        )
    }
}

#[derive(new)]
struct KmerSubset<'a> {
    high_mod_table: FxHashMap<KmerRef<'a>, u32>,
    low_mod_table: FxHashMap<KmerRef<'a>, u32>,
}

impl<'a> KmerSubset<'a> {
    #[inline]
    fn get_bools(
        &self,
        high_mod: bool,
        focus_base: usize,
    ) -> FxHashMap<(usize, u8), BitVec<usize, Lsb0>> {
        let table =
            if high_mod { &self.high_mod_table } else { &self.low_mod_table };

        table
            .par_iter()
            .fold(
                || FxHashMap::default(),
                |mut agg, (kmer, count)| {
                    for (idx, seq_base) in kmer
                        .iter()
                        .enumerate()
                        .filter(|(idx, _)| *idx != focus_base)
                    {
                        for nt in BASES {
                            let key = (idx, nt);
                            let val = nt == *seq_base;
                            let entry =
                                agg.entry(key)
                                    .or_insert(BitVec::<usize, Lsb0>::new());
                            for _ in 0..*count {
                                entry.push(val);
                            }
                        }
                    }
                    agg
                },
            )
            .reduce(
                || FxHashMap::default(),
                |mut a, mut b| {
                    for (k, bs) in b.iter_mut() {
                        a.entry(*k).or_insert(BitVec::new()).append(bs);
                    }
                    a
                },
            )
    }

    fn get_position_bools(&self, focus_position: usize) -> PositionBools {
        let now = std::time::Instant::now();
        let high_bools = self.get_bools(true, focus_position);
        let low_bools = self.get_bools(false, focus_position);

        debug_assert_eq!(
            high_bools.values().map(|bs| bs.len()).unique().count(),
            1
        );
        debug_assert_eq!(
            low_bools.values().map(|bs| bs.len()).unique().count(),
            1
        );

        let n_high = high_bools.values().next().map(|b| b.len()).unwrap_or(0);
        let n_low = low_bools.values().next().map(|b| b.len()).unwrap_or(0);

        let took = now.elapsed().as_micros();
        debug!("creating bool tables took {took} ms");
        let high_total = self.high_mod_table.values().sum::<u32>();
        let low_total = self.low_mod_table.values().sum::<u32>();

        PositionBools::new(
            high_bools, low_bools, n_high, n_low, high_total, low_total,
        )
    }

    #[inline]
    fn get_mod_contexts(
        &self,
        motif: &EnrichedMotif,
        high_mod: bool,
        focus_position: usize,
    ) -> Vec<(KmerRef, u32)> {
        let table =
            if high_mod { &self.high_mod_table } else { &self.low_mod_table };

        table
            .par_iter()
            .filter_map(|(kmer, count)| {
                if motif.multi_sequence.matches(kmer, focus_position) {
                    Some((*kmer, *count))
                } else {
                    None
                }
            })
            .collect()
    }

    fn get_high_mod_contexts(
        &self,
        motif: &EnrichedMotif,
        focus_position: usize,
    ) -> Vec<(KmerRef, u32)> {
        self.get_mod_contexts(motif, true, focus_position)
    }

    fn get_low_mod_contexts(
        &self,
        motif: &EnrichedMotif,
        focus_position: usize,
    ) -> Vec<(KmerRef, u32)> {
        self.get_mod_contexts(motif, false, focus_position)
    }
}

#[derive(new)]
struct ContextsCounts {
    high_match: u32,
    hig_not_match: u32,
    low_match: u32,
    low_not_match: u32,
}

impl ContextsCounts {
    fn log_odds(&self) -> f32 {
        util::log_odds(
            self.low_match,
            self.low_not_match,
            self.high_match,
            self.hig_not_match,
        )
    }
}

struct KmerModificationDb {
    context_bases: [usize; 2],
    low_mod_table: KmerTable,
    high_mod_table: KmerTable,
    mid_mod_table: KmerTable,
}

impl KmerModificationDb {
    fn get_fraction_modified_cached(
        &self,
        motif: &EnrichedMotif,
        cache: &RwLock<FxHashMap<String, (f32, u64)>>,
    ) -> f32 {
        let cache_ = cache.read().expect("cache poisoned");
        let key = motif.to_string();
        if let Some((f, _)) = cache_.get(&key) {
            *f
        } else {
            drop(cache_);
            let (high_counts, low_counts) = self.get_total_mod_counts(motif);
            let frac_mod =
                high_counts as f32 / (high_counts + low_counts) as f32;
            let mut w_cache = cache.write().expect("cache poisoned for writes");
            w_cache.insert(key, (frac_mod, high_counts));
            frac_mod
        }
    }

    fn get_inferred_mod_code_associations(
        &self,
        // todo allow overrides/definition of associations
        force_specification: bool,
    ) -> anyhow::Result<HashMap<ModCodeRepr, DnaBase>> {
        let context_iter = self
            .high_mod_table
            .counts
            .iter()
            .chain(self.mid_mod_table.counts.iter())
            .chain(self.low_mod_table.counts.iter());
        let mut focus_base_counter =
            FxHashMap::<ModCodeRepr, FxHashMap<DnaBase, usize>>::default();
        for (context, counts) in context_iter {
            let canonical_base =
                DnaBase::parse(context[self.focus_position()] as char)?;
            for mod_code in counts.keys() {
                *focus_base_counter
                    .entry(*mod_code)
                    .or_insert(FxHashMap::default())
                    .entry(canonical_base)
                    .or_insert(0) += 1;
            }
        }

        let associations = focus_base_counter
            .into_iter()
            .map(|(mod_code, counts)| {
                let inferred_primary_base = counts
                    .into_iter()
                    .max_by(|(_, a), (_, b)| a.cmp(b))
                    .map(|(b, _)| b)
                    .expect("should not be empty");
                (mod_code, inferred_primary_base)
            })
            .collect::<HashMap<ModCodeRepr, DnaBase>>();

        for (mod_code, prim_base) in associations.iter() {
            debug!(
                "inferred {prim_base:?} associated with modification code \
                 {mod_code}"
            );
            if let Some(expected_prim_base) = MOD_CODE_TO_DNA_BASE.get(mod_code)
            {
                if expected_prim_base != prim_base && force_specification {
                    bail!(
                        "modification code {mod_code} should be associated \
                         with {expected_prim_base:?}, use \
                         --force-override-spec to override."
                    )
                } else if expected_prim_base != prim_base {
                    warn!(
                        "modification code {mod_code} is normally associated \
                         with {expected_prim_base:?}, inferred to be \
                         associated with {prim_base:?}."
                    )
                }
            }
        }

        Ok(associations)
    }

    fn get_counts_and_frac_cached(
        &self,
        motif: &EnrichedMotif,
        cache: &RwLock<FxHashMap<String, (f32, u64)>>,
    ) -> (f32, u64) {
        let cache_ = cache.read().expect("cache poisoned");
        let key = motif.to_string();
        if let Some((f, c)) = cache_.get(&key) {
            (*f, *c)
        } else {
            drop(cache_);
            let (high_counts, low_counts) = self.get_total_mod_counts(motif);
            let frac_mod =
                high_counts as f32 / (high_counts + low_counts) as f32;
            let mut w_cache = cache.write().expect("cache poisoned for writes");
            w_cache.insert(key, (frac_mod, high_counts));
            (frac_mod, high_counts)
        }
    }

    fn focus_position(&self) -> usize {
        self.context_bases[0]
    }

    fn get_total_mod_counts(&self, motif: &EnrichedMotif) -> (u64, u64) {
        let low_counts = self.low_mod_table.count_total_matches(
            motif,
            self.context_bases[0],
            motif.canonical_base.as_byte(),
        );
        let high_counts = self.high_mod_table.count_total_matches(
            motif,
            self.context_bases[0],
            motif.canonical_base.as_byte(),
        );
        (high_counts, low_counts)
    }

    fn get_total_not_matching(
        &self,
        motif: &EnrichedMotif,
        n_high_matching: u64,
        n_low_matching: u64,
    ) -> (u64, u64) {
        let low_not_matching = self.low_mod_table.count_total_not_matching(
            motif,
            self.context_bases[0],
            n_low_matching,
        );
        let high_not_matching = self.high_mod_table.count_total_not_matching(
            motif,
            self.context_bases[0],
            n_high_matching,
        );
        (high_not_matching, low_not_matching)
    }

    fn get_mid_counts(&self, motif: &EnrichedMotif) -> u64 {
        self.mid_mod_table.count_total_matches(
            motif,
            self.context_bases[0],
            motif.canonical_base.as_byte(),
        )
    }

    fn get_mod_counts(
        &self,
        motif: &EnrichedMotif,
        mask: &KmerMask,
    ) -> (u64, u64) {
        let low_counts = self.low_mod_table.count_matches(
            motif,
            &mask.low_mod_mask,
            self.context_bases[0],
            motif.canonical_base.as_byte(),
        );
        let high_counts = self.high_mod_table.count_matches(
            motif,
            &mask.high_mod_mask,
            self.context_bases[0],
            motif.canonical_base.as_byte(),
        );
        (high_counts, low_counts)
    }

    fn get_high_mod_contexts(
        &self,
        motif: &EnrichedMotif,
        mask: &FxHashSet<KmerRef>,
    ) -> Vec<(KmerRef, u32)> {
        self.high_mod_table.get_matching_contexts(
            motif,
            mask,
            self.context_bases[0],
            motif.canonical_base.as_byte(),
        )
    }

    fn get_low_mod_contexts(
        &self,
        motif: &EnrichedMotif,
        mask: &FxHashSet<KmerRef>,
    ) -> Vec<(KmerRef, u32)> {
        self.low_mod_table.get_matching_contexts(
            motif,
            mask,
            self.context_bases[0],
            motif.canonical_base.as_byte(),
        )
    }

    fn count_mod_contexts(
        &self,
        mod_code: &ModCodeRepr,
        kmer_mask: &FxHashSet<KmerRef>,
        high_mod: bool,
    ) -> u64 {
        let table =
            if high_mod { &self.high_mod_table } else { &self.low_mod_table };

        table
            .counts
            .iter()
            .filter_map(|(kmer, counts)| {
                if !kmer_mask.contains(kmer.as_slice()) {
                    Some(*counts.get(&mod_code).unwrap_or(&0) as u64)
                } else {
                    None
                }
            })
            .sum::<u64>()
    }

    fn count_high_mod_contexts(
        &self,
        mod_code: &ModCodeRepr,
        kmer_mask: &KmerMask,
    ) -> u64 {
        self.count_mod_contexts(mod_code, &kmer_mask.high_mod_mask, true)
    }

    fn count_low_mod_contexts(
        &self,
        mod_code: &ModCodeRepr,
        kmer_mask: &KmerMask,
    ) -> u64 {
        self.count_mod_contexts(mod_code, &kmer_mask.low_mod_mask, false)
    }

    fn get_kmer_subset<'a>(
        &'a self,
        canonical_base: DnaBase,
        kmer_mask: &KmerMask,
        mod_code: ModCodeRepr,
    ) -> KmerSubset<'a> {
        let focus_position = self.context_bases[0];
        let canonical_base = canonical_base.as_byte();

        let get_table = |kmer_table: &'a KmerTable,
                         mask: &FxHashSet<KmerRef>|
         -> FxHashMap<KmerRef<'a>, u32> {
            kmer_table
                .counts
                .par_iter()
                // only contexts that have correct canonical base and NOT in the
                // mask
                .filter(|(kmer, _)| {
                    kmer[focus_position] == canonical_base
                        && !mask.contains(kmer.as_slice())
                })
                // filter to only contexts for which we have counts for this
                // modification code
                .filter_map(|(kmer, counts)| {
                    counts.get(&mod_code).map(|count| (kmer.as_slice(), *count))
                })
                .collect::<FxHashMap<KmerRef, u32>>()
        };

        let high_mod_table =
            get_table(&self.high_mod_table, &kmer_mask.high_mod_mask);
        let low_mod_table =
            get_table(&self.low_mod_table, &kmer_mask.low_mod_mask);

        KmerSubset { high_mod_table, low_mod_table }
    }

    fn get_enriched_motif_data(
        &self,
        motif: &EnrichedMotif,
    ) -> EnrichedMotifData {
        let (total_high_count, total_low_count) =
            self.get_total_mod_counts(&motif);
        let total_mid_count = self.get_mid_counts(&motif);
        let (total_high_not_match, total_low_not_match) = self
            .get_total_not_matching(&motif, total_high_count, total_low_count);
        EnrichedMotifData::new(
            motif.clone(),
            total_high_count,
            total_low_count,
            total_mid_count,
            total_high_not_match,
            total_low_not_match,
        )
    }
}

fn load_references_from_fasta(
    reference_fasta: &PathBuf,
    mpb: &MultiProgress,
) -> anyhow::Result<HashMap<String, Vec<u8>>> {
    info!("loading references from {:?}", reference_fasta);
    let pb = mpb.add(get_ticker());
    pb.set_message("sequences read");
    let reader = HtsLibFastaRecords::from_file(&reference_fasta)?;

    let (contigs, n_fails) = reader.into_iter().fold(
        (HashMap::new(), 0usize),
        |(mut agg, fails), record| match record {
            Ok((read_id, record_seq)) => {
                let seq = record_seq
                    .chars()
                    .map(|nt| nt.to_ascii_uppercase() as u8)
                    .collect::<Vec<u8>>();
                agg.insert(read_id, seq);
                pb.inc(1);
                (agg, fails)
            }
            Err(_) => (agg, fails + 1),
        },
    );

    if n_fails > 0 {
        info!("failed to load {n_fails} record(s)");
    }

    if contigs.is_empty() {
        bail!("failed to read any reference sequences");
    } else {
        pb.finish_and_clear();
        info!("loaded {} sequence(s)", contigs.len());
        Ok(contigs)
    }
}

fn parse_raw_known_motifs(
    raw_motifs: &[String],
    context: [usize; 2],
    mod_code_lookup: &HashMap<ModCodeRepr, DnaBase>,
) -> anyhow::Result<Vec<EnrichedMotif>> {
    let all_motifs = raw_motifs
        .chunks(3)
        .map(|parts| {
            EnrichedMotif::new_from_parts(
                &parts[0],
                &parts[2],
                &parts[1],
                context,
                mod_code_lookup,
            )
        })
        .collect::<anyhow::Result<Vec<EnrichedMotif>>>();

    all_motifs.map(|xs| xs.into_iter().unique().collect::<Vec<_>>())
}

fn parse_motifs_from_table(
    table_fp: &PathBuf,
    context: [usize; 2],
    mod_code_lookup: &HashMap<ModCodeRepr, DnaBase>,
) -> anyhow::Result<Vec<EnrichedMotif>> {
    fn parse_raw_motif_parts(
        l: &str,
    ) -> IResult<&str, (String, String, String)> {
        let (rest, raw_mod_code) = consume_string(l)?;
        let (rest, raw_motif_seq) =
            multispace1(rest).and_then(|(r, _)| consume_string(r))?;
        let (rest, raw_offset) =
            multispace1(rest).and_then(|(r, _)| consume_string(r))?;
        Ok((rest, (raw_mod_code, raw_motif_seq, raw_offset)))
    }

    let reader = BufReader::new(std::fs::File::open(table_fp)?);

    reader
        .lines()
        .skip_while(|r| {
            r.as_ref().map(|l| l.starts_with("mod_code")).unwrap_or(true)
        })
        .filter_map(|r| match r {
            Ok(l) => Some(l),
            Err(e) => {
                debug!("failed to read, {e}");
                None
            }
        })
        .map(|l| {
            parse_raw_motif_parts(&l)
                .map_err(|e| anyhow!("failed to parse line {l}, {e}"))
                .and_then(|(_, (raw_mod_code, raw_motif_seq, raw_offset))| {
                    EnrichedMotif::new_from_parts(
                        &raw_motif_seq,
                        &raw_mod_code,
                        &raw_offset,
                        context,
                        mod_code_lookup,
                    )
                })
        })
        .collect::<anyhow::Result<Vec<EnrichedMotif>>>()
}

fn load_bedmethyl_and_references(
    reference_fasta_fp: &PathBuf,
    bedmethyl_fp: &PathBuf,
    contig: Option<String>,
    min_coverage: u64,
    context_bases: [u64; 2],
    low_modification_threshold: f32,
    high_modification_threshold: f32,
    multi_progress: &MultiProgress,
    io_threads: usize,
    thread_pool: &rayon::ThreadPool,
) -> anyhow::Result<KmerModificationDb> {
    let reference_sequences =
        load_references_from_fasta(reference_fasta_fp, multi_progress)?;
    for x in context_bases {
        if x > 127u64 {
            bail!("context cannot be larger than 127x2 (255) bases")
        }
    }

    thread_pool.install(|| {
        load_bedmethyl(
            bedmethyl_fp,
            contig,
            min_coverage,
            context_bases,
            low_modification_threshold,
            high_modification_threshold,
            Arc::new(reference_sequences),
            &multi_progress,
            io_threads,
        )
    })
}

fn load_bedmethyl(
    bedmethyl_fp: &PathBuf,
    contig: Option<String>,
    min_coverage: u64,
    context_bases: [u64; 2],
    low_modification_threshold: f32,
    high_modification_threshold: f32,
    reference_sequences: Arc<HashMap<String, Vec<u8>>>,
    multi_progress: &MultiProgress,
    threads: usize,
) -> anyhow::Result<KmerModificationDb> {
    enum ModLevel {
        High,
        Middle,
        Low,
    }

    let file_fh = std::fs::File::open(bedmethyl_fp)
        .context(format!("failed to open bedMethyl at {bedmethyl_fp:?}"))?;

    let tbx_reader = crate::tabix::HtsTabixHandler::<BedMethylLine>::from_path(
        &bedmethyl_fp,
    );
    let use_tabix = tbx_reader.is_ok();
    if use_tabix {
        info!("using tabix/bgzf reader");
    }
    if !use_tabix && bedmethyl_fp.ends_with(".gz") {
        warn!(
            "failed to use tabix/bgzf reader, but file indicates compressed \
             input"
        );
    }
    if let (Err(e), Some(_)) = (&tbx_reader, contig.as_ref()) {
        bail!(
            "failed to use tabix index (required when --contig provided). \
             Error was {e}."
        )
    }
    if let Some(ctg) = contig.as_ref() {
        if !reference_sequences.contains_key(ctg.as_str()) {
            bail!("contig {ctg} not found in reference")
        }
    }

    let pb = multi_progress.add(get_ticker());
    pb.set_message("parsing bedMethyl records");

    let (snd, rcv) = crossbeam::channel::bounded(1000);
    let ref_seqs_reader = reference_sequences.clone();
    std::thread::spawn(move || {
        if let Ok(handler) = tbx_reader {
            let contigs = if let Some(ctg) = contig {
                let contig_length = ref_seqs_reader
                    .get(ctg.as_str())
                    .map(|x| x.len() as u64)
                    .unwrap(); // safe because of above check
                vec![(ctg.to_owned(), contig_length)]
            } else {
                ref_seqs_reader
                    .iter()
                    .map(|(name, s)| (name.to_owned(), s.len() as u64))
                    .collect::<Vec<_>>()
            };
            contigs.par_iter().for_each(|(name, len)| {
                match handler.read_bedmethyl(name, &(0..(*len)), threads) {
                    Ok(records) => {
                        records.into_iter().for_each(|r| match snd.send(r) {
                            Ok(_) => {}
                            Err(e) => {
                                error!("failed to send on channel, {e}")
                            }
                        });
                    }
                    Err(e) => {
                        debug!(
                            "failed to fetch bedMethyl records for {name}, {e}"
                        );
                    }
                }
            });
        } else {
            BufReader::new(file_fh)
                .lines()
                .par_bridge()
                .map(|r| {
                    r.map_err(|_| MkError::InvalidIO)
                        .and_then(|s| BedMethylLine::parse(&s))
                })
                .for_each(|r| match snd.send(r) {
                    Ok(_) => {}
                    Err(e) => {
                        error!("failed to send on channel, {e}");
                    }
                })
        }
    });

    let (bedmethyl_records, n_fails, n_discard) =
        rcv.into_iter().progress_with(pb).fold(
            (Vec::new(), 0usize, 0usize),
            |(mut records, mut fails, mut n_discard), next| {
                match next {
                    Ok(record) => {
                        if record.valid_coverage >= min_coverage {
                            if record.frac_modified()
                                <= low_modification_threshold
                            {
                                records.push((ModLevel::Low, record));
                            } else if record.frac_modified()
                                > high_modification_threshold
                            {
                                records.push((ModLevel::High, record));
                            } else {
                                records.push((ModLevel::Middle, record));
                            }
                        } else {
                            n_discard += 1;
                        }
                    }
                    Err(_) => fails += 1,
                };
                (records, fails, n_discard)
            },
        );

    if bedmethyl_records.is_empty() {
        bail!("failed to parse any bedmethyl records")
    }
    if n_fails > 0 {
        bail!("failed to parse {n_fails} bedmethyl records")
    }

    info!(
        "parsed {} bedmethyl records, discarded {n_discard} for low coverage, \
         {} total",
        bedmethyl_records.len(),
        bedmethyl_records.len() + n_discard
    );

    let mut low_mod_table = KmerTable::new_empty(&context_bases);
    let mut high_mod_table = KmerTable::new_empty(&context_bases);
    let mut mid_mod_table = KmerTable::new_empty(&context_bases);
    let mut discarded_contexts = 0usize;
    let mut mod_codes = FxHashSet::default();

    let check_kmer = |kmer: KmerRef| -> bool {
        kmer.iter().all(|b| match *b {
            iupac::nt_bytes::A
            | iupac::nt_bytes::C
            | iupac::nt_bytes::G
            | iupac::nt_bytes::T => true,
            _ => false,
        })
    };

    let pb = multi_progress
        .add(get_subroutine_progress_bar(bedmethyl_records.len()));
    pb.set_message("preparing contexts");
    for (mod_level, record) in bedmethyl_records {
        if let Some(seq) = reference_sequences.get(&record.chrom) {
            match record.strand {
                StrandRule::Both | StrandRule::Positive => {
                    let before_start =
                        record.start().checked_sub(context_bases[0]).is_none();
                    let after_end = ((record.start() + context_bases[1] + 1)
                        as usize)
                        > seq.len();
                    if before_start || after_end {
                        continue;
                    }
                    let start = record.start() - context_bases[0];
                    let end = record.start() + context_bases[1] + 1;
                    let interval = (start as usize)..(end as usize);
                    let kmer = &seq[interval];
                    if check_kmer(kmer) {
                        let table = match mod_level {
                            ModLevel::High => &mut high_mod_table,
                            ModLevel::Middle => &mut mid_mod_table,
                            ModLevel::Low => &mut low_mod_table,
                        };
                        table.add(kmer, record.raw_mod_code);
                        mod_codes.insert(record.raw_mod_code);
                    } else {
                        discarded_contexts += 1;
                    }
                }
                StrandRule::Negative => {
                    let before_start =
                        record.start().checked_sub(context_bases[1]).is_none();
                    let after_end = ((record.start() + context_bases[0] + 1)
                        as usize)
                        > seq.len();
                    if before_start || after_end {
                        continue;
                    }
                    let start = record.start() - context_bases[1];
                    let end = record.start() + context_bases[0] + 1;
                    let interval = (start as usize)..(end as usize);
                    let kmer = &seq[interval];
                    let kmer = bio::alphabets::dna::revcomp(kmer);
                    if check_kmer(&kmer) {
                        let table = match mod_level {
                            ModLevel::High => &mut high_mod_table,
                            ModLevel::Middle => &mut mid_mod_table,
                            ModLevel::Low => &mut low_mod_table,
                        };
                        table.add(&kmer, record.raw_mod_code);
                        mod_codes.insert(record.raw_mod_code);
                    } else {
                        discarded_contexts += 1;
                    }
                }
            }
        }
        pb.inc(1);
    }

    pb.finish_and_clear();

    info!(
        "loaded {} low-frequency, {} middle-frequency and {} high-frequency \
         contexts, {} modification codes, discarded {discarded_contexts} \
         contexts.",
        low_mod_table.counts.len(),
        mid_mod_table.counts.len(),
        high_mod_table.counts.len(),
        mod_codes.len(),
    );

    let context_bases = [context_bases[0] as usize, context_bases[1] as usize];
    Ok(KmerModificationDb {
        low_mod_table,
        high_mod_table,
        mid_mod_table,
        context_bases,
    })
}

#[derive(Debug, Eq, PartialEq, PartialOrd, Ord)]
pub(super) enum MotifRelationship {
    Equal,
    Subset,
    Superset,
    Disjoint { edit_distance: usize },
}

impl Display for MotifRelationship {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            MotifRelationship::Equal => "Equal",
            MotifRelationship::Subset => "Subset",
            MotifRelationship::Superset => "Superset",
            MotifRelationship::Disjoint { .. } => "Disjoint",
        };
        write!(f, "{}", s)
    }
}

#[derive(new)]
struct EnrichedMotifData {
    motif: EnrichedMotif,
    total_high_count: u64,
    total_low_count: u64,
    total_mid_count: u64,
    total_high_not_matching: u64,
    total_low_not_matching: u64,
}

impl EnrichedMotifData {
    fn frac_modified(&self) -> f32 {
        self.total_high_count as f32
            / (self.total_high_count + self.total_low_count) as f32
    }

    fn enough_sites(&self, min_sites: u64) -> bool {
        self.total_high_count >= min_sites || self.total_low_count >= min_sites
    }

    fn log_odds(&self) -> f32 {
        util::log_odds(
            self.total_low_count,
            self.total_low_not_matching,
            self.total_high_count,
            self.total_high_not_matching,
        )
    }
}

#[derive(new, Debug, Clone, Eq, PartialEq, Hash)]
struct EnrichedMotif {
    canonical_base: DnaBase,
    multi_sequence: MultiSequence,
}

impl EnrichedMotif {
    fn new_empty(mod_code: ModCodeRepr, canonical_base: DnaBase) -> Self {
        let multi_sequence = MultiSequence::new(mod_code, BTreeMap::new());
        Self { canonical_base, multi_sequence }
    }

    fn format_seq(&self) -> String {
        self.multi_sequence.format_seq(self.canonical_base)
    }

    fn new_from_parts(
        raw_seq: &str,
        raw_mod_code: &str,
        raw_offset: &str,
        context_size: [usize; 2],
        mod_code_lookup: &HashMap<ModCodeRepr, DnaBase>,
    ) -> anyhow::Result<Self> {
        let offset = raw_offset.parse::<usize>().context(format!(
            "failed to parse offset {raw_offset}, should be integer >= 0"
        ))?;
        if offset >= raw_seq.len() {
            bail!(
                "offset ({offset}) cannot be off the end of the sequence \
                 ({raw_seq}). Must be < length of the sequence (0-indexed)."
            )
        }
        let mod_code = ModCodeRepr::parse(raw_mod_code).context(format!(
            "failed to parse modification code {raw_mod_code}"
        ))?;

        let (before, after) = {
            let (b, a) = raw_seq.split_at(offset);
            (b.chars().collect::<Vec<char>>(), a.chars().collect::<Vec<char>>())
        };

        let canonical_base = DnaBase::parse(after[0]).context(format!(
            "failed to parse primary sequence base {}",
            after[0]
        ))?;

        match mod_code_lookup.get(&mod_code) {
            Some(dna_base) if dna_base != &canonical_base => {
                bail!(
                    "association of {mod_code} with primary sequence base \
                     {canonical_base:?} is different from association in \
                     bedMethyl ({dna_base:?}). {raw_seq} {raw_mod_code} \
                     {raw_offset} parsing failed."
                )
            }
            _ => {}
        }

        if before.len() > context_size[0]
            || (after.len() - 1usize) > context_size[1]
        {
            bail!(
                "known motif {raw_seq} is too large for context, [{},{}]",
                context_size[0],
                context_size[1]
            )
        }

        let left_context = before.len() as i8;
        let before = before
            .into_iter()
            .enumerate()
            .map(|(i, b)| {
                let offset = i as i8 - left_context;
                IupacBase::parse_char(b).map(|x| (offset, x))
            })
            .collect::<anyhow::Result<Vec<(i8, IupacBase)>>>()?;
        let after = after
            .into_iter()
            .skip(1)
            .enumerate()
            .map(|(offset, b)| {
                IupacBase::parse_char(b).map(|x| ((offset + 1usize) as i8, x))
            })
            .collect::<anyhow::Result<Vec<(i8, IupacBase)>>>()?;

        let seq = before
            .into_iter()
            .chain(after.into_iter())
            .filter(|(_, base)| !base.is_n())
            .collect::<BTreeMap<i8, IupacBase>>();
        let multi_sequence = MultiSequence::new(mod_code, seq);
        Ok(Self { canonical_base, multi_sequence })
    }

    fn extend_motif(
        mut self,
        kmer_subset: &KmerSubset,
        mod_db: &KmerModificationDb,
        cache: &RwLock<FxHashMap<String, (f32, u64)>>,
        min_sites: u64,
        frac_sites_thresh: f32,
        min_log_odds: f32,
        extend_iters: usize,
    ) -> Self {
        let focus_position = mod_db.focus_position();
        let count_other_bases =
            |counts: &FxHashMap<RawBase, FxHashMap<usize, u32>>,
             base: &RawBase,
             position: usize|
             -> u32 {
                counts
                    .iter()
                    .filter_map(|(b, counts)| {
                        if b != base {
                            Some(*counts.get(&position).unwrap_or(&0u32))
                        } else {
                            None
                        }
                    })
                    .sum::<u32>()
            };

        'extend_loop: for _it in 0..extend_iters {
            let high_matches =
                kmer_subset.get_high_mod_contexts(&self, focus_position);
            let low_matches =
                kmer_subset.get_low_mod_contexts(&self, focus_position);
            let high_count =
                high_matches.iter().map(|(_, c)| *c as u64).sum::<u64>();
            let low_count =
                low_matches.iter().map(|(_, c)| *c as u64).sum::<u64>();

            if high_count < min_sites || low_count < min_sites {
                // debug!("num sites too low {high_count}, {low_count}");
                break 'extend_loop;
            }

            let high_mod_frac =
                mod_db.get_fraction_modified_cached(&self, cache);
            if high_mod_frac > frac_sites_thresh {
                // debug!("fraction of sites satisfied {high_mod_frac}");
                break 'extend_loop;
            } else {
                // debug!("at iteration {it} fraction of sites
                // {high_mod_frac}");
            }

            let high_base_counts =
                util::aggregate_base_counts_on_position(&high_matches);
            let low_base_counts =
                util::aggregate_base_counts_on_position(&low_matches);

            let mut log_odds_per_base_per_pos = (0..=(mod_db.context_bases[0]
                + mod_db.context_bases[1]))
                .into_par_iter()
                .map(|idx| {
                    let motif_relative_pos =
                        (idx as i8) - mod_db.context_bases[0] as i8;
                    (idx, motif_relative_pos)
                })
                .filter(|(_idx, mot_pos)| *mot_pos != 0i8)
                .flat_map(|(idx, motif_position)| {
                    let log_odds_for_bases = BASES
                        .iter()
                        .filter(|base| {
                            !self.contains_base(base, motif_position)
                                && !self.is_fixed_position(motif_position)
                        })
                        .map(|base| {
                            let high_count = high_base_counts
                                .get(base)
                                .and_then(|counts| counts.get(&idx))
                                .unwrap_or(&0u32);
                            (base, high_count)
                        })
                        .map(|(base, high_count)| {
                            let low_count = low_base_counts
                                .get(base)
                                .and_then(|counts| counts.get(&idx))
                                .unwrap_or(&0u32);
                            // .map(|low_count| (base, high_count, low_count))
                            (base, high_count, low_count)
                        })
                        .map(|(base, &high_pos, &low_pos)| {
                            let high_neg =
                                count_other_bases(&high_base_counts, base, idx);
                            let low_neg =
                                count_other_bases(&low_base_counts, base, idx);
                            let log_odds = util::log_odds(
                                low_pos, low_neg, high_pos, high_neg,
                            );
                            (*base, log_odds, motif_position)
                        })
                        .collect::<Vec<(RawBase, f32, i8)>>();
                    log_odds_for_bases
                })
                .collect::<Vec<(RawBase, f32, i8)>>();
            log_odds_per_base_per_pos.par_sort_by(
                |(_raw_base_a, log_odds_a, _motif_position_a),
                 (_raw_base_b, log_odds_b, _motif_position_b)| {
                    let abs_log_odds_b = log_odds_b.abs();
                    let abs_log_odds_a = log_odds_a.abs();
                    abs_log_odds_b
                        .partial_cmp(&abs_log_odds_a)
                        .expect("should_compare")
                },
            );
            // debug!("before update {}", &self);
            assert!(
                !log_odds_per_base_per_pos.is_empty(),
                "log-odds per base should not be empty"
            );
            let (max_base, max_log_odds, max_motif_pos) =
                log_odds_per_base_per_pos[0];
            if max_log_odds.abs() < min_log_odds {
                // debug!("Log odds below threshold: {}", max_log_odds.abs());
                break 'extend_loop;
            }

            let allowed_bases = if max_log_odds > 0f32 {
                vec![max_base]
            } else {
                log_odds_per_base_per_pos
                    // .clone()
                    .iter()
                    .filter_map(|(base, log_odds, position)| {
                        if (*log_odds > 0f32) && *position == max_motif_pos {
                            Some(*base)
                        } else {
                            None
                        }
                    })
                    .unique()
                    .collect::<Vec<RawBase>>()
            };
            if allowed_bases.is_empty() {
                debug!("zero allowed bases, stopping early with {self}");
                break 'extend_loop;
            }
            assert!(
                !allowed_bases.is_empty(),
                "allowed bases cannot be empty, working on {self} \
                 {log_odds_per_base_per_pos:?}"
            );
            // debug!(
            //     "allowed bases {:?}",
            //     allowed_bases.iter().map(|c| *c as
            // char).collect::<Vec<char>>() );
            let allowed_base = allowed_bases
                .into_iter()
                .fold(IupacBase::Hole, |acc, next| acc.add(next));
            // debug!(
            //     "Iupac allowed base {allowed_base:?}, position
            // {max_motif_pos}" );
            self.update_motif_sequence(max_motif_pos, allowed_base, true);
            // debug!("after update {}", &self);
        }

        self
    }

    fn update_motif_sequence(
        &mut self,
        motif_position: i8,
        new_base: IupacBase,
        intersect: bool,
    ) {
        let base = if let Some(current_base) =
            self.multi_sequence.seq.remove(&motif_position)
        {
            if intersect {
                current_base.intersect(new_base)
            } else {
                current_base.union(new_base)
            }
        } else {
            new_base
        };
        assert!(
            self.multi_sequence.seq.insert(motif_position, base).is_none(),
            "should have added entry"
        );
        self.multi_sequence.clean();
    }

    fn contains_base(&self, base: &RawBase, motif_position: i8) -> bool {
        self.multi_sequence
            .seq
            .get(&motif_position)
            .map(|seq_base| seq_base.matches(*base))
            .unwrap_or(false)
    }

    fn is_fixed_position(&self, motif_position: i8) -> bool {
        self.multi_sequence
            .seq
            .get(&motif_position)
            .map(|base| base.is_fixed())
            .unwrap_or(false)
    }

    // panics if pos not in multiseq
    #[inline]
    fn exchange_base(&self, pos: &i8, base: RawBase) -> Self {
        let mut alt = Self {
            canonical_base: self.canonical_base,
            multi_sequence: self.multi_sequence.clone(),
        };
        let new_base = IupacBase::from_base_unchecked(base);
        assert!(
            alt.multi_sequence.seq.insert(*pos, new_base).is_some(),
            "should be replacing value"
        );
        alt
    }

    #[inline]
    fn propose_additional_bases(
        &self,
        mod_db: &KmerModificationDb,
        cache: &RwLock<FxHashMap<String, (f32, u64)>>,
        frac_sites_thresh: f32,
    ) -> FxHashMap<i8, IupacBase> {
        let pos_bases = iproduct!(self.multi_sequence.seq.iter(), BASES)
            .filter_map(|((pos, curr_base), base)| {
                if curr_base.matches(base) {
                    None
                } else {
                    Some((pos, base))
                }
            });

        let candidate_positions = pos_bases
            .par_bridge()
            .filter_map(|(pos, new_base)| {
                let alt = self.exchange_base(pos, new_base);
                let frac_mod = mod_db.get_fraction_modified_cached(&alt, cache);
                if frac_mod > frac_sites_thresh {
                    Some((*pos, new_base))
                } else {
                    None
                }
            })
            .collect::<Vec<(i8, u8)>>();

        candidate_positions.into_iter().fold(
            FxHashMap::default(),
            |mut agg, (pos, next_base)| {
                let base = agg.entry(pos).or_insert(IupacBase::Hole);
                *base = base.add(next_base);
                agg
            },
        )
    }

    fn add_bases_to_motif(
        mut self,
        mod_db: &KmerModificationDb,
        cache: &RwLock<FxHashMap<String, (f32, u64)>>,
        frac_sites_thresh: f32,
    ) -> Self {
        loop {
            let proposed_additions =
                self.propose_additional_bases(mod_db, cache, frac_sites_thresh);
            if proposed_additions.is_empty() {
                break;
            } else {
                for (pos, base) in proposed_additions {
                    self.update_motif_sequence(pos, base, false);
                }
            }
        }

        self
    }

    fn check_remove_base(
        &self,
        mod_db: &KmerModificationDb,
        cache: &RwLock<FxHashMap<String, (f32, u64)>>,
        pos_to_remove: i8,
        base_to_remove: RawBase,
    ) -> (f32, u64) {
        let new_seq = self
            .multi_sequence
            .seq
            .iter()
            .map(|(pos, base)| {
                if *pos == pos_to_remove {
                    (*pos, (*base).remove_to_n(base_to_remove))
                } else {
                    (*pos, *base)
                }
            })
            .filter(|(_, base)| !base.is_n())
            .collect::<BTreeMap<i8, IupacBase>>();
        let ms = MultiSequence::new(self.multi_sequence.mod_code, new_seq);

        let alt =
            Self { canonical_base: self.canonical_base, multi_sequence: ms };
        mod_db.get_counts_and_frac_cached(&alt, cache)
    }

    fn contract_motif(
        mut self,
        mod_db: &KmerModificationDb,
        cache: &RwLock<FxHashMap<String, (f32, u64)>>,
        frac_sites_thresh: f32,
        debug: bool,
    ) -> Self {
        loop {
            let candidate_positions =
                iproduct!(self.multi_sequence.seq.iter(), BASES.iter())
                    .filter(|((_position, motif_base), base)| {
                        motif_base.matches(**base)
                    })
                    .filter_map(|((pos_to_remove, _), base_to_remove)| {
                        let (frac_sites_high, alt_high_count) = self
                            .check_remove_base(
                                &mod_db,
                                cache,
                                *pos_to_remove,
                                *base_to_remove,
                            );
                        if frac_sites_high > frac_sites_thresh {
                            Some((
                                frac_sites_high,
                                alt_high_count,
                                *pos_to_remove,
                                *base_to_remove,
                            ))
                        } else {
                            None
                        }
                    })
                    .collect::<Vec<(f32, u64, i8, u8)>>();
            let contract_position = candidate_positions.iter().max_by(
                |(f, num_high_sites_a, _, a), (g, num_high_sites_b, _, b)| {
                    match f.partial_cmp(g).expect("should compare") {
                        Ordering::Equal => {
                            match num_high_sites_a.cmp(num_high_sites_b) {
                                Ordering::Equal => {
                                    let mut message = format!(
                                        "same number of high-mod sites \
                                         ({num_high_sites_a}) and frac ({f}) \
                                         mod for motif {self}\n \
                                         {candidate_positions:?}, picking \
                                         arbitrary position?"
                                    );
                                    if a == b {
                                        message.push_str(
                                            " Bases are also the same.",
                                        );
                                    }
                                    debug!("{message}");
                                    (*a as char).cmp(&(*b as char))
                                }
                                o @ _ => o,
                            }
                        }
                        o @ _ => o,
                    }
                },
            );

            if let Some((
                _max_frac_mod,
                _max_sites,
                max_position,
                max_base_to_remove,
            )) = contract_position
            {
                let mut base =
                    self.multi_sequence.seq.remove(max_position).unwrap();
                base = base.remove_to_n(*max_base_to_remove);
                if !base.is_n() {
                    assert!(self
                        .multi_sequence
                        .seq
                        .insert(*max_position, base)
                        .is_none());
                }
            } else {
                if debug {
                    info!("||> finished contracting motif: {self}")
                }
                break;
            }
        }
        self
    }

    fn refine(
        self,
        mod_db: &KmerModificationDb,
        cache: &RwLock<FxHashMap<String, (f32, u64)>>,
        kmer_subset: &KmerSubset,
        min_sites: u64,
        frac_sites_thresh: f32,
        min_log_odds: f32,
        stage: Stage,
    ) -> Self {
        let starting_pattern = self.to_string();
        let debug = false;
        let motif = self.extend_motif(
            &kmer_subset,
            mod_db,
            cache,
            min_sites,
            frac_sites_thresh,
            min_log_odds,
            24, // todo
        );
        let mut motif =
            motif.add_bases_to_motif(&mod_db, cache, frac_sites_thresh);
        let mut last_motif = motif.clone();
        let mut n_iters = 0usize;
        loop {
            motif =
                motif.contract_motif(&mod_db, cache, frac_sites_thresh, debug);
            motif = motif.add_bases_to_motif(&mod_db, cache, frac_sites_thresh);
            if motif == last_motif {
                debug!(
                    from_motif = starting_pattern,
                    motif = motif.to_string(),
                    mod_code = motif.mod_code().to_string(),
                    action = Action::Refined.to_string(),
                    stage = stage.to_string(),
                    "refined {starting_pattern} into {motif}, took {n_iters} \
                     iterations",
                );
                break;
            } else {
                last_motif = motif.clone();
                n_iters += 1;
            }
        }
        motif
    }

    fn is_superset(&self, other: &Self) -> bool {
        self.multi_sequence.is_superset(&other.multi_sequence)
    }

    fn is_subset(&self, other: &Self) -> bool {
        if self.multi_sequence.mod_code != other.multi_sequence.mod_code {
            return false;
        }
        if self.multi_sequence == other.multi_sequence {
            return true;
        }

        for (pos, base) in &self.multi_sequence.seq {
            match other.multi_sequence.seq.get(&pos) {
                // other has this position
                Some(other_base) => {
                    // if other_base is not a superset, this
                    // base cannot be a subset of it (kinda backwards)
                    if !other_base.is_superset(base) {
                        return false;
                    }
                }
                // if other does not have this position, other
                // self cannot be a subset of it
                None => {
                    return false;
                }
            }
        }

        true
    }

    pub(super) fn compare(
        &self,
        other: &Self,
        context_size: [usize; 2],
    ) -> MotifRelationship {
        if self == other {
            MotifRelationship::Equal
        } else if self.is_subset(other) {
            MotifRelationship::Subset
        } else if self.is_superset(other) {
            MotifRelationship::Superset
        } else {
            let edit_distance = self
                .multi_sequence
                .edit_distance(&other.multi_sequence, context_size);
            MotifRelationship::Disjoint { edit_distance }
        }
    }

    pub(super) fn mod_code(&self) -> ModCodeRepr {
        self.multi_sequence.mod_code
    }
}

impl Display for EnrichedMotif {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.multi_sequence)
    }
}

fn merge_motif(
    matches: FxHashMap<usize, FxHashSet<usize>>,
    enriched_motifs: Vec<EnrichedMotif>,
) -> Vec<EnrichedMotif> {
    let motif_idxs_to_discard = matches
        .iter()
        .max_by(|(_, a), (_, b)| a.len().cmp(&b.len()))
        .map(|(_, xs)| xs)
        .unwrap();

    enriched_motifs
        .into_iter()
        .enumerate()
        .filter(|(i, _)| !motif_idxs_to_discard.contains(i))
        .map(|(_, m)| m)
        .collect()
}

fn merge_motifs(mut enriched_motifs: Vec<EnrichedMotif>) -> Vec<EnrichedMotif> {
    // might be the case that you only need to do this loop once..
    loop {
        let matches = iproduct!(
            enriched_motifs.iter().enumerate(),
            enriched_motifs.iter().enumerate()
        )
        .filter_map(|((i, a), (j, b))| {
            if i == j {
                None
            } else {
                let is_superset = a.is_superset(b);
                if is_superset {
                    Some((i, j))
                } else {
                    None
                }
            }
        })
        .fold(FxHashMap::default(), |mut agg, (i, j)| {
            agg.entry(i).or_insert(FxHashSet::default()).insert(j);
            agg
        });

        if matches.iter().all(|(_, subsets)| subsets.is_empty()) {
            // debug!("no more merging to do");
            break;
        }

        enriched_motifs = merge_motif(matches, enriched_motifs);
    }

    enriched_motifs
}

fn get_fixed_length_motifs(
    canonical_base: DnaBase,
    mod_code: ModCodeRepr,
    initial_context: [usize; 2],
    min_log_odds_thresh: f32,
    kmer_mod_db: &KmerModificationDb,
    kmer_mask: &KmerMask,
) -> Vec<EnrichedMotif> {
    let sliced_low = kmer_mod_db.low_mod_table.sliced_kmers(
        canonical_base,
        &mod_code,
        &kmer_mask.low_mod_mask,
        &initial_context,
    );
    let sliced_high = kmer_mod_db.high_mod_table.sliced_kmers(
        canonical_base,
        &mod_code,
        &kmer_mask.high_mod_mask,
        &initial_context,
    );
    let low_total = sliced_low.values().sum::<u32>() as f32;
    let high_total = sliced_high.values().sum::<u32>() as f32;

    // debug!("high total {high_total}, low total {low_total}");
    // debug!("sliced_high has {} kmers", sliced_high.len());
    // debug!("sliced_low has {} kmers", sliced_low.len());

    let enriched_kmers = sliced_low
        .iter()
        .filter_map(|(kmer, low_count)| {
            // requires that we do the intersection of the two
            sliced_high
                .get(kmer)
                .map(|high_count| (kmer, low_count, high_count))
        })
        .filter_map(|(kmer, lo_count, hi_count)| {
            let numer = *hi_count as f32 * low_total;
            let denom = *lo_count as f32 * high_total;
            // debug!("kmer {}, {numer} / {denom}", string_kmer(kmer));
            let log_odds = (numer / denom).log2();
            if log_odds >= min_log_odds_thresh {
                Some(*kmer)
            } else {
                None
            }
        })
        .sorted()
        .collect::<Vec<_>>();

    // stopping condition
    if enriched_kmers.is_empty() {
        return Vec::new();
    } // otherwise join them up

    debug!(
        mod_code = mod_code.to_string(),
        stage = Stage::Seeded.to_string(),
        "there are {} enriched kmers for mod {mod_code}",
        enriched_kmers.len()
    );

    let joined_kmers = {
        let mut joined_kmers =
            FxHashMap::<KmerRef, FxHashSet<KmerRef>>::default();
        for (idx, kmer1) in
            enriched_kmers.iter().take(enriched_kmers.len() - 1).enumerate()
        {
            for kmer2 in enriched_kmers.iter().skip(idx + 1) {
                let ed = bio::alignment::distance::simd::hamming(kmer1, kmer2);
                assert!(ed >= 1, "should not have the same kmers!?");
                if ed == 1 {
                    joined_kmers
                        .entry(kmer1)
                        .or_insert(FxHashSet::default())
                        .insert(kmer2);
                    joined_kmers
                        .entry(kmer2)
                        .or_insert(FxHashSet::default())
                        .insert(kmer1);
                }
            }
        }
        joined_kmers
    };

    let core_kmers = {
        let mut core_kmers = Vec::<FxHashSet<KmerRef>>::new();
        let mut kmer_to_core = FxHashMap::<KmerRef, usize>::default();
        for kmer in enriched_kmers.iter() {
            let (motif_kmers, idx) = if let Some(idx) = kmer_to_core.get(kmer) {
                (core_kmers.get_mut(*idx).unwrap(), *idx)
            } else {
                let tmp_motif: FxHashSet<KmerRef> =
                    [*kmer].into_iter().collect();
                core_kmers.push(tmp_motif);
                let idx = core_kmers.len() - 1;
                kmer_to_core.insert(kmer, idx);
                (core_kmers.last_mut().unwrap(), idx)
            };
            if let Some(kmer2s) = joined_kmers.get(kmer) {
                for kmer2 in kmer2s {
                    if !motif_kmers.contains(kmer2) {
                        motif_kmers.insert(*kmer2);
                        kmer_to_core.insert(kmer2, idx);
                    }
                }
            }
        }
        core_kmers
    };

    // info!("core kmer sets");
    // for core_set in core_kmers.iter().sorted_by(|a, b| a.len().cmp(&b.len()))
    // {     let stringed = core_set.iter().map(|k|
    // string_kmer(k)).sorted().join(", ");     info!("{stringed} {}",
    // core_set.len()); }

    core_kmers
        .into_iter()
        .sorted_by(|a, b| a.len().cmp(&b.len()))
        .map(|kmers| {
            EnrichedMotif::new(
                canonical_base,
                MultiSequence::from_iter(
                    kmers.into_iter(),
                    initial_context,
                    mod_code,
                ),
            )
        })
        .collect::<Vec<EnrichedMotif>>()
}

fn get_seeded_motifs<'a>(
    canonical_base: DnaBase,
    mod_code: ModCodeRepr,
    mod_db: &'a KmerModificationDb,
    cache: &RwLock<FxHashMap<String, (f32, u64)>>,
    initial_context_size: &[usize],
    min_log_odds: f32,
    min_sites: u64,
    frac_sites_thresh: f32,
    multi_progress: &MultiProgress,
) -> (Vec<EnrichedMotif>, KmerMask<'a>) {
    let mut mod_motifs = Vec::<EnrichedMotif>::new();
    let mut kmer_mask = KmerMask::default();
    let mut kmer_subset =
        mod_db.get_kmer_subset(canonical_base, &kmer_mask, mod_code);
    let seed_pb = multi_progress.add(get_ticker());
    seed_pb
        .set_message(format!("iterations searching with seeds ({mod_code})"));
    loop {
        seed_pb.inc(1u64);
        let fixed_length_motifs = get_fixed_length_motifs(
            canonical_base,
            mod_code,
            [initial_context_size[0], initial_context_size[1]],
            min_log_odds,
            &mod_db,
            &kmer_mask,
        )
        .into_iter()
        .collect::<Vec<EnrichedMotif>>();
        if fixed_length_motifs.is_empty() {
            debug!(
                mod_code = mod_code.to_string(),
                stage = Stage::Seeded.to_string(),
                "zero fixed length motifs, finished in {} iterations",
                seed_pb.position()
            );
            break;
        }

        let refined_motifs = fixed_length_motifs
            .into_par_iter()
            .map(|motif| {
                motif.refine(
                    &mod_db,
                    cache,
                    &kmer_subset,
                    min_sites,
                    frac_sites_thresh,
                    min_log_odds,
                    Stage::Seeded,
                )
            })
            .collect::<FxHashSet<EnrichedMotif>>();

        let mut refined_motifs = refined_motifs
            .into_par_iter()
            .filter(|refined_motif| {
                // use the "total" method here.
                let (high_count, low_count) =
                    mod_db.get_mod_counts(refined_motif, &KmerMask::default());
                let frac_high =
                    high_count as f32 / (high_count + low_count) as f32;
                if high_count < min_sites {
                    debug!(
                        action = Action::Discard.to_string(),
                        motif = refined_motif.to_string(),
                        mod_code = refined_motif.mod_code().to_string(),
                        stage = Stage::Seeded.to_string(),
                        require = min_sites,
                        value = high_count,
                        "discarding {refined_motif}, high-modified sites too \
                         low, {high_count} requires {min_sites}"
                    );
                    false
                } else if frac_high <= frac_sites_thresh {
                    debug!(
                        action = Action::Discard.to_string(),
                        motif = refined_motif.to_string(),
                        mod_code = refined_motif.mod_code().to_string(),
                        stage = Stage::Seeded.to_string(),
                        require = frac_sites_thresh,
                        value = frac_high,
                        "discarding {refined_motif}, fraction modified too \
                         low, {frac_high}, requires {frac_sites_thresh}"
                    );
                    false
                } else {
                    debug!(
                        action = Action::Found.to_string(),
                        motif = refined_motif.to_string(),
                        mod_code = refined_motif.mod_code().to_string(),
                        stage = Stage::Seeded.to_string(),
                        "keeping refined motif {refined_motif}, it has \
                         {high_count} sites and {frac_high}% high"
                    );
                    true
                }
            })
            .filter(|refined_motif| {
                match mod_motifs
                    .iter()
                    .find(|mot| refined_motif.is_superset(mot))
                {
                    Some(mot) => {
                        debug!(
                            action = Action::Discard.to_string(),
                            motif = refined_motif.to_string(),
                            mod_code = refined_motif.mod_code().to_string(),
                            stage = Stage::Seeded.to_string(),
                            "discarding {refined_motif} as a superset of or \
                             equivalent to a previously found motif, {mot}"
                        );
                        false
                    }
                    None => true, // keep it
                }
            })
            .collect::<Vec<EnrichedMotif>>();

        if refined_motifs.is_empty() {
            debug!(
                "zero refined motifs, finished in {} iterations",
                seed_pb.position()
            );
            break;
        }

        // remove the contexts from these motifs
        let new_mask = match kmer_mask.update_with_check(
            &refined_motifs,
            mod_db,
            mod_code,
            Stage::Seeded,
        ) {
            Either::Left(new_mask) => {
                mod_motifs.append(&mut refined_motifs);
                mod_motifs = merge_motifs(mod_motifs);
                debug!(
                    stage = Stage::Seeded.to_string(),
                    mod_code = mod_code.to_string(),
                    "should always remove contexts, finishing seeded search \
                     after {} iterations",
                    seed_pb.position()
                );
                kmer_mask = new_mask;
                break;
            }
            Either::Right(new_mask) => new_mask,
        };

        mod_motifs.append(&mut refined_motifs);
        mod_motifs = merge_motifs(mod_motifs);

        let n_high_contexts =
            mod_db.count_high_mod_contexts(&mod_code, &new_mask);
        let n_low_contexts =
            mod_db.count_low_mod_contexts(&mod_code, &new_mask);

        debug!(
            stage = Stage::Seeded.to_string(),
            mod_code = mod_code.to_string(),
            "at iter {} there are {n_high_contexts} remaining high contexts \
             and {n_low_contexts} remaining low contexts",
            seed_pb.position()
        );

        kmer_subset =
            mod_db.get_kmer_subset(canonical_base, &new_mask, mod_code);
        kmer_mask = new_mask;
    }
    seed_pb.finish_and_clear();

    (mod_motifs, kmer_mask)
}

fn find_motifs_for_mod(
    search_config: SearchConfig,
    canonical_base: DnaBase,
    mod_code: ModCodeRepr,
    mod_db: &KmerModificationDb,
    initial_context_size: &[usize],
    min_log_odds: f32,
    min_sites: u64,
    frac_sites_thresh: f32,
    skip_search: bool,
    exhaustive_search_kmer_size: usize,
    exhaustive_search_min_log_odds: f32,
    multi_progress: &MultiProgress,
) -> Vec<EnrichedMotifData> {
    let cache = RwLock::new(FxHashMap::<String, (f32, u64)>::default());
    let n_high_contexts =
        mod_db.count_high_mod_contexts(&mod_code, &KmerMask::default());
    let n_low_contexts =
        mod_db.count_low_mod_contexts(&mod_code, &KmerMask::default());

    debug!(
        mod_code = mod_code.to_string(),
        "At start, there are {} high contexts and {} low contexts",
        n_high_contexts,
        n_low_contexts
    );

    let start = std::time::Instant::now();
    let (mut seeded_motifs, mut kmer_mask) = get_seeded_motifs(
        canonical_base,
        mod_code,
        mod_db,
        &cache,
        initial_context_size,
        min_log_odds,
        min_sites,
        frac_sites_thresh,
        multi_progress,
    );
    let took = start.elapsed();
    debug!(
        mod_code = mod_code.to_string(),
        stage = Stage::Seeded.to_string(),
        "{}",
        util::log_motifs(mod_code, &seeded_motifs, "seeded motifs", took)
    );

    let mut kmer_subset =
        mod_db.get_kmer_subset(canonical_base, &kmer_mask, mod_code);
    let seedless_pb = multi_progress.add(get_ticker());
    seedless_pb
        .set_message(format!("({mod_code}) finding motifs without seed"));
    'seedless_loop: loop {
        seedless_pb.inc(1);
        let motif = EnrichedMotif::new_empty(mod_code, canonical_base).refine(
            &mod_db,
            &cache,
            &kmer_subset,
            min_sites,
            frac_sites_thresh,
            min_log_odds,
            Stage::Seedless,
        );
        let (high_count, low_count) =
            mod_db.get_mod_counts(&motif, &KmerMask::default());
        let frac_high = high_count as f32 / (high_count + low_count) as f32;
        if high_count >= min_sites && frac_high > frac_sites_thresh {
            debug!(
                mod_code = mod_code.to_string(),
                stage = Stage::Seedless.to_string(),
                motif = motif.to_string(),
                action = Action::Found.to_string(),
                "({mod_code}) found motif {motif} in seedless stage"
            );
            seeded_motifs.push(motif);
            seeded_motifs = merge_motifs(seeded_motifs);
            debug!(
                mod_code = mod_code.to_string(),
                stage = Stage::Seedless.to_string(),
                "({mod_code}) seeded motifs is now {}",
                seeded_motifs.iter().map(|m| m.to_string()).join(",")
            );
            match kmer_mask.update_with_check(
                &seeded_motifs,
                &mod_db,
                mod_code,
                Stage::Seedless,
            ) {
                Either::Left(new_mask) => {
                    debug!(
                        "({mod_code}) done searching for seedless motifs, \
                         didn't remove any additional contexts , found {}",
                        seedless_pb.position()
                    );
                    kmer_mask = new_mask;
                    break 'seedless_loop;
                }
                Either::Right(new_mask) => {
                    kmer_mask = new_mask;
                    kmer_subset = mod_db.get_kmer_subset(
                        canonical_base,
                        &kmer_mask,
                        mod_code,
                    );
                    seedless_pb.inc(1);
                    continue 'seedless_loop;
                }
            }
        } else {
            if high_count < min_sites {
                debug!(
                    action = Action::Discard.to_string(),
                    motif = motif.to_string(),
                    mod_code = motif.mod_code().to_string(),
                    stage = Stage::Seedless.to_string(),
                    require = min_sites,
                    value = high_count,
                    "discarding {motif}, high-modified sites too low, \
                     {high_count} requires {min_sites}"
                );
            }
            if frac_high <= frac_sites_thresh {
                debug!(
                    action = Action::Discard.to_string(),
                    motif = motif.to_string(),
                    mod_code = motif.mod_code().to_string(),
                    stage = Stage::Seedless.to_string(),
                    require = frac_sites_thresh,
                    value = frac_high,
                    "discarding {motif}, fraction modified too low, \
                     {frac_high}, requires {frac_sites_thresh}"
                );
            }
            debug!(
                mod_code = mod_code.to_string(),
                stage = Stage::Seedless.to_string(),
                "({mod_code}) done searching for seedless motifs, found {}",
                seedless_pb.position()
            );
            break 'seedless_loop;
        }
    }
    seedless_pb.finish_and_clear();

    let final_motifs = if !skip_search {
        debug!(
            mod_code = mod_code.to_string(),
            stage = Stage::Search.to_string(),
            "performing search"
        );
        let start = std::time::Instant::now();
        let exhaustive_seed_motifs = find_exhaustive_seed_motifs(
            search_config,
            canonical_base,
            mod_code,
            exhaustive_search_kmer_size,
            exhaustive_search_min_log_odds,
            min_log_odds,
            min_sites,
            frac_sites_thresh,
            mod_db,
            &cache,
            kmer_mask,
            multi_progress,
        );
        let end = start.elapsed();

        let (exhaustive_seed_motifs, stopped_early) =
            match exhaustive_seed_motifs {
                Ok(x) => (x, false),
                Err(x) => (x, true),
            };

        debug!(
            mod_code = mod_code.to_string(),
            stage = Stage::Search.to_string(),
            "{}",
            util::log_motifs(
                mod_code,
                &exhaustive_seed_motifs,
                "exhaustive search",
                end,
            )
        );

        if stopped_early {
            multi_progress.suspend(|| {
                error!("stopped early without finishing search for {mod_code}")
            })
        };

        let exhaustive_seed_filt = exhaustive_seed_motifs
            .into_iter()
            .filter(|motif| {
                let subset_of_motifs = seeded_motifs
                    .iter()
                    .find(|seeded_mot| motif.is_subset(seeded_mot));
                if let Some(mot) = subset_of_motifs {
                    debug!(
                        action = Action::Discard.to_string(),
                        motif = motif.to_string(),
                        mod_code = motif.mod_code().to_string(),
                        stage = Stage::Search.to_string(),
                        "discarding {motif} as subset of previously found \
                         motif {mot}"
                    );
                    false
                } else {
                    debug!(
                        action = Action::Found.to_string(),
                        motif = motif.to_string(),
                        mod_code = motif.mod_code().to_string(),
                        stage = Stage::Search.to_string(),
                        "non-redundant motif from search {motif}"
                    );
                    true
                }
            })
            .collect::<Vec<EnrichedMotif>>();
        merge_motifs(
            exhaustive_seed_filt
                .into_iter()
                .chain(seeded_motifs.into_iter())
                .collect(),
        )
        .into_par_iter()
        .map(|motif| mod_db.get_enriched_motif_data(&motif))
        .collect::<Vec<EnrichedMotifData>>()
    } else {
        debug!(mod_code = mod_code.to_string(), "skipping search");
        seeded_motifs
            .into_par_iter()
            .map(|motif| mod_db.get_enriched_motif_data(&motif))
            .collect::<Vec<EnrichedMotifData>>()
    };

    final_motifs
}

fn find_exhaustive_seed_motifs<'a>(
    optim_config: SearchConfig,
    canonical_base: DnaBase,
    mod_code: ModCodeRepr,
    kmer_length: usize,
    search_min_log_odds: f32,
    refine_log_odds: f32,
    refine_min_sites: u64,
    refine_sites_thresh: f32,
    mod_db: &'a KmerModificationDb,
    cache: &RwLock<FxHashMap<String, (f32, u64)>>,
    mut kmer_mask: KmerMask<'a>,
    multi_progress: &MultiProgress,
) -> Result<Vec<EnrichedMotif>, Vec<EnrichedMotif>> {
    let start_time = std::time::Instant::now();
    let kmers =
        itertools::repeat_n(BASES, kmer_length).multi_cartesian_product();
    let positions = (0usize
        ..=(mod_db.context_bases[0] + mod_db.context_bases[1]))
        .filter(|idx| *idx != mod_db.focus_position());
    let position_mers = positions.combinations(kmer_length);
    let combs = iproduct!(kmers, position_mers).collect::<Vec<_>>();

    let mut kmer_subset =
        mod_db.get_kmer_subset(canonical_base, &kmer_mask, mod_code);
    let mut position_bools =
        kmer_subset.get_position_bools(mod_db.focus_position());

    let n_combs = combs.len();
    let get_scoring_pb = || {
        let pb = multi_progress.add(get_subroutine_progress_bar(n_combs));
        pb.set_message(format!("({mod_code}) scoring seeds"));
        pb
    };

    let mut kmer_positions = combs
        .par_iter()
        .progress_with(get_scoring_pb())
        .filter_map(|(kmer, positions)| {
            let lo = position_bools
                .get_counts(kmer.as_slice(), positions.as_slice())
                .log_odds();
            if lo >= search_min_log_odds {
                Some((kmer, positions, lo))
            } else {
                None
            }
        })
        .collect::<Vec<_>>()
        .into_iter()
        .sorted_by(|(_, _, a), (_, _, b)| {
            a.partial_cmp(b).unwrap_or(Ordering::Equal)
        })
        .collect::<Vec<_>>();

    // if !kmer_positions.is_empty() {
    //     let log_oddss =
    //         kmer_positions.iter().map(|(_, _, lo)|
    // *lo).collect::<Vec<f32>>();     multi_progress.suspend(|| {
    //         debug!(
    //         mod_code = mod_code.to_string(),
    //         stage = Stage::Search.to_string(),
    //         log_odds = ?log_oddss,
    //         );
    //     });
    // };

    let mut n_iter = 1usize;
    let mut results = Vec::new();
    'search_loop: loop {
        if kmer_positions.is_empty() {
            debug!(
                mod_code = mod_code.to_string(),
                stage = Stage::Search.to_string(),
                "zero seeds to search"
            );
            break;
        }
        multi_progress.suspend(|| {
            debug!(
                mod_code = mod_code.to_string(),
                stage = Stage::Search.to_string(),
                "there are {} total seeds to search at {n_iter}",
                kmer_positions.len()
            );
            info!(
                "there are {} seeds to search at {n_iter} for {mod_code}",
                kmer_positions.len()
            );
        });

        let search_batch = match optim_config {
            SearchConfig::FullSearch => kmer_positions.split_off(0),
            SearchConfig::TopFrac { frac, min_seeds, max_seeds }
            | SearchConfig::BatchAndNarrow {
                frac, min_seeds, max_seeds, ..
            }
            | SearchConfig::TimeoutAndNarrow {
                frac,
                min_seeds,
                max_seeds,
                ..
            } => {
                let head_n =
                    (kmer_positions.len() as f32 * frac).ceil() as usize;
                let head_n = std::cmp::min(max_seeds, head_n);
                let head_n = std::cmp::max(head_n, min_seeds);
                if head_n >= kmer_positions.len() {
                    kmer_positions.split_off(0)
                } else {
                    kmer_positions.split_off(kmer_positions.len() - head_n)
                }
            }
            SearchConfig::Timeout { batch_size, .. } => {
                if batch_size >= kmer_positions.len() {
                    kmer_positions.split_off(0)
                } else {
                    kmer_positions.split_off(kmer_positions.len() - batch_size)
                }
            }
        };

        let pb =
            multi_progress.add(get_subroutine_progress_bar(search_batch.len()));
        pb.set_message(format!("({mod_code}, {n_iter}) seeds searched"));
        debug!(
            mod_code = mod_code.to_string(),
            stage = Stage::Search.to_string(),
            "there are {} seeds in batch to search at {n_iter}",
            search_batch.len()
        );

        let enriched_motifs = search_batch
            .into_par_iter()
            .progress_with(pb)
            .filter(|(kmer, positions, _)| {
                position_bools
                    .get_counts(kmer.as_slice(), positions.as_slice())
                    .log_odds()
                    >= search_min_log_odds
            })
            .map(|(kmer, positions, _)| {
                assert_eq!(kmer.len(), positions.len());
                let seq = positions
                    .into_iter()
                    .zip(kmer.into_iter())
                    .map(|(pos, raw_base)| {
                        let offset = *pos as i8 - mod_db.focus_position() as i8;
                        (offset, IupacBase::from_base_unchecked(*raw_base))
                    })
                    .collect::<BTreeMap<i8, IupacBase>>();
                let multi_sequence = MultiSequence::new(mod_code, seq);
                EnrichedMotif::new(canonical_base, multi_sequence)
            })
            .map(|motif| {
                motif.refine(
                    mod_db,
                    cache,
                    &kmer_subset,
                    refine_min_sites,
                    refine_sites_thresh,
                    refine_log_odds,
                    Stage::Search,
                )
            })
            .filter(|refined_motif| {
                let (high_count, low_count) =
                    mod_db.get_mod_counts(refined_motif, &KmerMask::default());
                let frac_high =
                    high_count as f32 / (high_count + low_count) as f32;
                if high_count < refine_min_sites {
                    debug!(
                        action = Action::Discard.to_string(),
                        motif = refined_motif.to_string(),
                        mod_code = refined_motif.mod_code().to_string(),
                        stage = Stage::Search.to_string(),
                        require = refine_min_sites,
                        value = high_count,
                        "discarding {refined_motif}, high-modified sites too \
                         low, {high_count} requires {refine_min_sites}"
                    );
                    false
                } else if frac_high <= refine_sites_thresh {
                    debug!(
                        action = Action::Discard.to_string(),
                        motif = refined_motif.to_string(),
                        mod_code = refined_motif.mod_code().to_string(),
                        stage = Stage::Search.to_string(),
                        require = refine_sites_thresh,
                        value = frac_high,
                        "discarding {refined_motif}, fraction modified too \
                         low, {frac_high}, requires {refine_sites_thresh}"
                    );
                    false
                } else {
                    debug!(
                        action = Action::Found.to_string(),
                        motif = refined_motif.to_string(),
                        mod_code = refined_motif.mod_code().to_string(),
                        stage = Stage::Search.to_string(),
                        "found {refined_motif} during search"
                    );
                    true
                }
            })
            .collect::<FxHashSet<_>>()
            .into_iter()
            .collect::<Vec<EnrichedMotif>>();
        debug!(
            mod_code = mod_code.to_string(),
            stage = Stage::Search.to_string(),
            "found {} enriched motifs at iteration {}",
            enriched_motifs.len(),
            n_iter
        );
        match optim_config {
            // first stopping condition, we're only doing one loop
            SearchConfig::FullSearch | SearchConfig::TopFrac { .. } => {
                assert!(results.is_empty());
                return Ok(enriched_motifs);
            }
            SearchConfig::Timeout { total_time, .. } => {
                let so_far =
                    std::time::Instant::now().duration_since(start_time);
                if so_far >= total_time {
                    multi_progress.suspend(|| {
                        warn!(
                            "stopping search after {}",
                            humantime::format_duration(so_far)
                        )
                    });
                    results.extend(enriched_motifs.into_iter());
                    return if kmer_positions.is_empty() {
                        Ok(results)
                    } else {
                        Err(results)
                    };
                } else {
                    let time_left = total_time - so_far;
                    debug!(
                        mod_code = mod_code.to_string(),
                        stage = Stage::Search.to_string(),
                        "exhaustive search has found {} motifs, have {} \
                         remaining time",
                        results.len(),
                        humantime::format_duration(time_left)
                    );
                    n_iter += 1;
                    continue 'search_loop;
                }
            }
            _ => {}
        }
        // update
        kmer_mask = match kmer_mask.update_with_check(
            &enriched_motifs,
            mod_db,
            mod_code,
            Stage::Search,
        ) {
            // second stopping condition, we haven't updated the mask
            Either::Left(_new_mask) => {
                results.extend(enriched_motifs.into_iter());
                let so_far =
                    std::time::Instant::now().duration_since(start_time);
                multi_progress.suspend(|| {
                    info!(
                        "({mod_code}) didn't remove any contexts, after \
                         {n_iter} iteration(s) stopping, found {} candidate \
                         motifs, took {}",
                        results.len(),
                        humantime::format_duration(so_far)
                    );
                });
                debug!(
                    mod_code = mod_code.to_string(),
                    stage = Stage::Search.to_string(),
                    "didn't remove any contexts, stopping, found {} motifs, \
                     took {}",
                    results.len(),
                    humantime::format_duration(so_far)
                );
                return Ok(results);
            }
            Either::Right(new_mask) => new_mask,
        };
        kmer_subset =
            mod_db.get_kmer_subset(canonical_base, &kmer_mask, mod_code);
        position_bools =
            kmer_subset.get_position_bools(mod_db.focus_position());
        let pb = get_scoring_pb();
        kmer_positions = combs
            .par_iter()
            .progress_with(pb)
            .filter_map(|(kmer, positions)| {
                let lo = position_bools
                    .get_counts(kmer.as_slice(), positions.as_slice())
                    .log_odds();
                if lo >= search_min_log_odds {
                    Some((kmer, positions, lo))
                } else {
                    None
                }
            })
            .collect::<Vec<_>>()
            .into_iter()
            .sorted_by(|(_, _, a), (_, _, b)| {
                a.partial_cmp(b).unwrap_or(Ordering::Equal)
            })
            .collect::<Vec<_>>();

        results.extend(enriched_motifs.into_iter());
        match optim_config {
            SearchConfig::TimeoutAndNarrow {
                max_iters, total_time, ..
            } => {
                let so_far =
                    std::time::Instant::now().duration_since(start_time);
                if so_far >= total_time {
                    multi_progress.suspend(|| {
                        warn!(
                            "({mod_code}) stopping search after {}",
                            humantime::format_duration(so_far)
                        )
                    });
                    return if kmer_positions.is_empty() {
                        debug!(
                            mod_code = mod_code.to_string(),
                            stage = Stage::Search.to_string(),
                            "zero seeds to search at timeout"
                        );
                        Ok(results)
                    } else {
                        debug!(
                            mod_code = mod_code.to_string(),
                            stage = Stage::Search.to_string(),
                            "timed out when {} seeds to search on the next \
                             round",
                            kmer_positions.len()
                        );
                        Err(results)
                    };
                } else if max_iters.map(|i| n_iter >= i).unwrap_or(false) {
                    return if kmer_positions.is_empty() {
                        debug!(
                            mod_code = mod_code.to_string(),
                            stage = Stage::Search.to_string(),
                            "max iterations reached, zero seeds remaining",
                        );
                        Ok(results)
                    } else {
                        debug!(
                            mod_code = mod_code.to_string(),
                            stage = Stage::Search.to_string(),
                            "max iterations reached, {} seeds remaining",
                            kmer_positions.len()
                        );
                        Err(results)
                    };
                } else {
                    let time_left = total_time - so_far;
                    debug!(
                        mod_code = mod_code.to_string(),
                        stage = Stage::Search.to_string(),
                        "Batch and Narrow has found {} motifs, have {} \
                         remaining time",
                        results.len(),
                        humantime::format_duration(time_left)
                    );
                }
            }
            SearchConfig::BatchAndNarrow { max_iters, .. } => {
                if max_iters.map(|i| n_iter >= i).unwrap_or(false) {
                    return if kmer_positions.is_empty() {
                        debug!(
                            mod_code = mod_code.to_string(),
                            stage = Stage::Search.to_string(),
                            "max iterations reached, zero seeds remaining"
                        );
                        Ok(results)
                    } else {
                        debug!(
                            mod_code = mod_code.to_string(),
                            stage = Stage::Search.to_string(),
                            "max iterations reached, {} seeds remaining",
                            kmer_positions.len()
                        );
                        Err(results)
                    };
                }
            }
            SearchConfig::Timeout { .. }
            | SearchConfig::TopFrac { .. }
            | SearchConfig::FullSearch => unreachable!(),
        }
        n_iter += 1;
    }

    Ok(results)
}

fn parse_known_motifs(
    known_motifs_args: &KnownMotifsArgs,
    context_bases: [usize; 2],
    mod_code_lookup: &HashMap<ModCodeRepr, DnaBase>,
) -> anyhow::Result<Vec<EnrichedMotif>> {
    if known_motifs_args.known_motifs_table.is_none()
        && known_motifs_args.known_motifs.is_none()
    {
        bail!("must provide --known-motifs or --known-motifs-table")
    }

    let mut motifs_to_evaluate = Vec::new();
    if let Some(raw_known_motifs) = known_motifs_args.known_motifs.as_ref() {
        let command_line_motifs = parse_raw_known_motifs(
            raw_known_motifs,
            context_bases,
            mod_code_lookup,
        )?;
        motifs_to_evaluate.extend(command_line_motifs.into_iter());
    }

    if let Some(tab_fp) = known_motifs_args.known_motifs_table.as_ref() {
        let table_motifs =
            parse_motifs_from_table(tab_fp, context_bases, mod_code_lookup)?;
        debug!("parsed {} motifs from {tab_fp:?}", table_motifs.len());
        motifs_to_evaluate.extend(table_motifs.into_iter());
    }

    Ok(motifs_to_evaluate)
}

fn make_tables(motifs: &[EnrichedMotifData]) -> (pt::Table, pt::Table) {
    let human_header = row![
        "motif",
        "frac_mod",
        "high_count",
        "low_count",
        "mid_count",
        "log_odds"
    ];
    let mch_header = row![
        "mod_code",
        "motif",
        "offset",
        "frac_mod",
        "high_count",
        "low_count",
        "mid_count",
        "log_odds",
    ];
    let (mut human_table, mut machine_table) = motifs
        .iter()
        .sorted_by(|a, b| {
            b.frac_modified()
                .partial_cmp(&a.frac_modified())
                .unwrap_or(Ordering::Equal)
        })
        .fold(
            (pt::Table::new(), pt::Table::new()),
            |(mut hu, mut mch), next| {
                let human_row = row![
                    next.motif.to_string(),
                    next.frac_modified(),
                    next.total_high_count,
                    next.total_low_count,
                    next.total_mid_count,
                    next.log_odds(),
                ];
                let mach_row = row![
                    next.motif.multi_sequence.mod_code,
                    next.motif.format_seq(),
                    next.motif.multi_sequence.get_offset(),
                    next.frac_modified(),
                    next.total_high_count,
                    next.total_low_count,
                    next.total_mid_count,
                    next.log_odds(),
                ];
                hu.add_row(human_row);
                mch.add_row(mach_row);
                (hu, mch)
            },
        );
    human_table.set_titles(human_header);
    machine_table.set_titles(mch_header);
    (human_table, machine_table)
}

#[cfg(test)]
mod find_motifs_mod_tests {
    use std::collections::BTreeMap;

    use crate::mod_base_code::{DnaBase, ModCodeRepr, SIX_METHYL_ADENINE};
    use crate::motifs::iupac::IupacBase;
    use crate::motifs::{
        merge_motifs, EnrichedMotif, MotifRelationship, MultiSequence,
    };
    use bitvec::prelude::*;
    use common_macros::hash_map;

    fn make_enriched_motif(
        m: &[(i8, IupacBase)],
        canonical_base: DnaBase,
    ) -> EnrichedMotif {
        let mot = m.into_iter().copied().collect::<BTreeMap<_, _>>();
        EnrichedMotif::new(
            canonical_base,
            MultiSequence::new(ModCodeRepr::Code('m'), mot),
        )
    }

    #[test]
    fn test_multi_sequence() {
        use super::iupac::nt_bytes as nt;
        let kmers = vec![
            vec![nt::A, nt::A, nt::C, nt::G, nt::A],
            vec![nt::C, nt::A, nt::C, nt::G, nt::C],
            vec![nt::G, nt::A, nt::C, nt::G, nt::G],
            vec![nt::A, nt::A, nt::C, nt::G, nt::T],
        ];
        let multi_seq = MultiSequence::from_iter(
            kmers.iter().map(|x| x.as_slice()),
            [2, 2],
            ModCodeRepr::Code('m'),
        );
        assert_eq!("VA[m]G", &multi_seq.to_string());
        let seq = "ACTGACGATAC".chars().map(|x| x as u8).collect::<Vec<u8>>();
        let focus_position = 5usize;
        assert!(multi_seq.matches(&seq, focus_position));
        let seq = "AAGTACGATAA".chars().map(|x| x as u8).collect::<Vec<u8>>();
        assert!(!multi_seq.matches(&seq, focus_position));
    }

    #[test]
    fn test_multi_sequence_matches() {
        use super::iupac::nt_bytes as nt;
        let kmers = vec![
            vec![nt::A, nt::A, nt::C, nt::A, nt::T],
            vec![nt::C, nt::A, nt::C, nt::A, nt::T],
            vec![nt::G, nt::C, nt::C, nt::A, nt::T],
            vec![nt::T, nt::C, nt::C, nt::A, nt::T],
        ];
        let multi_seq = MultiSequence::from_iter(
            kmers.iter().map(|x| x.as_slice()),
            [2, 2],
            ModCodeRepr::ChEbi(21839),
        );

        assert_eq!(multi_seq.to_string(), "M[21839]AT".to_string());
        let context = "AAAAAAACCAAGCCTGGGCAAAAAC";
        assert!(!multi_seq.matches(context.as_bytes(), 12));
        let context = "AAAAAACGCCACCATGAGCGCATAC";
        assert!(multi_seq.matches(context.as_bytes(), 12));
    }

    #[test]
    fn test_contains_base() {
        use crate::motifs::iupac::nt_bytes;
        let seq = [(1i8, IupacBase::S), (2i8, IupacBase::G)]
            .into_iter()
            .collect::<BTreeMap<_, _>>();
        let ms = MultiSequence::new(ModCodeRepr::ChEbi(21839), seq);
        let motif = EnrichedMotif::new(DnaBase::C, ms);
        assert_eq!(motif.to_string(), "[21839]SG".to_string());
        assert!(motif.contains_base(&nt_bytes::G, 1))
    }

    #[test]
    fn test_create_motifs_from_fixed_length_kmers() {
        let kmers = vec!["ACCGG", "CCCGG", "TCCGG", "GCCGG"];
        let motif = MultiSequence::from_iter(
            kmers.iter().map(|x| x.as_bytes()),
            [2, 2],
            ModCodeRepr::ChEbi(21839),
        );
        assert_eq!(&motif.to_string(), "C[21839]GG");
    }

    #[test]
    fn test_is_superset() {
        let x = [(1i8, IupacBase::S), (2i8, IupacBase::G), (3i8, IupacBase::W)]
            .into_iter()
            .collect::<BTreeMap<_, _>>();
        let x = MultiSequence::new(ModCodeRepr::Code('m'), x);
        let y = [(1i8, IupacBase::G), (2i8, IupacBase::G), (3i8, IupacBase::W)]
            .into_iter()
            .collect::<BTreeMap<_, _>>();
        let y = MultiSequence::new(ModCodeRepr::Code('m'), y);
        assert!(x.is_superset(&y));
        let y = [(1i8, IupacBase::G), (2i8, IupacBase::G), (3i8, IupacBase::G)]
            .into_iter()
            .collect::<BTreeMap<_, _>>();
        let y = MultiSequence::new(ModCodeRepr::Code('m'), y);
        assert!(!x.is_superset(&y));

        let x = [(1i8, IupacBase::B), (2i8, IupacBase::D)]
            .into_iter()
            .collect::<BTreeMap<_, _>>();
        let x = MultiSequence::new(ModCodeRepr::Code('m'), x);
        let y = [(1i8, IupacBase::S), (2i8, IupacBase::W)]
            .into_iter()
            .collect::<BTreeMap<_, _>>();
        let y = MultiSequence::new(ModCodeRepr::Code('m'), y);
        assert!(x.is_superset(&y));
        assert!(!y.is_superset(&x));
    }

    #[test]
    fn test_is_subset() {
        // RNGA[21839]AY
        let x = {
            let mp = [
                (-4i8, IupacBase::R),
                (-2i8, IupacBase::G),
                (-1i8, IupacBase::A),
                (1i8, IupacBase::A),
                (2i8, IupacBase::Y),
            ]
            .into_iter()
            .collect::<BTreeMap<_, _>>();
            let ms = MultiSequence::new(ModCodeRepr::ChEbi(21839), mp);
            EnrichedMotif::new(DnaBase::C, ms)
        };
        assert_eq!(&x.to_string(), "RNGA[21839]AY");
        assert!(x.is_subset(&x.clone()));
        // GA[21839]AC
        let y = {
            let mp = [
                (-2i8, IupacBase::G),
                (-1i8, IupacBase::A),
                (1i8, IupacBase::A),
                (2i8, IupacBase::C),
            ]
            .into_iter()
            .collect::<BTreeMap<_, _>>();
            let ms = MultiSequence::new(ModCodeRepr::ChEbi(21839), mp);
            EnrichedMotif::new(DnaBase::C, ms)
        };
        assert_eq!(&y.to_string(), "GA[21839]AC");
        assert!(y.is_subset(&x));
    }

    #[test]
    fn test_merge_motif() {
        let dna_base = DnaBase::C;
        let w = make_enriched_motif(
            &[(1i8, IupacBase::A), (2i8, IupacBase::G), (3i8, IupacBase::T)],
            dna_base,
        );
        let x = make_enriched_motif(
            &[(1i8, IupacBase::S), (2i8, IupacBase::G), (3i8, IupacBase::W)],
            dna_base,
        );
        let y = make_enriched_motif(
            &[(1i8, IupacBase::G), (2i8, IupacBase::G), (3i8, IupacBase::A)],
            dna_base,
        );
        let z = make_enriched_motif(
            &[(1i8, IupacBase::C), (2i8, IupacBase::G), (3i8, IupacBase::T)],
            dna_base,
        );
        let merged = merge_motifs(vec![x, y, z, w]);
        assert_eq!(merged.len(), 2);
        assert!(merged.iter().find(|m| &m.to_string() == "[m]AGT").is_some());
        assert!(merged.iter().find(|m| &m.to_string() == "[m]SGW").is_some());
    }

    #[test]
    fn test_known_motifs() {
        let dna_base = DnaBase::C;
        let context_size = [3, 3];
        // G[]WSC
        let w = make_enriched_motif(
            &[
                (-1i8, IupacBase::G),
                (1i8, IupacBase::W),
                (2i8, IupacBase::S),
                (3i8, IupacBase::C),
            ],
            dna_base,
        );
        let bts = w.multi_sequence.to_bits([3, 3]);
        let expected = bits![
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1,
            0, 0
        ];
        assert_eq!(bts, expected);
        let t = make_enriched_motif(
            &[
                (-1i8, IupacBase::G),
                (1i8, IupacBase::A),
                (2i8, IupacBase::C),
                (3i8, IupacBase::C),
            ],
            dna_base,
        );

        let rel = w.compare(&t, context_size);
        assert_eq!(rel, MotifRelationship::Superset);

        let a = make_enriched_motif(
            &[(1i8, IupacBase::C), (2i8, IupacBase::G)],
            dna_base,
        );
        let b = make_enriched_motif(
            &[
                (-1i8, IupacBase::A),
                (1i8, IupacBase::W),
                (2i8, IupacBase::S),
                (3i8, IupacBase::C),
            ],
            dna_base,
        );
        let rel = w.compare(&b, context_size);
        assert_eq!(rel, MotifRelationship::Disjoint { edit_distance: 2 });
        let known_motifs = vec![a, b];
        let (closest, rel) = known_motifs
            .iter()
            .map(|m| w.compare(m, context_size))
            .enumerate()
            .min_by(|(_, a), (_, b)| a.cmp(b))
            .expect("should get min");
        assert_eq!(rel, MotifRelationship::Disjoint { edit_distance: 2 });
        assert_eq!(closest, 1usize);
    }

    #[test]
    fn test_motif_relationship_ord() {
        let equal = MotifRelationship::Equal;
        let subset = MotifRelationship::Subset;
        assert!(equal < subset);
        // 43210
        // 01234
        // GGCCAY
        // GGCCANNY
        let a = {
            let mp = [
                (-4i8, IupacBase::G),
                (-3i8, IupacBase::G),
                (-2i8, IupacBase::C),
                (-1i8, IupacBase::C),
                (1i8, IupacBase::Y),
            ]
            .into_iter()
            .collect::<BTreeMap<_, _>>();
            let ms = MultiSequence::new(ModCodeRepr::Code('a'), mp);
            EnrichedMotif::new(DnaBase::C, ms)
        };
        let b = {
            let mp = [
                (-4i8, IupacBase::G),
                (-3i8, IupacBase::G),
                (-2i8, IupacBase::C),
                (-1i8, IupacBase::C),
                (4i8, IupacBase::Y),
            ]
            .into_iter()
            .collect::<BTreeMap<_, _>>();
            let ms = MultiSequence::new(ModCodeRepr::Code('a'), mp);
            EnrichedMotif::new(DnaBase::C, ms)
        };
        assert_eq!(
            a.compare(&b, [4, 4]),
            MotifRelationship::Disjoint { edit_distance: 4 }
        );
        assert_eq!(
            b.compare(&a, [4, 4]),
            MotifRelationship::Disjoint { edit_distance: 4 }
        );
    }

    #[test]
    fn test_motif_subset_gh() {
        let codelookup = hash_map!(
            SIX_METHYL_ADENINE => DnaBase::A
        );
        let a = EnrichedMotif::new_from_parts(
            "GSATC",
            "a",
            "2",
            [12, 12],
            &codelookup,
        )
        .unwrap();
        let b = EnrichedMotif::new_from_parts(
            "GATC",
            "a",
            "1",
            [12, 12],
            &codelookup,
        )
        .unwrap();
        // dbg!(a.is_superset(&b));
        // dbg!(b.is_subset(&a));
        assert_eq!(
            a.compare(&b, [12, 12]),
            MotifRelationship::Disjoint { edit_distance: 2 }
        );
        // let x = vec![a, b];
        // dbg!(x.iter().map(|x| x.to_string()).collect::<Vec<String>>());
        // let x_merged = merge_motifs(x);
        // dbg!(x_merged.iter().map(|x|
        // x.to_string()).collect::<Vec<String>>());
    }
}
