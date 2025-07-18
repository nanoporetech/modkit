use itertools::Itertools;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{self, Read};
use rustc_hash::FxHashMap;
use std::collections::BTreeMap;

use crate::mod_bam::{BaseModCall, ModBaseInfo, WithModBaseInfos};
use crate::mod_base_code::{DnaBase, ModCodeRepr};
use crate::pileup::PileupFeatureCounts;
use crate::threshold_mod_caller::MultipleThresholdModCaller;
use crate::util::{Strand, TAB};

type ChromId = u32;
type RefPosition = u32;
type PositionStrand = (RefPosition, Strand);
type ChromIdToPositionCounts =
    BTreeMap<ChromId, BTreeMap<PositionStrand, Vec<PileupFeatureCounts>>>;

enum Call {
    Delete,
    Filtered(DnaBase),
    NoCall(DnaBase),
    CanonicalCall {
        primary_base: DnaBase,
        mod_codes: Vec<ModCodeRepr>,
    },
    ModifiedCall {
        primary_base: DnaBase,
        mod_code: ModCodeRepr,
        mod_codes: Vec<ModCodeRepr>,
    },
}

pub(super) struct ModPileup {
    counts: ChromIdToPositionCounts,
}

impl ModPileup {
    pub(super) fn make_records(&self, header: &bam::HeaderView) -> Vec<String> {
        let mut records = Vec::new();
        for (chrom_id, pos_counts) in self.counts.iter() {
            let chrom = header
                .tid2name(*chrom_id)
                .iter()
                .map(|x| *x as char)
                .collect::<String>();
            for ((pos, _strand), pileup_feature_counts) in pos_counts.iter() {
                for counts in pileup_feature_counts {
                    let row = format!(
                        "{}{TAB}\
                         {}{TAB}\
                         {}{TAB}\
                         {}{TAB}\
                         {}{TAB}\
                         {}{TAB}\
                         {}{TAB}\
                         {}{TAB}\
                         {}{TAB}\
                         {}{TAB}\
                         {}{TAB}\
                         {}{TAB}\
                         {}{TAB}\
                         {}{TAB}\
                         {}{TAB}\
                         {}{TAB}\
                         {}{TAB}\
                         {}\n",
                            chrom,
                            pos,
                            pos + 1,
                            counts.raw_mod_code.to_string(),
                            counts.filtered_coverage,
                            counts.raw_strand,
                            pos,
                            pos + 1,
                            "255,0,0",
                            counts.filtered_coverage,
                            format!("{:.2}",counts.fraction_modified * 100f32),
                            counts.n_modified,
                            counts.n_canonical,
                            counts.n_other_modified,
                            counts.n_delete,
                            counts.n_filtered,
                            counts.n_diff,
                            counts.n_nocall,
                    );
                    records.push(row);
                }
            }
        }
        records
    }
}

#[derive(Debug, Hash, Eq, PartialEq)]
enum ModifiedBase {
    UnModified,
    Modified(ModCodeRepr),
}

#[derive(Debug, Default)]
pub(super) struct Tally2 {
    n_delete: u32,
    n_filtered: FxHashMap<DnaBase, u32>,
    basecall_counts: FxHashMap<DnaBase, u32>,
    modcall_counts: FxHashMap<DnaBase, FxHashMap<ModifiedBase, u32>>,
}

#[derive(Default)]
struct Tallies {
    inner: FxHashMap<ChromId, FxHashMap<PositionStrand, Tally2>>,
}

impl Tallies {
    fn add_features(
        &mut self,
        chrom_id: ChromId,
        strand: Strand,
        features: FxHashMap<RefPosition, Call>,
    ) {
        let chrom_tallies =
            self.inner.entry(chrom_id).or_insert_with(FxHashMap::default);
        features.into_iter().for_each(|(ref_pos, call_feature)| {
            let tally = chrom_tallies
                .entry((ref_pos, strand))
                .or_insert_with(Tally2::default);
            tally.add_call(call_feature)
        });
    }

    fn finalize(self, combine_mods: bool) -> ModPileup {
        let counts = self.inner
            .into_iter()
            // for every chrom
            .map(|(chrom_id, tallies)| {
                let counts = tallies
                    // for every strand/position
                    .into_iter()
                    .map(|((ref_pos, strand), tally)| {
                        let counts = tally.into_counts(combine_mods, strand);
                        ((ref_pos, strand), counts)
                    })
                    .collect::<BTreeMap<PositionStrand, Vec<PileupFeatureCounts>>>();
                (chrom_id, counts)
            })
            .collect::<BTreeMap<ChromId, BTreeMap<PositionStrand, Vec<PileupFeatureCounts>>>>();
        ModPileup { counts }
    }
}

impl Tally2 {
    fn add_call(&mut self, call: Call) {
        match call {
            Call::Delete => {
                self.n_delete = self.n_delete.saturating_add(1);
            }
            Call::Filtered(read_base) => {
                self.n_filtered
                    .entry(read_base)
                    .and_modify(|x| *x = x.saturating_add(1))
                    .or_insert(1);
            }
            Call::NoCall(read_base) => {
                self.basecall_counts
                    .entry(read_base)
                    .and_modify(|x| *x = x.saturating_add(1))
                    .or_insert(1);
            }
            Call::CanonicalCall { primary_base, mod_codes } => {
                self.modcall_counts
                    .entry(primary_base)
                    .or_insert_with(FxHashMap::default)
                    .entry(ModifiedBase::UnModified)
                    .and_modify(|x| *x = x.saturating_add(1))
                    .or_insert(1);
                for mod_code in mod_codes {
                    self.modcall_counts
                        .entry(primary_base)
                        .or_insert_with(FxHashMap::default)
                        .entry(ModifiedBase::Modified(mod_code))
                        .or_insert(0);
                }
            }
            Call::ModifiedCall { primary_base, mod_code, mod_codes } => {
                self.modcall_counts
                    .entry(primary_base)
                    .or_insert_with(FxHashMap::default)
                    .entry(ModifiedBase::Modified(mod_code))
                    .and_modify(|x| *x = x.saturating_add(1))
                    .or_insert(1);
                for mod_code in mod_codes {
                    self.modcall_counts
                        .entry(primary_base)
                        .or_insert_with(FxHashMap::default)
                        .entry(ModifiedBase::Modified(mod_code))
                        .or_insert(0);
                }
            }
        }
    }

    fn into_counts(
        self,
        combine_mods: bool,
        strand: Strand,
    ) -> Vec<PileupFeatureCounts> {
        self.modcall_counts
            .into_iter()
            .flat_map(|(base, calls)| {
                let n_modified_all = calls
                    .iter()
                    .filter_map(|(mb, c)| match mb {
                        ModifiedBase::UnModified => None,
                        ModifiedBase::Modified(_) => Some(*c),
                    })
                    .sum::<u32>();

                let (n_diff, n_nocall) = self.basecall_counts.iter().fold(
                    (0u32, 0u32),
                    |(n_diff, n_nocall), (b, c)| {
                        if *b == base {
                            (n_diff, n_nocall.saturating_add(*c))
                        } else {
                            (n_diff.saturating_add(*c), n_nocall)
                        }
                    },
                );

                let n_filt =
                    self.n_filtered.get(&base).copied().unwrap_or(0u32);
                let n_canonical =
                    calls.get(&ModifiedBase::UnModified).copied().unwrap_or(0);
                let filtered_coverage = calls.values().sum::<u32>();
                assert!(filtered_coverage > 0);
                let calls = if combine_mods {
                    let mut tmp = FxHashMap::default();
                    tmp.insert(
                        ModifiedBase::Modified(ModCodeRepr::Code(base.char())),
                        n_modified_all,
                    );
                    tmp
                } else {
                    calls
                };
                let counts = calls
                    .into_iter()
                    .filter_map(|(mod_base, count)| match mod_base {
                        ModifiedBase::UnModified => None,
                        ModifiedBase::Modified(mod_code) => {
                            Some(PileupFeatureCounts::new(
                                strand.to_char(),
                                filtered_coverage,
                                mod_code,
                                count as f32 / filtered_coverage as f32,
                                n_canonical,
                                count,
                                n_modified_all.saturating_sub(count),
                                self.n_delete,
                                n_filt,
                                n_diff,
                                n_nocall,
                                None,
                            ))
                        }
                    })
                    .sorted_by(|a, b| a.raw_mod_code.cmp(&b.raw_mod_code))
                    .collect::<Vec<PileupFeatureCounts>>();
                counts
            })
            .collect()
    }
}

fn read_info_to_features(
    mod_base_info: &ModBaseInfo,
    caller: &MultipleThresholdModCaller,
    record: &bam::Record,
) -> FxHashMap<RefPosition, Call> {
    let seq_length = record.seq_len();
    let read_pos_to_calls = mod_base_info
        .pos_seq_base_mod_probs
        .iter()
        .flat_map(|(dna_base, pos_probs)| {
            pos_probs.pos_to_base_mod_probs.keys().map(|pos| (*pos, *dna_base))
        })
        .fold(FxHashMap::default(), |mut acc, (pos, base)| {
            let seen = acc.insert(pos, base).is_some();
            assert!(!seen);
            acc
        });

    let features = record
        .aligned_pairs_full()
        .filter_map(|pair| {
            let q_pos = pair[0];
            let r_pos = pair[1];
            match (q_pos, r_pos) {
                (Some(q_pos), Some(r_pos)) => {
                    if r_pos < 0 {
                        None
                    } else {
                        let r_pos = r_pos as u32;
                        assert!(q_pos >= 0);
                        let q_pos = q_pos as usize;
                        let forward_position = if record.is_reverse() {
                            seq_length
                                .checked_sub(q_pos)
                                .and_then(|x| x.checked_sub(1))
                                .expect("invalid aligned pair")
                        } else {
                            q_pos
                        };
                        let call_feature = if let Some(primary_base) =
                            read_pos_to_calls.get(&forward_position)
                        {
                            let threshold_base = if record.is_reverse() {
                                primary_base.complement()
                            } else {
                                *primary_base
                            };
                            let base_mod_probs = mod_base_info
                                .pos_seq_base_mod_probs
                                .get(primary_base)
                                .and_then(|x| {
                                    x.pos_to_base_mod_probs
                                        .get(&forward_position)
                                })
                                .unwrap();
                            let base_mod_call =
                                caller.call(&threshold_base, base_mod_probs);
                            let mod_codes = base_mod_probs
                                .iter_probs()
                                .map(|(code, _)| *code)
                                .filter(|c| match base_mod_call {
                                    BaseModCall::Canonical(_)
                                    | BaseModCall::Filtered => true,
                                    BaseModCall::Modified(_, code) => {
                                        code != *c
                                    }
                                })
                                .collect::<Vec<ModCodeRepr>>();
                            let call_feature = match base_mod_call {
                                BaseModCall::Canonical(_) => {
                                    Call::CanonicalCall {
                                        primary_base: *primary_base,
                                        mod_codes,
                                    }
                                }
                                BaseModCall::Modified(_, mod_code) => {
                                    Call::ModifiedCall {
                                        primary_base: *primary_base,
                                        mod_code,
                                        mod_codes,
                                    }
                                }
                                BaseModCall::Filtered => {
                                    Call::Filtered(*primary_base)
                                }
                            };
                            call_feature
                        } else {
                            let read_base = {
                                let b =
                                    DnaBase::parse(record.seq()[q_pos] as char)
                                        .unwrap();
                                if record.is_reverse() {
                                    b.complement()
                                } else {
                                    b
                                }
                            };
                            Call::NoCall(read_base)
                        };

                        Some((r_pos, call_feature))
                    }
                }
                (None, Some(r_pos)) => {
                    assert!(r_pos >= 0);
                    let call = Call::Delete;
                    Some((r_pos as u32, call))
                }
                _ => None,
            }
        })
        .fold(FxHashMap::default(), |mut acc, (r_pos, call)| {
            let seen = acc.insert(r_pos, call).is_some();
            assert!(!seen);
            acc
        });
    features
}

pub(super) fn process_stream(
    mut reader: bam::Reader,
    caller: &MultipleThresholdModCaller,
    combine_mods: bool,
) -> anyhow::Result<ModPileup> {
    let modbase_infos = reader.records().with_mod_base_info().filter_map(
        |(record, modbase_info)| {
            if record.tid() < 0i32 {
                None
            } else {
                Some((record, modbase_info))
            }
        },
    );

    let features = modbase_infos.map(|(record, modbase_info)| {
        let features = read_info_to_features(&modbase_info, caller, &record);
        let chrom_id = record.tid() as u32;
        let strand = if record.is_reverse() {
            Strand::Negative
        } else {
            Strand::Positive
        };
        ((chrom_id, strand), features)
    });

    let tallies = features.fold(
        Tallies::default(),
        |mut tallies, ((chrom_id, strand), features)| {
            tallies.add_features(chrom_id, strand, features);
            tallies
        },
    );

    let pileup = tallies.finalize(combine_mods);
    Ok(pileup)
}
