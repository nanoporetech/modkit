use crate::mod_bam::{
    filter_records_iter, BaseModCall, BaseModProbs, CollapseMethod,
};
use crate::mod_base_code::DnaBase;
use crate::monoid::Moniod;
use indicatif::ParallelProgressIterator;
use log::{debug, error};
use rust_htslib::bam::{self};
use std::collections::HashMap;

use crate::reads_sampler::record_sampler::{Indicator, RecordSampler};
use crate::record_processor::{RecordProcessor, WithRecords};
use crate::util::{get_master_progress_bar, get_query_name_string, Strand};
use rayon::prelude::*;

/// Read IDs mapped to their base modification probabilities, organized
/// by the canonical base. This data structure contains essentially all
/// of the same data as in the records themselves, but with the query
/// position and the alternative probabilities removed (i.e. it only has
/// the probability of the called modification).
pub(crate) struct ReadIdsToBaseModProbs {
    // mapping of read id to canonical base mapped to a vec
    // of base mod calls on that canonical base
    pub(crate) inner: HashMap<String, HashMap<DnaBase, Vec<BaseModProbs>>>,
}

impl ReadIdsToBaseModProbs {
    fn add_read_without_probs(&mut self, read_id: &str) {
        self.inner
            .entry(read_id.to_owned())
            .or_insert(HashMap::new());
    }

    fn add_mod_probs_for_read(
        &mut self,
        read_id: &str,
        canonical_base: DnaBase,
        mod_probs: Vec<BaseModProbs>,
    ) {
        let read_id_entry = self
            .inner
            .entry(read_id.to_owned())
            .or_insert(HashMap::new());
        let added = read_id_entry.insert(canonical_base, mod_probs);
        if added.is_some() {
            error!(
                "double added base mod calls for base {} and read {},\
            potentially a logic error!",
                canonical_base.char(),
                read_id
            );
        }
    }

    #[inline]
    /// Returns most likely probabilities for base modifications predicted for
    /// each canonical base.
    pub(crate) fn mle_probs_per_base(&self) -> HashMap<DnaBase, Vec<f32>> {
        let pb = get_master_progress_bar(self.inner.len());
        pb.set_message("aggregating per-base modification probabilities");
        self.inner
            .par_iter()
            .progress_with(pb)
            .map(|(_, can_base_to_base_mod_probs)| {
                can_base_to_base_mod_probs
                    .iter()
                    .map(|(canonical_base, base_mod_probs)| {
                        let probs = base_mod_probs
                            .iter()
                            .filter_map(|bmc| match bmc.argmax_base_mod_call() {
                                Ok(BaseModCall::Modified(f, _)) => Some(f),
                                Ok(BaseModCall::Canonical(f)) => Some(f),
                                Ok(BaseModCall::Filtered) => {
                                    unreachable!("argmax base mod call should not return Filtered")
                                }
                                Err(e) => {
                                    debug!("{}", e.to_string());
                                    None
                                }
                            })
                            .collect::<Vec<f32>>();
                        (*canonical_base, probs)
                    })
                    .collect::<HashMap<DnaBase, Vec<f32>>>()
            })
            .reduce(|| HashMap::zero(), |a, b| a.op(b))
    }

    /// return argmax probs for each mod-code
    pub(crate) fn mle_probs_per_base_mod(&self) -> HashMap<char, Vec<f64>> {
        // todo(arand) should really aggregate per mod-code
        let pb = get_master_progress_bar(self.inner.len());
        pb.set_message("aggregating per-mod probabilities");
        self.inner
            .par_iter()
            .progress_with(pb)
            .filter_map(|(_, base_mod_probs)| {
                let grouped = base_mod_probs
                    .iter()
                    .map(|(base, base_mod_probs)| {
                        // let canonical_code = base
                        //     .canonical_mod_code()
                        //     .map(|c| c.char())
                        //     .unwrap_or(base.char());
                        base_mod_probs
                            .iter()
                            // can make this .base_mod_call
                            .filter_map(|bmc| {
                                match bmc.argmax_base_mod_call() {
                                    Ok(BaseModCall::Modified(p, code)) => {
                                        Some((code.char(), p as f64))
                                    }
                                    Ok(BaseModCall::Canonical(p)) => {
                                        Some((base.char(), p as f64))
                                    }
                                    Ok(BaseModCall::Filtered) => {
                                        unreachable!("argmax base mod call should not return Filtered")
                                    }
                                    Err(e) => {
                                        debug!("{}", e.to_string());
                                        None
                                    }
                                }
                            })
                            .fold(
                                HashMap::<char, Vec<f64>>::new(),
                                |mut acc, (base, p)| {
                                    acc.entry(base)
                                        .or_insert(Vec::new())
                                        .push(p);
                                    acc
                                },
                            )
                    })
                    .reduce(|a, b| a.op(b));
                grouped
            })
            .reduce(|| HashMap::zero(), |a, b| a.op(b))
    }
}

impl Moniod for ReadIdsToBaseModProbs {
    fn zero() -> Self {
        Self {
            inner: HashMap::new(),
        }
    }

    fn op(self, other: Self) -> Self {
        let mut acc = self.inner;
        for (read_id, base_mod_calls) in other.inner {
            if acc.contains_key(&read_id) {
                continue;
            } else {
                acc.insert(read_id, base_mod_calls);
            }
        }

        Self { inner: acc }
    }

    fn op_mut(&mut self, other: Self) {
        for (read_id, base_mod_calls) in other.inner {
            if self.inner.contains_key(&read_id) {
                continue;
            } else {
                self.inner.insert(read_id, base_mod_calls);
            }
        }
    }

    fn len(&self) -> usize {
        self.inner.len()
    }
}

impl RecordProcessor for ReadIdsToBaseModProbs {
    type Output = Self;

    fn process_records<T: bam::Read>(
        records: bam::Records<T>,
        with_progress: bool,
        mut record_sampler: RecordSampler,
        collapse_method: Option<&CollapseMethod>,
    ) -> anyhow::Result<Self::Output> {
        let spinner = if with_progress {
            Some(record_sampler.get_progress_bar())
        } else {
            None
        };
        let mod_base_info_iter = filter_records_iter(records);
        let mut read_ids_to_mod_base_probs = Self::zero();
        for (record, mod_base_info) in mod_base_info_iter {
            match record_sampler.ask() {
                Indicator::Use => {
                    let record_name = get_query_name_string(&record)
                        .unwrap_or("FAILED_UTF_DECODE".to_string());
                    if mod_base_info.is_empty() {
                        // add count of unused/no calls
                        read_ids_to_mod_base_probs
                            .add_read_without_probs(&record_name);
                        continue;
                    }

                    let (_, base_mod_probs_iter) =
                        mod_base_info.into_iter_base_mod_probs();
                    for (raw_canonical_base, strand, seq_pos_base_mod_probs) in
                        base_mod_probs_iter
                    {
                        let canonical_base = match (
                            DnaBase::parse(raw_canonical_base),
                            strand,
                        ) {
                            (Err(_), _) => continue,
                            (Ok(dna_base), Strand::Positive) => dna_base,
                            (Ok(dna_base), Strand::Negative) => {
                                dna_base.complement()
                            }
                        };
                        let mod_probs = seq_pos_base_mod_probs
                            .pos_to_base_mod_probs
                            .into_iter()
                            .map(|(_q_pos, base_mod_probs)| {
                                if let Some(method) = collapse_method {
                                    base_mod_probs.into_collapsed(method)
                                } else {
                                    base_mod_probs
                                }
                            })
                            .collect::<Vec<BaseModProbs>>();
                        read_ids_to_mod_base_probs.add_mod_probs_for_read(
                            &record_name,
                            canonical_base,
                            mod_probs,
                        );
                    }
                    if let Some(pb) = &spinner {
                        pb.inc(1);
                    }
                }
                Indicator::Skip => continue,
                Indicator::Done => break,
            }
        }

        if let Some(pb) = &spinner {
            pb.finish_and_clear();
        }

        Ok(read_ids_to_mod_base_probs)
    }
}

impl WithRecords for ReadIdsToBaseModProbs {
    fn size(&self) -> u64 {
        let s = self
            .inner
            .iter()
            .map(|(_, base_mod_calls)| {
                base_mod_calls.values().map(|vs| vs.len()).sum::<usize>()
            })
            .sum::<usize>();
        s as u64
    }

    fn num_reads(&self) -> usize {
        self.inner.len()
    }
}
