use std::{
    collections::HashMap, marker::PhantomData, path::PathBuf, sync::Arc, usize,
};

use anyhow::{anyhow, bail};
use bitvec::{order::Lsb0, view::BitView};
use derive_new::new;
use indicatif::MultiProgress;
use itertools::Itertools;
use log::{debug, info, warn};
use log_once::{debug_once, warn_once};
use prettytable::row;
use rand::{rngs::SmallRng, Rng, SeedableRng};
use rust_htslib::{
    bam::{self, ext::BamRecordExtensions, FetchDefinition, HeaderView, Read},
    htslib,
};
use rustc_hash::FxHashMap;
use tracing::error;

use crate::{
    command_utils::get_motif_lookup_from_parts,
    errs::MkResult,
    interval_chunks::{
        ChromCoordinates, ChromCoordinatesFeeder, FocusPositions2, TotalLength,
    },
    mod_bam::{prob_to_qual, validate_mn_tag_on_record, EdgeFilter},
    mod_base_code::ModCodeRepr,
    modbam_util::subcommands::SampleModBaseProbs,
    motifs::motif_bed::RegexMotif,
    pileup::base_mods_adapter::{
        AlignedBaseModsIterator, BaseModsAdapter, ModState,
    },
    position_filter::StrandedPositionFilter,
    util::{
        get_human_readable_table, get_master_progress_bar, get_targets,
        get_ticker, record_is_not_primary, CheckedAddArr, Region, Strand,
    },
};
use crate::{
    mod_base_code::DnaBase,
    util::{qual_to_prob, reader_is_cram},
};

#[derive(new, Debug)]
pub(crate) struct PercentileThresholdQual {
    pub percentile: f32,
    pub threshold: f32,
    pub qual: u8,
}

#[derive(Debug, Copy, Clone)]
pub(crate) struct ModHist {
    pub total: u64,
    pub mod_code: ModCodeRepr,
    pub dna_base: DnaBase,
    pub hist: [u64; 256],
}

impl ModHist {
    fn new_one(mod_state: &ModState) -> Self {
        let total = 1u64;
        let mut hist = [0u64; 256];
        hist[mod_state.mod_qual as usize] = 1u64;
        Self {
            total,
            mod_code: mod_state.mod_code,
            dna_base: mod_state.primary_base,
            hist,
        }
    }

    fn clear(&mut self) {
        self.total = 0u64;
        self.hist.iter_mut().for_each(|x| *x = 0u64);
    }

    fn approx_checked_add(&mut self, other: &Self) -> Result<(), ()> {
        assert_eq!(
            self.mod_code, other.mod_code,
            "should not add two ModHists with different mod codes"
        );
        assert_eq!(
            self.dna_base, other.dna_base,
            "should not add two ModHists with different dna bases"
        );
        self.hist.approx_checked_add(&other.hist, &mut self.total, other.total)
    }

    // fn check(&self) -> anyhow::Result<()> {
    //     let hist_total = self.hist.iter().sum::<u64>();
    //     if hist_total == self.total {
    //         Ok(())
    //     } else {
    //         bail!("mod hist for {}:{} invalid", self.dna_base, self.mod_code)
    //     }
    // }
}

#[derive(Debug)]
pub(crate) struct QualHist {
    pub hist: [[u64; 256]; 4],
    pub base_totals: [u64; 4],
    pub num_records_with_base_mods: [u32; 4],
    pub explicit_canonical_probs: [u8; 4],
    pub mods_hists: Vec<ModHist>,
    pub interval_length: u32,
    pub erred_records: usize,
    pub ok_records: usize,
}

impl QualHist {
    pub(crate) fn clear(&mut self) {
        for ar in self.hist.iter_mut() {
            ar.iter_mut().for_each(|x| *x = 0u64);
        }
        self.explicit_canonical_probs.iter_mut().for_each(|x| *x = 0u8);
        self.base_totals.iter_mut().for_each(|x| *x = 0u64);
        self.mods_hists.iter_mut().for_each(|h| h.clear());
        self.num_records_with_base_mods.iter_mut().for_each(|x| *x = 0u32);
        self.erred_records = 0usize;
        self.ok_records = 0usize;
        self.interval_length = 0u32;
        // self.check_all("clear").expect("should work at clear");
    }

    pub(crate) fn combine(
        &mut self,
        other: &Self,
        multi_progress: &MultiProgress,
    ) -> anyhow::Result<()> {
        // self.check_all("combine start")?;
        // other.check_all("combine start - other")?;
        let Some(count) = self.erred_records.checked_add(other.erred_records)
        else {
            bail!("over-sampled failed records, subsample bam")
        };
        self.erred_records = count;
        self.ok_records = self.ok_records.saturating_add(other.ok_records);
        for (i, (ag, t)) in self
            .base_totals
            .iter_mut()
            .zip(other.base_totals.iter())
            .enumerate()
        {
            let sat = self.hist[i].approx_checked_add(&other.hist[i], ag, *t);
            if sat.is_err() {
                multi_progress.suspend(|| {
                    let b = DnaBase::try_from(i).unwrap();
                    warn_once!(
                        "Saturated counts for base {b}, counts will be \
                         approximate"
                    );
                })
            }
        }
        // self.check_all("combine middle")?;
        for (x, y) in self
            .explicit_canonical_probs
            .iter_mut()
            .zip(other.explicit_canonical_probs.iter())
        {
            *x = std::cmp::max(*x, *y);
        }
        for (x, y) in self
            .num_records_with_base_mods
            .iter_mut()
            .zip(other.num_records_with_base_mods)
        {
            *x = x.saturating_add(y);
        }

        for mod_hist in other.mods_hists.iter() {
            let this_idx = self
                .mods_hists
                .iter()
                .enumerate()
                .find(|(_, mh)| {
                    mh.mod_code == mod_hist.mod_code
                        && mh.dna_base == mod_hist.dna_base
                })
                .map(|(i, _)| i);
            if let Some(this_idx) = this_idx {
                let sat =
                    self.mods_hists[this_idx].approx_checked_add(&mod_hist);
                if sat.is_err() {
                    multi_progress.suspend(|| {
                        warn_once!(
                            "Saturated counts for modified base code {}, \
                             counts will be approximate",
                            mod_hist.mod_code
                        );
                    })
                }
            } else {
                self.mods_hists.push(*mod_hist);
            }
        }
        // self.check_all("combine end")
        Ok(())
    }

    fn add_bases_from_record(&mut self, bases: u8) {
        for (i, _) in
            bases.view_bits::<Lsb0>().iter().enumerate().filter(|(_, x)| **x)
        {
            let count = self.num_records_with_base_mods[i];
            self.num_records_with_base_mods[i] = count.saturating_add(1u32);
        }
    }

    fn calc_base_percentile_ranks(
        desired_percentiles: &[f32],
        total: u64,
    ) -> Vec<u64> {
        desired_percentiles
            .iter()
            .map(|&x| {
                let rank = if x <= 0f32 {
                    1u64
                } else {
                    ((x as f64 * total as f64) / 100f64).ceil() as u64
                };
                rank.min(total)
            })
            .collect()
    }

    fn calc_percentile_values(
        percentile_ranks: &[u64],
        hist: &[u64; 256],
        total: u64,
    ) -> Vec<Option<u8>> {
        let mut cum_sum = [0u64; 256];
        let mut sum = 0u64;
        for i in 0usize..256 {
            sum = sum.saturating_add(hist[i]);
            cum_sum[i] = sum;
        }
        assert_eq!(cum_sum[255], total, "{} =/= {total}", cum_sum[255]);
        let mut out = vec![None; percentile_ranks.len()];
        let mut qi = 0usize;

        for value in 0u16..256 {
            let cum = cum_sum[value as usize];
            while qi < percentile_ranks.len() && percentile_ranks[qi] <= cum {
                out[qi] = Some(value as u8);
                qi = qi.saturating_add(1);
            }

            if qi >= percentile_ranks.len() {
                break;
            }
        }

        out
    }

    pub(crate) fn percentiles(
        hist: &[[u64; 256]; 4],
        base_totals: &[u64; 4],
        desired_percentiles: &[f32],
    ) -> HashMap<DnaBase, Vec<PercentileThresholdQual>> {
        // debug!("{desired_percentiles:?}");
        let mut percentiles = HashMap::new();
        // debug!("base_totals={base_totals:?}");
        for b in 0usize..4 {
            if base_totals[b] > 0 {
                let dna_base = DnaBase::try_from(b).unwrap();
                // debug!("calculating for {dna_base} ({b})");
                let percentile_ranks = Self::calc_base_percentile_ranks(
                    desired_percentiles,
                    base_totals[b],
                );
                // debug!("percentile_ranks={percentile_ranks:?}");
                let percentile_values = Self::calc_percentile_values(
                    &percentile_ranks,
                    &hist[b],
                    base_totals[b],
                );
                // debug!("percentile_value={percentile_values:?}");
                let base_percentiles = desired_percentiles
                    .iter()
                    .zip(percentile_values)
                    .map(|(x, y)| {
                        let qual = y.unwrap();
                        PercentileThresholdQual::new(
                            *x,
                            qual_to_prob(qual),
                            qual,
                        )
                    })
                    .collect::<Vec<PercentileThresholdQual>>();
                percentiles.insert(dna_base, base_percentiles);
            }
        }
        percentiles
    }

    pub(crate) fn base_level_percentiles(
        &self,
        desired_percentiles: &[f32],
        multi_progress: &MultiProgress,
        verbose: bool,
    ) -> HashMap<DnaBase, Vec<PercentileThresholdQual>> {
        // self.check_all("base level percentiles")
        //     .expect("shold work at base level percentiles");
        let mut agg = Self::default();
        agg.base_totals = self.base_totals;
        agg.hist = self.hist;
        for mod_hist in self.mods_hists.iter() {
            // mod_hist.check().expect("mod_hist should be valid");
            let primary_base = mod_hist.dna_base;
            let mod_total = mod_hist.total;
            let sat = agg.hist[primary_base as usize].approx_checked_add(
                &mod_hist.hist,
                &mut agg.base_totals[primary_base as usize],
                mod_total,
            );
            if sat.is_err() {
                multi_progress.suspend(|| {
                    warn_once!(
                        "Saturated counts for modified base base {} for {} \
                         counts will be approximate",
                        mod_hist.dna_base,
                        mod_hist.mod_code
                    );
                })
            }
        }
        let mut total_counts_table = get_human_readable_table();
        let mut show_table = false;
        total_counts_table
            .set_titles(row!["primary_base", "probabilities_count"]);
        for base in (0..4usize).map(|x| DnaBase::try_from(x).unwrap()) {
            let total = agg.base_totals[base as usize];
            if total > 0 {
                total_counts_table
                    .add_row(row![base, agg.base_totals[base as usize]]);
                show_table = true
            }
        }
        if show_table && verbose {
            multi_progress.suspend(|| {
                info!(
                    "Total number of probabilities collected per primary \
                     sequence base:\n{total_counts_table}"
                );
            });
        } else {
            debug!(
                "Total number of probabilities collected per primary sequence \
                 base:\n{total_counts_table}"
            );
        }

        Self::percentiles(&agg.hist, &agg.base_totals, desired_percentiles)
    }

    pub(crate) fn combine_mods(&mut self, multi_progress: &MultiProgress) {
        let mut mod_hists = (0..4usize)
            .map(|x| {
                let primary_base = DnaBase::try_from(x).unwrap();
                ModHist {
                    total: 0,
                    mod_code: ModCodeRepr::Code(primary_base.char()),
                    dna_base: primary_base,
                    hist: [0u64; 256],
                }
            })
            .collect::<Vec<ModHist>>();
        for mod_hist in self.mods_hists.iter() {
            let idx = mod_hist.dna_base as usize;
            let agg = &mut mod_hists[idx];
            let sat = agg.hist.approx_checked_add(
                &mod_hist.hist,
                &mut agg.total,
                mod_hist.total,
            );
            if sat.is_err() {
                multi_progress.suspend(|| {
                    warn_once!(
                        "Saturated counts for base modifications to primary \
                         base {}, counts will be approximate",
                        mod_hist.dna_base
                    );
                })
            }
        }
        mod_hists.retain(|x| x.total > 0);
        let _ = std::mem::replace(&mut self.mods_hists, mod_hists);
    }

    pub(crate) fn from_records<R: bam::Read>(
        records: bam::Records<R>,
        stranded_position_filter: Option<StrandedPositionFilter<()>>,
        num_reads: Option<usize>,
        sampling_frac: Option<f64>,
        seed: Option<u64>,
        edge_filter: Option<&EdgeFilter>,
        collect_mod_probs: bool,
        allow_non_primary: bool,
        mapped_only: bool,
        multi_progress: &MultiProgress,
    ) -> anyhow::Result<Self> {
        let edge_filter_start =
            edge_filter.map(|ef| ef.edge_filter_start).unwrap_or(0usize);
        let edge_filter_end =
            edge_filter.map(|ef| ef.edge_filter_end).unwrap_or(0usize);
        let record_counter = multi_progress.add(get_ticker());
        let mut rng = SmallRng::seed_from_u64(seed.unwrap_or(42u64));
        record_counter.set_message("records processed");
        let mut mem = Self::default();
        'records: for res in records {
            if let Some(num_records) = num_reads {
                if mem.ok_records >= num_records {
                    // debug!(
                    //     "sampled {} records, done. {} records failed",
                    //     mem.ok_records, mem.erred_records
                    // );
                    return Ok(mem);
                }
            }
            let Ok(record) = res else {
                mem.erred_records = mem.erred_records.saturating_add(1);
                continue 'records;
            };
            if record_is_not_primary(&record) {
                if !allow_non_primary {
                    continue 'records;
                }
                if validate_mn_tag_on_record(&record).is_err() {
                    mem.erred_records = mem.erred_records.saturating_add(1);
                    continue 'records;
                }
            }
            if mapped_only && record.is_unmapped() {
                continue 'records;
            }
            let read_length = record.seq_len();
            let (left_edge_filter, right_edge_filter) = if record.is_reverse() {
                (edge_filter_end, read_length.saturating_sub(edge_filter_start))
            } else {
                (edge_filter_start, read_length.saturating_sub(edge_filter_end))
            };
            if let Some(sampling_frac) = sampling_frac {
                if !rng.gen_bool(sampling_frac) {
                    continue 'records;
                }
            }
            if let Some(spf) = stranded_position_filter.as_ref() {
                if record.is_unmapped() {
                    continue 'records;
                }
                let chrom_id = record.tid();
                let strand = if record.is_reverse() {
                    Strand::Negative
                } else {
                    Strand::Positive
                };
                if chrom_id < 0 {
                    continue 'records;
                }
                let Ok(aligned_pairs_iter) =
                    AlignedBaseModsIterator::new(&record)
                else {
                    mem.erred_records = mem.erred_records.saturating_add(1);
                    continue 'records;
                };

                let bases_in_record =
                    aligned_pairs_iter.primary_bases_in_record();
                mem.add_bases_from_record(bases_in_record);

                let aligned_pairs_iter = aligned_pairs_iter
                    .into_iter()
                    .filter_ok(move |(qpos, _rpos, _mod_state)| {
                        (*qpos as usize) >= left_edge_filter
                            && (*qpos as usize) < right_edge_filter
                    })
                    .filter_ok(|(_, rpos, _)| {
                        spf.contains(chrom_id, *rpos as u64, strand)
                    });

                for res in aligned_pairs_iter {
                    match res {
                        Ok((_, _, mod_state)) => {
                            if collect_mod_probs {
                                increment_mods_counts(
                                    mod_state,
                                    &mut rng,
                                    &mut mem.explicit_canonical_probs,
                                    &mut mem.hist,
                                    &mut mem.base_totals,
                                    &mut mem.mods_hists,
                                );
                            } else {
                                increment_base_counts(
                                    mod_state,
                                    &mut mem.explicit_canonical_probs,
                                    &mut mem.hist,
                                    &mut mem.base_totals,
                                    &mut rng,
                                );
                            }
                        }
                        Err(_) => {
                            mem.erred_records =
                                mem.erred_records.saturating_add(1);
                            continue 'records;
                        }
                    }
                }
                mem.ok_records = mem.ok_records.saturating_add(1);
                record_counter.inc(1);
            } else {
                let Ok(mut modbase_iter) = BaseModsAdapter::<16>::new(&record)
                else {
                    mem.erred_records = mem.erred_records.saturating_add(1);
                    continue 'records;
                };
                mem.add_bases_from_record(
                    modbase_iter.primary_bases_in_record(),
                );
                'mods: loop {
                    match modbase_iter.next_modified_position_no_thresh() {
                        Ok(Some(mod_state)) => {
                            if mod_state.mod_position >= left_edge_filter
                                && mod_state.mod_position < right_edge_filter
                            {
                                if collect_mod_probs {
                                    increment_mods_counts(
                                        mod_state,
                                        &mut rng,
                                        &mut mem.explicit_canonical_probs,
                                        &mut mem.hist,
                                        &mut mem.base_totals,
                                        &mut mem.mods_hists,
                                    );
                                } else {
                                    increment_base_counts(
                                        mod_state,
                                        &mut mem.explicit_canonical_probs,
                                        &mut mem.hist,
                                        &mut mem.base_totals,
                                        &mut rng,
                                    );
                                }
                            }
                        }
                        Ok(None) => break 'mods,
                        Err(_) => {
                            mem.erred_records =
                                mem.erred_records.saturating_add(1);
                            continue 'records;
                        }
                    }
                }
                mem.ok_records = mem.ok_records.saturating_add(1);
                record_counter.inc(1);
            }
        }
        Ok(mem)
    }

    pub(crate) fn get_base_thresholds(
        &self,
        filter_percentile: f32,
        max_thresholds_per_base: Option<[f32; 4]>,
        multi_progress: &MultiProgress,
    ) -> anyhow::Result<[f32; 4]> {
        if self.ok_records == 0 {
            bail!("Failed to sample any records to estimate threshold.")
        }
        let mut base_thresholds = [0f32; 4];
        let filter_percentile = if filter_percentile >= 1.0f32 {
            100f32
        } else if filter_percentile < 0f32 {
            bail!("filter percentile must be positive or zero")
        } else {
            filter_percentile * 100f32
        };
        for (base, vals) in QualHist::percentiles(
            &self.hist,
            &self.base_totals,
            &vec![filter_percentile],
        ) {
            assert_eq!(vals.len(), 1);
            let (t, q) = (vals[0].threshold, vals[0].qual);
            let user_specified_max_threshold =
                max_thresholds_per_base.map(|x| x[base as usize]);
            if q == 255u8
                || user_specified_max_threshold
                    .map(|x| qual_to_prob(q) > x)
                    .unwrap_or(false)
            {
                if let Some(user_specified_threshold) =
                    user_specified_max_threshold
                {
                    let user_specified_threshold_qual =
                        prob_to_qual(user_specified_threshold);
                    multi_progress.suspend(|| {
                        warn!(
                            "estimated threshold ({t}, qual: {q}) for {base} \
                             higher than specified threshold, \
                             {user_specified_threshold}. Setting threshold to \
                             {user_specified_threshold} (qual: \
                             {user_specified_threshold_qual})",
                        );
                    });
                    base_thresholds[base as usize] = user_specified_threshold;
                } else {
                    let q_fb = self.explicit_canonical_probs[base as usize];
                    let t_fb = qual_to_prob(q_fb);
                    multi_progress.suspend(|| {
                        warn!(
                            "estimated threshold for {base} too high ({t}, \
                             qual: {q}), setting threshold {t_fb} (qual: \
                             {q_fb}), highest discovered explicit canonical \
                             value",
                        );
                    });
                    base_thresholds[base as usize] = t_fb;
                }
            } else {
                multi_progress.suspend(|| {
                    info!(
                        "setting threshold {t} (qual: {q}) for base {base}, \
                         percentile {}, {} total probabilites",
                        filter_percentile, self.base_totals[base as usize]
                    );
                });
                base_thresholds[base as usize] = t;
            }
        }
        Ok(base_thresholds)
    }

    // fn check_all(&self, place: &str) -> anyhow::Result<()> {
    //     for i in 0..4 {
    //         let hist_total = self.hist[i].iter().sum::<u64>();
    //         assert_eq!(hist_total, self.base_totals[i]);
    //         if hist_total != self.base_totals[i] {
    //             let dna_base = DnaBase::try_from(i).unwrap();
    //             bail!("invalid dna counts for {dna_base} at {place}");
    //         }
    //     }
    //     for mod_hist in self.mods_hists.iter() {
    //         mod_hist.check()?;
    //     }
    //     debug!("check all good at {place}");
    //     Ok(())
    // }
}

impl Default for QualHist {
    fn default() -> Self {
        Self {
            hist: [[0u64; 256]; 4],
            explicit_canonical_probs: [0u8; 4],
            base_totals: Default::default(),
            num_records_with_base_mods: [0u32; 4],
            mods_hists: Default::default(),
            erred_records: Default::default(),
            interval_length: Default::default(),
            ok_records: Default::default(),
        }
    }
}

pub(crate) trait ExtractsMleProbs<T> {
    fn new(
        seed: u64,
        sample: bool,
        sample_frac: f64,
        motif_bases: [DnaBase; 4],
        edge_filter: Option<&EdgeFilter>,
    ) -> Self;
    fn process_record(
        &mut self,
        record: &bam::Record,
        chrom_corrds: &ChromCoordinates,
        explicit_canonical_probs: &mut [u8; 4],
        base_hist: &mut [[u64; 256]; 4],
        base_totals: &mut [u64; 4],
        mods_hists: &mut Vec<ModHist>,
        num_records_with_base_mods: &mut [u32; 4],
    ) -> anyhow::Result<bool>;
}

pub(crate) trait ExtractProbsWorker
where
    Self: Send,
{
    fn process(
        &mut self,
        item: ChromCoordinates,
        mem: QualHist,
    ) -> anyhow::Result<QualHist>;
}

pub(crate) struct AlignedBaseAndModArgmaxProbs {}
pub(crate) struct AlignedBaseArgmaxProbs {}
pub(crate) struct BaseArgmaxProbs {}
pub(crate) struct BaseAndModArgmaxProbs {}

pub(crate) struct RegionMleProbs<T, H> {
    reader: bam::IndexedReader,
    allow_non_primary: bool,
    hist: H,
    t: PhantomData<T>,
    chrom_to_counts: Option<Arc<FxHashMap<u32, usize>>>,
}

impl<T: Send + Sync, H: ExtractsMleProbs<T>> RegionMleProbs<T, H> {
    pub(crate) fn new(
        bam_fp: &PathBuf,
        reference_fp: Option<&PathBuf>,
        allow_non_primary: bool,
        motif_bases: [DnaBase; 4],
        edge_filter: Option<&EdgeFilter>,
        thread_pool: &rust_htslib::tpool::ThreadPool,
        worker_no: u64,
        sample: bool,
        sample_frac: f64,
        chrom_to_counts: Option<Arc<FxHashMap<u32, usize>>>,
    ) -> anyhow::Result<Self> {
        let mut reader = bam::IndexedReader::from_path(bam_fp)?;
        reader.set_thread_pool(&thread_pool)?;
        if reader_is_cram(&reader) {
            if let Some(reference_fp) = reference_fp {
                reader.set_reference(reference_fp)?;
                reader.set_cram_options(
                    htslib::hts_fmt_option_CRAM_OPT_REQUIRED_FIELDS,
                    htslib::sam_fields_SAM_FLAG
                        | htslib::sam_fields_SAM_RNAME
                        | htslib::sam_fields_SAM_POS
                        | htslib::sam_fields_SAM_MAPQ
                        | htslib::sam_fields_SAM_CIGAR
                        | htslib::sam_fields_SAM_SEQ
                        | htslib::sam_fields_SAM_AUX,
                )?;
            } else {
                bail!("CRAM input requires reference")
            }
        };
        let hist =
            H::new(worker_no, sample, sample_frac, motif_bases, edge_filter);
        Ok(Self {
            reader,
            allow_non_primary,
            hist,
            t: PhantomData::<T>,
            chrom_to_counts,
        })
    }
}

impl<T: Send + Sync, H: ExtractsMleProbs<T> + Send> ExtractProbsWorker
    for RegionMleProbs<T, H>
{
    fn process(
        &mut self,
        item: ChromCoordinates,
        mut mem: QualHist,
    ) -> anyhow::Result<QualHist> {
        let chrom_tid = item.chrom_tid;
        let start_pos = item.start_pos;
        let end_pos = item.end_pos;
        mem.interval_length = item.len();
        // debug!("worker starting on {chrom_tid}:{start_pos}-{end_pos}");

        self.reader.fetch(FetchDefinition::Region(
            chrom_tid as i32,
            start_pos as i64,
            end_pos as i64,
        ))?;

        let num_records =
            if let Some(chrom_to_counts) = self.chrom_to_counts.as_ref() {
                Some(*chrom_to_counts.get(&item.chrom_tid).unwrap_or(&0usize))
            } else {
                None
            };

        'records: for res in self.reader.records() {
            if let Some(num_records) = num_records {
                if mem.ok_records >= num_records {
                    // debug!(
                    //     "worker {chrom_tid}:{start_pos}-{end_pos} sampled {}
                    // \      records, done. {} records
                    // failed",     mem.ok_records,
                    // mem.erred_records );
                    return Ok(mem);
                }
            }
            let Ok(record) = res else {
                mem.erred_records = mem.erred_records.saturating_add(1);
                continue 'records;
            };

            if record_is_not_primary(&record) {
                if !self.allow_non_primary {
                    continue 'records;
                }
                if validate_mn_tag_on_record(&record).is_err() {
                    mem.erred_records = mem.erred_records.saturating_add(1);
                    continue 'records;
                }
            }
            match self.hist.process_record(
                &record,
                &item,
                &mut mem.explicit_canonical_probs,
                &mut mem.hist,
                &mut mem.base_totals,
                &mut mem.mods_hists,
                &mut mem.num_records_with_base_mods,
            ) {
                Ok(true) => {
                    mem.ok_records = mem.ok_records.saturating_add(1);
                }
                Ok(false) => {}
                Err(_) => {
                    mem.erred_records = mem.erred_records.saturating_add(1);
                    continue 'records;
                }
            }
        }
        // debug!(
        //     "worker {chrom_tid}:{start_pos}-{end_pos} processed {} records, \
        //      done. {} records failed",
        //     mem.ok_records, mem.erred_records
        // );
        Ok(mem)
    }
}

pub(crate) struct ProbsExtractor {
    rng: SmallRng,
    sample: bool,
    sample_frac: f64,
    motif_bases: [DnaBase; 4],
    edge_filter: Option<EdgeFilter>,
}

#[inline]
fn indexed_query_position_is_kept(
    edge_filter: Option<&EdgeFilter>,
    query_position: usize,
    read_length: usize,
    reverse: bool,
) -> bool {
    if query_position >= read_length {
        return false;
    }
    let forward_position =
        if reverse { read_length - 1 - query_position } else { query_position };
    edge_filter
        .map(|filter| {
            filter.keep_position(forward_position, read_length).unwrap_or(false)
        })
        .unwrap_or(true)
}

impl ProbsExtractor {
    fn get_aligned_mod_state_iterator<'a>(
        record: &'a bam::Record,
        chrom_coords: &'a ChromCoordinates,
        motif_bases: &'a [DnaBase; 4],
        records_with_base_mods: &'a mut [u32; 4],
        should_count_record: bool,
        start_pos: u32,
        end_pos: u32,
        edge_filter: Option<EdgeFilter>,
    ) -> anyhow::Result<impl Iterator<Item = MkResult<ModState>> + use<'a>>
    {
        let reverse = record.is_reverse();
        let reverse_offset = reverse as usize;
        let aligned_pairs_iter = AlignedBaseModsIterator::new(record)?;
        if should_count_record {
            let bases = aligned_pairs_iter.primary_bases_in_record();
            byte_to_bool_positions(bases, records_with_base_mods);
        }
        let read_length = record.seq_len();
        let mod_state_iter = aligned_pairs_iter
            .into_iter()
            .filter_ok(move |(_qpos, rpos, _mod_state)| {
                *rpos >= start_pos && *rpos < end_pos
            })
            .filter_ok(move |(qpos, _rpos, _mod_state)| {
                indexed_query_position_is_kept(
                    edge_filter.as_ref(),
                    *qpos as usize,
                    read_length,
                    reverse,
                )
            })
            .filter_map_ok(move |(qpos, rpos, mod_state)| {
                let rpos = rpos
                    .checked_sub(start_pos)
                    .expect("rpos should be larger than start pos");

                let mod_primary_base = mod_state.primary_base;
                match &chrom_coords.focus_positions {
                    FocusPositions2::MotifMask { mask, num_motifs } => {
                        let st = (rpos * 2 * (*num_motifs as u32))
                            + (*num_motifs * reverse as u8) as u32;
                        let st = st as usize;
                        let end = st + (*num_motifs as usize);
                        let bs = &mask[st..end];
                        debug_assert_eq!(
                            bs.len(),
                            *num_motifs as usize,
                            "bs:{bs:?} st:{st} end:{end} \
                             num_motifs:{num_motifs}"
                        );
                        #[cfg(debug_assertions)]
                        {
                            assert!(bs.count_ones() <= 1);
                        }
                        if bs[0] && mod_primary_base == motif_bases[0] {
                            Some((qpos, rpos, mod_state))
                        } else if *num_motifs >= 2u8
                            && bs[1]
                            && mod_primary_base == motif_bases[1]
                        {
                            Some((qpos, rpos, mod_state))
                        } else if *num_motifs >= 3u8
                            && bs[2]
                            && mod_primary_base == motif_bases[2]
                        {
                            Some((qpos, rpos, mod_state))
                        } else if *num_motifs >= 4u8
                            && bs[3]
                            && mod_primary_base == motif_bases[3]
                        {
                            Some((qpos, rpos, mod_state))
                        } else {
                            None
                        }
                    }
                    FocusPositions2::MaskedPositions { mask } => {
                        let keep =
                            mask[(rpos as usize * 2usize) + reverse_offset];
                        if keep {
                            Some((qpos, rpos, mod_state))
                        } else {
                            None
                        }
                    }
                    FocusPositions2::AllPositions => unreachable!(),
                }
            })
            .map(|res| res.map(|(_qpos, _rpos, mod_state)| mod_state));
        Ok(mod_state_iter)
    }
}

#[inline]
fn increment_mods_counts<R: rand::Rng>(
    mod_state: ModState,
    rng: &mut R,
    explicit_canonical_probs: &mut [u8; 4],
    base_hist: &mut [[u64; 256]; 4],
    base_totals: &mut [u64; 4],
    mods_hists: &mut Vec<ModHist>,
) {
    let qual_idx = mod_state.mod_qual as usize;
    if mod_state.modified {
        let mod_idx = mods_hists.iter_mut().find(|x| {
            x.mod_code == mod_state.mod_code
                && x.dna_base == mod_state.primary_base
        });

        if let Some(mod_hist) = mod_idx {
            match mod_hist.total.checked_add(1) {
                Some(val) => mod_hist.total = val,
                None => {
                    debug_once!(
                        "worker over-sampled mod probs for {}, performing \
                         sampling",
                        mod_state.mod_code
                    );
                    debug_assert_eq!(mod_hist.total, u64::MAX);
                    let target = rng.gen_range(0..u64::MAX);
                    let mut cum_sum = 0u64;
                    'decrement: for i in 0..256usize {
                        cum_sum = cum_sum.saturating_add(mod_hist.hist[i]);
                        if cum_sum >= target {
                            mod_hist.hist[i] =
                                mod_hist.hist[i].saturating_sub(1);
                            break 'decrement;
                        }
                    }
                }
            }
            let count = mod_hist.hist[qual_idx].saturating_add(1);
            mod_hist.hist[qual_idx] = count;
        } else {
            let mod_hist = ModHist::new_one(&mod_state);
            // mod_hist.check().expect("new mod hist should be valid");
            mods_hists.push(mod_hist);
        }
    } else {
        let i = mod_state.mod_qual as usize;
        let primary_base_idx = mod_state.primary_base as usize;
        let this_base_hist = &mut base_hist[primary_base_idx];
        match base_totals[primary_base_idx].checked_add(1) {
            Some(val) => base_totals[primary_base_idx] = val,
            None => {
                debug_once!(
                    "worker over-sampled base probs for {}, performing \
                     sampling",
                    mod_state.primary_base
                );
                debug_assert_eq!(base_totals[primary_base_idx], u64::MAX);
                let target = rng.gen_range(0..u64::MAX);
                let mut cum_sum = 0u64;
                'decrement: for i in 0..256usize {
                    cum_sum = cum_sum.saturating_add(this_base_hist[i]);
                    if cum_sum >= target {
                        this_base_hist[i] = this_base_hist[i].saturating_sub(1);
                        break 'decrement;
                    }
                }
            }
        }
        let count = base_hist[primary_base_idx][i].saturating_add(1);
        base_hist[primary_base_idx][i] = count;
        if !mod_state.modified && !mod_state.inferred {
            let p = explicit_canonical_probs[mod_state.primary_base as usize];
            explicit_canonical_probs[mod_state.primary_base as usize] =
                std::cmp::max(p, mod_state.mod_qual);
        }
    }
}

fn increment_base_counts<R: rand::Rng>(
    mod_state: ModState,
    explicit_canonical_probs: &mut [u8; 4],
    base_hist: &mut [[u64; 256]; 4],
    base_totals: &mut [u64; 4],
    rng: &mut R,
) {
    let i = mod_state.mod_qual as usize;
    let primary_base_idx = mod_state.primary_base as usize;
    let this_base_hist = &mut base_hist[primary_base_idx];
    match base_totals[primary_base_idx].checked_add(1) {
        Some(val) => base_totals[primary_base_idx] = val,
        None => {
            debug_once!(
                "worker over-sampled base probs for {}, performing sampling",
                mod_state.primary_base
            );
            debug_assert_eq!(base_totals[primary_base_idx], u64::MAX);
            let target = rng.gen_range(0..u64::MAX);
            let mut cum_sum = 0u64;
            'decrement: for i in 0..256usize {
                cum_sum = cum_sum.saturating_add(this_base_hist[i]);
                if cum_sum >= target {
                    this_base_hist[i] = this_base_hist[i].saturating_sub(1);
                    break 'decrement;
                }
            }
        }
    }
    let count = base_hist[primary_base_idx][i].saturating_add(1);
    base_hist[primary_base_idx][i] = count;
    if !mod_state.modified && !mod_state.inferred {
        let p = explicit_canonical_probs[mod_state.primary_base as usize];
        explicit_canonical_probs[mod_state.primary_base as usize] =
            std::cmp::max(p, mod_state.mod_qual);
    }
}

impl ExtractsMleProbs<BaseArgmaxProbs> for ProbsExtractor {
    fn new(
        seed: u64,
        sample: bool,
        sample_frac: f64,
        motif_bases: [DnaBase; 4],
        edge_filter: Option<&EdgeFilter>,
    ) -> Self {
        let rng = SmallRng::seed_from_u64(seed);
        Self {
            rng,
            sample,
            sample_frac,
            motif_bases,
            edge_filter: edge_filter.copied(),
        }
    }

    #[inline]
    fn process_record(
        &mut self,
        record: &bam::Record,
        chrom_coords: &ChromCoordinates,
        explicit_canonical_probs: &mut [u8; 4],
        base_hist: &mut [[u64; 256]; 4],
        base_totals: &mut [u64; 4],
        _mods_hists: &mut Vec<ModHist>,
        records_with_base_mods: &mut [u32; 4],
    ) -> anyhow::Result<bool> {
        if self.sample {
            if !self.rng.gen_bool(self.sample_frac) {
                return Ok(false);
            }
        }

        if !chrom_coords.final_interval {
            let aln_end = record.reference_end();
            if aln_end < 0i64 || (aln_end as u32) > chrom_coords.end_pos {
                return Ok(false);
            }
        }

        let mut modbase_iter = BaseModsAdapter::<16>::new(record)?;
        let bases = modbase_iter.primary_bases_in_record();
        byte_to_bool_positions(bases, records_with_base_mods);
        let read_length = record.seq_len();
        let reverse = record.is_reverse();
        loop {
            match modbase_iter.next_modified_position_no_thresh() {
                Ok(Some(mod_state)) => {
                    if indexed_query_position_is_kept(
                        self.edge_filter.as_ref(),
                        mod_state.mod_position,
                        read_length,
                        reverse,
                    ) {
                        increment_base_counts(
                            mod_state,
                            explicit_canonical_probs,
                            base_hist,
                            base_totals,
                            &mut self.rng,
                        );
                    }
                }
                Ok(None) => break,
                Err(e) => {
                    return Err(e.into());
                }
            }
        }
        Ok(true)
    }
}

impl ExtractsMleProbs<BaseAndModArgmaxProbs> for ProbsExtractor {
    fn new(
        seed: u64,
        sample: bool,
        sample_frac: f64,
        motif_bases: [DnaBase; 4],
        edge_filter: Option<&EdgeFilter>,
    ) -> Self {
        let rng = SmallRng::seed_from_u64(seed);
        Self {
            rng,
            sample,
            sample_frac,
            motif_bases,
            edge_filter: edge_filter.copied(),
        }
    }

    #[inline]
    fn process_record(
        &mut self,
        record: &bam::Record,
        chrom_coords: &ChromCoordinates,
        explicit_canonical_probs: &mut [u8; 4],
        base_hist: &mut [[u64; 256]; 4],
        base_totals: &mut [u64; 4],
        mods_hists: &mut Vec<ModHist>,
        records_with_base_mods: &mut [u32; 4],
    ) -> anyhow::Result<bool> {
        if self.sample {
            if !self.rng.gen_bool(self.sample_frac) {
                return Ok(false);
            }
        }
        if !chrom_coords.final_interval {
            let aln_end = record.reference_end();
            if aln_end < 0i64 || (aln_end as u32) > chrom_coords.end_pos {
                return Ok(false);
            }
        }

        let mut modbase_iter = BaseModsAdapter::<16>::new(record)?;
        let bases = modbase_iter.primary_bases_in_record();
        byte_to_bool_positions(bases, records_with_base_mods);
        let read_length = record.seq_len();
        let reverse = record.is_reverse();
        loop {
            match modbase_iter.next_modified_position_no_thresh() {
                Ok(Some(mod_state)) => {
                    if indexed_query_position_is_kept(
                        self.edge_filter.as_ref(),
                        mod_state.mod_position,
                        read_length,
                        reverse,
                    ) {
                        increment_mods_counts(
                            mod_state,
                            &mut self.rng,
                            explicit_canonical_probs,
                            base_hist,
                            base_totals,
                            mods_hists,
                        );
                    }
                }
                Ok(None) => break,
                Err(e) => {
                    return Err(e.into());
                }
            }
        }
        Ok(true)
    }
}

impl ExtractsMleProbs<AlignedBaseAndModArgmaxProbs> for ProbsExtractor {
    fn new(
        seed: u64,
        sample: bool,
        sample_frac: f64,
        motif_bases: [DnaBase; 4],
        edge_filter: Option<&EdgeFilter>,
    ) -> Self {
        let rng = SmallRng::seed_from_u64(seed);
        Self {
            rng,
            sample,
            sample_frac,
            motif_bases,
            edge_filter: edge_filter.copied(),
        }
    }

    #[inline]
    fn process_record(
        &mut self,
        record: &bam::Record,
        chrom_coords: &ChromCoordinates,
        explicit_canonical_probs: &mut [u8; 4],
        base_hist: &mut [[u64; 256]; 4],
        base_totals: &mut [u64; 4],
        mods_hists: &mut Vec<ModHist>,
        records_with_base_mods: &mut [u32; 4],
    ) -> anyhow::Result<bool> {
        if self.sample {
            if !self.rng.gen_bool(self.sample_frac) {
                return Ok(false);
            }
        }
        let start_pos = chrom_coords.start_pos;
        let end_pos = chrom_coords.end_pos;
        let should_count_record = if !chrom_coords.final_interval {
            let aln_end = record.reference_end();
            if aln_end < 0i64 || (aln_end as u32) > chrom_coords.end_pos {
                false
            } else {
                true
            }
        } else {
            // is the final interval
            true
        };
        let modstate_iter = Self::get_aligned_mod_state_iterator(
            record,
            chrom_coords,
            &self.motif_bases,
            records_with_base_mods,
            should_count_record,
            start_pos,
            end_pos,
            self.edge_filter,
        )?;

        for res in modstate_iter {
            match res {
                Ok(mod_state) => {
                    increment_mods_counts(
                        mod_state,
                        &mut self.rng,
                        explicit_canonical_probs,
                        base_hist,
                        base_totals,
                        mods_hists,
                    );
                }
                Err(e) => {
                    return Err(e.into());
                }
            }
        }
        Ok(should_count_record)
    }
}

impl ExtractsMleProbs<AlignedBaseArgmaxProbs> for ProbsExtractor {
    fn new(
        seed: u64,
        sample: bool,
        sample_frac: f64,
        motif_bases: [DnaBase; 4],
        edge_filter: Option<&EdgeFilter>,
    ) -> Self {
        let rng = SmallRng::seed_from_u64(seed);
        Self {
            rng,
            sample,
            sample_frac,
            motif_bases,
            edge_filter: edge_filter.copied(),
        }
    }

    #[inline]
    fn process_record(
        &mut self,
        record: &bam::Record,
        chrom_coords: &ChromCoordinates,
        explicit_canonical_probs: &mut [u8; 4],
        base_hist: &mut [[u64; 256]; 4],
        base_totals: &mut [u64; 4],
        _mods_hists: &mut Vec<ModHist>,
        records_with_base_mods: &mut [u32; 4],
    ) -> anyhow::Result<bool> {
        if self.sample {
            if !self.rng.gen_bool(self.sample_frac) {
                return Ok(false);
            }
        }
        let should_count_record = if !chrom_coords.final_interval {
            let aln_end = record.reference_end();
            if aln_end < 0i64 || (aln_end as u32) > chrom_coords.end_pos {
                false
            } else {
                true
            }
        } else {
            // is the final interval
            true
        };
        let start_pos = chrom_coords.start_pos;
        let end_pos = chrom_coords.end_pos;
        let modstate_iter = Self::get_aligned_mod_state_iterator(
            record,
            chrom_coords,
            &self.motif_bases,
            records_with_base_mods,
            should_count_record,
            start_pos,
            end_pos,
            self.edge_filter,
        )?;
        for res in modstate_iter {
            match res {
                Ok(mod_state) => {
                    increment_base_counts(
                        mod_state,
                        explicit_canonical_probs,
                        base_hist,
                        base_totals,
                        &mut self.rng,
                    );
                }
                Err(e) => {
                    return Err(e.into());
                }
            }
        }
        Ok(should_count_record)
    }
}

pub(crate) fn run_extract_probs_workers(
    workers: Vec<Box<dyn ExtractProbsWorker>>,
    multi_progress: MultiProgress,
    feeder: ChromCoordinatesFeeder,
) -> anyhow::Result<QualHist> {
    let tid_progress =
        multi_progress.add(get_master_progress_bar(feeder.total_length()));
    let (snd_work, rcv_work) = crossbeam::channel::unbounded();
    let (snd_results, rcv_results) = crossbeam::channel::unbounded();
    let (return_mem, rcv_mem) = crossbeam::channel::unbounded();
    for _ in 0..workers.len() {
        return_mem.send(QualHist::default()).expect("should stage memory 1");
        return_mem.send(QualHist::default()).expect("should stage memory 2");
    }

    let source_thread = std::thread::spawn({
        let results_tx = snd_results.clone();
        let mpb_handle = multi_progress.clone();
        let chrom_coord_feeder = feeder
            .into_iter()
            .inspect(move |r| match r {
                Ok(_) => {}
                Err(e) => {
                    mpb_handle
                        .suspend(|| error!("failed to fetch sequence, {e}"));
                }
            })
            .filter_ok(|cc| match &cc.focus_positions {
                FocusPositions2::MotifMask { mask, num_motifs: _ }
                | FocusPositions2::MaskedPositions { mask } => mask.any(),
                FocusPositions2::AllPositions => true,
            })
            .filter_map(|r| r.ok());
        move || {
            let get_mem = || -> Result<QualHist, ()> {
                match rcv_mem.recv() {
                    Ok(qh) => Ok(qh),
                    Err(_) => Err(()),
                }
            };

            for chrom_coords in chrom_coord_feeder {
                match get_mem() {
                    Ok(qh) => {
                        if snd_work.send((chrom_coords, qh)).is_err() {
                            break;
                        }
                    }
                    Err(_) => break,
                }
            }
            drop(snd_work);
            drop(results_tx);
        }
    });

    let mut worker_handles = Vec::with_capacity(workers.len());
    for mut worker in workers {
        let rcv_work = rcv_work.clone();
        let snd_results = snd_results.clone();
        worker_handles.push(std::thread::spawn(move || {
            while let Ok((chrom_coords, mem)) = rcv_work.recv() {
                let result = worker.process(chrom_coords, mem);
                if snd_results.send(result).is_err() {
                    break;
                }
            }
        }));
    }
    drop(snd_results);
    drop(rcv_work);

    let mut qual_hist = QualHist::default();
    for result in rcv_results {
        match result {
            Ok(mut qh) => {
                qual_hist.combine(&qh, &multi_progress)?;
                tid_progress.inc(qh.interval_length as u64);
                qh.clear();
                let _ = return_mem.send(qh);
            }
            Err(e) => {
                debug!("received {e} from rcv_results, stopping");
                break;
            }
        }
    }

    source_thread.join().expect("panic on source thread");
    for (i, worker_thread) in worker_handles.into_iter().enumerate() {
        worker_thread.join().expect(&format!("worker {i} paniced"));
    }

    Ok(qual_hist)
}

pub(crate) fn calc_per_base_thresholds_from_indexed_hts_file(
    bam_fp: &PathBuf,
    reference_fasta: Option<&PathBuf>,
    filter_percentile: f32,
    allow_non_primary: bool,
    preload_references: bool,
    mask: bool,
    sample_frac: Option<f64>,
    num_reads: usize,
    max_thresholds_per_base: Option<[f32; 4]>,
    stranded_position_filter: Option<StrandedPositionFilter<()>>,
    raw_motifs: Option<&Vec<String>>,
    cpg: bool,
    n_workers: usize,
    interval_size: u32,
    seed: Option<u64>,
    sampling_region: Option<&Region>,
    edge_filter: Option<&EdgeFilter>,
    io_threadpool: &rust_htslib::tpool::ThreadPool,
    multi_progress: MultiProgress,
) -> anyhow::Result<[f32; 4]> {
    let motifs = if let Some(raw_motif_parts) = raw_motifs {
        Some(RegexMotif::from_raw_parts(&raw_motif_parts, cpg)?)
    } else if cpg {
        Some(vec![RegexMotif::parse_string("CG", 0).unwrap()])
    } else {
        None
    };
    let mut qual_hist = get_base_mods_quals_from_indexed_hts_file(
        bam_fp,
        reference_fasta,
        false,
        allow_non_primary,
        preload_references,
        mask,
        sample_frac,
        Some(num_reads),
        stranded_position_filter,
        motifs,
        n_workers,
        interval_size,
        seed,
        sampling_region,
        edge_filter,
        io_threadpool,
        multi_progress.clone(),
    )?;
    if qual_hist.ok_records < 100 {
        multi_progress
            .suspend(|| info!("collecting probabilities from unmapped reads"));
        let mut records = bam::IndexedReader::from_path(bam_fp)?;
        records.fetch(FetchDefinition::Unmapped)?;
        let unmapped_qual_hist = QualHist::from_records(
            records.records(),
            None,
            Some(num_reads.saturating_sub(qual_hist.ok_records)),
            None,
            seed,
            edge_filter,
            false,
            allow_non_primary,
            false,
            &multi_progress,
        )?;
        qual_hist.combine(&unmapped_qual_hist, &multi_progress)?;
    }

    qual_hist.get_base_thresholds(
        filter_percentile,
        max_thresholds_per_base,
        &multi_progress,
    )
}

pub(crate) fn get_base_mods_quals_from_indexed_hts_file(
    bam_fp: &PathBuf,
    reference_fasta: Option<&PathBuf>,
    collect_mod_histograms: bool,
    allow_non_primary: bool,
    preload_references: bool,
    mask: bool,
    sample_frac: Option<f64>,
    num_reads: Option<usize>,
    stranded_position_filter: Option<StrandedPositionFilter<()>>,
    motifs: Option<Vec<RegexMotif>>,
    n_workers: usize,
    interval_size: u32,
    seed: Option<u64>,
    sampling_region: Option<&Region>,
    edge_filter: Option<&EdgeFilter>,
    io_threadpool: &rust_htslib::tpool::ThreadPool,
    multi_progress: MultiProgress,
) -> anyhow::Result<QualHist> {
    let bam_reader = bam::IndexedReader::from_path(bam_fp)?;
    let reference_records = get_targets(bam_reader.header(), sampling_region);
    let mut workers: Vec<Box<dyn ExtractProbsWorker>> = Vec::new();
    let (chrom_to_counts, rng_sample, sample_frac) =
        calculate_reads_per_contig(
            bam_reader.header(),
            sample_frac,
            num_reads,
            &multi_progress,
            interval_size,
            sampling_region,
        )?;
    let chrom_to_counts = chrom_to_counts.map(|x| Arc::new(x));
    let motif_lookup = get_motif_lookup_from_parts(
        motifs,
        reference_fasta,
        false,
        mask,
        preload_references,
    )?;
    let feeder = ChromCoordinatesFeeder::new(
        &reference_records,
        interval_size,
        motif_lookup,
        false,
        stranded_position_filter.clone(),
    )?;
    if let Some(motif_bases) = feeder.get_motif_bases() {
        for i in 0..n_workers {
            let worker: Box<dyn ExtractProbsWorker> = if collect_mod_histograms
            {
                Box::new(RegionMleProbs::<
                    AlignedBaseAndModArgmaxProbs,
                    ProbsExtractor,
                >::new(
                    bam_fp,
                    reference_fasta,
                    allow_non_primary,
                    motif_bases,
                    edge_filter,
                    io_threadpool,
                    seed.unwrap_or(i as u64),
                    rng_sample,
                    sample_frac,
                    chrom_to_counts.clone(),
                )?)
            } else {
                Box::new(RegionMleProbs::<
                    AlignedBaseArgmaxProbs,
                    ProbsExtractor,
                >::new(
                    bam_fp,
                    reference_fasta,
                    allow_non_primary,
                    motif_bases,
                    edge_filter,
                    io_threadpool,
                    seed.unwrap_or(i as u64),
                    rng_sample,
                    sample_frac,
                    chrom_to_counts.clone(),
                )?)
            };
            workers.push(worker);
        }
    } else if feeder.has_position_filter() {
        for i in 0..n_workers {
            let worker: Box<dyn ExtractProbsWorker> = if collect_mod_histograms
            {
                Box::new(RegionMleProbs::<
                    AlignedBaseAndModArgmaxProbs,
                    ProbsExtractor,
                >::new(
                    bam_fp,
                    reference_fasta,
                    allow_non_primary,
                    [DnaBase::A; 4],
                    edge_filter,
                    io_threadpool,
                    seed.unwrap_or(i as u64),
                    rng_sample,
                    sample_frac,
                    chrom_to_counts.clone(),
                )?)
            } else {
                Box::new(RegionMleProbs::<
                    AlignedBaseArgmaxProbs,
                    ProbsExtractor,
                >::new(
                    bam_fp,
                    reference_fasta,
                    allow_non_primary,
                    [DnaBase::A; 4],
                    edge_filter,
                    io_threadpool,
                    seed.unwrap_or(i as u64),
                    rng_sample,
                    sample_frac,
                    chrom_to_counts.clone(),
                )?)
            };
            workers.push(worker);
        }
    } else {
        for i in 0..n_workers {
            let worker: Box<dyn ExtractProbsWorker> = if collect_mod_histograms
            {
                Box::new(RegionMleProbs::<BaseAndModArgmaxProbs, ProbsExtractor>::new(
                                    bam_fp,
                                    reference_fasta,
                                    allow_non_primary,
                                    [DnaBase::A; 4],
                                    edge_filter,
                                    io_threadpool,
                                    seed.unwrap_or(i as u64),
                                    rng_sample,
                                    sample_frac,
                                    chrom_to_counts.clone(),
                                )?)
            } else {
                Box::new(
                    RegionMleProbs::<BaseArgmaxProbs, ProbsExtractor>::new(
                        bam_fp,
                        reference_fasta,
                        allow_non_primary,
                        [DnaBase::A; 4],
                        edge_filter,
                        io_threadpool,
                        seed.unwrap_or(i as u64),
                        rng_sample,
                        sample_frac,
                        chrom_to_counts.clone(),
                    )?,
                )
            };
            workers.push(worker);
        }
    }
    run_extract_probs_workers(workers, multi_progress.clone(), feeder.clone())
}

pub(crate) fn calculate_reads_per_contig(
    bam_header: &HeaderView,
    sampling_frac: Option<f64>,
    num_reads: Option<usize>,
    multi_progress: &MultiProgress,
    interval_size: u32,
    sampling_region: Option<&Region>,
) -> anyhow::Result<(Option<FxHashMap<u32, usize>>, bool, f64)> {
    let (chrom_to_counts, rng_sample, sample_frac) =
        if let Some(user_sampling_frac) = sampling_frac {
            multi_progress.suspend(|| {
                info!("sampling {}% of reads", user_sampling_frac * 100f64);
            });
            (None, true, user_sampling_frac)
        } else if let Some(num_reads) = num_reads {
            multi_progress.suspend(|| {
                info!("attempting to sample at least {num_reads} reads");
            });
            let counts = SampleModBaseProbs::calc_counts_per_chrom(
                interval_size,
                bam_header,
                num_reads,
                sampling_region,
            )?;
            let mut sample_schedule_table = get_human_readable_table();
            sample_schedule_table.set_titles(row![
                "contig",
                "num_reads_per_interval",
                "contig_length",
                "intervals_over_contig",
                "total_reads_per_contig"
            ]);
            for (tid, count) in
                counts.iter().sorted_by(|(a, _), (b, _)| a.cmp(b))
            {
                let contig_name = bam_header
                    .tid2name(*tid)
                    .into_iter()
                    .map(|b| *b as char)
                    .collect::<String>();
                let contig_length = bam_header
                    .target_len(*tid)
                    .ok_or(anyhow!("header missing length for {tid}"))?;
                let intervals_per_contig =
                    std::cmp::max(contig_length / interval_size as u64, 1u64);
                let total_reads_for_contig =
                    count.saturating_mul(intervals_per_contig as usize);
                sample_schedule_table.add_row(row![
                    contig_name,
                    count,
                    contig_length,
                    intervals_per_contig,
                    total_reads_for_contig
                ]);
            }
            debug!(
                "sample schedule for {num_reads} \
                 reads\n{sample_schedule_table}",
            );
            (Some(counts), false, 0f64)
        } else {
            (None, false, 0f64)
        };
    Ok((chrom_to_counts, rng_sample, sample_frac))
}

#[inline]
fn byte_to_bool_positions<const SIZE: usize>(b: u8, agg: &mut [u32; SIZE]) {
    for (i, _) in b.view_bits::<Lsb0>().iter().enumerate().filter(|(_, x)| **x)
    {
        let count = agg[i];
        agg[i] = count.saturating_add(1u32);
    }
}

#[cfg(test)]
mod tests {
    use bitvec::bitvec;
    use rust_htslib::bam::record::{Aux, Cigar, CigarString};

    use super::*;

    fn make_five_base_record(reverse: bool) -> bam::Record {
        let mut record = bam::Record::new();
        let cigar = CigarString(vec![Cigar::Match(5)]);
        let seq = if reverse { b"TTTTT" } else { b"AAAAA" };
        record.set(b"read", Some(&cigar), seq, &[255; 5]);
        record.set_tid(0);
        record.set_pos(0);
        record.push_aux(b"MM", Aux::String("A+a?,0,0,0,0,0;")).unwrap();
        let ml = vec![201u8, 202, 203, 204, 205];
        record.push_aux(b"ML", Aux::ArrayU8((&ml).into())).unwrap();
        if reverse {
            record.set_reverse();
        }
        record
    }

    fn retained_qualities<T>(
        record: &bam::Record,
        edge_filter: Option<&EdgeFilter>,
    ) -> Vec<u8>
    where
        ProbsExtractor: ExtractsMleProbs<T>,
    {
        let mut extractor = <ProbsExtractor as ExtractsMleProbs<T>>::new(
            42,
            false,
            0.0,
            [DnaBase::A; 4],
            edge_filter,
        );
        let chrom_coords = ChromCoordinates::new(
            0,
            0,
            5,
            FocusPositions2::MaskedPositions { mask: bitvec![1; 10] },
            true,
        );
        let mut hist = QualHist::default();
        let processed =
            <ProbsExtractor as ExtractsMleProbs<T>>::process_record(
                &mut extractor,
                record,
                &chrom_coords,
                &mut hist.explicit_canonical_probs,
                &mut hist.hist,
                &mut hist.base_totals,
                &mut hist.mods_hists,
                &mut hist.num_records_with_base_mods,
            )
            .unwrap();
        assert!(processed);

        let mut qualities = hist
            .hist
            .iter()
            .chain(hist.mods_hists.iter().map(|mod_hist| &mod_hist.hist))
            .flat_map(|counts| {
                counts.iter().enumerate().flat_map(|(qual, count)| {
                    std::iter::repeat(qual as u8).take(*count as usize)
                })
            })
            .collect::<Vec<_>>();
        qualities.sort_unstable();
        qualities
    }

    macro_rules! handler_results {
        ($marker:ty, $forward:expr, $reverse:expr, $filter:expr) => {
            (
                retained_qualities::<$marker>($forward, None),
                retained_qualities::<$marker>($forward, Some($filter)),
                retained_qualities::<$marker>($reverse, None),
                retained_qualities::<$marker>($reverse, Some($filter)),
            )
        };
    }

    #[test]
    fn indexed_handlers_apply_asymmetric_edge_filter() {
        let forward = make_five_base_record(false);
        let reverse = make_five_base_record(true);
        let edge_filter = EdgeFilter::new(1, 2, false);
        let actual = [
            handler_results!(
                AlignedBaseArgmaxProbs,
                &forward,
                &reverse,
                &edge_filter
            ),
            handler_results!(
                AlignedBaseAndModArgmaxProbs,
                &forward,
                &reverse,
                &edge_filter
            ),
            handler_results!(BaseArgmaxProbs, &forward, &reverse, &edge_filter),
            handler_results!(
                BaseAndModArgmaxProbs,
                &forward,
                &reverse,
                &edge_filter
            ),
        ];
        let all_qualities = vec![201, 202, 203, 204, 205];
        let retained_qualities = vec![202, 203];

        for (forward_all, forward_filtered, reverse_all, reverse_filtered) in
            actual
        {
            assert_eq!(forward_all, all_qualities);
            assert_eq!(forward_filtered, retained_qualities);
            assert_eq!(reverse_all, all_qualities);
            assert_eq!(reverse_filtered, retained_qualities);
        }
    }

    #[test]
    fn indexed_handlers_treat_overlong_edge_filter_as_empty() {
        let record = make_five_base_record(false);
        let edge_filter = EdgeFilter::new(1, 6, false);
        let actual = [
            retained_qualities::<BaseArgmaxProbs>(&record, Some(&edge_filter)),
            retained_qualities::<BaseAndModArgmaxProbs>(
                &record,
                Some(&edge_filter),
            ),
            retained_qualities::<AlignedBaseArgmaxProbs>(
                &record,
                Some(&edge_filter),
            ),
            retained_qualities::<AlignedBaseAndModArgmaxProbs>(
                &record,
                Some(&edge_filter),
            ),
        ];

        for qualities in actual {
            assert!(qualities.is_empty());
        }
    }

    #[test]
    fn indexed_query_filter_preserves_inversion() {
        let edge_filter = EdgeFilter::new(1, 2, true);
        let retained = |reverse| {
            (0..5)
                .filter(|qpos| {
                    indexed_query_position_is_kept(
                        Some(&edge_filter),
                        *qpos,
                        5,
                        reverse,
                    )
                })
                .collect::<Vec<_>>()
        };

        assert_eq!(retained(false), vec![0, 3, 4]);
        assert_eq!(retained(true), vec![0, 1, 4]);
    }
}
