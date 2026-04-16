use std::{
    collections::HashMap, marker::PhantomData, path::PathBuf, sync::Arc, usize,
};

use anyhow::{anyhow, bail};
use derive_new::new;
use indicatif::MultiProgress;
use itertools::Itertools;
use log::{debug, info};
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
    mod_bam::{validate_mn_tag_on_record, EdgeFilter},
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
        self.base_totals.iter_mut().for_each(|x| *x = 0u64);
        self.mods_hists.iter_mut().for_each(|h| h.clear());
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
        if show_table {
            multi_progress.suspend(|| {
                info!(
                    "Total number of probabilities collected per primary \
                     sequence base:\n{total_counts_table}"
                );
            });
        }

        Self::percentiles(&agg.hist, &agg.base_totals, desired_percentiles)
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
                    debug!(
                        "sampled {} records, done. {} records failed",
                        mem.ok_records, mem.erred_records
                    );
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
                                    &mut mem.hist,
                                    &mut mem.base_totals,
                                    &mut mem.mods_hists,
                                );
                            } else {
                                increment_base_counts(
                                    mod_state,
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
            } else {
                let Ok(mut modbase_iter) = BaseModsAdapter::<16>::new(&record)
                else {
                    mem.erred_records = mem.erred_records.saturating_add(1);
                    continue 'records;
                };
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
                                        &mut mem.hist,
                                        &mut mem.base_totals,
                                        &mut mem.mods_hists,
                                    );
                                } else {
                                    increment_base_counts(
                                        mod_state,
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
            multi_progress.suspend(|| {
                info!(
                    "setting threshold {t} (qual: {q}) for base {base}, \
                     percentile {}, {} total probabilites",
                    filter_percentile, self.base_totals[base as usize]
                );
            });
            base_thresholds[base as usize] = t;
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
            base_totals: Default::default(),
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
        base_hist: &mut [[u64; 256]; 4],
        base_totals: &mut [u64; 4],
        mods_hists: &mut Vec<ModHist>,
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
        debug!("worker starting on {chrom_tid}:{start_pos}-{end_pos}");

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
                    debug!(
                        "worker {chrom_tid}:{start_pos}-{end_pos} sampled {} \
                         records, done. {} records failed",
                        mem.ok_records, mem.erred_records
                    );
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
                debug!("non primary record encountered");
                if validate_mn_tag_on_record(&record).is_err() {
                    mem.erred_records = mem.erred_records.saturating_add(1);
                    continue 'records;
                }
            }
            match self.hist.process_record(
                &record,
                &item,
                &mut mem.hist,
                &mut mem.base_totals,
                &mut mem.mods_hists,
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
        debug!(
            "worker {chrom_tid}:{start_pos}-{end_pos} processed {} records, \
             done. {} records failed",
            mem.ok_records, mem.erred_records
        );
        Ok(mem)
    }
}

pub(crate) struct ProbsExtractor {
    rng: SmallRng,
    sample: bool,
    sample_frac: f64,
    motif_bases: [DnaBase; 4],
    edge_filter_start: usize,
    edge_filter_end: usize,
}

impl ProbsExtractor {
    fn get_aligned_mod_state_iterator<'a>(
        record: &'a bam::Record,
        chrom_coords: &'a ChromCoordinates,
        motif_bases: &'a [DnaBase; 4],
        start_pos: u32,
        end_pos: u32,
        edge_filter_start: usize,
        edge_filter_end: usize,
    ) -> anyhow::Result<impl Iterator<Item = MkResult<ModState>> + use<'a>>
    {
        let reverse = record.is_reverse();
        let reverse_offset = reverse as usize;
        let aligned_pairs_iter = AlignedBaseModsIterator::new(record)?;
        let read_length = record.seq_len();
        let mod_state_iter = aligned_pairs_iter
            .into_iter()
            .filter_ok(move |(_qpos, rpos, _mod_state)| {
                *rpos >= start_pos && *rpos < end_pos
            })
            .filter_ok(move |(qpos, _rpos, _mod_state)| {
                (*qpos as usize) >= edge_filter_start
                    && (*qpos as usize) < (read_length - edge_filter_end)
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
    }
}

fn increment_base_counts<R: rand::Rng>(
    mod_state: ModState,
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
}

impl ExtractsMleProbs<BaseArgmaxProbs> for ProbsExtractor {
    fn new(
        seed: u64,
        sample: bool,
        sample_frac: f64,
        motif_bases: [DnaBase; 4],
        _edge_filter: Option<&EdgeFilter>,
    ) -> Self {
        let rng = SmallRng::seed_from_u64(seed);
        Self {
            rng,
            sample,
            sample_frac,
            motif_bases,
            edge_filter_start: 0usize,
            edge_filter_end: 0usize,
        }
    }

    #[inline]
    fn process_record(
        &mut self,
        record: &bam::Record,
        chrom_corrds: &ChromCoordinates,
        base_hist: &mut [[u64; 256]; 4],
        base_totals: &mut [u64; 4],
        _mods_hists: &mut Vec<ModHist>,
    ) -> anyhow::Result<bool> {
        if self.sample {
            if !self.rng.gen_bool(self.sample_frac) {
                return Ok(false);
            }
        }

        if chrom_corrds.final_interval {
            let aln_end = record.reference_end();
            if aln_end < 0i64 {
                return Ok(false);
            } else if aln_end as u32 >= chrom_corrds.end_pos {
                return Ok(false);
            }
        }

        let mut modbase_iter = BaseModsAdapter::<16>::new(record)?;
        loop {
            match modbase_iter.next_modified_position_no_thresh() {
                Ok(Some(mod_state)) => {
                    increment_base_counts(
                        mod_state,
                        base_hist,
                        base_totals,
                        &mut self.rng,
                    );
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
        _edge_filter: Option<&EdgeFilter>,
    ) -> Self {
        let rng = SmallRng::seed_from_u64(seed);
        Self {
            rng,
            sample,
            sample_frac,
            motif_bases,
            edge_filter_start: 0usize,
            edge_filter_end: 0usize,
        }
    }

    #[inline]
    fn process_record(
        &mut self,
        record: &bam::Record,
        chrom_corrds: &ChromCoordinates,
        base_hist: &mut [[u64; 256]; 4],
        base_totals: &mut [u64; 4],
        mods_hists: &mut Vec<ModHist>,
    ) -> anyhow::Result<bool> {
        if self.sample {
            if !self.rng.gen_bool(self.sample_frac) {
                return Ok(false);
            }
        }

        if chrom_corrds.final_interval {
            let aln_end = record.reference_end();
            if aln_end < 0i64 {
                return Ok(false);
            } else if aln_end as u32 >= chrom_corrds.end_pos {
                return Ok(false);
            }
        }

        let mut modbase_iter = BaseModsAdapter::<16>::new(record)?;
        loop {
            match modbase_iter.next_modified_position_no_thresh() {
                Ok(Some(mod_state)) => {
                    increment_mods_counts(
                        mod_state,
                        &mut self.rng,
                        base_hist,
                        base_totals,
                        mods_hists,
                    );
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
        let edge_filter_start =
            edge_filter.map(|ef| ef.edge_filter_start).unwrap_or(0usize);
        let edge_filter_end =
            edge_filter.map(|ef| ef.edge_filter_end).unwrap_or(0usize);
        Self {
            rng,
            sample,
            sample_frac,
            motif_bases,
            edge_filter_start,
            edge_filter_end,
        }
    }

    #[inline]
    fn process_record(
        &mut self,
        record: &bam::Record,
        chrom_coords: &ChromCoordinates,
        base_hist: &mut [[u64; 256]; 4],
        base_totals: &mut [u64; 4],
        mods_hists: &mut Vec<ModHist>,
    ) -> anyhow::Result<bool> {
        if self.sample {
            if !self.rng.gen_bool(self.sample_frac) {
                return Ok(false);
            }
        }
        let start_pos = chrom_coords.start_pos;
        let end_pos = chrom_coords.end_pos;
        let modstate_iter = Self::get_aligned_mod_state_iterator(
            record,
            chrom_coords,
            &self.motif_bases,
            start_pos,
            end_pos,
            self.edge_filter_start,
            self.edge_filter_end,
        )?;
        for res in modstate_iter {
            match res {
                Ok(mod_state) => {
                    increment_mods_counts(
                        mod_state,
                        &mut self.rng,
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
        Ok(true)
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
        let edge_filter_start =
            edge_filter.map(|ef| ef.edge_filter_start).unwrap_or(0usize);
        let edge_filter_end =
            edge_filter.map(|ef| ef.edge_filter_end).unwrap_or(0usize);
        Self {
            rng,
            sample,
            sample_frac,
            motif_bases,
            edge_filter_start,
            edge_filter_end,
        }
    }

    #[inline]
    fn process_record(
        &mut self,
        record: &bam::Record,
        chrom_coords: &ChromCoordinates,
        base_hist: &mut [[u64; 256]; 4],
        base_totals: &mut [u64; 4],
        _mods_hists: &mut Vec<ModHist>,
    ) -> anyhow::Result<bool> {
        if self.sample {
            if !self.rng.gen_bool(self.sample_frac) {
                return Ok(false);
            }
        }
        let start_pos = chrom_coords.start_pos;
        let end_pos = chrom_coords.end_pos;
        let modstate_iter = Self::get_aligned_mod_state_iterator(
            record,
            chrom_coords,
            &self.motif_bases,
            start_pos,
            end_pos,
            self.edge_filter_start,
            self.edge_filter_end,
        )?;
        for res in modstate_iter {
            match res {
                Ok(mod_state) => {
                    increment_base_counts(
                        mod_state,
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
        Ok(true)
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
            .filter_ok(|cc| {
                let keep = match &cc.focus_positions {
                    FocusPositions2::MotifMask { mask, num_motifs: _ }
                    | FocusPositions2::MaskedPositions { mask } => mask.any(),
                    FocusPositions2::AllPositions => true,
                };
                if !keep {
                    debug!(
                        "discarding {}:{}-{}, zero positions",
                        cc.chrom_tid, cc.start_pos, cc.end_pos
                    );
                }
                keep
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

pub(crate) fn calc_per_base_thresholds_from_stream(
    bam_fp: &PathBuf,
    num_reads: usize,
    allow_non_primary: bool,
    stranded_position_filter: Option<StrandedPositionFilter<()>>,
    edge_filter: Option<&EdgeFilter>,
    filter_percentile: f32,
    io_threads: usize,
    multi_progress: &MultiProgress,
) -> anyhow::Result<[f32; 4]> {
    let mut records = bam::Reader::from_path(bam_fp)?;
    records.set_threads(io_threads)?;
    QualHist::from_records(
        records.records(),
        stranded_position_filter,
        Some(num_reads),
        None,
        None,
        edge_filter,
        false,
        allow_non_primary,
        multi_progress,
    )?
    .get_base_thresholds(filter_percentile, multi_progress)
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
    let motifs = if let Some(raw_motif_parts) = raw_motifs {
        Some(RegexMotif::from_raw_parts(&raw_motif_parts, cpg)?)
    } else if cpg {
        Some(vec![RegexMotif::parse_string("CG", 0).unwrap()])
    } else {
        None
    };
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
        multi_progress.suspend(|| {
            info!(
                "calculating threshold value with probabilities aligned to \
                 reference bases {}",
                motif_bases.iter().unique().join(",")
            )
        });
        for i in 0..n_workers {
            let worker =
                RegionMleProbs::<AlignedBaseArgmaxProbs, ProbsExtractor>::new(
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
                )?;
            workers.push(Box::new(worker));
        }
    } else if feeder.has_position_filter() {
        for i in 0..n_workers {
            let worker =
                RegionMleProbs::<AlignedBaseArgmaxProbs, ProbsExtractor>::new(
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
                )?;
            workers.push(Box::new(worker));
        }
    } else {
        for i in 0..n_workers {
            let worker =
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
                )?;
            workers.push(Box::new(worker));
        }
    }
    let mut qual_hist = run_extract_probs_workers(
        workers,
        multi_progress.clone(),
        feeder.clone(),
    )?;
    if qual_hist.ok_records < 100 {
        multi_progress
            .suspend(|| info!("collecting probabilities from unmapped reads"));
        let mut records = bam::IndexedReader::from_path(bam_fp)?;
        records.fetch(FetchDefinition::Unmapped)?;
        let unmapped_qual_hist = QualHist::from_records(
            records.records(),
            stranded_position_filter,
            Some(num_reads.saturating_sub(qual_hist.ok_records)),
            None,
            seed,
            edge_filter,
            false,
            allow_non_primary,
            &multi_progress,
        )?;
        qual_hist.combine(&unmapped_qual_hist, &multi_progress)?;
    }

    qual_hist.get_base_thresholds(filter_percentile, &multi_progress)
}

pub(crate) fn calculate_reads_per_contig(
    bam_header: &HeaderView,
    sampling_frac: Option<f64>,
    num_reads: usize,
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
        } else {
            multi_progress.suspend(|| {
                info!("attempting to sample at least {num_reads} reads");
            });
            let counts = SampleModBaseProbs::calc_counts_per_chrom(
                interval_size,
                bam_header,
                num_reads,
                sampling_region,
            )?;
            (Some(counts), false, 0f64)
        };
    if let Some(chrom_to_counts) = chrom_to_counts.as_ref() {
        let mut sample_schedule_table = get_human_readable_table();
        sample_schedule_table.set_titles(row![
            "contig",
            "num_reads_per_interval",
            "contig_length",
            "intervals_over_contig",
            "total_reads_per_contig"
        ]);
        for (tid, count) in
            chrom_to_counts.iter().sorted_by(|(a, _), (b, _)| a.cmp(b))
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
            "sample schedule for {num_reads} reads\n{sample_schedule_table}",
        );
        debug!("chrom_to_counts={chrom_to_counts:?}",);
    }
    Ok((chrom_to_counts, rng_sample, sample_frac))
}
