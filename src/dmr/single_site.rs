use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap, VecDeque};
use std::fs::File;
use std::io::Write;
use std::ops::Range;
use std::path::PathBuf;
use std::sync::Arc;

use anyhow::{bail, Context};
use derive_new::new;
use indicatif::{MultiProgress, ProgressBar};
use itertools::Itertools;
use log::{debug, error, info};
use rayon::prelude::*;

use crate::dmr::beta_diff::{BetaParams, PMapEstimator};
use crate::dmr::llr_model::{llk_ratio, AggregatedCounts};
use crate::dmr::tabix::{
    MultiSampleIndex, SampleToBedMethyLines, SingleSiteSampleIndex,
};
use crate::dmr::util::DmrBatchOfPositions;
use crate::genome_positions::GenomePositions;
use crate::hmm::{HmmModel, States};
use crate::mod_base_code::ModCodeRepr;
use crate::monoid::BorrowingMoniod;
use crate::thresholds::percentile_linear_interp;
use crate::util::{get_subroutine_progress_bar, get_ticker, Region};
use crate::writers::TsvWriter;

pub(super) struct SingleSiteDmrAnalysis {
    sample_index: Arc<SingleSiteSampleIndex>,
    genome_positions: Arc<GenomePositions>,
    pmap_estimator: Arc<PMapEstimator>,
    batch_size: usize,
    interval_size: u64,
    header: bool,
    segmentation_fp: Option<PathBuf>,
}

impl SingleSiteDmrAnalysis {
    pub(super) fn new(
        sample_index: MultiSampleIndex,
        genome_positions: GenomePositions,
        cap_coverages: bool,
        num_a: usize,
        num_b: usize,
        batch_size: usize,
        interval_size: u64,
        prior: Option<&Vec<f64>>,
        max_coverages: Option<&Vec<usize>>,
        rope: f64,
        sample_n: usize,
        header: bool,
        segmentation_fp: Option<&PathBuf>,
        progress: &MultiProgress,
    ) -> anyhow::Result<Self> {
        let sample_index =
            SingleSiteSampleIndex::new(sample_index, num_a, num_b)
                .map(|x| Arc::new(x))?;
        let genome_positions = Arc::new(genome_positions);
        if cap_coverages {
            info!("capping coverages when combining samples");
        }
        let prior = if let Some(raw_prior_params) = prior {
            if raw_prior_params[0] + raw_prior_params[1] < 1.0 {
                bail!("alpha + beta must be > 1.0 for numerical stability")
            }
            let prior =
                BetaParams::new(raw_prior_params[0], raw_prior_params[1])?;
            info!("using user-specified prior values {prior}");
            prior
        } else {
            let prior = BetaParams::new(0.55, 0.55).unwrap();
            info!("using default prior, {prior}");
            prior
        };

        let max_coverages = if let Some(raw_max_coverages) = max_coverages {
            let max_covs = [raw_max_coverages[0], raw_max_coverages[1]];
            info!(
                "using specified max coverage values control {} and \
                 experiment {}",
                max_covs[0], max_covs[1]
            );
            max_covs
        } else {
            calculate_max_coverages(
                sample_index.clone(),
                genome_positions.clone(),
                batch_size,
                interval_size,
                sample_n,
                progress,
            )?
        };
        if sample_index.min_valid_coverage() > 0 {
            info!(
                "min valid coverage set to {}",
                sample_index.min_valid_coverage()
            );
        }

        let pmap_estimator = Arc::new(PMapEstimator::new(
            max_coverages,
            num_a,
            num_b,
            prior,
            rope,
            cap_coverages,
        ));

        Ok(Self {
            sample_index,
            genome_positions,
            pmap_estimator,
            batch_size,
            interval_size,
            header,
            segmentation_fp: segmentation_fp.cloned(),
        })
    }

    pub(super) fn run(
        &self,
        multi_progress_bar: MultiProgress,
        pool: rayon::ThreadPool,
        max_gap_size: u64,
        dmr_prior: f64,
        diff_stay: f64,
        significance_factor: f64,
        decay_distance: u32,
        linear_transitions: bool,
        mut writer: Box<dyn Write>,
    ) -> anyhow::Result<()> {
        let matched_samples = self.sample_index.matched_replicate_samples();
        let multiple_samples = self.sample_index.multiple_samples();
        if matched_samples {
            info!("running with replicates and matched samples");
        } else if multiple_samples {
            info!("running with replicates, but not matched samples");
        }

        if self.header {
            writer.write(
                SingleSiteDmrScore::header(multiple_samples, matched_samples)
                    .as_bytes(),
            )?;
        }

        let mut segmenter: Box<dyn DmrSegmenter> =
            if let Some(segmentation_fp) = &self.segmentation_fp {
                Box::new(HmmDmrSegmenter::new(
                    segmentation_fp,
                    max_gap_size,
                    dmr_prior,
                    diff_stay,
                    0.3f64,
                    -0.1f64,
                    significance_factor,
                    linear_transitions,
                    decay_distance,
                    &multi_progress_bar,
                )?)
            } else {
                Box::new(DummySegmenter::new())
            };

        let (scores_snd, scores_rcv) = crossbeam::channel::bounded(1000);
        // let (segment_snd, segment_rcv) = crossbeam::channel::bounded(1000);
        let processed_batches = multi_progress_bar.add(get_ticker());
        let failure_counter = multi_progress_bar.add(get_ticker());
        let success_counter = multi_progress_bar.add(get_ticker());

        processed_batches.set_message("batches processed");
        failure_counter.set_message("sites failed");
        success_counter.set_message("sites processed successfully");

        let batch_iter = SingleSiteBatches::new(
            self.sample_index.clone(),
            self.genome_positions.clone(),
            self.batch_size,
            self.interval_size,
        )?;

        let sample_index = self.sample_index.clone();
        let pmap_estimator = self.pmap_estimator.clone();
        pool.spawn(move || {
            for super_batch in batch_iter.filter_map(|r| match r {
                Ok(super_batch) => Some(super_batch),
                Err(e) => {
                    debug!("batch failed, {e}");
                    None
                }
            }) {
                let mut results = Vec::new();
                let (super_batch_results, _) = rayon::join(
                    || {
                        super_batch
                            .into_par_iter()
                            .filter_map(|batch_of_positions| {
                                match process_batch_of_positions(
                                    batch_of_positions,
                                    sample_index.clone(),
                                    pmap_estimator.clone(),
                                ) {
                                    Ok(chrom_to_scores) => {
                                        Some(chrom_to_scores)
                                    }
                                    Err(_e) => {
                                        // error!("batch failed, {e}");
                                        None
                                    }
                                }
                            })
                            .collect::<Vec<Vec<ChromToSingleScores>>>()
                    },
                    || {
                        results.into_iter().for_each(
                            |chrom_to_scores: Vec<ChromToSingleScores>| {
                                match scores_snd.send(chrom_to_scores) {
                                    Ok(_) => processed_batches.inc(1),
                                    Err(e) => {
                                        error!(
                                            "failed to send on channel, {e}"
                                        );
                                    }
                                }
                            },
                        )
                    },
                );
                results = super_batch_results;
                results.into_iter().for_each(
                    |chrom_to_scores| match scores_snd.send(chrom_to_scores) {
                        Ok(_) => processed_batches.inc(1),
                        Err(e) => {
                            error!("failed to send on channel, {e}");
                        }
                    },
                );
            }
        });

        let mut success_count = 0usize;
        for batch_result in scores_rcv {
            if let Err(e) = segmenter.add(&batch_result) {
                debug!("segmentation error, {e}");
            }
            for (chrom, results) in batch_result {
                for result in results {
                    match result {
                        Ok(scores) => {
                            writer.write(
                                scores
                                    .to_row(
                                        multiple_samples,
                                        matched_samples,
                                        &chrom,
                                    )
                                    .as_bytes(),
                            )?;
                            success_counter.inc(1);
                            success_count += 1;
                        }
                        Err(e) => {
                            debug!("score error, {e}");
                            failure_counter.inc(1);
                        }
                    };
                }
            }
        }

        if let Err(e) = segmenter.run_current_chunk() {
            debug!("segmentation error, {e}")
        }
        success_counter.finish_and_clear();
        failure_counter.finish_and_clear();
        segmenter.clean_up()?;

        info!(
            "finished, processed {} sites successfully, {} failed",
            success_count,
            failure_counter.position(),
        );
        Ok(())
    }
}

struct SingleSiteBatches {
    /// map of contig ID to its size
    interval_queue: VecDeque<(String, Range<u64>)>,
    sample_index: Arc<SingleSiteSampleIndex>,
    genome_positions: Arc<GenomePositions>,
    batch_size: usize,
    interval_size: u64,
    curr_contig: String,
    curr_contig_end: u64,
    curr_pos: u64,
    done: bool,
}

impl SingleSiteBatches {
    fn new(
        sample_index: Arc<SingleSiteSampleIndex>,
        genome_positions: Arc<GenomePositions>,
        batch_size: usize,
        interval_size: u64,
    ) -> anyhow::Result<Self> {
        let mut interval_queue = genome_positions
            .contig_sizes()
            .filter(|(name, _)| sample_index.has_contig(name))
            .map(|(name, length)| (name.to_owned(), (0u64..(length as u64))))
            .sorted_by(|(a, _), (b, _)| a.cmp(b))
            .collect::<VecDeque<(String, Range<u64>)>>();

        if let Some((curr_contig, curr_contig_range)) =
            interval_queue.pop_front()
        {
            let curr_pos = curr_contig_range.start;
            let curr_contig_size = curr_contig_range.end;
            Ok(Self {
                interval_queue,
                sample_index,
                genome_positions,
                batch_size,
                interval_size,
                curr_contig,
                curr_contig_end: curr_contig_size,
                curr_pos,
                done: false,
            })
        } else {
            bail!("need at least 1 contig")
        }
    }

    fn update(&mut self) {
        if let Some((contig_name, range)) = self.interval_queue.pop_front() {
            self.curr_contig = contig_name;
            self.curr_pos = range.start;
            self.curr_contig_end = range.end;
        } else {
            self.done = true
        }
    }

    fn get_next_batch(
        &mut self,
    ) -> anyhow::Result<Option<Vec<DmrBatchOfPositions>>> {
        let mut batches = Vec::<DmrBatchOfPositions>::new();
        let mut current_batch = DmrBatchOfPositions::empty();
        let mut current_batch_length = 0u64;

        loop {
            if batches.len() >= self.batch_size || self.done {
                break;
            }

            let end = std::cmp::min(
                self.curr_pos.saturating_add(self.interval_size),
                self.curr_contig_end,
            );
            let interval = self.curr_pos..end;
            if let Some(positions) = self
                .genome_positions
                .get_positions(&self.curr_contig, &interval)
            {
                let interval_length = positions.len() as u64;
                let control_chunks = self
                    .sample_index
                    .query_control_chunks(&self.curr_contig, &interval)
                    .with_context(|| {
                        format!(
                            "failed to get control chunks at {}:{}-{}",
                            &self.curr_contig, &self.curr_pos, end
                        )
                    })?;
                let exp_chunks = self
                    .sample_index
                    .query_exp_chunks(&self.curr_contig, &interval)
                    .with_context(|| {
                        format!(
                            "failed to get exp chunks at {}:{}-{}",
                            &self.curr_contig, &self.curr_pos, end
                        )
                    })?;
                current_batch.add_chunks(
                    &self.curr_contig,
                    positions,
                    control_chunks,
                    exp_chunks,
                );
                current_batch_length += interval_length;
            }

            if current_batch_length >= self.interval_size {
                let finished_batch = std::mem::replace(
                    &mut current_batch,
                    DmrBatchOfPositions::empty(),
                );
                batches.push(finished_batch);
                current_batch_length = 0;
            }
            if end >= self.curr_contig_end {
                self.update()
            } else {
                self.curr_pos = end;
            }
        }

        if current_batch.num_chunks() > 0 {
            batches.push(current_batch);
        }

        if batches.is_empty() {
            Ok(None)
        } else {
            Ok(Some(batches))
        }
    }
}

impl Iterator for SingleSiteBatches {
    type Item = anyhow::Result<Vec<DmrBatchOfPositions>>;

    fn next(&mut self) -> Option<Self::Item> {
        self.get_next_batch().transpose()
    }
}

type SingleSiteDmrScoreResult = anyhow::Result<SingleSiteDmrScore>;
struct SingleSiteDmrScore {
    counts_a: AggregatedCounts,
    counts_b: AggregatedCounts,
    position: u64,
    score: f64,
    map_pval: f64,
    effect_size: f64,
    balanced_map_pval: f64,
    balanced_effect_size: f64,
    _balanced_score: f64,
    replicate_map_pval: Vec<f64>,
    replicate_effect_sizes: Vec<f64>,
    pct_a_samples: usize,
    pct_b_samples: usize,
}

impl SingleSiteDmrScore {
    fn header(multiple_samples: bool, matched_samples: bool) -> String {
        let mut fields = vec![
            "chrom",
            "start",
            "end",
            "name",
            "score",
            "a_counts",
            "a_total",
            "b_counts",
            "b_total",
            "a_mod_percentages",
            "b_mod_percentages",
            "a_pct_modified",
            "b_pct_modified",
            "map_pvalue",
            "effect_size",
        ];
        if multiple_samples {
            for field in [
                "balanced_map_pvalue",
                "balanced_effect_size",
                "pct_a_samples",
                "pct_b_samples",
            ] {
                fields.push(field)
            }
        }

        if matched_samples {
            for field in ["replicate_map_pvalues", "replicate_effect_sizes"] {
                fields.push(field)
            }
        }

        let mut s = fields.join("\t");
        s.push('\n');
        s
    }

    fn new_multi(
        counts_a: &[AggregatedCounts],
        counts_b: &[AggregatedCounts],
        sample_index: &SingleSiteSampleIndex,
        position: u64,
        estimator: &PMapEstimator,
    ) -> anyhow::Result<Self> {
        let (replicate_epmap, replicate_effect_sizes) = if sample_index
            .matched_replicate_samples()
            && counts_a.len() == counts_b.len()
        {
            let n_samples = counts_a.len();
            let mut replicate_epmap = Vec::with_capacity(n_samples);
            let mut replicate_effect_size = Vec::with_capacity(n_samples);
            for (a, b) in counts_a.iter().zip(counts_b) {
                let epmap = estimator.predict(a, b)?;
                replicate_epmap.push(epmap.e_pmap);
                replicate_effect_size.push(epmap.effect_size);
            }
            (replicate_epmap, replicate_effect_size)
        } else {
            (Vec::new(), Vec::new())
        };
        let pct_a_samples = ((counts_a.len() as f32
            / sample_index.num_a_samples() as f32)
            * 100f32)
            .floor() as usize;
        let pct_b_samples = ((counts_b.len() as f32
            / sample_index.num_b_samples() as f32)
            * 100f32)
            .floor() as usize;
        let balanced_counts_a = collapse_counts(counts_a, true);
        let balanced_counts_b = collapse_counts(counts_b, true);
        let epmap_balanced =
            estimator.predict(&balanced_counts_a, &balanced_counts_b)?;
        let balanced_llr_score =
            llk_ratio(&balanced_counts_a, &balanced_counts_b)?;
        let collapsed_a = collapse_counts(counts_a, false);
        let collapsed_b = collapse_counts(counts_b, false);
        let epmap = estimator.predict(&collapsed_a, &collapsed_b)?;
        let llr_score = llk_ratio(&collapsed_a, &collapsed_b)?;
        Ok(Self {
            counts_a: collapsed_a,
            counts_b: collapsed_b,
            position,
            score: llr_score,
            map_pval: epmap.e_pmap,
            effect_size: epmap.effect_size,
            balanced_map_pval: epmap_balanced.e_pmap,
            balanced_effect_size: epmap_balanced.effect_size,
            _balanced_score: balanced_llr_score,
            replicate_map_pval: replicate_epmap,
            replicate_effect_sizes,
            pct_a_samples,
            pct_b_samples,
        })
    }

    fn to_row(
        &self,
        multiple_samples: bool,
        matched_samples: bool,
        chrom: &str,
    ) -> String {
        let sep = '\t';
        if matched_samples {
            let replicate_map_pvals = if self.replicate_map_pval.is_empty() {
                "-".to_string()
            } else {
                self.replicate_map_pval.iter().join(",")
            };
            let replicate_effect_sizes =
                if self.replicate_effect_sizes.is_empty() {
                    "-".to_string()
                } else {
                    self.replicate_effect_sizes.iter().join(",")
                };
            format!(
                "\
            {}{sep}\
            {}{sep}\
            {}{sep}\
            {}{sep}\
            {}{sep}\
            {}{sep}\
            {}{sep}\
            {}{sep}\
            {}{sep}\
            {}{sep}\
            {}{sep}\
            {}{sep}\
            {}{sep}\
            {}{sep}\
            {}{sep}\
            {}{sep}\
            {}{sep}\
            {}{sep}\
            {}{sep}\
            {}{sep}\
            {}\n",
                chrom,
                self.position,
                self.position.saturating_add(1),
                '.',
                self.score,
                self.counts_a.string_counts(),
                self.counts_a.total,
                self.counts_b.string_counts(),
                self.counts_b.total,
                self.counts_a.string_percentages(),
                self.counts_b.string_percentages(),
                self.counts_a.pct_modified(),
                self.counts_b.pct_modified(),
                self.map_pval,
                self.effect_size,
                self.balanced_map_pval,
                self.balanced_effect_size,
                self.pct_a_samples,
                self.pct_b_samples,
                replicate_map_pvals,
                replicate_effect_sizes,
            )
        } else {
            self.to_row_pair(multiple_samples, chrom)
        }
    }

    fn to_row_pair(&self, multiple_samples: bool, chrom: &str) -> String {
        let sep = '\t';
        let row = format!(
            "\
            {}{sep}\
            {}{sep}\
            {}{sep}\
            {}{sep}\
            {}{sep}\
            {}{sep}\
            {}{sep}\
            {}{sep}\
            {}{sep}\
            {}{sep}\
            {}{sep}\
            {}{sep}\
            {}{sep}\
            {}{sep}\
            {}",
            chrom,
            self.position,
            self.position.saturating_add(1),
            '.',
            self.score,
            self.counts_a.string_counts(),
            self.counts_a.total,
            self.counts_b.string_counts(),
            self.counts_b.total,
            self.counts_a.string_percentages(),
            self.counts_b.string_percentages(),
            self.counts_a.pct_modified(),
            self.counts_b.pct_modified(),
            self.map_pval,
            self.effect_size,
        );
        let rest = if multiple_samples {
            format!(
                "\
                {sep}{}{sep}{}{sep}{}{sep}{}\n",
                self.balanced_map_pval,
                self.balanced_effect_size,
                self.pct_a_samples,
                self.pct_b_samples
            )
        } else {
            format!("\n")
        };

        format!("{row}{rest}")
    }
}

fn collapse_counts(
    counts: &[AggregatedCounts],
    balance: bool,
) -> AggregatedCounts {
    if counts.len() == 1 {
        counts[0].clone()
    } else if balance {
        let total_cov = counts.iter().map(|ac| ac.total).sum::<usize>();
        let n = counts.len();
        let target_cov = total_cov as f32 / n as f32;
        counts.iter().fold(AggregatedCounts::zero(), |agg, next| {
            let counts = next
                .iter_mod_fractions()
                .map(|(code, frac)| {
                    (code, (frac * target_cov).floor() as usize)
                })
                .collect::<HashMap<ModCodeRepr, usize>>();
            let total = target_cov as usize;
            let ac = AggregatedCounts::try_new(counts, total).unwrap();
            agg.op(&ac)
        })
    } else {
        counts.iter().fold(AggregatedCounts::zero(), |agg, next| agg.op(next))
    }
}

type ChromToSingleScores = (String, Vec<anyhow::Result<SingleSiteDmrScore>>);
fn process_batch_of_positions(
    batch: DmrBatchOfPositions,
    sample_index: Arc<SingleSiteSampleIndex>,
    pmap_estimator: Arc<PMapEstimator>,
) -> anyhow::Result<Vec<ChromToSingleScores>> {
    let (a_lines, b_lines) =
        sample_index.read_bedmethyl_lines_organized_by_position(batch)?;

    let counts = a_lines
        .into_iter()
        // intersect a_lines and b_lines on contig/chrom
        .filter_map(|(chrom, xs)| b_lines.get(&chrom).map(|ys| (chrom, xs, ys)))
        .map(|(chrom, xs, ys)| {
            // par iter over the positions
            let scores = xs
                .par_iter()
                .filter_map(|(pos, a_counts)| {
                    ys.get(&pos).map(|b_counts| (pos, a_counts, b_counts))
                })
                .map(|(pos, a_counts, b_counts)| {
                    // todo refactor this to be part of PMapEstimator
                    SingleSiteDmrScore::new_multi(
                        &a_counts,
                        &b_counts,
                        &sample_index,
                        pos.position,
                        &pmap_estimator,
                    )
                })
                .collect::<Vec<SingleSiteDmrScoreResult>>();
            (chrom, scores)
        })
        .collect::<Vec<_>>();

    Ok(counts)
}

struct Coverages {
    a_coverages: Vec<u64>,
    b_coverages: Vec<u64>,
}
fn get_coverage_from_batch(
    sample_index: &SingleSiteSampleIndex,
    batch: DmrBatchOfPositions,
) -> anyhow::Result<Coverages> {
    let (a_lines, b_lines) =
        sample_index.read_bedmethyl_lines_filtered_by_position(&batch)?;
    let get_covs = |x: SampleToBedMethyLines| -> Vec<u64> {
        x.into_iter()
            .flat_map(|(_sample_id, chrom_to_lines)| {
                chrom_to_lines
                    .values()
                    .flat_map(|ls| ls.iter().map(|l| l.valid_coverage))
                    .collect::<Vec<u64>>()
            })
            .collect::<Vec<u64>>()
    };
    let a_coverages = get_covs(a_lines);
    let b_coverages = get_covs(b_lines);
    Ok(Coverages { a_coverages, b_coverages })
}

fn calculate_max_coverages(
    sample_index: Arc<SingleSiteSampleIndex>,
    genome_positions: Arc<GenomePositions>,
    batch_size: usize,
    interval_size: u64,
    sample_n: usize,
    progress: &MultiProgress,
) -> anyhow::Result<[usize; 2]> {
    info!("estimating max coverages from data");
    let batch_iter = SingleSiteBatches::new(
        sample_index.clone(),
        genome_positions,
        batch_size,
        interval_size,
    )?;

    let pb = progress.add(get_subroutine_progress_bar(sample_n));
    pb.set_message("sampled bedMethyl records");
    let mut a_agg = Vec::with_capacity(sample_n);
    let mut b_agg = Vec::with_capacity(sample_n);
    for super_batch in batch_iter.filter_map(|r| match r {
        Ok(super_batch) => Some(super_batch),
        Err(e) => {
            debug!("batch failed during sampling, {e}");
            None
        }
    }) {
        let coverages = super_batch
            .into_par_iter()
            .filter_map(|batch_of_positions| {
                match get_coverage_from_batch(&sample_index, batch_of_positions)
                {
                    Ok(coverages) => Some(coverages),
                    Err(e) => {
                        debug!("failed to get coverages, {e}");
                        None
                    }
                }
            })
            .collect::<Vec<Coverages>>();
        for mut cov in coverages.into_iter() {
            a_agg.append(&mut cov.a_coverages);
            b_agg.append(&mut cov.b_coverages);
        }
        let n_so_far = std::cmp::min(a_agg.len(), b_agg.len());
        pb.set_position(n_so_far as u64);
        if n_so_far >= sample_n {
            break;
        }
    }

    pb.finish_and_clear();
    info!(
        "sampled {} a records and {} b records, calculating max coverages for \
         95th percentile",
        a_agg.len(),
        b_agg.len()
    );

    let a_coverages = a_agg
        .into_iter()
        .map(|x| x as f32)
        .sorted_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal))
        .collect::<Vec<f32>>();

    let a_max_cov =
        (percentile_linear_interp(&a_coverages, 0.95)?).floor() as usize;
    drop(a_coverages);

    let b_coverages = b_agg
        .into_iter()
        .map(|x| x as f32)
        .sorted_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal))
        .collect::<Vec<f32>>();
    let b_max_cov =
        (percentile_linear_interp(&b_coverages, 0.95)?).floor() as usize;

    info!("calculated max coverage for a: {a_max_cov} and b: {b_max_cov}");
    Ok([a_max_cov, b_max_cov])
}

trait DmrSegmenter {
    fn add(&mut self, dmr_scores: &[ChromToSingleScores])
        -> anyhow::Result<()>;
    fn run_current_chunk(&mut self) -> anyhow::Result<()>;
    fn clean_up(&mut self) -> anyhow::Result<()>;
}

#[derive(new)]
struct DummySegmenter {}

impl DmrSegmenter for DummySegmenter {
    fn add(
        &mut self,
        _dmr_scores: &[ChromToSingleScores],
    ) -> anyhow::Result<()> {
        Ok(())
    }

    fn run_current_chunk(&mut self) -> anyhow::Result<()> {
        Ok(())
    }

    fn clean_up(&mut self) -> anyhow::Result<()> {
        Ok(())
    }
}

struct HmmDmrSegmenter {
    writer: TsvWriter<File>,
    hmm: HmmModel,
    curr_region_scores: Vec<f64>,
    curr_region_positions: Vec<u64>,
    curr_counts_a: BTreeMap<u64, AggregatedCounts>,
    curr_counts_b: BTreeMap<u64, AggregatedCounts>,
    curr_chrom: Option<String>,
    curr_end: Option<u64>,
    max_gap_size: u64,
    size_gauge: ProgressBar,
    segments_written: ProgressBar,
}

impl DmrSegmenter for HmmDmrSegmenter {
    fn add(
        &mut self,
        dmr_scores: &[ChromToSingleScores],
    ) -> anyhow::Result<()> {
        for (chrom, scores) in dmr_scores.iter() {
            if let Some(curr_chrom) = self.curr_chrom.as_ref() {
                if chrom == curr_chrom {
                    let min_pos = scores.iter().find_map(|r| match r {
                        Ok(score) => Some(score.position),
                        Err(_) => None,
                    });
                    match (min_pos, self.curr_end) {
                        (Some(pos), Some(end)) => {
                            if pos
                                .checked_sub(end)
                                .map(|x| x < self.max_gap_size)
                                .unwrap_or(false)
                            {
                                // within limits, add to current
                                self.append_scores(&scores);
                            } else {
                                // next chunk is too far away, run current and
                                // reset
                                self.run_current_chunk()?;
                                self.append_scores(&scores);
                            }
                        }
                        (Some(_pos), None) => {
                            // maybe this never happens?
                            // don't have any data, append
                            self.append_scores(&scores);
                        }
                        (None, _) => {
                            // nothing to do
                            debug!("no valid results..");
                        }
                    }
                } else {
                    // finish current chunk and add this chunk to current
                    self.run_current_chunk()?;
                    // update chrom
                    self.curr_chrom = Some(chrom.to_string());
                    // update scores
                    assert_eq!(
                        self.curr_chrom.as_ref(),
                        Some(chrom),
                        "chroms arent' the same?"
                    );
                    self.append_scores(&scores);
                }
            } else {
                self.curr_chrom = Some(chrom.to_string());
                assert_eq!(
                    self.curr_chrom.as_ref(),
                    Some(chrom),
                    "chroms arent' the same?"
                );
                self.append_scores(&scores);
            }
        }
        Ok(())
    }

    #[inline]
    fn run_current_chunk(&mut self) -> anyhow::Result<()> {
        if self.curr_region_scores.is_empty() {
            debug!("no scores to run");
            assert!(
                self.curr_region_positions.is_empty(),
                "should not have positions and no scores"
            );
            return Ok(());
        }
        assert_eq!(
            self.curr_region_positions.len(),
            self.curr_region_scores.len(),
            "scores and positions should be the same length"
        );

        // these expects and asserts are safe because this method is only called
        // when self.curr_chrom is some
        let region =
            self.current_chunk_region().expect("region should not be None");
        assert!(self.curr_chrom.is_some());
        let start_time = std::time::Instant::now();
        let path = self.hmm.viterbi_path(
            &self.curr_region_scores,
            &self.curr_region_positions,
        );
        let took = start_time.elapsed();
        debug!(
            "segmenting {} ({} scores), took {took:?}",
            region.to_string(),
            self.curr_region_scores.len()
        );
        let integrated_path =
            path_to_region_labels(&path, &self.curr_region_positions);
        for (start, end, state) in integrated_path.iter() {
            let counts_a = self.get_counts_a(*start, *end);
            let counts_b = self.get_counts_b(*start, *end);
            let score = llk_ratio(&counts_a, &counts_b)?;
            let frac_mod_a = counts_a.pct_modified();
            let frac_mod_b = counts_b.pct_modified();
            let effect_size = frac_mod_a - frac_mod_b;
            let num_sites = self.curr_counts_a.range(*start..*end).count();

            let sep = '\t';
            let row = format!(
                "{}{sep}\
                {start}{sep}\
                {end}{sep}\
                {state}{sep}\
                {score}{sep}\
                {num_sites}{sep}\
                {}{sep}\
                {}{sep}\
                {}{sep}\
                {}{sep}\
                {frac_mod_a}{sep}\
                {frac_mod_b}{sep}\
                {effect_size}\n",
                self.curr_chrom.as_ref().unwrap(),
                counts_a.string_counts(),
                counts_b.string_counts(),
                counts_a.string_percentages(),
                counts_b.string_percentages(),
            );
            self.writer.write(row.as_bytes())?;
        }
        debug!("wrote {} segments", integrated_path.len());

        // reset everything
        self.curr_region_positions = Vec::new();
        self.curr_region_scores = Vec::new();
        self.curr_counts_a = BTreeMap::new();
        self.curr_counts_b = BTreeMap::new();
        self.curr_end = None;
        self.segments_written.inc(integrated_path.len() as u64);
        self.size_gauge.set_position(0u64);
        Ok(())
    }

    fn clean_up(&mut self) -> anyhow::Result<()> {
        self.size_gauge.finish_and_clear();
        self.segments_written.finish_and_clear();
        debug!(
            "HMM segmenter finished, wrote {} segments",
            self.segments_written.position()
        );
        Ok(())
    }
}

impl HmmDmrSegmenter {
    fn new(
        out_fp: &PathBuf,
        max_gap_size: u64,
        dmr_prior: f64,
        diff_stay: f64,
        same_state_factor: f64,
        diff_state_factor: f64,
        significance_factor: f64,
        linear_transitions: bool,
        decay_distance: u32,
        multi_progress: &MultiProgress,
    ) -> anyhow::Result<Self> {
        let hmm = HmmModel::new(
            dmr_prior,
            diff_stay,
            same_state_factor,
            diff_state_factor,
            significance_factor,
            decay_distance,
            linear_transitions,
        )?;
        let writer = TsvWriter::new_path(out_fp, true, None)?;
        let size_gauge = multi_progress.add(get_ticker());
        let segments_written = multi_progress.add(get_ticker());
        size_gauge.set_message("[segmenter] current region size");
        segments_written.set_message("[segmenter] segments finished");

        Ok(Self {
            writer,
            hmm,
            max_gap_size,
            curr_region_scores: Vec::new(),
            curr_region_positions: Vec::new(),
            curr_counts_a: BTreeMap::new(),
            curr_counts_b: BTreeMap::new(),
            curr_chrom: None,
            curr_end: None,
            size_gauge,
            segments_written,
        })
    }

    fn append_scores(&mut self, scores: &[anyhow::Result<SingleSiteDmrScore>]) {
        let mut rightmost = 0u64;
        for score in scores.iter().filter_map(|r| r.as_ref().ok()) {
            self.curr_region_scores.push(score.score);
            self.curr_region_positions.push(score.position);
            let check = self
                .curr_counts_a
                .insert(score.position, score.counts_a.clone());
            assert!(check.is_none());
            let check = self
                .curr_counts_b
                .insert(score.position, score.counts_b.clone());
            assert!(check.is_none());
            rightmost = std::cmp::max(rightmost, score.position);
        }
        // check, todo remove after testing
        if let Some(end) = self.curr_end {
            if rightmost > 0u64 {
                assert!(
                    end < rightmost,
                    "results were not sorted? {end} {rightmost}",
                );
            }
        }
        self.curr_end = Some(rightmost);
        self.size_gauge.set_position(self.curr_region_positions.len() as u64);
    }

    #[inline]
    fn current_chunk_start(&self) -> Option<&u64> {
        // todo can make this a .first()
        self.curr_region_positions.iter().min()
    }

    #[inline]
    fn current_chunk_region(&self) -> Option<Region> {
        match (
            self.curr_chrom.as_ref(),
            self.current_chunk_start(),
            self.curr_end,
        ) {
            (Some(chrom), Some(&start), Some(end)) => {
                Some(Region::new(chrom.to_string(), start as u32, end as u32))
            }
            _ => None,
        }
    }

    fn get_counts_a(&self, start: u64, stop: u64) -> AggregatedCounts {
        Self::get_counts_range(start..stop, &self.curr_counts_a)
    }

    fn get_counts_b(&self, start: u64, stop: u64) -> AggregatedCounts {
        Self::get_counts_range(start..stop, &self.curr_counts_b)
    }

    fn get_counts_range(
        r: Range<u64>,
        counts: &BTreeMap<u64, AggregatedCounts>,
    ) -> AggregatedCounts {
        counts.range(r).map(|(_, counts)| counts).fold(
            AggregatedCounts::zero(),
            |mut agg, x| {
                agg.op_mut(x);
                agg
            },
        )
    }
}

fn path_to_region_labels(
    path: &[States],
    positions: &[u64],
) -> Vec<(u64, u64, States)> {
    assert_eq!(path.len(), positions.len() - 1);
    if path.is_empty() {
        return Vec::new();
    } else {
        let mut curr_state = *path.first().unwrap();
        let mut curr_position = *positions.first().unwrap();
        let mut last_position = curr_position + 1;
        let mut agg = Vec::new();
        for (state, &pos) in path.iter().zip(positions).skip(1) {
            let position = pos;
            if state != &curr_state {
                let bedline = (curr_position, last_position, curr_state);
                agg.push(bedline);
                curr_position = position;
                last_position = position + 1;
                curr_state = *state;
            } else {
                last_position = position + 1;
            }
        }
        let final_bedline = (curr_position, last_position, curr_state);
        agg.push(final_bedline);

        agg
    }
}
