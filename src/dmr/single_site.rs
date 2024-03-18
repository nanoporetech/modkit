use std::cmp::Ordering;
use std::collections::{HashMap, VecDeque};
use std::io::Write;
use std::ops::Range;
use std::sync::Arc;

use anyhow::{bail, Context};
use indicatif::MultiProgress;
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
use crate::mod_base_code::ModCodeRepr;
use crate::monoid::BorrowingMoniod;
use crate::thresholds::percentile_linear_interp;
use crate::util::{get_subroutine_progress_bar, get_ticker};

pub(super) struct SingleSiteDmrAnalysis {
    sample_index: Arc<SingleSiteSampleIndex>,
    genome_positions: Arc<GenomePositions>,
    pmap_estimator: Arc<PMapEstimator>,
    batch_size: usize,
    interval_size: u64,
    header: bool,
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
        })
    }

    pub(super) fn run(
        &self,
        multi_progress_bar: MultiProgress,
        pool: rayon::ThreadPool,
        mut writer: Box<dyn Write>,
    ) -> anyhow::Result<()> {
        let matched_samples = self.sample_index.matched_replicate_samples();
        if self.header {
            writer.write(
                SingleSiteDmrScore::header(matched_samples).as_bytes(),
            )?;
        }

        let batch_iter = SingleSiteBatches::new(
            self.sample_index.clone(),
            self.genome_positions.clone(),
            self.batch_size,
            self.interval_size,
        )?;

        let (snd, rcv) = crossbeam::channel::bounded(1000);
        let processed_batches = multi_progress_bar.add(get_ticker());
        let failure_counter = multi_progress_bar.add(get_ticker());
        let success_counter = multi_progress_bar.add(get_ticker());

        processed_batches.set_message("batches processed");
        failure_counter.set_message("sites failed");
        success_counter.set_message("sites processed successfully");

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
                                match snd.send(chrom_to_scores) {
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
                results.into_iter().for_each(|chrom_to_scores| {
                    match snd.send(chrom_to_scores) {
                        Ok(_) => processed_batches.inc(1),
                        Err(e) => {
                            error!("failed to send on channel, {e}");
                        }
                    }
                });
            }
        });

        let mut success_count = 0usize;
        for batch_result in rcv {
            for (chrom, results) in batch_result {
                for result in results {
                    match result {
                        Ok(scores) => {
                            writer.write(scores.to_row(&chrom).as_bytes())?;
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

        success_counter.finish_and_clear();
        failure_counter.finish_and_clear();

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
    replicate_map_pval: Vec<f64>,
    replicate_effect_sizes: Vec<f64>,
}

impl SingleSiteDmrScore {
    fn header(matched_samples: bool) -> String {
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
        if matched_samples {
            for field in [
                "balanced_map_pvalue",
                "balanced_effect_size",
                "replicate_map_pvalues",
                "replicate_effect_sizes",
            ] {
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
        position: u64,
        estimator: &PMapEstimator,
        matched_samples: bool,
    ) -> anyhow::Result<Self> {
        let (replicate_epmap, replicate_effect_sizes) = if matched_samples {
            assert_eq!(
                counts_a.len(),
                counts_b.len(),
                "matched samples need to be the same"
            );
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
        let balanced_counts_a = collapse_counts(counts_a, true);
        let balanced_counts_b = collapse_counts(counts_b, true);
        let epmap_balanced =
            estimator.predict(&balanced_counts_a, &balanced_counts_b)?;
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
            replicate_map_pval: replicate_epmap,
            replicate_effect_sizes,
        })
    }

    fn to_row(&self, chrom: &str) -> String {
        let sep = '\t';
        if self.replicate_map_pval.is_empty() {
            debug_assert!(
                self.replicate_effect_sizes.is_empty(),
                "shouldn't have effect sizes and not map-pvalues"
            );
            self.to_row_pair(chrom)
        } else {
            let replicate_map_pvals = self.replicate_map_pval.iter().join(",");
            let replicate_effect_sizes =
                self.replicate_effect_sizes.iter().join(",");

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
                replicate_map_pvals,
                replicate_effect_sizes,
            )
        }
    }

    fn to_row_pair(&self, chrom: &str) -> String {
        let sep = '\t';
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
        )
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
    let matched_samples = sample_index.matched_replicate_samples();
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
                    if a_counts.len() != sample_index.num_a_samples()
                        || b_counts.len() != sample_index.num_b_samples()
                    {
                        bail!(
                            "don't have counts from all samples, a counts: \
                             {}, b counts: {}",
                            a_counts.len(),
                            b_counts.len()
                        )
                    } else {
                        SingleSiteDmrScore::new_multi(
                            &a_counts,
                            &b_counts,
                            pos.position,
                            &pmap_estimator,
                            matched_samples,
                        )
                    }
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
