use crate::interval_chunks::IntervalChunks;
use crate::monoid::Moniod;
use crate::util::{
    get_master_progress_bar, get_spinner, get_subroutine_progress_bar,
    get_targets, Region,
};
use anyhow::Result as AnyhowResult;
use indicatif::{MultiProgress, ParallelProgressIterator, ProgressBar};
use log::{debug, info};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use rust_htslib::bam::{self, Read};
use std::path::PathBuf;

pub struct RecordSampler {
    pub(crate) num_reads: Option<usize>,
    pub(crate) sample_frac: Option<f64>,
    pub(crate) seed: Option<u64>,
    rng: StdRng,
    reads_sampled: usize,
}

impl RecordSampler {
    pub(crate) fn new_num_reads(num_reads: usize) -> Self {
        Self {
            num_reads: Some(num_reads),
            sample_frac: None,
            seed: None,
            rng: StdRng::from_entropy(),
            reads_sampled: 0,
        }
    }

    pub(crate) fn new_sample_frac(sample_frac: f64, seed: Option<u64>) -> Self {
        let rng = seed
            .map(|s| StdRng::seed_from_u64(s))
            .unwrap_or(StdRng::from_entropy());
        Self {
            num_reads: None,
            sample_frac: Some(sample_frac),
            seed,
            rng,
            reads_sampled: 0,
        }
    }

    pub(crate) fn new_passthrough() -> Self {
        Self {
            num_reads: None,
            sample_frac: None,
            seed: None,
            rng: StdRng::from_entropy(),
            reads_sampled: 0,
        }
    }

    pub(crate) fn new_from_options(
        sample_frac: Option<f64>,
        num_reads: Option<usize>,
        seed: Option<u64>,
    ) -> Self {
        match (sample_frac, num_reads) {
            (_, Some(num_reads)) => RecordSampler::new_num_reads(num_reads),
            (Some(f), _) => RecordSampler::new_sample_frac(f, seed),
            (None, None) => RecordSampler::new_passthrough(),
        }
    }

    pub(crate) fn get_progress_bar(&self) -> ProgressBar {
        let spinner = if let Some(num) = self.num_reads {
            get_master_progress_bar(num)
        } else {
            get_spinner()
        };
        spinner.set_message("records sampled");
        spinner
    }

    fn check_num_reads(&mut self) -> Indicator {
        let indicator = if self.reads_sampled >= self.num_reads.unwrap() {
            Indicator::Done
        } else {
            Indicator::Use
        };
        self.reads_sampled += 1;
        indicator
    }

    fn check_sample_frac(&mut self) -> Indicator {
        if self.rng.gen_bool(self.sample_frac.unwrap()) {
            Indicator::Use
        } else {
            Indicator::Skip
        }
    }

    pub(crate) fn ask(&mut self) -> Indicator {
        match (self.num_reads, self.sample_frac) {
            (Some(_nr), _) => self.check_num_reads(),
            (_, Some(_sample_frac)) => self.check_sample_frac(),
            (None, None) => Indicator::Use,
        }
    }
}

pub enum Indicator {
    Use,
    Skip,
    Done,
}

pub trait IntervalProcessor {
    type Output;
    fn new(
        sample_frac: Option<f64>,
        num_reads: Option<usize>,
        seed: Option<u64>,
    ) -> Self;

    fn process_records<T: Read>(
        &mut self,
        records: bam::Records<T>,
        chrom_length: u32,
        total_length: u64,
        start: u64,
        end: u64,
    ) -> AnyhowResult<Self::Output>;

    fn label() -> &'static str;
}

#[inline]
fn process_interval_records<T: IntervalProcessor>(
    bam_fp: &PathBuf,
    chrom_tid: u32,
    start: u32,
    end: u32,
    total_length: u64,
    mut processor: T,
) -> AnyhowResult<T::Output> {
    let mut bam_reader = bam::IndexedReader::from_path(bam_fp)?;
    bam_reader.fetch(bam::FetchDefinition::Region(
        chrom_tid as i32,
        start as i64,
        end as i64,
    ))?;
    let chrom_length = end - start;
    processor.process_records(
        bam_reader.records(),
        chrom_length,
        total_length,
        start as u64,
        end as u64,
    )
}

pub fn run_sampled_region_processor<T: IntervalProcessor>(
    bam_fp: &PathBuf,
    threads: usize,
    interval_size: u32,
    sample_frac: Option<f64>,
    num_reads: Option<usize>,
    seed: Option<u64>,
    region: Option<&Region>,
) -> AnyhowResult<T::Output>
where
    T::Output: Send + Moniod,
{
    let reader = bam::IndexedReader::from_path(bam_fp)?;
    let header = reader.header();
    let references = get_targets(header, region);
    let total_length = references.iter().map(|r| r.length as u64).sum::<u64>();

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()?;

    let mut output_aggregator = T::Output::zero();
    let master_progress = MultiProgress::new();
    let tid_progress =
        master_progress.add(get_master_progress_bar(references.len()));
    tid_progress.set_message("contigs");
    let sampled_items = master_progress.add(get_spinner());
    sampled_items.set_message(format!("{} sampled", T::label()));

    pool.install(|| {
        for reference in references {
            let intervals = IntervalChunks::new(
                reference.start,
                reference.length,
                interval_size,
                reference.tid,
                None,
            )
            .collect::<Vec<(u32, u32)>>();
            let num_reads = num_reads.map(|nr| {
                let f = reference.length as f64 / total_length as f64;
                let nr = nr as f64 * f;
                std::cmp::max(nr.floor() as usize, 1usize)
            });

            let interval_progress = master_progress
                .add(get_subroutine_progress_bar(intervals.len()));
            interval_progress
                .set_message(format!("processing {}", &reference.name));
            let outputs = intervals
                .into_par_iter()
                .progress_with(interval_progress)
                .filter_map(|(start, end)| {
                    let processor = T::new(sample_frac, num_reads, seed);
                    match process_interval_records(
                        bam_fp,
                        reference.tid,
                        start,
                        end,
                        reference.length as u64,
                        processor,
                    ) {
                        Ok(ps) => Some(ps),
                        Err(er) => {
                            debug!(
                                // todo improve error message
                                "error sampling probs for region: {}",
                                er.to_string()
                            );
                            None
                        }
                    }
                })
                .reduce(|| T::Output::zero(), |a, b| a.op(b));
            tid_progress.inc(1);
            sampled_items.inc(outputs.len() as u64);
            output_aggregator.op_mut(outputs);
        }
    });

    tid_progress.finish_and_clear();
    info!("sampled {} {}", output_aggregator.len(), T::label());
    Ok(output_aggregator)
}
