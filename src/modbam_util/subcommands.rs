use std::path::{Path, PathBuf};

use anyhow::bail;
use clap::{Args, Subcommand};
use indicatif::{MultiProgress, ParallelProgressIterator, ProgressDrawTarget};
use log::{debug, info, warn};
use rayon::prelude::*;
use rust_htslib::bam;
use rust_htslib::bam::{FetchDefinition, Read};

use crate::command_utils::{get_serial_reader, using_stream};
use crate::interval_chunks::{
    ReferenceIntervalsFeeder, TotalLength, WithPrevEnd,
};
use crate::logging::init_logging;
use crate::modbam_util::check_tags::ModTagViews;
use crate::monoid::Moniod;
use crate::reads_sampler::record_sampler::RecordSampler;
use crate::reads_sampler::sample_reads_from_interval;
use crate::reads_sampler::sampling_schedule::SamplingSchedule;
use crate::record_processor::{RecordProcessor, WithRecords};
use crate::util::{
    get_master_progress_bar, get_subroutine_progress_bar, get_targets, Region,
};

#[derive(Subcommand)]
pub enum EntryModBam {
    #[command(name = "check-tags")]
    CheckTags(EntryCheckTags),
    // SampleReads(EntrySampleReads),
}

impl EntryModBam {
    pub fn run(&self) -> anyhow::Result<()> {
        match self {
            EntryModBam::CheckTags(x) => x.run(),
            // EntryModBam::SampleReads(x) => x.run(),
        }
    }
}

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct EntryCheckTags {
    /// Input modBam, can be a path to a file or one of `-` or
    /// `stdin` to specify a stream from standard input.
    in_bam: String,
    /// Don't exit 1 when invalid records are found in the input.
    #[arg(long, default_value_t = false)]
    permissive: bool,
    /// Write output tables into this directory. The directory will be created
    /// if it doesn't exist.
    #[clap(help_heading = "IO Options")]
    #[arg(short = 'o', long)]
    out_dir: Option<PathBuf>,
    /// Force overwrite of previous output
    #[clap(help_heading = "IO Options")]
    #[arg(short = 'f', long, default_value_t = false)]
    force: bool,
    /// Prefix output files with this string.
    #[clap(help_heading = "IO Options")]
    #[arg(long)]
    prefix: Option<String>,
    /// Number of threads to use.
    #[clap(help_heading = "Compute Options")]
    #[arg(short, long, default_value_t = 4)]
    threads: usize,
    /// Perform a linear scan of the modBAM even if the index is found.
    #[clap(help_heading = "Compute Options")]
    #[arg(long, default_value_t = false, hide_short_help = true)]
    ignore_index: bool,
    /// When using regions, interval chunk size in base pairs to process
    /// concurrently. Smaller interval chunk sizes will use less memory but
    /// incur more overhead.
    #[clap(help_heading = "Compute Options")]
    #[arg(short = 'i', long, default_value_t = 5_000_000)]
    interval_size: u32,
    /// Specify a file for debug logs to be written to, otherwise ignore them.
    /// Setting a file is recommended.
    #[clap(help_heading = "Logging Options")]
    #[arg(long, alias = "log")]
    log_filepath: Option<PathBuf>,
    /// Hide the progress bar.
    #[clap(help_heading = "Logging Options")]
    #[arg(long, default_value_t = false, hide_short_help = true)]
    suppress_progress: bool,

    /// Approximate maximum number of reads to use, especially recommended when
    /// using a large BAM without an index. If an indexed BAM is provided, the
    /// reads will be sampled evenly over the length of the aligned reference.
    /// If a region is passed with the --region option, they will be sampled
    /// over the genomic region. Actual number of reads used may deviate
    /// slightly from this number.
    #[clap(help_heading = "Selection Options")]
    #[arg(short = 'n', long)]
    num_reads: Option<usize>,
    /// Check tags on non-primary alignments as well. Keep in mind this
    /// may incur a double-counting of the read with its primary mapping.
    #[clap(help_heading = "Selection Options")]
    #[arg(long, default_value_t = false)]
    allow_non_primary: bool,
    /// Only check alignments that are mapped.
    #[clap(help_heading = "Selection Options")]
    #[arg(long, default_value_t = false)]
    only_mapped: bool,

    /// Process only the specified region of the BAM when collecting
    /// probabilities. Format should be <chrom_name>:<start>-<end> or
    /// <chrom_name>.
    #[clap(help_heading = "Selection Options")]
    #[arg(long)]
    region: Option<String>,
}

impl EntryCheckTags {
    pub fn run(&self) -> anyhow::Result<()> {
        use super::check_tags::output_filenames as ofn;
        let _handle = init_logging(self.log_filepath.as_ref());
        let mut reader = get_serial_reader(&self.in_bam)?;
        if let Some(out_d) = self.out_dir.as_ref() {
            if !out_d.exists() {
                info!("creating directory at {out_d:?}");
                std::fs::create_dir_all(out_d)?;
            }
            for fp in ofn::filenames {
                let out_fn = if let Some(p) = self.prefix.as_ref() {
                    out_d.join(format!("{p}_{fp}"))
                } else {
                    out_d.join(fp)
                };
                if out_fn.exists() && !self.force {
                    bail!("refusing to overwrite {fp:?}, use --force");
                }
            }
        }

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()?;

        let region = self
            .region
            .as_ref()
            .map(|raw_region| Region::parse_str(raw_region, reader.header()))
            .transpose()?;
        if region.is_some() && using_stream(&self.in_bam) {
            bail!("cannot use --region and stream from stdin");
        }
        let schedule = match (self.num_reads, using_stream(&self.in_bam)) {
            (_, true) | (None, false) => None,
            (Some(num_reads), false) => {
                match bam::IndexedReader::from_path(&self.in_bam) {
                    Ok(_) => {
                        let bam_fp = Path::new(&self.in_bam).to_path_buf();
                        if !bam_fp.exists() {
                            bail!("failed to find ${bam_fp:?}");
                        }
                        Some(SamplingSchedule::from_num_reads(
                            &self.in_bam,
                            num_reads,
                            region.as_ref(),
                            None,
                            !self.only_mapped,
                        )?)
                    }
                    Err(_) => {
                        warn!(
                            "didn't find BAM index, using first (head) \
                             {num_reads} reads"
                        );
                        None
                    }
                }
            }
        };

        let linear_scan = self.ignore_index
            || using_stream(&self.in_bam)
            || bam::IndexedReader::from_path(&self.in_bam).is_err();
        let tag_views = pool.install(|| {
            if linear_scan {
                reader.set_threads(self.threads)?;
                let record_sampler = if let Some(n_reads) =
                    self.num_reads.as_ref()
                {
                    info!("checking tags on {n_reads} randomly sampled reads");
                    RecordSampler::new_num_reads(*n_reads)
                } else {
                    RecordSampler::new_passthrough()
                };
                ModTagViews::process_records(
                    reader.records(),
                    !self.suppress_progress,
                    record_sampler,
                    None,
                    None,
                    None,
                    self.only_mapped,
                    self.allow_non_primary,
                    None,
                    None,
                )
            } else {
                let bam_fp = Path::new(&self.in_bam);
                if !bam_fp.exists() {
                    bail!("failed to find modBAM at {bam_fp:?}");
                }

                let reference_records =
                    get_targets(reader.header(), region.as_ref());
                let feeder = ReferenceIntervalsFeeder::new(
                    reference_records,
                    (self.threads as f32 * 1.5f32).floor() as usize,
                    self.interval_size,
                    false,
                    None,
                    None,
                )?;
                pool.install(|| {
                    self.run_check_tags_indexed(
                        bam_fp.to_path_buf(),
                        feeder,
                        schedule,
                    )
                })
            }
        })?;

        tag_views.report(
            self.out_dir.as_ref(),
            self.prefix.as_ref(),
            self.force,
            self.permissive,
        )?;
        Ok(())
    }

    fn run_check_tags_indexed(
        &self,
        bam_fp: PathBuf,
        feeder: ReferenceIntervalsFeeder,
        schedule: Option<SamplingSchedule>,
    ) -> anyhow::Result<ModTagViews> {
        let mpb = MultiProgress::new();
        if self.suppress_progress {
            mpb.set_draw_target(ProgressDrawTarget::hidden());
        }
        let prog_length = feeder.total_length() as usize;
        let master_progress = mpb.add(get_master_progress_bar(prog_length));
        master_progress.set_message("genome positions");

        let feeder = feeder.map(|x| x.unwrap()).with_prev_end();
        let (snd, rcv) = crossbeam::channel::bounded(1000);
        let allow_non_primary = self.allow_non_primary;
        let only_mapped = self.only_mapped;
        let bam_fp1 = bam_fp.clone();
        let with_progress = !self.suppress_progress;
        let n_reads = self.num_reads;
        let threads = self.threads;
        std::thread::spawn(move || {
            let mut records_used = 0usize;
            for super_batch in feeder {
                let total_batch_length =
                    super_batch.iter().map(|c| c.total_length()).sum::<u64>();
                let batch_progress =
                    mpb.add(get_subroutine_progress_bar(super_batch.len()));
                batch_progress.set_message("batch progress");
                let results = super_batch
                    .into_par_iter()
                    .progress_with(batch_progress)
                    .map(|batch| {
                        batch
                            .0
                            .into_par_iter()
                            .filter(|cc| {
                                schedule
                                    .as_ref()
                                    .map(|s| s.chrom_has_reads(cc.chrom_tid()))
                                    .unwrap_or(true)
                            })
                            .map(|cc| {
                                let record_sampler = schedule
                                    .as_ref()
                                    .map(|s| {
                                        s.get_record_sampler(
                                            cc.chrom_tid(),
                                            total_batch_length as u32,
                                            cc.start_pos(),
                                            cc.end_pos(),
                                        )
                                    })
                                    .unwrap_or_else(|| {
                                        RecordSampler::new_passthrough()
                                    });
                                sample_reads_from_interval::<ModTagViews>(
                                    &bam_fp1,
                                    cc.chrom_tid(),
                                    cc.start_pos(),
                                    cc.end_pos(),
                                    cc.prev_end(),
                                    record_sampler,
                                    None,
                                    None,
                                    None,
                                    true,
                                    allow_non_primary,
                                    None,
                                )
                            })
                            .collect::<anyhow::Result<Vec<ModTagViews>>>()
                    })
                    .collect::<Vec<anyhow::Result<Vec<ModTagViews>>>>();
                let records_accumulated = results
                    .par_iter()
                    .filter_map(|res| {
                        if let Ok(views) = res {
                            Some(
                                views
                                    .iter()
                                    .map(|x| x.num_reads())
                                    .sum::<usize>(),
                            )
                        } else {
                            None
                        }
                    })
                    .sum::<usize>();

                records_used = records_used.saturating_add(records_accumulated);

                match snd.send(results) {
                    Ok(_) => {}
                    Err(e) => {
                        debug!("failed to send on channel, {e}");
                    }
                }
                master_progress.inc(total_batch_length);
            }
            if !only_mapped {
                debug!("processing unmapped reads");
                let sampler = n_reads
                    .map(|nr| nr.checked_sub(records_used).unwrap_or(0))
                    .map(|n| {
                        debug!("processing {n} unmapped reads");
                        RecordSampler::new_num_reads(n)
                    })
                    .unwrap_or_else(|| RecordSampler::new_passthrough());
                let reader = bam::IndexedReader::from_path(&bam_fp)
                    .and_then(|mut reader| {
                        reader.fetch(FetchDefinition::Unmapped).map(|_| reader)
                    })
                    .and_then(|mut reader| {
                        reader.set_threads(threads).map(|_| reader)
                    });
                match reader {
                    Ok(mut reader) => {
                        let unmapped_tag_views = ModTagViews::process_records(
                            reader.records(),
                            with_progress,
                            sampler,
                            None,
                            None,
                            None,
                            false,
                            true,
                            None,
                            None,
                        );
                        if let Ok(unmapped) = unmapped_tag_views {
                            match snd.send(vec![Ok(vec![unmapped])]) {
                                Ok(_) => {}
                                Err(e) => {
                                    debug!("failed to send on channel, {e}");
                                }
                            }
                        }
                    }
                    Err(e) => {
                        debug!("failed to collect unmapped records, {e}");
                    }
                }
            }
        });

        let tag_views = rcv
            .iter()
            .par_bridge()
            .inspect(|x| {
                let n_errs = x
                    .iter()
                    .filter(|y| y.is_err())
                    // .inspect(|x| {
                    //     dbg!(&x);
                    // })
                    .count();

                if n_errs > 0 {
                    debug!("{n_errs} batches failed");
                }
            })
            .map(|xs| {
                xs.into_iter()
                    .filter_map(|x| x.ok())
                    .flatten()
                    .collect::<Vec<ModTagViews>>()
            })
            .fold(
                || ModTagViews::default(),
                |mut a, b| {
                    b.into_iter().for_each(|x| a.op_mut(x));
                    a
                },
            )
            .reduce(|| ModTagViews::default(), |a, b| a.op(b));
        Ok(tag_views)
    }
}

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct EntrySampleReads {
    in_bam: String,
}

impl EntrySampleReads {
    pub fn run(&self) -> anyhow::Result<()> {
        todo!()
    }
}
