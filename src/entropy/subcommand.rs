use std::collections::HashMap;
use std::io::Write;
use std::path::PathBuf;

use anyhow::{anyhow, bail, Context};
use clap::Args;
use indicatif::{MultiProgress, ProgressBar};
use log::{debug, error, info};
use rayon::prelude::*;
use rust_htslib::bam::{self, Read};
use rustc_hash::FxHashMap;

use crate::command_utils::{
    get_threshold_from_options, parse_per_mod_thresholds,
};
use crate::entropy::{
    process_entropy_window, EntropyCalculation, SlidingWindows,
};
use crate::errs::RunError;
use crate::logging::init_logging;
use crate::mod_base_code::DnaBase;
use crate::motif_bed::{get_masked_sequences, RegexMotif};
use crate::reads_sampler::sampling_schedule::IdxStats;
use crate::threshold_mod_caller::MultipleThresholdModCaller;
use crate::util::{get_targets, get_ticker, Strand, TAB};
use crate::writers::TsvWriter;

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct MethylationEntropy {
    /// Input mod-BAM
    in_bam: PathBuf,
    /// Output BED
    #[arg(short = 'o', long)]
    out_bed: Option<PathBuf>,
    /// Number of modified positions to consider at a time
    #[arg(short = 'n', long, default_value_t = 4)]
    num_positions: usize,
    /// Maximum length interval that "num_positions" modified bases can occur
    /// in. The maximum window size decides how dense the positions are
    /// packed. For example, consider that the num_positions is equal to 4, the
    /// motif is CpG, and the window_size is equal to 8, this configuration
    /// would require that the modified positions are immediately adjacent
    /// to each other, "CGCGCGCG". On the other hand, if the window_size
    /// was set to 12, then multiple sequences with various patterns of
    /// other bases can be used CGACGATCGGCG.
    #[arg(short = 'w', long, default_value_t = 50)]
    window_size: usize,
    /// Do not perform any filtering, include all mod base calls in output.
    #[arg(group = "thresholds", long, default_value_t = false)]
    no_filtering: bool,
    /// Filter out modified base calls where the probability of the predicted
    /// variant is below this confidence percentile. For example, 0.1 will
    /// filter out the 10% lowest confidence modification calls.
    #[arg(
        group = "thresholds",
        short = 'p',
        long,
        default_value_t = 0.1,
        hide_short_help = true
    )]
    filter_percentile: f32,
    /// Specify the filter threshold globally or for the canonical calls.
    /// When specified, base modification call probabilities will be required
    /// to be greater than or equal to this number. If `--mod-thresholds`
    /// is also specified, _this_ value will be used for canonical calls.
    #[arg(
        long,
        group = "thresholds",
        action = clap::ArgAction::Append,
        alias = "pass_threshold"
    )]
    filter_threshold: Option<f32>,
    /// Specify a passing threshold to use for a base modification, independent
    /// of the threshold for the primary sequence base or the default. For
    /// example, to set the pass threshold for 5hmC to 0.8 use
    /// `--mod-threshold h:0.8`. The pass threshold will still be estimated
    /// as usual and used for canonical cytosine and other modifications
    /// unless the `--filter-threshold` option is also passed.
    /// See the online documentation for more details.
    #[arg(
        long,
        alias = "mod-threshold",
        action = clap::ArgAction::Append
    )]
    mod_thresholds: Option<Vec<String>>,
    /// Number of threads to use.
    #[arg(short = 't', long, default_value_t = 4)]
    threads: usize,
    /// Number of BAM-reading threads to use.
    #[arg(long, hide_short_help = true)]
    io_threads: Option<usize>,
    /// Reference sequence in FASTA format.
    #[arg(long = "ref", alias = "reference")]
    reference_fasta: PathBuf,
    /// Motif to use for entropy calculation, default will be CpG. The motif
    /// must be reverse-complement palindromic.
    #[arg(long, num_args = 2)]
    motif: Option<Vec<String>>,
    /// Use CpG motifs. Short hand for --motif CG 0 --combine-strands
    #[arg(long, default_value_t = false)]
    cpg: bool,
    /// Primary sequence base to calculate modification entropy on.
    #[arg(long, conflicts_with = "motif")]
    base: Option<DnaBase>,
    /// Regions over which to calculate descriptive statistics
    #[arg(long = "regions")]
    regions_fp: Option<PathBuf>,
    /// Combine modification counts on the positive and negative strands and
    /// report entropy on just the positive strand.
    #[arg(long, conflicts_with_all=["base", "cpg"], default_value_t=false)]
    combine_strands: bool,
    /// Minimum coverage required at each position in the window. Windows
    /// without at least this many valid reads will be skipped, but
    /// positions within the window with enough coverage can be used by
    /// neighboring windows.
    #[arg(long = "min-coverage", default_value_t = 3)]
    min_valid_coverage: u32,
    /// Respect soft masking in the reference FASTA.
    #[arg(
        long,
        short = 'k',
        requires = "reference_fasta",
        default_value_t = false,
        hide_short_help = true
    )]
    mask: bool,
    /// Send debug logs to this file, setting this file is recommended.
    #[arg(long, alias = "log")]
    log_filepath: Option<PathBuf>,
    /// Hide progress bars
    #[arg(long, hide_short_help = true, default_value_t = false)]
    suppress_progress: bool,
    /// Force overwrite output
    #[arg(long, default_value_t = false)]
    force: bool,
    /// Write a header line
    #[arg(long, default_value_t = false)]
    header: bool,
    /// Omit windows with zero entropy
    #[arg(long, conflicts_with = "regions_fp", default_value_t = false)]
    drop_zeros: bool,
    /// Maximum number of filtered positions a read is allowed to have in a
    /// window, more than this number and the read will be discarded. Default
    /// will be 50% of `num_positions`.
    #[arg(long)]
    max_filtered_positions: Option<usize>,
}

impl MethylationEntropy {
    pub fn run(&self) -> anyhow::Result<()> {
        let _handle = init_logging(self.log_filepath.as_ref());
        if self.num_positions == 0 {
            bail!("num-positions must be at least 1")
        }
        if self.min_valid_coverage < 1 {
            bail!("min-valid-coverage must be at least 1")
        }
        IdxStats::check_any_mapped_reads(&self.in_bam, None, None).context(
            "did not find any mapped reads in mod-BAM, perform alignment first",
        )?;

        let mut writer: Box<dyn EntropyWriter> = if let Some(out_path) =
            self.out_bed.as_ref()
        {
            if let Some(p) = out_path.parent() {
                if !p.exists() && !(p == std::path::Path::new("")) {
                    info!("creating output directory {p:?}");
                    std::fs::create_dir_all(p)?;
                }
            }
            let tsv_writer = TsvWriter::new_path(out_path, self.force, None)?;
            Box::new(tsv_writer)
        } else {
            Box::new(TsvWriter::new_stdout(None))
        };
        if self.header {
            if self.regions_fp.is_some() {
                writer.write_line(&format!("\
                    chrom{TAB}\
                    start{TAB}\
                    end{TAB}\
                    region_name{TAB}\
                    mean_entropy{TAB}\
                    strand{TAB}\
                    median_entropy{TAB}\
                    max_entropy{TAB}\
                    min_entropy{TAB}\
                    mean_num_reads{TAB}\
                    min_num_reads{TAB}\
                    max_num_reads{TAB}\
                    successful_window_count{TAB}\
                    failed_window_count\n"))?;
            } else {
                writer.write_line(&format!("\
                    chrom{TAB}\
                    start{TAB}\
                    end{TAB}\
                    entropy{TAB}\
                    strand{TAB}\
                    num_reads\n"))?;
            }
        }

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()?;
        let multi_pb = MultiProgress::new();
        if self.suppress_progress {
            multi_pb.set_draw_target(indicatif::ProgressDrawTarget::hidden());
        }
        let bam_reader = bam::IndexedReader::from_path(&self.in_bam)?;
        let chrom_id_to_name = bam_reader
            .header()
            .target_names()
            .iter()
            .enumerate()
            .map(|(id, raw_name)| {
                let chrom_id = id as u32;
                String::from_utf8(raw_name.to_vec())
                    .map(|name| (chrom_id, name))
            })
            .collect::<Result<Vec<(u32, String)>, _>>()?
            .into_iter()
            .collect::<HashMap<u32, String>>();

        let bam_header_records = get_targets(bam_reader.header(), None);

        let (motif, combine_strands) =
            match (self.cpg, self.motif.as_ref(), self.base) {
                (true, _, _) => {
                    info!("using CpG motif");
                    (RegexMotif::parse_string("CG", 0).unwrap(), true)
                }
                (false, Some(raw_motif_parts), _) => {
                    let mut motifs =
                        RegexMotif::from_raw_parts(raw_motif_parts, false)?;
                    if motifs.len() == 1 {
                        let motif = motifs.remove(0);
                        if self.combine_strands && !motif.is_palendrome() {
                            bail!(
                                "motif must be reverse-complement palindromic \
                                 to combine strands"
                            )
                        }
                        info!("using user-specified motif {}", motif);
                        (motif, self.combine_strands)
                    } else {
                        bail!("cannot have more than 1 motif")
                    }
                }
                (false, None, Some(dna_base)) => {
                    if self.combine_strands {
                        bail!(
                            "cannot combine strands with single base \
                             modifications"
                        )
                    }
                    (
                        RegexMotif::parse_string(
                            &format!("{}", dna_base.char()),
                            0,
                        )?,
                        false,
                    )
                }
                _ => bail!(
                    "invalid input options, must provide --motif, --base, or \
                     specify --cpg"
                ),
            };

        let batch_size = (self.threads as f32 * 1.5f32).floor() as usize;
        let window_size = self.window_size;
        info!(
            "window size is set to {}, motif ({motif:?}) length is {}",
            window_size,
            motif.length()
        );
        if combine_strands {
            info!("combining (+)-strand and (-)-strand modification calls");
        }

        let sliding_windows = pool.install(|| {
            let names_to_tid = bam_header_records
                .iter()
                .map(|ref_record| (ref_record.name.as_str(), ref_record.tid))
                .collect::<HashMap<&str, u32>>();
            let reference_sequences = get_masked_sequences(
                &self.reference_fasta,
                &names_to_tid,
                self.mask,
                &multi_pb,
            )?
            .into_par_iter()
            .map(|(seq, tid)| (tid, seq.chars().collect::<Vec<char>>()))
            .collect::<HashMap<u32, Vec<char>>>();

            if let Some(regions_fp) = self.regions_fp.as_ref() {
                let tid_to_name = names_to_tid
                    .iter()
                    .map(|(name, tid)| (*tid, *name))
                    .collect::<FxHashMap<u32, &str>>();

                let reference_sequences = reference_sequences
                    .into_iter()
                    .filter_map(|(tid, seq)| {
                        tid_to_name.get(&tid).map(|name| (*name, seq))
                    })
                    .collect::<HashMap<&str, Vec<char>>>();

                SlidingWindows::new_with_regions(
                    &names_to_tid,
                    reference_sequences,
                    regions_fp,
                    motif,
                    combine_strands,
                    self.num_positions,
                    window_size,
                    batch_size,
                )
            } else {
                SlidingWindows::new(
                    bam_header_records,
                    reference_sequences,
                    motif,
                    combine_strands,
                    self.num_positions,
                    window_size,
                    batch_size,
                )
            }
        })?;

        let threshold_caller = self.get_threshold_caller(&pool)?;

        let (snd, rcv) = crossbeam::channel::bounded(10_000);
        let rows_written = multi_pb.add(get_ticker());
        let skipped_windows = multi_pb.add(get_ticker());
        let windows_failed = multi_pb.add(get_ticker());
        let batches_failed = multi_pb.add(get_ticker());

        let what =
            if self.regions_fp.is_some() { "regions" } else { "windows" };
        rows_written.set_message("rows written");
        windows_failed.set_message(format!("{what} failed"));
        skipped_windows.set_message(format!("{what} with zero coverage"));
        batches_failed.set_message("batches failed");

        let bam_fp = self.in_bam.clone();
        let min_coverage = self.min_valid_coverage;
        let threads = self.threads;
        let io_threads = self.io_threads.unwrap_or(threads);
        let max_filtered = self.max_filtered_positions.unwrap_or_else(|| {
            let max_filt_pos =
                (self.num_positions as f32 * 0.5f32).floor() as usize;
            info!("setting maximum filtered positions to {max_filt_pos}");
            max_filt_pos
        });
        pool.spawn(move || {
            for batch in sliding_windows {
                let mut results = Vec::new();
                let (entropies, _) = rayon::join(
                    || {
                        batch
                            .into_par_iter()
                            .map(|window| {
                                process_entropy_window(
                                    window,
                                    min_coverage,
                                    max_filtered,
                                    io_threads,
                                    &threshold_caller,
                                    &bam_fp,
                                )
                            })
                            .collect::<Vec<_>>()
                    },
                    || {
                        results.into_iter().for_each(|entropy| {
                            match snd.send(entropy) {
                                Ok(_) => {}
                                Err(e) => {
                                    error!("failed to send on channel, {e}");
                                }
                            }
                        })
                    },
                );
                results = entropies;
                results.into_iter().for_each(|entropy| {
                    match snd.send(entropy) {
                        Ok(_) => {}
                        Err(e) => {
                            error!("failed to send on channel, {e}");
                        }
                    }
                });
            }
        });

        for batch_result in rcv.iter() {
            match batch_result {
                Ok(entropy_calculation) => {
                    writer.write(
                        entropy_calculation,
                        &chrom_id_to_name,
                        self.drop_zeros,
                        &rows_written,
                        &skipped_windows,
                        &windows_failed,
                    )?;
                }
                Err(e) => {
                    debug!("batch failed, {e}");
                    batches_failed.inc(1);
                }
            }
        }

        multi_pb.clear()?;
        info!(
            "finished, {} {what} processed successfully, {} windows with zero \
             coverage, {} windows failed",
            rows_written.position(),
            skipped_windows.position(),
            windows_failed.position()
        );

        Ok(())
    }

    fn get_threshold_caller(
        &self,
        pool: &rayon::ThreadPool,
    ) -> anyhow::Result<MultipleThresholdModCaller> {
        let per_mod_thresholds = self
            .mod_thresholds
            .as_ref()
            .map(|raw_per_mod_thresholds| {
                parse_per_mod_thresholds(raw_per_mod_thresholds)
            })
            .transpose()?;
        if let Some(base_threshold) = self.filter_threshold {
            info!("using threshold {base_threshold}");
            if let Some(mod_thresholds) = per_mod_thresholds.as_ref() {
                mod_thresholds.iter().for_each(|(code, val)| {
                    info!("using threshold value {val} for mod-code {code}")
                });
            }
            Ok(MultipleThresholdModCaller::new(
                HashMap::new(),
                per_mod_thresholds.unwrap_or(HashMap::new()),
                base_threshold,
            ))
        } else {
            pool.install(|| {
                get_threshold_from_options(
                    &self.in_bam,
                    self.threads,
                    1_000_000,
                    None,
                    1042,
                    self.no_filtering,
                    self.filter_percentile,
                    None,
                    None,
                    per_mod_thresholds,
                    None,
                    None,
                    None,
                    true,
                    self.suppress_progress,
                )
            })
        }
    }
}

trait EntropyWriter {
    fn write(
        &mut self,
        entropy_calculation: EntropyCalculation,
        chrom_id_to_name: &HashMap<u32, String>,
        drop_zeros: bool,
        write_counter: &ProgressBar,
        skip_counter: &ProgressBar,
        failure_counter: &ProgressBar,
        // todo failure causes
    ) -> anyhow::Result<()>;

    fn write_line(&mut self, line: &str) -> anyhow::Result<()>;
}

impl<T: Write> EntropyWriter for TsvWriter<T> {
    fn write(
        &mut self,
        entropy_calculation: EntropyCalculation,
        chrom_id_to_name: &HashMap<u32, String>,
        drop_zeros: bool,
        write_counter: &ProgressBar,
        skip_counter: &ProgressBar,
        failure_counter: &ProgressBar,
    ) -> anyhow::Result<()> {
        match entropy_calculation {
            EntropyCalculation::Windows(entropy_windows) => {
                for entropy in entropy_windows {
                    let name = chrom_id_to_name
                        .get(&entropy.chrom_id)
                        .ok_or_else(|| {
                            anyhow!(
                                "missing chrom name for {}",
                                &entropy.chrom_id
                            )
                        })?;
                    match entropy.pos_me_entropy.as_ref() {
                        Ok(pos_entropy) => {
                            if (drop_zeros && !(pos_entropy.me_entropy == 0f32))
                                || !drop_zeros
                            {
                                let row = format!(
                                    "\
                                {name}\t{}\t{}\t{}\t{}\t{}\n",
                                    pos_entropy.interval.start,
                                    pos_entropy.interval.end,
                                    pos_entropy.me_entropy,
                                    Strand::Positive.to_char(),
                                    pos_entropy.num_reads
                                );
                                self.write(&row.as_bytes())?;
                                write_counter.inc(1);
                            }
                        }
                        Err(e) => {
                            match e {
                                RunError::Failed(e) => {
                                    debug!("(+) window failed, {e}");
                                    failure_counter.inc(1);
                                }
                                RunError::BadInput(reason) => {
                                    debug!(
                                        "(+) window bad input?, {}",
                                        &reason.0
                                    );
                                    failure_counter.inc(1);
                                }
                                RunError::Skipped(_e) => {
                                    skip_counter.inc(1);
                                    // debug!("window {}:{}-{} skipped, {e}",
                                    // name, entropy.interval.start,
                                    // entropy.interval.end);
                                }
                            }
                        }
                    }

                    match entropy.neg_me_entropy.as_ref() {
                        Some(Ok(neg_entropy)) => {
                            if (drop_zeros && !(neg_entropy.me_entropy == 0f32))
                                || !drop_zeros
                            {
                                let row = format!(
                                    "\
                                    {name}\t{}\t{}\t{}\t{}\t{}\n",
                                    neg_entropy.interval.start,
                                    neg_entropy.interval.end,
                                    neg_entropy.me_entropy,
                                    Strand::Negative.to_char(),
                                    neg_entropy.num_reads
                                );
                                self.write(&row.as_bytes())?;
                                write_counter.inc(1);
                            }
                        }
                        Some(Err(e)) => {
                            match e {
                                RunError::Failed(e) => {
                                    debug!("(-) window failed, {e}");
                                    failure_counter.inc(1);
                                }
                                RunError::BadInput(reason) => {
                                    debug!(
                                        "(-) window bad input?, {}",
                                        &reason.0
                                    );
                                    failure_counter.inc(1);
                                }
                                RunError::Skipped(_e) => {
                                    skip_counter.inc(1);
                                    // debug!("window {}:{}-{} skipped, {e}",
                                    // name, entropy.interval.start,
                                    // entropy.interval.end);
                                }
                            }
                        }
                        None => {}
                    }
                }
            }

            EntropyCalculation::Region(region_entropy) => {
                let chrom =
                    chrom_id_to_name.get(&region_entropy.chrom_id).expect(
                        "shouldn't have a result on a chrom without a chromId",
                    );
                let start = region_entropy.interval.start;
                let end = region_entropy.interval.end;
                let region_name = region_entropy.region_name;
                match region_entropy.pos_entropy_stats {
                    Ok(pos_entropy_stats) => {
                        let row = pos_entropy_stats.to_row(
                            &chrom,
                            start,
                            end,
                            Strand::Positive,
                            &region_name,
                        );
                        self.write(row.as_bytes())?;
                        write_counter.inc(1);
                    }
                    Err(e) => match e {
                        RunError::Failed(e) => {
                            debug!("(+) region failed, {e}");
                            failure_counter.inc(1);
                        }
                        RunError::BadInput(reason) => {
                            debug!("(+) region bad input?, {}", &reason.0);
                            failure_counter.inc(1);
                        }
                        RunError::Skipped(_e) => {
                            skip_counter.inc(1);
                            // debug!("window {}:{}-{} skipped, {e}", name,
                            // entropy.interval.start, entropy.interval.end);
                        }
                    },
                }
                match region_entropy.neg_entropy_stats {
                    Some(Ok(neg_entropy_stats)) => {
                        let row = neg_entropy_stats.to_row(
                            &chrom,
                            start,
                            end,
                            Strand::Negative,
                            &region_name,
                        );
                        self.write(row.as_bytes())?;
                        write_counter.inc(1);
                    }
                    Some(Err(e)) => match e {
                        RunError::Failed(e) => {
                            debug!("(-) region failed, {e}");
                            failure_counter.inc(1);
                        }
                        RunError::BadInput(reason) => {
                            debug!("(-) region bad input?, {}", &reason.0);
                            failure_counter.inc(1);
                        }
                        RunError::Skipped(_e) => {
                            skip_counter.inc(1);
                            // debug!("window {}:{}-{} skipped, {e}", name,
                            // entropy.interval.start, entropy.interval.end);
                        }
                    },
                    None => {}
                }
            }
        }

        // match (entropy.neg_me_entropy, entropy.neg_num_reads) {
        //     (Some(ne_entropy), Some(n)) => {
        //         if (drop_zeros && !(ne_entropy == 0f32)) || !drop_zeros {
        //             let row = format!(
        //                 "\
        //                 {name}\t{}\t{}\t{}\t{}\n",
        //                 entropy.interval.start,
        //                 entropy.interval.end,
        //                 ne_entropy,
        //                 n
        //             );
        //             self.write(&row.as_bytes())?;
        //         }
        //     },
        //     _ => {}
        // }

        Ok(())
    }

    fn write_line(&mut self, line: &str) -> anyhow::Result<()> {
        self.write(line.as_bytes())?;
        Ok(())
    }
}
