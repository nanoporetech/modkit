use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Arc;

use anyhow::{bail, Context};
use clap::Args;
use indicatif::MultiProgress;
use log::{debug, error, info};
use rayon::prelude::*;

use crate::command_utils::parse_per_mod_thresholds;
use crate::entropy::writers::{EntropyWriter, RegionsWriter, WindowsWriter};
use crate::entropy::{process_entropy_window, SlidingWindows};
use crate::find_motifs::motif_bed::RegexMotif;
use crate::logging::init_logging;
use crate::mod_base_code::DnaBase;
use crate::monoid::Moniod;
use crate::reads_sampler::sampling_schedule::{
    IdxStats, ReferenceSequencesLookup,
};
use crate::threshold_mod_caller::MultipleThresholdModCaller;
use crate::thresholds::{
    get_modbase_probs_from_bam, log_calculated_thresholds,
    percentile_linear_interp,
};
use crate::util::{get_master_progress_bar, get_ticker};

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct MethylationEntropy {
    /// Input mod-BAM, may be repeated multiple times to calculate entropy
    /// across all input mod-BAMs.
    #[arg(short = 's', long = "in-bam", required = true)]
    in_bams: Vec<PathBuf>,
    /// Output BED file, if using `--region` this must be a directory.
    #[arg(short = 'o', long)]
    out_bed: Option<PathBuf>,
    /// Only used with `--regions`, prefix files in output directory with this
    /// string.
    #[arg(long, requires = "regions_fp")]
    prefix: Option<String>,
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
    /// Sample this many reads when estimating the filtering threshold. Reads
    /// will be sampled evenly across aligned genome. If a region is
    /// specified, either with the --region option or the --sample-region
    /// option, then reads will be sampled evenly across the region given.
    /// This option is useful for large BAM files. In practice, 10-50
    /// thousand reads is sufficient to estimate the model output
    /// distribution and determine the filtering threshold.
    #[arg(long, default_value_t = 10_042)]
    num_reads: usize,
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
    /// Respect soft masking in the reference FASTA.
    #[arg(long, requires = "reference_fasta", default_value_t = false)]
    mask: bool,
    /// Motif to use for entropy calculation, default will be CpG.
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
    #[arg(long, default_value_t = false)]
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
        for bam_fp in self.in_bams.iter() {
            IdxStats::check_any_mapped_reads(&bam_fp, None, None)
                .with_context(|| {
                    format!(
                        "did not find any mapped reads in {bam_fp:?}, perform \
                         alignment first"
                    )
                })?;
        }

        let mut writer: Box<dyn EntropyWriter> =
            match (self.out_bed.as_ref(), self.regions_fp.is_some()) {
                (Some(out_fp), false) => Box::new(
                    WindowsWriter::new_file(out_fp, self.header)
                        .context("failed to make writer to file")?,
                ),
                (Some(out_dir), true) => Box::new(
                    RegionsWriter::new(
                        out_dir,
                        self.prefix.as_ref(),
                        self.header,
                    )
                    .context(
                        "failed to make regions writer, output must be a \
                         directory",
                    )?,
                ),
                (None, false) => Box::new(
                    WindowsWriter::new_stdout(self.header)
                        .context("failed to make writer to stdout")?,
                ),
                (None, true) => {
                    bail!("must provide output directory with regions")
                }
            };

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()?;
        let multi_pb = MultiProgress::new();
        if self.suppress_progress {
            multi_pb.set_draw_target(indicatif::ProgressDrawTarget::hidden());
        }

        let (motif, combine_strands) =
            match (self.cpg, self.motif.as_ref(), self.base) {
                (true, _, _) => {
                    info!("using CpG motif and combining strands");
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

        let reference_sequence_lookup = ReferenceSequencesLookup::new(
            &self.in_bams,
            &self.reference_fasta,
            self.mask,
            &multi_pb,
        )?;
        let chrom_id_to_name =
            reference_sequence_lookup.get_chrom_id_to_name_lookup();

        let sliding_windows = pool.install(|| {
            if let Some(regions_fp) = self.regions_fp.as_ref() {
                SlidingWindows::new_with_regions(
                    reference_sequence_lookup,
                    regions_fp,
                    motif,
                    combine_strands,
                    self.num_positions,
                    window_size,
                    batch_size,
                )
            } else {
                SlidingWindows::new(
                    reference_sequence_lookup,
                    motif,
                    combine_strands,
                    self.num_positions,
                    window_size,
                    batch_size,
                )
            }
        })?;

        let threshold_caller =
            self.get_threshold_caller(&pool).map(|c| Arc::new(c))?;

        let (snd, rcv) = crossbeam::channel::bounded(10_000);

        let bam_fps = self.in_bams.clone();
        let min_coverage = self.min_valid_coverage;
        let threads = self.threads;
        let io_threads = self.io_threads.unwrap_or(threads);
        let max_filtered = self.max_filtered_positions.unwrap_or_else(|| {
            let max_filt_pos =
                (self.num_positions as f32 * 0.5f32).floor() as usize;
            info!("setting maximum filtered positions to {max_filt_pos}");
            max_filt_pos
        });

        let genome_prog = multi_pb
            .add(get_master_progress_bar(sliding_windows.total_length()));
        let rows_written = multi_pb.add(get_ticker());
        let skipped_windows = multi_pb.add(get_ticker());
        let windows_failed = multi_pb.add(get_ticker());
        let batches_failed = multi_pb.add(get_ticker());

        let what =
            if self.regions_fp.is_some() { "regions" } else { "windows" };

        genome_prog.set_message("genome positions processed");
        rows_written.set_message("rows written");
        windows_failed.set_message(format!("{what} failed"));
        skipped_windows.set_message(format!("{what} with zero coverage"));
        batches_failed.set_message("batches failed");

        pool.spawn(move || {
            for batch in sliding_windows {
                let n_pos = batch
                    .iter()
                    .map(|gw| {
                        let r = gw.get_range();
                        r.end - r.start
                    })
                    .sum::<u64>();
                let mut results = Vec::new();
                let (entropies, _) = rayon::join(
                    || {
                        let rs = batch
                            .into_par_iter()
                            .map(|window| {
                                process_entropy_window(
                                    window,
                                    min_coverage,
                                    max_filtered,
                                    io_threads,
                                    threshold_caller.clone(),
                                    &bam_fps,
                                )
                            })
                            .collect::<Vec<_>>();
                        genome_prog.inc(n_pos);
                        rs
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
                let num_reads = self.num_reads / self.in_bams.len();
                let mut agg = HashMap::new();
                for in_bam in self.in_bams.iter() {
                    let per_base_thresholds = get_modbase_probs_from_bam(
                        in_bam,
                        self.threads,
                        1_000_000,
                        None,
                        Some(num_reads),
                        None,
                        None,
                        None,
                        None,
                        None,
                        true,
                        self.suppress_progress,
                    )?;
                    agg.op_mut(per_base_thresholds);
                }
                let per_base_thresholds = agg
                    .iter_mut()
                    .map(|(dna_base, mod_base_probs)| {
                        mod_base_probs
                            .par_sort_by(|x, y| x.partial_cmp(y).unwrap());
                        let threshold = percentile_linear_interp(
                            &mod_base_probs,
                            self.filter_percentile,
                        )?;
                        Ok((*dna_base, threshold))
                    })
                    .collect::<anyhow::Result<HashMap<DnaBase, f32>>>()?;
                log_calculated_thresholds(&per_base_thresholds);
                Ok(MultipleThresholdModCaller::new(
                    per_base_thresholds,
                    per_mod_thresholds.unwrap_or(HashMap::new()),
                    0f32,
                ))
            })
        }
    }
}
