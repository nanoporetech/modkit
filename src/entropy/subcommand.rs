use std::collections::HashMap;
use std::path::PathBuf;

use anyhow::{bail, Context};
use clap::Args;
use indicatif::MultiProgress;
use log::{debug, error, info};
use rayon::prelude::*;
use rust_htslib::bam::{self, Read};
use rustc_hash::FxHashMap;

use crate::command_utils::{
    get_threshold_from_options, parse_per_mod_thresholds,
};
use crate::entropy::writers::{EntropyWriter, RegionsWriter, WindowsWriter};
use crate::entropy::{process_entropy_window, SlidingWindows};
use crate::logging::init_logging;
use crate::mod_base_code::DnaBase;
use crate::motif_bed::{get_masked_sequences, RegexMotif};
use crate::reads_sampler::sampling_schedule::IdxStats;
use crate::threshold_mod_caller::MultipleThresholdModCaller;
use crate::util::{get_targets, get_ticker};

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct MethylationEntropy {
    /// Input mod-BAM
    in_bam: PathBuf,
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

                let idx_stats =
                    IdxStats::new_from_path(&self.in_bam, None, None)?;

                SlidingWindows::new_with_regions(
                    &names_to_tid,
                    reference_sequences,
                    &idx_stats,
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
