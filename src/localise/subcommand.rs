use std::collections::HashMap;
use std::fs::File;
use std::io::{stdout, BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::sync::Arc;

use anyhow::{anyhow, bail};
use clap::Args;
use indicatif::{MultiProgress, ParallelProgressIterator, ProgressIterator};
use itertools::Itertools;
use log::{debug, info};
use rayon::prelude::*;
use rustc_hash::FxHashMap;

use crate::dmr::bedmethyl::BedMethylLine;
use crate::localise::util::{LocalizedModCounts, StrandedFeatures};
use crate::logging::init_logging;
use crate::monoid::Moniod;
use crate::tabix::HtsTabixHandler;
use crate::util::{
    get_master_progress_bar, get_ticker, read_sequence_lengths_file,
    GenomeRegion, StrandRule,
};

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct EntryLocalize {
    /// Input bedMethyl table. Should be bgzip-compressed and have an
    /// associated Tabix index. The tabix index will be assumed to be
    /// $this_file.tbi
    in_bedmethyl: PathBuf,
    /// BED file of regions to calculate enrichment around. These BED records
    /// serve as the points from which the `--window` number of bases is
    /// centered.
    #[arg(long)]
    regions: PathBuf,
    /// Create plots showing %-modification vs. offset. Argument should be a
    /// path to a file.
    #[clap(help_heading = "Output Options")]
    #[arg(long = "chart")]
    chart_filepath: Option<PathBuf>,
    /// Give the HTML document and chart a name.
    #[clap(help_heading = "Output Options")]
    #[arg(long = "name", requires = "chart_filepath")]
    chart_name: Option<String>,
    /// Number of base pairs to search around, for example if your BED region
    /// records are single positions, a window of 500 will look 500 base
    /// pairs upstream and downstream of that position. If your region BED
    /// records are larger regions, this will expand from the midpoint of
    /// that region.
    #[arg(short = 'w', long = "window", default_value_t = 2000)]
    expand_window: u64,
    // todo
    // /// Expand the aggregation window `expand-window` base pairs from the
    // start /// and end of each region instead of the midpoint.
    // #[arg(long, alias = "widen")]
    // expand_from_edge: bool,
    // /// Only aggregate for a single modification code, default is to
    // aggregate /// for all modification codes at once.
    // #[arg(long)]
    // mod_code: Option<String>,
    /// Whether to only keep bedMethyl records on the "same" strand or
    /// "opposite" strand.
    #[arg(short = 's', long)]
    stranded: Option<StrandedFeatures>,
    /// Force use bedMethyl records from a particular strand, default is to use
    /// the strand as given in the BED file (will use BOTH for BED3).
    #[arg(long)]
    stranded_features: Option<StrandRule>,
    /// Minimum valid coverage to use a bedMethyl record
    #[arg(long, default_value_t = 3)]
    min_coverage: u64,
    /// TSV of genome sizes, should be <chrom>\t<size_in_bp>
    #[arg(long, short = 'r')]
    genome_sizes: PathBuf,
    /// Optionally specify a file to write output to, default is stdout.
    #[arg(long, short = 'o')]
    out_file: Option<PathBuf>,
    /// Specify a file to write debug logs to.
    #[clap(help_heading = "Logging Options")]
    #[arg(long, alias = "log")]
    log_filepath: Option<PathBuf>,
    /// Number of threads to use.
    #[clap(help_heading = "Compute Options")]
    #[arg(short = 't', long, default_value_t = 4)]
    threads: usize,
    /// Number of tabix/bgzf IO threads to use.
    #[clap(help_heading = "Compute Options")]
    #[arg(long, default_value_t = 2)]
    io_threads: usize,
    /// Force overwrite of existing output file.
    #[clap(help_heading = "Output Options")]
    #[arg(long, short = 'f', default_value_t = false)]
    force: bool,
    #[clap(help_heading = "Compute Options")]
    #[arg(
        long = "batch-size",
        hide_short_help = true,
        default_value_t = 500_000
    )]
    batch_size_bp: u64,
}

impl EntryLocalize {
    fn load_focus_regions(
        &self,
        sequence_lengths: &FxHashMap<String, u64>,
        index: &HtsTabixHandler<BedMethylLine>,
        multi_progress: &MultiProgress,
    ) -> anyhow::Result<Vec<GenomeRegion>> {
        let pb = multi_progress.add(get_ticker());
        pb.set_message("regions parsed");
        let mut reader = BufReader::new(File::open(&self.regions)?)
            .lines()
            .skip_while(|r| {
                r.as_ref().map(|l| l.starts_with('#')).unwrap_or(true)
            })
            .peekable();
        let parser = match reader.peek() {
            Some(Ok(l)) => {
                let num_fields = l.split_whitespace().count();
                if num_fields <= 4 {
                    |l: &str| GenomeRegion::parse_unstranded_bed_line(l)
                } else {
                    |l: &str| GenomeRegion::parse_stranded_bed_line(l)
                }
            }
            Some(Err(e)) => bail!("failed to inspect regions BED, {e}"),
            None => bail!("failed to inspect regions BED, no valid lines"),
        };

        let (regions, errs) = BufReader::new(File::open(&self.regions)?)
            .lines()
            .progress_with(pb)
            .map(|r| {
                r.map_err(|e| {
                    anyhow!("failed to read from regions BED file, {e}")
                })
                .and_then(|l| parser(&l))
            })
            .fold((Vec::new(), HashMap::new()), |(mut acc, mut errs), next| {
                match next {
                    Ok(gr) => {
                        acc.push(gr);
                        (acc, errs)
                    }
                    Err(e) => {
                        *errs.entry(e.to_string()).or_insert(0u32) += 1u32;
                        (acc, errs)
                    }
                }
            });

        if regions.is_empty() {
            let mut err_msg = "failure reasons: ".to_string();
            for (err, count) in
                errs.into_iter().sorted_by(|(_, x), (_, y)| x.cmp(y))
            {
                err_msg.push_str(&format!("\t{err}: {count}\n"));
            }
            bail!("failed to load any regions, {err_msg}")
        }
        let (regions, missing_from_sizes, missing_from_index) =
            regions.into_iter().fold(
                (Vec::new(), 0usize, 0usize),
                |(mut acc, no_size, no_index), mut next| {
                    if !sequence_lengths.contains_key(&next.chrom) {
                        (acc, no_size + 1, no_index)
                    } else if !index.has_contig(&next.chrom) {
                        (acc, no_size, no_index + 1)
                    } else {
                        let mp = next.midpoint();
                        let start = mp
                            .checked_sub(self.expand_window + 1)
                            .unwrap_or(0u64);
                        // safe because of conditional above
                        let end = std::cmp::min(
                            mp.saturating_add(self.expand_window),
                            *sequence_lengths.get(&next.chrom).unwrap(),
                        );
                        next.start = start;
                        next.end = end;
                        acc.push(next);
                        (acc, no_size, no_index)
                    }
                },
            );

        if missing_from_sizes > 0usize {
            debug!(
                "{missing_from_sizes} contigs in regions that are missing \
                 from contig sizes file"
            );
        }
        if missing_from_index > 0usize {
            debug!(
                "{missing_from_index} contigs in regions that are missing \
                 from tabix header"
            );
        }
        if regions.is_empty() {
            bail!("failed to find any valid regions")
        } else {
            debug!("{} valid regions parsed", regions.len());
        }

        info!("loaded {} regions", regions.len());
        Ok(regions)
    }

    pub fn run(&self) -> anyhow::Result<()> {
        let _ = init_logging(self.log_filepath.as_ref());

        // log out some of the settings...
        let min_cov = self.min_coverage;
        debug!("requiring min valid coverage of {min_cov}");
        if let Some(feature_strand) = self.stranded_features.as_ref() {
            info!("only using features from {feature_strand}");
        }
        if let Some(stranded) = self.stranded {
            info!(
                "only using features that are on {stranded:?} strand as the \
                 region of interest"
            );
        }

        let multi_progress = MultiProgress::new();
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()?;
        let writer: Box<dyn Write> =
            if let Some(out_fp) = self.out_file.as_ref() {
                if self.force {
                    Box::new(BufWriter::new(File::create(out_fp)?))
                } else {
                    Box::new(BufWriter::new(File::create_new(out_fp)?))
                }
            } else {
                Box::new(BufWriter::new(stdout()))
            };

        let writer =
            csv::WriterBuilder::new().delimiter('\t' as u8).from_writer(writer);

        info!("loading sequence lengths from {:?}", &self.genome_sizes);

        let sequence_lengths = read_sequence_lengths_file(&self.genome_sizes)?
            .into_iter()
            .collect::<FxHashMap<String, u64>>();
        let tabix_index = HtsTabixHandler::from_path(&self.in_bedmethyl)
            .map(|x| Arc::new(x))?;
        let genome_regions = self.load_focus_regions(
            &sequence_lengths,
            &tabix_index,
            &multi_progress,
        )?;
        info!("loaded {} regions", genome_regions.len());

        let successes =
            multi_progress.add(get_master_progress_bar(genome_regions.len()));

        let stranded_features = self.stranded;
        let counts = pool.install(|| {
            genome_regions
                .into_par_iter()
                .progress_with(successes)
                .map(|gr| {
                    gr.into_localized_mod_counts(
                        &tabix_index,
                        self.stranded_features,
                        stranded_features,
                        self.io_threads,
                    )
                })
                .fold(
                    || LocalizedModCounts::zero(),
                    |counts, next| match next {
                        Ok(lc) => counts.op(lc),
                        Err(e) => {
                            debug!("region failed, {e}");
                            counts
                        }
                    },
                )
                .reduce(|| LocalizedModCounts::zero(), |a, b| a.op(b))
        });

        let table = counts.get_table(&multi_progress);
        table.to_csv_writer(writer)?;

        if let Some(p) = self.chart_filepath.as_ref() {
            let fh = if self.force {
                File::create(p)?
            } else {
                File::create_new(p)?
            };
            let mut writer = BufWriter::new(fh);
            let blob = counts.get_plot(self.chart_name.as_ref())?;
            writer.write(blob.as_bytes())?;
        }

        multi_progress.clear()?;

        Ok(())
    }
}
