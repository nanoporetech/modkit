use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use anyhow::{anyhow, bail};
use clap::Args;
use indicatif::{ParallelProgressIterator, ProgressIterator};
use itertools::Itertools;
use log::{debug, error, info};
use log_once::info_once;
use rayon::prelude::*;
use rustc_hash::FxHashSet;

use crate::dmr::bedmethyl::BedMethylLine;
use crate::logging::init_logging;
use crate::mod_base_code::ModCodeRepr;
use crate::monoid::Moniod;
use crate::stats::MethylationStats;
use crate::tabix::HtsTabixHandler;
use crate::util::{get_subroutine_progress_bar, get_ticker, GenomeRegion};

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct EntryStats {
    /// Input bedMethyl table. Should be bgzip-compressed and have an
    /// associated Tabix index. The tabix index will be assumed to be
    /// $this_file.tbi
    in_bedmethyl: PathBuf,
    /// BED file of regions to aggregate base modification over.
    #[arg(long)]
    regions: PathBuf,
    /// Specify which base modification codes to use. Default will report
    /// information on all base modification codes encountered.
    #[arg(long, short = 'c', alias = "codes", value_delimiter = ',', action=clap::ArgAction::Append)]
    mod_codes: Option<Vec<String>>,
    /// Only use records with at least this much valid coverage.
    #[arg(short = 'm', long, alias = "min-cov", default_value_t = 1)]
    min_coverage: u64,
    /// Specify the output file to write the results table.
    #[arg(long, short = 'o', alias = "out")]
    out_table: PathBuf,
    /// Force overwrite the output file.
    #[arg(long, default_value_t = false)]
    force: bool,
    /// Don't add the header describing the columns to the output
    #[arg(long, default_value_t = false)]
    no_header: bool,
    /// Specify a file to write debug logs to.
    #[arg(long, alias = "log")]
    log_filepath: Option<PathBuf>,
    /// Number of threads to use
    #[arg(short = 't', long, default_value_t = 4)]
    threads: usize,
}

impl EntryStats {
    pub fn run(&self) -> anyhow::Result<()> {
        let _ = init_logging(self.log_filepath.as_ref());
        let index: HtsTabixHandler<BedMethylLine> =
            HtsTabixHandler::from_path(&self.in_bedmethyl)?;
        let fh = if self.force {
            File::create(&self.out_table)?
        } else {
            File::create_new(&self.out_table)?
        };
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()?;
        let mpb = indicatif::MultiProgress::new();
        let mut reader = BufReader::new(File::open(&self.regions)?)
            .lines()
            .skip_while(|r| {
                r.as_ref().map(|l| l.starts_with('#')).unwrap_or(true)
            })
            .peekable();
        let parser = match reader.peek() {
            Some(Ok(l)) => {
                let num_fields = l.split('\t').count();
                if num_fields <= 4 {
                    |l: &str| GenomeRegion::parse_unstranded_bed_line(l)
                } else {
                    |l: &str| GenomeRegion::parse_stranded_bed_line(l)
                }
            }
            Some(Err(e)) => bail!("failed to inspect regions BED, {e}"),
            None => bail!("failed to inspect regions BED, no valid lines"),
        };
        let parse_pb = mpb.add(get_ticker());
        parse_pb.set_message("parsing regions");
        let genome_regions = BufReader::new(File::open(&self.regions)?)
            .lines()
            .progress_with(parse_pb)
            .map(|r| {
                r.map_err(|e| anyhow!("failed to read from regions file, {e}"))
                    .and_then(|raw| parser(&raw))
            })
            .collect::<anyhow::Result<Vec<GenomeRegion>>>()?;

        if genome_regions.is_empty() {
            bail!("failed to load any regions")
        }

        let mod_codes = self
            .mod_codes
            .as_ref()
            .map(|codes| {
                codes
                    .iter()
                    .map(|raw| ModCodeRepr::parse(raw))
                    .collect::<anyhow::Result<FxHashSet<ModCodeRepr>>>()
            })
            .transpose()?;
        if let Some(codes) = mod_codes.as_ref() {
            let formatted = codes.iter().join(",");
            info!("narrowing calculation to modification codes: {formatted}");
        }

        let genome_regions = genome_regions
            .into_iter()
            .filter_map(|r| {
                if index.has_contig(&r.chrom) {
                    Some(r)
                } else {
                    info_once!(
                        "bedmethyl does not have contig {}, skipping.",
                        &r.chrom
                    );
                    None
                }
            })
            .collect::<Vec<GenomeRegion>>();

        debug!("calculating stats on {} regions", genome_regions.len());

        let stats_pb =
            mpb.add(get_subroutine_progress_bar(genome_regions.len()));
        stats_pb.set_message("calculating stats");
        let (stats, obs_codes) = pool.install(|| {
            genome_regions
                .into_par_iter()
                .progress_with(stats_pb)
                .filter_map(|gr| {
                    match gr.into_stats(
                        &index,
                        self.min_coverage,
                        mod_codes.as_ref(),
                    ) {
                        Ok(stats) => Some(stats),
                        Err(e) => {
                            debug!("failed to get stats, {e}");
                            None
                        }
                    }
                })
                .fold(
                    || (Vec::new(), FxHashSet::default()),
                    |(mut agg, mut codes), next| {
                        if mod_codes.is_none() {
                            codes.extend(
                                next.per_mod_methylation.keys().copied(),
                            );
                        }
                        agg.push(next);
                        (agg, codes)
                    },
                )
                .reduce(
                    || (Vec::new(), FxHashSet::default()),
                    |(a, b), (c, d)| (a.op(c), b.op(d)),
                )
        });

        let mod_codes =
            if let Some(codes) = mod_codes { codes } else { obs_codes }
                .into_iter()
                .sorted()
                .collect::<Vec<ModCodeRepr>>();

        if mod_codes.is_empty() && !stats.is_empty() {
            error!("zero modification codes detected or provided");
        }

        let mut table = prettytable::Table::new();
        if !self.no_header {
            table.set_titles(MethylationStats::header(&mod_codes));
        }
        stats.into_iter().for_each(|x| {
            table.add_row(x.into_row(&mod_codes));
        });

        let csv_writer =
            csv::WriterBuilder::new().delimiter('\t' as u8).from_writer(fh);
        table.to_csv_writer(csv_writer)?;

        Ok(())
    }
}
