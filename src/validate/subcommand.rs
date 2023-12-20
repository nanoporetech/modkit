use anyhow::bail;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use clap::Args;
use log::info;
use rust_htslib::bam::{self, Read};

use crate::logging::init_logging;
use crate::mod_base_code::ModCodeRepr;
use crate::util::{get_targets, get_ticker, reader_is_bam, Region};

#[derive(Args)]
pub struct ValidateFromModbam {
    // running args
    // convert to list of bam bed inputs
    /// Argument accepts 2 values. The first value is the BAM file path with
    /// modified base tags. The second is a bed file with ground truth
    /// reference positions. The name field in the ground truth bed file
    /// should be the short name (single letter code or ChEBI ID) for a
    /// modified base or the corresponding canonical base. This argument
    /// can be provided more than once for multiple samples.
    #[arg(
	long,
	action = clap::ArgAction::Append,
	num_args = 2,
	value_names = ["BAM", "BED"]
    )]
    bam_and_bed: Vec<PathBuf>,
    /// Hide the progress bar.
    #[arg(long, default_value_t = false, hide_short_help = true)]
    suppress_progress: bool,
    /// Specify a file for debug logs to be written to, otherwise ignore them.
    /// Setting a file is recommended. (alias: log)
    #[arg(long, alias = "log")]
    log_filepath: Option<PathBuf>,
}

impl ValidateFromModbam {
    pub fn run(&self) -> anyhow::Result<()> {
        let _handle = init_logging(self.log_filepath.as_ref());
        for bam_and_bed in self.bam_and_bed.chunks(2) {
            let bam = &bam_and_bed[0];
            let bed = &bam_and_bed[1];

            let header = bam::IndexedReader::from_path(&bam).map(|reader| {
                if !reader_is_bam(&reader) {
                    info!(
                        "detected non-BAM input format, please consider using \
                         BAM, CRAM may be unstable"
                    );
                }
                reader.header().to_owned()
            })?;
            let tids = get_targets(&header, Option::<&Region>::None);
            let chrom_to_tid = tids
                .iter()
                .map(|reference_record| {
                    (reference_record.name.as_str(), reference_record.tid)
                })
                .collect::<HashMap<&str, u32>>();
            let _mod_positions = Self::parse_mods_from_bed(
                bed,
                &chrom_to_tid,
                self.suppress_progress,
            );
        }
        Ok(())
    }

    pub fn parse_mods_from_bed(
        bed_fp: &PathBuf,
        chrom_to_target_id: &HashMap<&str, u32>,
        suppress_pb: bool,
    ) -> anyhow::Result<Vec<(u32, bool, u64, u64, ModCodeRepr)>> {
        info!("parsing BED at {}", bed_fp.to_str().unwrap_or("invalid-UTF-8"));

        let fh = File::open(bed_fp)?;
        let mut mod_positions = Vec::new();
        let lines_processed = get_ticker();
        if suppress_pb {
            lines_processed
                .set_draw_target(indicatif::ProgressDrawTarget::hidden());
        }
        lines_processed.set_message("rows processed");
        let mut warned = HashSet::new();

        let reader = BufReader::new(fh);
        for line in
            reader.lines().filter_map(|l| l.ok()).filter(|l| !l.is_empty())
        {
            let parts = line.split_ascii_whitespace().collect::<Vec<&str>>();
            let chrom_name = parts[0];
            if warned.contains(chrom_name) {
                continue;
            }
            if parts.len() < 6 {
                info!("improperly formatted BED line {line}");
                continue;
            }
            let raw_start = &parts[1].parse::<u64>();
            let raw_end = &parts[2].parse::<u64>();
            let (start, end) = match (raw_start, raw_end) {
                (Ok(start), Ok(end)) => (*start, *end),
                _ => {
                    info!("improperly formatted BED line {line}");
                    continue;
                }
            };
            let mod_code = ModCodeRepr::parse(parts[3])?;
            let (fwd_strand, rev_strand) = match parts[5] {
                "+" => (true, false),
                "-" => (false, true),
                "." => (true, true),
                _ => {
                    info!("improperly formatted strand field {}", &parts[5]);
                    continue;
                }
            };
            if let Some(chrom_id) = chrom_to_target_id.get(chrom_name) {
                if fwd_strand {
                    mod_positions.push((*chrom_id, true, start, end, mod_code))
                }
                if rev_strand {
                    mod_positions.push((*chrom_id, false, start, end, mod_code))
                }
            } else {
                info!("skipping chrom {chrom_name}, not present in BAM header");
                warned.insert(chrom_name.to_owned());
                continue;
            }
            lines_processed.inc(1);
        }
        if mod_positions.is_empty() {
            bail!("zero valid positions parsed from BED file")
        }
        lines_processed.finish_and_clear();
        info!("processed {} BED lines", lines_processed.position());

        Ok(mod_positions)
    }
}
