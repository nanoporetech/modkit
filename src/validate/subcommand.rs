use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use anyhow::{anyhow, Result};
use clap::Args;
use log::info;
use rust_htslib::bam::{Read, Reader, Record};

use crate::command_utils::parse_edge_filter_input;
use crate::logging::init_logging;
use crate::mod_bam::{CollapseMethod, EdgeFilter, ModBaseInfo};
use crate::mod_base_code::ModCodeRepr;
use crate::read_ids_to_base_mod_probs::ReadBaseModProfile;
use crate::util::{
    get_reference_mod_strand, get_ticker, record_is_secondary, Strand,
};

fn parse_ground_truth_bed_line(
    line: &str,
) -> Result<(String, Strand, ModCodeRepr, Vec<i64>), String> {
    let fields: Vec<&str> = line.split_whitespace().collect();
    if fields.len() >= 6 {
        let chromosome = fields[0].to_string();
        let strand_char = fields[5]
            .chars()
            .next()
            .ok_or_else(|| format!("Error parsing strand {}", fields[5]))?;
        let strand = Strand::parse_char(strand_char)?;
        let mod_code = fields[3].to_string();
        let start: i64 = fields[1]
            .parse()
            .map_err(|e| format!("Error parsing start: {}", e))?;
        let end: i64 = fields[2]
            .parse()
            .map_err(|e| format!("Error parsing end: {}", e))?;
        let positions = (start..end + 1).collect();
        Ok((
            chromosome,
            strand,
            ModCodeRepr::parse(&mod_code).expect(concat!(
                "Bad modified base code: ",
                stringify!(mod_code)
            )),
            positions,
        ))
    } else {
        Err("Invalid number of fields in BED line".to_string())
    }
}

type ChromToTid = HashMap<String, u32>;

type TidStrandNamePositions =
    HashMap<u32, HashMap<Strand, HashMap<ModCodeRepr, Vec<i64>>>>;

fn get_chrom_to_tid(reader: &Reader) -> ChromToTid {
    let header = reader.header().to_owned();
    let chrom_to_tid = (0..header.target_count())
        .map(|tid| {
            (String::from_utf8(header.tid2name(tid).to_vec()).unwrap(), tid)
        })
        .collect::<ChromToTid>();
    chrom_to_tid
}

fn parse_ground_truth_bed_file(
    file_path: &PathBuf,
    chrom_to_tid: ChromToTid,
    suppress_pb: bool,
) -> anyhow::Result<TidStrandNamePositions> {
    info!("parsing BED at {}", file_path.to_str().unwrap_or("invalid-UTF-8"));
    let mut result = HashMap::new();
    let lines_processed = get_ticker();
    if suppress_pb {
        lines_processed
            .set_draw_target(indicatif::ProgressDrawTarget::hidden());
    }
    lines_processed.set_message("rows processed");

    if let Ok(file) = File::open(file_path) {
        let reader = BufReader::new(file);
        for line in reader.lines() {
            if let Ok((chrom_name, strand, mod_code, positions)) =
                parse_ground_truth_bed_line(&line.as_ref().unwrap())
            {
                let Some(tid) = chrom_to_tid.get(&chrom_name) else {
                    continue;
                };
                result
                    .entry(*tid)
                    .or_insert_with(HashMap::new)
                    .entry(strand)
                    .or_insert_with(HashMap::new)
                    .entry(mod_code)
                    .or_insert_with(Vec::new)
                    .extend(positions);
                lines_processed.inc(1);
            } else {
                return Err(anyhow!(
                    "Error parsing BED line: {:?}",
                    line.unwrap()
                ));
            }
        }
        if result.is_empty() {
            return Err(anyhow!(
                "zero valid positions parsed from BED file".to_string()
            ));
        }
        lines_processed.finish_and_clear();
        info!("processed {} BED lines", lines_processed.position());
        Ok(result)
    } else {
        Err(anyhow!("Error opening BED file: {:?}", file_path))
    }
}

// ground truth and observed mod base codes pointing to vector of mod qualities
type ModBaseQuals = HashMap<(ModCodeRepr, ModCodeRepr), Vec<f32>>;

fn process_bam_record(
    record: &Record,
    mod_positions: &TidStrandNamePositions,
    collapse_method: Option<&CollapseMethod>,
    edge_filter: Option<&EdgeFilter>,
) -> anyhow::Result<ModBaseQuals> {
    let mut result = HashMap::new();
    let mbi = ModBaseInfo::new_from_record(&record)?;
    let record_name = String::from_utf8(record.qname().to_vec())
        .unwrap_or("utf-decode-failed".to_string());
    if record.is_unmapped() || record_is_secondary(&record) {
        return Err(anyhow!("Unmapped or secondary"));
    }

    let cgt_mod_pos =
        mod_positions.get(&(record.tid() as u32)).ok_or_else(|| {
            anyhow!(
                "No ground truth on this contig. {} not in {:?}",
                record_name,
                mod_positions.keys(),
            )
        })?;
    let mbp = ReadBaseModProfile::process_record(
        &record,
        &record_name,
        mbi,
        collapse_method,
        edge_filter,
        1,
    )?;
    for mod_call in mbp.profile {
        let Some(ref_pos) = mod_call.ref_position else {
            continue;
        };
        if ref_pos < 0 {
            continue;
        }
        let Some(alignment_strand) = mod_call.alignment_strand else {
            continue;
        };
        let ref_mod_strand =
            get_reference_mod_strand(mod_call.mod_strand, alignment_strand);
        let cgt_strand_mod_pos = cgt_mod_pos
            .get(&ref_mod_strand)
            .ok_or_else(|| anyhow!("No ground truth on this strand"))?;
        for (mod_code, ground_truth_pos) in cgt_strand_mod_pos.iter() {
            if ground_truth_pos.contains(&ref_pos) {
                result
                    .entry((*mod_code, mod_call.raw_mod_code))
                    .or_insert_with(Vec::new)
                    .push(mod_call.q_mod);
            }
        }
    }
    Ok(result)
}

fn process_bam_file(
    reader: &mut Reader,
    mod_positions: TidStrandNamePositions,
    collapse_method: Option<&CollapseMethod>,
    edge_filter: Option<&EdgeFilter>,
    suppress_pb: bool,
) {
    let lines_processed = get_ticker();
    if suppress_pb {
        lines_processed
            .set_draw_target(indicatif::ProgressDrawTarget::hidden());
    }
    lines_processed.set_message("BAM records processed");
    let mut errs: HashMap<String, i32> = HashMap::new();
    let mut gt_mod_quals = HashMap::new();
    for record in reader.records() {
        let record = match record {
            Ok(record) => record,
            Err(err) => {
                let err_counter = errs.entry(err.to_string()).or_insert(0);
                *err_counter += 1;
                continue;
            }
        };
        match process_bam_record(
            &record,
            &mod_positions,
            collapse_method,
            edge_filter,
        ) {
            Ok(mod_base_quals) => {
                for ((gt_code, call_code), mod_quals) in mod_base_quals.iter() {
                    gt_mod_quals
                        .entry((gt_code.clone(), call_code.clone()))
                        .or_insert_with(Vec::new)
                        .extend(mod_quals.clone());
                }
                lines_processed.inc(1);
            }
            Err(err) => {
                let err_counter = errs.entry(err.to_string()).or_insert(0);
                *err_counter += 1;
            }
        }
    }
    lines_processed.finish_and_clear();
    info!("processed {} BAM recrods", lines_processed.position());
    if !errs.is_empty() {
        let mut sorted_errs: Vec<_> = errs.iter().collect();
        sorted_errs.sort_by_key(|&(_, count)| std::cmp::Reverse(count));
        info!("Failed read reasons:");
        for (key, err) in sorted_errs {
            info!("\t{}: {}", key, err);
        }
    }
}

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
    /// Ignore a modified base class  _in_situ_ by redistributing base
    /// modification probability equally across other options. For example,
    /// if collapsing 'h', with 'm' and canonical options, half of the
    /// probability of 'h' will be added to both 'm' and 'C'. A full
    /// description of the methods can be found in collapse.md.
    #[arg(long, hide_short_help = true)]
    ignore: Option<String>,

    /// Discard base modification calls that are this many bases from the start
    /// or the end of the read. Two comma-separated values may be provided
    /// to asymmetrically filter out base modification calls from the start
    /// and end of the reads. For example, 4,8 will filter out base
    /// modification calls in the first 4 and last 8 bases of the read.
    #[arg(long)]
    edge_filter: Option<String>,
    /// Invert the edge filter, instead of filtering out base modification
    /// calls at the ends of reads, only _keep_ base modification calls at
    /// the ends of reads. E.g. if usually, "4,8" would remove (i.e. filter
    /// out) base modification calls in the first 4 and last 8 bases of the
    /// read, using this flag will keep only base modification calls in the
    /// first 4 and last 8 bases.
    #[arg(long, requires = "edge_filter", default_value_t = false)]
    invert_edge_filter: bool,

    /// Number of threads to use
    #[arg(short = 't', long, default_value_t = 4)]
    threads: usize,
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
        let collapse_method = match &self.ignore {
            Some(raw_mod_code) => {
                let mod_code = ModCodeRepr::parse(raw_mod_code)?;
                Some(CollapseMethod::ReDistribute(mod_code))
            }
            None => None,
        };
        let edge_filter = self
            .edge_filter
            .as_ref()
            .map(|raw| parse_edge_filter_input(raw, self.invert_edge_filter))
            .transpose()?;

        for bam_and_bed in self.bam_and_bed.chunks(2) {
            let bam_path = &bam_and_bed[0];
            let bed_path = &bam_and_bed[1];

            let mut reader = Reader::from_path(bam_path.as_path())?;
            reader.set_threads(self.threads)?;
            let chrom_to_tid = get_chrom_to_tid(&reader);
            let mod_positions = parse_ground_truth_bed_file(
                bed_path,
                chrom_to_tid,
                self.suppress_progress,
            )?;
            info!(
                "parsing BAM at {}",
                bam_path.to_str().unwrap_or("invalid-UTF-8")
            );
            let _gt_mod_quals = process_bam_file(
                &mut reader,
                mod_positions,
                collapse_method.as_ref(),
                edge_filter.as_ref(),
                self.suppress_progress,
            );
        }
        Ok(())
    }
}
