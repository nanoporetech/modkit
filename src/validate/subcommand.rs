use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::string::FromUtf8Error;

use anyhow::{anyhow, bail};
use clap::Args;
use log::info;
use rust_htslib::bam::{Read, Reader, Record};

use crate::command_utils::parse_edge_filter_input;
use crate::extract::writer::PositionModCalls;
use crate::logging::init_logging;
use crate::mod_bam::BaseModCall;
use crate::mod_bam::{CollapseMethod, EdgeFilter, ModBaseInfo};
use crate::mod_base_code::ModCodeRepr;
use crate::read_ids_to_base_mod_probs::ReadBaseModProfile;
use crate::thresholds::Percentiles;
use crate::util::{
    get_reference_mod_strand, get_ticker, record_is_secondary, Strand,
};

pub struct GroundTruthSite {
    pub chrom: String,
    pub strand: Strand,
    pub mod_code: ModCodeRepr,
    pub positions: Vec<i64>,
}

fn parse_ground_truth_bed_line(line: &str) -> anyhow::Result<GroundTruthSite> {
    // todo convert this function to use nom
    let fields: Vec<&str> = line.split_whitespace().collect();
    if fields.len() < 6 {
        bail!("Invalid number of fields in BED line");
    }
    let chrom = fields[0].to_string();
    let start: i64 =
        fields[1].parse().map_err(|e| anyhow!("Error parsing start: {}", e))?;
    let end: i64 =
        fields[2].parse().map_err(|e| anyhow!("Error parsing end: {}", e))?;
    let raw_mod_code = fields[3];
    let strand_char = fields[5]
        .chars()
        .next()
        .ok_or_else(|| anyhow!("Error parsing strand {}", fields[5]))?;

    let strand = Strand::parse_char(strand_char)?;
    let mod_code = ModCodeRepr::parse(&raw_mod_code)
        .map_err(|e| anyhow!("Error parsing modified base code: {}", e))?;
    let positions = (start..end + 1).collect();

    Ok(GroundTruthSite { chrom, strand, mod_code, positions })
}

type ChromToTid = HashMap<String, u32>;

type TidStrandNamePositions =
    HashMap<u32, HashMap<Strand, HashMap<ModCodeRepr, Vec<i64>>>>;

fn get_chrom_to_tid(reader: &Reader) -> anyhow::Result<ChromToTid> {
    let header = reader.header().to_owned();
    let chrom_to_tid_result = (0..header.target_count())
        .map(|tid: u32| {
            let raw_name = header.tid2name(tid);
            let stringed_name: Result<String, FromUtf8Error> =
                String::from_utf8(raw_name.to_vec());
            stringed_name.map(|name: String| (name, tid))
        })
        .collect::<Result<HashMap<String, u32>, _>>();
    let chrom_to_tid = chrom_to_tid_result?;
    Ok(chrom_to_tid)
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

    let reader = BufReader::new(File::open(file_path)?);
    for ground_truth_site in reader
        .lines()
        .filter_map(|r| r.ok())
        .filter_map(|line| {
            Some(parse_ground_truth_bed_line(&line).map(Some).transpose())
                .flatten()
        })
        .filter_map(|r| r.ok())
    {
        let Some(tid) = chrom_to_tid.get(&ground_truth_site.chrom) else {
            continue;
        };
        result
            .entry(*tid)
            .or_insert_with(HashMap::new)
            .entry(ground_truth_site.strand)
            .or_insert_with(HashMap::new)
            .entry(ground_truth_site.mod_code)
            .or_insert_with(Vec::new)
            .extend(ground_truth_site.positions);
        lines_processed.inc(1);
    }
    if result.is_empty() {
        bail!("zero valid positions parsed from BED file".to_string());
    }
    lines_processed.finish_and_clear();
    info!("processed {} BED lines", lines_processed.position());
    Ok(result)
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
        bail!("Unmapped or secondary");
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
    let position_calls = PositionModCalls::from_profile(&mbp);
    for mod_call in position_calls {
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
                let (mod_call_prob, mod_call_code) = match mod_call
                    .base_mod_probs
                    .argmax_base_mod_call()
                {
                    BaseModCall::Canonical(p) => (
                        p,
                        ModCodeRepr::parse(
                            &mod_call.canonical_base.char().to_string(),
                        )?,
                    ),
                    BaseModCall::Modified(p, code) => (p, code),
                    BaseModCall::Filtered => {
                        unreachable!("argmax should not output filtered calls")
                    }
                };
                result
                    .entry((*mod_code, mod_call_code))
                    .or_insert_with(Vec::new)
                    .push(mod_call_prob);
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
) -> anyhow::Result<ModBaseQuals> {
    let lines_processed = get_ticker();
    if suppress_pb {
        lines_processed
            .set_draw_target(indicatif::ProgressDrawTarget::hidden());
    }
    lines_processed.set_message("BAM records processed");
    let (gt_mod_quals, errs) = reader.records().fold(
        (HashMap::new(), HashMap::new()),
        |(mut gt_mod_quals, mut errs), record| {
            match record {
                Ok(record) => {
                    match process_bam_record(
                        &record,
                        &mod_positions,
                        collapse_method,
                        edge_filter,
                    ) {
                        Ok(mod_base_quals) => {
                            for ((gt_code, call_code), mod_quals) in
                                mod_base_quals.into_iter()
                            {
                                gt_mod_quals
                                    .entry((gt_code.clone(), call_code.clone()))
                                    .or_insert_with(Vec::new)
                                    .extend(mod_quals);
                            }
                            lines_processed.inc(1)
                        }
                        Err(err) => {
                            let err_counter =
                                errs.entry(err.to_string()).or_insert(0);
                            *err_counter += 1;
                        }
                    }
                }
                Err(err) => {
                    let err_counter = errs.entry(err.to_string()).or_insert(0);
                    *err_counter += 1;
                }
            }
            (gt_mod_quals, errs)
        },
    );
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
    Ok(gt_mod_quals)
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

    // todo make argument groups
    // args for BAM record manipulation
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

    // threshold args
    // todo add direct thresholds support
    /// Filter out modified base calls where the probability of the predicted
    /// variant is below this confidence percentile. For example, 0.1 will
    /// filter out the 10% lowest confidence modification calls.
    #[arg(short = 'q', long, default_value_t = 0.1)]
    filter_quantile: f32,

    // misc args
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

        let mut all_gt_mod_quals = HashMap::new();
        for bam_and_bed in self.bam_and_bed.chunks(2) {
            let bam_path = &bam_and_bed[0];
            let bed_path = &bam_and_bed[1];

            let mut reader = Reader::from_path(bam_path.as_path())?;
            reader.set_threads(self.threads)?;
            let chrom_to_tid = get_chrom_to_tid(&reader)?;
            let mod_positions = parse_ground_truth_bed_file(
                bed_path,
                chrom_to_tid,
                self.suppress_progress,
            )?;
            info!(
                "parsing BAM at {}",
                bam_path.to_str().unwrap_or("invalid-UTF-8")
            );
            let gt_mod_quals = process_bam_file(
                &mut reader,
                mod_positions,
                collapse_method.as_ref(),
                edge_filter.as_ref(),
                self.suppress_progress,
            )?;
            for ((gt_code, call_code), mod_quals) in gt_mod_quals.iter() {
                all_gt_mod_quals
                    .entry((gt_code.clone(), call_code.clone()))
                    .or_insert_with(Vec::new)
                    .extend(mod_quals.clone());
            }
        }

        let mut all_quals = Vec::<f32>::new();
        for ((gt_code, call_code), mod_quals) in all_gt_mod_quals.iter() {
            all_quals.extend(mod_quals);
            //println!("gt:{} call:{}", gt_code, call_code);
            //println!("\t{:?}", mod_quals);
        }
        let thresh = Percentiles::new(&mut all_quals, &[self.filter_quantile])?;
        info!("{}", thresh.report());

        Ok(())
    }
}
