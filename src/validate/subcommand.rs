use rustc_hash::FxHashMap;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fmt::{Debug, Display, Formatter};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;
use std::string::FromUtf8Error;

use ansi_term::Style;
use anyhow::{anyhow, bail};
use clap::Args;
use itertools::Itertools;
use lazy_static::lazy_static;
use log::info;
use log_once::warn_once;
use ndarray::Array1;
use prettytable::format::{
    FormatBuilder, LinePosition, LineSeparator, TableFormat,
};
use prettytable::{cell, row, Row, Table};
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{Read, Reader, Record};
use std::cmp::Ordering;

use crate::command_utils::parse_edge_filter_input;
use crate::logging::init_logging;
use crate::mod_bam::BaseModCall;
use crate::mod_bam::{CollapseMethod, EdgeFilter, ModBaseInfo};
use crate::mod_base_code::{
    DnaBase, ModCodeRepr, ANY_MOD_CODES, MOD_CODE_TO_DNA_BASE,
};
use crate::read_ids_to_base_mod_probs::{PositionModCalls, ReadBaseModProfile};
use crate::thresholds::percentile_linear_interp;
use crate::util::{
    format_int_with_commas, get_reference_mod_strand, get_ticker,
    record_is_not_primary, Strand,
};

/// todo investigate using this type in BaseModCall
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BaseStatus {
    Canonical,
    Modified(ModCodeRepr),
    NoCall,
    Mismatch(DnaBase),
    Deletion,
}

impl Display for BaseStatus {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match *self {
            BaseStatus::Canonical => write!(f, "-"),
            BaseStatus::Modified(b) => write!(f, "{}", b),
            BaseStatus::NoCall => write!(f, "No Call"),
            BaseStatus::Mismatch(b) => write!(f, "{}", b.char()),
            BaseStatus::Deletion => write!(f, "Deletion"),
        }
    }
}

impl BaseStatus {
    pub fn parse(raw: &str) -> anyhow::Result<Self> {
        if let Ok(code) = raw.parse::<char>() {
            if code == '-' {
                Ok(Self::Canonical)
            } else {
                Ok(Self::Modified(ModCodeRepr::Code(code)))
            }
        } else {
            if let Ok(chebi) = raw.parse::<u32>() {
                Ok(Self::Modified(ModCodeRepr::ChEbi(chebi)))
            } else {
                Err(anyhow!("failed to parse mod code {raw}"))
            }
        }
    }
}

lazy_static! {
    static ref TBL_FMT: TableFormat = {
        FormatBuilder::new()
            .column_separator('│')
            .borders('│')
            .separator(
                LinePosition::Title,
                LineSeparator::new('─', '┼', '├', '┤'),
            )
            .separator(
                LinePosition::Bottom,
                LineSeparator::new('─', '┴', '└', '┘'),
            )
            .separator(
                LinePosition::Top,
                LineSeparator::new('─', '┬', '┌', '┐'),
            )
            .padding(1, 1)
            .build()
    };
}

pub struct GroundTruthSite {
    pub chrom: String,
    pub strand: Strand,
    pub base_status: BaseStatus,
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
    let base_status = BaseStatus::parse(&raw_mod_code)
        .map_err(|e| anyhow!("Error parsing base status code: {}", e))?;
    if let BaseStatus::Modified(mod_code) = base_status {
        if ANY_MOD_CODES.contains(&mod_code) {
            warn_once!(
                "Ground truth file contains `any mod` code '{}'. If the \
                 intent was to represent a canonical position replace with a \
                 `-`",
                base_status
            )
        }
    }
    let positions = (start..end).collect();

    Ok(GroundTruthSite { chrom, strand, base_status, positions })
}

type TidToChrom = HashMap<u32, String>;

type ChromStrandPositionNames =
    HashMap<String, HashMap<Strand, BTreeMap<i64, BaseStatus>>>;

fn get_tid_to_chrom(reader: &Reader) -> anyhow::Result<TidToChrom> {
    let header = reader.header().to_owned();
    let tid_to_chrom_result = (0..header.target_count())
        .map(|tid: u32| {
            let raw_name = header.tid2name(tid);
            let stringed_name: Result<String, FromUtf8Error> =
                String::from_utf8(raw_name.to_vec());
            stringed_name.map(|name: String| (tid, name))
        })
        .collect::<Result<HashMap<u32, String>, _>>();
    let tid_to_chrom = tid_to_chrom_result?;
    Ok(tid_to_chrom)
}

fn parse_ground_truth_bed_file(
    file_path: &PathBuf,
    suppress_pb: bool,
) -> anyhow::Result<ChromStrandPositionNames> {
    info!("Parsing BED at {}", file_path.to_str().unwrap_or("invalid-UTF-8"));
    let mut result = HashMap::new();
    let lines_processed = get_ticker();
    if suppress_pb {
        lines_processed
            .set_draw_target(indicatif::ProgressDrawTarget::hidden());
    }
    lines_processed.set_message("rows processed");

    let reader = BufReader::new(File::open(file_path)?);
    for ground_truth_site in reader.lines().filter_map(|r| {
        r.map_err(|e| anyhow!("failed to read, {}", e.to_string()))
            .and_then(|line| parse_ground_truth_bed_line(&line))
            .ok()
    }) {
        let cs_res = result
            .entry(ground_truth_site.chrom)
            .or_insert_with(HashMap::new)
            .entry(ground_truth_site.strand)
            .or_insert_with(BTreeMap::new);
        for pos in ground_truth_site.positions {
            cs_res.insert(pos, ground_truth_site.base_status);
        }
        lines_processed.inc(1);
    }
    if result.is_empty() {
        bail!("zero valid positions parsed from BED file".to_string());
    }
    lines_processed.finish_and_clear();
    info!("Processed {} BED lines", lines_processed.position());
    Ok(result)
}

fn derive_canonical_base(
    gt_positions: &Vec<ChromStrandPositionNames>,
    can_base: Option<DnaBase>,
) -> anyhow::Result<DnaBase> {
    let mut base_statuses = HashSet::new();
    for chrom_strand_position_names in gt_positions {
        for (_, strand_map) in chrom_strand_position_names.iter() {
            for (_, base_status_map) in strand_map.iter() {
                for (_, base_status) in base_status_map.iter() {
                    base_statuses.insert(base_status);
                }
            }
        }
    }
    let mut can_base = can_base.clone();
    for base_status in base_statuses {
        match base_status {
            BaseStatus::Modified(mod_code) => {
                if let Some(existing_can_base) = &can_base {
                    if *existing_can_base
                        != *MOD_CODE_TO_DNA_BASE
                            .get(&mod_code)
                            .unwrap_or(existing_can_base)
                    {
                        return Err(anyhow::anyhow!(
                            "Multiple canonical bases represented in ground \
                             truth BED files: {} {}",
                            existing_can_base.char(),
                            MOD_CODE_TO_DNA_BASE
                                .get(&mod_code)
                                .unwrap_or(existing_can_base)
                                .char()
                        ));
                    }
                } else {
                    match MOD_CODE_TO_DNA_BASE.get(&mod_code) {
                        Some(extracted_can_base) => {
                            can_base = Some(*extracted_can_base);
                        }
                        None => {
                            continue;
                        }
                    }
                }
            }
            _ => {
                continue;
            }
        }
    }
    can_base.ok_or(anyhow::anyhow!(
        "Could not derive canonical base from ground truth."
    ))
}

// ground truth and observed base status pointing to vector of mod probabilities
type StatusProbs = HashMap<(BaseStatus, BaseStatus), Vec<f32>>;

fn process_bam_record(
    record: &Record,
    mod_positions: &ChromStrandPositionNames,
    tid_to_chrom: &TidToChrom,
    can_base: DnaBase,
    collapse_method: Option<&CollapseMethod>,
    edge_filter: Option<&EdgeFilter>,
) -> anyhow::Result<StatusProbs> {
    let mut result = HashMap::new();
    let mbi = ModBaseInfo::new_from_record(&record)?;
    let record_name = String::from_utf8(record.qname().to_vec())
        .unwrap_or("utf-decode-failed".to_string());
    if record.is_unmapped() || record_is_not_primary(&record) {
        bail!("Unmapped or secondary");
    }

    let chrom = tid_to_chrom
        .get(&(record.tid() as u32))
        .ok_or_else(|| anyhow!("Invalid record TID: {}", record.tid(),))?;
    let cgt_mod_pos = mod_positions
        .get(chrom)
        .ok_or_else(|| anyhow!("No ground truth on this contig.",))?;
    let mbp = ReadBaseModProfile::process_record(
        &record,
        &record_name,
        mbi,
        collapse_method,
        edge_filter,
        1,
    )?;
    let mut called_ref_pos = HashMap::new();
    for mod_call in PositionModCalls::from_profile(&mbp) {
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
        let Some(cs_mod_pos) = cgt_mod_pos.get(&ref_mod_strand) else {
            continue;
        };
        let strand_called_ref_pos =
            called_ref_pos.entry(ref_mod_strand).or_insert_with(HashSet::new);
        let Some(gt_code) = cs_mod_pos.get(&ref_pos) else {
            continue;
        };
        let (mod_call_prob, call_code) = match mod_call
            .base_mod_probs
            .argmax_base_mod_call()
        {
            BaseModCall::Canonical(p) => (p, BaseStatus::Canonical),
            BaseModCall::Modified(p, code) => (p, BaseStatus::Modified(code)),
            BaseModCall::Filtered => {
                unreachable!("argmax should not output filtered calls")
            }
        };
        result
            .entry((*gt_code, call_code))
            .or_insert_with(Vec::new)
            .push(mod_call_prob);
        strand_called_ref_pos.insert(ref_pos);
    }
    // add in no call, mismatch and deletion positions
    let r_st = record.reference_start();
    let r_en = record.reference_end();
    let q_seq = record.seq();
    let ref_to_query: FxHashMap<i64, i64> = record
        .aligned_pairs()
        .filter_map(|p| Some(p))
        .map(|pos| (pos[1], pos[0]))
        .collect();
    for (strand, positions) in called_ref_pos.iter() {
        let Some(cs_mod_pos) = cgt_mod_pos.get(&strand) else {
            continue;
        };
        for (pos, gt_code) in cs_mod_pos.range(r_st..r_en) {
            if positions.contains(pos) {
                // already recorded in result above
                continue;
            };
            let Some(q_pos) = ref_to_query.get(pos) else {
                result
                    .entry((*gt_code, BaseStatus::Deletion))
                    .or_insert_with(Vec::new)
                    .push(f32::NAN);
                continue;
            };
            let mut base = DnaBase::parse(q_seq[*q_pos as usize] as char)?;
            if record.is_reverse() {
                base = base.complement();
            }
            if base == can_base {
                result
                    .entry((*gt_code, BaseStatus::NoCall))
                    .or_insert_with(Vec::new)
                    .push(f32::NAN);
            } else {
                result
                    .entry((*gt_code, BaseStatus::Mismatch(base)))
                    .or_insert_with(Vec::new)
                    .push(f32::NAN);
            }
        }
    }

    Ok(result)
}

fn process_bam_file(
    reader: &mut Reader,
    mod_positions: &ChromStrandPositionNames,
    tid_to_chrom: &TidToChrom,
    can_base: DnaBase,
    collapse_method: Option<&CollapseMethod>,
    edge_filter: Option<&EdgeFilter>,
    suppress_pb: bool,
) -> anyhow::Result<StatusProbs> {
    let lines_processed = get_ticker();
    if suppress_pb {
        lines_processed
            .set_draw_target(indicatif::ProgressDrawTarget::hidden());
    }
    lines_processed.set_message("Records processed");
    let (status_probs, errs) = reader.records().fold(
        (HashMap::new(), HashMap::new()),
        |(mut status_probs, mut errs), record| {
            match record {
                Ok(record) => {
                    match process_bam_record(
                        &record,
                        &mod_positions,
                        tid_to_chrom,
                        can_base,
                        collapse_method,
                        edge_filter,
                    ) {
                        Ok(read_probs) => {
                            for ((gt_code, call_code), probs) in
                                read_probs.into_iter()
                            {
                                status_probs
                                    .entry((gt_code.clone(), call_code.clone()))
                                    .or_insert_with(Vec::new)
                                    .extend(probs);
                            }
                            lines_processed.inc(1)
                        }
                        Err(err) => {
                            *errs.entry(err.to_string()).or_insert(0) += 1;
                        }
                    }
                }
                Err(err) => {
                    let err_counter = errs.entry(err.to_string()).or_insert(0);
                    *err_counter += 1;
                }
            }
            (status_probs, errs)
        },
    );
    lines_processed.finish_and_clear();
    info!("Processed {} mapping records", lines_processed.position());
    if !errs.is_empty() {
        let mut sorted_errs: Vec<_> = errs.iter().collect();
        sorted_errs.sort_by_key(|&(_, count)| std::cmp::Reverse(count));
        info!("Failed read reasons:");
        for (key, err) in sorted_errs {
            info!("\t{}: {}", key, err);
        }
    }
    Ok(status_probs)
}

fn balance_ground_truth(status_probs: &mut StatusProbs) -> anyhow::Result<()> {
    // Calculate the total number of elements in each row
    let gt_totals: HashMap<_, _> = status_probs
        .iter()
        .map(|((gt_mod, _), values)| (*gt_mod, values.len()))
        .fold(HashMap::new(), |mut acc, (gt_mod, len)| {
            *acc.entry(gt_mod).or_insert(0) += len;
            acc
        });
    // Determine the target size for each row
    let gt_target_size = *gt_totals
        .values()
        .min()
        .ok_or_else(|| anyhow!("No minimum value found"))?;
    let totals_and_limits = gt_totals
        .into_iter()
        .map(|(gt_mod_repr, gt_total)| {
            let elements_to_remove =
                gt_total.checked_sub(gt_target_size).unwrap_or(0);
            (gt_mod_repr, (gt_total, elements_to_remove))
        })
        .collect::<HashMap<BaseStatus, (usize, usize)>>();

    for ((gt_mod, _), probs) in status_probs.iter_mut() {
        if probs.len() <= gt_target_size {
            continue;
        }
        if let Some((gt_total, elements_to_remove)) =
            totals_and_limits.get(gt_mod)
        {
            let n_obs = probs.len();
            let ratio = n_obs as f32 / *gt_total as f32;
            let samp_target_size =
                n_obs - (ratio * (*elements_to_remove) as f32).round() as usize;
            let keepers =
                Array1::linspace(0.0, (n_obs - 1) as f64, samp_target_size + 2)
                    .into_iter()
                    .skip(1)
                    .take(samp_target_size)
                    .map(|x| x.round() as usize)
                    .filter_map(|idx| probs.get(idx).copied())
                    .collect::<Vec<f32>>();
            *probs = keepers;
        }
    }
    Ok(())
}

fn machine_parseable_table(status_probs: &StatusProbs) -> String {
    let mut gt_codes: Vec<_> =
        status_probs.keys().map(|&(k, _)| k).unique().collect();
    gt_codes.sort();
    let call_codes: Vec<_> = status_probs.keys().map(|&(_, k)| k).collect();

    let mut all_codes: Vec<_> =
        gt_codes.iter().chain(call_codes.iter()).unique().collect();
    all_codes.sort();

    let mut out_str = "[[\"ground_truth_label\",\"".to_string();
    out_str.push_str(&all_codes.iter().map(|x| x.to_string()).join("\",\""));
    out_str.push_str("\"]");
    for gt_code in &gt_codes {
        out_str.push_str(",[\"");
        out_str.push_str(&gt_code.to_string());
        out_str.push_str("\"");
        for &call_code in &all_codes {
            let vector_length = status_probs
                .get(&(*gt_code, *call_code))
                .map(|v| v.len())
                .unwrap_or(0);
            out_str.push_str(",");
            out_str.push_str(&vector_length.to_string());
        }
        out_str.push_str("]");
    }
    out_str.push_str("]");
    out_str
}

fn print_table(status_probs: &StatusProbs, show_percentages: bool) {
    let mut gt_codes: Vec<_> =
        status_probs.keys().map(|&(k, _)| k).unique().collect();
    gt_codes.sort();
    let call_codes: Vec<_> = status_probs.keys().map(|&(_, k)| k).collect();

    let mut all_codes: Vec<_> =
        gt_codes.iter().chain(call_codes.iter()).unique().collect();
    all_codes.sort();

    let mut gt_totals: HashMap<BaseStatus, usize> = HashMap::new();
    if show_percentages {
        gt_totals = status_probs
            .iter()
            .map(|((gt_mod, _), values)| (*gt_mod, values.len()))
            .fold(HashMap::new(), |mut acc, (gt_mod, len)| {
                *acc.entry(gt_mod).or_insert(0) += len;
                acc
            });
    }

    let mut count_tbl = Table::new();
    count_tbl.set_format(*TBL_FMT);

    // Create a header row
    let mut header = Row::empty();
    header.add_cell(cell!(""));
    for &call_code in &all_codes {
        header.add_cell(cell!(&call_code.to_string()));
    }
    count_tbl.set_titles(header);

    // Create table rows
    for gt_code in &gt_codes {
        let mut row = Row::empty();
        row.add_cell(cell!(&gt_code.to_string()));
        for &call_code in &all_codes {
            let vector_length = status_probs
                .get(&(*gt_code, *call_code))
                .map(|v| v.len())
                .unwrap_or(0);
            if show_percentages {
                let gt_total = gt_totals.get(gt_code).unwrap();
                row.add_cell(
                    cell!(&format!(
                        "{:.2}%",
                        100.0 * vector_length as f32 / *gt_total as f32
                    ))
                    .style_spec("r"),
                );
            } else {
                row.add_cell(
                    cell!(&format!(
                        "{}",
                        format_int_with_commas(vector_length as isize)
                    ))
                    .style_spec("r"),
                );
            }
        }
        count_tbl.add_row(row);
    }

    let mut longest_gt_code_len: usize = 0;
    if let Some(longest_gt_code) = gt_codes
        .iter()
        .max_by(|&a, &b| a.to_string().len().cmp(&b.to_string().len()))
    {
        longest_gt_code_len = longest_gt_code.to_string().len();
    }
    let mut metatable = Table::new();
    metatable.set_format(*prettytable::format::consts::FORMAT_CLEAN);
    metatable.add_row(row!(
        "",
        format!("{:1$}Called Base", "", longest_gt_code_len + 4)
    ));
    metatable.add_row(row!("\n\n\nGround\nTruth", count_tbl));

    // Print the table
    metatable.printstd();
}

#[derive(Args)]
#[command(arg_required_else_help = true)]
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
    /// Canonical base to evaluate. By default, this will be derived from mod
    /// codes in ground truth BED files. For ground truth with only canonical
    /// sites and/or ChEBI codes this values must be set.
    #[clap(short = 'c', long)]
    canonical_base: Option<DnaBase>,

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
    /// Specify a file for machine parseable output.
    #[arg(short = 'o', long, alias = "out")]
    out_filepath: Option<PathBuf>,
    /// Specify a file for debug logs to be written to, otherwise ignore them.
    /// Setting a file is recommended. (alias: log)
    #[arg(long, alias = "log")]
    log_filepath: Option<PathBuf>,
}

impl ValidateFromModbam {
    pub fn run(&self) -> anyhow::Result<()> {
        let _handle = init_logging(self.log_filepath.as_ref());
        let mut out_handle: Option<File> = None;
        if let Some(file_path) = self.out_filepath.clone() {
            out_handle = Some(File::create(&file_path)?);
        }
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

        // parse bed files and determine canonical base
        let mut bed_paths = Vec::new();
        let mut bam_path_to_bed_indices: HashMap<PathBuf, Vec<usize>> =
            HashMap::new();
        for bam_and_bed in self.bam_and_bed.chunks(2) {
            let bam_path = &bam_and_bed[0]
                .canonicalize()
                .map_err(|e| anyhow::anyhow!("Cannot resolve path: {}", e))?;
            let bed_path = &bam_and_bed[1]
                .canonicalize()
                .map_err(|e| anyhow::anyhow!("Cannot resolve path: {}", e))?;

            let bed_idx = if let Some(bed_idx) =
                bed_paths.iter().position(|s| s == bed_path)
            {
                bed_idx
            } else {
                bed_paths.push(bed_path.to_path_buf());
                bed_paths.len() - 1
            };
            bam_path_to_bed_indices
                .entry(bam_path.clone())
                .or_insert_with(Vec::new)
                .push(bed_idx);
        }

        let gt_positions: Vec<_> = bed_paths
            .iter()
            .map(|bed_path| {
                parse_ground_truth_bed_file(bed_path, self.suppress_progress)
            })
            .collect::<Result<Vec<_>, _>>()?;
        let can_base =
            derive_canonical_base(&gt_positions, self.canonical_base)?;
        info!("Canonical base: {}", can_base.char());

        let mut all_probs = HashMap::new();
        for (bam_path, bed_indices) in bam_path_to_bed_indices {
            let mut reader = Reader::from_path(bam_path.as_path())?;
            reader.set_threads(self.threads)?;
            let tid_to_chrom = get_tid_to_chrom(&reader)?;
            info!(
                "Parsing mapping at {}",
                bam_path.to_str().unwrap_or("invalid-UTF-8")
            );

            for bed_idx in bed_indices {
                let status_probs = process_bam_file(
                    &mut reader,
                    &gt_positions[bed_idx],
                    &tid_to_chrom,
                    can_base,
                    collapse_method.as_ref(),
                    edge_filter.as_ref(),
                    self.suppress_progress,
                )?;
                for ((gt_code, call_code), probs) in status_probs.into_iter() {
                    all_probs
                        .entry((gt_code, call_code))
                        .or_insert_with(Vec::new)
                        .extend(probs);
                }
            }
        }

        // sort prob vectors
        for ((_, _), probs) in all_probs.iter_mut() {
            probs.sort_by_key(|&x| x.to_bits());
        }
        info!("Raw counts summary");
        print_table(&all_probs, false);
        if let Some(valid_out_handle) = &mut out_handle {
            valid_out_handle
                .write_all(
                    &format!(
                        "full_contingency_table: {}\n",
                        machine_parseable_table(&all_probs)
                    )
                    .into_bytes(),
                )
                .map_err(|e| anyhow::anyhow!("Error writing to file: {}", e))?;
        }

        // filter to only modified base calls
        all_probs.retain(|&(_, call_code), _| {
            call_code == BaseStatus::Canonical
                || matches!(call_code, BaseStatus::Modified(_))
        });

        info!("Balancing ground truth call totals");
        balance_ground_truth(&mut all_probs)?;
        let total_calls =
            all_probs.iter().map(|(_, values)| values.len()).sum::<usize>();
        let correct_calls = all_probs
            .iter()
            .filter(|&((gt_code, call_code), _)| gt_code == call_code)
            .map(|(_, values)| values.len())
            .sum::<usize>();
        let raw_acc = 100.0 * correct_calls as f32 / total_calls as f32;
        info!("Raw accuracy: {:.2}%", raw_acc);
        info!("Raw modified base calls contingency table");
        print_table(&all_probs, true);
        if let Some(valid_out_handle) = &mut out_handle {
            valid_out_handle
                .write_all(&format!("raw_accuracy: {}\n", raw_acc).into_bytes())
                .map_err(|e| anyhow::anyhow!("Error writing to file: {}", e))?;
            valid_out_handle
                .write_all(
                    &format!(
                        "raw_contingency_table: {}\n",
                        machine_parseable_table(&all_probs)
                    )
                    .into_bytes(),
                )
                .map_err(|e| anyhow::anyhow!("Error writing to file: {}", e))?;
        }

        let mut flat_probs = Vec::<f32>::new();
        for (_, probs) in all_probs.iter() {
            flat_probs.extend(probs);
        }
        flat_probs.sort_by(|x, y| x.partial_cmp(y).unwrap_or(Ordering::Equal));
        if flat_probs.iter().any(|v| v.is_nan()) {
            bail!("Failed to compare values");
        }
        let thresh =
            percentile_linear_interp(&flat_probs, self.filter_quantile)?;
        info!("Call probability threshold: {:.4}", thresh);

        // apply threshold and print filtered table
        all_probs.values_mut().for_each(|probs| {
            probs.retain(|&p| p > thresh);
        });
        let filt_calls =
            all_probs.iter().map(|(_, values)| values.len()).sum::<usize>();
        let correct_filt_calls = all_probs
            .iter()
            .filter(|&((gt_code, call_code), _)| gt_code == call_code)
            .map(|(_, values)| values.len())
            .sum::<usize>();
        let filt_acc = 100.0 * correct_filt_calls as f32 / filt_calls as f32;
        info!(
            "{}",
            Style::new()
                .bold()
                .paint(format!("Filtered accuracy: {:.2}%", filt_acc)),
        );
        info!("Filtered modified base calls contingency table");
        print_table(&all_probs, true);
        if let Some(valid_out_handle) = &mut out_handle {
            valid_out_handle
                .write_all(
                    &format!("filtered_accuracy: {}\n", filt_acc).into_bytes(),
                )
                .map_err(|e| anyhow::anyhow!("Error writing to file: {}", e))?;
            valid_out_handle
                .write_all(
                    &format!(
                        "filtered_contingency_table: {}\n",
                        machine_parseable_table(&all_probs)
                    )
                    .into_bytes(),
                )
                .map_err(|e| anyhow::anyhow!("Error writing to file: {}", e))?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod validate_tests {
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    use std::path::PathBuf;
    use tempfile::NamedTempFile;

    use super::ValidateFromModbam;

    #[test]
    fn test_validate() {
        let mut bam_bed_opts = Vec::new();
        bam_bed_opts.push(PathBuf::from("tests/resources/input_5mC.bam"));
        bam_bed_opts.push(PathBuf::from(
            "tests/resources/CGI_ladder_3.6kb_ref_CG_5mC.bed",
        ));
        bam_bed_opts.push(PathBuf::from("tests/resources/input_C.bam"));
        bam_bed_opts.push(PathBuf::from(
            "tests/resources/CGI_ladder_3.6kb_ref_CG_C.bed",
        ));

        let temp_file =
            NamedTempFile::new().expect("Output file creation failed");
        let file_path_buf = Some(temp_file.path().to_path_buf());
        println!(
            "Temp output file location: {}",
            temp_file.path().to_str().unwrap()
        );
        let log_temp_file =
            NamedTempFile::new().expect("Output file creation failed");
        let log_file_path_buf = Some(log_temp_file.path().to_path_buf());

        let val_opts = ValidateFromModbam {
            bam_and_bed: bam_bed_opts,
            ignore: None,
            edge_filter: None,
            invert_edge_filter: false,
            canonical_base: None,
            filter_quantile: 0.1,
            threads: 4,
            suppress_progress: false,
            out_filepath: file_path_buf,
            log_filepath: log_file_path_buf,
        };
        val_opts.run().unwrap();
        let file_path = temp_file.path().to_path_buf();
        let file = File::open(&file_path).unwrap();
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let line_content = line.unwrap();
            if line_content.starts_with("raw_accuracy") {
                let accuracy: f32 = line_content
                    .split_whitespace()
                    .skip(1)
                    .next()
                    .and_then(|s| s.parse().ok())
                    .unwrap_or_default();
                assert_eq!(accuracy, 79.97674);
            }
            if line_content.starts_with("filtered_accuracy") {
                let accuracy: f32 = line_content
                    .split_whitespace()
                    .skip(1)
                    .next()
                    .and_then(|s| s.parse().ok())
                    .unwrap_or_default();
                assert_eq!(accuracy, 84.24458);
            }
        }
    }
}
