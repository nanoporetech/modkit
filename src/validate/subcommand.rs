use std::collections::HashMap;
use std::fmt::{Debug, Display, Formatter};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::string::FromUtf8Error;

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
use rust_htslib::bam::{Read, Reader, Record};
use std::cmp::Ordering;

use crate::command_utils::parse_edge_filter_input;
use crate::extract::writer::PositionModCalls;
use crate::logging::init_logging;
use crate::mod_bam::BaseModCall;
use crate::mod_bam::{CollapseMethod, EdgeFilter, ModBaseInfo};
use crate::mod_base_code::{DnaBase, ModCodeRepr, ANY_MOD_CODES};
use crate::read_ids_to_base_mod_probs::ReadBaseModProfile;
use crate::thresholds::percentile_linear_interp;
use crate::util::{
    get_reference_mod_strand, get_ticker, record_is_secondary, Strand,
};

/// todo investigate using this type in BaseModCall
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum BaseStatus {
    Canonical,
    Modified(ModCodeRepr),
    Mismatch(DnaBase),
    Deletion,
}

impl Display for BaseStatus {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match *self {
            BaseStatus::Canonical => write!(f, "-"),
            BaseStatus::Modified(b) => write!(f, "{}", b),
            BaseStatus::Mismatch(b) => write!(f, "{}", b.char()),
            BaseStatus::Deletion => write!(f, "Deletion"),
        }
    }
}

impl Ord for BaseStatus {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match (self, other) {
            (BaseStatus::Canonical, BaseStatus::Canonical) => {
                std::cmp::Ordering::Equal
            }
            (BaseStatus::Canonical, _) => std::cmp::Ordering::Less,
            (_, BaseStatus::Canonical) => std::cmp::Ordering::Greater,
            (BaseStatus::Modified(mod1), BaseStatus::Modified(mod2)) => {
                mod1.cmp(mod2)
            }
            (BaseStatus::Modified(_), _) => std::cmp::Ordering::Less,
            (_, BaseStatus::Modified(_)) => std::cmp::Ordering::Greater,
            (BaseStatus::Mismatch(b1), BaseStatus::Mismatch(b2)) => b1.cmp(b2),
            (BaseStatus::Deletion, _) => std::cmp::Ordering::Greater,
            (_, BaseStatus::Deletion) => std::cmp::Ordering::Less,
        }
    }
}

impl PartialOrd for BaseStatus {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
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
    let positions = (start..end + 1).collect();

    Ok(GroundTruthSite { chrom, strand, base_status, positions })
}

type TidToChrom = HashMap<u32, String>;

type ChromStrandNamePositions =
    HashMap<String, HashMap<Strand, HashMap<BaseStatus, Vec<i64>>>>;

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
) -> anyhow::Result<ChromStrandNamePositions> {
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
        result
            .entry(ground_truth_site.chrom)
            .or_insert_with(HashMap::new)
            .entry(ground_truth_site.strand)
            .or_insert_with(HashMap::new)
            .entry(ground_truth_site.base_status)
            .or_insert_with(Vec::new)
            .extend(ground_truth_site.positions);
        lines_processed.inc(1);
    }
    if result.is_empty() {
        bail!("zero valid positions parsed from BED file".to_string());
    }
    lines_processed.finish_and_clear();
    info!("Processed {} BED lines", lines_processed.position());
    Ok(result)
}

// ground truth and observed base status pointing to vector of mod probabilities
type GtProbs = HashMap<(BaseStatus, BaseStatus), Vec<f32>>;

fn process_bam_record(
    record: &Record,
    mod_positions: &ChromStrandNamePositions,
    tid_to_chrom: &TidToChrom,
    collapse_method: Option<&CollapseMethod>,
    edge_filter: Option<&EdgeFilter>,
) -> anyhow::Result<GtProbs> {
    let mut result = HashMap::new();
    let mbi = ModBaseInfo::new_from_record(&record)?;
    let record_name = String::from_utf8(record.qname().to_vec())
        .unwrap_or("utf-decode-failed".to_string());
    if record.is_unmapped() || record_is_secondary(&record) {
        bail!("Unmapped or secondary");
    }

    let chrom = tid_to_chrom
        .get(&(record.tid() as u32))
        .ok_or_else(|| anyhow!("Invalid record TID: {}", record.tid(),))?;
    let cgt_mod_pos = mod_positions.get(chrom).ok_or_else(|| {
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
                    BaseModCall::Canonical(p) => (p, BaseStatus::Canonical),
                    BaseModCall::Modified(p, code) => {
                        (p, BaseStatus::Modified(code))
                    }
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
    mod_positions: &ChromStrandNamePositions,
    tid_to_chrom: &TidToChrom,
    collapse_method: Option<&CollapseMethod>,
    edge_filter: Option<&EdgeFilter>,
    suppress_pb: bool,
) -> anyhow::Result<GtProbs> {
    let lines_processed = get_ticker();
    if suppress_pb {
        lines_processed
            .set_draw_target(indicatif::ProgressDrawTarget::hidden());
    }
    lines_processed.set_message("Records processed");
    let (gt_mod_quals, errs) = reader.records().fold(
        (HashMap::new(), HashMap::new()),
        |(mut gt_mod_quals, mut errs), record| {
            match record {
                Ok(record) => {
                    match process_bam_record(
                        &record,
                        &mod_positions,
                        tid_to_chrom,
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
    info!("Processed {} mapping recrods", lines_processed.position());
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

fn balance_ground_truth(gt_mod_quals: &mut GtProbs) -> anyhow::Result<()> {
    // Calculate the total number of elements in each row
    let gt_totals: HashMap<_, _> = gt_mod_quals
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

    for ((gt_mod, _calls), probs) in gt_mod_quals.iter_mut() {
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
                    .filter_map(|idx| probs.get(idx).map(|x| *x))
                    .collect::<Vec<f32>>();
            *probs = keepers;
        }
    }
    Ok(())
}

fn print_table(gt_mod_quals: &GtProbs) {
    let mut gt_codes: Vec<_> = gt_mod_quals.keys().map(|&(k, _)| k).collect();
    gt_codes.sort();
    let call_codes: Vec<_> = gt_mod_quals.keys().map(|&(_, k)| k).collect();

    let mut all_codes: Vec<_> =
        gt_codes.iter().chain(call_codes.iter()).unique().collect();
    all_codes.sort();

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
    for &gt_code in &all_codes {
        let mut row = Row::empty();
        row.add_cell(cell!(&gt_code.to_string()));
        for &call_code in &all_codes {
            let vector_length = gt_mod_quals
                .get(&(*gt_code, *call_code))
                .map(|v| v.len())
                .unwrap_or(0);
            row.add_cell(cell!(&format!("{}", vector_length)));
        }
        count_tbl.add_row(row);
    }

    let mut metatable = Table::new();
    metatable.set_format(*prettytable::format::consts::FORMAT_CLEAN);
    metatable.add_row(row!("", cb->"Called Base"));
    metatable.add_row(row!(b->"\n\n\nGround\nTruth", count_tbl));

    // Print the table
    metatable.printstd();

    // todo also print percentages table
    // also print key metrics
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
        // todo determine canonnical base of interest

        let mut all_gt_probs = HashMap::new();
        for (bam_path, bed_indices) in bam_path_to_bed_indices {
            let mut reader = Reader::from_path(bam_path.as_path())?;
            reader.set_threads(self.threads)?;
            let tid_to_chrom = get_tid_to_chrom(&reader)?;
            info!(
                "Parsing mapping at {}",
                bam_path.to_str().unwrap_or("invalid-UTF-8")
            );

            for bed_idx in bed_indices {
                let gt_mod_probs = process_bam_file(
                    &mut reader,
                    &gt_positions[bed_idx],
                    &tid_to_chrom,
                    collapse_method.as_ref(),
                    edge_filter.as_ref(),
                    self.suppress_progress,
                )?;
                for ((gt_code, call_code), probs) in gt_mod_probs.iter() {
                    all_gt_probs
                        .entry((gt_code.clone(), call_code.clone()))
                        .or_insert_with(Vec::new)
                        .extend(probs.clone());
                }
            }
        }

        // sort prob vectors
        for ((_, _), probs) in all_gt_probs.iter_mut() {
            probs.sort_by_key(|&x| x.to_bits());
        }
        balance_ground_truth(&mut all_gt_probs)?;

        print_table(&all_gt_probs);

        let mut all_quals = Vec::<f32>::new();
        for (_, mod_quals) in all_gt_probs.iter() {
            all_quals.extend(mod_quals);
        }
        all_quals.sort_by(|x, y| x.partial_cmp(y).unwrap_or(Ordering::Equal));
        if all_quals.iter().any(|v| v.is_nan()) {
            bail!("Failed to compare values");
        }
        let thresh =
            percentile_linear_interp(&all_quals, self.filter_quantile)?;
        info!("Threshold: {}", thresh);

        // todo apply threshold and print table again

        Ok(())
    }
}
