use std::collections::HashMap;
use std::fs::File;
use std::io::BufRead;
use std::path::PathBuf;

use anyhow::bail;
use indicatif::ProgressBar;
use itertools::{Itertools, MinMaxResult};
use log::{debug, error};
use log_once::debug_once;
use noodles::bgzf;
use noodles::csi::index::reference_sequence::bin::Chunk as IndexChunk;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::dmr::bedmethyl::BedMethylLine;
use crate::dmr::model::{AggregatedCounts, ModificationCounts};
use crate::dmr::util::{DmrBatch, DmrIntervalIter};
use crate::mod_base_code::DnaBase;
use crate::position_filter::StrandedPositionFilter;
use crate::util::{Strand, StrandRule};

fn aggregate_counts(
    bm_lines: &[&BedMethylLine],
    chrom_id: u32,
    position_filter: &StrandedPositionFilter<DnaBase>,
) -> anyhow::Result<AggregatedCounts> {
    let grouped_by_position: FxHashMap<u64, Vec<&BedMethylLine>> = bm_lines
        .iter()
        .filter_map(|&bm_line| {
            let base = match bm_line.strand {
                StrandRule::Positive | StrandRule::Both => position_filter
                    .get_base_at_position_stranded(
                        chrom_id as i32,
                        bm_line.start(),
                        Strand::Positive,
                    ),
                StrandRule::Negative => position_filter
                    .get_base_at_position_stranded(
                        chrom_id as i32,
                        bm_line.start(),
                        Strand::Negative,
                    )
                    .map(|b| b.complement()),
            };
            base.and_then(|b| {
                if bm_line.check_base(b, None) {
                    // todo add user mappings
                    Some(bm_line)
                } else {
                    // this will eventually change when user-input is accepted
                    if !bm_line.check_mod_code_supported() {
                        debug_once!(
                            "encountered modification code {} in bedMethyl record, \
                             not currently supported",
                        bm_line.raw_mod_code
                    );

                    }
                    None
                }
            })
        })
        .fold(FxHashMap::default(), |mut acc, bm_line| {
            acc.entry(bm_line.start())
                .or_insert(Vec::new())
                .push(bm_line);
            acc
        });
    let (counts_per_code, total) = grouped_by_position.into_iter().try_fold(
        (HashMap::new(), 0),
        |(mut acc, mut total_so_far), (_pos, grouped)| {
            let valid_covs = grouped
                .iter()
                .map(|bml| bml.valid_coverage)
                .collect::<FxHashSet<u64>>();
            let chroms = grouped
                .iter()
                .map(|bml| &bml.chrom)
                .collect::<FxHashSet<&String>>();
            let valid_coverage = grouped[0].valid_coverage as usize;
            if valid_covs.len() != 1 {
                let mut message = format!("invalid data found, should not have more than 1 score per position for a \
                    base.");
                match grouped.iter().minmax_by(|a, b| a.start().cmp(&b.start())) {
                    MinMaxResult::NoElements => {},
                    MinMaxResult::MinMax(s, t) => {
                        message.push_str(&format!("starting at {}, ending at {}", s.start(), t.stop()))
                    }
                    MinMaxResult::OneElement(s) => {
                        message.push_str(&format!("starting at {}", s.start()))
                    }
                }
                bail!(message)
            }

            if chroms.len() != 1 {
                bail!(format!("should only get one chrom, got {}", chroms.len()))
            }
            for x in grouped {
                *acc.entry(x.raw_mod_code).or_insert(0) +=
                    x.count_methylated as usize;
            }
            total_so_far += valid_coverage;
            Ok((acc, total_so_far))
        },
    )?;
    match AggregatedCounts::try_new(counts_per_code, total) {
        Ok(x) => Ok(x),
        Err(e) => {
            debug!("{chrom_id}, {:?}\n{}", bm_lines, e.to_string());
            Err(e)
        }
    }
}

fn read_bedmethyl(
    fp: &PathBuf,
    chunks: &[IndexChunk],
) -> anyhow::Result<Vec<BedMethylLine>> {
    let mut reader = File::open(fp).map(bgzf::Reader::new)?;
    let mut bedmethyl_lines = Vec::new();
    let mut failed_to_parse = 0;
    let mut successfully_parsed = 0usize;
    // todo could be a fold instead
    for chunk in chunks {
        reader.seek(chunk.start())?;
        let mut lines = Vec::new();
        // todo come back and make this one loop
        'readloop: loop {
            let mut buf = String::new();
            let _byts = reader.read_line(&mut buf)?;
            lines.push(buf);
            let cur_pos = reader.virtual_position();
            if cur_pos >= chunk.end() {
                break 'readloop;
            }
        }
        for line in lines.iter() {
            match BedMethylLine::parse(line) {
                Ok(bm_line) => {
                    bedmethyl_lines.push(bm_line);
                    successfully_parsed += 1;
                }
                Err(_e) => {
                    // trace!("failed to parse line {line}");
                    failed_to_parse += 1
                }
            }
        }
    }

    if failed_to_parse > 0 {
        debug!("failed to parse {failed_to_parse} lines from {fp:?}");
    }

    if successfully_parsed == 0 {
        bail!("failed to parse any bedMethyl lines from {fp:?}",);
    }

    Ok(bedmethyl_lines)
}

pub(super) fn get_modification_counts(
    control_bedmethyl: &PathBuf,
    exp_bedmethyl: &PathBuf,
    dmr_batch: DmrBatch,
    position_filter: &StrandedPositionFilter<DnaBase>,
) -> anyhow::Result<Vec<anyhow::Result<ModificationCounts>>> {
    let control_bedmethyl_lines =
        read_bedmethyl(control_bedmethyl, &dmr_batch.get_control_chunks())?
            .into_iter()
            .fold(FxHashMap::default(), |mut acc, bm_line| {
                acc.entry(bm_line.chrom.clone())
                    .or_insert(Vec::new())
                    .push(bm_line);
                acc
            });
    let exp_bedmethyl_lines =
        read_bedmethyl(exp_bedmethyl, &dmr_batch.get_exp_chunks())?
            .into_iter()
            .fold(FxHashMap::default(), |mut acc, bm_line| {
                acc.entry(bm_line.chrom.clone())
                    .or_insert(Vec::new())
                    .push(bm_line);
                acc
            });

    let modification_counts_results = dmr_batch.dmr_chunks
        .into_par_iter()
        .map(|dmr_chunk| {
            let control = control_bedmethyl_lines
                .get(&dmr_chunk.dmr_interval.chrom)
                .map(|lines| {
                    lines
                        .iter()
                        .filter(|l| dmr_chunk.dmr_interval.interval.overlap(l.start(), l.stop()))
                        .collect::<Vec<&BedMethylLine>>()
                })
                .unwrap_or_else(|| Vec::new());
            let exp = exp_bedmethyl_lines
                .get(&dmr_chunk.dmr_interval.chrom)
                .map(|lines| {
                    lines
                        .iter()
                        .filter(|l| dmr_chunk.dmr_interval.interval.overlap(l.start(), l.stop()))
                        .collect::<Vec<&BedMethylLine>>()
                })
                .unwrap_or_else(|| Vec::new());
            let control_counts = aggregate_counts(&control, dmr_chunk.chrom_id, &position_filter);
            let exp_counts = aggregate_counts(&exp, dmr_chunk.chrom_id, &position_filter);
            match (control_counts, exp_counts) {
                (Ok(control_counts), Ok(exp_counts)) => {
                    ModificationCounts::new(
                        dmr_chunk.dmr_interval.start(),
                        dmr_chunk.dmr_interval.stop(),
                        control_counts,
                        exp_counts,
                        dmr_chunk.dmr_interval.clone())
                },
                (Err(e), Err(f)) => {
                    bail!("failed to aggregate control counts, {} and experimental counts, {}",
                        e.to_string(), f.to_string())
                }
                (Err(e), _) => {
                    bail!("failed to aggregate control counts, {}",
                        e.to_string())
                }
                (_, Err(e)) => {
                    bail!("failed to aggregate experiment counts, {}",
                        e.to_string())
                }
            }

        }).collect::<Vec<Result<ModificationCounts, _>>>();

    Ok(modification_counts_results)
}

pub(super) fn run_pairwise_dmr(
    control_bed_fp: &PathBuf,
    exp_bed_fp: &PathBuf,
    dmr_interval_iter: DmrIntervalIter,
    position_filter: StrandedPositionFilter<DnaBase>,
    mut writer: Box<dyn std::io::Write>,
    pb: ProgressBar,
    failure_counter: ProgressBar,
) -> anyhow::Result<usize> {
    let (snd, rcv) = crossbeam_channel::bounded(1000);
    let control_bedmethyl_fp = control_bed_fp.clone();
    let exp_bedmethyl_fp = exp_bed_fp.clone();

    std::thread::spawn(move || {
        for batch in dmr_interval_iter {
            match get_modification_counts(
                &control_bedmethyl_fp,
                &exp_bedmethyl_fp,
                batch,
                &position_filter,
            ) {
                Ok(results) => match snd.send(results) {
                    Ok(_) => {}
                    Err(e) => {
                        error!("failed to send results, {}", e.to_string())
                    }
                },
                Err(e) => {
                    debug!("failed entire dmr batch, {}", e.to_string())
                }
            }
        }
    });

    let mut success_count = 0;
    for results in rcv {
        for result in results {
            match result {
                Ok(counts) => {
                    writer.write(counts.to_row()?.as_bytes())?;
                    success_count += 1;
                    pb.inc(1);
                }
                Err(e) => {
                    debug!("unexpected error, {}", e.to_string());
                    failure_counter.inc(1);
                }
            }
        }
    }

    pb.finish_and_clear();

    Ok(success_count)
}
