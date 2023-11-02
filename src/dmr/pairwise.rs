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
use crate::dmr::util::{DmrInterval, DmrIntervalIter};
use crate::mod_base_code::DnaBase;
use crate::position_filter::{Iv, StrandedPositionFilter};
use crate::util::{Strand, StrandRule};

fn aggregate_counts(
    bm_lines: &[BedMethylLine],
    chrom_id: u32,
    position_filter: &StrandedPositionFilter<DnaBase>,
) -> anyhow::Result<AggregatedCounts> {
    let grouped_by_position: FxHashMap<u64, Vec<&BedMethylLine>> = bm_lines
        .iter()
        .filter_map(|bm_line| {
            let base = match bm_line.strand {
                StrandRule::Positive => position_filter
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
                StrandRule::Both => position_filter
                    .get_base_at_position_stranded(
                        chrom_id as i32,
                        bm_line.start(),
                        Strand::Positive,
                    ),
            };
            base.and_then(|b| {
                if bm_line.check_base(b, None) {
                    // todo add user mappings
                    Some(bm_line)
                } else {
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
    AggregatedCounts::try_new(counts_per_code, total)
}

fn get_mod_counts_for_condition(
    reader: &mut bgzf::Reader<File>,
    chunks: &[IndexChunk],
    interval: &Iv,
    chrom_id: u32,
    position_filter: &StrandedPositionFilter<DnaBase>,
    filename: &PathBuf,
) -> anyhow::Result<AggregatedCounts> {
    let mut bedmethyl_lines = Vec::new();
    let mut failed_to_parse = 0;
    let mut successfully_parsed = 0usize;
    // todo could be a fold instead
    for chunk in chunks {
        reader.seek(chunk.start())?;
        let mut lines = Vec::new();
        'readloop: loop {
            let mut buf = String::new();
            let _byts = reader.read_line(&mut buf)?;
            lines.push(buf);
            let cur_pos = reader.virtual_position();
            if cur_pos >= chunk.end() {
                break 'readloop;
            }
        }
        let (n_fail, lines) = lines.into_iter().fold(
            (0usize, Vec::new()),
            |(n_fail, mut acc), line| match BedMethylLine::parse(line.as_str())
            {
                Ok(bm_line) => {
                    acc.push(bm_line);
                    (n_fail, acc)
                }
                Err(_) => (n_fail + 1, acc),
            },
        );
        failed_to_parse += n_fail;
        successfully_parsed += lines.len();

        bedmethyl_lines.extend(
            lines
                .into_iter()
                .filter(|bml| interval.overlap(bml.start(), bml.stop())),
        );
    }

    if successfully_parsed == 0 {
        bail!("failed to parse any bedMethyl lines from {:?}", filename);
    }

    if failed_to_parse > 0 {
        debug!(
            "failed to parse {} lines from {:?}",
            failed_to_parse, filename
        );
    }

    aggregate_counts(&bedmethyl_lines, chrom_id, position_filter)
}

pub(super) fn get_modification_counts(
    control_bedmethyl: &PathBuf,
    exp_bedmethyl: &PathBuf,
    control_chunks: &[IndexChunk],
    exp_chunks: &[IndexChunk],
    dmr_interval: DmrInterval,
    position_filter: &StrandedPositionFilter<DnaBase>,
    chrom_id: u32,
) -> anyhow::Result<ModificationCounts> {
    let mut control_reader =
        File::open(control_bedmethyl).map(bgzf::Reader::new)?;
    let mut exp_bed_reader =
        File::open(exp_bedmethyl).map(bgzf::Reader::new)?;
    if control_chunks.len() != 1 {
        debug!("more than 1 control chunk?, got {}", control_chunks.len());
    }
    if exp_chunks.len() != 1 {
        debug!("more than 1 exp chunk?, got {}", exp_chunks.len());
    }
    let control_counts = get_mod_counts_for_condition(
        &mut control_reader,
        &control_chunks,
        &dmr_interval.interval,
        chrom_id,
        &position_filter,
        &control_bedmethyl,
    )?;
    let experimental_counts = get_mod_counts_for_condition(
        &mut exp_bed_reader,
        &exp_chunks,
        &dmr_interval.interval,
        chrom_id,
        &position_filter,
        &exp_bedmethyl,
    )?;

    ModificationCounts::new(
        dmr_interval.start(),
        dmr_interval.stop(),
        control_counts,
        experimental_counts,
        dmr_interval,
    )
}

pub(super) fn run_pairwise_dmr(
    control_bed_fp: &PathBuf,
    exp_bed_fp: &PathBuf,
    dmr_interval_iter: DmrIntervalIter,
    position_filter: StrandedPositionFilter<DnaBase>,
    mut writer: Box<dyn std::io::Write>,
    pb: ProgressBar,
) -> anyhow::Result<usize> {
    let (snd, rcv) = crossbeam_channel::bounded(1000);
    let control_bedmethyl_fp = control_bed_fp.clone();
    let exp_bedmethyl_fp = exp_bed_fp.clone();

    std::thread::spawn(move || {
        for chunks in dmr_interval_iter {
            let mut result: Vec<anyhow::Result<ModificationCounts>> = vec![];
            let (res, _) = rayon::join(
                || {
                    chunks
                        .into_par_iter()
                        .map(|dmr_chunk| {
                            get_modification_counts(
                                &control_bedmethyl_fp,
                                &exp_bedmethyl_fp,
                                &dmr_chunk.control_chunks,
                                &dmr_chunk.exp_chunks,
                                dmr_chunk.dmr_interval,
                                &position_filter,
                                dmr_chunk.chrom_id,
                            )
                        })
                        .collect::<Vec<_>>()
                },
                || {
                    result.into_iter().for_each(|counts| {
                        pb.inc(1);
                        match snd.send(counts) {
                            Ok(_) => pb.inc(1),
                            Err(e) => error!(
                                "failed to send results, {}",
                                e.to_string()
                            ),
                        }
                    })
                },
            );
            result = res;
            result.into_iter().for_each(|counts| {
                pb.inc(1);
                match snd.send(counts) {
                    Ok(_) => {}
                    Err(e) => {
                        error!("failed to send results, {}", e.to_string())
                    }
                }
            });
        }
        pb.finish_and_clear();
    });

    let mut success_count = 0;
    for result in rcv {
        match result {
            Ok(counts) => {
                writer.write(counts.to_row()?.as_bytes())?;
                success_count += 1;
            }
            Err(e) => {
                debug!("unexpected error, {}", e.to_string());
            }
        }
    }

    Ok(success_count)
}
