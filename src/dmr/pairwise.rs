use std::collections::HashMap;
use std::fs::File;
use std::io::BufRead;
use std::path::PathBuf;

use indicatif::ProgressBar;
use log::{debug, error};
use noodles::bgzf;
use noodles::csi::index::reference_sequence::bin::Chunk as IndexChunk;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::dmr::model::{AggregatedCounts, ModificationCounts};
use crate::dmr::util::{BedMethylLine, DmrInterval, DmrIntervalIter};
use crate::position_filter::{Iv, StrandedPositionFilter};
use crate::util::Strand;

fn aggregate_counts(
    bm_lines: &[BedMethylLine],
    chrom_id: u32,
    position_filter: &StrandedPositionFilter,
) -> anyhow::Result<AggregatedCounts> {
    let grouped_by_position: FxHashMap<u64, Vec<&BedMethylLine>> = bm_lines
        .iter()
        .filter(|bm_line| match bm_line.strand {
            '+' => position_filter.contains(
                chrom_id as i32,
                bm_line.start(),
                Strand::Positive,
            ),
            '-' => position_filter.contains(
                chrom_id as i32,
                bm_line.start(),
                Strand::Negative,
            ),
            '.' => position_filter.overlaps_not_stranded(
                chrom_id,
                bm_line.start(),
                bm_line.stop(),
            ),
            _ => {
                debug!(
                    "encountered illegal strand in bedmethyl {}",
                    bm_line.strand
                );
                false
            }
        })
        .fold(FxHashMap::default(), |mut acc, bm_line| {
            acc.entry(bm_line.start())
                .or_insert(Vec::new())
                .push(bm_line);
            acc
        });
    let (counts_per_code, total) = grouped_by_position.into_iter().fold(
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
            assert_eq!(valid_covs.len(), 1);
            assert_eq!(
                chroms.len(),
                1,
                "should only get 1 chrom, got {} {:?}, {:?}",
                chroms.len(),
                &chroms,
                &grouped
            );
            for x in grouped {
                *acc.entry(x.raw_mod_code).or_insert(0) +=
                    x.count_methylated as usize;
            }
            total_so_far += valid_coverage;
            (acc, total_so_far)
        },
    );
    AggregatedCounts::try_new(counts_per_code, total)
}

fn get_mod_counts_for_condition(
    reader: &mut bgzf::Reader<File>,
    chunks: &[IndexChunk],
    interval: &Iv,
    chrom_id: u32,
    position_filter: &StrandedPositionFilter,
) -> anyhow::Result<AggregatedCounts> {
    let mut bedmethyl_lines = Vec::new();
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
        bedmethyl_lines.extend(
            lines
                .into_iter()
                .filter_map(|l| BedMethylLine::parse(l.as_str()).ok())
                .filter(|bml| interval.overlap(bml.start(), bml.stop())),
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
    position_filter: &StrandedPositionFilter,
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
    )?;
    let experimental_counts = get_mod_counts_for_condition(
        &mut exp_bed_reader,
        &exp_chunks,
        &dmr_interval.interval,
        chrom_id,
        &position_filter,
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
    position_filter: StrandedPositionFilter,
    mut writer: Box<dyn std::io::Write>,
    pb: ProgressBar,
    successes: ProgressBar,
) -> anyhow::Result<()> {
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
                            Ok(_) => {}
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

    for result in rcv {
        match result {
            Ok(counts) => {
                writer.write(counts.to_row()?.as_bytes())?;
                successes.inc(1);
            }
            Err(e) => {
                debug!("unexpected error, {}", e.to_string());
            }
        }
    }

    Ok(())
}
