use std::collections::{HashMap, HashSet, VecDeque};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use anyhow::bail;
use itertools::Itertools;
use log::info;
use log_once::info_once;
use rust_htslib::bam::{self, Read};
use rust_lapper as lapper;
use rustc_hash::FxHashMap;

use crate::mod_base_code::DnaBase;
use crate::util::{get_targets, get_ticker, ReferenceRecord, Strand};

pub(crate) type Iv = lapper::Interval<u64, ()>;
pub(crate) type GenomeIntervals<T> = lapper::Lapper<u64, T>;

#[derive(Debug, Clone)]
pub struct StrandedPositionFilter<T: Send + Sync + Eq + Clone> {
    pub(crate) pos_positions: FxHashMap<u32, GenomeIntervals<T>>,
    pub(crate) neg_positions: FxHashMap<u32, GenomeIntervals<T>>,
}

impl<T: Send + Sync + Eq + Clone> StrandedPositionFilter<T> {
    pub fn contains(
        &self,
        chrom_id: i32,
        position: u64,
        strand: Strand,
    ) -> bool {
        let positions = match strand {
            Strand::Positive => &self.pos_positions,
            Strand::Negative => &self.neg_positions,
        };
        positions
            // todo(arand) chromId should really be an enum.. encoding things as
            // missing by making them  negative numbers is so.. C
            .get(&(chrom_id as u32))
            .map(|lp| lp.find(position, position + 1).count() > 0)
            .unwrap_or(false)
    }

    pub fn overlaps_not_stranded(
        &self,
        chrom_id: u32,
        start: u64,
        end: u64,
    ) -> bool {
        // check pos positions first, if overlaps with positive positions
        // eagerly return true otherwise check negative overlaps
        let pos_overlaps = self
            .pos_positions
            .get(&chrom_id)
            .map(|lp| lp.find(start, end).count() > 0)
            .unwrap_or(false);
        if !pos_overlaps {
            self.neg_positions
                .get(&chrom_id)
                .map(|lp| lp.find(start, end).count() > 0)
                .unwrap_or(false)
        } else {
            true
        }
    }

    pub fn contains_chrom_id(&self, chrom_id: &i64) -> bool {
        if *chrom_id < 0 {
            false
        } else {
            let chrom_id = *chrom_id as u32;
            self.pos_positions.contains_key(&chrom_id)
                || self.neg_positions.contains_key(&chrom_id)
        }
    }

    pub fn contig_ends(&self, contig_id: &u32) -> Option<(u64, u64)> {
        let get_start_end = |positions: &FxHashMap<u32, GenomeIntervals<T>>| -> Option<(u64, u64)> {
            positions.get(&contig_id)
                .and_then(|lp| {
                    let start = lp.intervals.first().map(|iv| iv.start);
                    let stop = lp.intervals.last().map(|iv| iv.stop);
                    match (start, stop) {
                        (Some(s), Some(t)) => Some((s, t)),
                        _ => None
                    }
                })
        };
        let pos_ends = get_start_end(&self.pos_positions);
        let neg_ends = get_start_end(&self.neg_positions);

        match (pos_ends, neg_ends) {
            (Some((a, b)), Some((c, d))) => {
                let start = std::cmp::min(a, c);
                let end = std::cmp::max(b, d);
                Some((start, end))
            }
            _ => None,
        }
    }

    fn group_genome_intervals(
        genome_intervals: GenomeIntervals<T>,
        reference_record: &ReferenceRecord,
        interval_size: u64,
    ) -> Vec<ReferenceRecord> {
        let mut intervals =
            genome_intervals.intervals.into_iter().collect::<VecDeque<_>>();
        if intervals.is_empty() {
            return Vec::new();
        }

        let mut current = intervals.pop_front().unwrap();
        let mut agg = Vec::new();
        while let Some(iv) = intervals.pop_front() {
            if current.stop.checked_sub(current.start).expect(&format!(
                "invalid interval coordinates, {}:{}",
                current.start, current.stop
            )) > interval_size
            {
                let finished = std::mem::replace(&mut current, iv);
                agg.push(finished);
                continue;
            }
            current.stop = iv.stop
        }
        agg.push(current);
        agg.into_iter()
            .map(|iv| {
                let start = iv.start as u32;
                let length = iv
                    .stop
                    .checked_sub(iv.start)
                    .expect("invalid final interval")
                    as u32;
                ReferenceRecord::new(
                    reference_record.tid,
                    start,
                    length,
                    reference_record.name.clone(),
                )
            })
            .collect()
    }

    pub(crate) fn optimize_reference_records(
        &self,
        reference_records: Vec<ReferenceRecord>,
        interval_size: u32,
    ) -> Vec<ReferenceRecord> {
        let lut = reference_records
            .into_iter()
            .map(|rec| (rec.tid, rec))
            .collect::<HashMap<u32, ReferenceRecord>>();

        let contig_ids = self
            .pos_positions
            .keys()
            .chain(self.neg_positions.keys())
            .unique()
            .copied()
            .collect::<Vec<u32>>();

        contig_ids
            .iter()
            // shouldn't really need this filter
            .filter_map(|tid| lut.get(tid))
            .flat_map(|ref_record| {
                let tid = &ref_record.tid;
                let mut pos = self
                    .pos_positions
                    .get(tid)
                    .map(|ivs| ivs.intervals.clone())
                    .unwrap_or_else(Vec::new);
                let mut neg = self
                    .neg_positions
                    .get(tid)
                    .map(|ivs| ivs.intervals.clone())
                    .unwrap_or_else(Vec::new);
                pos.append(&mut neg);
                #[cfg(debug_assertions)]
                {
                    for iv in pos.iter() {
                        assert!(iv.start <= iv.stop);
                    }
                }
                let mut genome_intervals = GenomeIntervals::new(pos);
                #[cfg(debug_assertions)]
                {
                    for iv in genome_intervals.intervals.iter() {
                        assert!(iv.start <= iv.stop);
                    }
                }
                genome_intervals.merge_overlaps();
                #[cfg(debug_assertions)]
                {
                    for iv in genome_intervals.intervals.iter() {
                        assert!(iv.start <= iv.stop);
                    }
                }
                Self::group_genome_intervals(
                    genome_intervals,
                    ref_record,
                    interval_size as u64,
                )
            })
            .collect()
    }
}

impl StrandedPositionFilter<()> {
    pub fn from_bam_and_bed(
        bam_fp: &PathBuf,
        bed_fp: &PathBuf,
        suppress_pb: bool,
    ) -> anyhow::Result<Self> {
        let bam_reader = bam::Reader::from_path(bam_fp)?;
        let targets = get_targets(bam_reader.header(), None);
        let chrom_to_tid = targets
            .iter()
            .map(|reference_record| {
                (reference_record.name.as_str(), reference_record.tid)
            })
            .collect::<HashMap<&str, u32>>();
        Self::from_bed_file(bed_fp, &chrom_to_tid, suppress_pb)
    }

    pub fn from_bed_file(
        bed_fp: &PathBuf,
        chrom_to_target_id: &HashMap<&str, u32>,
        suppress_pb: bool,
    ) -> anyhow::Result<Self> {
        info!("parsing BED at {}", bed_fp.to_str().unwrap_or("invalid-UTF-8"));

        let fh = File::open(bed_fp)?;
        let mut pos_positions = FxHashMap::default();
        let mut neg_positions = FxHashMap::default();
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
            let raw_start = parts[1].parse::<u64>();
            let raw_end = parts[2].parse::<u64>();
            let (start, stop) = match (raw_start, raw_end) {
                (Ok(start), Ok(end)) => (start, end),
                _ => {
                    info!(
                        "improperly formatted BED line, failed to parse start \
                         and/or stop, {line}"
                    );
                    continue;
                }
            };
            let (pos_strand, neg_strand) = match parts.len() {
                3 => {
                    // BED3 use both strands
                    info_once!(
                        "a BED3 line encountered, strand with be set to BOTH"
                    );
                    (true, true)
                }
                6.. => {
                    // BED6+ use specified strand
                    match parts[5] {
                        "+" => (true, false),
                        "-" => (false, true),
                        "." => (true, true),
                        _ => {
                            info!(
                                "improperly formatted strand field {}",
                                &parts[5]
                            );
                            continue;
                        }
                    }
                }
                _ => {
                    info!(
                        "improperly formatted BED line, must be BED3 or BED6 \
                         {line}"
                    );
                    continue;
                }
            };
            debug_assert!(start <= stop, "start should be before stop");
            if let Some(chrom_id) = chrom_to_target_id.get(chrom_name) {
                if pos_strand {
                    pos_positions
                        .entry(*chrom_id)
                        .or_insert(Vec::new())
                        .push(Iv { start, stop, val: () })
                }
                if neg_strand {
                    neg_positions
                        .entry(*chrom_id)
                        .or_insert(Vec::new())
                        .push(Iv { start, stop, val: () })
                }
                lines_processed.inc(1);
            } else {
                info!("skipping chrom {chrom_name}, not present in BAM header");
                warned.insert(chrom_name.to_owned());
                continue;
            }
        }
        if pos_positions.is_empty() && neg_positions.is_empty() {
            bail!("zero valid positions parsed from BED file")
        }

        let pos_intervals = pos_positions
            .into_iter()
            .map(|(chrom_id, intervals)| {
                let mut lp = lapper::Lapper::new(intervals);
                lp.merge_overlaps();
                (chrom_id, lp)
            })
            .collect::<FxHashMap<u32, GenomeIntervals<()>>>();

        let neg_intervals = neg_positions
            .into_iter()
            .map(|(chrom_id, intervals)| {
                let mut lp = lapper::Lapper::new(intervals);
                lp.merge_overlaps();
                (chrom_id, lp)
            })
            .collect::<FxHashMap<u32, GenomeIntervals<()>>>();

        lines_processed.finish_and_clear();
        info!("processed {} BED lines", lines_processed.position());

        Ok(Self { pos_positions: pos_intervals, neg_positions: neg_intervals })
    }
}

impl StrandedPositionFilter<DnaBase> {
    pub fn get_base_at_position_stranded(
        &self,
        chrom_id: i32,
        position: u64,
        strand: Strand,
    ) -> Option<DnaBase> {
        let positions = match strand {
            Strand::Positive => &self.pos_positions,
            Strand::Negative => &self.neg_positions,
        };
        positions.get(&(chrom_id as u32)).and_then(|lp| {
            lp.find(position, position + 1).map(|x| x.val).next()
        })
    }
}
