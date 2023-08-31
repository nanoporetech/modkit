use crate::util::{get_ticker, Strand};
use anyhow::bail;
use log::info;
use rust_lapper as lapper;
use rustc_hash::FxHashMap;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

type Iv = lapper::Interval<u64, ()>;
type GenomeLapper = lapper::Lapper<u64, ()>;

#[derive(Debug)]
pub struct StrandedPositionFilter {
    pos_positions: FxHashMap<u32, GenomeLapper>,
    neg_positions: FxHashMap<u32, GenomeLapper>,
}

impl StrandedPositionFilter {
    pub fn from_bed_file(
        bed_fp: &PathBuf,
        chrom_to_target_id: &HashMap<&str, u32>,
        suppress_pb: bool,
    ) -> anyhow::Result<Self> {
        info!(
            "parsing BED at {}",
            bed_fp.to_str().unwrap_or("invalid-UTF-8")
        );

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
        for line in reader
            .lines()
            .filter_map(|l| l.ok())
            .filter(|l| !l.is_empty())
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
            let (start, stop) = match (raw_start, raw_end) {
                (Ok(start), Ok(end)) => (*start, *end),
                _ => {
                    info!("improperly formatted BED line {line}");
                    continue;
                }
            };
            let (pos_strand, neg_strand) = match parts[5] {
                "+" => (true, false),
                "-" => (false, true),
                "." => (true, true),
                _ => {
                    info!("improperly formatted strand field {}", &parts[5]);
                    continue;
                }
            };
            if let Some(chrom_id) = chrom_to_target_id.get(chrom_name) {
                if pos_strand {
                    pos_positions.entry(*chrom_id).or_insert(Vec::new()).push(
                        Iv {
                            start,
                            stop,
                            val: (),
                        },
                    )
                }
                if neg_strand {
                    neg_positions.entry(*chrom_id).or_insert(Vec::new()).push(
                        Iv {
                            start,
                            stop,
                            val: (),
                        },
                    )
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

        let pos_lapper = pos_positions
            .into_iter()
            .map(|(chrom_id, intervals)| {
                let mut lp = lapper::Lapper::new(intervals);
                lp.merge_overlaps();
                (chrom_id, lp)
            })
            .collect::<FxHashMap<u32, GenomeLapper>>();

        let neg_lapper = neg_positions
            .into_iter()
            .map(|(chrom_id, intervals)| {
                let mut lp = lapper::Lapper::new(intervals);
                lp.merge_overlaps();
                (chrom_id, lp)
            })
            .collect::<FxHashMap<u32, GenomeLapper>>();

        lines_processed.finish_and_clear();
        info!("processed {} BED lines", lines_processed.position());

        Ok(Self {
            pos_positions: pos_lapper,
            neg_positions: neg_lapper,
        })
    }

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
            // todo(arand) chromId should really be an enum.. encoding things as missing by making them
            //  negative numbers is so.. C
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
        // check pos positions first, if overlaps with positive positions eagerly return true
        // otherwise check negative overlaps
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
}
