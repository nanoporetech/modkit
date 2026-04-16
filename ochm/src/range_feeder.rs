use mod_kit::util::{GenomeRegion, Region};
use rust_htslib::bam::{self, Read};
use std::collections::{HashMap, VecDeque};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::ops::Range;
use std::path::PathBuf;

use crate::util::InferenceMessage;
use anyhow::{anyhow, bail, Context};
use itertools::Itertools;
use log::{debug, info, warn};

// TODO change everything to u64
const MAX_CHR_LENGTH: u64 = u32::MAX as u64;

type Tid = u32;
type ChrLength = u32;
pub(crate) struct RangeFeeder {
    work_queue: VecDeque<(Tid, Range<ChrLength>)>,
    curr_tid: Tid,
    curr_range: Range<ChrLength>,
    region_width: ChrLength,
    super_batch_size: usize,
    tid_to_name: HashMap<Tid, String>,
    curr_super_batch_size: usize,
    break_batch: bool,
    done: bool,
}

impl RangeFeeder {
    pub(crate) fn new(
        in_bam: &PathBuf,
        region: Option<&Region>,
        include_bed: Option<&PathBuf>,
        chunk_width: u32,
        step_size: u32,
        batch_size: u32,
        super_batch_size: u32,
    ) -> anyhow::Result<Self> {
        let reader = bam::Reader::from_path(in_bam)?;
        let header = reader.header();
        let name_to_tid = header
            .target_names()
            .iter()
            .enumerate()
            .map(|(tid, raw_name)| {
                (tid as u32, header.target_len(tid as u32).unwrap(), raw_name)
            })
            .map(|(tid, length, raw_name)| {
                let name = String::from_utf8_lossy(raw_name).to_string();
                let length = std::cmp::min(MAX_CHR_LENGTH, length);
                (name, (tid, length as u32))
            })
            .collect::<HashMap<String, (u32, u32)>>();

        let mut work_queue = if let Some(region) = region {
            if let Some((tid, length)) = name_to_tid.get(&region.name) {
                let end = if region.end > *length {
                    warn!(
                        "contig is not long enough for region {region:?}, \
                         only {length}"
                    );
                    *length
                } else {
                    region.end
                };
                let start = region.start;
                VecDeque::from([(*tid, start..end)])
            } else {
                bail!("didn't find contig {region} in header")
            }
        } else if let Some(bed_fp) = include_bed {
            let reader = BufReader::new(File::open(bed_fp).context(
                format!("failed to open BED regions at {bed_fp:?}"),
            )?);
            reader
                .lines()
                .map(|res| {
                    res.map_err(|e| anyhow!("failed to read, {e}")).and_then(
                        |line| GenomeRegion::parse_unstranded_bed_line(&line),
                    )
                })
                .filter_map_ok(|region| {
                    name_to_tid.get(&region.chrom).map(|(tid, _length)| {
                        (*tid, (region.start as u32)..(region.end as u32))
                    })
                })
                .collect::<anyhow::Result<VecDeque<(u32, Range<u32>)>>>()?
        } else {
            name_to_tid
                .iter()
                .map(|(_name, (tid, l))| (*tid, 0..(*l)))
                .collect()
        };
        // calculate total genome positions that will be in one batch after
        // chunked up
        let region_width = ((step_size * batch_size) - step_size) + chunk_width;
        let super_batch_size_bp = region_width * super_batch_size;
        info!(
            "collecting regions of {region_width}bp ({chunk_width} bp \
             chunks), super batches of {super_batch_size} \
             ({super_batch_size_bp}bp). Stepping {step_size} bp at a time."
        );

        let (curr_tid, curr_range) = work_queue
            .pop_front()
            .ok_or(anyhow!("must be at least one contig to work on"))?;

        let tid_to_name = name_to_tid
            .into_iter()
            .map(|(name, (tid, _))| (tid, name))
            .collect::<HashMap<u32, String>>();

        Ok(Self {
            work_queue,
            region_width,
            curr_tid,
            curr_range,
            tid_to_name,
            super_batch_size: super_batch_size as usize,
            done: false,
            break_batch: false,
            curr_super_batch_size: 0usize,
        })
    }

    #[inline]
    fn update(&mut self) {
        if let Some((tid, range)) = self.work_queue.pop_front() {
            // self.break_batch = self.curr_super_batch_size >=
            // self.super_batch_size;
            self.curr_tid = tid;
            self.curr_range = range;
            self.done = false;
        } else {
            self.done = true
        }
    }
}

impl Iterator for RangeFeeder {
    type Item = InferenceMessage<Vec<Region>>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }

        if self.break_batch {
            self.break_batch = false;
            self.curr_super_batch_size = 0;
            return Some(InferenceMessage::Done);
        }

        let mut regions = Vec::new();
        let mut front = self.curr_range.start;
        loop {
            if self.curr_super_batch_size >= self.super_batch_size {
                self.break_batch = true;
                break;
            }
            let end =
                std::cmp::min(front + self.region_width, self.curr_range.end);
            let chrom = self.tid_to_name.get(&self.curr_tid).unwrap().clone();
            let region = Region::new(chrom, front, end); //todo might want a +1 here?
            regions.push(region);
            self.curr_super_batch_size += 1;
            if end >= self.curr_range.end {
                self.update();
                if self.done {
                    break;
                } else {
                    front = self.curr_range.start;
                    continue;
                }
            } else if self.curr_super_batch_size >= self.super_batch_size {
                self.curr_range.start = end;
                self.break_batch = true;
                debug!(
                    "breaking super batch, {} regions, starting at {front} \
                     next time",
                    regions.len()
                );
                break;
            } else {
                front = end;
            }
        }

        assert!(!regions.is_empty(), "regions should not be empty?");
        Some(InferenceMessage::Work(regions))
    }
}

#[cfg(test)]
mod range_feeder_tests {
    use crate::range_feeder::RangeFeeder;
    use common_macros::hash_map;
    use std::collections::VecDeque;

    #[test]
    fn test_range_feeder_regions() {
        let mut work_queue = vec![
            (1u32, 0..20),
            (1u32, 50..62),
            (1u32, 80..100),
            (2u32, 100..112),
        ]
        .into_iter()
        .collect::<VecDeque<_>>();
        // let mut work_queue = vec![
        //     (1u32, 0..234),
        // ].into_iter().collect::<VecDeque<_>>();
        let (curr_tid, curr_range) = work_queue.pop_front().unwrap();
        let chunk_width = 10u32;
        let overlap = 5u32;
        let step_size = chunk_width - overlap;
        assert_eq!(step_size, 5);
        let batch_size = 5u32;
        let region_width = ((step_size * batch_size) - step_size) + chunk_width;
        // dbg!(region_width);
        let tid_to_name = hash_map! {
            1u32 => "chr1".to_string(),
            2u32 => "chr2".to_string()
        };
        let range_feeder = RangeFeeder {
            work_queue,
            curr_tid,
            curr_range,
            region_width,
            super_batch_size: 2,
            tid_to_name,
            curr_super_batch_size: 0,
            break_batch: false,
            done: false,
        };

        // for batch in range_feeder {
        //     dbg!(batch);
        // }
    }
}
