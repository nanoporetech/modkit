use std::cmp::Ordering;
use std::collections::{HashMap, VecDeque};
use std::fmt::{Debug, Display, Formatter};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{anyhow, bail};
use derive_new::new;
use indicatif::ProgressBar;
use log::{debug, error};
use nom::character::complete::one_of;
use nom::multi::many0;
use nom::IResult;
use noodles::csi::index::{
    reference_sequence::bin::Chunk as IndexChunk, Index as CsiIndex,
};

use crate::parsing_utils::{
    consume_digit, consume_string, consume_string_spaces,
};
use crate::position_filter::Iv;

#[derive(new, Clone, Debug, Eq, PartialEq)]
pub(super) struct DmrInterval {
    pub(super) interval: Iv,
    pub(super) chrom: String,
    pub(super) name: String,
}

impl DmrInterval {
    pub(super) fn parse_bed_line(line: &str) -> IResult<&str, Self> {
        let (rest, chrom) = consume_string(line)?;
        let (rest, start) = consume_digit(rest)?;
        let (rest, stop) = consume_digit(rest)?;

        let (rest, interval, name) = many0(one_of(" \t\r\n"))(rest)
            .and_then(|(rest, _)| consume_string_spaces(rest))
            .map(|(rest, name)| {
                let interval = Iv {
                    start,
                    stop,
                    val: (),
                };
                (rest, interval, name)
            })
            .unwrap_or_else(|_| {
                let interval = Iv {
                    start,
                    stop,
                    val: (),
                };
                let name = format!("{}:{}-{}", chrom, start, stop);
                (rest, interval, name)
            });

        Ok((
            rest,
            Self {
                interval,
                chrom,
                name,
            },
        ))
    }

    pub(super) fn parse_str(line: &str) -> anyhow::Result<Self> {
        Self::parse_bed_line(line)
            .map(|(_, this)| this)
            .map_err(|e| anyhow!("{}", e.to_string()))
    }

    pub(super) fn start(&self) -> u64 {
        self.interval.start
    }

    pub(super) fn stop(&self) -> u64 {
        self.interval.stop
    }

    pub(super) fn get_index_chunks(
        &self,
        index: &CsiIndex,
        chrom_id: usize,
    ) -> std::io::Result<Vec<IndexChunk>> {
        let start =
            noodles::core::Position::new((self.start() + 1) as usize).unwrap();
        let end =
            noodles::core::Position::new((self.stop() + 1) as usize).unwrap();
        let interval = noodles::core::region::Interval::from(start..=end);
        index.query(chrom_id, interval)
    }
}

impl PartialOrd for DmrInterval {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(&other))
    }
}

impl Ord for DmrInterval {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.chrom.cmp(&other.chrom) {
            Ordering::Equal => self.interval.cmp(&other.interval),
            o @ _ => o,
        }
    }
}

impl Display for DmrInterval {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}:{}-{}", self.chrom, self.start(), self.stop())
    }
}

#[derive(new)]
pub(super) struct DmrChunk {
    pub(super) chrom_id: u32,
    pub(super) control_chunks: Vec<IndexChunk>,
    pub(super) exp_chunks: Vec<IndexChunk>,
    pub(super) dmr_interval: DmrInterval,
}

pub(super) struct DmrIntervalIter {
    control_fn: String,
    exp_fn: String,
    control_contig_lookup: Arc<HashMap<String, usize>>,
    exp_contig_lookup: HashMap<String, usize>,
    control_index: CsiIndex,
    exp_index: CsiIndex,
    regions_of_interest: VecDeque<DmrInterval>,
    chunk_size: usize,
    failures: ProgressBar,
}

impl DmrIntervalIter {
    pub(super) fn new(
        control_path: &PathBuf,
        exp_path: &PathBuf,
        control_contig_lookup: Arc<HashMap<String, usize>>,
        exp_contig_lookup: HashMap<String, usize>,
        control_index: CsiIndex,
        exp_index: CsiIndex,
        rois: VecDeque<DmrInterval>,
        chunk_size: usize,
        failure_counter: ProgressBar,
    ) -> Self {
        let control_fn = control_path
            .to_str()
            .map(|s| s.to_owned())
            .unwrap_or_else(|| format!("'a' failed path decode"));
        let exp_fn = exp_path
            .to_str()
            .map(|s| s.to_owned())
            .unwrap_or_else(|| format!("'a' failed path decode"));

        Self {
            control_fn,
            exp_fn,
            control_contig_lookup,
            exp_contig_lookup,
            control_index,
            exp_index,
            regions_of_interest: rois,
            chunk_size,
            failures: failure_counter,
        }
    }
}

impl Iterator for DmrIntervalIter {
    type Item = Vec<DmrChunk>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut chunks = Self::Item::with_capacity(self.chunk_size);
        loop {
            if let Some(dmr_interval) = self.regions_of_interest.pop_front() {
                let (control_chr_id, exp_chr_id) = match (
                    self.control_contig_lookup.get(&dmr_interval.chrom),
                    self.exp_contig_lookup.get(&dmr_interval.chrom),
                ) {
                    (Some(control_chr_id), Some(exp_chr_id)) => {
                        (*control_chr_id, *exp_chr_id)
                    }
                    (None, _) => {
                        self.failures.inc(1);
                        debug!(
                            "didn't find chrom id for {} in {} tabix header",
                            &self.control_fn, &dmr_interval.chrom
                        );
                        continue;
                    }
                    (_, None) => {
                        self.failures.inc(1);
                        debug!(
                            "didn't find chrom id for {} in {} tabix header",
                            &self.exp_fn, &dmr_interval.chrom
                        );
                        continue;
                    }
                };
                let control_chunks = dmr_interval
                    .get_index_chunks(&self.control_index, control_chr_id);
                let exp_chunks =
                    dmr_interval.get_index_chunks(&self.exp_index, exp_chr_id);
                let (control_chunks, exp_chunks) =
                    match (control_chunks, exp_chunks) {
                        (Ok(control_chunks), Ok(exp_chunks)) => {
                            (control_chunks, exp_chunks)
                        }
                        (Err(e), _) => {
                            self.failures.inc(1);
                            debug!(
                                "failed to index into {} bedMethyl \
                            for region {}, {}",
                                &self.control_fn,
                                dmr_interval,
                                e.to_string()
                            );
                            continue;
                        }
                        (_, Err(e)) => {
                            self.failures.inc(1);
                            debug!(
                                "failed to index into {} bedMethyl \
                            for chrom id {}, {}",
                                &self.exp_fn,
                                exp_chr_id,
                                e.to_string()
                            );
                            continue;
                        }
                    };
                let chunk = DmrChunk::new(
                    control_chr_id as u32,
                    control_chunks,
                    exp_chunks,
                    dmr_interval,
                );
                chunks.push(chunk);
                if chunks.len() >= self.chunk_size {
                    break;
                } else {
                    continue;
                }
            } else {
                break;
            }
        }
        if chunks.is_empty() {
            None
        } else {
            Some(chunks)
        }
    }
}

pub(super) fn parse_roi_bed<P: AsRef<Path>>(
    fp: P,
) -> anyhow::Result<Vec<DmrInterval>> {
    let intervals = BufReader::new(File::open(fp)?)
        .lines()
        .filter_map(|r| match r {
            Ok(l) => Some(l),
            Err(e) => {
                error!(
                    "error fetching line from regions BED, {}",
                    e.to_string()
                );
                None
            }
        })
        .map(|line| DmrInterval::parse_str(&line))
        .collect::<anyhow::Result<Vec<DmrInterval>>>()?;
    if intervals.is_empty() {
        bail!("didn't parse any regions")
    } else {
        Ok(intervals)
    }
}

#[cfg(test)]
mod dmr_util_tests {
    use crate::dmr::util::{parse_roi_bed, DmrInterval};
    use crate::position_filter::Iv;

    #[test]
    fn test_parse_rois() {
        let obs = DmrInterval::parse_str(
            "chr20\t279148\t279507\tCpG: 39 359\t39\t260\t21.7\t72.4\t0.83",
        )
        .unwrap();
        let expected = DmrInterval::new(
            Iv {
                start: 279148,
                stop: 279507,
                val: (),
            },
            "chr20".to_string(),
            "CpG: 39 359".to_string(),
        );
        assert_eq!(obs, expected);
        let obs = DmrInterval::parse_str("chr20\t279148\t279507\tCpGby_any_other_name\t39\t260\t21.7\t72.4\t0.83").unwrap();
        let expected = DmrInterval::new(
            Iv {
                start: 279148,
                stop: 279507,
                val: (),
            },
            "chr20".to_string(),
            "CpGby_any_other_name".to_string(),
        );
        assert_eq!(obs, expected);
        let obs = DmrInterval::parse_str("chr20\t279148\t279507\t").unwrap();
        let expected = DmrInterval::new(
            Iv {
                start: 279148,
                stop: 279507,
                val: (),
            },
            "chr20".to_string(),
            "chr20:279148-279507".to_string(),
        );
        assert_eq!(obs, expected);
        let obs = DmrInterval::parse_str("chr20\t279148\t279507 ").unwrap();
        assert_eq!(obs, expected);
    }

    #[test]
    fn test_roi_parsing() {
        let fp = "tests/resources/sim_cpg_regions.bed";
        let rois = parse_roi_bed(fp).unwrap();
        let expected = [
            DmrInterval {
                interval: Iv {
                    start: 10172120,
                    stop: 10172545,
                    val: (),
                },
                chrom: "chr20".to_string(),
                name: "r1".to_string(),
            },
            DmrInterval {
                interval: Iv {
                    start: 10217487,
                    stop: 10218336,
                    val: (),
                },
                chrom: "chr20".to_string(),
                name: "r2".to_string(),
            },
            DmrInterval {
                interval: Iv {
                    start: 10034963,
                    stop: 10035266,
                    val: (),
                },
                chrom: "chr20".to_string(),
                name: "r3".to_string(),
            },
        ]
        .to_vec();
        assert_eq!(rois, expected);
    }

    #[test]
    fn test_roi_parsing_noname() {
        let fp = "tests/resources/sim_cpg_regions_noname.bed";
        let rois = parse_roi_bed(fp).unwrap();
        let expected = [
            DmrInterval {
                interval: Iv {
                    start: 10172120,
                    stop: 10172545,
                    val: (),
                },
                chrom: "chr20".to_string(),
                name: "chr20:10172120-10172545".to_string(),
            },
            DmrInterval {
                interval: Iv {
                    start: 10217487,
                    stop: 10218336,
                    val: (),
                },
                chrom: "chr20".to_string(),
                name: "chr20:10217487-10218336".to_string(),
            },
            DmrInterval {
                interval: Iv {
                    start: 10034963,
                    stop: 10035266,
                    val: (),
                },
                chrom: "chr20".to_string(),
                name: "chr20:10034963-10035266".to_string(),
            },
        ]
        .to_vec();
        assert_eq!(rois, expected);
    }
}
