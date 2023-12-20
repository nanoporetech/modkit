use std::cmp::Ordering;
use std::collections::{BTreeSet, HashMap, VecDeque};
use std::fmt::{Debug, Display, Formatter};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{anyhow, bail};
use clap::ValueEnum;
use derive_new::new;
use indicatif::ProgressBar;

use log::{debug, error};
use log_once::debug_once;
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

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
#[allow(non_camel_case_types)]
pub(super) enum HandleMissing {
    quiet,
    warn,
    fail,
}

impl Display for HandleMissing {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            HandleMissing::quiet => write!(f, "quiet"),
            HandleMissing::warn => write!(f, "warn"),
            HandleMissing::fail => write!(f, "fail"),
        }
    }
}

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
                let interval = Iv { start, stop, val: () };
                (rest, interval, name)
            })
            .unwrap_or_else(|_| {
                let interval = Iv { start, stop, val: () };
                let name = format!("{}:{}-{}", chrom, start, stop);
                (rest, interval, name)
            });

        Ok((rest, Self { interval, chrom, name }))
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

#[derive(new, Debug, Eq, PartialEq, PartialOrd)]
pub(super) struct DmrChunk {
    pub(super) chrom_id: u32,
    pub(super) dmr_interval: DmrInterval,
}

impl Ord for DmrChunk {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.chrom_id.cmp(&other.chrom_id) {
            Ordering::Equal => self.dmr_interval.cmp(&other.dmr_interval),
            o @ _ => o,
        }
    }
}

impl Display for DmrChunk {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}:{}-{}",
            self.dmr_interval.chrom,
            self.dmr_interval.start(),
            self.dmr_interval.stop()
        )
    }
}

pub(super) struct DmrIntervalIter {
    control_fn: String,
    exp_fn: String,
    control_contig_lookup: Arc<ContigLookup>,
    exp_contig_lookup: ContigLookup,
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
        control_contig_lookup: Arc<ContigLookup>,
        exp_contig_lookup: ContigLookup,
        control_index: CsiIndex,
        exp_index: CsiIndex,
        rois: VecDeque<DmrInterval>,
        chunk_size: usize,
        failure_counter: ProgressBar,
        handle_missing: HandleMissing,
    ) -> anyhow::Result<Self> {
        let control_fn = control_path
            .to_str()
            .map(|s| s.to_owned())
            .unwrap_or_else(|| format!("'a' failed path decode"));
        let exp_fn = exp_path
            .to_str()
            .map(|s| s.to_owned())
            .unwrap_or_else(|| format!("'a' failed path decode"));
        let regions_of_interest =
            rois.into_iter().try_fold(Vec::new(), |mut acc, roi| {
                let a_found =
                    control_contig_lookup.inner.contains_key(&roi.chrom);
                let b_found = exp_contig_lookup.inner.contains_key(&roi.chrom);
                if a_found && b_found {
                    acc.push(roi);
                    Ok(acc)
                } else {
                    let which = if b_found {
                        if let Some(a_name) =
                            control_contig_lookup.sample_name.as_ref()
                        {
                            format!("'{a_name}'")
                        } else {
                            format!("{:?}", &control_contig_lookup.file_path)
                        }
                    } else {
                        match (
                            control_contig_lookup.sample_name.as_ref(),
                            exp_contig_lookup.sample_name.as_ref(),
                        ) {
                            (Some(a_name), Some(b_name)) => {
                                format!("'{a_name}' and '{b_name}'")
                            }
                            _ => {
                                format!(
                                    "{:?} and {:?}",
                                    &control_contig_lookup.file_path,
                                    &exp_contig_lookup.file_path
                                )
                            }
                        }
                    };
                    match handle_missing {
                        HandleMissing::quiet => Ok(acc),
                        HandleMissing::warn => {
                            debug_once!(
                                "chrom {} is missing from {which} bedMethyl \
                                 index, discarding",
                                &roi.chrom
                            );
                            Ok(acc)
                        }
                        HandleMissing::fail => {
                            let message = format!(
                                "chrom {} is missing from {which} bedMethyl \
                                 index, fatal error",
                                &roi.chrom
                            );
                            error!("{message}");
                            bail!(message)
                        }
                    }
                }
            })?;

        if regions_of_interest.is_empty() {
            bail!("no valid regions in input")
        }

        Ok(Self {
            control_fn,
            exp_fn,
            control_contig_lookup,
            exp_contig_lookup,
            control_index,
            exp_index,
            regions_of_interest: regions_of_interest.into_iter().collect(),
            chunk_size,
            failures: failure_counter,
        })
    }
}

#[derive(new, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
struct ProtoIndexChunk {
    start: u64,
    stop: u64,
}

impl Into<IndexChunk> for ProtoIndexChunk {
    fn into(self) -> IndexChunk {
        IndexChunk::new(self.start.into(), self.stop.into())
    }
}

impl From<IndexChunk> for ProtoIndexChunk {
    fn from(value: IndexChunk) -> Self {
        let start: u64 = value.start().into();
        let stop: u64 = value.end().into();
        Self { start, stop }
    }
}

#[derive(Default)]
pub(super) struct DmrBatch {
    pub(super) dmr_chunks: Vec<DmrChunk>, // todo could make this DmrInterval?
    control_chunks: BTreeSet<ProtoIndexChunk>,
    experiment_chunks: BTreeSet<ProtoIndexChunk>,
}

impl DmrBatch {
    fn append(
        &mut self,
        dmr_chunk: DmrChunk,
        control_index_chunks: Vec<IndexChunk>,
        experiment_index_chunks: Vec<IndexChunk>,
    ) {
        self.dmr_chunks.push(dmr_chunk);
        for x in control_index_chunks {
            self.control_chunks.insert(x.into());
        }
        for x in experiment_index_chunks {
            self.experiment_chunks.insert(x.into());
        }
    }

    fn size(&self) -> usize {
        self.dmr_chunks.len()
    }

    fn is_empty(&self) -> bool {
        self.dmr_chunks.is_empty()
    }

    #[inline]
    fn proto_iter(chunks: &BTreeSet<ProtoIndexChunk>) -> Vec<IndexChunk> {
        if let Some(&first) = chunks.first() {
            let (mut chunks, last) = chunks.iter().skip(1).fold(
                (Vec::new(), first),
                |(mut acc, mut front), chunk| {
                    if chunk.start <= front.stop {
                        front = ProtoIndexChunk::new(front.start, chunk.stop);
                        (acc, front)
                    } else {
                        acc.push(front);
                        (acc, *chunk)
                    }
                },
            );
            chunks.push(last);
            chunks.into_iter().map(|x| x.into()).collect()
        } else {
            Vec::new()
        }
    }

    pub(super) fn get_control_chunks(&self) -> Vec<IndexChunk> {
        Self::proto_iter(&self.control_chunks)
    }

    pub(super) fn get_exp_chunks(&self) -> Vec<IndexChunk> {
        Self::proto_iter(&self.experiment_chunks)
    }
}

impl Iterator for DmrIntervalIter {
    type Item = DmrBatch;

    fn next(&mut self) -> Option<Self::Item> {
        let mut batch = DmrBatch::default();
        // let mut chunks = Self::Item::with_capacity(self.chunk_size);
        loop {
            if let Some(dmr_interval) = self.regions_of_interest.pop_front() {
                let (control_chr_id, exp_chr_id) = match (
                    self.control_contig_lookup.inner.get(&dmr_interval.chrom),
                    self.exp_contig_lookup.inner.get(&dmr_interval.chrom),
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
                match (control_chunks, exp_chunks) {
                    (Ok(control_chunks), Ok(exp_chunks)) => {
                        let dmr_chunk =
                            DmrChunk::new(control_chr_id as u32, dmr_interval);
                        batch.append(dmr_chunk, control_chunks, exp_chunks);
                    }
                    (Err(e), _) => {
                        self.failures.inc(1);
                        debug!(
                            "failed to index into {} bedMethyl for region {}, \
                             {}",
                            &self.control_fn,
                            dmr_interval,
                            e.to_string()
                        );
                        continue;
                    }
                    (_, Err(e)) => {
                        self.failures.inc(1);
                        debug!(
                            "failed to index into {} bedMethyl for chrom id \
                             {}, {}",
                            &self.exp_fn,
                            exp_chr_id,
                            e.to_string()
                        );
                        continue;
                    }
                };

                if batch.size() >= self.chunk_size {
                    break;
                } else {
                    continue;
                }
            } else {
                break;
            }
        }
        if batch.is_empty() {
            None
        } else {
            Some(batch)
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
        // todo check that regions do not overlap
        .map(|line| DmrInterval::parse_str(&line))
        .collect::<anyhow::Result<Vec<DmrInterval>>>()?;
    if intervals.is_empty() {
        bail!("didn't parse any regions")
    } else {
        Ok(intervals)
    }
}

pub(super) struct ContigLookup {
    pub(super) sample_name: Option<String>,
    pub(super) file_path: PathBuf,
    pub(super) inner: HashMap<String, usize>,
}

impl ContigLookup {
    pub(super) fn new(
        index: &CsiIndex,
        index_fp: &PathBuf,
        sample_name: Option<&String>,
    ) -> anyhow::Result<Self> {
        let inner = index
            .header()
            .ok_or_else(|| anyhow!("failed to get control tabix header"))?
            .reference_sequence_names()
            .iter()
            .enumerate()
            .map(|(idx, r)| (r.to_owned(), idx))
            .collect::<HashMap<String, usize>>();

        Ok(Self {
            sample_name: sample_name.map(|x| x.to_owned()),
            file_path: index_fp.to_owned(),
            inner,
        })
    }
}

#[cfg(test)]
mod dmr_util_tests {
    use crate::dmr::util::{parse_roi_bed, DmrInterval, ProtoIndexChunk};
    use crate::position_filter::Iv;
    use noodles::bgzf::VirtualPosition;
    use noodles::csi::index::reference_sequence::bin::Chunk as IndexChunk;
    use std::collections::BTreeSet;

    #[test]
    #[rustfmt::skip]
    fn test_parse_rois() {
        let obs = DmrInterval::parse_str(
            "chr20\t279148\t279507\tCpG: 39 359\t39\t260\t21.7\t72.4\t0.83",
        )
        .unwrap();
        let expected = DmrInterval::new(
            Iv { start: 279148, stop: 279507, val: () },
            "chr20".to_string(),
            "CpG: 39 359".to_string(),
        );
        assert_eq!(obs, expected);
        let obs = DmrInterval::parse_str(
            "chr20\t279148\t279507\tCpGby_any_other_name\t39\t260\t21.7\t72.4\t0.83",
        )
        .unwrap();
        let expected = DmrInterval::new(
            Iv { start: 279148, stop: 279507, val: () },
            "chr20".to_string(),
            "CpGby_any_other_name".to_string(),
        );
        assert_eq!(obs, expected);
        let obs = DmrInterval::parse_str("chr20\t279148\t279507\t").unwrap();
        let expected = DmrInterval::new(
            Iv { start: 279148, stop: 279507, val: () },
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
                interval: Iv { start: 10172120, stop: 10172545, val: () },
                chrom: "chr20".to_string(),
                name: "r1".to_string(),
            },
            DmrInterval {
                interval: Iv { start: 10217487, stop: 10218336, val: () },
                chrom: "chr20".to_string(),
                name: "r2".to_string(),
            },
            DmrInterval {
                interval: Iv { start: 10034963, stop: 10035266, val: () },
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
                interval: Iv { start: 10172120, stop: 10172545, val: () },
                chrom: "chr20".to_string(),
                name: "chr20:10172120-10172545".to_string(),
            },
            DmrInterval {
                interval: Iv { start: 10217487, stop: 10218336, val: () },
                chrom: "chr20".to_string(),
                name: "chr20:10217487-10218336".to_string(),
            },
            DmrInterval {
                interval: Iv { start: 10034963, stop: 10035266, val: () },
                chrom: "chr20".to_string(),
                name: "chr20:10034963-10035266".to_string(),
            },
        ]
        .to_vec();
        assert_eq!(rois, expected);
    }

    #[test]
    fn test_position_offset_in_set() {
        let start = VirtualPosition::from(18176487209);
        let end = VirtualPosition::from(18762153665);
        let chunk1: IndexChunk = IndexChunk::from(start..end);
        let chunk2: IndexChunk = IndexChunk::from(start..end);
        let mut chunks: BTreeSet<ProtoIndexChunk> = BTreeSet::new();
        let added = chunks.insert(chunk1.into());
        assert!(added);
        let added = chunks.insert(chunk2.into());
        assert!(!added);
    }
}
