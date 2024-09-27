use std::cmp::Ordering;
use std::collections::VecDeque;
use std::fmt::{Debug, Display, Formatter};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::sync::Arc;

use anyhow::bail;
use clap::ValueEnum;
use derive_new::new;
use log::{debug, error};
use log_once::debug_once;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::dmr::tabix::MultiSampleIndex;
use crate::genome_positions::{GenomePositions, StrandedPosition};
use crate::mod_base_code::DnaBase;
use crate::position_filter::Iv;
use crate::util::{GenomeRegion, StrandRule};

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

// todo rename to ROI
#[derive(new, Clone, Debug, Eq, PartialEq)]
pub(super) struct DmrInterval {
    // todo refacter out Iv and lapper-things
    pub(super) interval: Iv,
    pub(super) chrom: String,
    pub(super) name: String,
    pub(super) strand: StrandRule,
}

impl DmrInterval {
    pub(super) fn parse_unstranded_bed_line(
        line: &str,
    ) -> anyhow::Result<Self> {
        let genome_region = GenomeRegion::parse_unstranded_bed_line(line)?;
        let name = genome_region.name.unwrap_or(format!(
            "{}:{}-{}",
            genome_region.chrom, genome_region.start, genome_region.end
        ));
        Ok(Self {
            interval: Iv {
                start: genome_region.start,
                stop: genome_region.end,
                val: (),
            },
            chrom: genome_region.chrom,
            name,
            strand: genome_region.strand,
        })
    }

    pub(super) fn parse_stranded_bed_line(line: &str) -> anyhow::Result<Self> {
        let genome_region = GenomeRegion::parse_stranded_bed_line(line)?;
        let name = genome_region.name.unwrap_or(format!(
            "{}:{}-{}",
            genome_region.chrom, genome_region.start, genome_region.end
        ));
        Ok(Self {
            interval: Iv {
                start: genome_region.start,
                stop: genome_region.end,
                val: (),
            },
            chrom: genome_region.chrom,
            name,
            strand: genome_region.strand,
        })
    }

    pub(super) fn start(&self) -> u64 {
        self.interval.start
    }

    pub(super) fn stop(&self) -> u64 {
        self.interval.stop
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

#[derive(Debug, Eq, PartialEq)]
pub(super) struct RegionOfInterest {
    // pub(super) chrom_id: u32,
    pub(super) positions: FxHashSet<StrandedPosition<DnaBase>>,
    pub(super) dmr_interval: DmrInterval,
}

impl RegionOfInterest {
    /// Create a new region of interest with the positions selected
    /// for just that region. None if there aren't any positions within
    /// that region.
    pub(super) fn new_from_interval(
        dmr_interval: &DmrInterval,
        genome_positions: &GenomePositions,
    ) -> Option<Self> {
        genome_positions
            .get_positions(
                &dmr_interval.chrom,
                &(dmr_interval.start()..dmr_interval.stop()),
                dmr_interval.strand,
            )
            .map(|positions| Self {
                positions: positions.into_iter().collect(),
                dmr_interval: dmr_interval.to_owned(),
            })
    }
}

impl PartialOrd for RegionOfInterest {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.dmr_interval.cmp(&other.dmr_interval))
    }
}

impl Ord for RegionOfInterest {
    fn cmp(&self, other: &Self) -> Ordering {
        self.dmr_interval.cmp(&other.dmr_interval)
    }
}

impl Display for RegionOfInterest {
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

pub(super) struct RoiIter {
    sample_index_a: Vec<usize>,
    sample_index_b: Vec<usize>,
    regions_of_interest: VecDeque<DmrInterval>,
    chunk_size: usize,
    genome_positions: Arc<GenomePositions>,
}

impl RoiIter {
    pub(super) fn new(
        sample_a_idxs: &[usize],
        sample_b_idxs: &[usize],
        sample_name_a: &str,
        sample_name_b: &str,
        sample_index: Arc<MultiSampleIndex>,
        rois: Vec<DmrInterval>,
        chunk_size: usize,
        handle_missing: HandleMissing,
        genome_positions: Arc<GenomePositions>,
    ) -> anyhow::Result<Self> {
        // there is a lot of lines below, but, this is really just a bunch of
        // input checking, depending on the "handle_missing" enum,
        // we might warn, fail, or do nothing
        let regions_of_interest =
            rois.into_iter().try_fold(Vec::new(), |mut acc, roi| {
                let a_found = sample_a_idxs
                    .iter()
                    .any(|i| sample_index.has_contig(*i, &roi.chrom));
                let b_found = sample_b_idxs
                    .iter()
                    .any(|i| sample_index.has_contig(*i, &roi.chrom));
                if a_found && b_found {
                    // happy path
                    acc.push(roi);
                    Ok(acc)
                } else {
                    // missing from one or the other..
                    let which = if b_found {
                        sample_name_a.to_string()
                    } else {
                        format!("'{sample_name_a}' and '{sample_name_b}'")
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
            sample_index_a: sample_a_idxs.to_vec(),
            sample_index_b: sample_b_idxs.to_vec(),
            regions_of_interest: regions_of_interest.into_iter().collect(),
            chunk_size,
            genome_positions,
        })
    }
}

#[derive(Default)]
pub(super) struct DmrBatch<T: Default> {
    // these contain the filtering info
    pub(super) dmr_chunks: T,
    // mapping of the sample index (order provided on the command line)
    // to the set of chunks to use
    pub(super) idxs_a: FxHashSet<usize>,
    pub(super) idxs_b: FxHashSet<usize>,
    pub(super) regions: FxHashMap<String, std::ops::Range<u64>>,
}

/// A batch of single positions, mapping of chrom to stranded position
pub type DmrBatchOfPositions =
    DmrBatch<FxHashMap<String, FxHashSet<StrandedPosition<DnaBase>>>>;
impl DmrBatchOfPositions {
    pub(super) fn num_chunks(&self) -> usize {
        self.dmr_chunks.values().map(|positions| positions.len()).sum::<usize>()
    }

    pub(super) fn add_chunks(
        &mut self,
        contig: &str,
        range: std::ops::Range<u64>,
        positions: Vec<StrandedPosition<DnaBase>>,
        control_idxs: &[usize],
        experiment_idxs: &[usize],
    ) {
        let contig_positions = self
            .dmr_chunks
            .entry(contig.to_owned())
            .or_insert_with(|| FxHashSet::default());
        for position in positions {
            contig_positions.insert(position);
        }
        for i in control_idxs {
            self.idxs_a.insert(*i);
        }
        for i in experiment_idxs {
            self.idxs_b.insert(*i);
        }
        if let Some(r) = self.regions.get_mut(contig) {
            r.start = std::cmp::min(r.start, range.start);
            r.end = std::cmp::max(r.end, range.end);
        } else {
            self.regions.insert(contig.to_owned(), range);
        }
    }

    pub(super) fn contains_position(
        &self,
        chrom: &str,
        position: &StrandedPosition<DnaBase>,
    ) -> bool {
        self.dmr_chunks
            .get(chrom)
            .map(|positions| positions.contains(position))
            .unwrap_or(false)
    }
}

impl DmrBatch<Vec<RegionOfInterest>> {
    fn size(&self) -> usize {
        self.dmr_chunks.len()
    }

    fn is_empty(&self) -> bool {
        self.dmr_chunks.is_empty()
    }

    fn add_chunks(
        &mut self,
        region_of_interest: RegionOfInterest,
        index_a: &[usize],
        index_b: &[usize],
    ) {
        let range = region_of_interest.dmr_interval.start()
            ..region_of_interest.dmr_interval.stop();
        let contig = &region_of_interest.dmr_interval.chrom;

        if let Some(r) = self.regions.get_mut(contig) {
            r.start = std::cmp::min(r.start, range.start);
            r.end = std::cmp::max(r.end, range.end);
        } else {
            self.regions.insert(contig.to_owned(), range);
        }

        self.dmr_chunks.push(region_of_interest);
        for i in index_a {
            self.idxs_a.insert(*i);
        }
        for i in index_b {
            self.idxs_b.insert(*i);
        }
    }
}

impl Iterator for RoiIter {
    type Item = DmrBatch<Vec<RegionOfInterest>>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut batch = DmrBatch::<Vec<RegionOfInterest>>::default();
        loop {
            if let Some(dmr_interval) = self.regions_of_interest.pop_front() {
                let region_of_interest = if let Some(roi) =
                    RegionOfInterest::new_from_interval(
                        &dmr_interval,
                        &self.genome_positions,
                    ) {
                    roi
                } else {
                    debug!(
                        "interval {dmr_interval} has zero comparative \
                         positions, skipping"
                    );
                    continue;
                };
                batch.add_chunks(
                    region_of_interest,
                    &self.sample_index_a,
                    &self.sample_index_b,
                );

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
    let mut reader = BufReader::new(File::open(fp)?)
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
        .skip_while(|l| l.starts_with('#'))
        .peekable();

    let num_fields = match reader.peek() {
        Some(l) => l.split('\t').count(),
        None => bail!("zero non-comment lines in regions"),
    };

    let parser = if num_fields <= 4 {
        |l: &str| DmrInterval::parse_unstranded_bed_line(l)
    } else {
        |l: &str| DmrInterval::parse_stranded_bed_line(l)
    };

    let intervals = reader
        // todo check that regions do not overlap
        .map(|line| parser(&line))
        .collect::<anyhow::Result<Vec<DmrInterval>>>()?;

    if intervals.is_empty() {
        bail!("didn't parse any regions")
    } else {
        Ok(intervals)
    }
}

pub(crate) fn n_choose_2(n: usize) -> anyhow::Result<usize> {
    match n {
        0 | 1 => bail!("n must be >= 2"),
        2 => Ok(1),
        _ => {
            let numerator = n * (n - 1usize);
            Ok(numerator / 2usize)
        }
    }
}

#[cfg(test)]
mod dmr_util_tests {
    use crate::dmr::util::{parse_roi_bed, DmrInterval};
    use crate::position_filter::Iv;
    use crate::util::StrandRule;

    #[test]
    #[rustfmt::skip]
    fn test_parse_rois() {
        let obs = DmrInterval::parse_stranded_bed_line(
            "chr20\t279148\t279507\tCpG: 39 359\t39\t+",
        )
        .unwrap();
        let expected = DmrInterval::new(
            Iv { start: 279148, stop: 279507, val: () },
            "chr20".to_string(),
            "CpG: 39 359".to_string(),
            StrandRule::Positive,
        );
        assert_eq!(obs, expected);
        let obs = DmrInterval::parse_stranded_bed_line(
            "chr20\t279148\t279507\tCpG: 39 359\t39\t.",
        )
            .unwrap();
        let expected = DmrInterval::new(
            Iv { start: 279148, stop: 279507, val: () },
            "chr20".to_string(),
            "CpG: 39 359".to_string(),
            StrandRule::Both,
        );
        assert_eq!(obs, expected);
        let obs = DmrInterval::parse_stranded_bed_line(
            "chr20\t279148\t279507\tCpGby_any_other_name\t39\t-\t",
        )
        .unwrap();
        let expected = DmrInterval::new(
            Iv { start: 279148, stop: 279507, val: () },
            "chr20".to_string(),
            "CpGby_any_other_name".to_string(),
            StrandRule::Negative,
        );
        assert_eq!(obs, expected);
        let obs = DmrInterval::parse_unstranded_bed_line("chr20\t279148\t279507\t").unwrap();
        let expected = DmrInterval::new(
            Iv { start: 279148, stop: 279507, val: () },
            "chr20".to_string(),
            "chr20:279148-279507".to_string(),
            StrandRule::Both,
        );
        assert_eq!(obs, expected);
        let obs = DmrInterval::parse_unstranded_bed_line("chr20\t279148\t279507 ").unwrap();
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
                strand: StrandRule::Both,
            },
            DmrInterval {
                interval: Iv { start: 10217487, stop: 10218336, val: () },
                chrom: "chr20".to_string(),
                name: "r2".to_string(),
                strand: StrandRule::Both,
            },
            DmrInterval {
                interval: Iv { start: 10034963, stop: 10035266, val: () },
                chrom: "chr20".to_string(),
                name: "r3".to_string(),
                strand: StrandRule::Both,
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
                strand: StrandRule::Both,
            },
            DmrInterval {
                interval: Iv { start: 10217487, stop: 10218336, val: () },
                chrom: "chr20".to_string(),
                name: "chr20:10217487-10218336".to_string(),
                strand: StrandRule::Both,
            },
            DmrInterval {
                interval: Iv { start: 10034963, stop: 10035266, val: () },
                chrom: "chr20".to_string(),
                name: "chr20:10034963-10035266".to_string(),
                strand: StrandRule::Both,
            },
        ]
        .to_vec();
        assert_eq!(rois, expected);
    }

    #[test]
    fn test_roi_parsing_motif_bed() {
        let fp = "tests/resources/test_motif_bed_drach.bed";
        let rois = parse_roi_bed(fp).unwrap();
        assert_eq!(rois.len(), 10);
    }
}
