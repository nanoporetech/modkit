use std::collections::{BTreeSet, HashMap, HashSet, VecDeque};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::ops::Range;
use std::path::PathBuf;

use anyhow::{anyhow, bail};
use derive_new::new;
use itertools::{Itertools, MinMaxResult};
use log::{debug, info};
use nom::character::complete::multispace1;
use nom::IResult;
use rayon::prelude::*;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{self, FetchDefinition, Read};
use rustc_hash::FxHashMap;

use crate::entropy::methylation_entropy::calc_me_entropy;
use crate::errs::RunError;
use crate::mod_bam::{BaseModCall, ModBaseInfo};
use crate::mod_base_code::ModCodeRepr;
use crate::motif_bed::RegexMotif;
use crate::read_ids_to_base_mod_probs::{PositionModCalls, ReadBaseModProfile};
use crate::threshold_mod_caller::MultipleThresholdModCaller;
use crate::thresholds::percentile_linear_interp;
use crate::util::{record_is_not_primary, ReferenceRecord, Strand};

mod methylation_entropy;
pub mod subcommand;

#[derive(Debug)]
pub(super) enum GenomeWindow {
    CombineStrands {
        interval: Range<u64>,
        neg_to_pos_positions: FxHashMap<u64, u64>,
        read_patterns: Vec<Vec<BaseModCall>>,
        position_valid_coverages: Vec<u32>,
    },
    Stranded {
        // todo instead of having pos/neg for everything, make one struct and
        // have  an optional for all of it
        pos_interval: Option<Range<u64>>,
        neg_interval: Option<Range<u64>>,
        pos_positions: Option<Vec<u64>>,
        neg_positions: Option<Vec<u64>>,
        pos_read_patterns: Vec<Vec<BaseModCall>>,
        neg_read_patterns: Vec<Vec<BaseModCall>>,
        pos_position_valid_coverages: Vec<u32>,
        neg_position_valid_coverages: Vec<u32>,
    },
}

impl GenomeWindow {
    fn new_combine_strands(
        interval: Range<u64>,
        num_positions: usize,
        neg_to_pos_positions: FxHashMap<u64, u64>,
    ) -> Self {
        let position_valid_coverages = vec![0u32; num_positions];
        Self::CombineStrands {
            interval,
            neg_to_pos_positions,
            read_patterns: Vec::new(),
            position_valid_coverages,
        }
    }

    fn new_stranded(
        pos_positions: Option<Vec<u64>>,
        neg_positions: Option<Vec<u64>>,
        num_positions: usize,
    ) -> Self {
        let pos_interval = pos_positions.as_ref().map(|positions| {
            match positions.iter().minmax() {
                MinMaxResult::MinMax(s, t) => *s..*t,
                MinMaxResult::OneElement(x) => *x..(*x + 1u64),
                MinMaxResult::NoElements => {
                    unreachable!("should have >0 elements")
                }
            }
        });
        let neg_interval = neg_positions.as_ref().map(|positions| {
            match positions.iter().minmax() {
                MinMaxResult::MinMax(s, t) => *s..*t,
                MinMaxResult::OneElement(x) => *x..(*x + 1u64),
                MinMaxResult::NoElements => {
                    unreachable!("should have >0 elements")
                }
            }
        });

        // TODO remove after testing
        let check = |positions: Option<&Vec<u64>>| {
            if let Some(ps) = positions {
                ps.iter().skip(1).fold(ps[0], |last, next| {
                    assert!(last < *next, "needs to be sorted");
                    *next
                });
            }
        };

        check(pos_positions.as_ref());
        check(neg_positions.as_ref());
        let pos_position_valid_coverages = vec![0u32; num_positions];
        let neg_position_valid_coverages = vec![0u32; num_positions];
        Self::Stranded {
            pos_interval,
            neg_interval,
            pos_positions,
            neg_positions,
            pos_read_patterns: Vec::new(),
            neg_read_patterns: Vec::new(),
            pos_position_valid_coverages,
            neg_position_valid_coverages,
        }
    }

    #[inline]
    fn inc_coverage(&mut self, pos: usize, strand: &Strand) {
        match self {
            Self::CombineStrands { position_valid_coverages, .. } => {
                assert!(
                    pos < position_valid_coverages.len(),
                    "pos is larger than the window size?"
                );
                position_valid_coverages[pos] += 1;
            }
            Self::Stranded {
                pos_position_valid_coverages,
                neg_position_valid_coverages,
                ..
            } => match strand {
                Strand::Positive => {
                    assert!(
                        pos < pos_position_valid_coverages.len(),
                        "pos is larger than the window size?"
                    );
                    pos_position_valid_coverages[pos] += 1;
                }
                Strand::Negative => {
                    assert!(
                        pos < neg_position_valid_coverages.len(),
                        "pos is larger than the window size?"
                    );
                    neg_position_valid_coverages[pos] += 1;
                }
            },
        };
    }

    fn add_pattern(&mut self, strand: &Strand, pattern: Vec<BaseModCall>) {
        match self {
            Self::Stranded { pos_read_patterns, neg_read_patterns, .. } => {
                match strand {
                    Strand::Positive => pos_read_patterns.push(pattern),
                    Strand::Negative => neg_read_patterns.push(pattern),
                }
            }
            Self::CombineStrands { read_patterns, .. } => {
                read_patterns.push(pattern);
            }
        }
    }

    fn leftmost(&self) -> u64 {
        std::cmp::min(
            self.start(&Strand::Positive).unwrap_or(0),
            self.start(&Strand::Negative).unwrap_or(0),
        )
    }

    fn rightmost(&self) -> u64 {
        std::cmp::max(
            self.end(&Strand::Positive).unwrap_or(0),
            self.end(&Strand::Negative).unwrap_or(0),
        )
    }

    fn start(&self, strand: &Strand) -> Option<u64> {
        match self {
            Self::CombineStrands { interval, .. } => Some(interval.start),
            Self::Stranded { pos_interval, neg_interval, .. } => match strand {
                Strand::Positive => pos_interval.as_ref().map(|x| x.start),
                Strand::Negative => neg_interval.as_ref().map(|x| x.start),
            },
        }
    }

    fn end(&self, strand: &Strand) -> Option<u64> {
        match self {
            Self::CombineStrands { interval, .. } => Some(interval.end),
            Self::Stranded { pos_interval, neg_interval, .. } => match strand {
                Strand::Positive => pos_interval.as_ref().map(|x| x.end),
                Strand::Negative => neg_interval.as_ref().map(|x| x.end),
            },
        }
    }

    fn add_read_to_patterns(
        &mut self,
        ref_pos_to_basemod_call: &FxHashMap<u64, BaseModCall>,
        record: &bam::Record,
        strand: &Strand,
        max_filtered_positions: usize,
    ) {
        // check that the read fully covers the interval
        let reference_start = if record.reference_start() >= 0 {
            Some(record.reference_start() as u64)
        } else {
            None
        };
        let reference_end = if reference_start
            .map(|x| record.reference_end() > x as i64)
            .unwrap_or(false)
        {
            Some(record.reference_end() as u64)
        } else {
            None
        };

        let overlaps = reference_start
            .and_then(|s| reference_end.map(|t| (s, t)))
            .map(|(s, t)| match (self.start(strand), self.end(strand)) {
                (Some(wind_start), Some(wind_end)) => {
                    s <= wind_start && t >= wind_end
                }
                _ => false,
            })
            // .map(|(s, t)| s <= self.start() && t >= self.end())
            .unwrap_or(false);
        if !overlaps {
            return;
        }

        let pattern: Vec<BaseModCall> = match strand {
            Strand::Positive => match &self {
                Self::Stranded { pos_positions: Some(positions), .. } => {
                    positions
                        .iter()
                        .map(|p| {
                            ref_pos_to_basemod_call
                                .get(p)
                                .copied()
                                .unwrap_or(BaseModCall::Filtered)
                        })
                        .collect()
                }
                Self::CombineStrands { neg_to_pos_positions, .. } => {
                    neg_to_pos_positions
                        .values()
                        .map(|p| {
                            let call = ref_pos_to_basemod_call
                                .get(p)
                                .copied()
                                .unwrap_or(BaseModCall::Filtered);
                            (p, call)
                        })
                        .sorted_by(|(a, _), (b, _)| a.cmp(b))
                        .map(|(_, call)| call)
                        .collect()
                }
                _ => return,
            },
            Strand::Negative => match &self {
                Self::Stranded { neg_positions: Some(positions), .. } => {
                    positions
                        .iter()
                        .map(|p| {
                            ref_pos_to_basemod_call
                                .get(p)
                                .copied()
                                .unwrap_or(BaseModCall::Filtered)
                        })
                        .collect()
                }
                Self::CombineStrands { neg_to_pos_positions, .. } => {
                    neg_to_pos_positions
                        .iter()
                        .map(|(neg_position, positive_position)| {
                            let call = ref_pos_to_basemod_call
                                .get(neg_position)
                                .copied()
                                .unwrap_or(BaseModCall::Filtered);
                            (positive_position, call)
                        })
                        .sorted_by(|(a, _), (b, _)| a.cmp(b))
                        .map(|(_, call)| call)
                        .collect()
                }
                _ => return,
            },
        };

        if pattern.iter().filter(|&bmc| bmc == &BaseModCall::Filtered).count()
            > max_filtered_positions
        {
            // skip when too many filtered positions
            return;
        }

        for (i, call) in pattern.iter().enumerate() {
            match call {
                BaseModCall::Filtered => {}
                _ => self.inc_coverage(i, strand),
            }
        }
        self.add_pattern(strand, pattern);
    }

    fn get_mod_code_lookup(&self) -> FxHashMap<ModCodeRepr, char> {
        let read_patterns: Box<dyn Iterator<Item = &Vec<BaseModCall>>> =
            match self {
                Self::Stranded {
                    pos_read_patterns, neg_read_patterns, ..
                } => {
                    Box::new(pos_read_patterns.iter().chain(neg_read_patterns))
                }
                Self::CombineStrands { read_patterns, .. } => {
                    Box::new(read_patterns.iter())
                }
            };
        read_patterns
            .flat_map(|pattern| {
                pattern.iter().filter_map(|call| match call {
                    BaseModCall::Modified(_, code) => Some(*code),
                    _ => None,
                })
            })
            .collect::<BTreeSet<ModCodeRepr>>()
            .into_iter()
            .enumerate()
            .map(|(id, code)| {
                // save 0 for canonical
                let id = id.saturating_add(1);
                let encoded = format!("{id}").parse::<char>().unwrap();
                (code, encoded)
            })
            .collect::<FxHashMap<ModCodeRepr, char>>()
    }

    fn encode_patterns(
        &self,
        chrom: &str,
        strand: Strand,
        patterns: &Vec<Vec<BaseModCall>>,
        mod_code_lookup: &FxHashMap<ModCodeRepr, char>,
        position_valid_coverages: &[u32],
        min_coverage: u32,
    ) -> Result<Vec<String>, RunError> {
        // todo remove these checks after testing
        assert!(
            self.start(&strand).is_some(),
            "start should be Some when encoding pattern"
        );
        assert!(
            self.end(&strand).is_some(),
            "end should be Some when encoding patttern"
        );

        if position_valid_coverages.iter().all(|x| *x >= min_coverage) {
            let encoded = patterns
                .iter()
                .map(|pat| {
                    let pattern = pat
                        .iter()
                        .map(|call| match call {
                            BaseModCall::Canonical(_) => '0',
                            BaseModCall::Modified(_, code) => {
                                *mod_code_lookup.get(code).unwrap()
                            }
                            BaseModCall::Filtered => '*',
                        })
                        .collect::<String>();
                    // todo remove after testing
                    assert_eq!(
                        pattern.len(),
                        position_valid_coverages.len(),
                        "pattern {pattern} is the wrong size? \
                         {position_valid_coverages:?}"
                    );
                    pattern
                })
                .collect();
            Ok(encoded)
        } else {
            let zero_coverage =
                position_valid_coverages.iter().all(|&cov| cov == 0);
            if zero_coverage {
                return Err(RunError::Skipped(format!("zero reads")));
            } else {
                let message = format!(
                    "window {}:{}-{} does not have enough coverage at all \
                     positions {:?}",
                    chrom,
                    self.start(&strand).expect(
                        "no start, should not be encoding patterns without \
                         start and end for strand"
                    ),
                    self.end(&strand).expect(
                        "no end, should not be encoding patterns without \
                         start and end for strand"
                    ),
                    position_valid_coverages,
                );
                return Err(RunError::Failed(message));
            }
        }
    }

    fn into_entropy(
        &self,
        chrom: &str,
        chrom_id: u32,
        min_valid_coverage: u32,
    ) -> WindowEntropy {
        let window_size = self.size();
        let constant = 1f32 / window_size as f32; // todo make this configurable

        let mod_code_lookup = self.get_mod_code_lookup();
        let positive_encoded_patterns = match &self {
            Self::CombineStrands {
                read_patterns,
                position_valid_coverages,
                ..
            } => self.encode_patterns(
                chrom,
                Strand::Positive,
                read_patterns,
                &mod_code_lookup,
                position_valid_coverages,
                min_valid_coverage,
            ),
            Self::Stranded {
                pos_read_patterns,
                pos_position_valid_coverages,
                ..
            } => self.encode_patterns(
                chrom,
                Strand::Positive,
                pos_read_patterns,
                &mod_code_lookup,
                &pos_position_valid_coverages,
                min_valid_coverage,
            ),
        };
        let negative_patterns = match &self {
            Self::Stranded {
                neg_read_patterns,
                neg_position_valid_coverages,
                ..
            } => Some(self.encode_patterns(
                chrom,
                Strand::Negative,
                neg_read_patterns,
                &mod_code_lookup,
                neg_position_valid_coverages,
                min_valid_coverage,
            )),
            _ => None,
        };

        // todo remove this after testing or make it a result/debug conditional
        if let Ok(patterns) = positive_encoded_patterns.as_ref() {
            assert!(
                patterns.iter().all(|x| x.len() == window_size),
                "patterns are the wrong size {positive_encoded_patterns:?}"
            );
        }
        if let Some(Ok(neg_patterns)) = negative_patterns.as_ref() {
            assert!(neg_patterns.iter().all(|x| x.len() == window_size));
        }

        let pos_me_entropy = positive_encoded_patterns.map(|patterns| {
            let me_entropy = calc_me_entropy(&patterns, window_size, constant);
            let num_reads = patterns.len();
            let interval = self.start(&Strand::Positive).unwrap()
                ..self.end(&Strand::Positive).unwrap();
            MethylationEntropy::new(me_entropy, num_reads, interval)
        });
        let neg_me_entropy = negative_patterns.map(|maybe_patterns| {
            maybe_patterns.map(|patterns| {
                let me_entropy =
                    calc_me_entropy(&patterns, window_size, constant);
                let num_reads = patterns.len();
                let interval = self.start(&Strand::Negative).unwrap()
                    ..self.end(&Strand::Negative).unwrap();
                MethylationEntropy::new(me_entropy, num_reads, interval)
            })
        });

        WindowEntropy::new(chrom_id, pos_me_entropy, neg_me_entropy)
    }

    #[inline]
    fn size(&self) -> usize {
        match self {
            Self::Stranded { pos_position_valid_coverages, .. } => {
                pos_position_valid_coverages.len()
            }
            Self::CombineStrands { position_valid_coverages, .. } => {
                position_valid_coverages.len()
            }
        }
    }
}

pub(super) struct GenomeWindows {
    chrom_id: u32,
    entropy_windows: Vec<GenomeWindow>,
    region_name: Option<String>,
}

pub(super) enum EntropyCalculation {
    Windows(Vec<WindowEntropy>),
    Region(RegionEntropy),
}

impl GenomeWindows {
    fn new(
        chrom_id: u32,
        entropy_windows: Vec<GenomeWindow>,
        region_name: Option<String>,
    ) -> Self {
        assert!(!entropy_windows.is_empty());
        Self { chrom_id, entropy_windows, region_name }
    }

    fn get_range(&self) -> Range<u64> {
        // these expects are checked in a few places, make them .unwrap()s
        let start = self
            .entropy_windows
            .first()
            .expect("self.entropy_windows should not be empty")
            .leftmost();
        let end = self
            .entropy_windows
            .last()
            .expect("self.entropy_windows should not be empty")
            .rightmost();
        start..end
    }

    fn get_fetch_definition(&self) -> FetchDefinition {
        let range = self.get_range();
        let start = range.start as i64;
        let end = range.end as i64;
        let chrom_id = self.chrom_id;
        FetchDefinition::Region(chrom_id as i32, start, end)
    }

    fn into_entropy_calculation(
        self,
        chrom: &str,
        chrom_id: u32,
        min_coverage: u32,
    ) -> EntropyCalculation {
        // to appease the bC we have to get the interval
        // here, but it's only used if we're summarizing a region
        let interval = self.get_range();
        let window_entropies = self
            .entropy_windows
            .par_iter()
            .map(|ew| ew.into_entropy(chrom, chrom_id, min_coverage))
            .collect::<Vec<_>>();
        let chrom_id = self.chrom_id;
        if let Some(region_name) = self.region_name {
            let mut pos_entropies = Vec::with_capacity(window_entropies.len());
            let mut pos_num_reads = Vec::with_capacity(window_entropies.len());
            let mut pos_num_fails = 0usize;
            let mut neg_entropies = Vec::with_capacity(window_entropies.len());
            let mut neg_num_reads = Vec::with_capacity(window_entropies.len());
            let mut neg_num_fails = 0usize;

            for window_entropy in window_entropies.iter() {
                match window_entropy.pos_me_entropy.as_ref() {
                    Ok(me_entropy) => {
                        pos_entropies.push(me_entropy.me_entropy);
                        pos_num_reads.push(me_entropy.num_reads);
                    }
                    Err(_e) => {
                        pos_num_fails += 1;
                    }
                }
                match window_entropy.neg_me_entropy.as_ref() {
                    Some(Ok(me_entropy)) => {
                        neg_entropies.push(me_entropy.me_entropy);
                        neg_num_reads.push(me_entropy.num_reads);
                    }
                    Some(Err(_e)) => {
                        neg_num_fails += 1;
                    }
                    // this means it was combine strands
                    None => {}
                }
            }

            let pos_entropy_stats = DescriptiveStats::new(
                &pos_entropies,
                &pos_num_reads,
                pos_num_fails,
            );
            // if neg_entropies is empty and there are no fails, we never saw
            // any negative strand me entropies
            let neg_entropy_stats = if neg_entropies.is_empty()
                && neg_num_fails == 0
            {
                assert!(
                    neg_num_reads.is_empty(),
                    "neg num reads and window entropies should both be empty"
                );
                None
            } else {
                // this will fail correctly if there are neg_entropies is empty
                // but there are fails
                Some(DescriptiveStats::new(
                    &neg_entropies,
                    &neg_num_reads,
                    neg_num_fails,
                ))
            };

            let region_entropy = RegionEntropy::new(
                chrom_id,
                interval,
                pos_entropy_stats,
                neg_entropy_stats,
                region_name,
            );
            EntropyCalculation::Region(region_entropy)
        } else {
            EntropyCalculation::Windows(window_entropies)
        }
    }
}

struct SlidingWindows {
    motif: RegexMotif,
    work_queue: VecDeque<(ReferenceRecord, Vec<char>)>,
    region_names: VecDeque<String>,
    window_size: usize,
    num_positions: usize,
    batch_size: usize,
    curr_position: usize,
    curr_contig: ReferenceRecord,
    curr_seq: Vec<char>,
    curr_region_name: Option<String>,
    combine_strands: bool,
    done: bool,
}

impl<'a> SlidingWindows {
    fn new_with_regions(
        names_to_tid: &HashMap<&str, u32>,
        reference_sequences: HashMap<&str, Vec<char>>,
        regions_bed_fp: &PathBuf,
        motif: RegexMotif,
        combine_strands: bool,
        num_positions: usize,
        window_size: usize,
        batch_size: usize,
    ) -> anyhow::Result<Self> {
        let regions_iter = BufReader::new(File::open(regions_bed_fp)?)
            .lines()
            .map(|r| r.map_err(|e| anyhow!("failed to read line, {e}")))
            .map(|r| r.and_then(|l| BedRegion::parse_str(&l)));

        let mut work_queue = VecDeque::new();
        let mut region_queue = VecDeque::new();
        let mut failures = HashMap::new();

        let mut add_failure = |cause: String| {
            *failures.entry(cause).or_insert(0) += 1;
        };

        for res in regions_iter {
            match res {
                Ok(region) => {
                    if let Some(seq) =
                        reference_sequences.get(region.chrom.as_str())
                    {
                        let length = region.length();
                        let region_name = region.name;
                        let chrom_name = region.chrom;
                        let start = region.interval.start as u32;
                        let end = region.interval.end;
                        if end >= seq.len() {
                            debug!(
                                "region {chrom_name}:{start}-{end} \
                                 ({region_name}) beyond length of contig {}",
                                seq.len()
                            );
                            add_failure(format!(
                                "region beyond length of contig"
                            ));
                            continue;
                        }

                        let sub_seq = seq[region.interval]
                            .iter()
                            .copied()
                            .collect::<Vec<char>>();
                        let tid = names_to_tid
                            .get(chrom_name.as_str())
                            .expect("should have tid");
                        let ref_record = ReferenceRecord::new(
                            *tid,
                            start,
                            length as u32,
                            chrom_name,
                        );
                        work_queue.push_back((ref_record, sub_seq));
                        region_queue.push_back(region_name);
                    } else {
                        add_failure(format!(
                            "contig {} not in BAM header",
                            &region.name
                        ));
                        continue;
                    }
                }
                Err(e) => {
                    add_failure(e.to_string());
                }
            }
        }

        if !failures.is_empty() {
            debug!("failure reasons while parsing regions BED file");
            for (cause, count) in
                failures.iter().sorted_by(|(_, a), (_, b)| a.cmp(b))
            {
                debug!("\t {cause}: {count}")
            }
        }

        if work_queue.is_empty() {
            bail!("no valid regions parsed");
        }

        assert_eq!(region_queue.len(), work_queue.len());
        let (curr_contig, curr_seq, curr_position, curr_region_name) = loop {
            let (ref_record, subseq, region_name) =
                match (work_queue.pop_front(), region_queue.pop_front()) {
                    (Some((rr, subseq)), Some(region_name)) => {
                        anyhow::Ok((rr, subseq, region_name))
                    }
                    _ => bail!(
                        "didn't find at least 1 sequence with valid start \
                         position"
                    ),
                }?;
            if let Some(start_position) =
                Self::find_start_position(&subseq, &motif)
            {
                info!(
                    "starting with region {region_name} at 0-based position \
                     {} on contig {}",
                    start_position + ref_record.start as usize,
                    &ref_record.name
                );
                break (ref_record, subseq, start_position, region_name);
            } else {
                info!("region {region_name} has no valid positions, skipping");
                continue;
            }
        };
        debug!(
            "parsed {} regions, starting with {} on contig {}",
            region_queue.len() + 1usize,
            &curr_region_name,
            curr_contig.name
        );

        Ok(Self {
            motif,
            work_queue,
            region_names: region_queue,
            window_size,
            num_positions,
            batch_size,
            curr_position,
            curr_contig,
            curr_seq,
            curr_region_name: Some(curr_region_name),
            combine_strands,
            done: false,
        })
    }

    fn new(
        bam_header_records: Vec<ReferenceRecord>,
        mut reference_sequences: HashMap<u32, Vec<char>>,
        motif: RegexMotif,
        combine_strands: bool,
        num_positions: usize,
        window_size: usize,
        batch_size: usize,
    ) -> anyhow::Result<Self> {
        let mut work_queue = bam_header_records
            .into_iter()
            .filter_map(|record| {
                match reference_sequences.remove(&record.tid) {
                    Some(seq) => Some((record, seq)),
                    None => None,
                }
            })
            .collect::<VecDeque<_>>();

        let (curr_contig, curr_seq, curr_position) = loop {
            let (curr_record, curr_seq) =
                work_queue.pop_front().ok_or_else(|| {
                    anyhow!(
                        "didn't find at least 1 sequence with a valid start \
                         position"
                    )
                })?;
            if let Some(pos) = Self::find_start_position(&curr_seq, &motif) {
                info!(
                    "starting with contig {} at 0-based position {pos}",
                    &curr_record.name
                );
                break (curr_record, curr_seq, pos);
            } else {
                info!(
                    "contig {} had no valid motif positions, skipping..",
                    curr_record.name
                );
            }
        };

        Ok(Self {
            motif,
            work_queue,
            region_names: VecDeque::new(),
            window_size,
            num_positions,
            batch_size,
            curr_position,
            curr_contig,
            curr_seq,
            curr_region_name: None,
            combine_strands,
            done: false,
        })
    }

    #[inline]
    fn take_hits_if_enough(
        &self,
        motif_hits: &[(usize, Strand)],
    ) -> Option<Vec<u64>> {
        let positions = motif_hits
            .into_iter()
            .take(self.num_positions)
            .map(|(pos, _strand)| *pos as u64)
            .sorted()
            .collect::<Vec<_>>();
        if positions.len() == self.num_positions {
            Some(positions)
        } else {
            None
        }
    }

    #[inline]
    fn enough_hits_for_window(
        &self,
        pos_hits: &[(usize, Strand)],
        neg_hits: &[(usize, Strand)],
    ) -> Option<GenomeWindow> {
        if self.combine_strands {
            let neg_to_pos = pos_hits
                .into_iter()
                .filter(|(_, strand)| strand == &Strand::Positive)
                .take(self.num_positions)
                .filter_map(|(pos_position, strand)| {
                    assert_eq!(strand, &Strand::Positive, "logic error!");
                    self.motif
                        .motif_info
                        .negative_strand_position(*pos_position as u32)
                        .map(|neg_position| {
                            (neg_position as u64, *pos_position as u64)
                        })
                })
                .collect::<FxHashMap<u64, u64>>();
            if neg_to_pos.len() < self.num_positions {
                None
            } else {
                let (start, end) = match neg_to_pos
                    .keys()
                    .chain(neg_to_pos.values())
                    .minmax()
                {
                    MinMaxResult::MinMax(s, t) => (*s, *t),
                    MinMaxResult::OneElement(x) => (*x, *x + 1), /* should probably fail here too? */
                    _ => unreachable!("there must be more than 1 element"),
                };
                let interval = start..end;
                Some(GenomeWindow::new_combine_strands(
                    interval,
                    self.num_positions,
                    neg_to_pos,
                ))
            }
        } else {
            if pos_hits.len() >= self.num_positions
                || neg_hits.len() >= self.num_positions
            {
                let pos_positions = self.take_hits_if_enough(pos_hits);
                let neg_positions = self.take_hits_if_enough(neg_hits);
                Some(GenomeWindow::new_stranded(
                    pos_positions,
                    neg_positions,
                    self.num_positions,
                ))
            } else {
                None
            }
        }
    }

    fn next_window(&mut self) -> Option<GenomeWindow> {
        while !self.at_end_of_contig() {
            // search forward for hits
            let end = std::cmp::min(
                self.curr_position.saturating_add(self.window_size),
                self.curr_seq.len(),
            );
            // todo optimize?
            let subseq = self.curr_seq[self.curr_position..end]
                .iter()
                .map(|x| *x)
                .collect::<String>();
            // N.B. the 'position' in these tuples are  _genome coordinates_!
            // this is because when we fetch reads we need to do it with the
            // proper genome coordinates. when we're using normal
            // sliding windows, the relative coordinates and the
            // genome coordinates _should_ be the same however when
            // using regions, we slice the reference genome, so the
            // relative (to the sequence) and genome coordinates will _not_ be
            // the same
            let (pos_hits, neg_hits): (
                Vec<(usize, Strand)>,
                Vec<(usize, Strand)>,
            ) = self
                .motif
                .find_hits(&subseq)
                .into_iter()
                .map(|(pos, strand)| {
                    let adjusted_position = pos
                        .saturating_add(self.curr_position)
                        .saturating_add(self.curr_contig.start as usize);
                    (adjusted_position, strand)
                })
                .partition(|(_pos, strand)| strand == &Strand::Positive);
            if let Some(entropy_window) =
                self.enough_hits_for_window(&pos_hits, &neg_hits)
            {
                let new_genome_space_position =
                    (entropy_window.leftmost() as usize).saturating_add(1usize);
                // need to re-adjust to relative coordinates instead of genome
                // coordinates
                self.curr_position = new_genome_space_position
                    .checked_sub(self.curr_contig.start as usize)
                    .expect(
                        "should be able to subtract contig start from position",
                    );

                return Some(entropy_window);
            } else {
                // not enough on (+) or (-)
                let hits = pos_hits
                    .into_iter()
                    .chain(neg_hits)
                    .map(|(p, _)| {
                        // need to re-adjust to relative coordinates instead of
                        // genome coordinates
                        p.checked_sub(self.curr_contig.start as usize)
                            .expect("should be able to re-adjust position")
                    })
                    .collect::<BTreeSet<usize>>();
                if let Some(&first) = hits.first() {
                    // at least 1
                    if self.curr_position == first {
                        match hits.iter().nth(1) {
                            Some(&second_hit) => {
                                self.curr_position = second_hit
                            }
                            None => {
                                // there was only 1
                                self.curr_position = end;
                            }
                        }
                    } else {
                        self.curr_position = first;
                    }
                } else {
                    // hits was empty, set to end
                    self.curr_position = end;
                }
                continue;
            }
        }
        None
    }

    fn find_start_position(seq: &[char], motif: &RegexMotif) -> Option<usize> {
        seq.par_chunks(10_000).find_map_first(|c| {
            let s = c.iter().collect::<String>();
            let min_pos =
                motif.find_hits(&s).into_iter().map(|(pos, _)| pos).nth(0);
            min_pos
        })
    }

    #[inline]
    fn at_end_of_contig(&self) -> bool {
        self.curr_position >= self.curr_contig.length as usize
    }

    fn update_current_contig(&mut self) {
        if let Some((record, seq, pos)) =
            self.work_queue.pop_front().and_then(|(ref_record, seq)| {
                Self::find_start_position(&seq, &self.motif)
                    .map(|pos| (ref_record, seq, pos))
            })
        {
            self.curr_contig = record;
            self.curr_position = pos;
            self.curr_seq = seq;
            let region_name = self.region_names.pop_front();
            if let Some(rn) = region_name.as_ref() {
                debug!("next region is named {rn}");
            }
            self.curr_region_name = region_name;
        } else {
            assert!(self.region_names.is_empty());
            self.done = true
        }
    }
}

impl Iterator for SlidingWindows {
    type Item = Vec<GenomeWindows>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut batch = Vec::with_capacity(self.batch_size);
        let mut windows = Vec::new();
        loop {
            // stopping conditions
            if self.done || batch.len() >= self.batch_size {
                break;
            }

            // grab the next window
            if let Some(entropy_window) = self.next_window() {
                windows.push(entropy_window);
            }

            // update conditions
            if self.at_end_of_contig() {
                // need to rotate the windows since we're moving on to another
                // contig
                let finished_windows =
                    std::mem::replace(&mut windows, Vec::new());
                let finished_region =
                    std::mem::replace(&mut self.curr_region_name, None);
                if !finished_windows.is_empty() {
                    let entropy_windows = GenomeWindows::new(
                        self.curr_contig.tid,
                        finished_windows,
                        finished_region,
                    );
                    batch.push(entropy_windows);
                }
                self.update_current_contig();
                continue;
            }

            // N.B. semantics, if the current region name is None, we're just
            // doing sliding windows over the genome and we can cut
            // this batch once the window size is the batch size.
            // otoh, if current_region is Some, we never cut the batch until
            // we've finished the contig so that an entire region ends up
            // in a single batch
            if self.curr_region_name.is_none()
                && windows.len() > self.batch_size
            {
                assert!(
                    self.region_names.is_empty(),
                    "region names should be empty here!"
                );
                let finished_windows =
                    std::mem::replace(&mut windows, Vec::new());
                if !finished_windows.is_empty() {
                    let entropy_windows = GenomeWindows::new(
                        self.curr_contig.tid,
                        finished_windows,
                        None,
                    );
                    batch.push(entropy_windows);
                }
            }
        }

        if !windows.is_empty() {
            assert!(
                self.region_names.is_empty(),
                "region names should be empty here also!"
            );
            let entropy_windows =
                GenomeWindows::new(self.curr_contig.tid, windows, None);
            batch.push(entropy_windows)
        }

        if batch.is_empty() {
            None
        } else {
            Some(batch)
        }
    }
}

#[derive(new, Debug)]
pub(super) struct MethylationEntropy {
    me_entropy: f32,
    num_reads: usize,
    interval: Range<u64>,
}

// todo make this an enum, one for regions
#[derive(new, Debug)]
pub(super) struct WindowEntropy {
    chrom_id: u32,
    pos_me_entropy: Result<MethylationEntropy, RunError>,
    neg_me_entropy: Option<Result<MethylationEntropy, RunError>>,
}

struct DescriptiveStats {
    mean_entropy: f32,
    median_entropy: f32,
    max_entropy: f32,
    min_entropy: f32,
    mean_num_reads: f32,
    max_num_reads: usize,
    min_num_reads: usize,
    failed_count: usize,
    successful_count: usize,
}

impl DescriptiveStats {
    fn mean(xs: &[f32]) -> f32 {
        xs.iter().sum::<f32>() / (xs.len() as f32)
    }

    fn new(
        measurements: &[f32],
        n_reads: &[usize],
        n_fails: usize,
    ) -> Result<Self, RunError> {
        if measurements.is_empty() {
            assert!(
                n_reads.is_empty(),
                "measurements and reads should be empty together"
            );
            Err(RunError::new_failed("all reads failed"))
        } else {
            assert_eq!(
                measurements.len(),
                n_reads.len(),
                "measurements and n_reads should be the same length"
            );
            let mean_entropy = Self::mean(measurements);
            let median_entropy = percentile_linear_interp(measurements, 0.5f32)
                .map_err(|e| RunError::new_failed(e.to_string()))?;
            // safe because of above check
            let (min_entropy, max_entropy) = match measurements.iter().minmax()
            {
                MinMaxResult::OneElement(x) => (*x, *x),
                MinMaxResult::MinMax(m, x) => (*m, *x),
                MinMaxResult::NoElements => {
                    unreachable!("checked for empty above")
                }
            };

            let mean_num_reads = Self::mean(
                &n_reads.iter().map(|&x| x as f32).collect::<Vec<_>>(),
            );
            let (min_num_reads, max_num_reads) = match n_reads.iter().minmax() {
                MinMaxResult::OneElement(x) => (*x, *x),
                MinMaxResult::MinMax(m, x) => (*m, *x),
                MinMaxResult::NoElements => {
                    unreachable!("checked for empty above")
                }
            };

            let success_count = measurements.len();

            Ok(Self {
                mean_entropy,
                median_entropy,
                max_entropy,
                min_entropy,
                mean_num_reads,
                max_num_reads,
                min_num_reads,
                successful_count: success_count,
                failed_count: n_fails,
            })
        }
    }

    pub(super) fn to_row(
        &self,
        chrom: &str,
        start: u64,
        end: u64,
        strand: Strand,
        region_name: &str,
    ) -> String {
        use crate::util::TAB;

        format!("\
            {chrom}{TAB}\
            {start}{TAB}\
            {end}{TAB}\
            {region_name}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}{TAB}\
            {}\n",
              self.mean_entropy,
              strand.to_char(),
              self.median_entropy,
              self.max_entropy,
              self.min_entropy,
              self.mean_num_reads,
              self.min_num_reads,
              self.max_num_reads,
              self.successful_count,
              self.failed_count
        )
    }
}

#[derive(new)]
pub(super) struct RegionEntropy {
    chrom_id: u32,
    interval: Range<u64>,
    pos_entropy_stats: Result<DescriptiveStats, RunError>,
    neg_entropy_stats: Option<Result<DescriptiveStats, RunError>>,
    region_name: String,
}

pub(super) fn process_entropy_window(
    mut entropy_windows: GenomeWindows,
    min_coverage: u32,
    max_filtered_positions: usize,
    io_threads: usize,
    caller: &MultipleThresholdModCaller,
    bam_fp: &PathBuf,
) -> anyhow::Result<EntropyCalculation> {
    let mut reader = bam::IndexedReader::from_path(bam_fp)?;
    reader.set_threads(io_threads)?;
    let header = reader.header();
    let chrom_id = entropy_windows.chrom_id;
    let chrom = String::from_utf8_lossy(header.tid2name(chrom_id)).to_string();
    let fd = entropy_windows.get_fetch_definition();
    reader.fetch(fd)?;

    let record_iter = reader
        .records()
        .filter_map(|r| r.ok())
        .filter(|record| {
            !record.is_unmapped()
                && !(record_is_not_primary(&record) || record.seq_len() == 0)
        })
        .filter_map(|record| {
            String::from_utf8(record.qname().to_vec())
                .ok()
                .map(|name| (record, name))
        })
        .filter_map(|(record, name)| {
            match ModBaseInfo::new_from_record(&record) {
                Ok(modbase_info) => Some((modbase_info, record, name)),
                Err(run_error) => {
                    debug!(
                        "read {name}, failed to parse modbase info, \
                         {run_error}"
                    );
                    None
                }
            }
        })
        .filter_map(|(modbase_info, record, name)| {
            match ReadBaseModProfile::process_record(
                &record,
                &name,
                modbase_info,
                None,
                None,
                1,
            ) {
                Ok(profile) => {
                    let position_calls =
                        PositionModCalls::from_profile(&profile);
                    let strands = position_calls
                        .iter()
                        .map(|p| p.mod_strand)
                        .collect::<HashSet<Strand>>();
                    if strands.len() > 1 {
                        debug!("duplex not yet supported");
                        None
                    } else {
                        let strand = if record.is_reverse() {
                            Strand::Negative
                        } else {
                            Strand::Positive
                        };
                        let mod_calls = position_calls
                            .into_iter()
                            .filter_map(|p| {
                                match (p.ref_position, p.alignment_strand) {
                                    (Some(ref_pos), Some(aln_strand)) => {
                                        Some((p, ref_pos, aln_strand))
                                    }
                                    _ => None,
                                }
                            })
                            .map(|(p, ref_pos, _alignment_strand)| {
                                let mod_base_call = caller
                                    .call(&p.canonical_base, &p.base_mod_probs);
                                (ref_pos as u64, mod_base_call)
                            })
                            .collect::<FxHashMap<u64, BaseModCall>>();
                        Some((mod_calls, strand, record, name))
                    }
                }
                Err(e) => {
                    debug!("read {name} failed to extract modbase info, {e}");
                    None
                }
            }
        });

    for (mod_calls, strand, record, _name) in record_iter {
        entropy_windows.entropy_windows.par_iter_mut().for_each(|window| {
            window.add_read_to_patterns(
                &mod_calls,
                &record,
                &strand,
                max_filtered_positions,
            )
        });
    }

    Ok(entropy_windows.into_entropy_calculation(&chrom, chrom_id, min_coverage))
}

#[derive(new, Debug)]
struct BedRegion {
    chrom: String,
    interval: Range<usize>,
    name: String,
}

impl BedRegion {
    fn length(&self) -> usize {
        self.interval.end - self.interval.start
    }

    fn parser(raw: &str) -> IResult<&str, Self> {
        let n_parts = raw.split('\t').count();
        let (rest, chrom) = crate::parsing_utils::consume_string(raw)?;
        let (rest, start) = crate::parsing_utils::consume_digit(rest)?;
        let (rest, stop) = crate::parsing_utils::consume_digit(rest)?;
        let (rest, name) = if n_parts == 3 {
            (rest, format!("{chrom}:{start}-{stop}"))
        } else {
            let (rest, _) = multispace1(rest)?;
            crate::parsing_utils::consume_string(rest)?
        };

        let interval = (start as usize)..(stop as usize);
        let this = Self { chrom, interval, name };
        Ok((rest, this))
    }

    fn parse_str(raw: &str) -> anyhow::Result<Self> {
        Self::parser(raw)
            .map_err(|e| anyhow!("failed to parse {raw} into BED3 line, {e}"))
            .and_then(|(_, this)| {
                if this.interval.end > this.interval.start {
                    Ok(this)
                } else {
                    bail!("end must be after start")
                }
            })
    }
}

#[cfg(test)]
mod entropy_mod_tests {
    use crate::entropy::BedRegion;

    #[test]
    fn test_bed_region_parsing() {
        let raw = "chr1\t100\t101\tfoo\n";
        let bed_region = BedRegion::parse_str(raw).expect("should parse");
        assert_eq!(&bed_region.chrom, "chr1");
        assert_eq!(bed_region.interval, 100usize..101);
        assert_eq!(&bed_region.name, "foo");
        let raw = "chr1\t100\t101\tfoo\t400\t.\tmorestuff\n";
        let bed_region = BedRegion::parse_str(raw).expect("should parse");
        assert_eq!(&bed_region.chrom, "chr1");
        assert_eq!(bed_region.interval, 100usize..101);
        assert_eq!(&bed_region.name, "foo");
    }
}
