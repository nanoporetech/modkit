use std::collections::{BTreeSet, HashMap, HashSet, VecDeque};
use std::ops::Range;
use std::path::PathBuf;

use anyhow::anyhow;
use derive_new::new;
use itertools::{Itertools, MinMaxResult};
use log::{debug, info};
use rayon::prelude::*;
use rust_htslib::bam::{self, FetchDefinition, Read};
use rust_htslib::bam::ext::BamRecordExtensions;
use rustc_hash::FxHashMap;

use crate::entropy::methylation_entropy::calc_me_entropy;
use crate::errs::RunError;
use crate::mod_bam::{BaseModCall, ModBaseInfo};
use crate::mod_base_code::ModCodeRepr;
use crate::motif_bed::RegexMotif;
use crate::read_ids_to_base_mod_probs::{PositionModCalls, ReadBaseModProfile};
use crate::threshold_mod_caller::MultipleThresholdModCaller;
use crate::util::{record_is_not_primary, ReferenceRecord, Strand};

mod methylation_entropy;
pub mod subcommand;

pub(super) struct EntropyWindow {
    interval: Range<u64>,
    // provides way to combine (-)-strand mappings to their (+)-strand
    // positions
    neg_to_pos_positions: FxHashMap<u64, u64>,
    // every read's contribution to this window
    // todo(n.b.) mapped back to the positive position
    read_patterns: Vec<Vec<BaseModCall>>,
    position_valid_coverages: Vec<u32>,
}

impl EntropyWindow {
    fn new(
        interval: Range<u64>,
        num_positions: usize,
        neg_to_pos_positions: FxHashMap<u64, u64>,
    ) -> Self {
        let position_valid_coverages = vec![0u32; num_positions];
        assert_eq!(neg_to_pos_positions.len(), num_positions);
        Self {
            interval,
            neg_to_pos_positions,
            position_valid_coverages,
            read_patterns: Vec::new(),
        }
    }

    fn add_read_to_patterns(
        &mut self,
        ref_pos_to_basemod_call: &FxHashMap<u64, BaseModCall>,
        record: &bam::Record,
        strand: &Strand,
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
            .map(|(s, t)| s <= self.interval.start && t >= self.interval.end)
            .unwrap_or(false);
        if !overlaps {
            return;
        }

        let pattern_iter: Box<dyn Iterator<Item = (u64, BaseModCall)>> =
            match strand {
                Strand::Negative => {
                    Box::new(self.neg_to_pos_positions.iter().map(
                        |(neg_pos, pos_pos)| {
                            let bmc = ref_pos_to_basemod_call
                                .get(neg_pos)
                                .copied()
                                .unwrap_or(BaseModCall::Filtered);
                            (*pos_pos, bmc)
                        },
                    ))
                }
                Strand::Positive => Box::new(
                    self.neg_to_pos_positions.values().map(|pos_position| {
                        let bmc = ref_pos_to_basemod_call
                            .get(pos_position)
                            .copied()
                            .unwrap_or(BaseModCall::Filtered);
                        (*pos_position, bmc)
                    }),
                ),
            };

        let pattern = pattern_iter
            .sorted_by(|(x, _), (y, _)| x.cmp(y))
            .map(|(_, c)| c)
            .collect::<Vec<BaseModCall>>();

        for (i, call) in pattern.iter().enumerate() {
            assert!(i < self.position_valid_coverages.len());
            match call {
                BaseModCall::Filtered => {}
                _ => self.position_valid_coverages[i] += 1u32,
            }
        }
        self.read_patterns.push(pattern);
    }

    fn into_entropy(
        self,
        chrom: &str,
        chrom_id: u32,
        min_valid_coverage: u32,
    ) -> Result<Entropy, RunError> {
        if !self
            .position_valid_coverages
            .iter()
            .all(|&cov| cov > min_valid_coverage)
        {
            let zero_coverage =
                self.position_valid_coverages.iter().all(|&cov| cov == 0);
            if zero_coverage {
                return Err(RunError::Skipped(format!("zero reads")));
            } else {
                let message = format!(
                    "window {}:{}-{} does not have enough coverage at all \
                     positions {:?}",
                    chrom,
                    self.interval.start,
                    self.interval.end,
                    self.position_valid_coverages,
                );
                return Err(RunError::Failed(message));
            }
        }

        let window_size = self.size();
        let mod_code_lookup = self
            .read_patterns
            .iter()
            .flat_map(|pat| {
                pat.iter().filter_map(|call| match call {
                    BaseModCall::Modified(_, code) => Some(*code),
                    _ => None,
                })
            })
            .collect::<BTreeSet<ModCodeRepr>>()
            .into_iter()
            .enumerate()
            .map(|(id, code)| {
                let id = id.saturating_add(1); // save 0 for canonical
                let encoded = format!("{id}").parse::<char>().unwrap();
                (code, encoded)
            })
            .collect::<FxHashMap<ModCodeRepr, char>>();

        let patterns = self
            .read_patterns
            .into_iter()
            .map(|pat| {
                pat.into_iter()
                    .map(|c| match c {
                        BaseModCall::Canonical(_) => '0',
                        BaseModCall::Modified(_, code) => {
                            *mod_code_lookup.get(&code).unwrap()
                        }
                        BaseModCall::Filtered => '*',
                    })
                    .collect::<String>()
            })
            .collect::<Vec<String>>();

        assert!(patterns.iter().all(|x| x.len() == window_size));

        let constant = 1f32 / window_size as f32;
        let me_entropy = calc_me_entropy(&patterns, window_size, constant);
        let num_reads = patterns.len();
        let entropy =
            Entropy::new(chrom_id, self.interval, me_entropy, num_reads);
        Ok(entropy)
    }

    #[inline]
    fn size(&self) -> usize {
        self.position_valid_coverages.len()
    }
}

pub(super) struct EntropyWindows {
    chrom_id: u32,
    entropy_windows: Vec<EntropyWindow>,
}

impl EntropyWindows {
    fn new(chrom_id: u32, entropy_windows: Vec<EntropyWindow>) -> Self {
        assert!(!entropy_windows.is_empty());
        Self { chrom_id, entropy_windows }
    }

    fn get_fetch_definition(&self) -> bam::FetchDefinition {
        let interval = &self
            .entropy_windows
            .first()
            .expect("entropy windows should not be empty")
            .interval;
        let start = interval.start as i64;
        let end = interval.end as i64;
        let chrom_id = self.chrom_id;
        FetchDefinition::Region(chrom_id as i32, start, end)
    }
}

struct SlidingWindows {
    motif: RegexMotif,
    work_queue: VecDeque<(ReferenceRecord, Vec<char>)>,
    window_size: usize, // need to change here for "unbounded windows"
    num_positions: usize,
    batch_size: usize,
    curr_position: usize,
    curr_contig: ReferenceRecord,
    curr_seq: Vec<char>,
    done: bool,
}

impl<'a> SlidingWindows {
    fn new(
        reference_records: Vec<ReferenceRecord>,
        mut reference_sequences: HashMap<u32, Vec<char>>,
        motif: RegexMotif,
        num_positions: usize,
        window_size: usize,
        batch_size: usize,
    ) -> anyhow::Result<Self> {
        let mut work_queue = reference_records
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
            window_size,
            num_positions,
            batch_size,
            curr_position,
            curr_contig,
            curr_seq,
            done: false,
        })
    }

    fn next_window(&mut self) -> Option<EntropyWindow> {
        while !self.at_end_of_contig() {
            // search forward for hits
            let end = std::cmp::min(
                self.curr_position.saturating_add(self.window_size),
                self.curr_contig.length as usize,
            );
            // todo optimize?
            let subseq = self.curr_seq[self.curr_position..end]
                .iter()
                .map(|x| *x)
                .collect::<String>();
            let hits = self
                .motif
                .find_hits(&subseq)
                .into_iter()
                .filter_map(|(pos, strand)| {
                    if strand == Strand::Positive {
                        Some(pos.saturating_add(self.curr_position))
                    } else {
                        None
                    }
                })
                .collect::<Vec<usize>>();
            if hits.len() >= self.num_positions {
                // we have hits and there are enough of them
                let neg_to_pos = hits
                    .into_iter()
                    .take(self.num_positions)
                    .filter_map(|pos_position| {
                        self.motif
                            .motif_info
                            .negative_strand_position(pos_position as u32)
                            .map(|neg_position| {
                                (neg_position as u64, pos_position as u64)
                            })
                    })
                    .collect::<FxHashMap<u64, u64>>();
                let (start, end) =
                    match neg_to_pos.keys().chain(neg_to_pos.values()).minmax()
                    {
                        MinMaxResult::MinMax(s, t) => (*s, *t),
                        MinMaxResult::OneElement(x) => (*x, *x + 1),
                        _ => unreachable!("there must be more than 1 element"),
                    };
                let interval = start..end;
                self.curr_position = (start as usize).saturating_add(1usize);
                let entropy_window = EntropyWindow::new(
                    interval,
                    self.num_positions,
                    neg_to_pos,
                );
                return Some(entropy_window);
            } else {
                if let Some(&first) = hits.first() {
                    // at least 1
                    if self.curr_position == first {
                        match hits.iter().nth(1) {
                            Some(&second_hit) => self.curr_position = second_hit,
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
        } else {
            self.done = true
        }
    }
}

impl Iterator for SlidingWindows {
    type Item = Vec<EntropyWindows>;

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
                if !finished_windows.is_empty() {
                    let entropy_windows = EntropyWindows::new(
                        self.curr_contig.tid,
                        finished_windows,
                    );
                    batch.push(entropy_windows);
                }
                self.update_current_contig();
                continue;
            }

            // todo batch size maybe not the correct size to use here
            if windows.len() > self.batch_size {
                let finished_windows =
                    std::mem::replace(&mut windows, Vec::new());
                if !finished_windows.is_empty() {
                    let entropy_windows = EntropyWindows::new(
                        self.curr_contig.tid,
                        finished_windows,
                    );
                    batch.push(entropy_windows);
                }
            }
        }

        if !windows.is_empty() {
            let entropy_windows =
                EntropyWindows::new(self.curr_contig.tid, windows);
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
pub(super) struct Entropy {
    chrom_id: u32,
    interval: Range<u64>,
    me_entropy: f32,
    num_reads: usize,
}

pub(super) fn process_entropy_window(
    mut entropy_windows: EntropyWindows,
    min_coverage: u32,
    io_threads: usize,
    caller: &MultipleThresholdModCaller,
    bam_fp: &PathBuf,
) -> anyhow::Result<Vec<Result<Entropy, RunError>>> {
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

    for (mod_calls, strand, record, name) in record_iter {
        entropy_windows.entropy_windows.par_iter_mut().for_each(|window| {
            window.add_read_to_patterns(&mod_calls, &record, &strand)
        });
    }
    let entropy_calcs = entropy_windows
        .entropy_windows
        .into_par_iter()
        .map(|ew| ew.into_entropy(&chrom, chrom_id, min_coverage))
        .collect::<Vec<_>>();
    Ok(entropy_calcs)
}
