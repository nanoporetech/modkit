use std::cmp::Ordering;
use std::collections::HashMap;
use std::path::Path;

use anyhow::bail;
use derive_new::new;
use log::debug;
use rust_htslib::bam::{self, FetchDefinition, Read};
use rustc_hash::FxHashMap;

use crate::interval_chunks::{FocusPositions, MultiChromCoordinates};
use crate::mod_bam::{DuplexModCall, DuplexPattern, EdgeFilter};
use crate::pileup::{get_forward_read_base, PileupIter, PileupNumericOptions};
use crate::read_cache::DuplexReadCache;
use crate::threshold_mod_caller::MultipleThresholdModCaller;
use crate::util::record_is_not_primary;

/// Summarizes the duplex (hemi) methylation patterns for
/// a genomic interval
#[derive(Debug)]
pub struct DuplexModBasePileup {
    /// name of the contig/chrom
    pub chrom_name: String,
    /// duplex pattern counts per genomic position
    pub pileup_counts: FxHashMap<u32, DuplexPileupFeatureCounts>,
    /// number of records used
    pub processed_records: usize,
    /// number of records skipped
    pub skipped_records: usize,
}

#[derive(new, Debug, Eq, PartialEq)]
pub struct DuplexPatternCounts {
    pattern: DuplexPattern,
    pub count: usize,
    pub n_other_pattern: usize,
    pub n_diff: usize,
    pub n_canonical: usize,
    pub n_fail: usize,
    pub n_nocall: usize,
}

impl DuplexPatternCounts {
    pub fn frac_pattern(&self) -> f32 {
        self.count as f32 / self.valid_coverage() as f32
    }

    pub fn valid_coverage(&self) -> usize {
        self.count + self.n_other_pattern
    }
}

impl DuplexPatternCounts {
    pub fn pattern_string(&self, primary_base: char) -> String {
        format!("{},{},{}", self.pattern[0], self.pattern[1], primary_base)
    }
}

impl Ord for DuplexPatternCounts {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.pattern[0].cmp(&other.pattern[0]) {
            Ordering::Equal => self.pattern[1].cmp(&other.pattern[1]),
            o @ _ => o,
        }
    }
}

impl PartialOrd for DuplexPatternCounts {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Summarizes the pileup for a single genomic position
#[derive(new, Debug)]
pub struct DuplexPileupFeatureCounts {
    /// pattern counts, per read base (so there will be a vec
    /// for A and C, for example).
    pub pattern_counts: HashMap<char, Vec<DuplexPatternCounts>>,
    /// Number of reads with DEL at this position
    pub n_delete: usize,
}

enum DuplexFeature {
    ModCall(DuplexModCall),
    Delete,
}

#[derive(Default, Debug)]
struct DuplexFeatureVector {
    duplex_mod_counts: FxHashMap<DuplexModCall, usize>,
    delete_counts: usize,
}

impl DuplexFeatureVector {
    fn add_feature(
        &mut self,
        feature: DuplexFeature,
        pileup_options: &PileupNumericOptions,
    ) {
        match feature {
            DuplexFeature::Delete => self.delete_counts += 1,
            DuplexFeature::ModCall(duplex_mod_call) => match pileup_options {
                PileupNumericOptions::Collapse(_)
                | PileupNumericOptions::Passthrough => {
                    *self
                        .duplex_mod_counts
                        .entry(duplex_mod_call)
                        .or_insert(0) += 1;
                }
                PileupNumericOptions::Combine => {
                    *self
                        .duplex_mod_counts
                        .entry(duplex_mod_call.into_combined())
                        .or_insert(0) += 1;
                }
            },
        }
    }

    fn decode(self) -> DuplexPileupFeatureCounts {
        let grouped_by_base = self.duplex_mod_counts.into_iter().fold(
            HashMap::new(),
            |mut acc, (mod_call, count)| {
                let primary_base = mod_call.primary_base();
                acc.entry(primary_base)
                    .or_insert(Vec::new())
                    .push((mod_call, count));
                acc
            },
        );

        let mut agg_pattern_counts = HashMap::new();
        for (primary_base, duplex_calls) in &grouped_by_base {
            let mut agg = Vec::new();
            let pattern_counts = duplex_calls
                .iter()
                .filter_map(|(x, c)| x.pattern().map(|p| (p, *c)))
                .collect::<HashMap<DuplexPattern, usize>>();

            let n_diff = grouped_by_base
                .iter()
                .filter(|(base, _)| *base != primary_base)
                .map(|(_, calls)| {
                    calls
                        .iter()
                        .filter_map(|(x, c)| {
                            if x.is_mod_call() || x.is_canonical() {
                                Some(*c)
                            } else {
                                None
                            }
                        })
                        .sum::<usize>()
                })
                .sum::<usize>();

            let n_canonical = duplex_calls
                .iter()
                .filter_map(
                    |(x, c)| if x.is_canonical() { Some(*c) } else { None },
                )
                .sum::<usize>();

            let n_filtered = duplex_calls
                .iter()
                .filter_map(
                    |(x, c)| if x.is_filtered() { Some(*c) } else { None },
                )
                .sum::<usize>();

            let n_nocall = duplex_calls
                .iter()
                .filter_map(
                    |(x, c)| if x.is_nocall() { Some(*c) } else { None },
                )
                .sum::<usize>();

            for (pattern, count) in duplex_calls
                .iter()
                .filter_map(|(x, count)| x.pattern().map(|p| (p, count)))
            {
                let n_other_pattern = pattern_counts
                    .iter()
                    .filter(|(p, _c)| *p != &pattern)
                    .map(|(_, c)| *c)
                    .sum::<usize>();
                let duplex_pattern_counts = DuplexPatternCounts::new(
                    pattern,
                    *count,
                    n_other_pattern,
                    n_diff,
                    n_canonical,
                    n_filtered,
                    n_nocall,
                );
                agg.push(duplex_pattern_counts);
            }

            agg_pattern_counts.insert(*primary_base, agg);
        }

        DuplexPileupFeatureCounts::new(agg_pattern_counts, self.delete_counts)
    }
}

// todo this function should be removed in favor of a more
//  generic version in pileup/mod.rs
pub fn process_region_duplex_batch<T: AsRef<Path> + Copy>(
    chromosome_coordintes: &MultiChromCoordinates,
    bam_fp: T,
    caller: &MultipleThresholdModCaller,
    pileup_numeric_options: &PileupNumericOptions,
    force_allow: bool,
    max_depth: u32,
    edge_filter: Option<&EdgeFilter>,
) -> Vec<anyhow::Result<DuplexModBasePileup>> {
    chromosome_coordintes
        .0
        .iter()
        .map(|chrom_coords| {
            process_region_duplex(
                bam_fp,
                chrom_coords.chrom_tid,
                chrom_coords.start_pos,
                chrom_coords.end_pos,
                caller,
                pileup_numeric_options,
                force_allow,
                max_depth,
                &chrom_coords.focus_positions,
                edge_filter,
            )
        })
        .collect()
}

fn process_region_duplex<T: AsRef<Path>>(
    bam_fp: T,
    chrom_tid: u32,
    start_pos: u32,
    end_pos: u32,
    caller: &MultipleThresholdModCaller,
    pileup_numeric_options: &PileupNumericOptions,
    force_allow: bool,
    max_depth: u32,
    focus_positions: &FocusPositions,
    edge_filter: Option<&EdgeFilter>,
) -> anyhow::Result<DuplexModBasePileup> {
    let positions_to_motifs = match focus_positions {
        FocusPositions::MotifCombineStrands { positive_motifs, .. } => {
            positive_motifs
        }
        _ => bail!("duplex requires a motif"),
    };

    let mut bam_reader = bam::IndexedReader::from_path(bam_fp)?;
    let chrom_name =
        String::from_utf8_lossy(bam_reader.header().tid2name(chrom_tid))
            .to_string();
    bam_reader.fetch(FetchDefinition::Region(
        chrom_tid as i32,
        start_pos as i64,
        end_pos as i64,
    ))?;
    let mut read_cache = DuplexReadCache::new(
        pileup_numeric_options.get_collapse_method(),
        caller,
        edge_filter,
        force_allow,
    );

    let mut position_feature_counts = FxHashMap::default();

    let hts_pileup = {
        let mut tmp_pileup = bam_reader.pileup();
        tmp_pileup.set_max_depth(max_depth);
        tmp_pileup
    };

    let pileup_iter =
        PileupIter::new(hts_pileup, start_pos, end_pos, focus_positions);

    for (pileup, motif) in pileup_iter.filter_map(|pileup| {
        let motifs = positions_to_motifs.get(&pileup.bam_pileup.pos())?;
        if motifs.len() > 1 {
            debug!("more than 1 motif not supported yet");
        };
        let (motif, _idx) = &motifs[0];
        Some((pileup, motif))
    }) {
        let pos = pileup.bam_pileup.pos();
        let mut feature_vector = DuplexFeatureVector::default();
        let alignment_iter =
            pileup.bam_pileup.alignments().filter(|alignment| {
                if alignment.is_refskip() {
                    false
                } else {
                    let record = alignment.record();
                    !(record_is_not_primary(&record) || record.seq_len() == 0)
                }
            });

        for alignment in alignment_iter {
            if alignment.is_del() {
                feature_vector.add_feature(
                    DuplexFeature::Delete,
                    &pileup_numeric_options,
                );
                continue;
            }
            let record = alignment.record();
            let read_base = get_forward_read_base(&alignment, &record);
            if read_base.is_none() {
                continue;
            }
            let read_base = read_base.unwrap();
            if let Some(duplex_mod_call) =
                read_cache.get_duplex_mod_call(&record, pos, read_base, motif)
            {
                feature_vector.add_feature(
                    DuplexFeature::ModCall(duplex_mod_call),
                    &pileup_numeric_options,
                );
            }
        }
        let pileup_counts = feature_vector.decode();
        position_feature_counts.insert(pos, pileup_counts);
    }

    let (processed_records, skipped_records) =
        read_cache.get_records_used_and_skipped();
    Ok(DuplexModBasePileup {
        chrom_name,
        pileup_counts: position_feature_counts,
        processed_records,
        skipped_records,
    })
}
