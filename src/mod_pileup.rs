use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::Path;

use derive_new::new;
use itertools::Itertools;
use log::{debug, error};
use rust_htslib::bam;
use rust_htslib::bam::{FetchDefinition, Read};

use crate::mod_bam::{BaseModCall, CollapseMethod, EdgeFilter};
use crate::mod_base_code::{DnaBase, ModCode};
use crate::motif_bed::{MotifLocations, RegexMotif};
use crate::read_cache::ReadCache;
use crate::threshold_mod_caller::MultipleThresholdModCaller;
use crate::util::{get_query_name_string, record_is_secondary, Strand};

#[derive(Debug, Copy, Clone)]
enum Feature {
    Delete,
    Filtered,
    NoCall(DnaBase),
    ModCall(ModCode),
}

impl Feature {
    fn from_base_mod_call(
        base_mod_call: BaseModCall,
        read_base: DnaBase,
    ) -> Self {
        match base_mod_call {
            BaseModCall::Canonical(_) => Feature::ModCall(
                read_base.canonical_mod_code().expect("should get base"),
            ),
            BaseModCall::Modified(_, mod_code) => Feature::ModCall(mod_code),
            BaseModCall::Filtered => Feature::Filtered,
        }
    }
}

#[derive(Debug, Copy, Clone, new, Default)]
pub struct PileupFeatureCounts {
    pub raw_strand: char,
    pub filtered_coverage: u32,
    pub raw_mod_code: char,
    pub fraction_modified: f32,
    pub n_canonical: u32,
    pub n_modified: u32,
    pub n_other_modified: u32,
    pub n_delete: u32,
    pub n_filtered: u32,
    pub n_diff: u32,
    pub n_nocall: u32,
}

impl PileupFeatureCounts {
    fn new_empty(raw_strand: char, raw_mod_code: char) -> Self {
        Self {
            raw_strand,
            raw_mod_code,
            ..Default::default()
        }
    }

    fn combine_counts_ignore_strand(self, other: Self) -> Self {
        if self.raw_mod_code != other.raw_mod_code {
            error!(
                "shouldn't be combining counts with different mod codes!\
            {} vs {}",
                self.raw_mod_code, other.raw_mod_code
            );
        }
        let n_modified = self.n_modified + other.n_modified;
        let n_canonical = self.n_canonical + other.n_canonical;
        let n_other_modified = self.n_other_modified + other.n_other_modified;
        let filtered_coverage =
            self.filtered_coverage + other.filtered_coverage;
        let n_delete = self.n_delete + other.n_delete;
        let n_filtered = self.n_filtered + other.n_filtered;
        let n_diff = self.n_diff + other.n_diff;
        let n_nocall = self.n_nocall + other.n_nocall;

        let fraction_modified = n_modified as f32 / filtered_coverage as f32;

        Self::new(
            self.raw_strand,
            filtered_coverage,
            self.raw_mod_code,
            fraction_modified,
            n_canonical,
            n_modified,
            n_other_modified,
            n_delete,
            n_filtered,
            n_diff,
            n_nocall,
        )
    }
}

#[allow(non_snake_case)]
#[derive(Debug, Default)]
struct Tally {
    n_delete: u32,
    n_filtered: u32,
    n_basecall_A: u32,
    n_basecall_C: u32,
    n_basecall_G: u32,
    n_basecall_T: u32,
    n_modcall_A: u32,
    n_modcall_C: u32,
    n_modcall_a: u32,
    n_modcall_h: u32,
    n_modcall_m: u32,
}

impl Tally {
    fn add_feature(&mut self, feature: Feature) {
        match feature {
            Feature::Filtered => self.n_filtered += 1,
            Feature::Delete => self.n_delete += 1,
            Feature::ModCall(mod_base) => match mod_base {
                ModCode::C => self.n_modcall_C += 1,
                ModCode::h => self.n_modcall_h += 1,
                ModCode::m => self.n_modcall_m += 1,
                ModCode::A => self.n_modcall_A += 1,
                ModCode::a => self.n_modcall_a += 1,
                _ => {}
            },
            Feature::NoCall(dna_base) => match dna_base {
                DnaBase::A => self.n_basecall_A += 1,
                DnaBase::C => self.n_basecall_C += 1,
                DnaBase::G => self.n_basecall_G += 1,
                DnaBase::T => self.n_basecall_T += 1,
            },
        }
    }
}

#[derive(Debug, Default)]
struct FeatureVector {
    pos_tally: Tally,
    neg_tally: Tally,
}

impl FeatureVector {
    pub(crate) fn new() -> Self {
        Self::default()
    }

    /// Add counts to the tally.
    pub(crate) fn add_feature(
        &mut self,
        alignment_strand: Strand,
        feature: Feature,
        read_strand: Strand,
        strand_rule: &StrandRule,
    ) {
        match strand_rule {
            StrandRule::Both => match (alignment_strand, read_strand) {
                (Strand::Positive, Strand::Positive) => {
                    self.pos_tally.add_feature(feature)
                }
                (Strand::Negative, Strand::Positive) => {
                    self.neg_tally.add_feature(feature)
                }

                (Strand::Positive, Strand::Negative) => {
                    self.neg_tally.add_feature(feature)
                }
                (Strand::Negative, Strand::Negative) => {
                    self.pos_tally.add_feature(feature)
                }
            },
            StrandRule::Positive => match (alignment_strand, read_strand) {
                (Strand::Positive, Strand::Positive) => {
                    self.pos_tally.add_feature(feature)
                }
                (Strand::Negative, Strand::Negative) => {
                    self.pos_tally.add_feature(feature)
                }
                _ => {}
            },
            StrandRule::Negative => match (alignment_strand, read_strand) {
                (Strand::Negative, Strand::Positive) => {
                    self.neg_tally.add_feature(feature)
                }

                (Strand::Positive, Strand::Negative) => {
                    self.neg_tally.add_feature(feature)
                }
                _ => {}
            },
        }
    }

    fn add_pileup_counts(
        pileup_options: &PileupNumericOptions,
        counts: &mut Vec<PileupFeatureCounts>,
        observed_mods: &HashSet<ModCode>,
        strand: Strand,
        filtered_coverage: u32,
        n_h: u32,
        n_m: u32,
        n_canonical: u32,
        n_delete: u32,
        n_filtered: u32,
        n_diff: u32,
        n_nocall: u32,
    ) {
        match pileup_options {
            PileupNumericOptions::Passthrough
            | PileupNumericOptions::Collapse(_) => {
                for (mod_code, (n_modified, n_other_modified)) in
                    [(ModCode::h, (n_h, n_m)), (ModCode::m, (n_m, n_h))]
                {
                    if observed_mods.contains(&mod_code) {
                        let percent_modified =
                            n_modified as f32 / filtered_coverage as f32;
                        counts.push(PileupFeatureCounts {
                            raw_strand: strand.to_char(),
                            filtered_coverage,
                            raw_mod_code: mod_code.char(),
                            fraction_modified: percent_modified,
                            n_canonical,
                            n_modified,
                            n_other_modified,
                            n_delete,
                            n_filtered,
                            n_diff,
                            n_nocall,
                        })
                    }
                }
            }
            PileupNumericOptions::Combine => {
                // todo maybe want a check for seen the mod here?
                let n_modified = n_h + n_m;
                let percent_modified =
                    n_modified as f32 / filtered_coverage as f32;
                counts.push(PileupFeatureCounts {
                    raw_strand: strand.to_char(),
                    filtered_coverage,
                    raw_mod_code: ModCode::C.char(),
                    fraction_modified: percent_modified,
                    n_canonical,
                    n_modified,
                    n_other_modified: 0,
                    n_delete,
                    n_filtered,
                    n_diff,
                    n_nocall,
                })
            }
        }
    }

    fn add_tally_to_counts(
        counts: &mut Vec<PileupFeatureCounts>,
        tally: &Tally,
        strand: Strand,
        observed_mods: &HashSet<ModCode>,
        pileup_options: &PileupNumericOptions,
    ) {
        if (tally.n_modcall_A + tally.n_modcall_a) > 0 {
            let n_canonical = tally.n_modcall_A;
            let n_mod = tally.n_modcall_a;
            let filtered_coverage = n_canonical + n_mod;
            let raw_mod_code = ModCode::a.char();
            let n_nocall = tally.n_basecall_A;
            let percent_modified =
                n_mod as f32 / (n_mod as f32 + n_canonical as f32);
            let n_diff = tally.n_basecall_C
                + tally.n_basecall_T
                + tally.n_basecall_G
                + tally.n_modcall_C
                + tally.n_modcall_m
                + tally.n_modcall_h;
            counts.push(PileupFeatureCounts {
                raw_strand: strand.to_char(),
                filtered_coverage,
                raw_mod_code,
                fraction_modified: percent_modified,
                n_canonical,
                n_modified: n_mod,
                n_other_modified: 0,
                n_delete: tally.n_delete,
                n_filtered: tally.n_filtered,
                n_diff,
                n_nocall,
            });
        }

        // + and - strand C-mods
        if (tally.n_modcall_h + tally.n_modcall_m + tally.n_modcall_C) > 0 {
            let n_canonical = tally.n_modcall_C;
            let n_nocall = tally.n_basecall_C;

            let n_diff = tally.n_basecall_A
                + tally.n_basecall_G
                + tally.n_basecall_T
                + tally.n_modcall_A
                + tally.n_modcall_a;

            let n_h = tally.n_modcall_h;
            let n_m = tally.n_modcall_m;
            let filtered_coverage = n_canonical + n_h + n_m;
            Self::add_pileup_counts(
                pileup_options,
                counts,
                observed_mods,
                strand,
                filtered_coverage,
                n_h,
                n_m,
                n_canonical,
                tally.n_delete,
                tally.n_filtered,
                n_diff,
                n_nocall,
            );
        }
    }

    pub fn decode(
        self,
        pos_observed_mods: &HashSet<ModCode>,
        neg_observed_mods: &HashSet<ModCode>,
        pileup_options: &PileupNumericOptions,
    ) -> Vec<PileupFeatureCounts> {
        let mut counts = Vec::new();
        Self::add_tally_to_counts(
            &mut counts,
            &self.pos_tally,
            Strand::Positive,
            pos_observed_mods,
            pileup_options,
        );
        Self::add_tally_to_counts(
            &mut counts,
            &self.neg_tally,
            Strand::Negative,
            neg_observed_mods,
            pileup_options,
        );

        counts
    }
}

fn combine_strand_features(
    positions: Option<&HashMap<u32, Strand>>,
    motif: Option<&RegexMotif>,
    mut position_feature_counts: HashMap<u32, Vec<PileupFeatureCounts>>,
) -> HashMap<u32, Vec<PileupFeatureCounts>> {
    if positions.is_none() || motif.is_none() {
        error!(
            "asked to combine feature counts with no motif information present"
        );
        return position_feature_counts;
    }

    let positions = positions.unwrap();
    let motif = motif.unwrap();
    let positions_to_combine =
        positions.iter().filter_map(|(pos, strand)| match strand {
            Strand::Positive => {
                let offset =
                    (motif.reverse_offset - motif.forward_offset) as u32;
                Some((*pos, pos + offset))
            }
            Strand::Negative => None,
        });

    let mut combined_position_feature_counts = HashMap::new();
    for (pos_position, neg_position) in positions_to_combine {
        let pos_features = position_feature_counts.remove(&pos_position);
        let neg_features = position_feature_counts.remove(&neg_position);
        match (pos_features, neg_features) {
            (Some(pos_feats), Some(neg_feats)) => {
                // important to use b-tree map here so that ordering is deterministic
                let mut grouped_by_mod_code = BTreeMap::new();
                for feature_counts in
                    pos_feats.into_iter().chain(neg_feats.into_iter())
                {
                    let feats_for_mod_code = grouped_by_mod_code
                        .entry(feature_counts.raw_mod_code)
                        .or_insert(Vec::new());
                    feats_for_mod_code.push(feature_counts);
                }
                let combined = grouped_by_mod_code
                    .into_iter()
                    .map(|(mod_code, feature_counts)| {
                        feature_counts.into_iter().fold(
                            // use unknown/ambiguous strand because we're combining
                            PileupFeatureCounts::new_empty('.', mod_code),
                            |acc, next| acc.combine_counts_ignore_strand(next),
                        )
                    })
                    .collect::<Vec<PileupFeatureCounts>>();
                combined_position_feature_counts.insert(pos_position, combined);
            }
            (None, Some(neg_feats)) => {
                let strand_switched_feats = neg_feats
                    .into_iter()
                    .map(|f| {
                        PileupFeatureCounts::new(
                            // use unknown/ambiguous strand because we're combining
                            '.',
                            f.filtered_coverage,
                            f.raw_mod_code,
                            f.fraction_modified,
                            f.n_canonical,
                            f.n_modified,
                            f.n_other_modified,
                            f.n_delete,
                            f.n_filtered,
                            f.n_diff,
                            f.n_nocall,
                        )
                    })
                    .collect::<Vec<_>>();
                combined_position_feature_counts
                    .insert(pos_position, strand_switched_feats);
            }
            (Some(pos_feats), None) => {
                combined_position_feature_counts
                    .insert(pos_position, pos_feats);
            }
            (None, None) => {}
        }
    }

    combined_position_feature_counts
}

fn get_motif_locations_for_region(
    motif_locations: &MotifLocations,
    reference_id: u32,
    start_pos: u32,
    end_pos: u32,
) -> HashMap<u32, Strand> {
    motif_locations
        .get_locations_unchecked(reference_id)
        .iter()
        .filter_map(|(pos, strand)| {
            if pos >= &start_pos && pos < &end_pos {
                Some((*pos, *strand))
            } else {
                None
            }
        })
        .collect()
}

enum StrandRule {
    Positive,
    Negative,
    Both,
}

#[derive(new)]
struct StrandPileup {
    pub(crate) bam_pileup: bam::pileup::Pileup,
    strand_rule: StrandRule,
}

struct PileupIter<'a> {
    pileups: bam::pileup::Pileups<'a, bam::IndexedReader>,
    start_pos: u32,
    end_pos: u32,
    motif_locations: Option<&'a HashMap<u32, Strand>>,
}

impl<'a> PileupIter<'a> {
    fn new(
        pileups: bam::pileup::Pileups<'a, bam::IndexedReader>,
        start_pos: u32,
        end_pos: u32,
        motif_locations: Option<&'a HashMap<u32, Strand>>,
    ) -> Self {
        Self {
            pileups,
            start_pos,
            end_pos,
            motif_locations,
        }
    }
}

impl<'a> Iterator for PileupIter<'a> {
    type Item = StrandPileup;

    fn next(&mut self) -> Option<Self::Item> {
        let mut pileup: Option<Self::Item> = None;
        while let Some(Ok(plp)) = self.pileups.next() {
            let off_end = plp.pos() >= self.end_pos;
            if off_end {
                // we're done
                return None;
            } else if plp.pos() < self.start_pos {
                // advance into region we're looking at
                continue;
            } else {
                match &self.motif_locations {
                    Some(locations) => {
                        if let Some(strand) = locations.get(&plp.pos()) {
                            let strand_rule = match strand {
                                Strand::Positive => StrandRule::Positive,
                                Strand::Negative => StrandRule::Negative,
                            };
                            pileup = Some(StrandPileup::new(plp, strand_rule));
                            break;
                        } else {
                            continue;
                        }
                    }
                    None => {
                        pileup = Some(StrandPileup::new(plp, StrandRule::Both));
                        break;
                    }
                }
            }
        }
        pileup
    }
}

pub struct ModBasePileup {
    pub chrom_name: String,
    position_feature_counts: HashMap<u32, Vec<PileupFeatureCounts>>,
    pub(crate) skipped_records: usize,
    pub(crate) processed_records: usize,
}

impl ModBasePileup {
    pub fn num_results(&self) -> usize {
        self.position_feature_counts.len()
    }

    pub fn iter_counts(
        &self,
    ) -> impl Iterator<Item = (&u32, &Vec<PileupFeatureCounts>)> {
        self.position_feature_counts
            .iter()
            .sorted_by(|(x, _), (y, _)| x.cmp(y))
    }
}

pub enum PileupNumericOptions {
    Passthrough,
    Combine,
    Collapse(CollapseMethod),
}

impl PileupNumericOptions {
    fn get_collapse_method(&self) -> Option<&CollapseMethod> {
        match self {
            Self::Collapse(method) => Some(method),
            _ => None,
        }
    }
}

pub fn process_region<T: AsRef<Path>>(
    bam_fp: T,
    chrom_tid: u32,
    start_pos: u32,
    end_pos: u32,
    caller: &MultipleThresholdModCaller,
    pileup_numeric_options: &PileupNumericOptions,
    force_allow: bool,
    combine_strands: bool,
    motif_locations: Option<&MotifLocations>,
    edge_filter: Option<&EdgeFilter>,
) -> Result<ModBasePileup, String> {
    let mut bam_reader =
        bam::IndexedReader::from_path(bam_fp).map_err(|e| e.to_string())?;
    let chrom_name =
        String::from_utf8_lossy(bam_reader.header().tid2name(chrom_tid))
            .to_string();
    bam_reader
        .fetch(FetchDefinition::Region(
            chrom_tid as i32,
            start_pos as i64,
            end_pos as i64,
        ))
        .map_err(|e| e.to_string())?;

    let motif_positions = motif_locations.map(|mls| {
        get_motif_locations_for_region(mls, chrom_tid, start_pos, end_pos)
    });

    let mut read_cache = ReadCache::new(
        pileup_numeric_options.get_collapse_method(),
        caller,
        edge_filter,
        force_allow,
    );
    let mut position_feature_counts = HashMap::new();
    let pileup_iter = PileupIter::new(
        bam_reader.pileup(),
        start_pos,
        end_pos,
        motif_positions.as_ref(),
    );
    let mut dupe_reads = HashMap::new();
    for pileup in pileup_iter {
        let mut feature_vector = FeatureVector::new();
        let mut pos_strand_observed_mod_codes = HashSet::new();
        let mut neg_strand_observed_mod_codes = HashSet::new();
        let pos = pileup.bam_pileup.pos();

        // structures around warning about dupes
        let mut observed_read_ids_to_pos = HashMap::new();

        let alignment_iter =
            pileup.bam_pileup.alignments().filter(|alignment| {
                if alignment.is_refskip() {
                    false
                } else {
                    let record = alignment.record();
                    !(record_is_secondary(&record) || record.seq_len() == 0)
                }
            });
        for alignment in alignment_iter {
            assert!(!alignment.is_refskip());
            let record = alignment.record();
            read_cache.add_mod_codes_for_record(
                &record,
                &mut pos_strand_observed_mod_codes,
                &mut neg_strand_observed_mod_codes,
            );

            if let Ok(read_name) = get_query_name_string(&record) {
                (*observed_read_ids_to_pos
                    .entry(read_name)
                    .or_insert(0usize)) += 1
            }

            // alignment stand is the strand the read is aligned to
            let alignment_strand = if record.is_reverse() {
                Strand::Negative
            } else {
                Strand::Positive
            };

            if alignment.is_del() {
                feature_vector.add_feature(
                    alignment_strand,
                    Feature::Delete,
                    Strand::Positive,
                    &pileup.strand_rule,
                );
                continue;
            }

            // not delete or skip, add base
            let read_base = alignment.qpos().and_then(|pos| {
                if pos >= record.seq_len() {
                    debug!("Record position is not included in sequence?");
                    None
                } else {
                    DnaBase::parse(record.seq()[pos] as char).ok()
                }
            });

            let read_base = if let Some(base) = read_base {
                if record.is_reverse() {
                    base.complement()
                } else {
                    base
                }
            } else {
                continue;
            };

            match read_cache.get_mod_call(&record, pos, read_base.char()) {
                // a read can report on the read-positive or read-negative
                // strand (see the docs for .get_mod_call above) so the
                // pos_call and neg_call below are _read oriented_, the
                // `read_strand` in add_feature (see the docs there too)
                // is meant to pass along the information regarding which
                // strand of a read the feature belongs to. In almost all
                // cases this is Positive, because we sequence single
                // stranded DNA. However, for duplex and other double-
                // stranded techs, you could have a read with a mod call on
                // the negative strand. You must pass along the
                // `alignment_strand` so that everything can be oriented to
                // the positive strand of the reference.
                (Some(pos_call), Some(neg_call)) => {
                    let pos_feature =
                        Feature::from_base_mod_call(pos_call, read_base);
                    let neg_feature = Feature::from_base_mod_call(
                        neg_call,
                        read_base.complement(),
                    );
                    feature_vector.add_feature(
                        alignment_strand,
                        pos_feature,
                        Strand::Positive,
                        &pileup.strand_rule,
                    );
                    feature_vector.add_feature(
                        alignment_strand,
                        neg_feature,
                        Strand::Negative,
                        &pileup.strand_rule,
                    );
                }
                (Some(pos_call), None) => {
                    let pos_feature =
                        Feature::from_base_mod_call(pos_call, read_base);
                    feature_vector.add_feature(
                        alignment_strand,
                        pos_feature,
                        Strand::Positive,
                        &pileup.strand_rule,
                    );
                }
                (None, Some(neg_call)) => {
                    let neg_feature = Feature::from_base_mod_call(
                        neg_call,
                        read_base.complement(),
                    );

                    feature_vector.add_feature(
                        alignment_strand,
                        neg_feature,
                        Strand::Negative,
                        &pileup.strand_rule,
                    );
                }
                (None, None) => feature_vector.add_feature(
                    alignment_strand,
                    Feature::NoCall(read_base),
                    Strand::Positive,
                    &pileup.strand_rule,
                ),
            }
        } // alignment loop
        position_feature_counts.insert(
            pos,
            feature_vector.decode(
                &pos_strand_observed_mod_codes,
                &neg_strand_observed_mod_codes,
                &pileup_numeric_options,
            ),
        );
        observed_read_ids_to_pos
            .into_iter()
            .filter(|(_, count)| *count > 1usize)
            .for_each(|(read_id, count)| {
                dupe_reads.entry(read_id).or_insert(Vec::new()).push(count);
            })
    } // position loop

    let position_feature_counts = if combine_strands {
        combine_strand_features(
            motif_positions.as_ref(),
            motif_locations.map(|mls| mls.motif()),
            position_feature_counts,
        )
    } else {
        position_feature_counts
    };
    let (processed_records, skipped_records) =
        read_cache.get_records_used_and_skipped();

    let should_warn = !dupe_reads.is_empty();
    for (read_id, counts) in dupe_reads {
        let avg_times =
            counts.iter().map(|c| *c as f32).sum::<f32>() / counts.len() as f32;
        debug!(
            "read {read_id} was observed multiple times, \
            avg {avg_times} at {} positions on contig {chrom_name}, \
            between {start_pos} and {end_pos}",
            counts.len()
        );
    }
    if should_warn {
        debug!("consider marking duplicate alignments");
    }

    Ok(ModBasePileup {
        chrom_name,
        position_feature_counts,
        processed_records,
        skipped_records,
    })
}

#[cfg(test)]
mod mod_pileup_tests {
    use std::collections::HashSet;

    use crate::mod_pileup::{
        DnaBase, Feature, FeatureVector, ModCode, PileupNumericOptions,
        StrandRule,
    };
    use crate::util::Strand;

    #[test]
    fn test_feature_vector_basic() {
        let pos_observed_mods = HashSet::from([ModCode::m, ModCode::h]);
        let neg_observed_mods = HashSet::new();
        let mut fv = FeatureVector::new();
        fv.add_feature(
            Strand::Positive,
            Feature::NoCall(DnaBase::A),
            Strand::Positive,
            &StrandRule::Both,
        );
        fv.add_feature(
            Strand::Positive,
            Feature::ModCall(ModCode::C),
            Strand::Positive,
            &StrandRule::Both,
        );
        fv.add_feature(
            Strand::Positive,
            Feature::ModCall(ModCode::m),
            Strand::Positive,
            &StrandRule::Both,
        );
        fv.add_feature(
            Strand::Positive,
            Feature::ModCall(ModCode::m),
            Strand::Positive,
            &StrandRule::Both,
        );
        fv.add_feature(
            Strand::Positive,
            Feature::NoCall(DnaBase::C),
            Strand::Positive,
            &StrandRule::Both,
        );
        fv.add_feature(
            Strand::Negative,
            Feature::NoCall(DnaBase::G),
            Strand::Positive,
            &StrandRule::Both,
        );
        fv.add_feature(
            Strand::Negative,
            Feature::NoCall(DnaBase::G),
            Strand::Positive,
            &StrandRule::Both,
        );
        let counts = fv.decode(
            &pos_observed_mods,
            &neg_observed_mods,
            &PileupNumericOptions::Passthrough,
        );
        assert_eq!(counts.len(), 2); // h and m, negative strand should not be there
        for pileup_counts in counts {
            assert_eq!(pileup_counts.filtered_coverage, 3);
            assert_eq!(pileup_counts.n_nocall, 1);
            assert_eq!(pileup_counts.n_diff, 1);
            assert_eq!(pileup_counts.raw_strand, Strand::Positive.to_char());
        }
        let mut fv = FeatureVector::new();
        let neg_observed_mods = HashSet::from([ModCode::m, ModCode::h]);
        fv.add_feature(
            Strand::Positive,
            Feature::ModCall(ModCode::C),
            Strand::Positive,
            &StrandRule::Both,
        );
        fv.add_feature(
            Strand::Negative,
            Feature::ModCall(ModCode::m),
            Strand::Positive,
            &StrandRule::Both,
        );
        fv.add_feature(
            Strand::Negative,
            Feature::NoCall(DnaBase::G),
            Strand::Positive,
            &StrandRule::Both,
        );
        fv.add_feature(
            Strand::Negative,
            Feature::NoCall(DnaBase::G),
            Strand::Positive,
            &StrandRule::Both,
        );
        let counts = fv.decode(
            &pos_observed_mods,
            &neg_observed_mods,
            &PileupNumericOptions::Passthrough,
        );
        assert_eq!(counts.len(), 4);
        counts
            .iter()
            .filter(|c| c.raw_strand == Strand::Negative.to_char())
            .for_each(|c| assert_eq!(c.n_diff, 2));
    }

    #[test]
    fn test_feature_vector_with_strand_rules() {
        let mut fv = FeatureVector::new();
        let pos_observed_mods = HashSet::from([ModCode::m]);
        fv.add_feature(
            Strand::Positive,
            Feature::ModCall(ModCode::m),
            Strand::Positive,
            &StrandRule::Positive,
        );
        // this feature should be ignored because it's on the wrong
        // strand
        fv.add_feature(
            Strand::Negative,
            Feature::ModCall(ModCode::m),
            Strand::Positive,
            &StrandRule::Positive,
        );
        let counts = fv.decode(
            &pos_observed_mods,
            &HashSet::new(),
            &PileupNumericOptions::Passthrough,
        );
        assert_eq!(counts.len(), 1);
        let count = &counts[0];
        // change alignment strand to Positive and this will be 2
        assert_eq!(count.n_modified, 1);
    }
}
