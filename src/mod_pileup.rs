use crate::mod_bam::{BaseModCall, CollapseMethod};
use crate::mod_base_code::{DnaBase, ModCode};
use crate::read_cache::ReadCache;
use crate::util::{record_is_secondary, Strand};
use itertools::Itertools;
use log::debug;
use rust_htslib::bam;
use rust_htslib::bam::{FetchDefinition, Read};
use std::collections::{HashMap, HashSet};
use std::path::Path;

#[derive(Debug, Copy, Clone)]
enum Feature {
    Delete,
    Filtered,
    NoCall(DnaBase),
    ModCall(ModCode),
}

#[derive(Debug)]
pub struct PileupFeatureCounts {
    pub strand: Strand,
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
    fn new() -> Self {
        Self::default()
    }

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

    pub(crate) fn add_feature(&mut self, strand: Strand, feature: Feature) {
        match strand {
            Strand::Positive => self.pos_tally.add_feature(feature),
            Strand::Negative => self.neg_tally.add_feature(feature),
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
                            strand,
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
                let n_modified = n_h + n_m;
                let percent_modified =
                    n_modified as f32 / filtered_coverage as f32;
                counts.push(PileupFeatureCounts {
                    strand,
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
                strand,
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

        // + strand C-mods
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

struct PileupIter<'a> {
    pileups: bam::pileup::Pileups<'a, bam::IndexedReader>,
    start_pos: u32,
    end_pos: u32,
}

impl<'a> PileupIter<'a> {
    fn new(
        pileups: bam::pileup::Pileups<'a, bam::IndexedReader>,
        start_pos: u32,
        end_pos: u32,
    ) -> Self {
        Self {
            pileups,
            start_pos,
            end_pos,
        }
    }
}

impl<'a> Iterator for PileupIter<'a> {
    type Item = bam::pileup::Pileup;

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
                pileup = Some(plp);
                break;
            }
        }
        pileup
    }
}

pub struct ModBasePileup {
    pub chrom_name: String,
    position_feature_counts: HashMap<u32, Vec<PileupFeatureCounts>>,
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
    fn get_collapse_method(&self) -> Option<CollapseMethod> {
        match self {
            Self::Collapse(method) => Some(*method),
            _ => None,
        }
    }
}

pub fn process_region<T: AsRef<Path>>(
    bam_fp: T,
    chrom_tid: u32,
    start_pos: u32,
    end_pos: u32,
    threshold: f32,
    pileup_numeric_options: &PileupNumericOptions,
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

    let mut read_cache =
        ReadCache::new(pileup_numeric_options.get_collapse_method());
    let mut position_feature_counts = HashMap::new();
    let pileup_iter = PileupIter::new(bam_reader.pileup(), start_pos, end_pos);
    for pileup in pileup_iter {
        let mut feature_vector = FeatureVector::new();
        let mut pos_strand_observed_mod_codes = HashSet::new();
        let mut neg_strand_observed_mod_codes = HashSet::new();
        // let mut observed_mod_codes = HashSet::new();
        let pos = pileup.pos();

        let alignment_iter = pileup.alignments().filter_map(|alignment| {
            if alignment.is_refskip() {
                None
            } else {
                let record = alignment.record();
                if record_is_secondary(&record) || record.seq_len() == 0 {
                    None
                } else {
                    Some(alignment)
                }
            }
        });
        for alignment in alignment_iter {
            assert!(!alignment.is_refskip());
            let record = alignment.record();
            // observed_mod_codes
            //     .extend(read_cache.get_mod_codes_for_record(&record));
            let strand = if record.is_reverse() {
                neg_strand_observed_mod_codes
                    .extend(read_cache.get_mod_codes_for_record(&record));
                Strand::Negative
            } else {
                pos_strand_observed_mod_codes
                    .extend(read_cache.get_mod_codes_for_record(&record));
                Strand::Positive
            };

            if alignment.is_del() {
                feature_vector.add_feature(strand, Feature::Delete);
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

            let feature = if let Some(mod_call) = read_cache.get_mod_call(
                &record,
                pos,
                read_base.char(),
                threshold,
            ) {
                match mod_call {
                    BaseModCall::Canonical(_) => Feature::ModCall(
                        read_base.canonical_mod_code().unwrap(),
                    ),
                    BaseModCall::Filtered => Feature::Filtered,
                    BaseModCall::Modified(_, raw_code) => Feature::ModCall(
                        ModCode::parse_raw_mod_code(raw_code).unwrap(),
                    ),
                }
            } else {
                Feature::NoCall(read_base)
            };
            feature_vector.add_feature(strand, feature);
        } // alignment loop
        position_feature_counts.insert(
            pos,
            feature_vector.decode(
                &pos_strand_observed_mod_codes,
                &neg_strand_observed_mod_codes,
                &pileup_numeric_options,
            ),
        );
    } // position loop

    Ok(ModBasePileup {
        chrom_name,
        position_feature_counts,
    })
}

#[cfg(test)]
mod mod_pileup_tests {
    use crate::mod_pileup::{
        DnaBase, Feature, FeatureVector, ModCode, PileupNumericOptions,
    };
    use crate::util::Strand;
    use std::collections::HashSet;

    #[test]
    fn test_feature_vector() {
        let pos_observed_mods = HashSet::from([ModCode::m, ModCode::h]);
        let neg_observed_mods = HashSet::new();
        let mut fv = FeatureVector::new();
        fv.add_feature(Strand::Positive, Feature::NoCall(DnaBase::A));
        fv.add_feature(Strand::Positive, Feature::ModCall(ModCode::C));
        fv.add_feature(Strand::Positive, Feature::ModCall(ModCode::m));
        fv.add_feature(Strand::Positive, Feature::ModCall(ModCode::m));
        fv.add_feature(Strand::Positive, Feature::NoCall(DnaBase::C));
        fv.add_feature(Strand::Negative, Feature::NoCall(DnaBase::G));
        fv.add_feature(Strand::Negative, Feature::NoCall(DnaBase::G));
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
            assert_eq!(pileup_counts.strand, Strand::Positive);
        }
        let mut fv = FeatureVector::new();
        let neg_observed_mods = HashSet::from([ModCode::m, ModCode::h]);
        fv.add_feature(Strand::Positive, Feature::ModCall(ModCode::C));
        fv.add_feature(Strand::Negative, Feature::ModCall(ModCode::m));
        fv.add_feature(Strand::Negative, Feature::NoCall(DnaBase::G));
        fv.add_feature(Strand::Negative, Feature::NoCall(DnaBase::G));
        let counts = fv.decode(
            &pos_observed_mods,
            &neg_observed_mods,
            &PileupNumericOptions::Passthrough,
        );
        assert_eq!(counts.len(), 4);
        counts
            .iter()
            .filter(|c| c.strand == Strand::Negative)
            .for_each(|c| assert_eq!(c.n_diff, 2));
    }
}
