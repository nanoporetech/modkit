use std::{
    cmp::Ordering, collections::HashMap, marker::PhantomData, ops::Range,
    path::PathBuf,
};

use anyhow::{anyhow, bail};
use bitvec::slice::BitSlice;
use common_macros::hash_map;
use itertools::Itertools;
use rustc_hash::{FxHashMap, FxHashSet};
use std::ops::BitOrAssign;

use crate::{
    errs::{MkError, MkResult},
    interval_chunks::{ChromCoordinates, FocusPositions2},
    mod_bam::{
        validate_mn_tag_on_record, BaseModCall, BaseModProbs, EdgeFilter,
        MmTagInfo,
    },
    mod_base_code::{
        DnaBase, ModCodeRepr, ANY_CYTOSINE, METHYL_CYTOSINE, SIX_METHYL_ADENINE,
    },
    motifs::motif_bed::MotifInfo,
    pileup::{
        base_mods_adapter::{BaseModsAdapter, ModState},
        ModBasePileup2, PileupFeatureCounts2, PileupNumericOptions,
    },
    threshold_mod_caller::MultipleThresholdModCaller,
    util::{
        get_haplotype_tag, qual_to_prob, reader_is_cram, record_is_not_primary,
        record_is_primary, SamTag, Strand,
    },
};
use rust_htslib::{
    bam::{
        self, ext::BamRecordExtensions, record::BaseModificationsPositionIter,
        FetchDefinition, Read,
    },
    htslib::{self},
};

pub trait PileupWorker
where
    Self: Send,
{
    fn process(
        &mut self,
        chrom_coordinates: ChromCoordinates,
        pileup_space: ModBasePileup2,
    ) -> anyhow::Result<ModBasePileup2>;
}

pub(super) struct DnaPileupWorker<
    T,
    M,
    const STRANDS: usize,
    const CHECK_DEPTH: bool,
> {
    reader: bam::IndexedReader,
    phased: bool,
    dna_mod_option: DnaModOption,
    combine_mods: bool,
    allow_non_primary: bool,
    filter_thresholds: [f32; 4],
    mod_thresholds: Vec<(ModCodeRepr, f32)>,
    motif_bases: [DnaBase; 4],
    edge_filter_start: usize,
    edge_filter_end: usize,
    matrix: M,
    t: PhantomData<T>,
}

impl<
        T: Send + Sync,
        M: ACountsMatrix<T, STRANDS, CHECK_DEPTH> + Sync + Send,
        const STRANDS: usize,
        const CHECK_DEPTH: bool,
    > DnaPileupWorker<T, M, STRANDS, CHECK_DEPTH>
{
    pub(crate) fn new(
        bam_fp: &PathBuf,
        reference_fp: &PathBuf,
        phased: bool,
        filter_thresholds: [f32; 4],
        motif_bases: [DnaBase; 4],
        dna_mod_option: DnaModOption,
        allow_non_primary: bool,
        mod_codes: Vec<(DnaBase, ModCodeRepr)>,
        mod_thresholds: Vec<(ModCodeRepr, f32)>,
        edge_filter: Option<&EdgeFilter>,
        motif_offset: u32,
        max_depth: u16,
    ) -> anyhow::Result<Self> {
        let combine_mods = match dna_mod_option {
            DnaModOption::Combine => true,
            _ => false,
        };
        let mut reader = bam::IndexedReader::from_path(bam_fp)?;
        if reader_is_cram(&reader) {
            reader.set_reference(reference_fp)?;
            reader.set_cram_options(
                htslib::hts_fmt_option_CRAM_OPT_REQUIRED_FIELDS,
                htslib::sam_fields_SAM_FLAG
                    | htslib::sam_fields_SAM_RNAME
                    | htslib::sam_fields_SAM_POS
                    | htslib::sam_fields_SAM_MAPQ
                    | htslib::sam_fields_SAM_CIGAR
                    | htslib::sam_fields_SAM_SEQ
                    | htslib::sam_fields_SAM_AUX,
            )?;
        };
        let edge_filter_start =
            edge_filter.map(|ef| ef.edge_filter_start).unwrap_or(0usize);
        let edge_filter_end =
            edge_filter.map(|ef| ef.edge_filter_end).unwrap_or(0usize);
        let matrix = M::new(0usize, phased, mod_codes, motif_offset, max_depth);

        Ok(Self {
            reader,
            phased,
            dna_mod_option,
            combine_mods,
            allow_non_primary,
            filter_thresholds,
            mod_thresholds,
            edge_filter_start,
            edge_filter_end,
            motif_bases,
            matrix,
            t: PhantomData::<T>,
        })
    }
}

impl<
        T: Send + Sync,
        M: ACountsMatrix<T, STRANDS, CHECK_DEPTH> + Sync + Send,
        const STRANDS: usize,
        const CHECK_DEPTH: bool,
    > PileupWorker for DnaPileupWorker<T, M, STRANDS, CHECK_DEPTH>
{
    fn process(
        &mut self,
        item: ChromCoordinates,
        mut pileup_space: ModBasePileup2,
    ) -> anyhow::Result<ModBasePileup2> {
        let mut erred_records = 0usize;
        let chrom_tid = item.chrom_tid;
        let start_pos = item.start_pos;
        let end_pos = item.end_pos;

        let width = end_pos
            .checked_sub(start_pos)
            .ok_or_else(|| anyhow!("interval end before start"))?
            as usize;
        self.matrix.reset(width, self.phased);

        let chrom_name =
            String::from_utf8_lossy(self.reader.header().tid2name(chrom_tid))
                .to_string();

        self.reader.fetch(FetchDefinition::Region(
            chrom_tid as i32,
            start_pos as i64,
            end_pos as i64,
        ))?;
        let records = self
            .reader
            .records()
            .filter_ok(|record| record.tid() >= 0i32)
            .filter_ok(|record| {
                let read_length = record.seq_len();
                !(read_length <= self.edge_filter_start
                    || read_length <= self.edge_filter_end)
            })
            .filter_ok(|record| {
                record_is_primary(record) || self.allow_non_primary
            })
            .filter_map(|res| match res {
                Ok(record) => {
                    if self.phased {
                        let hp =
                            get_haplotype_tag(&record, &HPTAG).unwrap_or(0u8);
                        Some((record, hp))
                    } else {
                        Some((record, 0u8))
                    }
                }
                // TODO: tabulate this error
                Err(_) => None,
            });

        'records: for (record, hp) in records {
            if self.allow_non_primary && record_is_not_primary(&record) {
                if validate_mn_tag_on_record(&record).is_err() {
                    erred_records = erred_records.saturating_add(1);
                    continue 'records;
                }
            }

            let reverse = record.is_reverse();
            let read_length = record.seq_len();
            let Ok(mut modbase_iter) = BaseModsAdapter::<16>::new(&record)
            else {
                erred_records = erred_records.saturating_add(1);
                continue 'records;
            };

            let Ok(mut mod_state) = modbase_iter.next_modified_position(
                self.filter_thresholds,
                &self.mod_thresholds,
            ) else {
                erred_records = erred_records.saturating_add(1);
                continue 'records;
            };

            let aligned_pairs_iter = record
                .aligned_pairs_full()
                // TODO: invalid negative numbers?
                .filter_map(|pair| {
                    let qpos = pair[0].and_then(|x| {
                        if x < 0 {
                            erred_records = erred_records.saturating_add(1);
                            return None;
                        } else {
                            Some(x as usize)
                        }
                    });
                    let r_pos = pair[1];
                    match r_pos {
                        None => None,
                        Some(pos) if pos < 0 => None,
                        Some(pos) => {
                            let pos = pos as u32;
                            if pos >= start_pos && pos < end_pos {
                                Some((
                                    qpos,
                                    pos.checked_sub(start_pos).unwrap(),
                                ))
                            } else {
                                None
                            }
                        }
                    }
                })
                .filter_map(|(qpos, rpos)| match &item.focus_positions {
                    FocusPositions2::MotifMask { mask, num_motifs } => {
                        let st = (rpos * 2 * (*num_motifs as u32))
                            + (*num_motifs * reverse as u8) as u32;
                        let st = st as usize;
                        let end = st + (*num_motifs as usize);
                        let bs = &mask[st..end];
                        debug_assert_eq!(
                            bs.len(),
                            *num_motifs as usize,
                            "bs:{bs:?} st:{st} end:{end} \
                             num_motifs:{num_motifs}"
                        );
                        if bs[0] {
                            Some((qpos, rpos, self.motif_bases[0]))
                        } else if *num_motifs >= 2u8 && bs[1] {
                            Some((qpos, rpos, self.motif_bases[1]))
                        } else if *num_motifs >= 3u8 && bs[2] {
                            Some((qpos, rpos, self.motif_bases[2]))
                        } else if *num_motifs >= 4u8 && bs[3] {
                            Some((qpos, rpos, self.motif_bases[3]))
                        } else {
                            None
                        }
                    }
                    _ => unreachable!(),
                })
                .filter_map(|(qpos, rpos, ref_base)| match qpos {
                    Some(qpos)
                        if qpos >= self.edge_filter_start
                            && qpos < (read_length - self.edge_filter_end) =>
                    {
                        Some((Some(qpos), rpos, ref_base))
                    }
                    Some(_) => None,
                    None => Some((None, rpos, ref_base)),
                });
            'pileup: for (qpos, rpos, ref_base) in aligned_pairs_iter {
                if CHECK_DEPTH {
                    if self.matrix.reached_max_depth(rpos, reverse, hp) {
                        continue 'pileup;
                    }
                }
                match (qpos, mod_state) {
                    // at this aligned position we have a base modification call
                    (Some(q), Some(mp)) if q == mp.mod_position => {
                        self.matrix.incr_call(
                            mp,
                            rpos,
                            ref_base,
                            reverse,
                            hp,
                            self.combine_mods,
                        );

                        if let Ok(next_mod_state) = modbase_iter
                            .next_modified_position(
                                self.filter_thresholds,
                                &self.mod_thresholds,
                            )
                        {
                            mod_state = next_mod_state;
                        } else {
                            erred_records = erred_records.saturating_add(1);
                            continue 'records;
                        }
                        continue 'pileup;
                    }
                    // at this aligned position, but there is no modification
                    // call
                    (Some(q), Some(mp)) if q < mp.mod_position => {
                        let base = {
                            let Ok(tmp) = DnaBase::try_from(record.seq()[q])
                            else {
                                erred_records = erred_records.saturating_add(1);
                                continue 'records;
                            };
                            if record.is_reverse() {
                                tmp.complement()
                            } else {
                                tmp
                            }
                        };
                        self.matrix
                            .incr_diff_call(rpos, base, ref_base, reverse, hp);
                    }
                    // we're at an aligned pair to the right of the last
                    // modified position
                    (Some(q), Some(mp)) if q > mp.mod_position => {
                        'overran: loop {
                            match modbase_iter.next_modified_position(
                                self.filter_thresholds,
                                &self.mod_thresholds,
                            ) {
                                Ok(Some(ms)) => {
                                    let pos = ms.mod_position;
                                    if q > pos {
                                        // discard a modification call that is
                                        // not aligned
                                        continue 'overran;
                                    } else {
                                        if q == pos {
                                            self.matrix.incr_call(
                                                ms,
                                                rpos,
                                                ref_base,
                                                reverse,
                                                hp,
                                                self.combine_mods,
                                            );
                                            mod_state = Some(ms);
                                            break 'overran;
                                        } else {
                                            assert!(pos > q);
                                            let base = {
                                                let tmp = DnaBase::try_from(
                                                    record.seq()[q],
                                                )
                                                .unwrap();
                                                if record.is_reverse() {
                                                    tmp.complement()
                                                } else {
                                                    tmp
                                                }
                                            };
                                            self.matrix.incr_diff_call(
                                                rpos, base, ref_base, reverse,
                                                hp,
                                            );
                                            mod_state = Some(ms);
                                            break 'overran;
                                        }
                                    }
                                }
                                Err(_) => {
                                    erred_records =
                                        erred_records.saturating_add(1);
                                    continue 'records;
                                }
                                Ok(None) => {
                                    let base = {
                                        let tmp =
                                            DnaBase::try_from(record.seq()[q])
                                                .unwrap();
                                        if record.is_reverse() {
                                            tmp.complement()
                                        } else {
                                            tmp
                                        }
                                    };
                                    self.matrix.incr_diff_call(
                                        rpos, base, ref_base, reverse, hp,
                                    );
                                    mod_state = None;
                                    break 'overran;
                                }
                            }
                        }
                        continue 'pileup;
                    }
                    (None, _) => {
                        self.matrix.incr_delete(rpos, reverse, hp);
                    }
                    (Some(q), None) => {
                        let base = {
                            let tmp =
                                DnaBase::try_from(record.seq()[q]).unwrap();
                            if record.is_reverse() {
                                tmp.complement()
                            } else {
                                tmp
                            }
                        };
                        self.matrix
                            .incr_diff_call(rpos, base, ref_base, reverse, hp);
                    }
                    (Some(_), Some(_)) => unreachable!(),
                }
            }
        }

        let position_feature_counts =
            self.matrix.decode(start_pos, self.dna_mod_option, self.phased);

        if self.phased {
            let hp_width =
                width * STRANDS * self.matrix.stride(self.dna_mod_option);
            let mut combined = position_feature_counts;
            let mut hp1 = combined.split_off(hp_width);
            let hp2 = hp1.split_off(hp_width);
            debug_assert_eq!(combined.len(), hp1.len());
            debug_assert_eq!(hp2.len(), hp1.len());
            let combined_counts = combined
                .into_iter()
                .zip(hp2.iter())
                .zip(hp1.iter())
                .map(|((mut combined, hp2), hp1)| {
                    combined.add(hp1);
                    combined.add(hp2);
                    combined
                })
                .filter(|x| x.is_valid())
                .collect::<Vec<PileupFeatureCounts2>>();
            let hp1 =
                hp1.into_iter().filter(|x| x.is_valid()).collect::<Vec<_>>();
            let hp2 =
                hp2.into_iter().filter(|x| x.is_valid()).collect::<Vec<_>>();
            pileup_space.chrom_name = chrom_name;
            pileup_space.interval_width = width;
            pileup_space.stride = self.matrix.stride(self.dna_mod_option);
            pileup_space.position_feature_counts = combined_counts;
            pileup_space.phased_feature_counts = [hp1, hp2];
            pileup_space.failed_records = erred_records;
            Ok(pileup_space)
        } else {
            pileup_space.chrom_name = chrom_name;
            pileup_space.interval_width = width;
            pileup_space.stride = self.matrix.stride(self.dna_mod_option);
            pileup_space.position_feature_counts = position_feature_counts
                .into_iter()
                .filter(|x| x.is_valid())
                .collect();
            pileup_space.failed_records = erred_records;
            Ok(pileup_space)
        }
    }
}

#[derive(Debug, Copy, Clone)]
pub(crate) enum DnaModOption {
    Other(ModCodeRepr),
    Single,
    Combine,
}

pub struct GenericPileupWorker {
    reader: bam::IndexedReader,
    motif_bases: Vec<MotifInfo>,
    caller: MultipleThresholdModCaller,
    pileup_numeric_options: PileupNumericOptions,
    combine_strands: bool,
    max_depth: u16,
}

impl GenericPileupWorker {
    pub(crate) fn new(
        in_bam_fp: &PathBuf,
        reference_fp: Option<&PathBuf>,
        motif_bases: Vec<MotifInfo>,
        default_threshold: f32,
        base_thresholds: [f32; 4],
        mod_code_thresholds: &[(ModCodeRepr, f32)],
        pileup_numeric_options: PileupNumericOptions,
        combine_strands: bool,
        max_depth: u16,
    ) -> anyhow::Result<Self> {
        let mut reader = bam::IndexedReader::from_path(in_bam_fp)?;
        if reader_is_cram(&reader) {
            if let Some(reference_fp) = reference_fp {
                reader.set_reference(reference_fp)?;
                reader.set_cram_options(
                    htslib::hts_fmt_option_CRAM_OPT_REQUIRED_FIELDS,
                    htslib::sam_fields_SAM_FLAG
                        | htslib::sam_fields_SAM_RNAME
                        | htslib::sam_fields_SAM_POS
                        | htslib::sam_fields_SAM_MAPQ
                        | htslib::sam_fields_SAM_CIGAR
                        | htslib::sam_fields_SAM_SEQ
                        | htslib::sam_fields_SAM_AUX,
                )?;
            } else {
                bail!("CRAM input requires reference")
            }
        };

        let per_base_thresholds = hash_map! {
            DnaBase::A => base_thresholds[DnaBase::A as usize],
            DnaBase::C => base_thresholds[DnaBase::C as usize],
            DnaBase::G => base_thresholds[DnaBase::G as usize],
            DnaBase::T => base_thresholds[DnaBase::T as usize],
        };
        let per_mod_thresholds = mod_code_thresholds
            .iter()
            .copied()
            .collect::<HashMap<ModCodeRepr, f32>>();
        let caller = MultipleThresholdModCaller::new(
            per_base_thresholds,
            per_mod_thresholds,
            default_threshold,
        );
        Ok(Self {
            reader,
            motif_bases,
            caller,
            pileup_numeric_options,
            combine_strands,
            max_depth,
        })
    }
}
impl PileupWorker for GenericPileupWorker {
    fn process(
        &mut self,
        chrom_coordinates: ChromCoordinates,
        mut pileup_space: ModBasePileup2,
    ) -> anyhow::Result<ModBasePileup2> {
        let mut erred_records = 0usize;
        let chrom_tid = chrom_coordinates.chrom_tid;
        let start_pos = chrom_coordinates.start_pos;
        let end_pos = chrom_coordinates.end_pos;

        let width = end_pos
            .checked_sub(start_pos)
            .ok_or_else(|| anyhow!("interval end before start"))?
            as usize;

        let chrom_name =
            String::from_utf8_lossy(self.reader.header().tid2name(chrom_tid))
                .to_string();

        self.reader.fetch(FetchDefinition::Region(
            chrom_tid as i32,
            start_pos as i64,
            end_pos as i64,
        ))?;
        let records = self
            .reader
            .records()
            .filter_ok(|record| record.tid() >= 0i32)
            .filter_ok(|record| record_is_primary(record))
            // TODO: Capture errors here. Also check for partition-tag
            .filter_map(|res| res.ok());

        let mut chrom_features: FxHashMap<PositionStrand, Tally2> =
            FxHashMap::default();
        let mut add_to_tally = |position_strand: PositionStrand,
                                call: Call,
                                motif_info: Option<&MotifInfo>,
                                motif_idxs: u8| {
            let call = match motif_info.map(|x| x.primary_base) {
                Some(x) if call.matches_dna_base(&x) => call,
                Some(_x) => match call {
                    Call::Delete => call,
                    Call::Filtered(primary_base)
                    | Call::NoCall(primary_base)
                    | Call::CanonicalCall { primary_base, .. }
                    | Call::ModifiedCall { primary_base, .. } => {
                        Call::NoCall(primary_base)
                    }
                },
                None => call,
            };
            let position_strand = match motif_info {
                Some(motif_info) if self.combine_strands => {
                    let (rpos, strand) = position_strand;
                    if strand == Strand::Negative {
                        (
                            rpos.saturating_sub(
                                (motif_info.reverse_offset
                                    - motif_info.forward_offset)
                                    as u32,
                            ),
                            Strand::Positive,
                        )
                    } else {
                        position_strand
                    }
                }
                _ => position_strand,
            };

            chrom_features
                .entry(position_strand)
                .or_insert_with(Tally2::default)
                .add_call(call, motif_idxs, self.max_depth);
        };

        let mut mod_pos = Option::<usize>::None;
        let mut canonical_base = Option::<DnaBase>::None;
        let mut pos_base_mod_call = Option::<BaseModCall>::None;
        let mut mod_strand = Option::<Strand>::None;
        'records: for record in records {
            let reverse = record.is_reverse();
            let record_strand = if record.is_reverse() {
                Strand::Negative
            } else {
                Strand::Positive
            };
            let Ok(mut modbase_iter) = record.basemods_position_iter() else {
                erred_records = erred_records.saturating_add(1);
                continue 'records;
            };

            let mm_tag_infos = MmTagInfo::from_record(&record)?;
            let (bases_to_codes, implicit_bases) =
                get_base_codes_and_implicits(&mm_tag_infos);

            let Ok(_) = update_mods_iter2(
                &mut modbase_iter,
                &mut mod_pos,
                &mut canonical_base,
                &mut pos_base_mod_call,
                &mut mod_strand,
                reverse,
                &self.caller,
            ) else {
                erred_records = erred_records.saturating_add(1);
                continue 'records;
            };

            let aligned_pairs_iter = record
                .aligned_pairs_full()
                .filter_map(|pair| {
                    let qpos = pair[0].and_then(|x| {
                        if x < 0 {
                            erred_records = erred_records.saturating_add(1);
                            return None;
                        } else {
                            Some(x as usize)
                        }
                    });
                    let r_pos = pair[1];
                    match r_pos {
                        None => None,
                        Some(pos) if pos < 0 => None,
                        Some(pos) => {
                            let pos = pos as u32;
                            if pos >= start_pos && pos < end_pos {
                                Some((
                                    qpos,
                                    pos.checked_sub(start_pos).unwrap(),
                                ))
                            } else {
                                None
                            }
                        }
                    }
                })
                .filter_map(|(qpos, rpos)| {
                    match &chrom_coordinates.focus_positions {
                        FocusPositions2::MotifMask { mask, num_motifs } => {
                            let st = (rpos * 2 * (*num_motifs as u32))
                                + (*num_motifs * reverse as u8) as u32;
                            let st = st as usize;
                            let end = st + (*num_motifs as usize);
                            let bs = &mask[st..end];
                            debug_assert_eq!(
                                bs.len(),
                                *num_motifs as usize,
                                "bs:{bs:?} st:{st} end:{end} \
                                 num_motifs:{num_motifs}"
                            );
                            bs.first_one().map(|idx| {
                                (
                                    qpos,
                                    rpos,
                                    Some(&self.motif_bases[idx]),
                                    indices_to_byte(bs),
                                )
                            })
                        }
                        FocusPositions2::AllPositions => {
                            Some((qpos, rpos, Option::<&MotifInfo>::None, 0u8))
                        }
                        FocusPositions2::MaskedPositions { mask: _ } => {
                            unreachable!()
                        }
                    }
                });
            'pileup: for (qpos, rpos, ref_base, motif_idxs) in
                aligned_pairs_iter
            {
                match (qpos, mod_pos) {
                    (Some(q), Some(mp)) if q == mp => {
                        debug_assert!(canonical_base.is_some());
                        debug_assert!(pos_base_mod_call.is_some());
                        debug_assert!(mod_strand.is_some());

                        let read_primary_base = canonical_base.unwrap();
                        let call = Call::from_basemodcall(
                            pos_base_mod_call.unwrap(),
                            read_primary_base,
                            &bases_to_codes,
                        );

                        add_to_tally(
                            (rpos, mod_strand.unwrap()),
                            call,
                            ref_base,
                            motif_idxs,
                        );

                        let Ok(_) = update_mods_iter2(
                            &mut modbase_iter,
                            &mut mod_pos,
                            &mut canonical_base,
                            &mut pos_base_mod_call,
                            &mut mod_strand,
                            reverse,
                            &self.caller,
                        ) else {
                            erred_records = erred_records.saturating_add(1);
                            continue 'records;
                        };
                        continue 'pileup;
                    }
                    (Some(q), Some(mp)) if q < mp => {
                        let base = {
                            let Ok(tmp) = DnaBase::try_from(record.seq()[q])
                            else {
                                erred_records = erred_records.saturating_add(1);
                                continue 'records;
                            };
                            if record.is_reverse() {
                                tmp.complement()
                            } else {
                                tmp
                            }
                        };
                        if implicit_bases.contains(&base) {
                            add_to_tally(
                                (rpos, record_strand),
                                Call::CanonicalCall {
                                    primary_base: base,
                                    mod_codes: bases_to_codes
                                        .get(&base)
                                        .cloned()
                                        .unwrap_or_else(Vec::new),
                                },
                                ref_base,
                                motif_idxs,
                            );
                        } else {
                            add_to_tally(
                                (rpos, record_strand),
                                Call::NoCall(base),
                                ref_base,
                                motif_idxs,
                            );
                        }
                    }
                    (Some(q), Some(mp)) if q > mp => {
                        mod_pos = None;
                        canonical_base = None;
                        pos_base_mod_call = None;
                        mod_strand = None;
                        'overran: while let Some(r) = modbase_iter.next() {
                            let Ok((pos, codes)) = r else {
                                erred_records = erred_records.saturating_add(1);
                                continue 'records;
                            };
                            let pos = pos as usize;
                            if q > pos {
                                continue 'overran;
                            } else {
                                let probs = codes
                                    .iter()
                                    .map(|hts_code| {
                                        let mod_code = ModCodeRepr::from(
                                            hts_code.modified_base,
                                        );
                                        let prob = qual_to_prob(hts_code.qual);
                                        (mod_code, prob)
                                    })
                                    .collect();
                                let Ok(can_base) =
                                    DnaBase::try_from(codes[0].canonical_base)
                                else {
                                    erred_records =
                                        erred_records.saturating_add(1);
                                    continue 'records;
                                };
                                let Ok(mod_tag_strand) =
                                    Strand::try_from(codes[0].strand)
                                else {
                                    erred_records =
                                        erred_records.saturating_add(1);
                                    continue 'records;
                                };
                                let can_base = match mod_tag_strand {
                                    Strand::Positive => can_base,
                                    Strand::Negative => can_base.complement(),
                                };
                                let Ok(pos_mod_strand) = duplex_aware_strand(
                                    mod_tag_strand,
                                    reverse,
                                ) else {
                                    erred_records =
                                        erred_records.saturating_add(1);
                                    continue 'records;
                                };
                                let bmp = BaseModProbs::new(probs, false);
                                let call = self.caller.call(&can_base, &bmp);

                                if q == pos {
                                    add_to_tally(
                                        (rpos, pos_mod_strand),
                                        Call::from_basemodcall(
                                            call,
                                            can_base,
                                            &bases_to_codes,
                                        ),
                                        ref_base,
                                        motif_idxs,
                                    );
                                    mod_pos = Some(pos);
                                    canonical_base = Some(can_base);
                                    pos_base_mod_call = Some(call);
                                    mod_strand = Some(pos_mod_strand);
                                    break 'overran;
                                } else {
                                    assert!(pos > q);
                                    let base = {
                                        let tmp =
                                            DnaBase::try_from(record.seq()[q])
                                                .unwrap();
                                        if record.is_reverse() {
                                            tmp.complement()
                                        } else {
                                            tmp
                                        }
                                    };
                                    if implicit_bases.contains(&base) {
                                        add_to_tally(
                                            (rpos, record_strand),
                                            Call::CanonicalCall {
                                                primary_base: base,
                                                mod_codes: bases_to_codes
                                                    .get(&base)
                                                    .cloned()
                                                    .unwrap(),
                                            },
                                            ref_base,
                                            motif_idxs,
                                        );
                                    } else {
                                        add_to_tally(
                                            (rpos, record_strand),
                                            Call::NoCall(base),
                                            ref_base,
                                            motif_idxs,
                                        );
                                    }
                                    mod_pos = Some(pos);
                                    canonical_base = Some(can_base);
                                    pos_base_mod_call = Some(call);
                                    mod_strand = Some(pos_mod_strand);
                                    break 'overran;
                                }
                            }
                        }
                        continue 'pileup;
                    }
                    (None, _) => {
                        add_to_tally(
                            (rpos, record_strand),
                            Call::Delete,
                            ref_base,
                            motif_idxs,
                        );
                    }
                    (Some(q), None) => {
                        let base = {
                            let tmp =
                                DnaBase::try_from(record.seq()[q]).unwrap();
                            if record.is_reverse() {
                                tmp.complement()
                            } else {
                                tmp
                            }
                        };
                        if implicit_bases.contains(&base) {
                            add_to_tally(
                                (rpos, record_strand),
                                Call::CanonicalCall {
                                    primary_base: base,
                                    mod_codes: bases_to_codes
                                        .get(&base)
                                        .cloned()
                                        .unwrap(),
                                },
                                ref_base,
                                motif_idxs,
                            );
                        } else {
                            add_to_tally(
                                (rpos, record_strand),
                                Call::NoCall(base),
                                ref_base,
                                motif_idxs,
                            );
                        }
                    }
                    (Some(_), Some(_)) => unreachable!(),
                }
            }
        }

        let combine_mods = match self.pileup_numeric_options {
            PileupNumericOptions::Combine => true,
            _ => false,
        };
        let position_feature_counts = chrom_features
            .into_iter()
            .sorted_by(|((x, a), _), ((y, b), _)| match x.cmp(y) {
                Ordering::Equal => a.cmp(b),
                o @ _ => o,
            })
            .flat_map(|((ref_pos, strand), tally)| {
                tally.into_counts(
                    start_pos,
                    ref_pos,
                    combine_mods,
                    if self.combine_strands { '.' } else { strand.to_char() },
                )
            })
            .filter(|x| x.is_valid())
            .collect::<Vec<PileupFeatureCounts2>>();

        pileup_space.chrom_name = chrom_name;
        pileup_space.interval_width = width;
        pileup_space.position_feature_counts = position_feature_counts;
        pileup_space.failed_records = erred_records;
        Ok(pileup_space)
    }
}

pub(super) struct DnaAllContext {}
pub(super) struct DnaCytosineCombine {}
pub(super) struct DnaCpGCombineStrands {}
pub(super) struct Dynamic {}

pub(super) const DNA_N_FEATURES: usize = 10usize;
// G | T counts are offset 0
pub(super) const DELETE_OFFSET: usize = 1usize;
pub(super) const NO_CALL_OFFSET: usize = 2usize;
pub(super) const CAN_C_OFFSET: usize = 3usize;
pub(super) const FILT_C_OFFSET: usize = 4usize;
pub(super) const METH_C_OFFSET: usize = 5usize;
pub(super) const OTHER_C_OFFSET: usize = 6usize;
pub(super) const FILT_A_OFFSET: usize = 7usize;
pub(super) const CAN_A_OFFSET: usize = 8usize;
pub(super) const METH_A_OFFSET: usize = 9usize;
// You don't need a NO_CALL per base since we only count matches to the
// reference. pub(super) const NO_CALL_A_OFFSET: usize = 10usize;

pub(super) const DNA_N_FEATURES_C_COMBINE: usize = 6usize;
pub(super) const DNA_BASES_CYTOSINE_FIRST: [DnaBase; 4] =
    [DnaBase::C, DnaBase::A, DnaBase::T, DnaBase::G];

pub(super) const HPTAG: SamTag = SamTag { inner: [b'H', b'P'] };
const PHASED_OFFSETS: &'static [usize] = &[0usize, 1, 2];
const UNPHASED_OFFSET: &'static [usize] = &[0usize];

// Dyn Offsets
// DIff = 0
// Delete = 1
// NOCall = 2
const DYN_FILTERED: usize = 3usize;
const DYN_CAN_A: usize = 4usize;
const DYN_CAN_C: usize = 5usize;
const DYN_CAN_G: usize = 6usize;
const DYN_CAN_T: usize = 7usize;
const DYN_OTHER_MOD_A: usize = 8usize;
const DYN_OTHER_MOD_C: usize = 9usize;
const DYN_OTHER_MOD_G: usize = 10usize;
const DYN_OTHER_MOD_T: usize = 11usize;
const DYN_N_CONSTANT_COUNTS: usize = 12usize;

fn dynamic_canonical_offset(canonical_base: DnaBase) -> usize {
    match canonical_base {
        DnaBase::A => DYN_CAN_A,
        DnaBase::C => DYN_CAN_C,
        DnaBase::G => DYN_CAN_G,
        DnaBase::T => DYN_CAN_T,
    }
}

fn dynamic_other_mod_offset(canonical_base: DnaBase) -> usize {
    match canonical_base {
        DnaBase::A => DYN_OTHER_MOD_A,
        DnaBase::C => DYN_OTHER_MOD_C,
        DnaBase::G => DYN_OTHER_MOD_G,
        DnaBase::T => DYN_OTHER_MOD_T,
    }
}

#[inline]
fn dynamic_anchor_rpos(rpos: u32, reverse: bool, motif_offset: u32) -> u32 {
    if reverse {
        rpos.saturating_sub(motif_offset)
    } else {
        rpos
    }
}

fn dynamic_mod_offset(
    mod_codes: &[(DnaBase, ModCodeRepr)],
    canonical_base: DnaBase,
    mod_code: ModCodeRepr,
    combine_mods: bool,
) -> usize {
    let compact_slot = if combine_mods {
        mod_codes.iter().position(|(base, _)| *base == canonical_base)
    } else {
        mod_codes.iter().position(|(base, code)| {
            *base == canonical_base && *code == mod_code
        })
    };
    let mod_offset = compact_slot
        .map(|slot| DYN_N_CONSTANT_COUNTS + slot)
        .unwrap_or_else(|| dynamic_other_mod_offset(canonical_base));
    let row_stride = DYN_N_CONSTANT_COUNTS + mod_codes.len();
    assert!(
        mod_offset < row_stride,
        "dynamic modification offset {mod_offset} exceeds row stride \
         {row_stride}"
    );
    mod_offset
}

pub(super) trait ACountsMatrix<T, const STRANDS: usize, const CHECK_DEPTH: bool>
{
    fn new(
        width: usize,
        phased: bool,
        mod_codes: Vec<(DnaBase, ModCodeRepr)>,
        motif_offset: u32,
        max_depth: u16,
    ) -> Self;
    fn reset(&mut self, width: usize, phased: bool);
    #[inline]
    fn incr_call(
        &mut self,
        mod_state: ModState,
        rpos: u32,
        reference_base: DnaBase,
        reverse: bool,
        haplotype: u8,
        combine_mods: bool,
    ) {
        if mod_state.primary_base != reference_base {
            self.incr_diff_call(
                rpos,
                mod_state.primary_base,
                reference_base,
                reverse,
                haplotype,
            );
            return;
        }

        if CHECK_DEPTH {
            if self.reached_max_depth(rpos, reverse, haplotype) {
                return;
            }
        }

        if mod_state.filtered {
            self.incr_filtered_count(
                rpos,
                mod_state.primary_base,
                reverse,
                haplotype,
            );
        } else if mod_state.modified {
            debug_assert!(!mod_state.filtered);
            self.incr_modified_count(
                rpos,
                mod_state.primary_base,
                mod_state.mod_code,
                reverse,
                haplotype,
                combine_mods,
            )
        } else {
            debug_assert!(!mod_state.filtered);
            self.incr_canonical_count(
                rpos,
                mod_state.primary_base,
                reference_base,
                reverse,
                haplotype,
            );
        }
    }
    fn incr_delete(&mut self, rpos: u32, reverse: bool, haplotype: u8);
    fn incr_canonical_count(
        &mut self,
        rpos: u32,
        base: DnaBase,
        reference_base: DnaBase,
        reverse: bool,
        haplotype: u8,
    );
    fn incr_modified_count(
        &mut self,
        rpos: u32,
        canonical_base: DnaBase,
        mod_code: ModCodeRepr,
        reverse: bool,
        haplotype: u8,
        combine_mods: bool,
    );
    fn incr_filtered_count(
        &mut self,
        rpos: u32,
        canonical_base: DnaBase,
        reverse: bool,
        haplotype: u8,
    );
    fn incr_diff_call(
        &mut self,
        rpos: u32,
        base: DnaBase,
        ref_base: DnaBase,
        reverse: bool,
        haplotype: u8,
    );
    fn stride(&self, dna_mod_option: DnaModOption) -> usize;
    fn strand_width(&self) -> usize;
    fn reached_max_depth(
        &self,
        rpos: u32,
        reverse: bool,
        haplotype: u8,
    ) -> bool;
    fn decode(
        &self,
        start_pos: u32,
        mod_code_op: DnaModOption,
        phased: bool,
    ) -> Vec<PileupFeatureCounts2> {
        let offsets = if phased { PHASED_OFFSETS } else { UNPHASED_OFFSET };
        offsets
            .iter()
            .flat_map(|&offset| {
                let pos_start = offset * (self.strand_width() * STRANDS);
                let pos_end = pos_start + self.strand_width();
                let neg_start = pos_end;
                let neg_end = neg_start + self.strand_width();
                self.decode_inner(
                    start_pos,
                    mod_code_op,
                    pos_start..pos_end,
                    neg_start..neg_end,
                )
            })
            .collect()
    }

    fn decode_inner(
        &self,
        start_pos: u32,
        mod_code_op: DnaModOption,
        pos_range: Range<usize>,
        neg_range: Range<usize>,
    ) -> Vec<PileupFeatureCounts2>;
}

pub(super) struct CountsMatrix {
    strand_width: usize,
    inner: Vec<u16>,
    mod_codes: Vec<(DnaBase, ModCodeRepr)>,
    motif_offset: u32,
    max_depth: u16,
}

impl ACountsMatrix<DnaAllContext, 2, false> for CountsMatrix {
    fn new(
        width: usize,
        phased: bool,
        mod_codes: Vec<(DnaBase, ModCodeRepr)>,
        motif_offset: u32,
        max_depth: u16,
    ) -> Self {
        let strand_width = width * DNA_N_FEATURES;
        let size = if phased {
            strand_width * 2 * 3 // 2 for each strand, x3 for HP1, HP2, unphased
        } else {
            strand_width * 2 // 2 for each strand
        };
        let inner = vec![0u16; size];
        Self { strand_width, inner, mod_codes, motif_offset, max_depth }
    }

    fn reset(&mut self, width: usize, phased: bool) {
        self.inner.iter_mut().for_each(|x| *x = 0u16);
        let strand_width = width * DNA_N_FEATURES;
        let size = if phased {
            strand_width * 2 * 3 // 2 for each strand, x3 for HP1, HP2, unphased
        } else {
            strand_width * 2 // 2 for each strand
        };
        self.strand_width = strand_width;
        self.inner.resize(size, 0u16);
    }

    #[inline]
    fn incr_delete(&mut self, rpos: u32, reverse: bool, haplotype: u8) {
        let offset = calc_offset::<DNA_N_FEATURES, 2>(
            rpos,
            self.strand_width,
            reverse,
            haplotype,
        );
        self.inner[offset + DELETE_OFFSET] =
            self.inner[offset + DELETE_OFFSET].saturating_add(1);
    }

    #[inline]
    fn incr_canonical_count(
        &mut self,
        rpos: u32,
        base: DnaBase,
        reference_base: DnaBase,
        reverse: bool,
        haplotype: u8,
    ) {
        if base != reference_base {
            return;
        }
        let offset = calc_offset::<DNA_N_FEATURES, 2>(
            rpos,
            self.strand_width,
            reverse,
            haplotype,
        );
        match base {
            DnaBase::A => {
                self.inner[offset + CAN_A_OFFSET] =
                    self.inner[offset + CAN_A_OFFSET].saturating_add(1);
            }
            DnaBase::C => {
                self.inner[offset + CAN_C_OFFSET] =
                    self.inner[offset + CAN_C_OFFSET].saturating_add(1);
            }
            DnaBase::G => {}
            DnaBase::T => {}
        }
    }

    #[inline]
    fn incr_modified_count(
        &mut self,
        rpos: u32,
        base: DnaBase,
        mod_code_repr: ModCodeRepr,
        reverse: bool,
        haplotype: u8,
        combine_mods: bool,
    ) {
        let offset = calc_offset::<DNA_N_FEATURES, 2>(
            rpos,
            self.strand_width,
            reverse,
            haplotype,
        );
        match base {
            DnaBase::A => {
                if combine_mods || mod_code_repr == SIX_METHYL_ADENINE {
                    self.inner[offset + METH_A_OFFSET] =
                        self.inner[offset + METH_A_OFFSET].saturating_add(1);
                } else {
                    self.inner[offset + OTHER_C_OFFSET] =
                        self.inner[offset + OTHER_C_OFFSET].saturating_add(1);
                }
            }
            DnaBase::C => {
                if combine_mods || mod_code_repr == METHYL_CYTOSINE {
                    self.inner[offset + METH_C_OFFSET] =
                        self.inner[offset + METH_C_OFFSET].saturating_add(1);
                } else {
                    self.inner[offset + OTHER_C_OFFSET] =
                        self.inner[offset + OTHER_C_OFFSET].saturating_add(1);
                }
            }
            DnaBase::G => {}
            DnaBase::T => {}
        }
    }

    #[inline]
    fn incr_filtered_count(
        &mut self,
        rpos: u32,
        base: DnaBase,
        reverse: bool,
        haplotype: u8,
    ) {
        let offset = calc_offset::<DNA_N_FEATURES, 2>(
            rpos,
            self.strand_width,
            reverse,
            haplotype,
        );
        match base {
            DnaBase::A => {
                self.inner[offset + FILT_A_OFFSET] =
                    self.inner[offset + FILT_A_OFFSET].saturating_add(1);
            }
            DnaBase::C => {
                self.inner[offset + FILT_C_OFFSET] =
                    self.inner[offset + FILT_C_OFFSET].saturating_add(1);
            }
            DnaBase::G => {}
            DnaBase::T => {}
        }
    }

    #[inline]
    fn incr_diff_call(
        &mut self,
        rpos: u32,
        base: DnaBase,
        reference_base: DnaBase,
        reverse: bool,
        haplotype: u8,
    ) {
        let offset = calc_offset::<DNA_N_FEATURES, 2>(
            rpos,
            self.strand_width,
            reverse,
            haplotype,
        );
        if base == reference_base {
            self.inner[offset + NO_CALL_OFFSET] =
                self.inner[offset + NO_CALL_OFFSET].saturating_add(1);
        } else {
            self.inner[offset] = self.inner[offset].saturating_add(1);
        }
    }

    fn strand_width(&self) -> usize {
        self.strand_width
    }

    #[inline]
    fn reached_max_depth(
        &self,
        _rpos: u32,
        _reverse: bool,
        _haplotype: u8,
    ) -> bool {
        unreachable!()
    }

    fn stride(&self, dna_mod_option: DnaModOption) -> usize {
        match dna_mod_option {
            DnaModOption::Other(_) => 3,
            DnaModOption::Combine | DnaModOption::Single => 2,
        }
    }

    #[inline]
    fn decode_inner(
        &self,
        start_pos: u32,
        dna_mod_option: DnaModOption,
        pos_range: Range<usize>,
        neg_range: Range<usize>,
    ) -> Vec<PileupFeatureCounts2> {
        let pos_strand = &self.inner[pos_range];
        let neg_strand = &self.inner[neg_range];
        assert_eq!(pos_strand.len(), neg_strand.len());

        let (mod_code_op, methyl_code) = match dna_mod_option {
            DnaModOption::Other(mod_code_repr) => {
                (Some(mod_code_repr), METHYL_CYTOSINE)
            }
            DnaModOption::Combine => (None, ANY_CYTOSINE),
            DnaModOption::Single => (None, METHYL_CYTOSINE),
        };
        slice_to_counts_dna_all_context(
            pos_strand,
            start_pos,
            '+',
            methyl_code,
            mod_code_op,
        )
        .zip(slice_to_counts_dna_all_context(
            neg_strand,
            start_pos,
            '-',
            methyl_code,
            mod_code_op,
        ))
        .flat_map(|(pos, neg)| [pos, neg])
        .collect()
    }
}

impl ACountsMatrix<DnaCytosineCombine, 2, false> for CountsMatrix {
    fn new(
        width: usize,
        phased: bool,
        mod_codes: Vec<(DnaBase, ModCodeRepr)>,
        motif_offset: u32,
        max_depth: u16,
    ) -> Self {
        let strand_width = width * DNA_N_FEATURES_C_COMBINE;
        let size = if phased {
            strand_width * 2 * 3 // 2 for each strand, x3 for HP1, HP2, unphased
        } else {
            strand_width * 2 // 2 for each strand
        };
        let inner = vec![0u16; size];
        Self { strand_width, inner, mod_codes, motif_offset, max_depth }
    }

    fn reset(&mut self, width: usize, phased: bool) {
        let strand_width = width * DNA_N_FEATURES_C_COMBINE;
        self.inner.iter_mut().for_each(|x| *x = 0u16);
        let size = if phased {
            strand_width * 2 * 3 // 2 for each strand, x3 for HP1, HP2, unphased
        } else {
            strand_width * 2 // 2 for each strand
        };
        self.strand_width = strand_width;
        self.inner.resize(size, 0u16);
    }

    #[inline]
    fn incr_delete(&mut self, rpos: u32, reverse: bool, haplotype: u8) {
        let offset = calc_offset::<DNA_N_FEATURES_C_COMBINE, 2>(
            rpos,
            self.strand_width,
            reverse,
            haplotype,
        );
        self.inner[offset + DELETE_OFFSET] =
            self.inner[offset + DELETE_OFFSET].saturating_add(1);
    }

    #[inline]
    fn incr_canonical_count(
        &mut self,
        rpos: u32,
        base: DnaBase,
        _reference_base: DnaBase,
        reverse: bool,
        haplotype: u8,
    ) {
        let offset = calc_offset::<DNA_N_FEATURES_C_COMBINE, 2>(
            rpos,
            self.strand_width,
            reverse,
            haplotype,
        );
        match base {
            DnaBase::C => {
                self.inner[offset + CAN_C_OFFSET] =
                    self.inner[offset + CAN_C_OFFSET].saturating_add(1);
            }
            // Increase diff count
            _ => self.inner[offset] = self.inner[offset].saturating_add(1),
        }
    }

    #[inline]
    fn incr_modified_count(
        &mut self,
        rpos: u32,
        canonical_base: DnaBase,
        _mod_code: ModCodeRepr,
        reverse: bool,
        haplotype: u8,
        _combine_mods: bool,
    ) {
        let offset = calc_offset::<DNA_N_FEATURES_C_COMBINE, 2>(
            rpos,
            self.strand_width,
            reverse,
            haplotype,
        );
        match canonical_base {
            DnaBase::C => {
                self.inner[offset + METH_C_OFFSET] =
                    self.inner[offset + METH_C_OFFSET].saturating_add(1);
            }
            _ => {}
        }
    }

    #[inline]
    fn incr_filtered_count(
        &mut self,
        rpos: u32,
        canonical_base: DnaBase,
        reverse: bool,
        haplotype: u8,
    ) {
        let offset = calc_offset::<DNA_N_FEATURES_C_COMBINE, 2>(
            rpos,
            self.strand_width,
            reverse,
            haplotype,
        );
        match canonical_base {
            DnaBase::C => {
                self.inner[offset + FILT_C_OFFSET] =
                    self.inner[offset + FILT_C_OFFSET].saturating_add(1);
            }
            _ => {}
        }
    }

    fn strand_width(&self) -> usize {
        self.strand_width
    }

    fn reached_max_depth(
        &self,
        _rpos: u32,
        _reverse: bool,
        _haplotype: u8,
    ) -> bool {
        unreachable!()
    }

    fn stride(&self, _dna_mod_option: DnaModOption) -> usize {
        1
    }

    #[inline]
    fn incr_diff_call(
        &mut self,
        rpos: u32,
        base: DnaBase,
        _ref_base: DnaBase,
        reverse: bool,
        haplotype: u8,
    ) {
        let offset = calc_offset::<DNA_N_FEATURES_C_COMBINE, 2>(
            rpos,
            self.strand_width,
            reverse,
            haplotype,
        );
        match base {
            DnaBase::C => {
                self.inner[offset + NO_CALL_OFFSET] =
                    self.inner[offset + NO_CALL_OFFSET].saturating_add(1)
            }
            _ => self.inner[offset] = self.inner[offset].saturating_add(1),
        }
    }

    #[inline]
    fn decode_inner(
        &self,
        start_pos: u32,
        _mod_code_op: DnaModOption,
        pos_range: Range<usize>,
        neg_range: Range<usize>,
    ) -> Vec<PileupFeatureCounts2> {
        let pos_strand = &self.inner[pos_range];
        let neg_strand = &self.inner[neg_range];
        assert_eq!(pos_strand.len(), neg_strand.len());

        slice_to_counts_dna_c_combine(pos_strand, start_pos, '+')
            .zip(slice_to_counts_dna_c_combine(neg_strand, start_pos, '-'))
            .flat_map(|(pos, neg)| [pos, neg])
            .collect()
    }
}

impl ACountsMatrix<DnaCpGCombineStrands, 1, false> for CountsMatrix {
    fn new(
        width: usize,
        phased: bool,
        mod_codes: Vec<(DnaBase, ModCodeRepr)>,
        motif_offset: u32,
        max_depth: u16,
    ) -> Self {
        let strand_width = width * DNA_N_FEATURES;
        let size = if phased {
            strand_width * 3 // 2 for each strand, x3 for HP1, HP2, unphased
        } else {
            strand_width
        };
        let inner = vec![0u16; size];
        Self { strand_width, inner, mod_codes, motif_offset, max_depth }
    }

    fn reset(&mut self, width: usize, phased: bool) {
        self.inner.iter_mut().for_each(|x| *x = 0u16);
        let strand_width = width * DNA_N_FEATURES;
        let size = if phased {
            strand_width * 3 // x3 for HP1, HP2, unphased
        } else {
            strand_width
        };
        self.strand_width = strand_width;
        self.inner.resize(size, 0u16);
    }

    #[inline]
    fn incr_delete(&mut self, rpos: u32, reverse: bool, haplotype: u8) {
        let rpos = if reverse { rpos.saturating_sub(1) } else { rpos };
        let offset = calc_offset::<DNA_N_FEATURES, 1>(
            rpos,
            self.strand_width,
            false,
            haplotype,
        );
        self.inner[offset + DELETE_OFFSET] =
            self.inner[offset + DELETE_OFFSET].saturating_add(1);
    }

    #[inline]
    fn incr_canonical_count(
        &mut self,
        rpos: u32,
        base: DnaBase,
        _reference_base: DnaBase,
        reverse: bool,
        haplotype: u8,
    ) {
        let rpos = if reverse { rpos.saturating_sub(1) } else { rpos };
        let offset = calc_offset::<DNA_N_FEATURES, 1>(
            rpos,
            self.strand_width,
            false,
            haplotype,
        );
        match base {
            DnaBase::C => {
                assert!(
                    (offset + METH_C_OFFSET) < self.inner.len(),
                    "[incr_canonical_count] offset: {offset}, rpos: {rpos}, \
                     strand_width: {}, reverse: {reverse}, haplotype: \
                     {haplotype}",
                    self.strand_width
                );
                self.inner[offset + CAN_C_OFFSET] =
                    self.inner[offset + CAN_C_OFFSET].saturating_add(1);
            }
            _ => {}
        }
    }

    #[inline]
    fn incr_modified_count(
        &mut self,
        rpos: u32,
        canonical_base: DnaBase,
        mod_code: ModCodeRepr,
        reverse: bool,
        haplotype: u8,
        combine_mods: bool,
    ) {
        let rpos = if reverse { rpos.saturating_sub(1) } else { rpos };
        let offset = calc_offset::<DNA_N_FEATURES, 1>(
            rpos,
            self.strand_width,
            false,
            haplotype,
        );
        match canonical_base {
            DnaBase::C => {
                if mod_code == METHYL_CYTOSINE || combine_mods {
                    debug_assert!(
                        (offset + METH_C_OFFSET) < self.inner.len(),
                        "[incr_modified_count] offset: {offset}, rpos: \
                         {rpos}, strand_width: {}, reverse: {reverse}, \
                         haplotype: {haplotype}",
                        self.strand_width
                    );
                    self.inner[offset + METH_C_OFFSET] =
                        self.inner[offset + METH_C_OFFSET].saturating_add(1);
                } else {
                    self.inner[offset + OTHER_C_OFFSET] =
                        self.inner[offset + OTHER_C_OFFSET].saturating_add(1);
                }
            }
            _ => {}
        }
    }

    #[inline]
    fn incr_filtered_count(
        &mut self,
        rpos: u32,
        canonical_base: DnaBase,
        reverse: bool,
        haplotype: u8,
    ) {
        let rpos = if reverse { rpos.saturating_sub(1) } else { rpos };
        let offset = calc_offset::<DNA_N_FEATURES, 1>(
            rpos,
            self.strand_width,
            false,
            haplotype,
        );
        match canonical_base {
            DnaBase::C => {
                self.inner[offset + FILT_C_OFFSET] =
                    self.inner[offset + FILT_C_OFFSET].saturating_add(1);
            }
            _ => {}
        }
    }

    #[inline]
    fn incr_diff_call(
        &mut self,
        rpos: u32,
        base: DnaBase,
        _ref_base: DnaBase,
        reverse: bool,
        haplotype: u8,
    ) {
        let rpos = if reverse { rpos.saturating_sub(1) } else { rpos };
        let offset = calc_offset::<DNA_N_FEATURES, 1>(
            rpos,
            self.strand_width,
            false,
            haplotype,
        );
        match base {
            DnaBase::C => {
                self.inner[offset + NO_CALL_OFFSET] =
                    self.inner[offset + NO_CALL_OFFSET].saturating_add(1);
            }
            _ => {
                self.inner[offset] = self.inner[offset].saturating_add(1);
            }
        }
    }

    fn strand_width(&self) -> usize {
        self.strand_width
    }

    fn reached_max_depth(
        &self,
        _rpos: u32,
        _reverse: bool,
        _haplotype: u8,
    ) -> bool {
        unreachable!()
    }

    fn stride(&self, dna_mod_option: DnaModOption) -> usize {
        match dna_mod_option {
            DnaModOption::Other(_) => 3,
            DnaModOption::Combine | DnaModOption::Single => 2,
        }
    }

    #[inline]
    fn decode_inner(
        &self,
        start_pos: u32,
        mod_code_op: DnaModOption,
        pos_range: Range<usize>,
        _neg_range: Range<usize>,
    ) -> Vec<PileupFeatureCounts2> {
        let (mod_code_op, methyl_code) = match mod_code_op {
            DnaModOption::Other(mod_code_repr) => {
                (Some(mod_code_repr), METHYL_CYTOSINE)
            }
            DnaModOption::Single => (None, METHYL_CYTOSINE),
            DnaModOption::Combine => (None, ANY_CYTOSINE),
        };
        slice_to_counts_dna_all_context(
            &self.inner[pos_range],
            start_pos,
            '.',
            methyl_code,
            mod_code_op,
        )
        .collect()
    }
}

impl<const STRANDS: usize, const CHECK_DEPTH: bool>
    ACountsMatrix<Dynamic, STRANDS, CHECK_DEPTH> for CountsMatrix
{
    fn new(
        width: usize,
        phased: bool,
        mod_codes: Vec<(DnaBase, ModCodeRepr)>,
        motif_offset: u32,
        max_depth: u16,
    ) -> Self {
        // N <- n_modified_bases
        // [N_diff, N_delete, N_nocall, N_canonical{A,C,G,T},
        // N_filtered{A,C,G,T}, N_mod{**mod_codes},
        // N_nocall{len(mod_codes)}]
        let strand_width = (DYN_N_CONSTANT_COUNTS + mod_codes.len()) * width;
        let size = if phased {
            strand_width * STRANDS * 3 // *STRANDS for each strand, x3 for HP1,
                                       // HP2, unphased
        } else {
            strand_width * STRANDS // *STRANDS for each strand
        };
        let inner = vec![0u16; size];
        Self { strand_width, inner, mod_codes, motif_offset, max_depth }
    }

    fn reset(&mut self, width: usize, phased: bool) {
        self.inner.iter_mut().for_each(|x| *x = 0u16);
        let strand_width =
            (DYN_N_CONSTANT_COUNTS + self.mod_codes.len()) * width;
        let size = if phased {
            strand_width * STRANDS * 3 // *STRANDS for each strand, x3 for HP1,
                                       // HP2, unphased
        } else {
            strand_width * STRANDS // *STRANDS for each strand
        };
        self.strand_width = strand_width;
        self.inner.resize(size, 0u16);
    }

    fn incr_delete(&mut self, rpos: u32, reverse: bool, haplotype: u8) {
        if <CountsMatrix as ACountsMatrix<Dynamic, STRANDS, CHECK_DEPTH>>::reached_max_depth(self, rpos, reverse, haplotype) {
            return;
        }

        let rpos = dynamic_anchor_rpos(rpos, reverse, self.motif_offset);
        let offset = calc_offset_dyn::<STRANDS>(
            rpos,
            self.strand_width,
            DYN_N_CONSTANT_COUNTS + self.mod_codes.len(),
            reverse,
            haplotype,
        );
        assert!(offset + DELETE_OFFSET < self.inner.len());
        self.inner[offset + DELETE_OFFSET] =
            self.inner[offset + DELETE_OFFSET].saturating_add(1);
    }

    fn incr_canonical_count(
        &mut self,
        rpos: u32,
        base: DnaBase,
        reference_base: DnaBase,
        reverse: bool,
        haplotype: u8,
    ) {
        if base != reference_base {
            return;
        }
        let rpos = dynamic_anchor_rpos(rpos, reverse, self.motif_offset);
        let offset = calc_offset_dyn::<STRANDS>(
            rpos,
            self.strand_width,
            DYN_N_CONSTANT_COUNTS + self.mod_codes.len(),
            reverse,
            haplotype,
        );
        let base_offset = dynamic_canonical_offset(base);
        assert!(offset + base_offset < self.inner.len());
        self.inner[offset + base_offset] =
            self.inner[offset + base_offset].saturating_add(1);
    }

    fn incr_modified_count(
        &mut self,
        rpos: u32,
        canonical_base: DnaBase,
        mod_code: ModCodeRepr,
        reverse: bool,
        haplotype: u8,
        combine_mods: bool,
    ) {
        let rpos = dynamic_anchor_rpos(rpos, reverse, self.motif_offset);
        let offset = calc_offset_dyn::<STRANDS>(
            rpos,
            self.strand_width,
            DYN_N_CONSTANT_COUNTS + self.mod_codes.len(),
            reverse,
            haplotype,
        );
        let mod_offset = dynamic_mod_offset(
            &self.mod_codes,
            canonical_base,
            mod_code,
            combine_mods,
        );
        assert!(offset + mod_offset < self.inner.len(), "STRANDS {STRANDS}");
        self.inner[offset + mod_offset] =
            self.inner[offset + mod_offset].saturating_add(1);
    }

    fn incr_filtered_count(
        &mut self,
        rpos: u32,
        _canonical_base: DnaBase,
        reverse: bool,
        haplotype: u8,
    ) {
        let rpos = dynamic_anchor_rpos(rpos, reverse, self.motif_offset);
        let offset = calc_offset_dyn::<STRANDS>(
            rpos,
            self.strand_width,
            DYN_N_CONSTANT_COUNTS + self.mod_codes.len(),
            reverse,
            haplotype,
        );
        assert!(offset + DYN_FILTERED < self.inner.len());
        self.inner[offset + DYN_FILTERED] =
            self.inner[offset + DYN_FILTERED].saturating_add(1);
    }

    fn incr_diff_call(
        &mut self,
        rpos: u32,
        base: DnaBase,
        ref_base: DnaBase,
        reverse: bool,
        haplotype: u8,
    ) {
        if <CountsMatrix as ACountsMatrix<Dynamic, STRANDS, CHECK_DEPTH>>::reached_max_depth(self, rpos, reverse, haplotype) {
            return;
        }
        let rpos = dynamic_anchor_rpos(rpos, reverse, self.motif_offset);
        let offset = calc_offset_dyn::<STRANDS>(
            rpos,
            self.strand_width,
            DYN_N_CONSTANT_COUNTS + self.mod_codes.len(),
            reverse,
            haplotype,
        );
        if base == ref_base {
            self.inner[offset + NO_CALL_OFFSET] =
                self.inner[offset + NO_CALL_OFFSET].saturating_add(1);
        } else {
            self.inner[offset] = self.inner[offset].saturating_add(1);
        }
    }

    fn strand_width(&self) -> usize {
        self.strand_width
    }

    fn reached_max_depth(
        &self,
        rpos: u32,
        reverse: bool,
        haplotype: u8,
    ) -> bool {
        let rpos = dynamic_anchor_rpos(rpos, reverse, self.motif_offset);
        let offset = calc_offset_dyn::<STRANDS>(
            rpos,
            self.strand_width,
            DYN_N_CONSTANT_COUNTS + self.mod_codes.len(),
            reverse,
            haplotype,
        );
        let coverage = self.inner
            [offset..(offset + DYN_N_CONSTANT_COUNTS + self.mod_codes.len())]
            .iter()
            .copied()
            .reduce(|a, b| a.saturating_add(b))
            .unwrap_or(0u16);
        coverage >= self.max_depth
    }

    fn stride(&self, _dna_mod_option: DnaModOption) -> usize {
        self.mod_codes.len()
    }

    fn decode_inner(
        &self,
        start_pos: u32,
        _mod_code_op: DnaModOption,
        pos_range: Range<usize>,
        neg_range: Range<usize>,
    ) -> Vec<PileupFeatureCounts2> {
        if STRANDS == 2 {
            let pos_strand = &self.inner[pos_range];
            let neg_strand = &self.inner[neg_range];
            assert_eq!(pos_strand.len(), neg_strand.len());
            slice_to_counts_dynamic(pos_strand, start_pos, '+', &self.mod_codes)
                .zip(slice_to_counts_dynamic(
                    neg_strand,
                    start_pos,
                    '-',
                    &self.mod_codes,
                ))
                .flat_map(|(pos, neg)| [pos, neg])
                .collect()
        } else {
            let pos_strand = &self.inner[pos_range];
            slice_to_counts_dynamic(pos_strand, start_pos, '.', &self.mod_codes)
                .collect()
        }
    }
}

fn slice_to_counts_dynamic<'a>(
    sl: &'a [u16],
    start_pos: u32,
    strand: char,
    mod_codes: &'a [(DnaBase, ModCodeRepr)],
) -> impl Iterator<Item = PileupFeatureCounts2> + 'a {
    sl.chunks_exact(DYN_N_CONSTANT_COUNTS + mod_codes.len())
        .enumerate()
        .flat_map(move |(i, chunk)| {
            let position = (i as u32).saturating_add(start_pos);
            let n_diff = chunk[0];
            let n_delete = chunk[DELETE_OFFSET];
            let n_nocall = chunk[NO_CALL_OFFSET];
            let n_filtered = chunk[DYN_FILTERED];
            mod_codes.iter().enumerate().map(
                move |(j, (primary_base, code))| {
                    let n_modified = chunk[DYN_N_CONSTANT_COUNTS + j];
                    let n_canonical =
                        chunk[dynamic_canonical_offset(*primary_base)];
                    let n_anon_modified =
                        chunk[dynamic_other_mod_offset(*primary_base)];
                    let n_other_modified = mod_codes
                        .iter()
                        .enumerate()
                        .map(|(k, (b, _))| {
                            if b == primary_base {
                                chunk[DYN_N_CONSTANT_COUNTS + k]
                            } else {
                                0
                            }
                        })
                        .reduce(|a, b| a.saturating_add(b))
                        .unwrap_or(0u16)
                        .saturating_sub(n_modified)
                        .saturating_add(n_anon_modified);
                    let filtered_coverage = n_modified
                        .saturating_add(n_canonical)
                        .saturating_add(n_other_modified);
                    PileupFeatureCounts2::new(
                        position,
                        strand,
                        filtered_coverage,
                        *code,
                        n_canonical,
                        n_modified,
                        n_other_modified,
                        n_delete,
                        n_filtered,
                        n_diff,
                        n_nocall,
                        0u8,
                    )
                },
            )
        })
}

fn slice_to_counts_dna_all_context(
    sl: &[u16],
    start_pos: u32,
    strand: char,
    methyl_code: ModCodeRepr,
    mod_code_op: Option<ModCodeRepr>,
) -> impl Iterator<Item = PileupFeatureCounts2> + '_ {
    sl.chunks_exact(DNA_N_FEATURES).enumerate().flat_map(move |(i, chunk)| {
        let position = (i as u32).saturating_add(start_pos);
        let filtered_coverage_a =
            chunk[CAN_A_OFFSET].saturating_add(chunk[METH_A_OFFSET]);
        let filtered_coverage_c = chunk[CAN_C_OFFSET]
            .saturating_add(chunk[METH_C_OFFSET])
            .saturating_add(chunk[OTHER_C_OFFSET]);
        let a_rec = PileupFeatureCounts2::new(
            position,
            strand,
            filtered_coverage_a,
            SIX_METHYL_ADENINE,
            chunk[CAN_A_OFFSET],
            chunk[METH_A_OFFSET],
            0,
            chunk[DELETE_OFFSET],
            chunk[FILT_A_OFFSET],
            chunk[0]
                .saturating_add(chunk[FILT_C_OFFSET])
                .saturating_add(chunk[CAN_C_OFFSET])
                .saturating_add(chunk[METH_C_OFFSET])
                .saturating_add(chunk[OTHER_C_OFFSET]),
            chunk[NO_CALL_OFFSET],
            0u8,
        );

        let n_not_c = chunk[0]
            .saturating_add(chunk[FILT_A_OFFSET])
            .saturating_add(chunk[CAN_A_OFFSET])
            .saturating_add(chunk[METH_A_OFFSET]);
        let meth_c_record = PileupFeatureCounts2::new(
            position,
            strand,
            filtered_coverage_c,
            methyl_code,
            chunk[CAN_C_OFFSET],
            chunk[METH_C_OFFSET],
            chunk[OTHER_C_OFFSET],
            chunk[DELETE_OFFSET],
            chunk[FILT_C_OFFSET],
            n_not_c,
            chunk[NO_CALL_OFFSET],
            0u8,
        );
        if let Some(code) = mod_code_op {
            let third_mod_c_record = PileupFeatureCounts2::new(
                position,
                strand,
                filtered_coverage_c,
                code,
                chunk[CAN_C_OFFSET],
                chunk[OTHER_C_OFFSET],
                chunk[METH_C_OFFSET],
                chunk[DELETE_OFFSET],
                chunk[FILT_C_OFFSET],
                n_not_c,
                chunk[NO_CALL_OFFSET],
                0u8,
            );
            vec![a_rec, third_mod_c_record, meth_c_record]
        } else {
            vec![a_rec, meth_c_record]
        }
    })
}

fn slice_to_counts_dna_c_combine(
    sl: &[u16],
    start_pos: u32,
    strand: char,
) -> impl Iterator<Item = PileupFeatureCounts2> + '_ {
    sl.chunks_exact(DNA_N_FEATURES_C_COMBINE).enumerate().map(
        move |(idx, chunk)| {
            let filtered_coverage_c =
                chunk[CAN_C_OFFSET].saturating_add(chunk[METH_C_OFFSET]);
            let n_not_c = chunk[0];
            let position = start_pos.saturating_add(idx as u32);

            let meth_c_record = PileupFeatureCounts2::new(
                position,
                strand,
                filtered_coverage_c,
                ANY_CYTOSINE,
                chunk[CAN_C_OFFSET],
                chunk[METH_C_OFFSET],
                0,
                chunk[DELETE_OFFSET],
                chunk[FILT_C_OFFSET],
                n_not_c,
                chunk[NO_CALL_OFFSET],
                0u8,
            );
            meth_c_record
        },
    )
}

#[inline]
fn calc_offset<const DIM: usize, const STRANDS: usize>(
    rpos: u32,
    strand_width: usize,
    reverse: bool,
    haplotype: u8,
) -> usize {
    if reverse {
        rpos as usize * DIM
            + strand_width
            // TODO: should this be STRANDS here?
            + (strand_width * 2 * haplotype as usize)
    } else {
        rpos as usize * DIM + (strand_width * STRANDS * haplotype as usize)
    }
}

#[inline]
fn calc_offset_dyn<const STRANDS: usize>(
    rpos: u32,
    strand_width: usize,
    dim: usize,
    reverse: bool,
    haplotype: u8,
) -> usize {
    if reverse && STRANDS > 1 {
        rpos as usize * dim
            + strand_width
            // TODO: should this be STRANDS here?
            + (strand_width * 2 * haplotype as usize)
    } else {
        rpos as usize * dim + (strand_width * STRANDS * haplotype as usize)
    }
}

#[derive(Clone, Debug)]
pub enum Call {
    Delete,
    Filtered(DnaBase),
    NoCall(DnaBase),
    CanonicalCall {
        primary_base: DnaBase,
        mod_codes: Vec<ModCodeRepr>,
    },
    ModifiedCall {
        primary_base: DnaBase,
        mod_code: ModCodeRepr,
        mod_codes: Vec<ModCodeRepr>,
    },
}

impl Call {
    fn from_basemodcall(
        base_mod_call: BaseModCall,
        primary_base: DnaBase,
        bases_to_codes: &FxHashMap<DnaBase, Vec<ModCodeRepr>>,
    ) -> Self {
        match base_mod_call {
            BaseModCall::Canonical(_) => Call::CanonicalCall {
                primary_base,
                mod_codes: bases_to_codes
                    .get(&primary_base)
                    .cloned()
                    .unwrap_or_else(Vec::new),
            },
            BaseModCall::Modified(_, mod_code) => Call::ModifiedCall {
                primary_base,
                mod_code,
                mod_codes: bases_to_codes
                    .get(&primary_base)
                    .cloned()
                    .unwrap_or_else(Vec::new),
            },
            BaseModCall::Filtered => Call::Filtered(primary_base),
        }
    }

    fn matches_dna_base(&self, ref_base: &DnaBase) -> bool {
        match self {
            Call::CanonicalCall { primary_base, .. }
            | Call::ModifiedCall { primary_base, .. } => {
                ref_base == primary_base
            }
            _ => true,
        }
    }
}

#[derive(Debug, Hash, Eq, PartialEq)]
enum ModifiedBase {
    UnModified,
    Modified(ModCodeRepr),
}

#[derive(Debug, Default)]
pub(super) struct Tally2 {
    total_cov: u16,
    n_delete: u16,
    n_filtered: [u16; 4],
    basecall_counts: [u16; 4],
    modcall_counts: [FxHashMap<ModifiedBase, u16>; 4],
    motif_idxs: u8,
}

impl Tally2 {
    pub(super) fn add_call(
        &mut self,
        call: Call,
        motif_idxs: u8,
        max_coverage: u16,
    ) {
        self.motif_idxs.bitor_assign(motif_idxs);
        self.total_cov = self.total_cov.saturating_add(1);
        if self.total_cov > max_coverage {
            return;
        }
        match call {
            Call::Delete => {
                self.n_delete = self.n_delete.saturating_add(1);
            }
            Call::Filtered(read_base) => {
                self.n_filtered[read_base as usize] =
                    self.n_filtered[read_base as usize].saturating_add(1);
            }
            Call::NoCall(read_base) => {
                self.basecall_counts[read_base as usize] =
                    self.basecall_counts[read_base as usize].saturating_add(1);
            }
            Call::CanonicalCall { primary_base, mod_codes } => {
                let mod_calls = &mut self.modcall_counts[primary_base as usize];
                mod_calls
                    .entry(ModifiedBase::UnModified)
                    .and_modify(|x| *x = x.saturating_add(1))
                    .or_insert(1);
                for mod_code in mod_codes {
                    mod_calls
                        .entry(ModifiedBase::Modified(mod_code))
                        .or_insert(0);
                }
            }
            Call::ModifiedCall { primary_base, mod_code, mod_codes } => {
                let mod_calls = &mut self.modcall_counts[primary_base as usize];
                mod_calls
                    .entry(ModifiedBase::Modified(mod_code))
                    .and_modify(|x| *x = x.saturating_add(1))
                    .or_insert(1);
                for mod_code in mod_codes.into_iter().filter(|x| *x != mod_code)
                {
                    mod_calls
                        .entry(ModifiedBase::Modified(mod_code))
                        .or_insert(0);
                }
            }
        }
    }

    pub(super) fn into_counts(
        self,
        start_pos: u32,
        relative_position: u32,
        combine_mods: bool,
        strand: char,
    ) -> Vec<PileupFeatureCounts2> {
        self.modcall_counts
            .into_iter()
            .enumerate()
            .filter_map(|(i, calls)| {
                let valid_coverage = calls
                    .values()
                    .copied()
                    .reduce(|a, b| a.saturating_add(b))
                    .unwrap_or(0u16);
                if valid_coverage == 0 {
                    None
                } else {
                    Some((DnaBase::try_from(i).unwrap(), calls, valid_coverage))
                }
            })
            .flat_map(|(read_base, calls, filtered_coverage)| {
                let n_modified_all = calls
                    .iter()
                    .filter_map(|(mb, c)| match mb {
                        ModifiedBase::UnModified => None,
                        ModifiedBase::Modified(_) => Some(*c),
                    })
                    .reduce(|a, b| a.saturating_add(b))
                    .unwrap_or(0u16);

                let (n_diff, n_nocall) = self
                    .basecall_counts
                    .iter()
                    .zip(&self.n_filtered)
                    .enumerate()
                    .fold((0u16, 0u16), |(n_diff, n_nocall), (b, (c, f))| {
                        if b == (read_base as usize) {
                            (n_diff, n_nocall.saturating_add(*c))
                        } else {
                            (
                                n_diff.saturating_add(*c).saturating_add(*f),
                                n_nocall,
                            )
                        }
                    });

                let n_filt = self.n_filtered[read_base as usize];
                let n_canonical =
                    calls.get(&ModifiedBase::UnModified).copied().unwrap_or(0);
                assert!(filtered_coverage > 0, "{read_base} {calls:?}");
                let calls = if combine_mods {
                    let mut tmp = FxHashMap::default();
                    tmp.insert(
                        ModifiedBase::Modified(ModCodeRepr::Code(
                            read_base.char(),
                        )),
                        n_modified_all,
                    );
                    tmp
                } else {
                    calls
                };
                let counts = calls
                    .into_iter()
                    .filter_map(|(mod_base, count)| match mod_base {
                        ModifiedBase::UnModified => None,
                        ModifiedBase::Modified(mod_code) => {
                            Some(PileupFeatureCounts2::new(
                                relative_position.saturating_add(start_pos),
                                strand,
                                filtered_coverage,
                                mod_code,
                                n_canonical,
                                count,
                                n_modified_all.saturating_sub(count),
                                self.n_delete,
                                n_filt,
                                n_diff,
                                n_nocall,
                                self.motif_idxs,
                            ))
                        }
                    })
                    .sorted_by(|a, b| a.mod_code.cmp(&b.mod_code))
                    .collect::<Vec<PileupFeatureCounts2>>();
                counts
            })
            .collect()
    }
}

type PositionStrand = (u32, Strand);

fn get_base_codes_and_implicits(
    mm_tag_infos: &[MmTagInfo],
) -> (FxHashMap<DnaBase, Vec<ModCodeRepr>>, FxHashSet<DnaBase>) {
    mm_tag_infos
        .iter()
        .flat_map(|info| {
            info.fundamental_base.expand_bases().iter().map(|primary_base| {
                (*primary_base, info.is_implicit(), &info.mod_base_codes)
            })
        })
        .fold(
            (FxHashMap::default(), FxHashSet::default()),
            |(mut agg, mut bagg), (base, implicit, codes)| {
                agg.entry(base)
                    .or_insert_with(Vec::new)
                    .extend(codes.iter().copied());
                if implicit {
                    bagg.insert(base);
                }
                (agg, bagg)
            },
        )
}

#[inline]
fn duplex_aware_strand(tag_strand: Strand, reverse: bool) -> MkResult<Strand> {
    match (tag_strand, reverse) {
        (Strand::Positive, true) | (Strand::Negative, false) => {
            Ok(Strand::Negative)
        }
        (Strand::Positive, false) | (Strand::Negative, true) => {
            Ok(Strand::Positive)
        }
    }
}

#[inline]
fn update_mods_iter2<'a>(
    modbase_iter: &mut BaseModificationsPositionIter<'a>,
    mod_pos: &mut Option<usize>,
    canonical_base: &mut Option<DnaBase>,
    pos_base_mod_call: &mut Option<BaseModCall>,
    mod_strand: &mut Option<Strand>,
    reverse: bool,
    caller: &MultipleThresholdModCaller,
) -> MkResult<()> {
    if let Some(res) = modbase_iter.next() {
        match res {
            Ok((pos, codes)) => {
                let pos = pos as usize;
                let Ok(can_base) = DnaBase::try_from(codes[0].canonical_base)
                else {
                    return Err(MkError::InvalidMm(
                        "invalid DNA base".to_string(),
                    ));
                };
                let tag_mod_strand = Strand::try_from(codes[0].strand)?;
                let can_base = match tag_mod_strand {
                    Strand::Positive => can_base,
                    Strand::Negative => can_base.complement(),
                };
                let base_mod_probs = codes
                    .iter()
                    .map(|hts_code| {
                        (
                            ModCodeRepr::from(hts_code.modified_base),
                            qual_to_prob(hts_code.qual),
                        )
                    })
                    .collect();
                let base_mod_probs = BaseModProbs::new(base_mod_probs, false);

                *pos_base_mod_call =
                    Some(caller.call(&can_base, &base_mod_probs));
                *mod_pos = Some(pos);
                *canonical_base = Some(can_base);
                let strand = duplex_aware_strand(tag_mod_strand, reverse)?;
                *mod_strand = Some(strand)
            }
            Err(e) => {
                return Err(MkError::HtsLibError(e));
            }
        }
    } else {
        *mod_pos = None;
        *canonical_base = None;
        *pos_base_mod_call = None;
        *mod_strand = None;
    }
    Ok(())
}

fn indices_to_byte(motif_idxs: &BitSlice) -> u8 {
    let mut agg = 0u8;
    for (i, b) in motif_idxs.iter().enumerate() {
        let b = *b as u8;
        agg |= b << i;
    }
    agg
}

#[cfg(test)]
mod tests {
    use rust_htslib::bam::{
        self,
        header::HeaderRecord,
        record::{Aux, Cigar, CigarString},
    };
    use tempfile::tempdir;

    use crate::{
        mod_base_code::{
            DnaBase, ModCodeRepr, ANY_CYTOSINE, FORMYL_CYTOSINE,
            HYDROXY_METHYL_CYTOSINE, METHYL_CYTOSINE, TWO_OME_CYTOSINE,
        },
        pileup::base_mods_adapter::ModState,
    };

    use super::{
        ACountsMatrix, ChromCoordinates, CountsMatrix, DnaCytosineCombine,
        DnaModOption, Dynamic, FocusPositions2, GenericPileupWorker,
        ModBasePileup2, PileupNumericOptions, PileupWorker, SIX_METHYL_ADENINE,
    };

    fn make_aligned_mod_record(
        name: &[u8],
        mm: &str,
        ml: &[u8],
    ) -> bam::Record {
        let mut record = bam::Record::new();
        record.set(
            name,
            Some(&CigarString(vec![Cigar::Match(5)])),
            b"AAAAA",
            &[255; 5],
        );
        record.set_tid(0);
        record.set_pos(0);
        record.set_mapq(60);
        record.push_aux(b"MM", Aux::String(mm)).unwrap();
        record.push_aux(b"ML", Aux::ArrayU8(ml.into())).unwrap();
        record
    }

    #[test]
    fn generic_pileup_clears_mod_state_between_records() {
        let temp_dir = tempdir().unwrap();
        let bam_path = temp_dir.path().join("two-records.bam");
        let mut header = bam::Header::new();
        let mut hd = HeaderRecord::new(b"HD");
        hd.push_tag(b"VN", "1.6").push_tag(b"SO", "coordinate");
        header.push_record(&hd);
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1").push_tag(b"LN", 5);
        header.push_record(&sq);

        {
            let mut writer =
                bam::Writer::from_path(&bam_path, &header, bam::Format::Bam)
                    .unwrap();
            writer
                .write(&make_aligned_mod_record(b"modified", "A+a.,4;", &[255]))
                .unwrap();
            writer
                .write(&make_aligned_mod_record(b"canonical", "A+a.;", &[]))
                .unwrap();
        }
        bam::index::build(&bam_path, None, bam::index::Type::Bai, 1).unwrap();

        let mut worker = GenericPileupWorker::new(
            &bam_path,
            None,
            Vec::new(),
            0.0,
            [0.0; 4],
            &[],
            PileupNumericOptions::Passthrough,
            false,
            u16::MAX,
        )
        .unwrap();
        let coordinates =
            ChromCoordinates::new(0, 0, 5, FocusPositions2::AllPositions, true);
        let pileup =
            worker.process(coordinates, ModBasePileup2::new_empty()).unwrap();
        let counts = pileup
            .position_feature_counts
            .iter()
            .find(|counts| {
                counts.position == 4 && counts.mod_code == SIX_METHYL_ADENINE
            })
            .unwrap();

        assert_eq!(counts.filtered_coverage, 2);
        assert_eq!(counts.n_modified, 1);
        assert_eq!(counts.n_canonical, 1);
        assert_eq!(counts.n_other_modified, 0);
        assert_eq!(counts.n_delete, 0);
        assert_eq!(counts.n_filtered, 0);
        assert_eq!(counts.n_diff, 0);
        assert_eq!(counts.n_nocall, 0);
        assert_eq!(
            f32::from(counts.n_modified) / f32::from(counts.filtered_coverage),
            0.5
        );
    }

    type DynamicCounts =
        (u32, char, u16, ModCodeRepr, u16, u16, u16, u16, u16, u16, u16, u8);

    fn increment_dynamic_call(
        matrix: &mut CountsMatrix,
        rpos: u32,
        primary_base: DnaBase,
        reference_base: DnaBase,
        mod_code: ModCodeRepr,
        modified: bool,
        filtered: bool,
        reverse: bool,
        haplotype: u8,
    ) {
        <CountsMatrix as ACountsMatrix<Dynamic, 2, true>>::incr_call(
            matrix,
            ModState {
                mod_position: rpos as usize,
                modified,
                filtered,
                mod_code,
                primary_base,
                inferred: false,
                mod_qual: 255,
            },
            rpos,
            reference_base,
            reverse,
            haplotype,
            false,
        );
    }

    fn valid_dynamic_counts(
        counts: &[crate::pileup::PileupFeatureCounts2],
    ) -> Vec<DynamicCounts> {
        counts
            .iter()
            .filter(|counts| counts.is_valid())
            .map(|counts| {
                (
                    counts.position,
                    counts.raw_strand,
                    counts.filtered_coverage,
                    counts.mod_code,
                    counts.n_canonical,
                    counts.n_modified,
                    counts.n_other_modified,
                    counts.n_delete,
                    counts.n_filtered,
                    counts.n_diff,
                    counts.n_nocall,
                    counts.motif_idxs,
                )
            })
            .collect()
    }

    fn new_cytosine_combine_matrix() -> CountsMatrix {
        <CountsMatrix as ACountsMatrix<DnaCytosineCombine, 2, false>>::new(
            1,
            false,
            Vec::new(),
            0,
            u16::MAX,
        )
    }

    fn increment_cytosine_call(
        matrix: &mut CountsMatrix,
        modified: bool,
        filtered: bool,
        mod_code: ModCodeRepr,
    ) {
        <CountsMatrix as ACountsMatrix<
            DnaCytosineCombine,
            2,
            false,
        >>::incr_call(
            matrix,
            ModState {
                mod_position: 0,
                modified,
                filtered,
                mod_code,
                primary_base: DnaBase::C,
                inferred: false,
                mod_qual: 255,
            },
            0,
            DnaBase::C,
            false,
            0,
            true,
        );
    }

    fn decode_cytosine_combine(
        matrix: &CountsMatrix,
    ) -> Vec<crate::pileup::PileupFeatureCounts2> {
        <CountsMatrix as ACountsMatrix<DnaCytosineCombine, 2, false>>::decode(
            matrix,
            100,
            DnaModOption::Combine,
            false,
        )
    }

    #[test]
    fn dynamic_explicit_mod_offsets_are_keyed_by_base_and_code() {
        let mod_codes =
            [(DnaBase::A, METHYL_CYTOSINE), (DnaBase::C, METHYL_CYTOSINE)];

        assert_eq!(
            super::dynamic_mod_offset(
                &mod_codes,
                DnaBase::A,
                METHYL_CYTOSINE,
                false,
            ),
            super::DYN_N_CONSTANT_COUNTS
        );
        assert_eq!(
            super::dynamic_mod_offset(
                &mod_codes,
                DnaBase::C,
                METHYL_CYTOSINE,
                false,
            ),
            super::DYN_N_CONSTANT_COUNTS + 1
        );
        assert_eq!(
            super::dynamic_mod_offset(
                &mod_codes,
                DnaBase::T,
                METHYL_CYTOSINE,
                false,
            ),
            super::DYN_OTHER_MOD_T
        );
    }

    #[test]
    fn dynamic_explicit_slots_preserve_statuses_across_strand_and_phase() {
        let mod_codes =
            vec![(DnaBase::A, METHYL_CYTOSINE), (DnaBase::C, METHYL_CYTOSINE)];
        let mut matrix = <CountsMatrix as ACountsMatrix<Dynamic, 2, true>>::new(
            2,
            true,
            mod_codes,
            0,
            u16::MAX,
        );

        // The unphased partition proves that identical modification codes on
        // different primary bases retain independent requested slots.
        increment_dynamic_call(
            &mut matrix,
            0,
            DnaBase::A,
            DnaBase::A,
            METHYL_CYTOSINE,
            true,
            false,
            false,
            0,
        );
        increment_dynamic_call(
            &mut matrix,
            1,
            DnaBase::C,
            DnaBase::C,
            METHYL_CYTOSINE,
            true,
            false,
            false,
            0,
        );

        // Exercise every dynamic status in the reverse HP1 partition. A zero
        // motif offset intentionally isolates slot mapping from motif shifts.
        for (mod_code, modified, filtered) in [
            (METHYL_CYTOSINE, true, false),
            (METHYL_CYTOSINE, false, false),
            (HYDROXY_METHYL_CYTOSINE, true, false),
            (METHYL_CYTOSINE, false, true),
        ] {
            increment_dynamic_call(
                &mut matrix,
                1,
                DnaBase::C,
                DnaBase::C,
                mod_code,
                modified,
                filtered,
                true,
                1,
            );
        }
        <CountsMatrix as ACountsMatrix<Dynamic, 2, true>>::incr_diff_call(
            &mut matrix,
            1,
            DnaBase::A,
            DnaBase::C,
            true,
            1,
        );
        <CountsMatrix as ACountsMatrix<Dynamic, 2, true>>::incr_diff_call(
            &mut matrix,
            1,
            DnaBase::C,
            DnaBase::C,
            true,
            1,
        );
        <CountsMatrix as ACountsMatrix<Dynamic, 2, true>>::incr_delete(
            &mut matrix,
            1,
            true,
            1,
        );

        // Repeat the same status oracle for adenine in the forward HP2
        // partition, including an unrequested modification code.
        for (mod_code, modified, filtered) in [
            (METHYL_CYTOSINE, true, false),
            (METHYL_CYTOSINE, false, false),
            (ModCodeRepr::Code('x'), true, false),
            (METHYL_CYTOSINE, false, true),
        ] {
            increment_dynamic_call(
                &mut matrix,
                0,
                DnaBase::A,
                DnaBase::A,
                mod_code,
                modified,
                filtered,
                false,
                2,
            );
        }
        <CountsMatrix as ACountsMatrix<Dynamic, 2, true>>::incr_diff_call(
            &mut matrix,
            0,
            DnaBase::C,
            DnaBase::A,
            false,
            2,
        );
        <CountsMatrix as ACountsMatrix<Dynamic, 2, true>>::incr_diff_call(
            &mut matrix,
            0,
            DnaBase::A,
            DnaBase::A,
            false,
            2,
        );
        <CountsMatrix as ACountsMatrix<Dynamic, 2, true>>::incr_delete(
            &mut matrix,
            0,
            false,
            2,
        );

        let decoded = <CountsMatrix as ACountsMatrix<Dynamic, 2, true>>::decode(
            &matrix,
            100,
            DnaModOption::Other(ANY_CYTOSINE),
            true,
        );
        let partition_width = 2 * 2 * 2;
        assert_eq!(decoded.len(), partition_width * 3);
        assert_eq!(
            valid_dynamic_counts(&decoded[..partition_width]),
            vec![
                (100, '+', 1, METHYL_CYTOSINE, 0, 1, 0, 0, 0, 0, 0, 0),
                (101, '+', 1, METHYL_CYTOSINE, 0, 1, 0, 0, 0, 0, 0, 0),
            ]
        );
        assert_eq!(
            valid_dynamic_counts(
                &decoded[partition_width..partition_width * 2]
            ),
            vec![(101, '-', 3, METHYL_CYTOSINE, 1, 1, 1, 1, 1, 1, 1, 0,)]
        );
        assert_eq!(
            valid_dynamic_counts(&decoded[partition_width * 2..]),
            vec![(100, '+', 3, METHYL_CYTOSINE, 1, 1, 1, 1, 1, 1, 1, 0,)]
        );
    }

    #[test]
    fn dna_cytosine_combine_counts_every_c_modification() {
        let mut matrix = new_cytosine_combine_matrix();
        increment_cytosine_call(&mut matrix, false, false, ANY_CYTOSINE);
        for mod_code in [
            METHYL_CYTOSINE,
            HYDROXY_METHYL_CYTOSINE,
            FORMYL_CYTOSINE,
            TWO_OME_CYTOSINE,
        ] {
            increment_cytosine_call(&mut matrix, true, false, mod_code);
        }

        let decoded = decode_cytosine_combine(&matrix);
        let counts = &decoded[0];
        assert_eq!(
            (
                counts.position,
                counts.raw_strand,
                counts.filtered_coverage,
                counts.mod_code,
                counts.n_canonical,
                counts.n_modified,
                counts.n_other_modified,
                counts.n_delete,
                counts.n_filtered,
                counts.n_diff,
                counts.n_nocall,
                counts.motif_idxs,
            ),
            (100, '+', 5, ANY_CYTOSINE, 1, 4, 0, 0, 0, 0, 0, 0)
        );
        assert_eq!(
            counts.filtered_coverage,
            counts
                .n_canonical
                .saturating_add(counts.n_modified)
                .saturating_add(counts.n_other_modified)
        );
    }

    #[test]
    fn dna_cytosine_combine_counts_filtered_c_calls_saturating() {
        let mut matrix = new_cytosine_combine_matrix();
        increment_cytosine_call(&mut matrix, false, false, ANY_CYTOSINE);
        increment_cytosine_call(&mut matrix, false, true, ANY_CYTOSINE);
        let decoded = decode_cytosine_combine(&matrix);
        let counts = &decoded[0];
        assert_eq!(
            (
                counts.position,
                counts.raw_strand,
                counts.filtered_coverage,
                counts.mod_code,
                counts.n_canonical,
                counts.n_modified,
                counts.n_other_modified,
                counts.n_delete,
                counts.n_filtered,
                counts.n_diff,
                counts.n_nocall,
                counts.motif_idxs,
            ),
            (100, '+', 1, ANY_CYTOSINE, 1, 0, 0, 0, 1, 0, 0, 0)
        );
        assert_eq!(
            counts.filtered_coverage,
            counts
                .n_canonical
                .saturating_add(counts.n_modified)
                .saturating_add(counts.n_other_modified)
        );

        matrix.inner[super::FILT_C_OFFSET] = u16::MAX;
        increment_cytosine_call(&mut matrix, false, true, ANY_CYTOSINE);
        assert_eq!(decode_cytosine_combine(&matrix)[0].n_filtered, u16::MAX);
    }
}
