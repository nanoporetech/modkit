use std::collections::{BTreeMap, HashSet};
use std::fmt::Debug;
use std::ops::Range;

use anyhow::bail;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::dmr::bedmethyl::{aggregate_counts2, BedMethylLine};
use crate::dmr::llr_model::AggregatedCounts;
use crate::dmr::util::{n_choose_2, DmrBatch, DmrBatchOfPositions};
use crate::errs::{MkError, MkResult};
use crate::genome_positions::StrandedPosition;
use crate::mod_base_code::{DnaBase, ModCodeRepr};
use crate::monoid::Moniod;
use crate::tabix::BedMethylTbxIndex;

/// Chrom -> {Sample id -> <bedmethyl_records>}
pub(super) type ChromToSampleBMLines =
    FxHashMap<String, FxHashMap<usize, Vec<BedMethylLine>>>;
/// Sample id -> {Chrom -> <BedMethyl_records>}
pub(super) type SampleToChromBMLines =
    FxHashMap<usize, FxHashMap<String, Vec<BedMethylLine>>>;
#[derive(Debug)]
pub(super) struct SampleCount {
    pub(super) sample_id: usize,
    pub(super) counts: AggregatedCounts,
}
/// Chrom -> {StrandedPosition -> <X_i>} for all i in samples
pub(super) type ChromToPosAggregatedCounts =
    FxHashMap<String, BTreeMap<StrandedPosition<DnaBase>, Vec<SampleCount>>>;
/// Usually (control, experiment)
pub(super) type BedMethylLinesResult<T> = MkResult<(T, T)>;

pub(super) struct MultiSampleIndex {
    index_handlers: Vec<BedMethylTbxIndex>,
    pub code_lookup: FxHashMap<ModCodeRepr, DnaBase>,
    min_valid_coverage: u64,
    io_threads: usize,
}

impl MultiSampleIndex {
    pub(super) fn new(
        handlers: Vec<BedMethylTbxIndex>,
        code_lookup: FxHashMap<ModCodeRepr, DnaBase>,
        min_valid_coverage: u64,
        io_threads: usize,
    ) -> Self {
        Self {
            index_handlers: handlers,
            min_valid_coverage,
            code_lookup,
            io_threads,
        }
    }

    #[inline]
    fn read_bedmethyl_files(
        &self,
        idxs: &FxHashSet<usize>,
        chunks: &FxHashMap<String, Range<u64>>,
    ) -> MkResult<SampleToChromBMLines> {
        // take all the mappings of sample_id to chunks
        let groups =
            idxs.par_iter() // yah
                .filter_map(|id| {
                    // get the index handler for each
                    // shouldn't ever really get a miss here, but
                    // just in case do a filter_map
                    self.index_handlers
                        .get(*id)
                        .map(|handler| (*id, handler, chunks))
                })
                // chunks is a mapping of each chrom to the range in that chrom
                // to fetch
                .map(|(sample_id, handler, chunks)| {
                    // actually read the bedmethyl here
                    let grouped_by_chrom =
                        chunks
                            .par_iter()
                            // here we read the bedmethyl and have a mapping of
                            // chrom to records
                            .map(|(chrom, range)| {
                                let bm_lines = handler
                                    .read_bedmethyl_check_code(
                                        chrom,
                                        range,
                                        self.min_valid_coverage,
                                        &self.code_lookup,
                                        self.io_threads
                                    );
                                bm_lines.map(|lines| (chrom.to_owned(), lines))
                            })
                            .collect::<MkResult<
                                FxHashMap<String, Vec<BedMethylLine>>,
                            >>();
                    grouped_by_chrom.map(|grouped| (sample_id, grouped))
                })
                .collect::<MkResult<
                    FxHashMap<usize, FxHashMap<String, Vec<BedMethylLine>>>,
                >>()?;

        Ok(groups)
    }

    fn read_bedmethyl_lines<T: Default + Debug>(
        &self,
        dmr_batch: &DmrBatch<T>,
    ) -> BedMethylLinesResult<SampleToChromBMLines> {
        let bedmethyl_lines_a =
            self.read_bedmethyl_files(&dmr_batch.idxs_a, &dmr_batch.regions)?;
        let bedmethyl_lines_b =
            self.read_bedmethyl_files(&dmr_batch.idxs_b, &dmr_batch.regions)?;

        Ok((bedmethyl_lines_a, bedmethyl_lines_b))
    }

    pub(super) fn read_bedmethyl_group_by_chrom<T: Default + Debug>(
        &self,
        dmr_batch: &DmrBatch<T>,
    ) -> BedMethylLinesResult<ChromToSampleBMLines> {
        let (bedmethyl_lines_a, bedmethyl_lines_b) =
            self.read_bedmethyl_lines(dmr_batch)?;
        // todo I think this could be replaced by moniod
        let traverse_records =
            |sample_lines: SampleToChromBMLines| -> ChromToSampleBMLines {
                sample_lines.into_iter().fold(
                    FxHashMap::default(),
                    |mut agg, (sample, records)| {
                        for (chrom, lines) in records {
                            agg.entry(chrom)
                                .or_insert_with(FxHashMap::default)
                                .entry(sample)
                                .or_insert_with(Vec::new)
                                .extend(lines);
                        }
                        agg
                    },
                )
            };
        let a = traverse_records(bedmethyl_lines_a);
        let b = traverse_records(bedmethyl_lines_b);
        Ok((a, b))
    }

    #[allow(unused)]
    pub(super) fn num_combinations(&self) -> anyhow::Result<usize> {
        n_choose_2(self.index_handlers.len())
    }

    pub(super) fn has_contig(
        &self,
        sample_id: usize,
        contig_name: &str,
    ) -> bool {
        self.index_handlers
            .get(sample_id)
            .map(|handler| handler.has_contig(&contig_name))
            .unwrap_or(false)
    }

    // todo try and make this return &String
    pub(super) fn all_contigs(&self) -> HashSet<String> {
        self.index_handlers
            .iter()
            .flat_map(|handler| handler.get_contigs())
            .collect()
    }
}

/// Handles getting chunks for multiple indices and reading bedMethyl data

pub(super) struct SingleSiteSampleIndex {
    multi_sample_index: MultiSampleIndex,
    pub(super) control_idxs: Vec<usize>,
    pub(super) exp_idxs: Vec<usize>,
    single_code_op: Option<ModCodeRepr>,
}

impl SingleSiteSampleIndex {
    pub(super) fn new(
        multi_sample_index: MultiSampleIndex,
        num_a: usize,
        num_b: usize,
        single_code_op: Option<ModCodeRepr>,
    ) -> anyhow::Result<Self> {
        if num_a == 0 || num_b == 0 {
            bail!("must be at least 1 sample for 'a' and 'b'")
        }

        let control_idxs = (0..num_a).collect::<Vec<usize>>();
        let exp_idxs = (0..num_b).map(|x| x + num_a).collect::<Vec<usize>>();

        Ok(Self { multi_sample_index, control_idxs, exp_idxs, single_code_op })
    }

    pub(super) fn has_contig(&self, name: &str) -> bool {
        let control_has_contig = self
            .control_idxs
            .iter()
            .any(|idx| self.multi_sample_index.has_contig(*idx, &name));
        let exp_has_contig = || {
            self.exp_idxs
                .iter()
                .any(|idx| self.multi_sample_index.has_contig(*idx, &name))
        };
        control_has_contig && exp_has_contig()
    }

    fn organize_bedmethy_lines(
        &self,
        mut sample: SampleToChromBMLines,
        configured_sample_ids: &[usize],
        code_lookup: &FxHashMap<ModCodeRepr, DnaBase>,
    ) -> MkResult<ChromToPosAggregatedCounts> {
        let mut agg = FxHashMap::default();

        // samples should be length ~1-5
        for sample_id in configured_sample_ids {
            let chrom_to_filtered_bm_records =
                sample.remove(sample_id).ok_or_else(|| {
                    MkError::InvalidBedMethyl(format!(
                        "missing configured sample ID {sample_id} while \
                         organizing bedMethyl records"
                    ))
                })?;
            let chrom_to_counts = chrom_to_filtered_bm_records
                .into_iter()
                .map(|(chrom, lines)| {
                    let grouped_by_position = lines
                        .into_par_iter()
                        .fold(
                            || BTreeMap::new(),
                            |mut agg, l| {
                                agg.entry(
                                    l.get_stranded_position(&code_lookup),
                                )
                                .or_insert(Vec::new())
                                .push(l);
                                agg
                            },
                        )
                        .reduce(|| BTreeMap::new(), |a, b| a.op(b))
                        .into_par_iter()
                        .map(|(position, lines)| {
                            (
                                position,
                                aggregate_counts2(
                                    &lines,
                                    &code_lookup,
                                    self.single_code_op,
                                ),
                            )
                        })
                        .collect::<FxHashMap<
                            StrandedPosition<DnaBase>,
                            MkResult<AggregatedCounts>,
                        >>();
                    (chrom, grouped_by_position)
                })
                .collect::<Vec<(
                    String,
                    FxHashMap<
                        StrandedPosition<DnaBase>,
                        MkResult<AggregatedCounts>,
                    >,
                )>>();
            // There should be only a few chroms
            for (chrom, grouped_by_position) in chrom_to_counts {
                let chrom_agg = agg.entry(chrom).or_insert(BTreeMap::new());
                // there will be max interval-length positions. todo consider
                // making this parallel and dropping the error?
                for (position, res) in grouped_by_position {
                    match res {
                        Ok(aggregated_counts) => chrom_agg
                            .entry(position)
                            .or_insert(Vec::new())
                            .push(SampleCount {
                                sample_id: *sample_id,
                                counts: aggregated_counts,
                            }),
                        Err(e) => return Err(e),
                    }
                }
            }
        }
        if !sample.is_empty() {
            let mut unexpected_ids = sample.keys().copied().collect::<Vec<_>>();
            unexpected_ids.sort_unstable();
            return Err(MkError::InvalidBedMethyl(format!(
                "found unexpected sample IDs while organizing bedMethyl \
                 records: {unexpected_ids:?}"
            )));
        }
        Ok(agg)
    }

    #[inline]
    fn intersect_bedmethyl_lines_with_sites(
        dmr_batch: &DmrBatchOfPositions,
        bedmethyl_lines: SampleToChromBMLines,
        code_lookup: &FxHashMap<ModCodeRepr, DnaBase>,
    ) -> SampleToChromBMLines {
        bedmethyl_lines
            .into_iter()
            .map(|(sample, bm_lines)| {
                let chrom_to_filtered_lines = bm_lines
                    .into_iter()
                    .filter_map(|(chrom, lines)| {
                        let filtered_lines = lines
                            .into_par_iter()
                            .filter(|l| {
                                dmr_batch.contains_position(
                                    &chrom,
                                    &l.get_stranded_position(&code_lookup),
                                )
                            })
                            .collect::<Vec<BedMethylLine>>();
                        let filtered_lines = if filtered_lines.is_empty() {
                            None
                        } else {
                            Some(filtered_lines)
                        };
                        filtered_lines.map(|ls| (chrom, ls))
                    })
                    .collect::<FxHashMap<String, Vec<BedMethylLine>>>();
                (sample, chrom_to_filtered_lines)
            })
            .collect::<FxHashMap<usize, FxHashMap<String, Vec<BedMethylLine>>>>(
            )
    }

    pub(super) fn read_bedmethyl_lines_filtered_by_position(
        &self,
        dmr_batch: &DmrBatchOfPositions,
    ) -> BedMethylLinesResult<SampleToChromBMLines> {
        let (bedmethyl_lines_a, bedmethyl_lines_b) =
            self.multi_sample_index.read_bedmethyl_lines(&dmr_batch)?;

        // filter down to just the sites we are scoring
        let filt_lines_a = Self::intersect_bedmethyl_lines_with_sites(
            &dmr_batch,
            bedmethyl_lines_a,
            &self.multi_sample_index.code_lookup,
        );
        let filt_lines_b = Self::intersect_bedmethyl_lines_with_sites(
            &dmr_batch,
            bedmethyl_lines_b,
            &self.multi_sample_index.code_lookup,
        );

        Ok((filt_lines_a, filt_lines_b))
    }

    pub(super) fn read_bedmethyl_lines_organized_by_position(
        &self,
        dmr_batch: DmrBatchOfPositions,
    ) -> BedMethylLinesResult<ChromToPosAggregatedCounts> {
        // read the records for the two samples
        let (bedmethyl_lines_a, bedmethyl_lines_b) =
            self.read_bedmethyl_lines_filtered_by_position(&dmr_batch)?;

        // group by chrom, this can fail if the records are deemed invalid
        let counts_a = self.organize_bedmethy_lines(
            bedmethyl_lines_a,
            &self.control_idxs,
            &self.multi_sample_index.code_lookup,
        )?;
        let counts_b = self.organize_bedmethy_lines(
            bedmethyl_lines_b,
            &self.exp_idxs,
            &self.multi_sample_index.code_lookup,
        )?;
        Ok((counts_a, counts_b))
    }

    pub(super) fn num_a_samples(&self) -> usize {
        self.control_idxs.len()
    }

    pub(super) fn num_b_samples(&self) -> usize {
        self.exp_idxs.len()
    }

    pub(super) fn min_valid_coverage(&self) -> u64 {
        self.multi_sample_index.min_valid_coverage
    }

    #[inline]
    pub(super) fn matched_replicate_samples(&self) -> bool {
        self.num_a_samples() == self.num_b_samples() && self.num_a_samples() > 1
    }

    pub(super) fn has_complete_positive_matched_counts(
        &self,
        counts_a: &[SampleCount],
        counts_b: &[SampleCount],
    ) -> bool {
        let aligned_and_positive =
            |expected_ids: &[usize], counts: &[SampleCount]| {
                expected_ids.len() == counts.len()
                    && expected_ids.iter().zip(counts).all(
                        |(expected_id, sample_count)| {
                            expected_id == &sample_count.sample_id
                                && sample_count.counts.total > 0
                        },
                    )
            };
        self.matched_replicate_samples()
            && aligned_and_positive(&self.control_idxs, counts_a)
            && aligned_and_positive(&self.exp_idxs, counts_b)
    }

    #[inline]
    pub(super) fn multiple_samples(&self) -> bool {
        self.num_a_samples() > 1 || self.num_b_samples() > 1
    }
}

#[cfg(test)]
mod sample_identity_tests {
    use rustc_hash::FxHashMap;

    use super::{
        MultiSampleIndex, SampleToChromBMLines, SingleSiteSampleIndex,
    };
    use crate::dmr::bedmethyl::BedMethylLine;
    use crate::mod_base_code::{DnaBase, ModCodeRepr};

    fn line(modified: usize) -> BedMethylLine {
        BedMethylLine::parse(&format!(
            "chr1\t0\t1\tm\t10\t+\t0\t1\t255,0,0\t10\t0.00\t{}\t{}\t0\t0\t0\t0\t0",
            modified,
            10 - modified
        ))
        .unwrap()
    }

    #[test]
    fn organizer_preserves_configured_sample_order() {
        let mut code_lookup = FxHashMap::default();
        code_lookup.insert(ModCodeRepr::Code('m'), DnaBase::C);
        let sample_index = SingleSiteSampleIndex::new(
            MultiSampleIndex::new(Vec::new(), code_lookup.clone(), 0, 1),
            3,
            3,
            None,
        )
        .unwrap();
        let mut samples = SampleToChromBMLines::default();
        for (sample_id, modified) in [(5, 9), (3, 1), (4, 5)] {
            samples.insert(
                sample_id,
                FxHashMap::from_iter([(
                    "chr1".to_string(),
                    vec![line(modified)],
                )]),
            );
        }

        let organized = sample_index
            .organize_bedmethy_lines(
                samples,
                &sample_index.exp_idxs,
                &code_lookup,
            )
            .unwrap();
        let counts = organized.get("chr1").unwrap().values().next().unwrap();

        assert_eq!(
            counts
                .iter()
                .map(|sample_count| {
                    (
                        sample_count.sample_id,
                        sample_count.counts.modified_counts(),
                    )
                })
                .collect::<Vec<_>>(),
            vec![(3, 1), (4, 5), (5, 9)]
        );
    }
}
