use std::collections::{BTreeMap, HashSet};
use std::ops::Range;

use anyhow::bail;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::dmr::bedmethyl::{aggregate_counts2, BedMethylLine};
use crate::dmr::llr_model::AggregatedCounts;
use crate::dmr::util::{n_choose_2, DmrBatch, DmrBatchOfPositions};
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
/// Chrom -> {StrandedPosition -> <X_i>} for all i in samples
pub(super) type ChromToPosAggregatedCounts = FxHashMap<
    String,
    BTreeMap<StrandedPosition<DnaBase>, Vec<AggregatedCounts>>,
>;
/// Usually (control, experiment)
pub(super) type BedMethylLinesResult<T> = anyhow::Result<(T, T)>;

pub(super) struct MultiSampleIndex {
    index_handlers: Vec<BedMethylTbxIndex>,
    pub code_lookup: FxHashMap<ModCodeRepr, DnaBase>,
    min_valid_coverage: u64,
}

impl MultiSampleIndex {
    pub(super) fn new(
        handlers: Vec<BedMethylTbxIndex>,
        code_lookup: FxHashMap<ModCodeRepr, DnaBase>,
        min_valid_coverage: u64,
    ) -> Self {
        Self { index_handlers: handlers, min_valid_coverage, code_lookup }
    }

    #[inline]
    fn read_bedmethyl_files(
        &self,
        idxs: &FxHashSet<usize>,
        chunks: &FxHashMap<String, Range<u64>>,
    ) -> anyhow::Result<FxHashMap<usize, FxHashMap<String, Vec<BedMethylLine>>>>
    {
        // take all the mappings of sample_id to chunks
        let groups = idxs
            .par_iter() // yah
            .filter_map(|id| {
                // get the index handler for each
                // shouldn't ever really get a miss here, but
                // just in case do a filter_map
                self.index_handlers
                    .get(*id)
                    .map(|handler| (*id, handler, chunks))
            })
            // chunks is a mapping of each chrom to the range in that chrom to
            // fetch
            .map(|(sample_id, handler, chunks)| {
                // actually read the bedmethyl here
                let grouped_by_chrom =
                        chunks
                            .par_iter()
                            // here we read the bedmethyl and have a mapping of chrom to records
                            .map(|(chrom, range)| {
                                let bm_lines = handler.read_bedmethyl(
                                    chrom,
                                    range,
                                    self.min_valid_coverage,
                                    &self.code_lookup,
                                );
                                bm_lines.map(|lines| (chrom.to_owned(), lines))
                            })
                            .collect::<anyhow::Result<
                                FxHashMap<String, Vec<BedMethylLine>>,
                            >>();
                grouped_by_chrom.map(|grouped| (sample_id, grouped))
            })
            .collect::<anyhow::Result<
                FxHashMap<usize, FxHashMap<String, Vec<BedMethylLine>>>,
            >>()?;

        Ok(groups)
    }

    // remember this method _only_ reads the lines, it does not perform any
    // filtering
    fn read_bedmethyl_lines<T: Default>(
        &self,
        dmr_batch: &DmrBatch<T>,
    ) -> BedMethylLinesResult<SampleToChromBMLines> {
        let bedmethyl_lines_a =
            self.read_bedmethyl_files(&dmr_batch.idxs_a, &dmr_batch.regions)?;
        let bedmethyl_lines_b =
            self.read_bedmethyl_files(&dmr_batch.idxs_b, &dmr_batch.regions)?;

        Ok((bedmethyl_lines_a, bedmethyl_lines_b))
    }

    pub(super) fn read_bedmethyl_group_by_chrom<T: Default>(
        &self,
        dmr_batch: &DmrBatch<T>,
    ) -> BedMethylLinesResult<ChromToSampleBMLines> {
        let (bedmethyl_lines_a, bedmethyl_lines_b) =
            self.read_bedmethyl_lines(dmr_batch)?;
        // todo I think this could be replaced my moniod
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
}

impl SingleSiteSampleIndex {
    pub(super) fn new(
        multi_sample_index: MultiSampleIndex,
        num_a: usize,
        num_b: usize,
    ) -> anyhow::Result<Self> {
        if num_a == 0 || num_b == 0 {
            bail!("must be at least 1 sample for 'a' and 'b'")
        }

        let control_idxs = (0..num_a).collect::<Vec<usize>>();
        let exp_idxs = (0..num_b).map(|x| x + num_a).collect::<Vec<usize>>();

        Ok(Self { multi_sample_index, control_idxs, exp_idxs })
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
        sample: SampleToChromBMLines,
        code_lookup: &FxHashMap<ModCodeRepr, DnaBase>,
    ) -> anyhow::Result<ChromToPosAggregatedCounts> {
        let mut agg = FxHashMap::default();

        // samples should be length ~1-5
        for (_sample_id, filtered_chrom_to_lines) in sample {
            let chrom_to_counts = filtered_chrom_to_lines
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
                            (position, aggregate_counts2(&lines, &code_lookup))
                        })
                        .collect::<FxHashMap<
                            StrandedPosition<DnaBase>,
                            anyhow::Result<AggregatedCounts>,
                        >>();
                    (chrom, grouped_by_position)
                })
                .collect::<Vec<(
                    String,
                    FxHashMap<
                        StrandedPosition<DnaBase>,
                        anyhow::Result<AggregatedCounts>,
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
                            .push(aggregated_counts),
                        Err(e) => {
                            bail!("failed to aggregate counts, {e}")
                        }
                    }
                }
            }
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
        let (bedmethyl_lines_a, bedmethyl_lines_b) =
            self.read_bedmethyl_lines_filtered_by_position(&dmr_batch)?;

        let counts_a = Self::organize_bedmethy_lines(
            bedmethyl_lines_a,
            &self.multi_sample_index.code_lookup,
        )?;
        let counts_b = Self::organize_bedmethy_lines(
            bedmethyl_lines_b,
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

    #[inline]
    pub(super) fn multiple_samples(&self) -> bool {
        self.num_a_samples() > 1 || self.num_b_samples() > 1
    }
}
