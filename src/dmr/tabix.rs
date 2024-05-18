use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::fs::File;
use std::io::BufRead;
use std::ops::Range;
use std::path::{Path, PathBuf};

use anyhow::{anyhow, bail, Context};
use log::debug;
use noodles::csi::{
    binning_index::index::reference_sequence::bin::Chunk as IndexChunk,
    BinningIndex,
};
use noodles::tabix::Index as CsiIndex;
use rayon::prelude::*;
use rustc_hash::FxHashMap;

use crate::dmr::bedmethyl::{aggregate_counts2, BedMethylLine};
use crate::dmr::llr_model::AggregatedCounts;
use crate::dmr::util::{
    n_choose_2, DmrBatch, DmrBatchOfPositions, ProtoIndexChunk,
};
use crate::genome_positions::StrandedPosition;
use crate::mod_base_code::DnaBase;
use crate::monoid::Moniod;

pub(super) type ChromToBedmethylLines = FxHashMap<String, Vec<BedMethylLine>>;
pub(super) type SampleToBedMethyLines =
    FxHashMap<usize, FxHashMap<String, Vec<BedMethylLine>>>;
/// Chrom -> {StrandedPosition -> <X_i>} for all i in samples
pub(super) type ChromToAggregatedCountsByPosition = FxHashMap<
    String,
    BTreeMap<StrandedPosition<DnaBase>, Vec<AggregatedCounts>>,
>;
pub(super) type BedMethylLinesResult<T> = anyhow::Result<(T, T)>;

pub(super) struct IndexHandler {
    csi_index: CsiIndex,
    bedmethyl_fp: PathBuf,
    contig_to_id: FxHashMap<String, usize>,
}

impl IndexHandler {
    pub(super) fn new(
        bedmethyl_fp: &PathBuf,
        index_fp: &PathBuf,
    ) -> anyhow::Result<Self> {
        let csi_index = noodles::tabix::read(index_fp).with_context(|| {
            format!("failed to read tabix index {index_fp:?}")
        })?;
        if let Some(header) = csi_index.header() {
            if bedmethyl_fp.exists() {
                let contig_to_id = header
                    .reference_sequence_names()
                    .iter()
                    .enumerate()
                    .map(|(id, name)| (name.to_owned(), id))
                    .collect::<FxHashMap<String, usize>>();
                Ok(Self {
                    csi_index,
                    bedmethyl_fp: bedmethyl_fp.to_owned(),
                    contig_to_id,
                })
            } else {
                bail!("bedMethyl file not found at {bedmethyl_fp:?}")
            }
        } else {
            bail!("could not read tabix header from {index_fp:?}")
        }
    }

    pub(super) fn new_infer_index_filepath(
        bedmethyl_fp: &PathBuf,
        index_fp: Option<&PathBuf>,
    ) -> anyhow::Result<Self> {
        if let Some(index_fp) = index_fp {
            debug!("loading specified index at {}", index_fp.to_string_lossy());
            IndexHandler::new(bedmethyl_fp, index_fp)
        } else {
            let bedmethyl_root = bedmethyl_fp.to_str().ok_or_else(|| {
                anyhow!("could not format bedmethyl filepath, {bedmethyl_fp:?}")
            })?;
            let index_fp = format!("{}.tbi", bedmethyl_root);
            debug!(
                "looking for index associated with {} at {}",
                bedmethyl_root, &index_fp
            );
            let index_path = Path::new(&index_fp).to_path_buf();
            IndexHandler::new(bedmethyl_fp, &index_path).with_context(|| {
                format!(
                    "failed to read index inferred from bedMethyl file name \
                     at {index_fp}"
                )
            })
        }
    }

    /// Reads bedmethyl records from the file. Filters:
    /// 1. N_valid_cov >= min_valid_cov
    /// 2. mod_code must be in [`crate::mod_base_code::SUPPORTED_CODES`]
    fn read_bedmethyl(
        &self,
        chunks: &[IndexChunk],
        min_valid_coverage: u64,
    ) -> anyhow::Result<Vec<BedMethylLine>> {
        let mut reader =
            File::open(&self.bedmethyl_fp).map(noodles::bgzf::Reader::new)?;
        let mut bedmethyl_lines = Vec::new();
        let mut failed_to_parse = 0;
        // let mut successfully_parsed = 0usize;
        // todo could be a fold instead
        for chunk in chunks {
            reader.seek(chunk.start())?;
            let mut lines = Vec::new();
            // todo come back and make this one loop
            'readloop: loop {
                let mut buf = String::new();
                let _byts = reader.read_line(&mut buf)?;
                lines.push(buf);
                let cur_pos = reader.virtual_position();
                if cur_pos >= chunk.end() {
                    break 'readloop;
                }
            }
            for line in lines.iter() {
                match BedMethylLine::parse(line) {
                    Ok(bm_line) => {
                        bedmethyl_lines.push(bm_line);
                        // successfully_parsed += 1;
                    }
                    Err(_e) => {
                        // trace!("failed to parse line {line}");
                        failed_to_parse += 1
                    }
                }
            }
        }

        if failed_to_parse > 0 {
            debug!(
                "failed to parse {failed_to_parse} lines from {:?}",
                &self.bedmethyl_fp
            );
        }

        Ok(bedmethyl_lines
            .into_iter()
            .filter(|bml| bml.valid_coverage >= min_valid_coverage)
            .filter(|bml| {
                if bml.check_mod_code_supported() {
                    true
                } else {
                    debug!(
                        "encountered illegal mod-code for DMR, {}",
                        bml.raw_mod_code
                    );
                    false
                }
            })
            .collect())
    }

    fn has_contig(&self, contig_name: &str) -> bool {
        self.contig_to_id.contains_key(contig_name)
    }
}

pub(super) struct MultiSampleIndex {
    index_handlers: Vec<IndexHandler>,
    min_valid_coverage: u64,
}

impl MultiSampleIndex {
    pub(super) fn new(
        handlers: Vec<IndexHandler>,
        min_valid_coverage: u64,
    ) -> Self {
        Self { index_handlers: handlers, min_valid_coverage }
    }

    #[inline]
    fn read_bedmethyl_lines_from_batch(
        &self,
        // todo needs documentation of what this data structure is
        chunks: &FxHashMap<usize, BTreeSet<ProtoIndexChunk>>,
    ) -> anyhow::Result<FxHashMap<usize, FxHashMap<String, Vec<BedMethylLine>>>>
    {
        // take all the mappings of sample_id to chunks
        let groups = chunks
            .par_iter() // yah
            .filter_map(|(id, chunks)| {
                // get the index handler for each
                // shouldn't ever really get an miss here, but
                // just in case do a filter_map
                self.index_handlers.get(*id).map(|handler| {
                    // optimize/join overlapping chunks
                    let index_chunks = join_chunks_together(chunks);
                    // each pair is for one tabix index
                    (*id, handler, index_chunks)
                })
            })
            .map(|(sample_id, handler, chunks)| {
                // actually read the bedmethyl here
                let grouped_by_chrom = handler
                    .read_bedmethyl(&chunks, self.min_valid_coverage)
                    .map(|bedmethyl_lines| {
                        // group the bedmethyl lines by contig
                        // in case this batch has more than one contigs'
                        // worth in it
                        bedmethyl_lines.into_iter().fold(
                            FxHashMap::default(),
                            |mut agg, bml| {
                                // groupby chrom
                                agg.entry(bml.chrom.to_owned())
                                    .or_insert(Vec::new())
                                    .push(bml);
                                agg
                            },
                        )
                    });
                grouped_by_chrom.map(|grouped| (sample_id, grouped))
            })
            .collect::<anyhow::Result<
                FxHashMap<usize, FxHashMap<String, Vec<BedMethylLine>>>,
            >>()?;

        Ok(groups)
    }

    // remember this method _only_ reads the lines, it does not perform any
    // filtering
    fn read_bedmethyl_lines<T>(
        &self,
        dmr_batch: &DmrBatch<T>,
    ) -> BedMethylLinesResult<SampleToBedMethyLines> {
        // std::thread::scope(|sc| )
        let bedmethyl_lines_a =
            self.read_bedmethyl_lines_from_batch(&dmr_batch.chunks_a)?;
        let bedmethyl_lines_b =
            self.read_bedmethyl_lines_from_batch(&dmr_batch.chunks_b)?;

        Ok((bedmethyl_lines_a, bedmethyl_lines_b))
    }

    pub(super) fn read_bedmethyl_lines_collapse_on_chrom<T>(
        &self,
        dmr_batch: &DmrBatch<T>,
    ) -> BedMethylLinesResult<ChromToBedmethylLines> {
        let (bedmethyl_lines_a, bedmethyl_lines_b) =
            self.read_bedmethyl_lines(dmr_batch)?;
        let collapse_f =
            |sample_lines: SampleToBedMethyLines| -> ChromToBedmethylLines {
                sample_lines.into_iter().fold(
                    FxHashMap::default(),
                    |mut agg, (_sample_id, chrom_to_lines)| {
                        for (chrom, mut lines) in chrom_to_lines {
                            agg.entry(chrom)
                                .or_insert(Vec::new())
                                .append(&mut lines);
                        }
                        agg
                    },
                )
            };
        let a = collapse_f(bedmethyl_lines_a);
        let b = collapse_f(bedmethyl_lines_b);
        Ok((a, b))
    }

    pub(super) fn samples(&self) -> Vec<usize> {
        (0..self.index_handlers.len()).collect()
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

    pub(super) fn get_index_chunks(
        &self,
        sample_index: usize,
        chrom: &str,
        interval: &Range<u64>,
    ) -> anyhow::Result<Vec<IndexChunk>> {
        if let Some((tabix_index, chrom_id)) =
            self.index_handlers.get(sample_index).and_then(|index| {
                index.contig_to_id.get(chrom).map(|chrom_id| (index, chrom_id))
            })
        {
            let start =
                noodles::core::Position::new((interval.start + 1) as usize)
                    .unwrap();
            let end = noodles::core::Position::new((interval.end + 1) as usize)
                .unwrap();
            let interval = noodles::core::region::Interval::from(start..=end);
            tabix_index.csi_index.query(*chrom_id, interval).map_err(|e| {
                anyhow!(
                    "failed to query sample index {sample_index} for interval \
                     {interval:?}, {e}"
                )
            })
        } else {
            Ok(Vec::new())
        }
    }

    // todo try and make this return &String
    pub(super) fn all_contigs(&self) -> HashSet<String> {
        self.index_handlers
            .iter()
            .flat_map(|handler| {
                handler
                    .contig_to_id
                    .keys()
                    .map(|s| s.to_owned())
                    .collect::<HashSet<String>>()
            })
            .collect()
    }
}

/// Handles getting chunks for multiple indices and reading bedMethyl data

pub(super) struct SingleSiteSampleIndex {
    multi_sample_index: MultiSampleIndex,
    control_idxs: Vec<usize>,
    exp_idxs: Vec<usize>,
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

    #[inline]
    fn query_chunks(
        &self,
        contig_name: &str,
        interval: &Range<u64>,
        is_control: bool,
    ) -> anyhow::Result<FxHashMap<usize, Vec<IndexChunk>>> {
        let idxs = if is_control { &self.control_idxs } else { &self.exp_idxs };
        idxs.iter()
            .map(|idx| {
                let chunks = self.multi_sample_index.get_index_chunks(
                    *idx,
                    contig_name,
                    interval,
                );
                chunks.map(|x| (*idx, x))
            })
            .collect::<anyhow::Result<FxHashMap<usize, Vec<IndexChunk>>>>()
    }

    pub(super) fn query_control_chunks(
        &self,
        contig_name: &str,
        interval: &Range<u64>,
    ) -> anyhow::Result<FxHashMap<usize, Vec<IndexChunk>>> {
        self.query_chunks(contig_name, interval, true)
    }

    pub(super) fn query_exp_chunks(
        &self,
        contig_name: &str,
        interval: &Range<u64>,
    ) -> anyhow::Result<FxHashMap<usize, Vec<IndexChunk>>> {
        self.query_chunks(contig_name, interval, false)
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
        sample: SampleToBedMethyLines,
    ) -> anyhow::Result<ChromToAggregatedCountsByPosition> {
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
                                agg.entry(l.get_stranded_position())
                                    .or_insert(Vec::new())
                                    .push(l);
                                agg
                            },
                        )
                        .reduce(|| BTreeMap::new(), |a, b| a.op(b))
                        .into_par_iter()
                        .map(|(position, lines)| {
                            (position, aggregate_counts2(&lines))
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
        bedmethyl_lines: SampleToBedMethyLines,
    ) -> SampleToBedMethyLines {
        bedmethyl_lines
            .into_iter()
            .map(|(sample, bm_lines)| {
                let filtered_chrom_to_lines = bm_lines
                    .into_iter()
                    .filter_map(|(chrom, lines)| {
                        let filtered_lines = lines
                            .into_par_iter()
                            .filter(|l| {
                                dmr_batch.contains_position(
                                    &chrom,
                                    &l.get_stranded_position(),
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
                (sample, filtered_chrom_to_lines)
            })
            .collect::<FxHashMap<usize, FxHashMap<String, Vec<BedMethylLine>>>>(
            )
    }

    pub(super) fn read_bedmethyl_lines_filtered_by_position(
        &self,
        dmr_batch: &DmrBatchOfPositions,
    ) -> BedMethylLinesResult<SampleToBedMethyLines> {
        let (bedmethyl_lines_a, bedmethyl_lines_b) =
            self.multi_sample_index.read_bedmethyl_lines(&dmr_batch)?;

        let filt_lines_a = Self::intersect_bedmethyl_lines_with_sites(
            &dmr_batch,
            bedmethyl_lines_a,
        );
        let filt_lines_b = Self::intersect_bedmethyl_lines_with_sites(
            &dmr_batch,
            bedmethyl_lines_b,
        );

        Ok((filt_lines_a, filt_lines_b))
    }

    pub(super) fn read_bedmethyl_lines_organized_by_position(
        &self,
        dmr_batch: DmrBatchOfPositions,
    ) -> BedMethylLinesResult<ChromToAggregatedCountsByPosition> {
        let (bedmethyl_lines_a, bedmethyl_lines_b) =
            self.read_bedmethyl_lines_filtered_by_position(&dmr_batch)?;

        let counts_a = Self::organize_bedmethy_lines(bedmethyl_lines_a)?;
        let counts_b = Self::organize_bedmethy_lines(bedmethyl_lines_b)?;
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

#[inline]
fn join_chunks_together(chunks: &BTreeSet<ProtoIndexChunk>) -> Vec<IndexChunk> {
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
