use crate::position_filter::StrandedPositionFilter;
use crate::util::{ReferenceRecord, Region};
use anyhow::{anyhow, bail, Context};
use derive_new::new;
use itertools::Itertools;
use log::debug;
use rust_htslib::bam::{self, Read};
use rustc_hash::FxHashMap;
use std::path::{Path, PathBuf};

#[derive(Debug)]
pub(crate) struct SamplingSchedule {
    counts_for_chroms: FxHashMap<u32, usize>,
    index_stats: FxHashMap<i64, CountsFrac>,
    pub(crate) unmapped_count: usize,
}

#[derive(new, Debug)]
struct CountsFrac {
    counts: usize,
    frac: f32,
}

impl SamplingSchedule {
    fn get_index_stats(
        reader: &mut bam::IndexedReader,
        region: Option<&Region>,
        position_filter: Option<&StrandedPositionFilter>,
        include_unmapped: bool,
    ) -> anyhow::Result<FxHashMap<i64, CountsFrac>> {
        let stats = IdxStats::new_from_reader(reader, region, position_filter)?;
        let total = if include_unmapped {
            stats.total()
        } else {
            stats.mapped_read_count
        };
        if total == 0 {
            bail!("zero reads found in bam index")
        }
        let total = total as f32;
        let mut index_stats = stats
            .tid_to_mapped_read_count
            .into_iter()
            .map(|(target_id, n)| {
                (target_id, CountsFrac::new(n as usize, n as f32 / total))
            })
            .collect::<FxHashMap<i64, CountsFrac>>();
        if include_unmapped {
            let n = stats.unmapped_read_count;
            let unmapped_counts_frac =
                CountsFrac::new(n as usize, n as f32 / total);
            index_stats.insert(-1, unmapped_counts_frac);
        }
        Ok(index_stats)
    }

    fn log_schedule(
        counts_for_chroms: &FxHashMap<u32, usize>,
        unmapped_count: usize,
        total_to_sample: usize,
    ) {
        let contigs_with_reads = counts_for_chroms
            .iter()
            .filter(|(_, &counts)| counts > 0)
            .count();
        let noun = if contigs_with_reads > 1 {
            "contigs"
        } else {
            "contig"
        };
        debug!("derived sampling schedule, sampling total {total_to_sample} reads from \
            {} {noun}, {} unmapped reads", contigs_with_reads, unmapped_count);

        let report = counts_for_chroms
            .iter()
            .sorted_by(|(_, counts_a), (_, counts_b)| counts_b.cmp(counts_a))
            .fold(format!("schedule: "), |mut acc, (contig, counts)| {
                acc.push_str(&format!("SQ: {}, {} reads ", contig, counts));
                acc
            });
        debug!("{report}");
    }

    pub fn from_num_reads<T: AsRef<Path>>(
        bam_fp: T,
        num_reads: usize,
        region: Option<&Region>,
        position_filter: Option<&StrandedPositionFilter>,
        include_unmapped: bool,
    ) -> anyhow::Result<Self> {
        let mut reader = bam::IndexedReader::from_path(bam_fp)?;
        let index_stats = Self::get_index_stats(
            &mut reader,
            region,
            position_filter,
            include_unmapped,
        )?;
        let mut total_to_sample = 0usize;
        let counts_for_chroms = index_stats
            .iter()
            .filter(|(&chrom_id, counts_frac)| {
                chrom_id >= 0 && counts_frac.counts > 0
            })
            .map(|(&chrom_id, counts_frac)| {
                // use ceil here so that if there is at least 1 read aligned to this
                // contig we sample it.
                let num_reads_for_chrom = std::cmp::min(
                    (num_reads as f32 * counts_frac.frac).ceil() as usize,
                    counts_frac.counts,
                );
                total_to_sample += num_reads_for_chrom;
                (chrom_id as u32, num_reads_for_chrom)
            })
            .collect::<FxHashMap<u32, usize>>();
        let unmapped_count = if include_unmapped {
            index_stats
                .get(&-1)
                .map(|counts_frac| {
                    (num_reads as f32 * counts_frac.frac).ceil() as usize
                })
                .unwrap_or(0)
        } else {
            0
        };
        total_to_sample += unmapped_count;
        Self::log_schedule(&counts_for_chroms, unmapped_count, total_to_sample);

        Ok(Self {
            index_stats,
            counts_for_chroms,
            unmapped_count,
        })
    }

    pub fn from_sample_frac<T: AsRef<Path>>(
        bam_fp: T,
        sample_frac: f32,
        region: Option<&Region>,
        position_filter: Option<&StrandedPositionFilter>,
        include_unmapped: bool,
    ) -> anyhow::Result<Self> {
        if sample_frac > 1.0 {
            bail!("sample fraction must be <= 1")
        }
        let mut reader = bam::IndexedReader::from_path(bam_fp)?;
        let index_stats = Self::get_index_stats(
            &mut reader,
            region,
            position_filter,
            include_unmapped,
        )?;
        let mut total_to_sample = 0usize;
        let counts_for_chroms = index_stats
            .iter()
            .filter(|(&chrom_id, counts_frac)| {
                chrom_id >= 0 && counts_frac.counts > 0
            })
            .map(|(&chrom_id, counts_frac)| {
                let n_reads =
                    (counts_frac.counts as f32 * sample_frac).ceil() as usize;
                total_to_sample += n_reads;
                (chrom_id as u32, n_reads)
            })
            .collect::<FxHashMap<u32, usize>>();
        let unmapped_count = if include_unmapped {
            index_stats
                .get(&-1)
                .map(|couts_frac| {
                    (couts_frac.counts as f32 * sample_frac).ceil() as usize
                })
                .unwrap_or(0usize)
        } else {
            0usize
        };
        total_to_sample += unmapped_count;

        Self::log_schedule(&counts_for_chroms, unmapped_count, total_to_sample);
        Ok(Self {
            index_stats,
            counts_for_chroms,
            unmapped_count,
        })
    }

    pub(crate) fn get_num_reads(&self, chrom_id: u32) -> usize {
        *self.counts_for_chroms.get(&chrom_id).unwrap_or(&0)
    }

    pub(crate) fn get_num_reads_for_interval(
        &self,
        reference_record: &ReferenceRecord,
        total_interval_length: u32,
        start: u32,
        end: u32,
    ) -> usize {
        let reads_for_chrom = self.get_num_reads(reference_record.tid);
        let reads_aligned_to_chrom = self
            .index_stats
            .get(&(reference_record.tid as i64))
            .map(|counts_frac| counts_frac.counts)
            .unwrap_or(0);
        if reads_for_chrom == 0 || reads_aligned_to_chrom == 0 {
            return 0;
        }

        if reads_for_chrom == reads_aligned_to_chrom {
            reads_aligned_to_chrom
        } else {
            let f = (end - start) as f64 / total_interval_length as f64;
            let nr = reads_for_chrom as f64 * f;
            nr.ceil() as usize
        }
    }
}
pub(crate) struct IdxStats {
    tid_to_mapped_read_count: FxHashMap<i64, u64>,
    pub(crate) unmapped_read_count: u64,
    pub(crate) mapped_read_count: u64,
}

impl IdxStats {
    pub(crate) fn new_from_path(
        bam_fp: &PathBuf,
        region: Option<&Region>,
        position_filter: Option<&StrandedPositionFilter>,
    ) -> anyhow::Result<Self> {
        let mut reader = bam::IndexedReader::from_path(bam_fp)
            .context("could not create reader for getting mapping stats")?;
        Self::new_from_reader(&mut reader, region, position_filter)
    }

    pub(crate) fn new_from_reader(
        reader: &mut bam::IndexedReader,
        region: Option<&Region>,
        position_filter: Option<&StrandedPositionFilter>,
    ) -> anyhow::Result<Self> {
        let header = reader.header();
        let region_tid = region
            .map(|region| {
                let tid_for_region = (0..header.target_count())
                    .filter_map(|tid| {
                        String::from_utf8(header.tid2name(tid).to_vec())
                            .ok()
                            .map(|name| (tid, name))
                    })
                    .find_map(|(tid, name)| {
                        if name == region.name {
                            Some(tid)
                        } else {
                            None
                        }
                    });
                match tid_for_region {
                    None => Err(anyhow!(
                        "did not find target_id for region {} in header",
                        region.to_string()
                    )),
                    Some(tid_for_region) => Ok(tid_for_region),
                }
            })
            .transpose()?;

        let idx_stats =
            reader.index_stats().context("failed to get index stats")?;
        let mut mapped_read_count = 0u64;
        let mut unmapped_read_count = 0u64;
        let tid_to_mapped_read_count = idx_stats
            .into_iter()
            .filter(|(target_id, _, _, _)| {
                match (region_tid, position_filter) {
                    (Some(tid), _) => tid == *target_id as u32,
                    (None, Some(position_filter)) => {
                        position_filter.contains_chrom_id(target_id)
                    }
                    (None, None) => true,
                }
            })
            .filter_map(|(target_id, _length, n_mapped, n_unmapped)| {
                if target_id >= 0 {
                    mapped_read_count += n_mapped;
                    unmapped_read_count += n_unmapped;
                    Some((target_id, n_mapped))
                } else {
                    unmapped_read_count += n_unmapped;
                    None
                }
            })
            .collect::<FxHashMap<i64, u64>>();

        Ok(Self {
            tid_to_mapped_read_count,
            unmapped_read_count,
            mapped_read_count,
        })
    }

    pub(crate) fn total(&self) -> u64 {
        self.mapped_read_count + self.unmapped_read_count
    }
}

#[cfg(test)]
mod record_sampler_tests {
    use crate::reads_sampler::sampling_schedule::SamplingSchedule;

    #[test]
    fn test_record_sampler_sampling_schedule() {
        let sched = SamplingSchedule::from_num_reads(
            "tests/resources/bc_anchored_10_reads.sorted.bam",
            1023,
            None,
            None,
            false,
        )
        .unwrap();
        assert_eq!(sched.counts_for_chroms.get(&0), Some(&10));
        let sched = SamplingSchedule::from_num_reads(
            "tests/resources/bc_anchored_10_reads.sorted.bam",
            5,
            None,
            None,
            false,
        )
        .unwrap();
        assert_eq!(sched.counts_for_chroms.get(&0), Some(&5));
        let sched = SamplingSchedule::from_sample_frac(
            "tests/resources/bc_anchored_10_reads.sorted.bam",
            1f32,
            None,
            None,
            false,
        )
        .unwrap();
        assert_eq!(sched.counts_for_chroms.get(&0), Some(&10));
        let sched = SamplingSchedule::from_sample_frac(
            "tests/resources/bc_anchored_10_reads.sorted.bam",
            0.5f32,
            None,
            None,
            false,
        )
        .unwrap();
        assert_eq!(sched.counts_for_chroms.get(&0), Some(&5));
    }
}
