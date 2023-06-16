use crate::position_filter::StrandedPositionFilter;
use crate::util::{
    get_master_progress_bar, get_spinner, ReferenceRecord, Region,
};
use anyhow::bail;
use derive_new::new;
use indicatif::ProgressBar;
use itertools::Itertools;
use log::debug;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rustc_hash::FxHashMap;
use std::path::Path;

/// A utility data structure that when used in an interator allows
/// to randomly sample either a preset number of reads or a fraction
/// of reads. If sampling a preset number of reads, say N reads, then
/// the first N reads are taken.
pub struct RecordSampler {
    pub(crate) num_reads: Option<usize>,
    pub(crate) sample_frac: Option<f64>,
    rng: StdRng,
    reads_sampled: usize,
}

impl RecordSampler {
    pub(crate) fn new_num_reads(num_reads: usize) -> Self {
        Self {
            num_reads: Some(num_reads),
            sample_frac: None,
            rng: StdRng::from_entropy(),
            reads_sampled: 0,
        }
    }

    pub(crate) fn new_sample_frac(sample_frac: f64, seed: Option<u64>) -> Self {
        let rng = seed
            .map(|s| StdRng::seed_from_u64(s))
            .unwrap_or(StdRng::from_entropy());
        Self {
            num_reads: None,
            sample_frac: Some(sample_frac),
            rng,
            reads_sampled: 0,
        }
    }

    pub(crate) fn new_passthrough() -> Self {
        Self {
            num_reads: None,
            sample_frac: None,
            rng: StdRng::from_entropy(),
            reads_sampled: 0,
        }
    }

    pub(crate) fn new_from_options(
        sample_frac: Option<f64>,
        num_reads: Option<usize>,
        seed: Option<u64>,
    ) -> Self {
        match (sample_frac, num_reads) {
            (_, Some(num_reads)) => RecordSampler::new_num_reads(num_reads),
            (Some(f), _) => RecordSampler::new_sample_frac(f, seed),
            (None, None) => RecordSampler::new_passthrough(),
        }
    }

    pub(crate) fn get_progress_bar(&self) -> ProgressBar {
        let spinner = if let Some(num) = self.num_reads {
            get_master_progress_bar(num)
        } else {
            get_spinner()
        };
        spinner.set_message("records sampled");
        spinner
    }

    fn check_num_reads(&mut self) -> Indicator {
        if self.reads_sampled >= self.num_reads.unwrap() {
            Indicator::Done
        } else {
            Indicator::Use(Token)
        }
    }

    fn check_sample_frac(&mut self) -> Indicator {
        if self.rng.gen_bool(self.sample_frac.unwrap()) {
            Indicator::Use(Token)
        } else {
            Indicator::Skip
        }
    }

    pub(crate) fn ask(&mut self) -> Indicator {
        match (self.num_reads, self.sample_frac) {
            (Some(_nr), _) => self.check_num_reads(),
            (_, Some(_sample_frac)) => self.check_sample_frac(),
            (None, None) => Indicator::Use(Token),
        }
    }

    pub(crate) fn used(&mut self, _token: Token) {
        self.reads_sampled += 1;
    }
}

pub(crate) struct Token;

pub(crate) enum Indicator {
    Use(Token),
    Skip,
    Done,
}

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
        let header = reader.header();
        let region_tid = if let Some(region) = region {
            let tid_for_region =
                (0..header.target_count())
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
            if tid_for_region.is_none() {
                bail!(
                    "did not find target_id for region {} in header",
                    region.name
                )
            } else {
                tid_for_region // is always Some
            }
        } else {
            None
        };

        let idx_stats = reader.index_stats()?;
        let mut total = 0usize;
        let idx_stats = idx_stats
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
                    total += n_mapped as usize;
                    Some((target_id, n_mapped))
                } else if include_unmapped {
                    total += n_unmapped as usize;
                    Some((target_id, n_unmapped))
                } else {
                    None
                }
            })
            .collect::<Vec<(i64, u64)>>();
        if total == 0 {
            bail!("zero reads in bam index")
        }
        let total = total as f32;
        // could be sped up with pulp/simd
        Ok(idx_stats
            .into_iter()
            .map(|(target_id, n)| {
                (target_id, CountsFrac::new(n as usize, n as f32 / total))
            })
            .collect())
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

#[cfg(test)]
mod record_sampler_tests {
    use crate::reads_sampler::record_sampler::SamplingSchedule;

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
