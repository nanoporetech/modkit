use crate::position_filter::StrandedPositionFilter;
use crate::reads_sampler::record_sampler::RecordSampler;
use crate::util::{reader_is_bam, ReferenceRecord, Region};
use anyhow::{anyhow, bail, Context};
use derive_new::new;
use itertools::Itertools;
use log::debug;
use rust_htslib::bam::{self, FetchDefinition, Read};
use rustc_hash::FxHashMap;
use std::cmp::Ordering;
use std::path::{Path, PathBuf};

/// Count is an exact count, Sample is a fraction to sample
#[derive(Debug, PartialEq)]
enum CountOrSample {
    Count(usize),
    Sample(f32),
}

impl Eq for CountOrSample {}

impl PartialOrd for CountOrSample {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        match (self, other) {
            (CountOrSample::Count(x), CountOrSample::Count(y)) => {
                Some(x.cmp(y))
            }
            (CountOrSample::Sample(x), CountOrSample::Sample(y)) => {
                x.partial_cmp(y)
            }
            _ => None,
        }
    }
}

#[derive(Debug)]
pub(crate) struct SamplingSchedule {
    counts_for_chroms: FxHashMap<u32, CountOrSample>,
    unmapped_count: Option<CountOrSample>,
}

#[derive(new, Debug)]
struct CountsFrac {
    counts: usize,
    frac: f32,
}

impl SamplingSchedule {
    fn get_contig_to_counts_frac(
        index_stats: IdxStats,
        include_unmapped: bool,
    ) -> anyhow::Result<FxHashMap<i64, CountsFrac>> {
        debug_assert!(index_stats.is_bam, "should be BAM-index based");
        let total = if include_unmapped {
            index_stats.total()
        } else {
            index_stats.mapped_read_count
        };
        if total == 0 {
            bail!("zero reads found in bam index")
        }
        let total = total as f32;
        let mut contig_to_counts_frac = index_stats
            .tid_to_mapped_read_count
            .into_iter()
            .map(|(target_id, n)| {
                (target_id, CountsFrac::new(n as usize, n as f32 / total))
            })
            .collect::<FxHashMap<i64, CountsFrac>>();
        if include_unmapped {
            let n = index_stats.unmapped_read_count;
            let unmapped_counts_frac =
                CountsFrac::new(n as usize, n as f32 / total);
            contig_to_counts_frac.insert(-1, unmapped_counts_frac);
        }
        Ok(contig_to_counts_frac)
    }

    fn log_schedule(
        is_bam: bool,
        counts_for_chroms: &FxHashMap<u32, CountOrSample>,
        unmapped_count: Option<&CountOrSample>,
        total_to_sample: CountOrSample,
    ) {
        if !is_bam {
            debug!("using CRAM index, sampling schedule is approximate!");
        }
        let contigs_with_reads = counts_for_chroms.len();
        let noun = if contigs_with_reads > 1 {
            "contigs"
        } else {
            "contig"
        };
        let total_to_sample = match total_to_sample {
            CountOrSample::Count(count) => format!("{}", count),
            CountOrSample::Sample(frac) => format!("{}% of", frac * 100f32),
        };
        let unmapped = match unmapped_count {
            Some(CountOrSample::Count(0)) => format!("including"),
            Some(CountOrSample::Count(count)) => format!("{}", count),
            Some(CountOrSample::Sample(s)) => {
                format!("{}% of", (*s * 100f32).round())
            }
            None => format!("0"),
        };
        debug!("derived sampling schedule, sampling total {total_to_sample} reads from \
            {} {noun}, {} unmapped reads", contigs_with_reads, unmapped);

        let report = counts_for_chroms
            .iter()
            .sorted_by(|(_, counts_a), (_, counts_b)| {
                counts_b.partial_cmp(counts_a).unwrap_or(Ordering::Equal)
            })
            .fold(format!("schedule: "), |mut acc, (contig, counts)| {
                acc.push_str(&format!("SQ: {}, {:?} reads ", contig, counts));
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
        let header = reader.header().to_owned();
        let index_stats =
            IdxStats::new_from_reader(&mut reader, region, position_filter)?;
        drop(reader); // mostly as a safety because we leave the reader in a mutated state
        if index_stats.is_bam {
            let contig_to_counts_frac =
                Self::get_contig_to_counts_frac(index_stats, include_unmapped)?;
            let mut total_to_sample = 0usize;

            let counts_for_chroms = contig_to_counts_frac
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
                    (chrom_id as u32, CountOrSample::Count(num_reads_for_chrom))
                })
                .collect::<FxHashMap<u32, CountOrSample>>();

            let unmapped_count = if include_unmapped {
                let count = contig_to_counts_frac
                    .get(&-1)
                    .map(|counts_frac| {
                        (num_reads as f32 * counts_frac.frac).ceil() as usize
                    })
                    .unwrap_or(0);
                total_to_sample += count;
                Some(CountOrSample::Count(count))
            } else {
                None
            };
            Self::log_schedule(
                true,
                &counts_for_chroms,
                unmapped_count.as_ref(),
                CountOrSample::Count(total_to_sample),
            );

            Ok(Self {
                counts_for_chroms,
                unmapped_count,
            })
        } else {
            // using CRAM distribute num_reads over the contigs that we found at least 1
            // record for (N.B. that we assume the target IDs here are are >=0
            let contigs_with_records = index_stats
                .tid_to_mapped_read_count
                .keys()
                .copied()
                .collect::<Vec<i64>>();

            let total_length = contigs_with_records
                .iter()
                .filter_map(|&target_id| header.target_len(target_id as u32))
                .sum::<u64>() as f32;

            let counts_for_chroms = contigs_with_records
                .iter()
                .filter_map(|&target_id| {
                    header.target_len(target_id as u32).map(|length| {
                        // get the proportion of the total length contributed by this contig
                        let weight = length as f32 / total_length;
                        let count = (num_reads as f32 * weight).ceil() as usize;
                        (target_id as u32, CountOrSample::Count(count))
                    })
                })
                .collect::<FxHashMap<u32, CountOrSample>>();

            let unmapped_count = if include_unmapped {
                if index_stats.unmapped_read_count > 0 {
                    Some(CountOrSample::Count(0))
                } else {
                    None
                }
            } else {
                None
            };
            Self::log_schedule(
                false,
                &counts_for_chroms,
                unmapped_count.as_ref(),
                CountOrSample::Count(num_reads),
            );
            Ok(Self {
                counts_for_chroms,
                unmapped_count,
            })
        }
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
        let index_stats =
            IdxStats::new_from_reader(&mut reader, region, position_filter)?;
        drop(reader);
        if index_stats.is_bam {
            let contig_to_counts_frac =
                Self::get_contig_to_counts_frac(index_stats, include_unmapped)?;
            let mut total_to_sample = 0usize;
            let counts_for_chroms = contig_to_counts_frac
                .iter()
                .filter(|(&chrom_id, counts_frac)| {
                    chrom_id >= 0 && counts_frac.counts > 0
                })
                .map(|(&chrom_id, counts_frac)| {
                    if sample_frac == 1.0f32 {
                        total_to_sample += counts_frac.counts;
                        (chrom_id as u32, CountOrSample::Sample(1.0))
                    } else {
                        let n_reads = (counts_frac.counts as f32 * sample_frac)
                            .ceil()
                            as usize;
                        total_to_sample += n_reads;
                        (chrom_id as u32, CountOrSample::Count(n_reads))
                    }
                })
                .collect::<FxHashMap<u32, CountOrSample>>();
            let unmapped_count = if include_unmapped {
                let count = contig_to_counts_frac
                    .get(&-1)
                    .map(|couts_frac| {
                        (couts_frac.counts as f32 * sample_frac).ceil() as usize
                    })
                    .unwrap_or(0usize);
                total_to_sample += count;
                Some(CountOrSample::Count(count))
            } else {
                None
            };

            Self::log_schedule(
                true,
                &counts_for_chroms,
                unmapped_count.as_ref(),
                CountOrSample::Count(total_to_sample),
            );
            Ok(Self {
                counts_for_chroms,
                unmapped_count,
            })
        } else {
            let counts_for_chroms = index_stats
                .tid_to_mapped_read_count
                .iter()
                .map(|(&target_id, _)| {
                    (target_id as u32, CountOrSample::Sample(sample_frac))
                })
                .collect::<FxHashMap<u32, CountOrSample>>();

            let unmapped_count = if include_unmapped {
                if index_stats.unmapped_read_count > 0 {
                    Some(CountOrSample::Sample(sample_frac))
                } else {
                    None
                }
            } else {
                None
            };
            Self::log_schedule(
                false,
                &counts_for_chroms,
                unmapped_count.as_ref(),
                CountOrSample::Sample(sample_frac),
            );
            Ok(Self {
                counts_for_chroms,
                unmapped_count,
            })
        }
    }

    pub(crate) fn chrom_has_reads(&self, chrom_id: u32) -> bool {
        self.counts_for_chroms.contains_key(&chrom_id)
    }

    // pub(crate) fn get_num_reads_for_interval(
    //     &self,
    //     reference_record: &ReferenceRecord,
    //     total_interval_length: u32,
    //     start: u32,
    //     end: u32,
    // ) -> usize {
    //     let reads_for_chrom = self.get_num_reads(reference_record.tid);
    //     let reads_aligned_to_chrom = self
    //         .index_stats
    //         .get(&(reference_record.tid as i64))
    //         .map(|counts_frac| counts_frac.counts)
    //         .unwrap_or(0);
    //     if reads_for_chrom == 0 || reads_aligned_to_chrom == 0 {
    //         return 0;
    //     }
    //
    //     if reads_for_chrom == reads_aligned_to_chrom {
    //         reads_aligned_to_chrom
    //     } else {
    //         let f = (end - start) as f64 / total_interval_length as f64;
    //         let nr = reads_for_chrom as f64 * f;
    //         nr.ceil() as usize
    //     }
    // }

    pub(crate) fn get_record_sampler(
        &self,
        reference_record: &ReferenceRecord,
        total_interval_length: u32,
        start: u32,
        end: u32,
    ) -> RecordSampler {
        self.counts_for_chroms
            .get(&reference_record.tid)
            .map(|counts_or_sample| match counts_or_sample {
                CountOrSample::Count(count) => {
                    let f = (end - start) as f64 / total_interval_length as f64;
                    let nr = *count as f64 * f;
                    RecordSampler::new_num_reads(nr.ceil() as usize)
                }
                CountOrSample::Sample(frac) => {
                    RecordSampler::new_sample_frac(*frac as f64, None)
                }
            })
            .unwrap_or_else(|| RecordSampler::new_sample_frac(0.0, None))
    }

    pub(crate) fn has_unmapped(&self) -> bool {
        self.unmapped_count.is_some()
    }
}

pub(crate) struct IdxStats {
    tid_to_mapped_read_count: FxHashMap<i64, u64>,
    is_bam: bool,
    pub(crate) unmapped_read_count: u64,
    pub(crate) mapped_read_count: u64,
}

impl IdxStats {
    fn check_is_bam_format(reader: &bam::IndexedReader) -> bool {
        reader_is_bam(reader)
    }

    pub(crate) fn check_any_mapped_reads(
        bam_fp: &PathBuf,
        region: Option<&Region>,
        position_filter: Option<&StrandedPositionFilter>,
    ) -> anyhow::Result<bool> {
        Self::new_from_path(bam_fp, region, position_filter)
            .map(|idx_stats| idx_stats.mapped_read_count > 0)
    }

    fn new_from_path(
        bam_fp: &PathBuf,
        region: Option<&Region>,
        position_filter: Option<&StrandedPositionFilter>,
    ) -> anyhow::Result<Self> {
        let mut reader = bam::IndexedReader::from_path(bam_fp)
            .context("could not create reader for getting mapping stats")?;
        Self::new_from_reader(&mut reader, region, position_filter)
    }

    fn new_from_reader(
        reader: &mut bam::IndexedReader,
        region: Option<&Region>,
        position_filter: Option<&StrandedPositionFilter>,
    ) -> anyhow::Result<Self> {
        let header = reader.header();
        // get the tid of the region if we're using it
        let region_tid = region
            .map(|region| {
                // todo use header directly here
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

        let is_bam = Self::check_is_bam_format(&reader);
        if is_bam {
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
                is_bam,
                unmapped_read_count,
                mapped_read_count,
            })
        } else {
            // when using CRAM the semantics are that the target ids must have at least 1 found read
            let tid_to_mapped_read_count = (0..header.target_count())
                .filter(|&target_id| match (region_tid, position_filter) {
                    (Some(tid), _) => tid == target_id,
                    (None, Some(position_filter)) => {
                        position_filter.contains_chrom_id(&(target_id as i64))
                    }
                    (None, None) => true,
                })
                .filter_map(|target_id| {
                    if let Err(e) = reader.fetch(FetchDefinition::CompleteTid(target_id as i32)) {
                        debug!("failed to fetch contig {target_id}, {}", e.to_string());
                        None
                    } else {
                        loop {
                            match reader.records().next() {
                                Some(Ok(_)) => break Some(target_id),
                                None => break None,
                                Some(Err(e)) => {
                                    debug!("failed to parse record for target id {target_id}, {}", e.to_string());
                                    break None
                                }
                            }
                        }
                    }
                })
                .map(|target_id| (target_id as i64, 1u64)).collect::<FxHashMap<i64, u64>>();
            reader
                .fetch(FetchDefinition::Unmapped)
                .context("failed to fetch unmapped reads")?;
            let has_unmapped = loop {
                match reader.records().next() {
                    Some(Ok(_)) => break true,
                    None => break false,
                    Some(Err(e)) => {
                        debug!(
                            "failed to parse unmapped record, {}",
                            e.to_string()
                        );
                        break false;
                    }
                }
            };

            Ok(Self {
                tid_to_mapped_read_count,
                is_bam,
                unmapped_read_count: if has_unmapped { 1 } else { 0 },
                mapped_read_count: 1,
            })
        }
    }

    pub(crate) fn total(&self) -> u64 {
        self.mapped_read_count + self.unmapped_read_count
    }
}

#[cfg(test)]
mod record_sampler_tests {
    use crate::reads_sampler::sampling_schedule::{
        CountOrSample, SamplingSchedule,
    };

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
        assert_eq!(
            sched.counts_for_chroms.get(&0),
            Some(&CountOrSample::Count(10))
        );
        let sched = SamplingSchedule::from_num_reads(
            "tests/resources/bc_anchored_10_reads.sorted.bam",
            5,
            None,
            None,
            false,
        )
        .unwrap();
        assert_eq!(
            sched.counts_for_chroms.get(&0),
            Some(&CountOrSample::Count(5))
        );
        let sched = SamplingSchedule::from_sample_frac(
            "tests/resources/bc_anchored_10_reads.sorted.bam",
            1f32,
            None,
            None,
            false,
        )
        .unwrap();
        assert_eq!(
            sched.counts_for_chroms.get(&0),
            Some(&CountOrSample::Sample(1.0))
        );

        let sched = SamplingSchedule::from_sample_frac(
            "tests/resources/bc_anchored_10_reads.sorted.bam",
            0.5f32,
            None,
            None,
            false,
        )
        .unwrap();
        assert_eq!(
            sched.counts_for_chroms.get(&0),
            Some(&CountOrSample::Count(5))
        );
    }

    #[test]
    #[ignore = "no reference server, TODO: figure out workaround for CI"]
    fn test_record_sampler_sampling_schedule_cram() {
        let sched = SamplingSchedule::from_num_reads(
            "tests/resources/bc_anchored_10_reads.sorted.cram",
            1023,
            None,
            None,
            false,
        )
        .unwrap();
        assert_eq!(
            sched.counts_for_chroms.get(&0),
            // will over-ask because CRAM doesn't have all of the counts
            Some(&CountOrSample::Count(1023))
        );
        assert_eq!(sched.counts_for_chroms.len(), 1);
        let sched = SamplingSchedule::from_sample_frac(
            "tests/resources/bc_anchored_10_reads.sorted.cram",
            1.0,
            None,
            None,
            false,
        )
        .unwrap();
        assert_eq!(
            sched.counts_for_chroms.get(&0),
            Some(&CountOrSample::Sample(1.0))
        );
        assert_eq!(sched.counts_for_chroms.len(), 1);
    }

    #[test]
    #[ignore = "no reference server, TODO: figure out workaround for CI"]
    fn test_record_sampler_sampling_schedule_cram_unmapped() {
        let sched = SamplingSchedule::from_num_reads(
            "tests/resources/bc_anchored_10_reads_unmapped.sorted.cram",
            1023,
            None,
            None,
            true,
        )
        .unwrap();
        assert_eq!(sched.unmapped_count, Some(CountOrSample::Count(0)));
        let sched = SamplingSchedule::from_sample_frac(
            "tests/resources/bc_anchored_10_reads_unmapped.sorted.cram",
            0.05,
            None,
            None,
            true,
        )
        .unwrap();
        assert_eq!(sched.unmapped_count, Some(CountOrSample::Sample(0.05)));
    }
}
