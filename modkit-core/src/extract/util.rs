use crate::errs::{MkError, MkResult};
use crate::extract::args::InputArgs;
use crate::interval_chunks::{
    ReferenceIntervalBatchesFeeder, TotalLength, WithPrevEnd,
};
use crate::mod_bam::{CollapseMethod, EdgeFilter, TrackingModRecordIter};
use crate::mod_base_code::ModCodeRepr;
use crate::monoid::Moniod;
use crate::motifs::motif_bed::{
    find_motif_hits, MotifPositionLookup, RegexMotif,
};
use crate::pileup::base_mods_adapter::BaseModsAdapter;
use crate::position_filter::{GenomeIntervals, Iv, StrandedPositionFilter};
use crate::read_ids_to_base_mod_probs::{
    ModProfile, ReadBaseModProfile, ReadsBaseModProfile,
};
use crate::reads_sampler::record_sampler::RecordSampler;
use crate::reads_sampler::sample_reads_from_interval_with_reference;
use crate::reads_sampler::sampling_schedule::SamplingSchedule;
use crate::record_processor::WithRecords;
use crate::sample_probs::QualHist;
use crate::util::{
    get_guage, get_master_progress_bar, get_query_name_string,
    get_reference_mod_strand, get_subroutine_progress_bar, get_targets,
    get_ticker, record_is_primary, set_reference_for_cram_indexed_reader,
    set_reference_for_cram_reader, Region, Strand,
};
use derive_new::new;
use indicatif::{MultiProgress, ParallelProgressIterator, ProgressBar};
use itertools::Itertools;
use log::{debug, error, info};
use rayon::prelude::*;
use rayon::ThreadPool;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{self, FetchDefinition, Read};
use rustc_hash::FxHashMap;
use std::collections::HashMap;
use std::path::{Path, PathBuf};

#[derive(new)]
pub(super) struct ReferencePositionFilter {
    pub(super) include_pos: Option<StrandedPositionFilter<()>>,
    pub(super) exclude_pos: Option<StrandedPositionFilter<()>>,
    pub(super) include_unmapped_reads: bool,
    pub(super) include_unmapped_positions: bool,
}

impl ReferencePositionFilter {
    fn keep(
        &self,
        chrom_id: u32,
        position: u64,
        alignment_strand: Strand,
        mod_strand: Strand,
    ) -> bool {
        let reference_mod_strand =
            get_reference_mod_strand(mod_strand, alignment_strand);
        let include_hit = self
            .include_pos
            .as_ref()
            .map(|flt| {
                flt.contains(chrom_id as i32, position, reference_mod_strand)
            })
            .unwrap_or(true);
        let exclude_hit = self
            .exclude_pos
            .as_ref()
            .map(|filt| {
                filt.contains(chrom_id as i32, position, reference_mod_strand)
            })
            .unwrap_or(false);

        include_hit && !exclude_hit
    }

    pub(super) fn filter_read_base_mod_probs(
        &self,
        reads_base_mods_profile: ReadsBaseModProfile,
    ) -> ReadsBaseModProfile {
        let mut n_skipped = reads_base_mods_profile.num_skips;
        let n_failed = reads_base_mods_profile.num_fails;
        let profiles = reads_base_mods_profile
            .profiles
            .into_par_iter()
            .map(|read_base_mod_profile| {
                let read_name = read_base_mod_profile.record_name;
                let chrom_id = read_base_mod_profile.chrom_id;
                let flag = read_base_mod_profile.flag;
                let alignment_start = read_base_mod_profile.alignment_start;
                let alignment_end = read_base_mod_profile.alignment_end;
                let profile = read_base_mod_profile
                    .profile
                    .into_par_iter()
                    .filter(|mod_profile| {
                        match (
                            chrom_id,
                            mod_profile.ref_position,
                            mod_profile.alignment_strand,
                        ) {
                            (Some(chrom_id), Some(ref_pos), Some(strand)) => {
                                self.keep(
                                    chrom_id,
                                    ref_pos as u64,
                                    strand,
                                    mod_profile.mod_strand,
                                )
                            }
                            _ => self.include_unmapped_positions,
                        }
                    })
                    .collect::<Vec<ModProfile>>();
                ReadBaseModProfile::new(
                    read_name,
                    chrom_id,
                    flag,
                    alignment_start,
                    alignment_end,
                    profile,
                )
            })
            .collect::<Vec<ReadBaseModProfile>>();
        let empty = profiles
            .iter()
            .filter(|read_base_mod_profile| {
                read_base_mod_profile.profile.is_empty()
            })
            .count();
        n_skipped += empty;
        ReadsBaseModProfile::new(profiles, n_skipped, n_failed)
    }
}

pub(super) fn load_regions(
    input_args: &InputArgs,
    using_stdin: bool,
    name_to_tid: &HashMap<&str, u32>,
    region: Option<&Region>,
    contigs: &HashMap<String, Vec<u8>>,
    master_progress_bar: &MultiProgress,
    thread_pool: &ThreadPool,
) -> anyhow::Result<(
    Option<ReferenceIntervalBatchesFeeder>,
    ReferencePositionFilter,
    Option<MotifPositionLookup>,
)> {
    let (include_unmapped_reads, include_unmapped_positions) = if input_args
        .include_bed
        .is_some()
    {
        info!("specifying include-only BED outputs only mapped sites");
        (false, false)
    } else if input_args.motif.is_some() || input_args.cpg {
        info!("specifying a motif (including --cpg) outputs only mapped sites");
        (false, false)
    } else if region.is_some() {
        info!("specifying a region outputs only mapped reads");
        if input_args.mapped_only {
            info!("including only mapped positions");
        } else {
            info!("including unmapped positions within mapped reads");
        }
        (false, !input_args.mapped_only)
    } else {
        (!input_args.mapped_only, !input_args.mapped_only)
    };

    let motifs = if let Some(raw_motif_parts) = &input_args.motif {
        Some(RegexMotif::from_raw_parts(&raw_motif_parts, input_args.cpg)?)
    } else if input_args.cpg {
        Some(vec![RegexMotif::parse_string("CG", 0).unwrap()])
    } else {
        None
    };

    let include_positions = input_args
        .include_bed
        .as_ref()
        .map(|fp| {
            StrandedPositionFilter::from_bed_file(
                fp,
                name_to_tid,
                input_args.suppress_progress,
            )
        })
        .transpose()?;

    let exclude_positions = input_args
        .exclude_bed
        .as_ref()
        .map(|fp| {
            StrandedPositionFilter::from_bed_file(
                fp,
                name_to_tid,
                input_args.suppress_progress,
            )
        })
        .transpose()?;

    // extract the motif positions, if given
    let tid_motif_to_positions = motifs.as_ref().map(|motifs| {
        let pb =
            master_progress_bar.add(get_subroutine_progress_bar(contigs.len()));
        master_progress_bar.suspend(|| {
            info!("searching for {} motifs", motifs.len());
        });
        pb.set_message("contigs searched");
        let contigs_sorted_by_size = contigs
            .iter()
            .sorted_by(|(_, s), (_, p)| s.len().cmp(&p.len()))
            .collect::<Vec<(&String, &Vec<u8>)>>();
        thread_pool.install(|| {
            contigs_sorted_by_size
                .into_par_iter()
                .progress_with(pb)
                .filter_map(|(name, raw_seq)| {
                    name_to_tid.get(name.as_str()).map(|tid| (*tid, raw_seq))
                })
                .map(|(tid, raw_seq)| {
                    let seq =
                        raw_seq.iter().map(|&b| b as char).collect::<String>();
                    let seq = if input_args.mask {
                        seq
                    } else {
                        seq.to_ascii_uppercase()
                    };
                    motifs
                        .par_iter()
                        .enumerate()
                        .map(|(motif_idx, motif)| {
                            let positions = find_motif_hits(&seq, motif);
                            let positions = if let Some(filter) =
                                include_positions.as_ref()
                            {
                                positions
                                    .into_iter()
                                    .filter(|(pos, strand)| {
                                        filter.contains(
                                            tid as i32,
                                            *pos as u64,
                                            *strand,
                                        )
                                    })
                                    .collect::<Vec<(usize, Strand)>>()
                            } else {
                                positions
                            };
                            ((tid, motif_idx), positions)
                        })
                        .collect::<FxHashMap<(u32, usize), Vec<(usize, Strand)>>>(
                        )
                })
                .reduce(|| FxHashMap::zero(), |a, b| a.op(b))
        })
    });

    // intersect the motif positions with the include positions from the BED
    // file
    let to_lapper = |intervals: Vec<(Iv, Strand)>| -> GenomeIntervals<()> {
        let intervals = intervals.into_iter().map(|(iv, _)| iv).collect();
        GenomeIntervals::new(intervals)
    };

    let include_positions =
        match (tid_motif_to_positions.as_ref(), input_args.annotate_motifs) {
            (Some(tid_motif_to_pos), false) => {
                let tid_to_positions = tid_motif_to_pos.iter().fold(
                    FxHashMap::default(),
                    |mut agg, ((tid, _motif_idx), positions)| {
                        let ps = agg.entry(*tid).or_insert_with(Vec::new);
                        ps.extend(positions.iter().map(|(pos, strand)| {
                            let iv = Iv {
                                start: *pos as u64,
                                stop: (*pos).saturating_add(1) as u64,
                                val: (),
                            };
                            (iv, *strand)
                        }));
                        agg
                    },
                );

                let (pos_lappers, neg_lappers) =
                    tid_to_positions.into_iter().fold(
                        (FxHashMap::default(), FxHashMap::default()),
                        |(mut pos, mut neg), (tid, intervals)| {
                            let (pos_intervals, neg_intervals): (
                                Vec<(Iv, Strand)>,
                                Vec<(Iv, Strand)>,
                            ) = intervals.into_iter().partition(
                                |(_iv, strand)| *strand == Strand::Positive,
                            );
                            let mut pos_lapper = to_lapper(pos_intervals);
                            let mut neg_lapper = to_lapper(neg_intervals);
                            pos_lapper.merge_overlaps();
                            neg_lapper.merge_overlaps();
                            pos.insert(tid, pos_lapper);
                            neg.insert(tid, neg_lapper);
                            (pos, neg)
                        },
                    );

                Some(StrandedPositionFilter {
                    pos_positions: pos_lappers,
                    neg_positions: neg_lappers,
                })
            }
            (Some(_), true) => {
                master_progress_bar.suspend(|| {
                    info!("annotating motifs, but not restricting output");
                });
                include_positions
            }
            _ => include_positions,
        };

    let reference_and_intervals = if !using_stdin && !input_args.ignore_index {
        match bam::IndexedReader::from_path(&input_args.in_bam) {
            Ok(reader) => {
                info!(
                    "found BAM index, processing reads in {} base pair chunks",
                    input_args.interval_size
                );
                let reference_records = get_targets(reader.header(), region);
                let zero_reference_targets = reference_records.is_empty();
                let reference_records =
                    if let Some(pf) = include_positions.as_ref() {
                        pf.optimize_reference_records(
                            reference_records,
                            input_args.interval_size,
                        )
                    } else {
                        reference_records
                    };

                let batch_size =
                    (input_args.threads as f32 * 1.5f32).floor() as usize;
                let feeder = if include_unmapped_reads && zero_reference_targets
                {
                    ReferenceIntervalBatchesFeeder::
                        new_allowing_zero_reference_targets(
                            reference_records,
                            batch_size,
                            input_args.interval_size,
                            false,
                            None,
                            None,
                        )?
                } else {
                    ReferenceIntervalBatchesFeeder::new(
                        reference_records,
                        batch_size,
                        input_args.interval_size,
                        false,
                        None,
                        None,
                    )?
                };
                Some(feeder)
            }
            Err(_) => {
                info!(
                    "did not find index to modBAM, defaulting to serial scan"
                );
                None
            }
        }
    } else {
        None
    };

    let motif_lookup = match (tid_motif_to_positions, motifs) {
        (Some(tid_motif_positions), Some(motifs)) => {
            Some(MotifPositionLookup::new(tid_motif_positions, motifs))
        }
        _ => None,
    };

    let reference_position_filter = ReferencePositionFilter::new(
        include_positions,
        exclude_positions,
        include_unmapped_reads,
        include_unmapped_positions,
    );

    Ok((reference_and_intervals, reference_position_filter, motif_lookup))
}

pub(super) fn run_extract_reads(
    mut reader: bam::Reader,
    in_bam: String,
    reference_fasta: Option<PathBuf>,
    references_and_intervals: Option<ReferenceIntervalBatchesFeeder>,
    schedule: Option<SamplingSchedule>,
    collapse_method: Option<CollapseMethod>,
    edge_filter: Option<EdgeFilter>,
    allow_non_primary: bool,
    kmer_size: usize,
    remove_inferred: bool,
    reference_position_filter: ReferencePositionFilter,
    snd: crossbeam::channel::Sender<anyhow::Result<ReadsBaseModProfile>>,
    queue_size: usize,
    n_reads: Option<usize>,
    threads: usize,
    mapped_only: bool,
    multi_prog: MultiProgress,
) {
    let gauge = multi_prog.add(get_guage(queue_size));
    gauge.set_message("enqueued processed reads");
    gauge.set_position(snd.len() as u64);
    // references_and_intervals is only some when we have an index
    if let Some(feeder) = references_and_intervals {
        drop(reader);
        // should make this a method on this struct?
        let bam_fp = Path::new(&in_bam).to_path_buf();

        let prog_length = feeder.total_length() as usize;
        let master_progress =
            multi_prog.add(get_master_progress_bar(prog_length));
        master_progress.set_message("genome positions");

        let mut num_aligned_reads_used = 0usize;
        // TODO/WARN currently this cannot fail since we don't use the
        // FastaIndex here
        let feeder = feeder.map(|x| x.unwrap());
        for super_batch in feeder.with_prev_end() {
            let total_batch_length =
                super_batch.iter().map(|c| c.total_length()).sum::<u64>();
            let batch_progress =
                multi_prog.add(get_subroutine_progress_bar(super_batch.len()));
            batch_progress.set_message("batch progress");
            let n_reads_used = super_batch
                .into_par_iter()
                .progress_with(batch_progress)
                .map(|batch| {
                    let successful_reads_in_batch = batch
                        .0
                        .into_par_iter()
                        .filter(|cc| {
                            schedule
                                .as_ref()
                                .map(|s| s.chrom_has_reads(cc.chrom_tid()))
                                .unwrap_or(true)
                        })
                        .map(|cc| {
                            let record_sampler = schedule
                                .as_ref()
                                .map(|s| {
                                    s.get_record_sampler(
                                        cc.chrom_tid(),
                                        total_batch_length as u32,
                                        cc.start_pos(),
                                        cc.end_pos(),
                                    )
                                })
                                .unwrap_or_else(|| {
                                    RecordSampler::new_passthrough()
                                });
                            let batch_result =
                                sample_reads_from_interval_with_reference::<
                                    ReadsBaseModProfile,
                                >(
                                    &bam_fp,
                                    reference_fasta.as_ref(),
                                    cc.chrom_tid(),
                                    cc.start_pos(),
                                    cc.end_pos(),
                                    cc.prev_end(),
                                    record_sampler,
                                    collapse_method.as_ref(),
                                    edge_filter.as_ref(),
                                    None,
                                    false,
                                    allow_non_primary,
                                    Some(kmer_size),
                                );
                            let batch_result = batch_result
                                .map(|reads_base_mod_profile| {
                                    if remove_inferred {
                                        reads_base_mod_profile.remove_inferred()
                                    } else {
                                        reads_base_mod_profile
                                    }
                                })
                                .map(|reads_base_mod_profile| {
                                    reference_position_filter
                                        .filter_read_base_mod_probs(
                                            reads_base_mod_profile,
                                        )
                                });

                            let num_reads_success = batch_result
                                .as_ref()
                                .map(|r| r.num_reads())
                                .unwrap_or(0);

                            match snd.send(batch_result) {
                                Ok(_) => {
                                    let n_enqueued = snd.len() as u64;
                                    gauge.set_position(n_enqueued);
                                    num_reads_success
                                }
                                Err(e) => {
                                    error!(
                                        "failed to send result to writer, {}",
                                        e.to_string()
                                    );
                                    0
                                }
                            }
                        })
                        .sum::<usize>();
                    successful_reads_in_batch
                })
                .sum::<usize>();
            num_aligned_reads_used += n_reads_used;
            master_progress.inc(total_batch_length);
        }

        if reference_position_filter.include_unmapped_reads {
            let n_unmapped_reads = n_reads
                .map(|nr| nr.checked_sub(num_aligned_reads_used).unwrap_or(0));
            if let Some(n) = n_unmapped_reads {
                debug!("processing {n} unmapped reads");
            } else {
                debug!("processing unmapped reads");
            }
            let reader = (|| -> anyhow::Result<bam::IndexedReader> {
                let mut reader = bam::IndexedReader::from_path(&bam_fp)?;
                set_reference_for_cram_indexed_reader(
                    &mut reader,
                    reference_fasta.as_ref(),
                )?;
                reader.fetch(FetchDefinition::Unmapped)?;
                reader.set_threads(threads)?;
                Ok(reader)
            })();
            match reader {
                Ok(mut reader) => {
                    let (skip, fail) = process_records_to_chan(
                        reader.records(),
                        &multi_prog,
                        &reference_position_filter,
                        snd.clone(),
                        n_unmapped_reads,
                        collapse_method.as_ref(),
                        edge_filter.as_ref(),
                        false,
                        false,
                        "unmapped ",
                        kmer_size,
                    );
                    let _ = snd.send(Ok(ReadsBaseModProfile::new(
                        Vec::new(),
                        skip,
                        fail,
                    )));
                }
                Err(e) => {
                    error!(
                        "failed to get indexed reader for unmapped read \
                         processing, {}",
                        e.to_string()
                    );
                }
            }
        }
    } else {
        let (skip, fail) = process_records_to_chan(
            reader.records(),
            &multi_prog,
            &reference_position_filter,
            snd.clone(),
            n_reads,
            collapse_method.as_ref(),
            edge_filter.as_ref(),
            mapped_only,
            allow_non_primary,
            "",
            kmer_size,
        );
        let _ = snd.send(Ok(ReadsBaseModProfile::new(Vec::new(), skip, fail)));
    }
}

fn process_records_to_chan<'a, T: Read>(
    records: bam::Records<T>,
    multi_pb: &MultiProgress,
    reference_position_filter: &ReferencePositionFilter,
    snd: crossbeam::channel::Sender<anyhow::Result<ReadsBaseModProfile>>,
    n_reads: Option<usize>,
    collapse_method: Option<&CollapseMethod>,
    edge_filter: Option<&EdgeFilter>,
    only_mapped: bool,
    allow_non_primary: bool,
    message: &'static str,
    kmer_size: usize,
) -> (usize, usize) {
    let mut mod_iter =
        TrackingModRecordIter::new(records, false, allow_non_primary);
    let pb = multi_pb.add(get_ticker());
    pb.set_message(format!("{message}records processed"));
    for (record, read_id, mod_base_info) in &mut mod_iter {
        if record.is_unmapped() && only_mapped {
            continue;
        }
        let mod_profile = match ReadBaseModProfile::process_record(
            &record,
            &read_id,
            mod_base_info,
            collapse_method,
            edge_filter,
            kmer_size,
        ) {
            Ok(mod_profile) => {
                ReadsBaseModProfile::new(vec![mod_profile], 0, 0)
            }
            Err(_) => ReadsBaseModProfile::new(Vec::new(), 0, 1),
        };
        let mod_profile =
            reference_position_filter.filter_read_base_mod_probs(mod_profile);
        match snd.send(Ok(mod_profile)) {
            Ok(_) => {
                pb.inc(1);
            }
            Err(snd_error) => {
                error!(
                    "failed to send results to writer, {}",
                    snd_error.to_string()
                );
            }
        }
        let done =
            n_reads.map(|nr| pb.position() as usize >= nr).unwrap_or(false);
        if done {
            debug!("stopping after processing {} reads", pb.position());
            break;
        }
    }
    pb.finish_and_clear();
    (mod_iter.num_skipped, mod_iter.num_failed)
}

pub(super) struct ReadModsStatsRecord<const SIZE: usize> {
    pub read_id: String,
    pub chrom_id: i32,
    pub aln_start: i64,
    pub count_unmodified: [u32; 4],
    pub count_unmodified_fail: [u32; 4],
    pub other_modified: [u32; 4],
    pub other_modified_fail: [u32; 4],
    pub mod_counts: [u32; SIZE],
    pub filtered_mod_counts: [u32; SIZE],
    pub read_length: usize,
}

impl<const SIZE: usize> ReadModsStatsRecord<SIZE> {
    fn zero(&mut self) {
        self.count_unmodified = [0u32; 4];
        self.count_unmodified_fail = [0u32; 4];
        self.other_modified = [0u32; 4];
        self.other_modified_fail = [0u32; 4];
        self.mod_counts = [0u32; SIZE];
        self.filtered_mod_counts = [0u32; SIZE];
        self.read_length = 0usize;
    }

    pub(super) fn empty() -> Self {
        Self {
            read_id: "".to_string(),
            chrom_id: -1i32,
            aln_start: -1i64,
            count_unmodified: [0u32; 4],
            count_unmodified_fail: [0u32; 4],
            other_modified: [0u32; 4],
            other_modified_fail: [0u32; 4],
            mod_counts: [0u32; SIZE],
            filtered_mod_counts: [0u32; SIZE],
            read_length: 0usize,
        }
    }

    fn new(
        record: &bam::Record,
        mut mem: Self,
        mods_per_base: &[u8; 4],
        mods_codes_and_idx: &[(ModCodeRepr, usize)],
        filter_thresholds: [f32; 4],
        mod_thresholds: &[(ModCodeRepr, f32)],
    ) -> MkResult<Self> {
        mem.zero();
        let read_id = get_query_name_string(record)?;
        let ref_start = record.reference_start();
        let chrom_id = record.tid();
        let mut adapter = BaseModsAdapter::<SIZE>::new(record)
            .map_err(|_| MkError::InvalidTags)?;
        loop {
            match adapter
                .next_modified_position(filter_thresholds, mod_thresholds)
            {
                Ok(Some(mod_pos)) => {
                    let idx = mod_pos.primary_base as usize;
                    if mods_per_base[idx] > 0 {
                        if mod_pos.modified {
                            if let Some(i) =
                                mods_codes_and_idx.iter().find_map(|(c, x)| {
                                    if *c == mod_pos.mod_code {
                                        Some(*x)
                                    } else {
                                        None
                                    }
                                })
                            {
                                let ar = if mod_pos.filtered {
                                    &mut mem.filtered_mod_counts
                                } else {
                                    &mut mem.mod_counts
                                };
                                ar[i] = ar[i].saturating_add(1);
                            } else {
                                if mod_pos.filtered {
                                    let count = mem.other_modified_fail[idx]
                                        .saturating_add(1);
                                    mem.other_modified_fail[idx] = count;
                                } else {
                                    let count = mem.other_modified[idx]
                                        .saturating_add(1);
                                    mem.other_modified[idx] = count;
                                }
                            }
                        } else {
                            let ar = if mod_pos.filtered {
                                &mut mem.count_unmodified_fail
                            } else {
                                &mut mem.count_unmodified
                            };
                            ar[idx] = ar[idx].saturating_add(1);
                        }
                    } else {
                        continue;
                    }
                }
                Ok(None) => {
                    break;
                }
                Err(e) => return Err(e),
            }
        }
        mem.read_length = record.seq_len();
        mem.read_id = read_id;
        mem.chrom_id = chrom_id;
        mem.aln_start = ref_start;
        return Ok(mem);
    }
}

pub(super) struct ReadModStatsProcessor<const SIZE: usize> {
    mods_per_base: [u8; 4],
    mod_codes_and_idx: Vec<(ModCodeRepr, usize)>,
    reader: bam::Reader,
    filter_thresholds: [f32; 4],
    mod_thresholds: Vec<(ModCodeRepr, f32)>,
    failed_records: ProgressBar,
    skipped_records: ProgressBar,
    in_channel: crossbeam::channel::Receiver<ReadModsStatsRecord<SIZE>>,
    out_channel:
        crossbeam::channel::Sender<MkResult<ReadModsStatsRecord<SIZE>>>,
}

impl<const SIZE: usize> ReadModStatsProcessor<SIZE> {
    pub(super) fn new(
        in_channel: crossbeam::channel::Receiver<ReadModsStatsRecord<SIZE>>,
        out_channel: crossbeam::channel::Sender<
            MkResult<ReadModsStatsRecord<SIZE>>,
        >,
        reader: bam::Reader,
        mods_per_base: [u8; 4],
        mod_codes: &[ModCodeRepr],
        mod_code_idxs: &[usize],
        filter_thresholds: [f32; 4],
        mod_thresholds: Vec<(ModCodeRepr, f32)>,
        multi_prog: MultiProgress,
    ) -> Self {
        assert_eq!(mod_codes.len(), mod_code_idxs.len());
        let mod_codes_and_idx = mod_codes
            .iter()
            .zip(mod_code_idxs.iter())
            .map(|(code, idx)| (*code, *idx))
            .collect();
        let failed_records = multi_prog.add(get_ticker());
        failed_records.set_message("failed records");
        let skipped_records = multi_prog.add(get_ticker());
        skipped_records.set_message("skipped records");

        Self {
            in_channel,
            mods_per_base,
            mod_codes_and_idx,
            reader,
            filter_thresholds,
            mod_thresholds,
            out_channel,
            failed_records,
            skipped_records,
        }
    }

    pub(super) fn process_bam(&mut self) {
        let get_mem = || -> Result<ReadModsStatsRecord<SIZE>, ()> {
            match self.in_channel.recv() {
                Ok(mem) => Ok(mem),
                Err(_) => Err(()),
            }
        };
        let iter = self
            .reader
            .records()
            .map(|r| r.map_err(|e| MkError::from(e)))
            .filter_ok(|rec| {
                if record_is_primary(&rec) {
                    true
                } else {
                    self.skipped_records.inc(1);
                    false
                }
            });
        'records: for res in iter {
            match res {
                Ok(record) => match get_mem() {
                    Ok(mem) => {
                        let stats_res = ReadModsStatsRecord::new(
                            &record,
                            mem,
                            &self.mods_per_base,
                            &self.mod_codes_and_idx,
                            self.filter_thresholds,
                            &self.mod_thresholds,
                        );
                        if let Err(_) = self.out_channel.send(stats_res) {
                            break 'records;
                        }
                    }
                    Err(_) => {
                        self.failed_records.inc(1);
                        break 'records;
                    }
                },
                Err(e) => {
                    error!("reading from BAM failed, {e}");
                    break 'records;
                }
            }
        }
    }
}

pub(super) fn calc_per_base_thresholds_from_stream(
    bam_fp: &PathBuf,
    reference_fasta: Option<&PathBuf>,
    num_reads: usize,
    allow_non_primary: bool,
    stranded_position_filter: Option<StrandedPositionFilter<()>>,
    edge_filter: Option<&EdgeFilter>,
    filter_percentile: f32,
    io_threads: usize,
    max_thresholds_per_base: Option<[f32; 4]>,
    multi_progress: &MultiProgress,
) -> anyhow::Result<[f32; 4]> {
    multi_progress.suspend(|| {
        info!(
            "estimating base thresholds from modBAM, taking first {num_reads} \
             including any unmapped records"
        );
    });
    let mut records = bam::Reader::from_path(bam_fp)?;
    set_reference_for_cram_reader(&mut records, reference_fasta)?;
    records.set_threads(io_threads)?;
    QualHist::from_records(
        records.records(),
        stranded_position_filter,
        Some(num_reads),
        None,
        None,
        edge_filter,
        false,
        allow_non_primary,
        false,
        multi_progress,
    )?
    .get_base_thresholds(
        filter_percentile,
        max_thresholds_per_base,
        multi_progress,
    )
}
