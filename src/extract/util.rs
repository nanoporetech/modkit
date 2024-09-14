use crate::errs::RunError;
use crate::extract::args::InputArgs;
use crate::find_motifs::motif_bed::{find_motif_hits, RegexMotif};
use crate::interval_chunks::{ReferenceIntervalsFeeder, WithPrevEnd};
use crate::mod_bam::{CollapseMethod, EdgeFilter, TrackingModRecordIter};
use crate::monoid::Moniod;
use crate::position_filter::{GenomeIntervals, Iv, StrandedPositionFilter};
use crate::read_ids_to_base_mod_probs::{
    ModProfile, ReadBaseModProfile, ReadsBaseModProfile,
};
use crate::reads_sampler::record_sampler::RecordSampler;
use crate::reads_sampler::sample_reads_from_interval;
use crate::reads_sampler::sampling_schedule::SamplingSchedule;
use crate::record_processor::WithRecords;
use crate::util::{
    get_guage, get_master_progress_bar, get_reference_mod_strand,
    get_subroutine_progress_bar, get_targets, get_ticker, Region, Strand,
};
use derive_new::new;
use indicatif::{MultiProgress, ParallelProgressIterator};
use itertools::Itertools;
use log::{debug, error, info};
use rayon::prelude::*;
use rayon::ThreadPool;
use rust_htslib::bam::{self, FetchDefinition, Read};
use rustc_hash::FxHashMap;
use std::collections::HashMap;
use std::path::Path;

#[derive(new)]
pub(super) struct ReferencePositionFilter {
    pub(super) include_pos: Option<StrandedPositionFilter<()>>,
    pub(super) exclude_pos: Option<StrandedPositionFilter<()>>,
    pub(super) include_unmapped_reads: bool,
    pub(super) include_unmapped_positions: bool,
}

impl ReferencePositionFilter {
    pub(super) fn only_mapped_positions(&self) -> bool {
        !self.include_unmapped_positions
    }

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
) -> anyhow::Result<(Option<ReferenceIntervalsFeeder>, ReferencePositionFilter)>
{
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

    // intersect the motif positions with the include positions from the BED
    // file
    let include_positions = if let Some(motifs) = motifs {
        let pb =
            master_progress_bar.add(get_subroutine_progress_bar(contigs.len()));
        pb.set_message("contigs searched");
        let contigs_sorted_by_size = contigs
            .iter()
            .sorted_by(|(_, s), (_, p)| s.len().cmp(&p.len()))
            .collect::<Vec<(&String, &Vec<u8>)>>();
        let tid_to_positions = thread_pool.install(|| {
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
                        .map(|motif| {
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
                            (tid, positions)
                        })
                        .collect::<HashMap<u32, Vec<(usize, Strand)>>>()
                })
                .reduce(|| HashMap::zero(), |a, b| a.op(b))
        });
        let (pos_lappers, neg_lappers) = tid_to_positions.into_iter().fold(
            (FxHashMap::default(), FxHashMap::default()),
            |(mut pos, mut neg), (tid, positions)| {
                let to_lapper =
                    |intervals: Vec<(Iv, Strand)>| -> GenomeIntervals<()> {
                        let intervals =
                            intervals.into_iter().map(|(iv, _)| iv).collect();
                        GenomeIntervals::new(intervals)
                    };

                let (pos_positions, neg_positions): (
                    Vec<(Iv, Strand)>,
                    Vec<(Iv, Strand)>,
                ) = positions
                    .into_iter()
                    .map(|(position, strand)| {
                        let iv = Iv {
                            start: position as u64,
                            stop: (position + 1) as u64,
                            val: (),
                        };
                        (iv, strand)
                    })
                    .partition(|(_iv, strand)| *strand == Strand::Positive);
                let pos_lapper = to_lapper(pos_positions);
                let neg_lapper = to_lapper(neg_positions);
                pos.insert(tid, pos_lapper);
                neg.insert(tid, neg_lapper);
                (pos, neg)
            },
        );

        Some(StrandedPositionFilter {
            pos_positions: pos_lappers,
            neg_positions: neg_lappers,
        })
    } else {
        include_positions
    };

    let reference_and_intervals = if !using_stdin && !input_args.ignore_index {
        match bam::IndexedReader::from_path(&input_args.in_bam) {
            Ok(reader) => {
                info!(
                    "found BAM index, processing reads in {} base pair chunks",
                    input_args.interval_size
                );
                let reference_records = get_targets(reader.header(), region);
                let reference_records =
                    if let Some(pf) = include_positions.as_ref() {
                        pf.optimize_reference_records(
                            reference_records,
                            input_args.interval_size,
                        )
                    } else {
                        reference_records
                    };

                let feeder = ReferenceIntervalsFeeder::new(
                    reference_records,
                    (input_args.threads as f32 * 1.5f32).floor() as usize,
                    input_args.interval_size,
                    false,
                    None,
                    None,
                )?;
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

    let reference_position_filter = ReferencePositionFilter::new(
        include_positions,
        exclude_positions,
        include_unmapped_reads,
        include_unmapped_positions,
    );

    Ok((reference_and_intervals, reference_position_filter))
}

pub(super) fn run_extract_reads(
    mut reader: bam::Reader,
    in_bam: String,
    references_and_intervals: Option<ReferenceIntervalsFeeder>,
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

        let prog_length = feeder.total_length();
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
                            let batch_result = sample_reads_from_interval::<
                                ReadsBaseModProfile,
                            >(
                                &bam_fp,
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
                            )
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
            let reader = bam::IndexedReader::from_path(&bam_fp)
                .and_then(|mut reader| {
                    reader.fetch(FetchDefinition::Unmapped).map(|_| reader)
                })
                .and_then(|mut reader| {
                    reader.set_threads(threads).map(|_| reader)
                });
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
            Err(run_error) => match run_error {
                RunError::BadInput(_) | RunError::Failed(_) => {
                    ReadsBaseModProfile::new(Vec::new(), 0, 1)
                }
                RunError::Skipped(_) => {
                    ReadsBaseModProfile::new(Vec::new(), 1, 0)
                }
            },
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
