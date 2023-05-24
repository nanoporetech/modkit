use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::thread;

use bio::io::fasta::Reader as FastaReader;
use clap::Args;
use crossbeam_channel::{bounded, Sender};
use derive_new::new;
use indicatif::{MultiProgress, ParallelProgressIterator, ProgressIterator};
use log::{debug, error, info};
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use rust_htslib::bam::{self, FetchDefinition, Read};
use rust_lapper as lapper;

use crate::errs::RunError;
use crate::interval_chunks::IntervalChunks;
use crate::logging::init_logging;
use crate::mod_bam::{CollapseMethod, TrackingModRecordIter};
use crate::mod_base_code::ModCode;
use crate::read_ids_to_base_mod_probs::{
    ModProfile, ReadBaseModProfile, ReadsBaseModProfile,
};
use crate::reads_sampler::record_sampler::RecordSampler;
use crate::reads_sampler::sample_reads_from_interval;
use crate::record_processor::WithRecords;
use crate::util::{
    get_master_progress_bar, get_spinner, get_subroutine_progress_bar,
    get_targets, ReferenceRecord, Region, Strand,
};
use crate::writers::{
    OutwriterWithMemory, TsvWriter, TsvWriterWithContigNames,
};

#[derive(Args)]
pub struct ExtractMods {
    /// Path to modBAM file to extract read-level information from.
    in_bam: PathBuf,
    /// Path to output file, "stdout" or "-" will direct output to stdout.
    out_path: String,
    /// Number of threads to use
    #[arg(short = 't', long, default_value_t = 4)]
    threads: usize,
    /// Path to file to write run log.
    #[arg(long, alias = "log")]
    log_filepath: Option<PathBuf>,
    /// Include only mapped bases in output (alias: mapped).
    #[arg(long, alias = "mapped", default_value_t = false)]
    mapped_only: bool,
    /// Number of reads to use
    #[arg(long)]
    num_reads: Option<usize>,
    /// Process only reads that are aligned to a specified region of the BAM.
    /// Format should be <chrom_name>:<start>-<end> or <chrom_name>.
    #[arg(long)]
    region: Option<String>,
    /// Force overwrite of output file
    #[arg(long, default_value_t = false)]
    force: bool,

    /// Path to reference FASTA to extract reference context information from.
    /// If no reference is provided, `ref_kmer` column will be "." in the output.
    #[arg(long, alias = "ref")]
    reference: Option<PathBuf>,

    /// BED file with regions to include (alias: include).
    #[arg(long, alias = "include")]
    include_bed: Option<PathBuf>,
    /// BED file with regions to _exclude_ (alias: exclude).
    #[arg(long, alias = "exclude", short = 'v')]
    exclude_bed: Option<PathBuf>,

    /// Ignore a modified base class  _in_situ_ by redistributing base modification
    /// probability equally across other options. For example, if collapsing 'h',
    /// with 'm' and canonical options, half of the probability of 'h' will be added to
    /// both 'm' and 'C'. A full description of the methods can be found in
    /// collapse.md.
    #[arg(long, hide_short_help = true)]
    ignore: Option<char>,

    /// Interval chunk size in base pairs to process concurrently. Smaller interval
    /// chunk sizes will use less memory but incur more overhead. Only used when an
    /// indexed modBAM is provided.
    #[arg(
        short = 'i',
        long,
        default_value_t = 100_000,
        hide_short_help = true
    )]
    interval_size: u32,
}

type ReferenceAndIntervals = Vec<(ReferenceRecord, IntervalChunks)>;
type Iv = lapper::Interval<u64, ()>;
type GenomeLapper = lapper::Lapper<u64, ()>;

impl ExtractMods {
    fn load_regions(
        &self,
        name_to_tid: &HashMap<&str, u32>,
        region: Option<&Region>,
    ) -> anyhow::Result<(Option<ReferenceAndIntervals>, ReferencePositionFilter)>
    {
        let include_positions = self
            .include_bed
            .as_ref()
            .map(|fp| StrandedPositionFilter::from_bed_file(fp, name_to_tid))
            .transpose()?;
        let exclude_positions = self
            .exclude_bed
            .as_ref()
            .map(|fp| StrandedPositionFilter::from_bed_file(fp, name_to_tid))
            .transpose()?;

        let reference_and_intervals =
            match bam::IndexedReader::from_path(&self.in_bam) {
                Ok(reader) => {
                    let reference_records =
                        get_targets(reader.header(), region);
                    let reference_and_intervals = reference_records
                        .into_iter()
                        .map(|reference_record| {
                            let interval_chunks = IntervalChunks::new(
                                reference_record.start,
                                reference_record.length,
                                self.interval_size,
                                reference_record.tid,
                                None,
                            );
                            (reference_record, interval_chunks)
                        })
                        .collect::<ReferenceAndIntervals>();
                    Some(reference_and_intervals)
                }
                Err(_) => {
                    info!(
                    "did not find index to modBAM, defaulting to serial scan"
                );
                    None
                }
            };
        let reference_position_filter = ReferencePositionFilter::new(
            include_positions,
            exclude_positions,
            !self.mapped_only,
        );

        Ok((reference_and_intervals, reference_position_filter))
    }

    pub(crate) fn run(&self) -> anyhow::Result<()> {
        let _handle = init_logging(self.log_filepath.as_ref());

        let pool =
            ThreadPoolBuilder::new().num_threads(self.threads).build()?;

        let collapse_method = match &self.ignore {
            Some(raw_mod_code) => {
                let _ = ModCode::parse_raw_mod_code(*raw_mod_code)?;
                Some(CollapseMethod::ReDistribute(*raw_mod_code))
            }
            None => None,
        };
        let mut reader = bam::Reader::from_path(&self.in_bam)?;
        let header = reader.header().to_owned();

        let (snd, rcv) = bounded(100_000);

        let in_bam = self.in_bam.clone();
        let tid_to_name = (0..header.target_count())
            .filter_map(|tid| {
                match String::from_utf8(header.tid2name(tid).to_vec()) {
                    Ok(contig) => Some((tid, contig)),
                    Err(e) => {
                        error!(
                            "failed to parse contig {tid}, {}",
                            e.to_string()
                        );
                        None
                    }
                }
            })
            .collect::<HashMap<u32, String>>();
        let name_to_tid = tid_to_name
            .iter()
            .map(|(tid, name)| (name.as_str(), *tid))
            .collect::<HashMap<&str, u32>>();

        let chrom_to_seq = match self.reference.as_ref() {
            Some(fp) => {
                let reader = FastaReader::from_file(fp)?;
                let pb = get_spinner();
                pb.set_message("parsing FASTA records");
                reader
                    .records()
                    .progress_with(pb)
                    .filter_map(|r| r.ok())
                    .filter(|record| name_to_tid.get(record.id()).is_some())
                    .map(|record| {
                        (record.id().to_owned(), record.seq().to_vec())
                    })
                    .collect::<HashMap<String, Vec<u8>>>()
            }
            None => HashMap::new(),
        };

        let region = self
            .region
            .as_ref()
            .map(|raw_region| Region::parse_str(raw_region, reader.header()))
            .transpose()?;

        let (regions, reference_position_filter) =
            self.load_regions(&name_to_tid, region.as_ref())?;

        let multi_prog = MultiProgress::new();
        // multi_prog.set_draw_target(indicatif::ProgressDrawTarget::hidden());
        let n_failed = multi_prog.add(get_spinner());
        n_failed.set_message("~records failed");
        let n_skipped = multi_prog.add(get_spinner());
        n_skipped.set_message("~records skipped");
        let n_used = multi_prog.add(get_spinner());
        n_used.set_message("~records used");
        let n_rows = multi_prog.add(get_spinner());
        n_rows.set_message("rows written");
        reader.set_threads(self.threads)?;
        let n_reads = self.num_reads;
        let threads = self.threads;

        thread::spawn(move || {
            pool.install(|| {
                if let Some(interval_chunks) = regions {
                    drop(reader);
                    let master_progress = multi_prog.add(get_master_progress_bar(interval_chunks.len()));
                    master_progress.set_message("contigs");

                    let total_length = interval_chunks.iter().map(|(r, _)| r.length as usize).sum::<usize>();
                    let mut num_aligned_reads_used = 0usize;
                    for (reference_record, interval_chunks) in interval_chunks {
                        let interval_chunks =
                            interval_chunks.collect::<Vec<(u32, u32)>>();

                        // todo use the index api to pre-calculate the number of reads per ref
                        let n_reads_for_reference = if let Some(nr) = n_reads {
                            let f = reference_record.length as f64 / total_length as f64;
                            let nr = nr as f64 * f;
                            Some(std::cmp::max(nr.floor() as usize, 1usize))
                        } else {
                            None
                        };

                        let interval_pb = multi_prog.add(get_subroutine_progress_bar(interval_chunks.len()));
                        interval_pb.set_message(format!("processing {}", &reference_record.name));
                        let n_reads_used = interval_chunks.into_par_iter()
                            .progress_with(interval_pb)
                            .map(
                                |(start, end)| {
                                    let record_sampler = n_reads_for_reference.map(|nr| {
                                        let f = (end - start) as f64 / reference_record.length as f64;
                                        let nr = nr as f64 * f;
                                        let nr = std::cmp::max(nr.floor() as usize, 1usize);
                                        RecordSampler::new_num_reads(nr)
                                    }).unwrap_or(RecordSampler::new_passthrough());

                                    let batch_result = sample_reads_from_interval::<
                                        ReadsBaseModProfile,
                                    >(
                                        &in_bam,
                                        reference_record.tid,
                                        start,
                                        end,
                                        record_sampler,
                                        collapse_method.as_ref(),
                                    ).map(|reads_base_mod_profile| {
                                        reference_position_filter.filter_read_base_mod_probs(reads_base_mod_profile)
                                    });
                                    let num_reads_success = batch_result.as_ref().map(|batch| batch.num_reads()).unwrap_or(0);

                                    match snd.send(batch_result) {
                                        Ok(_) => {
                                            num_reads_success
                                        }
                                        Err(e) => {
                                            error!( "failed to send result to writer, {}", e.to_string() );
                                            0
                                        }
                                    }
                                }
                            ).sum::<usize>();
                        num_aligned_reads_used += n_reads_used;
                    }

                    if reference_position_filter.include_unmapped {
                        let n_unmapped_reads = n_reads.map(|nr| {
                            nr - num_aligned_reads_used
                        });
                        if let Some(n) = n_unmapped_reads {
                            debug!("processing {n} unmapped reads");
                        } else {
                            debug!("processing unmapped reads");
                        }
                        let reader = bam::IndexedReader::from_path(&in_bam)
                            .and_then(|mut reader| reader.fetch(FetchDefinition::Unmapped).map(|_| reader))
                            .and_then(|mut reader| reader.set_threads(threads).map(|_| reader));
                        match reader {
                            Ok(mut reader) => {
                                let (skip, fail) = Self::process_records_to_chan(
                                    reader.records(),
                                    &multi_prog,
                                    &reference_position_filter,
                                    snd.clone(),
                                    n_unmapped_reads,
                                    collapse_method.as_ref());
                                let _ = snd.send(Ok(ReadsBaseModProfile::new(Vec::new(), skip, fail)));
                            },
                            Err(e) => {
                                error!("failed to get indexed reader for unmapped read processing, {}", e.to_string());
                            }
                        }
                    }
                } else {
                    let (skip, fail) = Self::process_records_to_chan(
                        reader.records(),
                        &multi_prog,
                        &reference_position_filter,
                        snd.clone(),
                        n_reads,
                        collapse_method.as_ref()
                    );
                    let _ = snd.send(Ok(ReadsBaseModProfile::new(Vec::new(), skip, fail)));
                }
            })
        });

        let mut writer: Box<dyn OutwriterWithMemory<ReadsBaseModProfile>> =
            match self.out_path.as_str() {
                "stdout" | "-" => {
                    let tsv_writer =
                        TsvWriter::new_stdout(Some(ModProfile::header()));
                    let writer = TsvWriterWithContigNames::new(
                        tsv_writer,
                        tid_to_name,
                        chrom_to_seq,
                        HashSet::new(),
                    );
                    Box::new(writer)
                }
                _ => {
                    let tsv_writer = TsvWriter::new_file(
                        &self.out_path,
                        self.force,
                        Some(ModProfile::header()),
                    )?;
                    let writer = TsvWriterWithContigNames::new(
                        tsv_writer,
                        tid_to_name,
                        chrom_to_seq,
                        HashSet::new(),
                    );
                    Box::new(writer)
                }
            };

        for result in rcv {
            match result {
                Ok(mod_profile) => {
                    n_used.inc(mod_profile.num_reads() as u64);
                    n_failed.inc(mod_profile.num_fails as u64);
                    n_skipped.inc(mod_profile.num_skips as u64);
                    match writer.write(mod_profile) {
                        Ok(n) => n_rows.inc(n),
                        Err(e) => {
                            error!("failed to write {}", e.to_string());
                        }
                    }
                }
                Err(e) => {
                    debug!(
                        "failed to calculate read-level mod probs, {}",
                        e.to_string()
                    );
                }
            }
        }
        n_failed.finish_and_clear();
        n_skipped.finish_and_clear();
        n_used.finish_and_clear();
        n_rows.finish_and_clear();
        info!(
            "processed {} reads, {} rows, skipped {} reads, failed {} reads",
            writer.num_reads(),
            n_rows.position(),
            n_skipped.position(),
            n_failed.position()
        );
        Ok(())
    }

    fn process_records_to_chan<'a, T: Read>(
        records: bam::Records<T>,
        multi_pb: &MultiProgress,
        reference_position_filter: &ReferencePositionFilter,
        snd: Sender<anyhow::Result<ReadsBaseModProfile>>,
        n_reads: Option<usize>,
        collapse_method: Option<&CollapseMethod>,
    ) -> (usize, usize) {
        let mut mod_iter = TrackingModRecordIter::new(records);
        let pb = multi_pb.add(get_spinner());
        pb.set_message("records processed");
        for (record, read_id, mod_base_info) in &mut mod_iter {
            let mod_profile = match ReadBaseModProfile::process_record(
                &record,
                &read_id,
                mod_base_info,
                collapse_method,
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
            let mod_profile = reference_position_filter
                .filter_read_base_mod_probs(mod_profile);
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
            let done = n_reads
                .map(|nr| pb.position() as usize >= nr)
                .unwrap_or(false);
            if done {
                // pb.finish_and_clear();
                debug!("stopping after processing {} reads", pb.position());
                break;
            }
        }
        pb.finish_and_clear();
        (mod_iter.num_skipped, mod_iter.num_failed)
    }
}

struct StrandedPositionFilter {
    pos_positions: HashMap<u32, GenomeLapper>,
    neg_positions: HashMap<u32, GenomeLapper>,
}

impl StrandedPositionFilter {
    fn from_bed_file(
        bed_fp: &PathBuf,
        chrom_to_target_id: &HashMap<&str, u32>,
    ) -> anyhow::Result<Self> {
        info!(
            "parsing BED at {}",
            bed_fp.to_str().unwrap_or("invalid-UTF-8")
        );

        let fh = File::open(bed_fp)?;
        let mut pos_positions = HashMap::new();
        let mut neg_positions = HashMap::new();
        let lines_processed = get_spinner();
        let mut warned = HashSet::new();

        let reader = BufReader::new(fh);
        for line in reader.lines().filter_map(|l| l.ok()) {
            let parts = line.split_ascii_whitespace().collect::<Vec<&str>>();
            let chrom_name = parts[0];
            if warned.contains(chrom_name) {
                continue;
            }
            if parts.len() < 6 {
                info!("improperly formatted BED line {line}");
                continue;
            }
            let raw_start = &parts[1].parse::<u64>();
            let raw_end = &parts[2].parse::<u64>();
            let (start, stop) = match (raw_start, raw_end) {
                (Ok(start), Ok(end)) => (*start, *end),
                _ => {
                    info!("improperly formatted BED line {line}");
                    continue;
                }
            };
            let (pos_strand, neg_strand) = match parts[5] {
                "+" => (true, false),
                "-" => (false, true),
                "." => (true, true),
                _ => {
                    info!("improperly formatted strand field {}", &parts[5]);
                    continue;
                }
            };
            if let Some(chrom_id) = chrom_to_target_id.get(chrom_name) {
                if pos_strand {
                    pos_positions.entry(*chrom_id).or_insert(Vec::new()).push(
                        Iv {
                            start,
                            stop,
                            val: (),
                        },
                    )
                }
                if neg_strand {
                    neg_positions.entry(*chrom_id).or_insert(Vec::new()).push(
                        Iv {
                            start,
                            stop,
                            val: (),
                        },
                    )
                }
                lines_processed.inc(1);
            } else {
                info!("skipping chrom {chrom_name}, not present in BAM header");
                warned.insert(chrom_name.to_owned());
                continue;
            }
        }

        let pos_lapper = pos_positions
            .into_iter()
            .map(|(chrom_id, intervals)| {
                let mut lp = lapper::Lapper::new(intervals);
                lp.merge_overlaps();
                (chrom_id, lp)
            })
            .collect::<HashMap<u32, GenomeLapper>>();

        let neg_lapper = neg_positions
            .into_iter()
            .map(|(chrom_id, intervals)| {
                let mut lp = lapper::Lapper::new(intervals);
                lp.merge_overlaps();
                (chrom_id, lp)
            })
            .collect::<HashMap<u32, GenomeLapper>>();

        lines_processed.finish_and_clear();
        info!("processed {} BED lines", lines_processed.position());
        Ok(Self {
            pos_positions: pos_lapper,
            neg_positions: neg_lapper,
        })
    }

    fn contains(&self, chrom_id: u32, position: u64, strand: Strand) -> bool {
        let positions = match strand {
            Strand::Positive => &self.pos_positions,
            Strand::Negative => &self.neg_positions,
        };
        positions
            .get(&chrom_id)
            .map(|lp| lp.find(position, position + 1).count() > 0)
            .unwrap_or(false)
    }
}

#[derive(new)]
struct ReferencePositionFilter {
    include_pos: Option<StrandedPositionFilter>,
    exclude_pos: Option<StrandedPositionFilter>,
    include_unmapped: bool,
}

impl ReferencePositionFilter {
    fn keep(&self, chrom_id: u32, position: u64, strand: Strand) -> bool {
        let include_hit = self
            .include_pos
            .as_ref()
            .map(|flt| flt.contains(chrom_id, position, strand))
            .unwrap_or(true);
        let exclude_hit = self
            .exclude_pos
            .as_ref()
            .map(|filt| filt.contains(chrom_id, position, strand))
            .unwrap_or(false);

        include_hit && !exclude_hit
    }

    fn filter_read_base_mod_probs(
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
                                self.keep(chrom_id, ref_pos as u64, strand)
                            }
                            _ => self.include_unmapped,
                        }
                    })
                    .collect::<Vec<ModProfile>>();
                ReadBaseModProfile::new(read_name, chrom_id, profile)
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
