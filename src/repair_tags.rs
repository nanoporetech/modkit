use crate::errs::RunError;
use crate::logging::init_logging;
use crate::mod_bam::{
    format_mm_ml_tag, BaseModProbs, DeltaListConverter, ModBaseInfo,
    SeqPosBaseModProbs, ML_TAGS, MM_TAGS,
};
use crate::util::{
    get_forward_sequence, get_query_name_string, get_ticker,
    record_is_secondary,
};
use anyhow::{anyhow, bail, Context};
use clap::Args;
use derive_new::new;
use indicatif::{MultiProgress, ProgressBar};
use log::{debug, error, info, warn};
use rayon::prelude::*;
use rust_htslib::bam::record::{Aux, AuxArray};
use rust_htslib::bam::{self, Read};
use rustc_hash::FxHashMap;

use std::path::PathBuf;
use std::sync::Arc;

#[derive(Args)]
pub struct RepairTags {
    /// Donor modBAM with original MM/ML tags. Must be sorted by read name.
    #[arg(long, short = 'd', alias = "donor")]
    donor_bam: PathBuf,
    /// Acceptor modBAM with reads to have MM/ML base modification data
    /// projected on to. Must be sorted by read name.
    #[arg(long, short = 'a', alias = "acceptor")]
    acceptor_bam: PathBuf,
    /// output modBAM location.
    #[arg(long, short = 'o', alias = "output")]
    output_bam: PathBuf,
    /// File to write logs to, it is recommended to use this option as some
    /// reads may be rejected and logged here.
    #[arg(long)]
    log_filepath: Option<PathBuf>,
    /// The number of threads to use.
    #[arg(long, short = 't', default_value_t = 4)]
    threads: usize,
}

impl RepairTags {
    pub(crate) fn run(&self) -> anyhow::Result<()> {
        let _handle = init_logging(self.log_filepath.as_ref());

        let reader_threads = {
            let half = self.threads / 2;
            std::cmp::min(half, 16)
        };
        let threads_per_reader = std::cmp::max(reader_threads / 2, 1);
        let pool_threads =
            self.threads.checked_sub(reader_threads).unwrap_or(1);
        debug!(
            "assigning {threads_per_reader} to each reader and using \
             {pool_threads} to process records"
        );

        let (pair_snd, pair_rcv) = std::sync::mpsc::sync_channel(1000);
        let mut donor_records = bam::Reader::from_path(&self.donor_bam)?;
        donor_records.set_threads(threads_per_reader)?;
        let mut acceptor_records = bam::Reader::from_path(&self.acceptor_bam)?;
        acceptor_records.set_threads(threads_per_reader)?;
        let header = bam::Header::from_template(acceptor_records.header());
        let mut writer = bam::Writer::from_path(
            &self.output_bam,
            &header,
            bam::Format::Bam,
        )?;
        info!(
            "repairing records in {} with base modification information in {}",
            &self.acceptor_bam.to_str().unwrap_or_else(|| "??"),
            &self.donor_bam.to_str().unwrap_or_else(|| "??")
        );

        // pb stuff
        let master_progress = MultiProgress::new();
        let donor_ticker = master_progress.add(get_ticker());
        donor_ticker.set_message("~donor records processed");
        let acceptor_ticker = master_progress.add(get_ticker());
        acceptor_ticker.set_message("~acceptor records processed");
        let repaired_ticker = master_progress.add(get_ticker());
        repaired_ticker.set_message("~records repaired");
        let written_ticker =
            master_progress.add(get_ticker()).with_message("~records written");

        std::thread::spawn(move || {
            let pair_iter = ZipRecordsIter::new(
                donor_records.records(),
                acceptor_records.records(),
                donor_ticker,
                acceptor_ticker,
            );
            for pair in pair_iter {
                match pair_snd.send(pair) {
                    Ok(_) => {}
                    Err(e) => {
                        error!(
                            "failed to send record pair on channel, {}",
                            e.to_string()
                        );
                    }
                }
            }
        });

        let (repair_snd, repair_rcv) = std::sync::mpsc::sync_channel(1000);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(pool_threads)
            .build()
            .context("failed to make thread pool")?;
        std::thread::spawn(move || {
            pool.install(|| {
                pair_rcv.into_iter().par_bridge().for_each(|record_pair| {
                    let repaired = repair_record_pair(record_pair);
                    match repair_snd.send(repaired) {
                        Ok(_) => repaired_ticker.inc(1),
                        Err(e) => {
                            error!(
                                "failed to send repaired record on channel, {}",
                                e.to_string()
                            );
                        }
                    }
                })
            })
        });

        let mut n_repaired = 0usize;
        let mut n_failed = 0usize;
        for res in repair_rcv {
            match res {
                Ok(record) => {
                    if let Err(e) = writer.write(&record) {
                        error!("failed to write record {}", e.to_string());
                        n_failed += 1;
                    } else {
                        written_ticker.inc(1);
                        n_repaired += 1;
                    }
                }
                Err(e) => {
                    debug!("record failed to be repaired: {}", e.to_string());
                    n_failed += 1;
                }
            }
        }

        info!("finished, repaired {n_repaired} records, {n_failed} failed.");
        Ok(())
    }
}

#[derive(new)]
struct RecordPair {
    donor: Arc<bam::Record>,
    acceptor: bam::Record,
}

struct ZipRecordsIter<'a, T: Read> {
    donor_records: bam::Records<'a, T>,
    acceptor_records: bam::Records<'a, T>,
    cur_donor_record: Option<Arc<bam::Record>>,
    cur_acceptor_record: Option<bam::Record>,
    donor_ticker: ProgressBar,
    acceptor_ticker: ProgressBar,
}

impl<'a, T: Read> ZipRecordsIter<'a, T> {
    fn new(
        donor: bam::Records<'a, T>,
        acceptor: bam::Records<'a, T>,
        donor_ticker: ProgressBar,
        acceptor_ticker: ProgressBar,
    ) -> Self {
        Self {
            donor_records: donor,
            acceptor_records: acceptor,
            cur_donor_record: None,
            cur_acceptor_record: None,
            donor_ticker,
            acceptor_ticker,
        }
    }
}

fn get_next_record<T: Read>(
    records: &mut bam::Records<T>,
    donor: bool,
) -> Option<bam::Record> {
    loop {
        match records.next() {
            Some(Ok(record)) => {
                if record_is_secondary(&record) && donor {
                    continue;
                } else {
                    break Some(record);
                }
            }
            Some(Err(e)) => {
                let label = if donor { "donor" } else { "acceptor" };
                warn!(
                    "failed to parse record from {label} BAM, {}",
                    e.to_string()
                );
                continue;
            }
            None => break None,
        }
    }
}

impl<'a, T: Read> ZipRecordsIter<'a, T> {
    fn advance_donor_record(&mut self) {
        match self.cur_donor_record {
            Some(_) => {
                return;
            }
            None => {
                self.cur_donor_record =
                    get_next_record(&mut self.donor_records, true)
                        .map(|rec| Arc::new(rec))
            }
        }
    }
    fn advance_acceptor_record(&mut self) {
        match self.cur_acceptor_record {
            Some(_) => {
                return;
            }
            None => {
                self.cur_acceptor_record =
                    get_next_record(&mut self.acceptor_records, false);
            }
        }
    }
}

impl<'a, T: Read> Iterator for ZipRecordsIter<'a, T> {
    type Item = RecordPair;

    fn next(&mut self) -> Option<Self::Item> {
        // advance a to next record
        loop {
            self.advance_donor_record();
            self.advance_acceptor_record();

            return match (
                self.cur_donor_record.as_ref(),
                self.cur_acceptor_record.as_ref(),
            ) {
                (Some(donor), Some(acceptor)) => {
                    match donor.qname().eq(acceptor.qname()) {
                        true => {
                            // unwrap are safe because of the above match,
                            // advances acceptor on next
                            // call to .next
                            let acceptor_record = std::mem::replace(
                                &mut self.cur_acceptor_record,
                                None,
                            )
                            .unwrap();
                            self.acceptor_ticker.inc(1);
                            return Some(RecordPair::new(
                                donor.clone(),
                                acceptor_record,
                            ));
                        }
                        false => {
                            // advance donor record in attempt to find this
                            // acceptor
                            // todo consider logging?
                            let _ = std::mem::replace(
                                &mut self.cur_donor_record,
                                None,
                            );
                            self.donor_ticker.inc(1);
                            continue;
                        }
                    }
                }
                (None, Some(_)) => {
                    // no more donors, but still some acceptors.. error case
                    error!("ran out of donor records");
                    None
                }
                (Some(_), None) => {
                    debug!("exhausted acceptor BAM reader, finished.");
                    None
                }
                (None, None) => None,
            };
        }
    }
}

fn repair_record_pair(record_pair: RecordPair) -> anyhow::Result<bam::Record> {
    let read_name =
        get_query_name_string(&record_pair.donor).unwrap_or_else(|e| {
            format!("failed to parse query name, {}", e.to_string())
        });
    let modbase_info = ModBaseInfo::new_from_record(&record_pair.donor)
        .map_err(|e| anyhow!("record {read_name} failed, {}", e.to_string()))?;

    let donor_seq = get_forward_sequence(&record_pair.donor).map_err(|e| {
        anyhow!(
            "donor sequence for record {read_name} failed, {}",
            e.to_string()
        )
    })?;
    let acceptor_seq =
        get_forward_sequence(&record_pair.acceptor).map_err(|e| {
            anyhow!(
                "acceptor sequence for record {read_name} failed, {}",
                e.to_string()
            )
        })?;

    if donor_seq.len() < acceptor_seq.len() {
        bail!("donor sequence for {read_name} is longer than acceptor sequence")
    }
    let matches = donor_seq.match_indices(&acceptor_seq);

    let starts =
        matches.into_iter().map(|(start, _)| start).collect::<Vec<usize>>();
    if starts.len() > 1 {
        bail!("multiple potential corrections found for {read_name}")
    } else if starts.is_empty() {
        bail!("acceptor sequence is not a substring of the donor sequence")
    } else {
        let start = *starts.get(0).unwrap();
        let end = start + acceptor_seq.len();

        let mm_style = modbase_info.mm_style;
        let ml_style = modbase_info.ml_style;

        let mut mm_agg = String::new();
        let mut ml_agg = Vec::new();

        let (_, base_mod_probs_iter) = modbase_info.into_iter_base_mod_probs();
        for (primary_base, strand, seq_pos_base_mod_probs) in
            base_mod_probs_iter
        {
            let converter =
                DeltaListConverter::new(&acceptor_seq, primary_base);
            let skip_mode = seq_pos_base_mod_probs.skip_mode;
            let adjusted = seq_pos_base_mod_probs
                .pos_to_base_mod_probs
                .into_iter()
                .filter_map(|(pos, base_mod_probs)| {
                    if pos >= start && pos < end {
                        Some((pos - start, base_mod_probs))
                    } else {
                        None
                    }
                })
                .collect::<FxHashMap<usize, BaseModProbs>>();
            let repaired_seq_pos_base_mod_probs =
                SeqPosBaseModProbs::new(adjusted, skip_mode);
            let (mm, mut ml) = format_mm_ml_tag(
                repaired_seq_pos_base_mod_probs,
                strand,
                &converter,
            );
            mm_agg.push_str(&mm);
            ml_agg.extend_from_slice(&mut ml);
        }

        let mm = Aux::String(&mm_agg);
        let ml_arr: AuxArray<u8> = {
            let sl = &ml_agg;
            sl.into()
        };
        let ml = Aux::ArrayU8(ml_arr);

        let mut repaired_record = record_pair.acceptor;
        for tag in MM_TAGS.iter().chain(ML_TAGS.iter()) {
            let _ = repaired_record.remove_aux(tag.as_bytes());
        }
        repaired_record.push_aux(mm_style.as_bytes(), mm).map_err(|e| {
            RunError::new_failed(format!(
                "failed to add MM tag, {}",
                e.to_string()
            ))
        })?;
        repaired_record.push_aux(ml_style.as_bytes(), ml).map_err(|e| {
            RunError::new_failed(format!(
                "failed to add ML tag, {}",
                e.to_string()
            ))
        })?;
        Ok(repaired_record)
    }
}
