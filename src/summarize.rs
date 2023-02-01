use crate::errs::RunError;
use crate::mod_bam::{
    base_mod_probs_from_record, get_canonical_bases_with_mod_calls,
    BaseModCall, DeltaListConverter,
};
use crate::mod_base_code::{DnaBase, ModCode};
use crate::util::record_is_secondary;
use indicatif::{ProgressBar, ProgressStyle};
use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::collections::{HashMap, HashSet};
use std::path::Path;

#[derive(Debug)]
pub struct ModSummary {
    pub mod_called_bases: Vec<DnaBase>,
    pub reads_with_mod_calls: HashMap<DnaBase, u64>,
    pub mod_call_counts: HashMap<DnaBase, HashMap<ModCode, u64>>,
    pub total_reads_used: usize,
}

pub fn summarize_modbam<T: AsRef<Path>>(
    bam_fp: T,
    threads: usize,
) -> Result<ModSummary, RunError> {
    let mut reader = bam::Reader::from_path(bam_fp)
        .map_err(|e| RunError::new_input_error(e.to_string()))?;
    reader.set_threads(threads).map_err(|e| {
        RunError::new_failed(format!(
            "failed to set threads on reader, {}",
            e.to_string()
        ))
    })?;

    let record_iter = reader
        .records()
        // skip records that fail to parse htslib
        .filter_map(|res| res.ok())
        // skip non-primary
        .filter(|record| !record_is_secondary(&record))
        // skip records with empty sequences
        .filter(|record| record.seq_len() > 0)
        // pull out the canonical bases in the MM tags, drop records that fail to parse
        .filter_map(|record| {
            get_canonical_bases_with_mod_calls(&record)
                .map(|bases| (bases, record))
                .ok()
        });

    let spinner = ProgressBar::new_spinner();
    spinner.set_style(
        ProgressStyle::with_template(
            "{spinner:.blue} [{elapsed_precise}] {pos} {msg}",
        )
        .unwrap()
        .tick_strings(&[
            "▹▹▹▹▹",
            "▸▹▹▹▹",
            "▹▸▹▹▹",
            "▹▹▸▹▹",
            "▹▹▹▸▹",
            "▹▹▹▹▸",
            "▪▪▪▪▪",
        ]),
    );
    spinner.set_message("records processed");

    let mut total_reads_used = 0;
    let mut mod_called_bases = HashSet::new();
    let mut reads_with_mod_calls = HashMap::new();
    let mut mod_call_counts = HashMap::new();
    for (i, (canonical_bases, record)) in record_iter.enumerate() {
        for canonical_base in canonical_bases {
            mod_called_bases.insert(canonical_base);
            let mod_counts = mod_call_counts
                .entry(canonical_base)
                .or_insert(HashMap::new());

            let converter = DeltaListConverter::new_from_record(
                &record,
                canonical_base.char(),
            )?;
            match base_mod_probs_from_record(
                &record,
                &converter,
                canonical_base.char(),
            ) {
                Ok(seq_pos_base_mod_probs) => {
                    let count =
                        reads_with_mod_calls.entry(canonical_base).or_insert(0);
                    *count += 1;

                    for (_position, base_mod_probs) in seq_pos_base_mod_probs {
                        match base_mod_probs.base_mod_call() {
                            BaseModCall::Canonical(_p) => {
                                let count = mod_counts
                                    .entry(
                                        canonical_base
                                            .canonical_mod_code()
                                            .unwrap(),
                                    )
                                    .or_insert(0);
                                *count += 1;
                            }
                            BaseModCall::Modified(_p, mod_code) => {
                                let count =
                                    mod_counts.entry(mod_code).or_insert(0);
                                *count += 1;
                            }
                            BaseModCall::Filtered => {}
                        }
                    }
                }
                Err(_err) => {}
            }
        }
        total_reads_used = i;
        spinner.inc(1);
    }

    let mut mod_called_bases = mod_called_bases.into_iter().collect::<Vec<_>>();
    mod_called_bases.sort();

    Ok(ModSummary {
        mod_called_bases,
        reads_with_mod_calls,
        mod_call_counts,
        total_reads_used: total_reads_used + 1,
    })
}
