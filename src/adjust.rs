use std::collections::HashSet;

use anyhow::anyhow;
use log::{debug, info};
use rust_htslib::bam::record::{Aux, AuxArray};
use rust_htslib::bam::{self, Read};

use crate::errs::{InputError, RunError};
use crate::mod_bam::{
    collapse_mod_probs, format_mm_ml_tag, CollapseMethod, EdgeFilter,
    ModBaseInfo,
};
use crate::mod_base_code::DnaBase;
use crate::threshold_mod_caller::MultipleThresholdModCaller;
use crate::util::{
    get_forward_sequence, get_query_name_string, get_spinner,
    record_is_secondary,
};

pub fn record_is_valid(record: &bam::Record) -> Result<(), RunError> {
    if record_is_secondary(&record) {
        return Err(RunError::new_skipped("not primary"));
    }
    if record.seq_len() == 0 {
        return Err(RunError::new_failed("seq is zero length"));
    }
    Ok(())
}

pub fn adjust_mod_probs(
    mut record: bam::Record,
    methods: &[CollapseMethod],
    caller: Option<&MultipleThresholdModCaller>,
    edge_filter: Option<&EdgeFilter>,
) -> Result<bam::Record, RunError> {
    let _ok = record_is_valid(&record)?;

    let mod_base_info = ModBaseInfo::new_from_record(&record)?;
    let mm_style = mod_base_info.mm_style;
    let ml_style = mod_base_info.ml_style;

    let mut mm_agg = String::new();
    let mut ml_agg = Vec::new();

    let record_name = get_query_name_string(&record)
        .unwrap_or("FAILED-UTF8-DECODE".to_string());
    let (converters, mod_prob_iter) = mod_base_info.into_iter_base_mod_probs();
    for (base, strand, seq_pos_mod_probs) in mod_prob_iter {
        let converter = converters.get(&base).unwrap();
        let filtered_seq_pos_mod_probs = if let Some(edge_filter) = edge_filter
        {
            let forward_sequence = get_forward_sequence(&record)?;
            match seq_pos_mod_probs
                .edge_filter_positions(edge_filter, record.seq_len())
                .map(|mod_probs| {
                    let codes_to_remove = methods
                        .iter()
                        .flat_map(|method| method.get_codes_to_remove())
                        .collect::<HashSet<char>>();
                    mod_probs.add_implicit_mod_calls(
                        &forward_sequence,
                        base,
                        &codes_to_remove,
                        Some(edge_filter),
                    )
                }) {
                Some(x) => Some(x),
                None => {
                    debug!("all base mod positions for record {record_name} and canonical \
                        base {base} were filtered out");
                    None
                }
            }
        } else {
            Some(seq_pos_mod_probs)
        };
        if let Some(mut seq_pos_mod_probs) = filtered_seq_pos_mod_probs {
            for method in methods {
                seq_pos_mod_probs =
                    collapse_mod_probs(seq_pos_mod_probs, method);
            }
            match (caller, DnaBase::parse(base)) {
                (Some(caller), Ok(dna_base)) => {
                    seq_pos_mod_probs = caller
                        .call_seq_pos_mod_probs(&dna_base, seq_pos_mod_probs)
                        .map_err(|e| {
                            RunError::new_input_error(e.to_string())
                        })?;
                }
                (Some(_), Err(e)) => {
                    let e = e.context(format!(
                        "failed to parse DNA base, cannot use threshold."
                    ));
                    return Err(RunError::new_input_error(e.to_string()));
                }
                _ => {}
            }

            let (mm, mut ml) =
                format_mm_ml_tag(seq_pos_mod_probs, strand, converter);
            mm_agg.push_str(&mm);
            ml_agg.extend_from_slice(&mut ml);
        } else {
            continue;
        }
    }

    record.remove_aux(mm_style.as_bytes()).map_err(|e| {
        RunError::new_failed(format!(
            "failed to remove MM tag, {}",
            e.to_string()
        ))
    })?;
    record.remove_aux(ml_style.as_bytes()).map_err(|e| {
        RunError::new_failed(format!(
            "failed to remove ML tag, {}",
            e.to_string()
        ))
    })?;

    let mm = Aux::String(&mm_agg);
    let ml_arr: AuxArray<u8> = {
        let sl = &ml_agg;
        sl.into()
    };
    let ml = Aux::ArrayU8(ml_arr);
    record.push_aux(mm_style.as_bytes(), mm).map_err(|e| {
        RunError::new_failed(format!("failed to add MM tag, {}", e.to_string()))
    })?;
    record.push_aux(ml_style.as_bytes(), ml).map_err(|e| {
        RunError::new_failed(format!("failed to add ML tag, {}", e.to_string()))
    })?;

    Ok(record)
}

pub fn adjust_modbam(
    reader: &mut bam::Reader,
    writer: &mut bam::Writer,
    collapse_methods: &[CollapseMethod],
    threshold_caller: Option<&MultipleThresholdModCaller>,
    edge_filter: Option<&EdgeFilter>,
    fail_fast: bool,
    verb: &'static str,
    suppress_progress: bool,
) -> anyhow::Result<()> {
    let spinner = get_spinner();
    if suppress_progress {
        spinner.set_draw_target(indicatif::ProgressDrawTarget::hidden())
    }
    spinner.set_message(verb);
    let mut total = 0usize;
    let mut total_failed = 0usize;
    let mut total_skipped = 0usize;
    for (i, result) in reader.records().enumerate() {
        if let Ok(record) = result {
            let record_name =
                get_query_name_string(&record).unwrap_or("???".to_owned());
            match adjust_mod_probs(
                record,
                &collapse_methods,
                threshold_caller,
                edge_filter,
            ) {
                Err(RunError::BadInput(InputError(err)))
                | Err(RunError::Failed(err)) => {
                    if fail_fast {
                        return Err(anyhow!("{}", err.to_string()));
                    } else {
                        debug!("read {} failed, {}", record_name, err);
                        total_failed += 1;
                    }
                }
                Err(RunError::Skipped(_reason)) => {
                    total_skipped += 1;
                }
                Ok(record) => {
                    if let Err(err) = writer.write(&record) {
                        if fail_fast {
                            return Err(anyhow!(
                                "failed to write {}",
                                err.to_string()
                            ));
                        } else {
                            debug!("failed to write {}", err);
                            total_failed += 1;
                        }
                    } else {
                        spinner.inc(1);
                        total = i;
                    }
                }
            }
        } else {
            if fail_fast {
                let err = result.err().unwrap().to_string();
                return Err(anyhow!("{}", err));
            }
            total_failed += 1;
        }
    }
    spinner.finish_and_clear();

    info!(
        "done, {} records processed, {} failed, {} skipped",
        total + 1,
        total_failed,
        total_skipped
    );
    Ok(())
}
