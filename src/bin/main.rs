use std::error::Error;
use std::fmt::Formatter;
use std::ptr::write;

use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use rust_htslib::bam;
use rust_htslib::bam::record::{Aux, AuxArray};
use rust_htslib::bam::Read;

use mod_flatten::errs::InputError;
use mod_flatten::mod_bam::{
    collapse_mod_probs, extract_mod_probs, format_mm_ml_tag, get_mod_probs_for_query_positions,
    DeltaListConverter,
};

const MM_TAGS: [&str; 2] = ["MM", "Mm"];
const ML_TAGS: [&str; 2] = ["ML", "Ml"];

#[derive(Parser)]
#[command(about = "flattens multi-way base-modification calls into one", long_about = None)]
struct Cli {
    in_bam: String,
    out_bam: String,
    #[arg(
        short,
        long,
        help = "canonical base to flatten calls for",
        default_value_t = 'C'
    )]
    base: char,
    #[arg(
        short,
        long,
        help = "mod base code to flatten/remove",
        default_value_t = 'h'
    )]
    mod_base: char,
    #[arg(short, long, help = "number of threads to use", default_value_t = 1)]
    threads: usize,

    #[arg(
        short,
        long = "ff",
        help = "exit on bad reads, otherwise continue",
        default_value_t = false
    )]
    fail_fast: bool,
}

type CliResult<T> = Result<T, InputError>;

fn get_mm_tag(mm_aux: &Aux, tag_key: &str) -> Result<String, InputError> {
    match mm_aux {
        Aux::String(s) => Ok(s.to_string()),
        _ => Err(InputError::new(&format!(
            "incorrect {} tag, should be string",
            tag_key
        ))),
    }
}

fn get_ml_tag(ml_aux: &Aux, tag_key: &str) -> CliResult<Vec<u16>> {
    match ml_aux {
        Aux::ArrayU8(arr) => Ok(arr.iter().map(|x| x as u16).collect()),
        _ => Err(InputError::new(&format!(
            "invalid {} tag, expected array",
            tag_key
        ))),
    }
}

/// tag keys should be the new then old tags, for example ["MM", "Mm"].
fn get_tag<T>(
    record: &bam::Record,
    tag_keys: &[&str; 2],
    parser: &dyn Fn(&Aux, &str) -> CliResult<T>,
) -> CliResult<T> {
    let tag_new = record.aux(tag_keys[0].as_bytes());
    let tag_old = record.aux(tag_keys[1].as_bytes());

    match (tag_new, tag_old) {
        (Ok(aux), _) => parser(&aux, tag_keys[0]),
        (Err(_), Ok(aux)) => parser(&aux, tag_keys[1]),
        _ => Err(InputError::new("didn't find MM or Ml tag")),
    }
}

pub(crate) fn get_spinner() -> ProgressBar {
    let spinner = ProgressBar::new_spinner();
    spinner.set_style(
        ProgressStyle::with_template(
            "{spinner:.blue} {msg} [{elapsed_precise}] {pos} reads processed",
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
    spinner
}

fn main() -> Result<(), String> {
    let cli = Cli::parse();

    let fp: &str = &cli.in_bam;
    let out_fp: &str = &cli.out_bam;
    let threads: usize = cli.threads;
    let canonical_base: char = cli.base;
    let mod_base_to_remove: char = cli.mod_base;
    let fail_fast: bool = cli.fail_fast;

    let mut reader = bam::Reader::from_path(fp).map_err(|e| e.to_string())?;
    reader.set_threads(threads).map_err(|e| e.to_string())?;

    let header = bam::Header::from_template(reader.header());
    let mut out_bam =
        bam::Writer::from_path(out_fp, &header, bam::Format::Bam).map_err(|e| e.to_string())?;

    let spinner = get_spinner();
    let message = format!(
        "Removing mod base {} from {}, new bam {}",
        mod_base_to_remove, fp, out_fp
    );
    println!("> {}", message);
    spinner.set_message("Processing BAM");
    let mut total = 0usize;
    let mut total_failed = 0usize;
    for (i, result) in reader.records().enumerate() {
        if let Ok(mut record) = result {
            let raw_seq = if record.is_reverse() {
                bio::alphabets::dna::revcomp(record.seq().as_bytes())
            } else {
                record.seq().as_bytes()
            };
            let seq = String::from_utf8(raw_seq).map_err(|e| e.to_string())?;

            let converter = DeltaListConverter::new(&seq, canonical_base);

            let tags = {
                let mm = get_tag::<String>(&record, &MM_TAGS, &get_mm_tag);
                let ml = get_tag::<Vec<u16>>(&record, &ML_TAGS, &get_ml_tag);
                match (mm, ml) {
                    (Ok(mm), Ok(ml)) => Ok((mm, ml)),
                    _ => Err("no mod tags".to_string()),
                }
            };
            if let Err(e) = tags {
                if fail_fast {
                    return Err(e);
                } else {
                    total_failed += 1;
                    continue;
                }
            }

            let (mm, ml) = tags.unwrap();
            let probs_for_positions = extract_mod_probs(&mm, &ml, canonical_base, &converter)?;
            let collapsed_probs_for_positions =
                collapse_mod_probs(probs_for_positions, mod_base_to_remove);
            let (mm, ml) = format_mm_ml_tag(collapsed_probs_for_positions, 'C', &converter);

            record
                .remove_aux("MM".as_bytes())
                .expect("failed to remove MM tag");
            record
                .remove_aux("ML".as_bytes())
                .expect("failed to remove ML tag");
            let mm = Aux::String(&mm);
            let ml_arr: AuxArray<u16> = {
                let sl = &ml;
                sl.into()
            };
            let ml = Aux::ArrayU16(ml_arr);
            record
                .push_aux("MM".as_bytes(), mm)
                .expect("failed to add MM tag");
            record
                .push_aux("ML".as_bytes(), ml)
                .expect("failed to add ML tag");

            if let Err(err) = out_bam.write(&record) {
                if fail_fast {
                    return Err(err.to_string());
                } else {
                    total_failed += 1;
                }
            } else {
                spinner.inc(1);
                total = i;
            }
        } else {
            if fail_fast {
                let err = result.err().unwrap().to_string();
                return Err(err);
            }
            total_failed += 1;
        }
    }
    spinner.finish_and_clear();

    println!(
        "> done, {} records processed, {} failed",
        total + 1,
        total_failed
    );
    Ok(())
}
