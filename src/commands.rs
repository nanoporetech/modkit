use crate::errs::{InputError, RunError};
use crate::mod_bam::{
    base_mod_probs_from_record, collapse_mod_probs, format_mm_ml_tag, DeltaListConverter,
};
use std::io::BufWriter;
// use crate::mod_base_code::ModificationMotif;
use crate::interval_chunks::IntervalChunks;
use crate::mod_pileup::{process_region, ModBasePileup};
use crate::writers::{BEDWriter, OutWriter};
use clap::{Args, Subcommand};
use crossbeam_channel::bounded;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use rust_htslib::bam;
use rust_htslib::bam::record::{Aux, AuxArray};
use rust_htslib::bam::Read;
use std::path::PathBuf;
use std::thread;

#[derive(Subcommand)]
pub enum Commands {
    /// Collapse N-way base modification calls to (N-1)-way
    Collapse(Collapse),
    /// Pileup (combine) mod calls across genomic positions.
    Pileup(ModBamPileup),
}

impl Commands {
    pub fn run(&self) -> Result<(), String> {
        match self {
            Self::Collapse(x) => x.run(),
            Self::Pileup(x) => x.run(),
        }
    }
}

#[derive(Args)]
pub struct Collapse {
    /// BAM file to collapse mod call from
    in_bam: PathBuf,
    /// File path to new BAM file
    out_bam: PathBuf,
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

pub(crate) fn get_spinner() -> ProgressBar {
    let spinner = ProgressBar::new_spinner();
    spinner.set_style(
        ProgressStyle::with_template("{spinner:.blue} {msg} [{elapsed_precise}] {pos}")
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

type CliResult<T> = Result<T, RunError>;

fn flatten_mod_probs(
    mut record: bam::Record,
    canonical_base: char,
    mod_base_to_remove: char,
) -> CliResult<bam::Record> {
    if record.is_supplementary() || record.is_secondary() || record.is_duplicate() {
        return Err(RunError::new_skipped("not primary"));
    }
    if record.seq_len() == 0 {
        return Err(RunError::new_failed("seq is zero length"));
    }

    let converter = DeltaListConverter::new_from_record(&record, canonical_base)?;
    let probs_for_positions = base_mod_probs_from_record(&record, &converter, canonical_base)?;
    let collapsed_probs_for_positions = collapse_mod_probs(probs_for_positions, mod_base_to_remove);
    let (mm, ml) = format_mm_ml_tag(collapsed_probs_for_positions, canonical_base, &converter);

    record
        .remove_aux("MM".as_bytes())
        .expect("failed to remove MM tag");
    record
        .remove_aux("ML".as_bytes())
        .expect("failed to remove ML tag");
    let mm = Aux::String(&mm);
    let ml_arr: AuxArray<u8> = {
        let sl = &ml;
        sl.into()
    };
    let ml = Aux::ArrayU8(ml_arr);
    record
        .push_aux("MM".as_bytes(), mm)
        .expect("failed to add MM tag");
    record
        .push_aux("ML".as_bytes(), ml)
        .expect("failed to add ML tag");
    Ok(record)
}

impl Collapse {
    pub fn run(&self) -> Result<(), String> {
        let fp = &self.in_bam;
        let out_fp = &self.out_bam;
        let threads = self.threads;
        let canonical_base = self.base;
        let mod_base_to_remove = self.mod_base;
        let fail_fast = self.fail_fast;

        let mut reader = bam::Reader::from_path(fp).map_err(|e| e.to_string())?;
        reader.set_threads(threads).map_err(|e| e.to_string())?;

        let header = bam::Header::from_template(reader.header());
        let mut out_bam =
            bam::Writer::from_path(out_fp, &header, bam::Format::Bam).map_err(|e| e.to_string())?;

        let spinner = get_spinner();
        let message = format!(
            "Removing mod base {} from {}, new bam {}",
            mod_base_to_remove,
            fp.to_str().unwrap_or("???"),
            out_fp.to_str().unwrap_or("???")
        );
        eprintln!("> {}", message);
        spinner.set_message("Flattening ModBAM");
        let mut total = 0usize;
        let mut total_failed = 0usize;
        let mut total_skipped = 0usize;
        for (i, result) in reader.records().enumerate() {
            if let Ok(record) = result {
                match flatten_mod_probs(record, canonical_base, mod_base_to_remove) {
                    Err(RunError::BadInput(InputError(err))) | Err(RunError::Failed(err)) => {
                        if fail_fast {
                            return Err(err.to_string());
                        } else {
                            total_failed += 1;
                        }
                    }
                    Err(RunError::Skipped(_reason)) => {
                        total_skipped += 1;
                    }
                    Ok(record) => {
                        if let Err(err) = out_bam.write(&record) {
                            if fail_fast {
                                return Err(format!("failed to write {}", err.to_string()));
                            } else {
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
                    return Err(err);
                }
                total_failed += 1;
            }
        }
        spinner.finish_and_clear();

        eprintln!(
            "> done, {} records processed, {} failed, {} skipped",
            total, total_failed, total_skipped
        );
        Ok(())
    }
}

const ALLOWED_MOD_CODES: [char; 4] = ['h', 'm', 'a', 'c'];
fn check_raw_modbase_code(raw_code: &str) -> Result<String, String> {
    for raw_modbase_code in raw_code.chars() {
        if !ALLOWED_MOD_CODES.contains(&raw_modbase_code) {
            return Err(format!(
                "mod base code {raw_modbase_code} not allowed, options are {:?}",
                ALLOWED_MOD_CODES
            ));
        }
    }
    return Ok(raw_code.to_string());
}

#[derive(Args)]
pub struct ModBamPileup {
    /// Input BAM, should be sorted and have associated index
    in_bam: PathBuf,
    /// Output file
    out_bed: PathBuf,
    /// TODO, unused atm
    #[arg(
        group="mod-args",
        short='c',
        long,
        default_value_t=String::from("hm"),
        value_parser = check_raw_modbase_code)
    ]
    modbase_code: String,
}

impl ModBamPileup {
    fn run(&self) -> Result<(), String> {
        let header = bam::Reader::from_path(&self.in_bam)
            .map_err(|e| e.to_string())
            .map(|reader| reader.header().to_owned())?;
        let tids = (0..header.target_count())
            .map(|tid| {
                let size = header.target_len(tid).unwrap() as u32;
                (tid, size)
            })
            .collect::<Vec<(u32, u32)>>();
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap();

        let (snd, rx) = bounded(0);
        let in_bam_fp = self.in_bam.clone();
        thread::spawn(move || {
            pool.install(|| {
                for (tid, size) in tids {
                    let intervals = IntervalChunks::new(size, 50, 0).collect::<Vec<(u32, u32)>>();
                    let mut result: Vec<Result<ModBasePileup, String>> = vec![];
                    let (res, _) = rayon::join(
                        || {
                            intervals
                                .into_par_iter()
                                .map(|(start, end)| process_region(&in_bam_fp, tid, start, end))
                                .collect::<Vec<Result<ModBasePileup, String>>>()
                        },
                        || {
                            result
                                .into_iter()
                                .for_each(|mod_base_pileup| snd.send(mod_base_pileup).unwrap());
                        },
                    );
                    result = res;
                    result
                        .into_iter()
                        .for_each(|pileup| snd.send(pileup).unwrap());
                }
            });
        });

        let out_fp_str = self.out_bed.clone();
        let out_fp = std::fs::File::create(out_fp_str).unwrap();
        let mut writer = BEDWriter::new(BufWriter::new(out_fp));
        let spinner = get_spinner();
        for result in rx.into_iter() {
            match result {
                Ok(mod_base_pileup) => {
                    let rows_written = writer.write(mod_base_pileup)?;
                    spinner.inc(rows_written);
                }
                Err(message) => {
                    eprintln!("> unexpected error {message}");
                }
            }
        }
        Ok(())
    }
}
