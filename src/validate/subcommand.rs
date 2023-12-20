use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use clap::Args;
use log::info;
use rust_htslib::bam;

use crate::logging::init_logging;
use crate::mod_base_code::ModCodeRepr;
use crate::util::get_ticker;

#[derive(Debug, PartialEq, Eq, Hash)]
struct ChromosomeStrand {
    chromosome: String,
    fwd_strand: bool,
}

fn parse_mod_bed_line(
    line: &str,
) -> Result<(ChromosomeStrand, ModCodeRepr, Vec<u32>), String> {
    let fields: Vec<&str> = line.split_whitespace().collect();
    if fields.len() >= 6 {
        let chromosome = fields[0].to_string();
        let strand = fields[5].to_string();
        let mod_code = fields[3].to_string();
        let start: u32 = fields[1]
            .parse()
            .map_err(|e| format!("Error parsing start: {}", e))?;
        let end: u32 = fields[2]
            .parse()
            .map_err(|e| format!("Error parsing end: {}", e))?;
        let fwd_strand = match strand.as_str() {
            "+" => true,
            "-" => false,
            _ => {
                return Err(format!("Strand must be + or -. Found: {}", strand))
            }
        };
        let positions = (start..end + 1).collect();
        Ok((
            ChromosomeStrand { chromosome, fwd_strand },
            ModCodeRepr::parse(&mod_code).expect(concat!(
                "Bad modified base code: ",
                stringify!(mod_code)
            )),
            positions,
        ))
    } else {
        Err("Invalid number of fields in BED line".to_string())
    }
}

type ChromStrandNamePositions =
    HashMap<ChromosomeStrand, HashMap<ModCodeRepr, Vec<u32>>>;

fn parse_bed_file(
    file_path: &PathBuf,
    suppress_pb: bool,
) -> Result<ChromStrandNamePositions, String> {
    info!("parsing BED at {}", file_path.to_str().unwrap_or("invalid-UTF-8"));
    let mut result = HashMap::new();
    let lines_processed = get_ticker();
    if suppress_pb {
        lines_processed
            .set_draw_target(indicatif::ProgressDrawTarget::hidden());
    }
    lines_processed.set_message("rows processed");

    if let Ok(file) = File::open(file_path) {
        let reader = BufReader::new(file);
        for line in reader.lines() {
            if let Ok((chrom_strand, mod_code, positions)) =
                parse_mod_bed_line(&line.as_ref().unwrap())
            {
                result
                    .entry(chrom_strand)
                    .or_insert_with(HashMap::new)
                    .entry(mod_code)
                    .or_insert_with(Vec::new)
                    .extend(positions);
                lines_processed.inc(1);
            } else {
                return Err(format!(
                    "Error parsing BED line: {:?}",
                    line.unwrap()
                ));
            }
        }
        if result.is_empty() {
            return Err("zero valid positions parsed from BED file".to_string());
        }
        lines_processed.finish_and_clear();
        info!("processed {} BED lines", lines_processed.position());
        return Ok(result);
    } else {
        return Err(format!("Error opening BED file: {:?}", file_path));
    }
}

#[derive(Args)]
pub struct ValidateFromModbam {
    // running args
    // convert to list of bam bed inputs
    /// Argument accepts 2 values. The first value is the BAM file path with
    /// modified base tags. The second is a bed file with ground truth
    /// reference positions. The name field in the ground truth bed file
    /// should be the short name (single letter code or ChEBI ID) for a
    /// modified base or the corresponding canonical base. This argument
    /// can be provided more than once for multiple samples.
    #[arg(
	long,
	action = clap::ArgAction::Append,
	num_args = 2,
	value_names = ["BAM", "BED"]
    )]
    bam_and_bed: Vec<PathBuf>,
    /// Number of threads to use
    #[arg(short = 't', long, default_value_t = 4)]
    threads: usize,
    /// Hide the progress bar.
    #[arg(long, default_value_t = false, hide_short_help = true)]
    suppress_progress: bool,
    /// Specify a file for debug logs to be written to, otherwise ignore them.
    /// Setting a file is recommended. (alias: log)
    #[arg(long, alias = "log")]
    log_filepath: Option<PathBuf>,
}

impl ValidateFromModbam {
    pub fn run(&self) -> anyhow::Result<()> {
        let _handle = init_logging(self.log_filepath.as_ref());

        for bam_and_bed in self.bam_and_bed.chunks(2) {
            let bam = &bam_and_bed[0];
            let bed = &bam_and_bed[1];

            match parse_bed_file(bed, self.suppress_progress) {
                Ok(mod_positions) => {
                    // Print the result for demonstration
                    for (chrom_strand, mod_positions) in &mod_positions {
                        for (mod_code, positions) in mod_positions {
                            println!(
                                "{:?} - {:?}: {:?}",
                                chrom_strand, mod_code, positions
                            );
                        }
                    }
                }
                Err(error) => eprintln!("Error: {}", error),
            }
            let mut _bam_reader =
                bam::Reader::from_path(&bam).expect("Failed to load BAM");
            // loop over records
            //for record in bam_reader.records() {
            //    process_bam_record(
            //        &record.expect("Error reading BAM record"),
            //        mod_positions,
            //    );
            //}
        }
        Ok(())
    }

    /*pub fn process_bam_record(
        record: &bam::Record,
        mod_positions: ChromStrandNamePositions,
    ) {
        //mbi = ModBaseInfo::;
        //ReadModBaseProfile;
    }*/
}
