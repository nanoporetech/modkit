use crossbeam_channel::bounded;
use indicatif::{ProgressBar, ProgressStyle};
use mod_kit::interval_chunks::IntervalChunks;
use mod_kit::mod_base_code::{CpGModificationMotif, HydroxyMethylCytosineCode, ModificationMotif};
use mod_kit::mod_pileup::{ModBasePileup, ModBasePileupProcessor};
use rayon::prelude::*;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::faidx;
use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::ops::Range;
use std::path::Path;
use std::thread;

pub(crate) fn get_spinner() -> ProgressBar {
    let spinner = ProgressBar::new_spinner();
    spinner.set_style(
        ProgressStyle::with_template(
            "{spinner:.blue} {msg} [{elapsed_precise}] {pos} positions processed",
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
    let bam_fp = "data/bc_anchored_10_reads.sorted.bam";
    // let bam_fp = "data/bc_anchored_hac_synthetic_ground_truth_5hmC.sorted.bam";
    // let bam_fp = "tests/resources/fwd_rev_modbase_records.sorted.bam";
    let fasta_fp = "tests/resources/CGI_ladder_3.6kb_ref.fa";
    let header = bam::Reader::from_path(bam_fp)
        .map_err(|e| e.to_string())
        .map(|reader| reader.header().to_owned())?;
    let tids = (0..header.target_count())
        .map(|tid| {
            let size = header.target_len(tid).unwrap() as u32;
            let ref_name = String::from_utf8(header.tid2name(tid).to_vec()).unwrap();
            (tid, size, ref_name)
        })
        .collect::<Vec<(u32, u32, String)>>();

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap();

    let motif = CpGModificationMotif;
    let code = HydroxyMethylCytosineCode;

    let mut fasta_reader = faidx::Reader::from_path(fasta_fp).unwrap();
    let reference_sequences = tids
        .iter()
        .map(|(tid, size, name)| {
            let size = *size as usize;
            // todo need a fail here
            let seq = fasta_reader.fetch_seq_string(name, 0, size).unwrap();
            (*tid, seq)
        })
        .collect::<HashMap<u32, String>>();

    let (snd, rx) = bounded(0);
    thread::spawn(move || {
        pool.install(|| {
            for (tid, size, ref_name) in tids {
                let reference_seq = reference_sequences.get(&tid).unwrap();
                let intervals = IntervalChunks::new(size, 50, motif.required_overlap())
                    .collect::<Vec<(u32, u32)>>();
                let mut result: Vec<Result<ModBasePileup, String>> = vec![];
                let (res, _) = rayon::join(
                    || {
                        intervals
                            .into_par_iter()
                            .map(|(start, end)| {
                                ModBasePileupProcessor::new(bam_fp, tid, start, end).map(
                                    |mut processor| {
                                        let ref_seq = reference_seq
                                            .char_indices()
                                            .filter_map(|(pos, nt)| {
                                                if pos >= start as usize && pos <= end as usize {
                                                    Some(nt)
                                                } else {
                                                    None
                                                }
                                            })
                                            .collect::<String>();
                                        processor.process_region(&code, &motif, &ref_name, &ref_seq)
                                    },
                                )
                            })
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

                // intervals.into_par_iter().for_each(|(start, end)| {
                //     let mut processor =
                //         ModBasePileupProcessor::new(bam_fp, fasta_fp, tid, start, end).unwrap();
                //     let pileup = processor.process_region(&code, &motif, &ref_name, &reference_seq);
                //     snd.send(pileup).unwrap();
                // });
            }
        });
    });
    // thread::spawn(move || {
    //     for (tid, size, ref_name) in tids {
    //         let reference_seq = reference_sequences.get(&tid).unwrap();
    //         let intervals = IntervalChunks::new(size, 20, motif.required_overlap())
    //             .collect::<Vec<(u32, u32)>>();
    //         intervals.into_iter().for_each(|(start, end)| {
    //             let mut processor =
    //                 ModBasePileupProcessor::new(bam_fp, fasta_fp, tid, start, end).unwrap();
    //             let pileup = processor.process_region(&code, &motif, &ref_name, &reference_seq);
    //             snd.send(pileup).unwrap();
    //         });
    //     }
    // });

    let out_fp_str = "data/test-modbam2bed.bed";
    let out_fp = std::fs::File::create(out_fp_str).unwrap();
    let mut writer = BufWriter::new(out_fp);

    // let spinner = get_spinner();
    for result in rx.into_iter() {
        match result {
            Ok(mod_base_pileup) => {
                let mut total_inc = 0u64;
                if let Some((rows, n)) = mod_base_pileup.decode(&code, '\t') {
                    writer.write(rows.as_bytes()).unwrap();
                    total_inc += n;
                }
                // spinner.inc(total_inc);
            }
            Err(message) => {
                eprintln!("> unexpected error {message}");
            }
        }
    }

    eprintln!("done, results at {}", out_fp_str);
    Ok(())
}
