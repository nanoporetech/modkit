use crate::entropy::{EntropyCalculation, WindowEntropy};
use crate::errs::RunError;
use crate::util::{Strand, TAB};
use anyhow::{anyhow, bail};
use indicatif::ProgressBar;
use log::debug;
use std::collections::HashMap;
use std::fs::File;
use std::io::{stdout, BufWriter, Write};
use std::path::PathBuf;

#[inline(always)]
fn write_entropy_windows<T: Write>(
    writer: &mut BufWriter<T>,
    window_entropies: &[WindowEntropy],
    chrom_id_to_name: &HashMap<u32, String>,
    drop_zeros: bool,
    write_counter: &ProgressBar,
    skip_counter: &ProgressBar,
    failure_counter: &ProgressBar,
) -> anyhow::Result<()> {
    for entropy in window_entropies {
        let name =
            chrom_id_to_name.get(&entropy.chrom_id).ok_or_else(|| {
                anyhow!("missing chrom name for {}", &entropy.chrom_id)
            })?;
        match entropy.pos_me_entropy.as_ref() {
            Some(Ok(pos_entropy)) => {
                if (drop_zeros && !(pos_entropy.me_entropy == 0f32))
                    || !drop_zeros
                {
                    let row = format!(
                        "{name}\t{}\t{}\t{}\t{}\t{}\n",
                        pos_entropy.interval.start,
                        pos_entropy.interval.end,
                        pos_entropy.me_entropy,
                        Strand::Positive.to_char(),
                        pos_entropy.num_reads
                    );
                    writer.write(&row.as_bytes())?;
                    write_counter.inc(1);
                }
            }
            Some(Err(e)) => {
                match e {
                    RunError::Failed(e) => {
                        // debug!("(+) window failed, {e}");
                        failure_counter.inc(1);
                    }
                    RunError::BadInput(reason) => {
                        debug!("(+) window bad input, {}", &reason.0);
                        failure_counter.inc(1);
                    }
                    RunError::Skipped(_e) => {
                        skip_counter.inc(1);
                        // debug!("window {}:{}-{} skipped, {e}",
                        // name, entropy.interval.start,
                        // entropy.interval.end);
                    }
                }
            }
            None => {}
        }

        match entropy.neg_me_entropy.as_ref() {
            Some(Ok(neg_entropy)) => {
                if (drop_zeros && !(neg_entropy.me_entropy == 0f32))
                    || !drop_zeros
                {
                    let row = format!(
                        "{name}\t{}\t{}\t{}\t{}\t{}\n",
                        neg_entropy.interval.start,
                        neg_entropy.interval.end,
                        neg_entropy.me_entropy,
                        Strand::Negative.to_char(),
                        neg_entropy.num_reads
                    );
                    writer.write(&row.as_bytes())?;
                    write_counter.inc(1);
                }
            }
            Some(Err(e)) => {
                match e {
                    RunError::Failed(e) => {
                        // debug!("(-) window failed, {e}");
                        failure_counter.inc(1);
                    }
                    RunError::BadInput(reason) => {
                        debug!("(-) window bad input, {}", &reason.0);
                        failure_counter.inc(1);
                    }
                    RunError::Skipped(_e) => {
                        skip_counter.inc(1);
                        // debug!("window {}:{}-{} skipped, {e}",
                        // name, entropy.interval.start,
                        // entropy.interval.end);
                    }
                }
            }
            None => {}
        }
    }
    Ok(())
}

pub(super) trait EntropyWriter {
    fn write(
        &mut self,
        entropy_calculation: EntropyCalculation,
        chrom_id_to_name: &HashMap<u32, String>,
        drop_zeros: bool,
        write_counter: &ProgressBar,
        skip_counter: &ProgressBar,
        failure_counter: &ProgressBar,
        // todo failure causes
    ) -> anyhow::Result<()>;
}

const WINDOWS_HEADER: &'static str = "\
        chrom\tstart\tend\tentropy\tstrand\tnum_reads\n";

pub(super) struct WindowsWriter<T: Write> {
    output: BufWriter<T>,
}

impl WindowsWriter<File> {
    pub(super) fn new_file(
        out_fp: &PathBuf,
        header: bool,
    ) -> anyhow::Result<Self> {
        let mut output = BufWriter::new(File::create(out_fp)?);
        if header {
            output.write(WINDOWS_HEADER.as_bytes())?;
        }
        Ok(Self { output })
    }
}

impl WindowsWriter<std::io::Stdout> {
    pub(super) fn new_stdout(header: bool) -> anyhow::Result<Self> {
        let mut output = BufWriter::new(stdout());
        if header {
            output.write(WINDOWS_HEADER.as_bytes())?;
        }
        Ok(Self { output })
    }
}

pub(super) struct RegionsWriter {
    regions_bed_out: BufWriter<File>,
    windows_bed_out: BufWriter<File>,
}

impl RegionsWriter {
    pub(super) fn new(
        out_dir: &PathBuf,
        prefix: Option<&String>,
        header: bool,
    ) -> anyhow::Result<Self> {
        if out_dir.is_file() {
            bail!("regions output location must be a directory")
        }
        std::fs::create_dir_all(out_dir)?;
        debug_assert!(out_dir.exists(), "out_dir should exist now");
        let mut regions_bed_out = if let Some(p) = prefix {
            let fp = out_dir.join(format!("{p}_regions.bed"));
            BufWriter::new(File::create(fp)?)
        } else {
            let fp = out_dir.join("regions.bed");
            BufWriter::new(File::create(fp)?)
        };

        let mut windows_bed_out = if let Some(p) = prefix {
            let fp = out_dir.join(format!("{p}_windows.bedgraph"));
            BufWriter::new(File::create(fp)?)
        } else {
            let fp = out_dir.join("windows.bedgraph");
            BufWriter::new(File::create(fp)?)
        };

        if header {
            windows_bed_out.write(WINDOWS_HEADER.as_bytes())?;
            regions_bed_out.write(
                &format!(
                    "\
                chrom{TAB}\
                start{TAB}\
                end{TAB}\
                region_name{TAB}\
                mean_entropy{TAB}\
                strand{TAB}\
                median_entropy{TAB}\
                min_entropy{TAB}\
                max_entropy{TAB}\
                mean_num_reads{TAB}\
                min_num_reads{TAB}\
                max_num_reads{TAB}\
                successful_window_count{TAB}\
                failed_window_count\n"
                )
                .as_bytes(),
            )?;
        }

        Ok(Self { windows_bed_out, regions_bed_out })
    }
}

impl<T: Write> EntropyWriter for WindowsWriter<T> {
    fn write(
        &mut self,
        entropy_calculation: EntropyCalculation,
        chrom_id_to_name: &HashMap<u32, String>,
        drop_zeros: bool,
        write_counter: &ProgressBar,
        skip_counter: &ProgressBar,
        failure_counter: &ProgressBar,
    ) -> anyhow::Result<()> {
        match entropy_calculation {
            EntropyCalculation::Windows(entropy_windows) => {
                write_entropy_windows(
                    &mut self.output,
                    &entropy_windows,
                    chrom_id_to_name,
                    drop_zeros,
                    write_counter,
                    skip_counter,
                    failure_counter,
                )?;
            }
            EntropyCalculation::Region(_) => bail!("shouldn't have regions"),
        }
        Ok(())
    }
}

impl EntropyWriter for RegionsWriter {
    fn write(
        &mut self,
        entropy_calculation: EntropyCalculation,
        chrom_id_to_name: &HashMap<u32, String>,
        drop_zeros: bool,
        write_counter: &ProgressBar,
        skip_counter: &ProgressBar,
        failure_counter: &ProgressBar,
    ) -> anyhow::Result<()> {
        match entropy_calculation {
            EntropyCalculation::Region(region_entropy) => {
                let chrom =
                    chrom_id_to_name.get(&region_entropy.chrom_id).expect(
                        "shouldn't have a result on a chrom without a chromId",
                    );
                let start = region_entropy.interval.start;
                let end = region_entropy.interval.end;
                let region_name = region_entropy.region_name;
                match region_entropy.pos_entropy_stats {
                    Ok(pos_entropy_stats) => {
                        let row = pos_entropy_stats.to_row(
                            &chrom,
                            start,
                            end,
                            Strand::Positive,
                            &region_name,
                        );
                        self.regions_bed_out.write(row.as_bytes())?;
                        write_counter.inc(1);
                    }
                    Err(e) => match e {
                        RunError::Failed(_e) => {
                            // debug!("(+) region failed, {e}");
                            failure_counter.inc(1);
                        }
                        RunError::BadInput(reason) => {
                            debug!("(+) region bad input?, {}", &reason.0);
                            failure_counter.inc(1);
                        }
                        RunError::Skipped(_e) => {
                            skip_counter.inc(1);
                            // debug!("window {}:{}-{} skipped, {e}", name,
                            // entropy.interval.start, entropy.interval.end);
                        }
                    },
                }
                match region_entropy.neg_entropy_stats {
                    Some(Ok(neg_entropy_stats)) => {
                        let row = neg_entropy_stats.to_row(
                            &chrom,
                            start,
                            end,
                            Strand::Negative,
                            &region_name,
                        );
                        self.regions_bed_out.write(row.as_bytes())?;
                        write_counter.inc(1);
                    }
                    Some(Err(e)) => match e {
                        RunError::Failed(_e) => {
                            // debug!("(-) region failed, {e}");
                            failure_counter.inc(1);
                        }
                        RunError::BadInput(reason) => {
                            debug!("(-) region bad input?, {}", &reason.0);
                            failure_counter.inc(1);
                        }
                        RunError::Skipped(_e) => {
                            skip_counter.inc(1);
                            // debug!("window {}:{}-{} skipped, {e}", name,
                            // entropy.interval.start, entropy.interval.end);
                        }
                    },
                    None => {}
                }
                write_entropy_windows(
                    &mut self.windows_bed_out,
                    &region_entropy.window_entropies,
                    chrom_id_to_name,
                    drop_zeros,
                    write_counter,
                    skip_counter,
                    failure_counter,
                )?;
            }
            EntropyCalculation::Windows(_) => {
                bail!("shouldn't have windows with regions")
            }
        }

        Ok(())
    }
}
