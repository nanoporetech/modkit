use crate::entropy::{EntropyCalculation, WindowEntropy};
use crate::errs::MkError;
use crate::util::{Strand, TAB};
use anyhow::{anyhow, bail};
use indicatif::ProgressBar;
use log::debug;
use rustc_hash::FxHashMap;
use std::collections::HashMap;
use std::fs::{File, OpenOptions};
use std::io::{stdout, BufWriter, ErrorKind, Write};
use std::ops::AddAssign;
use std::path::{Component, Path, PathBuf};

#[inline(always)]
fn write_entropy_windows<T: Write>(
    writer: &mut BufWriter<T>,
    window_entropies: &[WindowEntropy],
    chrom_id_to_name: &HashMap<u32, String>,
    drop_zeros: bool,
    write_counter: &ProgressBar,
    failure_counter: &ProgressBar,
    failure_reasons: &mut FxHashMap<String, usize>,
    verbose: bool,
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
            Some(Err(e)) => match e {
                _ => {
                    if verbose {
                        match e {
                            MkError::EntropyZeroCoverage {
                                chrom_id,
                                start,
                                end,
                            } => {
                                if let Some(chrom) =
                                    chrom_id_to_name.get(chrom_id)
                                {
                                    debug!(
                                        "{chrom}:{start}-{end}: zero coverage"
                                    );
                                } else {
                                    debug!(
                                        "{chrom_id}:{start}-{end}: zero \
                                         coverage"
                                    );
                                }
                            }
                            MkError::EntropyInsufficientCoverage {
                                chrom_id,
                                start,
                                end,
                            } => {
                                if let Some(chrom) =
                                    chrom_id_to_name.get(chrom_id)
                                {
                                    debug!(
                                        "{chrom}:{start}-{end}: zero coverage"
                                    );
                                } else {
                                    debug!(
                                        "{chrom_id}:{start}-{end}: zero \
                                         coverage"
                                    );
                                }
                            }
                            _ => {}
                        }
                    }
                    failure_counter.inc(1);
                    failure_reasons
                        .entry(e.to_string())
                        .or_insert(0usize)
                        .add_assign(1usize);
                }
            },
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
                failure_counter.inc(1);
                failure_reasons
                    .entry(e.to_string())
                    .or_insert(0usize)
                    .add_assign(1usize);
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
        failure_counter: &ProgressBar,
        failure_reasons: &mut FxHashMap<String, usize>,
    ) -> anyhow::Result<()>;
}

const WINDOWS_HEADER: &'static str = "\
        #chrom\tstart\tend\tentropy\tstrand\tnum_reads\n";

fn preflight_output_file(out_fp: &Path) -> anyhow::Result<bool> {
    match std::fs::symlink_metadata(out_fp) {
        Ok(metadata) if metadata.file_type().is_file() => Ok(true),
        Ok(_) => bail!(
            "entropy output target must be a regular file and may not be a \
             symbolic link: {}",
            out_fp.display()
        ),
        Err(error) if error.kind() == ErrorKind::NotFound => Ok(false),
        Err(error) => Err(error.into()),
    }
}

fn open_output_file_untruncated(
    out_fp: &Path,
    existed: bool,
) -> std::io::Result<File> {
    let mut options = OpenOptions::new();
    options.write(true);
    if !existed {
        options.create_new(true);
    }
    #[cfg(unix)]
    {
        use std::os::unix::fs::OpenOptionsExt;
        options.custom_flags(libc::O_NOFOLLOW);
    }
    options.open(out_fp)
}

fn open_output_file(out_fp: &Path, force: bool) -> anyhow::Result<File> {
    let existed = preflight_output_file(out_fp)?;
    if existed && !force {
        bail!(
            "entropy output already exists, use --force to overwrite it: {}",
            out_fp.display()
        )
    }
    let file = open_output_file_untruncated(out_fp, existed)?;
    if force {
        file.set_len(0)?;
    }
    Ok(file)
}

fn validate_regions_prefix(
    prefix: Option<&String>,
) -> anyhow::Result<Option<&str>> {
    let Some(prefix) = prefix else {
        return Ok(None);
    };
    let path = Path::new(prefix);
    let mut components = path.components();
    match (components.next(), components.next()) {
        (Some(Component::Normal(component)), None)
            if !component.is_empty() && path.as_os_str() == component =>
        {
            Ok(Some(prefix))
        }
        _ => bail!(
            "entropy regions prefix must be exactly one non-empty filename \
             component"
        ),
    }
}

pub(super) struct WindowsWriter<T: Write> {
    output: BufWriter<T>,
    verbose: bool,
}

impl WindowsWriter<File> {
    pub(super) fn new_file(
        out_fp: &PathBuf,
        header: bool,
        verbose: bool,
        force: bool,
    ) -> anyhow::Result<Self> {
        let mut output = BufWriter::new(open_output_file(out_fp, force)?);
        if header {
            output.write(WINDOWS_HEADER.as_bytes())?;
        }
        Ok(Self { output, verbose })
    }
}

impl WindowsWriter<std::io::Stdout> {
    pub(super) fn new_stdout(
        header: bool,
        verbose: bool,
    ) -> anyhow::Result<Self> {
        let mut output = BufWriter::new(stdout());
        if header {
            output.write(WINDOWS_HEADER.as_bytes())?;
        }
        Ok(Self { output, verbose })
    }
}

pub(super) struct RegionsWriter {
    regions_bed_out: BufWriter<File>,
    windows_bed_out: BufWriter<File>,
    verbose: bool,
}

impl RegionsWriter {
    pub(super) fn new(
        out_dir: &PathBuf,
        prefix: Option<&String>,
        header: bool,
        verbose: bool,
        force: bool,
    ) -> anyhow::Result<Self> {
        let prefix = validate_regions_prefix(prefix)?;
        if out_dir.exists() && !out_dir.is_dir() {
            bail!("regions output location must be a directory")
        }
        if !out_dir.exists() {
            std::fs::create_dir_all(out_dir)?;
        }
        debug_assert!(out_dir.exists(), "out_dir should exist now");
        let regions_fp = if let Some(p) = prefix {
            out_dir.join(format!("{p}_regions.bed"))
        } else {
            out_dir.join("regions.bed")
        };
        let windows_fp = if let Some(p) = prefix {
            out_dir.join(format!("{p}_windows.bedgraph"))
        } else {
            out_dir.join("windows.bedgraph")
        };
        let regions_existed = preflight_output_file(&regions_fp)?;
        let windows_existed = preflight_output_file(&windows_fp)?;
        if !force && (regions_existed || windows_existed) {
            bail!(
                "entropy region output already exists, use --force to \
                 overwrite it"
            )
        }
        let regions_file =
            open_output_file_untruncated(&regions_fp, regions_existed)?;
        let windows_file =
            match open_output_file_untruncated(&windows_fp, windows_existed) {
                Ok(file) => file,
                Err(error) => {
                    drop(regions_file);
                    if !regions_existed {
                        std::fs::remove_file(&regions_fp)?;
                    }
                    return Err(error.into());
                }
            };
        if force {
            regions_file.set_len(0)?;
            windows_file.set_len(0)?;
        }

        let mut regions_bed_out = BufWriter::new(regions_file);
        let mut windows_bed_out = BufWriter::new(windows_file);

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

        Ok(Self { windows_bed_out, regions_bed_out, verbose })
    }
}

impl<T: Write> EntropyWriter for WindowsWriter<T> {
    fn write(
        &mut self,
        entropy_calculation: EntropyCalculation,
        chrom_id_to_name: &HashMap<u32, String>,
        drop_zeros: bool,
        write_counter: &ProgressBar,
        failure_counter: &ProgressBar,
        failure_reasons: &mut FxHashMap<String, usize>,
    ) -> anyhow::Result<()> {
        match entropy_calculation {
            EntropyCalculation::Windows(entropy_windows) => {
                write_entropy_windows(
                    &mut self.output,
                    &entropy_windows,
                    chrom_id_to_name,
                    drop_zeros,
                    write_counter,
                    failure_counter,
                    failure_reasons,
                    self.verbose,
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
        failure_counter: &ProgressBar,
        failure_reasons: &mut FxHashMap<String, usize>,
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
                    Err(e) => {
                        failure_counter.inc(1);
                        failure_reasons
                            .entry(e.to_string())
                            .or_insert(0usize)
                            .add_assign(1usize);
                    }
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
                    Some(Err(e)) => {
                        if self.verbose {
                            debug!("{chrom}:{start}-{end}, {e}");
                        }
                        failure_counter.inc(1);
                        failure_reasons
                            .entry(e.to_string())
                            .or_insert(0usize)
                            .add_assign(1usize);
                    }
                    None => {}
                }
                write_entropy_windows(
                    &mut self.windows_bed_out,
                    &region_entropy.window_entropies,
                    chrom_id_to_name,
                    drop_zeros,
                    write_counter,
                    failure_counter,
                    failure_reasons,
                    self.verbose,
                )?;
            }
            EntropyCalculation::Windows(_) => {
                bail!("shouldn't have windows with regions")
            }
        }

        Ok(())
    }
}
