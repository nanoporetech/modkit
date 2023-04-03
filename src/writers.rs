use crate::mod_pileup::ModBasePileup;
use crate::summarize::ModSummary;
use crate::util::Strand;
use anyhow::{anyhow, Context, Result as AnyhowResult};

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

pub trait OutWriter<T> {
    fn write(&mut self, item: T) -> AnyhowResult<u64>;
}

pub struct BedMethylWriter {
    buf_writer: BufWriter<File>,
    tabs_and_spaces: bool,
}

impl BedMethylWriter {
    pub fn new(buf_writer: BufWriter<File>, tabs_and_spaces: bool) -> Self {
        Self {
            buf_writer,
            tabs_and_spaces,
        }
    }
}

impl OutWriter<ModBasePileup> for BedMethylWriter {
    fn write(&mut self, item: ModBasePileup) -> AnyhowResult<u64> {
        let mut rows_written = 0;
        let tab = '\t';
        let space = if self.tabs_and_spaces { tab } else { ' ' };
        for (pos, feature_counts) in item.iter_counts() {
            for feature_count in feature_counts {
                let row = format!(
                    "{}{tab}\
                    {}{tab}\
                    {}{tab}\
                    {}{tab}\
                    {}{tab}\
                    {}{tab}\
                    {}{tab}\
                    {}{tab}\
                    {}{tab}\
                    {}{space}\
                    {}{space}\
                    {}{space}\
                    {}{space}\
                    {}{space}\
                    {}{space}\
                    {}{space}\
                    {}{space}\
                    {}\n",
                    item.chrom_name,
                    pos,
                    pos + 1,
                    feature_count.raw_mod_code,
                    feature_count.filtered_coverage,
                    feature_count.strand.to_char(),
                    pos,
                    pos + 1,
                    "255,0,0",
                    feature_count.filtered_coverage,
                    format!("{:.2}", feature_count.fraction_modified * 100f32),
                    feature_count.n_modified,
                    feature_count.n_canonical,
                    feature_count.n_other_modified,
                    feature_count.n_delete,
                    feature_count.n_filtered,
                    feature_count.n_diff,
                    feature_count.n_nocall,
                );
                self.buf_writer
                    .write(row.as_bytes())
                    .with_context(|| "failed to write row")?;
                rows_written += 1;
            }
        }
        Ok(rows_written)
    }
}

pub struct BedGraphWriter {
    out_dir: PathBuf,
    router: HashMap<(char, Strand), BufWriter<File>>,
}

impl BedGraphWriter {
    pub fn new(out_dir: PathBuf) -> AnyhowResult<Self> {
        if out_dir.is_file() {
            Err(anyhow!("out dir cannot be a file, needs to be a directory"))
        } else {
            if !out_dir.exists() {
                std::fs::create_dir_all(out_dir.clone())?;
            }
            Ok(Self {
                out_dir,
                router: HashMap::new(),
            })
        }
    }

    fn get_writer_for_modstrand(
        &mut self,
        strand: Strand,
        raw_mod_code: char,
    ) -> &mut BufWriter<File> {
        self.router
            .entry((raw_mod_code, strand))
            .or_insert_with(|| {
                let strand_label = match strand {
                    Strand::Positive => "positive",
                    Strand::Negative => "negative",
                };
                let fp = self.out_dir.join(format!(
                    "{}_{}.bedgraph",
                    raw_mod_code, strand_label
                ));
                let fh = File::create(fp).unwrap();
                BufWriter::new(fh)
            })
    }
}

impl OutWriter<ModBasePileup> for BedGraphWriter {
    fn write(&mut self, item: ModBasePileup) -> AnyhowResult<u64> {
        let mut rows_written = 0;
        let tab = '\t';
        for (pos, feature_counts) in item.iter_counts() {
            for feature_count in feature_counts {
                let fh = self.get_writer_for_modstrand(
                    feature_count.strand,
                    feature_count.raw_mod_code,
                );
                let row = format!(
                    "{}{tab}\
                     {}{tab}\
                     {}{tab}\
                     {}\n",
                    item.chrom_name,
                    pos,
                    pos + 1,
                    feature_count.fraction_modified
                );
                fh.write(row.as_bytes()).unwrap();
                rows_written += 1;
            }
        }

        Ok(rows_written)
    }
}

pub struct TsvWriter<W: Write> {
    buf_writer: BufWriter<W>,
}

impl TsvWriter<std::io::Stdout> {
    pub fn new() -> Self {
        let out = BufWriter::new(std::io::stdout());

        Self { buf_writer: out }
    }
}

impl<W: Write> OutWriter<ModSummary> for TsvWriter<W> {
    fn write(&mut self, item: ModSummary) -> AnyhowResult<u64> {
        let mut report = String::new();
        let mod_called_bases = item
            .mod_call_counts
            .keys()
            .map(|d| d.char().to_string())
            .collect::<Vec<String>>()
            .join(",");
        report.push_str(&format!("mod_bases\t{}\n", mod_called_bases));
        for (dna_base, read_count) in item.reads_with_mod_calls {
            report.push_str(&format!(
                "count_reads_{}\t{}\n",
                dna_base.char(),
                read_count
            ));
        }
        for (canonical_base, mod_counts) in item.mod_call_counts {
            let total_calls = mod_counts.values().sum::<u64>() as f64;
            let total_filtered_calls = item
                .filtered_mod_calls
                .get(&canonical_base)
                .map(|filtered_counts| filtered_counts.values().sum::<u64>())
                .unwrap_or(0);
            for (mod_code, counts) in mod_counts {
                let label = if mod_code.is_canonical() {
                    format!("unmodified")
                } else {
                    format!("modified_{}", mod_code.char())
                };
                let filtered = *item
                    .filtered_mod_calls
                    .get(&canonical_base)
                    .and_then(|filtered_counts| filtered_counts.get(&mod_code))
                    .unwrap_or(&0);
                report.push_str(&format!(
                    "{}_calls_{}\t{}\n",
                    canonical_base.char(),
                    label,
                    counts
                ));
                report.push_str(&format!(
                    "{}_frac_{}\t{}\n",
                    canonical_base.char(),
                    label,
                    counts as f64 / total_calls
                ));
                report.push_str(&format!(
                    "{}_filtered_{}\t{}\n",
                    canonical_base.char(),
                    label,
                    filtered
                ));
            }
            report.push_str(&format!(
                "{}_total_mod_calls\t{}\n",
                canonical_base.char(),
                total_calls as u64
            ));
            report.push_str(&format!(
                "{}_total_filtered_mod_calls\t{}\n",
                canonical_base.char(),
                total_filtered_calls
            ));
        }

        report.push_str(&format!(
            "total_reads_used\t{}\n",
            item.total_reads_used
        ));

        self.buf_writer.write(report.as_bytes())?;
        Ok(1)
    }
}
