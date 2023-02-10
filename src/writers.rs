use crate::mod_pileup::ModBasePileup;
use crate::summarize::ModSummary;
use anyhow::Context;
use std::fs::File;
use std::io::{BufWriter, Write};

pub trait OutWriter<T> {
    fn write(&mut self, item: T) -> Result<u64, String>;
}

pub struct BedMethylWriter {
    buf_writer: BufWriter<File>,
}

impl BedMethylWriter {
    pub fn new(buf_writer: BufWriter<File>) -> Self {
        Self { buf_writer }
    }
}

impl OutWriter<ModBasePileup> for BedMethylWriter {
    fn write(&mut self, item: ModBasePileup) -> Result<u64, String> {
        let mut rows_written = 0;
        let sep = '\t';
        for (pos, feature_counts) in item.iter_counts() {
            for feature_count in feature_counts {
                let row = format!(
                    "{}{sep}{}{sep}{}{sep}{}{sep}{}{sep}{}{sep}{}{sep}{}{sep}\
                    {}{sep}{}{sep}{}{sep}{}{sep}{}{sep}{}{sep}{}{sep}{}{sep}{}\
                    {sep}{}{sep}\n",
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
                    .with_context(|| "failed to write row")
                    .map_err(|e| e.to_string())?;
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
    fn write(&mut self, item: ModSummary) -> Result<u64, String> {
        let mut report = String::new();
        let mod_called_bases = item
            .mod_called_bases
            .iter()
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
            for (mod_code, counts) in mod_counts {
                let label = if mod_code.is_canonical() {
                    format!("unmodified")
                } else {
                    format!("modified_{}", mod_code.char())
                };
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
            }
            report.push_str(&format!(
                "total_mod_calls\t{}\n",
                total_calls as u64
            ));
        }

        report.push_str(&format!(
            "total_reads_used\t{}\n",
            item.total_reads_used
        ));

        self.buf_writer
            .write(report.as_bytes())
            .map_err(|e| e.to_string())?;
        Ok(1)
    }
}
