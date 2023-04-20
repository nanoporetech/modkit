use crate::mod_pileup::ModBasePileup;
use crate::summarize::ModSummary;
use anyhow::{anyhow, Context, Result as AnyhowResult};

use log::debug;
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
                    feature_count.raw_strand,
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
    prefix: Option<String>,
    out_dir: PathBuf,
    router: HashMap<(char, char), BufWriter<File>>,
}

impl BedGraphWriter {
    pub fn new(
        out_dir: PathBuf,
        prefix: Option<&String>,
    ) -> AnyhowResult<Self> {
        if out_dir.is_file() {
            Err(anyhow!("out dir cannot be a file, needs to be a directory"))
        } else {
            if !out_dir.exists() {
                std::fs::create_dir_all(out_dir.clone())?;
            }
            Ok(Self {
                prefix: prefix.map(|s| s.to_owned()),
                out_dir,
                router: HashMap::new(),
            })
        }
    }

    fn get_writer_for_modstrand(
        &mut self,
        strand: char,
        raw_mod_code: char,
    ) -> &mut BufWriter<File> {
        self.router
            .entry((raw_mod_code, strand))
            .or_insert_with(|| {
                let strand_label = match strand {
                    '+' => "positive",
                    '-' => "negative",
                    '.' => "combined",
                    _ => "_unknown",
                };
                let filename = if let Some(p) = &self.prefix {
                    format!("{}_{}_{}.bedgraph", p, raw_mod_code, strand_label)
                } else {
                    format!("{}_{}.bedgraph", raw_mod_code, strand_label)
                };
                let fp = self.out_dir.join(filename);
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
                    feature_count.raw_strand,
                    feature_count.raw_mod_code,
                );
                let row = format!(
                    "{}{tab}\
                     {}{tab}\
                     {}{tab}\
                     {}{tab}\
                     {}\n",
                    item.chrom_name,
                    pos,
                    pos + 1,
                    feature_count.fraction_modified,
                    feature_count.filtered_coverage,
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

impl<'a, W: Write> OutWriter<ModSummary<'a>> for TsvWriter<W> {
    fn write(&mut self, item: ModSummary) -> AnyhowResult<u64> {
        let tab = '\t';
        let newline = '\n';

        let mut header = String::new();
        let mod_called_bases = item
            .mod_call_counts
            .keys()
            .map(|d| d.char().to_string())
            .collect::<Vec<String>>()
            .join(",");
        header.push_str(&format!(
            "#   total_reads_used{tab}{}{newline}",
            item.total_reads_used
        ));
        header.push_str(&format!(
            "#              bases{tab}{}{newline}",
            mod_called_bases
        ));
        for (dna_base, reads_with_calls) in item.reads_with_mod_calls {
            header.push_str(&format!(
                "#      count_reads_{}{tab}{}{newline}",
                dna_base.char(),
                reads_with_calls
            ));
        }
        for (dna_base, threshold) in item.per_base_thresholds {
            header.push_str(&format!(
                "# filter_threshold_{}{tab}{}{newline}",
                dna_base.char(),
                threshold
            ));
        }

        if let Some(region) = item.region {
            header.push_str(&format!(
                "#             region{tab}{}{newline}",
                region.to_string()
            ));
        }

        let mut report = String::new();
        report.push_str(&format!(
            "base{tab}code{tab}count{tab}frac{tab}filt_count{tab}filt_frac{newline}"
        ));
        for (canonical_base, mod_counts) in item.mod_call_counts {
            // total calls here are filtered counts, (i.e. after filtering)
            let total_calls = mod_counts.values().sum::<u64>() as f64;
            let total_filtered_calls = item
                .filtered_mod_call_counts
                .get(&canonical_base)
                .map(|filtered_counts| filtered_counts.values().sum::<u64>())
                .unwrap_or(0);

            // counts here are _filtered_ counts
            for (mod_code, counts) in mod_counts {
                let label = if mod_code.is_canonical() {
                    format!("-")
                } else {
                    format!("{}", mod_code.char())
                };
                let filtered = *item
                    .filtered_mod_call_counts
                    .get(&canonical_base)
                    .and_then(|filtered_counts| filtered_counts.get(&mod_code))
                    .unwrap_or(&0);
                let call_frac = counts as f32 / total_calls as f32;
                let filt_frac = filtered as f32 / total_filtered_calls as f32;
                report.push_str(&format!("{}{tab}{label}{tab}{counts}{tab}{call_frac}{tab}{filtered}{tab}{filt_frac}{newline}", canonical_base.char()));
            }
        }
        self.buf_writer.write(header.as_bytes())?;
        self.buf_writer.write(report.as_bytes())?;
        Ok(1)
    }
}
