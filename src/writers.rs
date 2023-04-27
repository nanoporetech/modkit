use crate::mod_pileup::ModBasePileup;
use crate::summarize::ModSummary;
use anyhow::{anyhow, Context, Result as AnyhowResult};

use log::warn;
use prettytable::format::FormatBuilder;
use prettytable::{row, Table};
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

pub struct TableWriter<W: Write> {
    writer: BufWriter<W>,
}

impl TableWriter<std::io::Stdout> {
    pub fn new() -> Self {
        let out = BufWriter::new(std::io::stdout());
        Self { writer: out }
    }
}

impl<'a, W: Write> OutWriter<ModSummary<'a>> for TableWriter<W> {
    fn write(&mut self, item: ModSummary<'a>) -> AnyhowResult<u64> {
        let mut metadata_table = Table::new();
        let metadata_format =
            FormatBuilder::new().padding(1, 1).left_border('#').build();
        metadata_table.set_format(metadata_format);
        metadata_table.add_row(row!["bases", item.mod_bases()]);
        metadata_table.add_row(row!["total_reads_used", item.total_reads_used]);
        for (dna_base, reads_with_calls) in item.reads_with_mod_calls {
            metadata_table.add_row(row![
                format!("count_reads_{}", dna_base.char()),
                reads_with_calls
            ]);
        }
        for (dna_base, threshold) in item.per_base_thresholds {
            metadata_table.add_row(row![
                format!("filter_threshold_{}", dna_base.char()),
                threshold
            ]);
        }
        if let Some(region) = item.region {
            metadata_table.add_row(row!["region", region.to_string()]);
        }
        let emitted = metadata_table.print(&mut self.writer)?;

        let mut report_table = Table::new();
        report_table.set_format(*prettytable::format::consts::FORMAT_CLEAN);
        report_table.set_titles(row![
            "base",
            "code",
            "all_count",
            "all_frac",
            "pass_count",
            "pass_frac"
        ]);

        for (canonical_base, pass_mod_to_counts) in item.mod_call_counts {
            let total_pass_calls = pass_mod_to_counts.values().sum::<u64>();
            let total_filtered_calls = item
                .filtered_mod_call_counts
                .get(&canonical_base)
                .map(|filtered_counts| filtered_counts.values().sum::<u64>())
                .unwrap_or(0);
            let total_calls = total_filtered_calls + total_pass_calls;

            for (mod_code, pass_counts) in pass_mod_to_counts {
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
                let all_counts = pass_counts + filtered;
                let all_frac = all_counts as f32 / total_calls as f32;
                let pass_frac = pass_counts as f32 / total_pass_calls as f32;
                report_table.add_row(row![
                    canonical_base.char(),
                    label,
                    all_counts,
                    all_frac,
                    pass_counts,
                    pass_frac
                ]);
            }
        }
        let mut report_emitted = report_table.print(&mut self.writer)?;
        report_emitted += emitted;
        Ok(report_emitted as u64)
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
        warn!("this output format will not be default in the next version, the table output \
            (set with --table) will become default and this format will require the --tsv option"
        );
        let mut report = String::new();
        let mod_called_bases = item.mod_bases();
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
                .filtered_mod_call_counts
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
                    .filtered_mod_call_counts
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
