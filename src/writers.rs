use crate::pileup::{ModBasePileup, PartitionKey, PileupFeatureCounts};
use crate::summarize::ModSummary;
use anyhow::{anyhow, Context, Result as AnyhowResult};

use crate::read_ids_to_base_mod_probs::ReadsBaseModProfile;
use crate::thresholds::Percentiles;
use derive_new::new;
use histo_fp::Histogram;
use log::{debug, info, warn};
use prettytable::format::FormatBuilder;
use prettytable::{cell, row, Table};
use rustc_hash::FxHashMap;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Stdout, Write};
use std::path::{Path, PathBuf};

pub trait OutwriterWithMemory<T> {
    fn write(&mut self, item: T) -> AnyhowResult<u64>;
    fn num_reads(&self) -> usize;
}

pub trait OutWriter<T> {
    fn write(&mut self, item: T) -> AnyhowResult<u64>;
}

pub struct BedMethylWriter<T: Write> {
    buf_writer: BufWriter<T>,
    tabs_and_spaces: bool,
}

impl<T: Write + Sized> BedMethylWriter<T> {
    pub fn new(buf_writer: BufWriter<T>, tabs_and_spaces: bool) -> Self {
        Self {
            buf_writer,
            tabs_and_spaces,
        }
    }

    #[inline]
    fn write_feature_counts(
        pos: u32,
        chrom_name: &str,
        feature_counts: &[PileupFeatureCounts],
        writer: &mut BufWriter<T>,
        tabs_and_spaces: bool,
    ) -> AnyhowResult<u64> {
        let tab = '\t';
        let space = if tabs_and_spaces { ' ' } else { tab };
        let mut rows_written = 0u64;
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
                chrom_name,
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
            writer
                .write(row.as_bytes())
                .with_context(|| "failed to write row")?;
            rows_written += 1;
        }

        Ok(rows_written)
    }
}

impl<T: Write> OutWriter<ModBasePileup> for BedMethylWriter<T> {
    fn write(&mut self, item: ModBasePileup) -> AnyhowResult<u64> {
        let mut rows_written = 0;
        // let tab = '\t';
        // let space = if self.tabs_and_spaces { tab } else { ' ' };
        for (pos, feature_counts) in item.iter_counts_sorted() {
            match feature_counts.get(&PartitionKey::NoKey) {
                Some(feature_counts) => {
                    rows_written += BedMethylWriter::write_feature_counts(
                        *pos,
                        &item.chrom_name,
                        &feature_counts,
                        &mut self.buf_writer,
                        self.tabs_and_spaces,
                    )?;
                }
                None => {}
            }
        }
        Ok(rows_written)
    }
}

#[derive(new, Hash, Eq, PartialEq, Copy, Clone)]
struct BedGraphFileKey {
    partition_key: PartitionKey,
    strand: char,
    raw_mode_code: char,
}

pub struct BedGraphWriter {
    prefix: Option<String>,
    out_dir: PathBuf,
    router: HashMap<BedGraphFileKey, BufWriter<File>>,
    use_groupings: bool,
}

impl BedGraphWriter {
    pub fn new(
        out_dir: &str,
        prefix: Option<&String>,
        use_groupings: bool,
    ) -> AnyhowResult<Self> {
        let out_dir_fp = Path::new(out_dir).to_path_buf();
        if !out_dir_fp.exists() {
            info!("creating directory for bedgraph output at {out_dir}");
            std::fs::create_dir_all(out_dir_fp.clone())?;
        }
        Ok(Self {
            prefix: prefix.map(|s| s.to_owned()),
            out_dir: out_dir_fp,
            router: HashMap::new(),
            use_groupings,
        })
    }

    fn get_writer_for_modstrand(
        &mut self,
        key: BedGraphFileKey,
        key_name: &str,
    ) -> &mut BufWriter<File> {
        self.router
            .entry(key)
            .or_insert_with(|| {
                let strand = key.strand;
                let raw_mod_code = key.raw_mode_code;
                let delim = if key_name == "" {
                    ""
                } else {
                    "_"
                };
                let strand_label = match strand {
                    '+' => "positive",
                    '-' => "negative",
                    '.' => "combined",
                    _ => "_unknown",
                };
                let filename = if let Some(p) = &self.prefix {
                    format!("{p}{delim}{key_name}{delim}{raw_mod_code}_{strand_label}.bedgraph")
                } else {
                    format!("{key_name}{delim}{raw_mod_code}_{strand_label}.bedgraph")
                };
                let fp = self.out_dir.join(filename);
                // todo(arand) danger, should remove this unwrap
                let fh = File::create(fp).unwrap();
                BufWriter::new(fh)
            })
    }
}

impl OutWriter<ModBasePileup> for BedGraphWriter {
    fn write(&mut self, item: ModBasePileup) -> AnyhowResult<u64> {
        let mut rows_written = 0;
        let tab = '\t';
        for (pos, feature_counts) in item.iter_counts_sorted() {
            for (partition_key, pileup_feature_counts) in feature_counts {
                let key_name = match partition_key {
                    PartitionKey::NoKey => {
                        if self.use_groupings {
                            UNGROUPED
                        } else {
                            ""
                        }
                    }
                    PartitionKey::Key(idx) => item
                        .partition_keys
                        .get_index(*idx)
                        .map(|s| s.as_str())
                        .unwrap_or(NOT_FOUND),
                };
                for feature_count in pileup_feature_counts {
                    let key = BedGraphFileKey::new(
                        *partition_key,
                        feature_count.raw_strand,
                        feature_count.raw_mod_code,
                    );
                    let fh = self.get_writer_for_modstrand(key, key_name);
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
        }

        Ok(rows_written)
    }
}

pub struct TableWriter<W: Write> {
    writer: BufWriter<W>,
}

impl TableWriter<Stdout> {
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
                format!("pass_threshold_{}", dna_base.char()),
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
            "pass_count",
            "pass_frac",
            "all_count",
            "all_frac",
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
                    pass_counts,
                    pass_frac,
                    all_counts,
                    all_frac,
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

impl TsvWriter<Stdout> {
    pub fn new_stdout(header: Option<String>) -> Self {
        let out = BufWriter::new(std::io::stdout());
        if let Some(header) = header {
            println!("{header}");
        }

        Self { buf_writer: out }
    }
}

impl TsvWriter<File> {
    pub fn new_file(
        fp: &str,
        force: bool,
        header: Option<String>,
    ) -> AnyhowResult<Self> {
        let p = Path::new(fp);
        if p.exists() && !force {
            return Err(anyhow!("refusing to write over existing file {fp}"));
        }
        let fh = File::create(p)?;
        let mut buf_writer = BufWriter::new(fh);
        if let Some(header) = header {
            buf_writer.write(format!("{header}\n").as_bytes())?;
        }
        Ok(Self { buf_writer })
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
                    "{}_pass_calls_{}\t{}\n",
                    canonical_base.char(),
                    label,
                    counts
                ));
                report.push_str(&format!(
                    "{}_pass_frac_{}\t{}\n",
                    canonical_base.char(),
                    label,
                    counts as f64 / total_calls
                ));
                report.push_str(&format!(
                    "{}_fail_calls_{}\t{}\n",
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
                "{}_total_fail_mod_calls\t{}\n",
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

#[derive(new)]
pub(crate) struct MultiTableWriter {
    out_dir: PathBuf,
}

#[derive(new)]
pub(crate) struct SampledProbs {
    histograms: Option<HashMap<char, Histogram>>,
    percentiles: HashMap<char, Percentiles>,
    prefix: Option<String>,
}

impl SampledProbs {
    fn get_thresholds_filename(&self) -> String {
        if let Some(prefix) = &self.prefix {
            format!("{prefix}_thresholds.tsv")
        } else {
            format!("thresholds.tsv")
        }
    }

    fn get_probabilities_filenames(&self) -> (String, String) {
        if let Some(prefix) = &self.prefix {
            (
                format!("{prefix}_probabilities.tsv"),
                format!("{prefix}_probabilities.txt"),
            )
        } else {
            (format!("probabilities.tsv"), format!("probabilities.txt"))
        }
    }

    pub(crate) fn check_path(
        &self,
        p: &PathBuf,
        force: bool,
    ) -> AnyhowResult<()> {
        let filename = self.get_thresholds_filename();
        let fp = p.join(filename);
        if fp.exists() && !force {
            return Err(anyhow!("refusing to overwrite {:?}", fp));
        } else if fp.exists() && force {
            debug!("thresholds file at {:?} will be overwritten", fp);
        }
        if let Some(_) = &self.histograms {
            let (probs_table_fn, probs_plots_fn) =
                self.get_probabilities_filenames();
            let probs_table_fp = p.join(probs_table_fn);
            let probs_plots_fp = p.join(probs_plots_fn);
            for fp in [probs_table_fp, probs_plots_fp] {
                if fp.exists() && !force {
                    return Err(anyhow!("refusing to overwrite {:?}", fp));
                } else if fp.exists() && force {
                    debug!(
                        "probabilities file at {:?} will be overwritten",
                        fp
                    );
                }
            }
        }

        Ok(())
    }

    fn thresholds_table(&self) -> Table {
        let mut table = Table::new();
        table.set_format(*prettytable::format::consts::FORMAT_CLEAN);
        table.set_titles(row!["base", "percentile", "threshold"]);
        for (base, percentiles) in &self.percentiles {
            for (q, p) in percentiles.qs.iter() {
                let q = *q * 100f32;
                table.add_row(row![base, q, *p]);
            }
        }
        table
    }
}

impl OutWriter<SampledProbs> for MultiTableWriter {
    fn write(&mut self, item: SampledProbs) -> AnyhowResult<u64> {
        let mut rows_written = 0u64;
        let thresh_table = item.thresholds_table();

        let threshold_fn = self.out_dir.join(item.get_thresholds_filename());
        let mut fh = File::create(threshold_fn)?;
        let n_written = thresh_table.print(&mut fh)?;
        rows_written += n_written as u64;

        if let Some(histograms) = &item.histograms {
            let (probs_table_fn, probs_plots_fn) =
                item.get_probabilities_filenames();
            let mut probs_table_fh =
                File::create(self.out_dir.join(probs_table_fn))?;
            let mut probs_plots_fh =
                File::create(self.out_dir.join(probs_plots_fn))?;

            let mut histogram_table = Table::new();
            histogram_table
                .set_format(*prettytable::format::consts::FORMAT_CLEAN);
            histogram_table.set_titles(row![
                "code",
                "bucket",
                "range_start",
                "range_end",
                "count",
                "frac"
            ]);
            let mut total_rows = 0;
            for (raw_mod_base_code, histogram) in histograms {
                let mut row_count = Vec::new();
                for (i, bucket) in histogram.buckets().enumerate() {
                    histogram_table.add_row(row![
                        raw_mod_base_code,
                        i + 1,
                        format!("{:.3}", bucket.start()),
                        format!("{:.3}", bucket.end()),
                        bucket.count()
                    ]);
                    row_count.push(bucket.count());
                }
                let total = row_count.iter().sum::<u64>() as f32;
                for (i, count) in row_count.iter().enumerate() {
                    let frac = *count as f32 / total;
                    histogram_table
                        .get_mut_row(i + total_rows)
                        .unwrap()
                        .add_cell(cell!(frac));
                }
                total_rows += row_count.len();
            }
            let n_written = histogram_table.print(&mut probs_table_fh)?;
            rows_written += n_written as u64;

            for (raw_mod_code, hist) in histograms {
                probs_plots_fh
                    .write(format!("# code {raw_mod_code}\n").as_bytes())?;
                probs_plots_fh.write(format!("{hist}").as_bytes())?;
            }
        }

        Ok(rows_written)
    }
}

impl OutWriter<SampledProbs> for TsvWriter<Stdout> {
    fn write(&mut self, item: SampledProbs) -> AnyhowResult<u64> {
        let mut rows_written = 0u64;
        let thresholds_table = item.thresholds_table();
        let n_written = thresholds_table.print(&mut self.buf_writer)?;
        rows_written += n_written as u64;
        if let Some(histograms) = &item.histograms {
            for (raw_mod_code, hist) in histograms {
                println!("# code {raw_mod_code}");
                println!("{hist}");
            }
        }

        Ok(rows_written)
    }
}

#[derive(new)]
pub struct TsvWriterWithContigNames<W: Write> {
    tsv_writer: TsvWriter<W>,
    tid_to_name: HashMap<u32, String>,
    name_to_seq: HashMap<String, Vec<u8>>,
    written_reads: HashSet<String>,
}

impl<W: Write> OutwriterWithMemory<ReadsBaseModProfile>
    for TsvWriterWithContigNames<W>
{
    fn write(&mut self, item: ReadsBaseModProfile) -> AnyhowResult<u64> {
        let missing_chrom = ".".to_string();
        let mut rows_written = 0u64;
        for profile in item.profiles.iter() {
            if self.written_reads.contains(&profile.record_name) {
                continue;
            } else {
                let chrom_name = if let Some(chrom_id) = profile.chrom_id {
                    self.tid_to_name.get(&chrom_id)
                } else {
                    None
                };
                for mod_profile in profile.profile.iter() {
                    let row = mod_profile.to_row(
                        &profile.record_name,
                        chrom_name.unwrap_or(&missing_chrom),
                        &self.name_to_seq,
                    );
                    self.tsv_writer.buf_writer.write(row.as_bytes())?;
                    rows_written += 1;
                }
                self.written_reads.insert(profile.record_name.to_owned());
            }
        }
        Ok(rows_written)
    }

    fn num_reads(&self) -> usize {
        self.written_reads.len()
    }
}

pub struct PartitioningBedMethylWriter {
    prefix: Option<String>,
    out_dir: PathBuf,
    tabs_and_spaces: bool,
    router: FxHashMap<PartitionKey, BufWriter<File>>,
}

impl PartitioningBedMethylWriter {
    pub fn new(
        out_path: &String,
        only_tabs: bool,
        prefix: Option<&String>,
    ) -> anyhow::Result<Self> {
        let dir_path = Path::new(out_path);
        if !dir_path.is_dir() {
            info!("creating {out_path}");
            std::fs::create_dir_all(dir_path)?;
        }
        let out_dir = dir_path.to_path_buf();
        let prefix = prefix.cloned();
        let router = FxHashMap::default();
        Ok(Self {
            out_dir,
            prefix,
            router,
            tabs_and_spaces: !only_tabs,
        })
    }

    fn get_writer_for_key(
        &mut self,
        partition_key: PartitionKey,
        key_name: &str,
    ) -> &mut BufWriter<File> {
        self.router.entry(partition_key).or_insert_with(|| {
            let filename = if let Some(prefix) = self.prefix.as_ref() {
                format!("{prefix}_{key_name}.bed")
            } else {
                format!("{key_name}.bed")
            };
            let fp = self.out_dir.join(filename);
            let fh = File::create(fp).unwrap();
            BufWriter::new(fh)
        })
    }
}

const NOT_FOUND: &str = "not_found";
const UNGROUPED: &str = "ungrouped";

impl OutWriter<ModBasePileup> for PartitioningBedMethylWriter {
    fn write(&mut self, item: ModBasePileup) -> AnyhowResult<u64> {
        let tabs_and_spaces = self.tabs_and_spaces;
        let mut rows_written = 0u64;
        for (&pos, partitioned_feature_counts) in item.iter_counts_sorted() {
            for (&partition_key, pileup_feature_counts) in
                partitioned_feature_counts
            {
                let key_name = match partition_key {
                    PartitionKey::NoKey => UNGROUPED,
                    PartitionKey::Key(idx) => item
                        .partition_keys
                        .get_index(idx)
                        .map(|s| s.as_str())
                        .unwrap_or(NOT_FOUND),
                };
                let writer = self.get_writer_for_key(partition_key, key_name);
                rows_written += BedMethylWriter::write_feature_counts(
                    pos,
                    &item.chrom_name,
                    &pileup_feature_counts,
                    writer,
                    tabs_and_spaces,
                )?;
            }
        }

        Ok(rows_written)
    }
}
