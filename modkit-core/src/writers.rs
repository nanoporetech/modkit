use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{stdout, BufWriter, Cursor, Stdout, Write};
use std::path::{Path, PathBuf};

use anyhow::{anyhow, bail, Context, Result as AnyhowResult};
use bitvec::order::Lsb0;
use bitvec::view::BitView;
use charming::component::{
    Axis, DataZoom, DataZoomType, Feature, Legend, Restore, SaveAsImage, Title,
    Toolbox, ToolboxDataZoom,
};
use charming::element::{
    AxisPointer, AxisPointerType, AxisType, Color, Tooltip, Trigger,
};
use charming::series::Bar;
use charming::{Chart, HtmlRenderer};
use crossbeam_channel::Sender;
use derive_new::new;
use gzp::deflate::Bgzf;
use gzp::par::compress::{ParCompress, ParCompressBuilder};
use indicatif::{MultiProgress, ProgressBar};
use itertools::Itertools;
use log::{debug, info, warn};
use prettytable::format::FormatBuilder;
use prettytable::{row, Table};
use random_color::RandomColor;
use rust_htslib::bam::HeaderView;

use crate::mod_base_code::{
    BaseState, DnaBase, ModCodeRepr, ModifiedBasesOptions, ProbHistogram,
    DNA_BASE_COLORS, MOD_COLORS,
};
use crate::pileup::bedrmod::BedRModArgs;
use crate::pileup::duplex::DuplexModBasePileup;
use crate::pileup::{ModBasePileup2, PileupFeatureCounts2};
use crate::sample_probs::PercentileThresholdQual;
use crate::summarize::ModSummary;
use crate::util::{get_ticker_with_rate, TAB};

pub trait PileupWriter<T> {
    fn write(
        &mut self,
        item: T,
        motif_labels: &[String],
    ) -> anyhow::Result<u64>;
}

pub trait OutWriter<T> {
    fn write(&mut self, item: T) -> AnyhowResult<u64>;
}

pub struct BedMethylWriter<T: Write> {
    buf_writer: BufWriter<T>,
    tabs_and_spaces: bool,
}

pub struct BedMethylWriter2<T: Write> {
    buff: Cursor<Vec<u8>>,
    inner: RecordingWriter<T>,
    return_mem: Sender<ModBasePileup2>,
    bedrmod_spec: bool,
}

impl BedMethylWriter2<BufWriter<std::io::Stdout>> {
    pub(crate) fn new_stdout(
        with_header: bool,
        bed_rmod_args: &BedRModArgs,
        header: &HeaderView,
        modified_bases_options: Option<&Vec<ModifiedBasesOptions>>,
        multi_progress: MultiProgress,
        return_mem: Sender<ModBasePileup2>,
    ) -> anyhow::Result<Self> {
        let write_pb = multi_progress.add(get_ticker_with_rate());
        write_pb.set_message(format!("B written stdout"));
        write_pb.set_position(0);
        let mut writer = RecordingWriter::new_stdout(write_pb);
        if with_header {
            writer.write(bedmethyl_header().as_bytes())?;
        } else if bed_rmod_args.enabled() {
            let modified_bases_options =
                modified_bases_options.ok_or_else(|| {
                    anyhow!("--modified-bases required for --bedrmod")
                })?;
            let bedrmod_header =
                bed_rmod_args.header(&header, modified_bases_options)?;
            writer.write(bedrmod_header.as_bytes())?;
        }

        let buff = Cursor::new(vec![0u8; 1 << 20]);
        Ok(Self {
            buff,
            inner: writer,
            return_mem,
            bedrmod_spec: bed_rmod_args.enabled(),
        })
    }
}

impl BedMethylWriter2<ParCompress<Bgzf>> {
    pub(crate) fn new_bgzf(
        out_path: &PathBuf,
        with_header: bool,
        bed_rmod_args: &BedRModArgs,
        header: &HeaderView,
        modified_bases_options: Option<&Vec<ModifiedBasesOptions>>,
        multi_progress: MultiProgress,
        return_mem: Sender<ModBasePileup2>,
        bgzf_threads: usize,
    ) -> anyhow::Result<Self> {
        let fh = File::create(out_path)?;
        let write_pb = multi_progress.add(get_ticker_with_rate());
        let out_fn = out_path
            .file_name()
            .and_then(|x| x.to_str())
            .map(|x| x.to_string())
            .unwrap_or_else(|| "failed to parse filename".to_string());
        write_pb.set_message(format!("B written to output: {out_fn}"));
        write_pb.set_position(0);
        let mut writer = RecordingWriter::new_bgzf(fh, bgzf_threads, write_pb);
        if with_header {
            writer.write(bedmethyl_header().as_bytes())?;
        } else if bed_rmod_args.enabled() {
            let modified_bases_options =
                modified_bases_options.ok_or_else(|| {
                    anyhow!("--modified-bases required for --bedrmod")
                })?;
            let bedrmod_header =
                bed_rmod_args.header(&header, modified_bases_options)?;
            writer.write(bedrmod_header.as_bytes())?;
        }
        let buff = Cursor::new(vec![0u8; 1 << 20]);
        Ok(Self {
            buff,
            inner: writer,
            return_mem,
            bedrmod_spec: bed_rmod_args.enabled(),
        })
    }
}

impl BedMethylWriter2<BufWriter<File>> {
    pub(crate) fn new(
        out_path: &PathBuf,
        with_header: bool,
        bed_rmod_args: &BedRModArgs,
        header: &HeaderView,
        modified_bases_options: Option<&Vec<ModifiedBasesOptions>>,
        multi_progress: MultiProgress,
        return_mem: Sender<ModBasePileup2>,
    ) -> anyhow::Result<Self> {
        let fh = File::create(out_path)?;
        let write_pb = multi_progress.add(get_ticker_with_rate());
        let out_fn = out_path
            .file_name()
            .and_then(|x| x.to_str())
            .map(|x| x.to_string())
            .unwrap_or_else(|| "failed to parse filename".to_string());
        write_pb.set_message(format!("B written to output: {out_fn}"));
        write_pb.set_position(0);
        let mut writer = RecordingWriter::new_file(fh, write_pb);
        if with_header {
            writer.write(bedmethyl_header().as_bytes())?;
        } else if bed_rmod_args.enabled() {
            let modified_bases_options =
                modified_bases_options.ok_or_else(|| {
                    anyhow!("--modified-bases required for --bedrmod")
                })?;
            let bedrmod_header =
                bed_rmod_args.header(&header, modified_bases_options)?;
            writer.write(bedrmod_header.as_bytes())?;
        }

        let buff = Cursor::new(vec![0u8; 1 << 20]);
        Ok(Self {
            buff,
            inner: writer,
            return_mem,
            bedrmod_spec: bed_rmod_args.enabled(),
        })
    }
}

impl<T: Write> PileupWriter<ModBasePileup2> for BedMethylWriter2<T> {
    fn write(
        &mut self,
        item: ModBasePileup2,
        _motif_labels: &[String],
    ) -> anyhow::Result<u64> {
        self.buff.set_position(0);
        let n_rows = item.position_feature_counts.len() as u64;
        let chrom_name = &item.chrom_name;
        for pfc in item.position_feature_counts.iter() {
            format_feature_counts2(
                chrom_name,
                &mut self.buff,
                pfc,
                self.bedrmod_spec,
            )
            .unwrap();
            let pos = self.buff.position() as usize;
            if pos >= 1 << 20 {
                self.inner.write(&self.buff.get_ref()[..pos]).unwrap();
                self.buff.set_position(0);
            }
        }
        let pos = self.buff.position() as usize;
        self.inner.write(&self.buff.get_ref()[..pos]).unwrap();
        let _ = self.return_mem.send(item);
        Ok(n_rows)
    }
}

pub fn bedmethyl_header() -> String {
    bedmethyl_header_op("valid_coverage")
}

pub fn bedrmod_bedmethyl_header() -> String {
    bedmethyl_header_op("coverage")
}

fn bedmethyl_header_op(coverage_field: &'static str) -> String {
    let fields = [
        "chrom",
        "chromStart",
        "chromEnd",
        "name",
        "score",
        "strand",
        "thickStart",
        "thickEnd",
        "color",
        coverage_field,
        "percent_modified",
        "count_modified",
        "count_canonical",
        "count_other_mod",
        "count_delete",
        "count_fail",
        "count_diff",
        "count_nocall",
    ];
    let fields = fields.join("\t");
    format!("#{fields}\n")
}

impl<T: Write + Sized> BedMethylWriter<T> {
    fn header() -> String {
        bedmethyl_header()
    }

    pub fn new(
        mut buf_writer: BufWriter<T>,
        tabs_and_spaces: bool,
        with_header: bool,
    ) -> anyhow::Result<Self> {
        if with_header {
            buf_writer.write(Self::header().as_bytes())?;
        }

        Ok(Self { buf_writer, tabs_and_spaces })
    }

    #[inline]
    fn write_feature_counts2(
        chrom_name: &str,
        feature_count: &PileupFeatureCounts2,
        writer: &mut BufWriter<T>,
    ) -> anyhow::Result<()> {
        let pos = feature_count.position;
        let pp1 = pos + 1;
        let mut buff = Cursor::new([0u8; 512]);
        let fraction_modified = feature_count.n_modified as f32
            / feature_count.filtered_coverage as f32;
        write!(&mut buff, "{}{TAB}", chrom_name)?;
        write!(&mut buff, "{pos}{TAB}{pp1}{TAB}")?;
        write!(&mut buff, "{}{TAB}", feature_count.mod_code)?;
        write!(&mut buff, "{}{TAB}", feature_count.filtered_coverage)?;
        write!(&mut buff, "{}{TAB}", feature_count.raw_strand)?;
        write!(&mut buff, "{pos}{TAB}{pp1}{TAB}")?;
        write!(&mut buff, "255,0,0{TAB}")?;
        write!(&mut buff, "{}{TAB}", feature_count.filtered_coverage)?;
        write!(&mut buff, "{:.2}{TAB}", fraction_modified * 100f32)?;
        write!(&mut buff, "{}{TAB}", feature_count.n_modified)?;
        write!(&mut buff, "{}{TAB}", feature_count.n_canonical)?;
        write!(&mut buff, "{}{TAB}", feature_count.n_other_modified)?;
        write!(&mut buff, "{}{TAB}", feature_count.n_delete)?;
        write!(&mut buff, "{}{TAB}", feature_count.n_filtered)?;
        write!(&mut buff, "{}{TAB}", feature_count.n_diff)?;
        write!(&mut buff, "{}\n", feature_count.n_nocall)?;
        let pos = buff.position() as usize;

        writer
            .write(&buff.get_ref()[..pos])
            .with_context(|| "failed to write row")?;

        Ok(())
    }
}

impl<T: Write> PileupWriter<ModBasePileup2> for BedMethylWriter<T> {
    fn write(
        &mut self,
        item: ModBasePileup2,
        _motif_labels: &[String],
    ) -> anyhow::Result<u64> {
        let mut rows_written = 0;
        for feature_counts in
            item.position_feature_counts.iter().filter(|x| x.is_valid())
        {
            BedMethylWriter::write_feature_counts2(
                &item.chrom_name,
                feature_counts,
                &mut self.buf_writer,
            )?;
            rows_written += 1;
        }
        std::thread::spawn(|| drop(item));
        Ok(rows_written)
    }
}

impl<T: Write> PileupWriter<DuplexModBasePileup> for BedMethylWriter<T> {
    fn write(
        &mut self,
        item: DuplexModBasePileup,
        _motif_labels: &[String],
    ) -> AnyhowResult<u64> {
        let tab = '\t';
        let space = if !self.tabs_and_spaces { tab } else { ' ' };
        let mut rows_written = 0;
        for (pos, duplex_pileup_counts) in item
            .pileup_counts
            .iter()
            // sort by position
            .sorted_by(|(a, _), (b, _)| a.cmp(b))
        {
            // sort by base
            for (base, patterns) in duplex_pileup_counts
                .pattern_counts
                .iter()
                .sorted_by(|(a, _), (b, _)| a.cmp(b))
            {
                for pattern in patterns.iter().sorted() {
                    let name = pattern.pattern_string(*base);
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
                        name,
                        pattern.valid_coverage(),
                        '.',
                        pos,
                        pos + 1,
                        "255,0,0",
                        pattern.valid_coverage(),
                        format!("{:.2}", pattern.frac_pattern() * 100f32),
                        pattern.count,
                        pattern.n_canonical,
                        pattern.n_other_pattern,
                        duplex_pileup_counts.n_delete,
                        pattern.n_fail,
                        pattern.n_diff,
                        pattern.n_nocall,
                    );
                    self.buf_writer
                        .write(row.as_bytes())
                        .with_context(|| "failed to write row")?;
                    rows_written += 1;
                }
            }
        }
        Ok(rows_written)
    }
}

pub struct MultipleMotifBedmethylWriter<T: Write> {
    writer: RecordingWriter<T>,
    return_mem: Sender<ModBasePileup2>,
    bedrmod_spec: bool,
}

impl MultipleMotifBedmethylWriter<BufWriter<std::io::Stdout>> {
    pub(crate) fn new_stdout(
        with_header: bool,
        bed_rmod_args: &BedRModArgs,
        header: &HeaderView,
        modified_bases_options: Option<&Vec<ModifiedBasesOptions>>,
        return_mem: Sender<ModBasePileup2>,
        multi_progress: MultiProgress,
    ) -> anyhow::Result<Self> {
        let mut writer = BufWriter::new(stdout());
        if with_header {
            writer.write(bedmethyl_header().as_bytes())?;
        } else if bed_rmod_args.enabled() {
            let modified_bases_options =
                modified_bases_options.ok_or_else(|| {
                    anyhow!("--modified-bases required for --bedrmod")
                })?;
            let bedrmod_header =
                bed_rmod_args.header(&header, modified_bases_options)?;
            writer.write(bedrmod_header.as_bytes())?;
        }

        let write_pb = multi_progress.add(get_ticker_with_rate());
        write_pb.set_message(format!("B written to stdout"));
        write_pb.set_position(0);
        let writer = RecordingWriter { inner: writer, pb: write_pb };

        Ok(Self { writer, return_mem, bedrmod_spec: bed_rmod_args.enabled() })
    }
}

impl MultipleMotifBedmethylWriter<BufWriter<File>> {
    pub(crate) fn new_file(
        out_path: &PathBuf,
        with_header: bool,
        bed_rmod_args: &BedRModArgs,
        header: &HeaderView,
        modified_bases_options: Option<&Vec<ModifiedBasesOptions>>,
        multi_progress: MultiProgress,
        return_mem: Sender<ModBasePileup2>,
    ) -> anyhow::Result<Self> {
        let fh = File::create(out_path)?;
        let write_pb = multi_progress.add(get_ticker_with_rate());
        let out_fn = out_path
            .file_name()
            .and_then(|x| x.to_str())
            .map(|x| x.to_string())
            .unwrap_or_else(|| "failed to parse filename".to_string());
        write_pb.set_message(format!("B written to output: {out_fn}"));
        write_pb.set_position(0);

        let mut writer = RecordingWriter::new_file(fh, write_pb);
        if with_header {
            writer.write(bedmethyl_header().as_bytes())?;
        } else if bed_rmod_args.enabled() {
            let modified_bases_options =
                modified_bases_options.ok_or_else(|| {
                    anyhow!("--modified-bases required for --bedrmod")
                })?;
            let bedrmod_header =
                bed_rmod_args.header(&header, modified_bases_options)?;
            writer.write(bedrmod_header.as_bytes())?;
        }

        Ok(Self { writer, return_mem, bedrmod_spec: bed_rmod_args.enabled() })
    }
}

impl MultipleMotifBedmethylWriter<ParCompress<Bgzf>> {
    pub(crate) fn new_bgzf(
        out_path: &PathBuf,
        with_header: bool,
        bed_rmod_args: &BedRModArgs,
        header: &HeaderView,
        modified_bases_options: Option<&Vec<ModifiedBasesOptions>>,
        multi_progress: MultiProgress,
        return_mem: Sender<ModBasePileup2>,
        bgzf_threads: usize,
    ) -> anyhow::Result<Self> {
        let fh = File::create(out_path)?;
        let write_pb = multi_progress.add(get_ticker_with_rate());
        let out_fn = out_path
            .file_name()
            .and_then(|x| x.to_str())
            .map(|x| x.to_string())
            .unwrap_or_else(|| "failed to parse filename".to_string());
        write_pb.set_message(format!("B written to output: {out_fn}"));
        write_pb.set_position(0);

        let mut writer = RecordingWriter::new_bgzf(fh, bgzf_threads, write_pb);
        if with_header {
            writer.write(bedmethyl_header().as_bytes())?;
        } else if bed_rmod_args.enabled() {
            let modified_bases_options =
                modified_bases_options.ok_or_else(|| {
                    anyhow!("--modified-bases required for --bedrmod")
                })?;
            let bedrmod_header =
                bed_rmod_args.header(&header, modified_bases_options)?;
            writer.write(bedrmod_header.as_bytes())?;
        }

        Ok(Self { writer, return_mem, bedrmod_spec: bed_rmod_args.enabled() })
    }
}

impl<T: Write> PileupWriter<ModBasePileup2>
    for MultipleMotifBedmethylWriter<T>
{
    fn write(
        &mut self,
        item: ModBasePileup2,
        motif_labels: &[String],
    ) -> anyhow::Result<u64> {
        let mut rows_written = 0u64;
        let chrom_name = &item.chrom_name;
        let mut buff = Cursor::new([0u8; 512]);
        for feature_count in item.position_feature_counts.iter() {
            let idxs = feature_count.motif_idxs.view_bits::<Lsb0>().iter_ones();
            for idx in idxs {
                let motif_label = &motif_labels[idx];
                let pos = feature_count.position;
                let pp1 = pos + 1;
                let fraction_modified = feature_count.n_modified as f32
                    / feature_count.filtered_coverage as f32;
                write!(&mut buff, "{}{TAB}", chrom_name)?;
                write!(&mut buff, "{pos}{TAB}{pp1}{TAB}")?;
                write!(
                    &mut buff,
                    "{},{motif_label}{TAB}",
                    feature_count.mod_code
                )?;
                write!(&mut buff, "{}{TAB}", feature_count.filtered_coverage)?;
                write!(&mut buff, "{}{TAB}", feature_count.raw_strand)?;
                write!(&mut buff, "{pos}{TAB}{pp1}{TAB}")?;
                write!(&mut buff, "255,0,0{TAB}")?;
                let cov_val = if self.bedrmod_spec {
                    feature_count.filtered_coverage + feature_count.n_filtered
                } else {
                    feature_count.filtered_coverage
                };
                write!(&mut buff, "{}{TAB}", cov_val)?;
                write!(&mut buff, "{:.2}{TAB}", fraction_modified * 100f32)?;
                write!(&mut buff, "{}{TAB}", feature_count.n_modified)?;
                write!(&mut buff, "{}{TAB}", feature_count.n_canonical)?;
                write!(&mut buff, "{}{TAB}", feature_count.n_other_modified)?;
                write!(&mut buff, "{}{TAB}", feature_count.n_delete)?;
                write!(&mut buff, "{}{TAB}", feature_count.n_filtered)?;
                write!(&mut buff, "{}{TAB}", feature_count.n_diff)?;
                write!(&mut buff, "{}\n", feature_count.n_nocall)?;
                let pos = buff.position() as usize;

                self.writer
                    .write(&buff.get_ref()[..pos])
                    .with_context(|| "failed to write row")?;
                buff.set_position(0);
                rows_written = rows_written.saturating_add(1);
            }
        }
        let _ = self.return_mem.send(item);
        Ok(rows_written)
    }
}

pub struct TableWriter<W: Write> {
    writer: BufWriter<W>,
}

/// Stable writer-local ordering: letter codes lexicographically, followed by
/// ChEBI identifiers numerically.
fn compare_summary_mod_codes(
    left: &ModCodeRepr,
    right: &ModCodeRepr,
) -> Ordering {
    match (left, right) {
        (ModCodeRepr::Code(left), ModCodeRepr::Code(right)) => left.cmp(right),
        (ModCodeRepr::Code(_), ModCodeRepr::ChEbi(_)) => Ordering::Less,
        (ModCodeRepr::ChEbi(_), ModCodeRepr::Code(_)) => Ordering::Greater,
        (ModCodeRepr::ChEbi(left), ModCodeRepr::ChEbi(right)) => {
            left.cmp(right)
        }
    }
}

fn summary_output_bases(item: &ModSummary<'_>) -> Vec<DnaBase> {
    [DnaBase::A, DnaBase::C, DnaBase::G, DnaBase::T]
        .into_iter()
        .filter(|base| {
            item.mod_call_counts.contains_key(base)
                || item.filtered_mod_call_counts.contains_key(base)
                || item.per_base_mod_codes.contains_key(base)
        })
        .collect()
}

fn summary_mod_codes(item: &ModSummary<'_>, base: DnaBase) -> Vec<ModCodeRepr> {
    let mut codes = item
        .per_base_mod_codes
        .get(&base)
        .into_iter()
        .flatten()
        .copied()
        .collect::<Vec<_>>();
    for counts in [
        item.mod_call_counts.get(&base),
        item.filtered_mod_call_counts.get(&base),
    ]
    .into_iter()
    .flatten()
    {
        codes.extend(counts.keys().filter_map(|state| match state {
            BaseState::Canonical(_) => None,
            BaseState::Modified(code) => Some(*code),
        }));
    }
    codes.sort_by(compare_summary_mod_codes);
    codes.dedup();
    codes
}

fn summary_fraction(numerator: u64, denominator: u64) -> f32 {
    if denominator == 0 {
        0.0
    } else {
        numerator as f32 / denominator as f32
    }
}

impl TableWriter<Stdout> {
    pub fn new() -> Self {
        let out = BufWriter::new(std::io::stdout());
        Self { writer: out }
    }
}

impl<'a, W: Write> OutWriter<ModSummary<'a>> for TableWriter<W> {
    fn write(&mut self, item: ModSummary<'a>) -> AnyhowResult<u64> {
        let output_bases = summary_output_bases(&item);
        let mut metadata_table = Table::new();
        let metadata_format =
            FormatBuilder::new().padding(1, 1).left_border('#').build();
        metadata_table.set_format(metadata_format);
        metadata_table.add_row(row!["bases", item.mod_bases()]);
        metadata_table.add_row(row!["total_reads_used", item.total_reads_used]);
        for (dna_base, reads_with_calls) in item
            .reads_with_mod_calls
            .iter()
            .sorted_by(|(a, _), (b, _)| a.cmp(b))
        {
            metadata_table.add_row(row![
                format!("count_reads_{}", dna_base.char()),
                reads_with_calls
            ]);
        }
        for (dna_base, threshold) in
            item.per_base_thresholds.iter().sorted_by(|(a, _), (b, _)| a.cmp(b))
        {
            metadata_table.add_row(row![
                format!("pass_threshold_{}", dna_base.char()),
                threshold
            ]);
        }
        if let Some(region) = item.region {
            metadata_table.add_row(row!["region", region.to_string()]);
        }
        for base in output_bases.iter().copied() {
            let codes = summary_mod_codes(&item, base);
            if !codes.is_empty() {
                metadata_table.add_row(row![
                    format!("modification_codes_for_{base}"),
                    codes.into_iter().join(",")
                ]);
            }
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

        for canonical_base in output_bases {
            let mod_codes = summary_mod_codes(&item, canonical_base);
            let pass_mod_to_counts = item.mod_call_counts.get(&canonical_base);
            let filtered_counts =
                item.filtered_mod_call_counts.get(&canonical_base);
            let total_pass_calls = pass_mod_to_counts
                .map(|counts| counts.values().sum::<u64>())
                .unwrap_or(0);
            let total_filtered_calls = filtered_counts
                .map(|counts| counts.values().sum::<u64>())
                .unwrap_or(0);
            let total_calls = total_filtered_calls + total_pass_calls;

            let mut add_row = |base_state: BaseState, label: String| {
                let pass_count = pass_mod_to_counts
                    .and_then(|counts| counts.get(&base_state))
                    .copied()
                    .unwrap_or(0);
                let filtered = filtered_counts
                    .and_then(|counts| counts.get(&base_state))
                    .copied()
                    .unwrap_or(0);
                let all_counts = pass_count.saturating_add(filtered);
                let all_frac = summary_fraction(all_counts, total_calls);
                let pass_frac = summary_fraction(pass_count, total_pass_calls);
                report_table.add_row(row![
                    canonical_base.char(),
                    label,
                    pass_count,
                    pass_frac,
                    all_counts,
                    all_frac,
                ]);
            };

            add_row(BaseState::Canonical(canonical_base), "-".to_string());
            for mod_code in mod_codes {
                add_row(BaseState::Modified(mod_code), mod_code.to_string());
            }
        }
        let mut report_emitted = report_table.print(&mut self.writer)?;
        report_emitted += emitted;
        Ok(report_emitted as u64)
    }
}

pub struct TsvWriter<W> {
    writer: W,
}

impl<T: Write> TsvWriter<T> {
    pub fn write(&mut self, raw: &[u8]) -> std::io::Result<usize> {
        self.writer.write(raw)
    }
}

impl TsvWriter<BufWriter<std::io::Sink>> {
    pub fn new_null() -> Self {
        let out = BufWriter::new(std::io::sink());
        Self { writer: out }
    }
}

impl TsvWriter<BufWriter<Stdout>> {
    pub fn new_stdout(header: Option<String>) -> Self {
        let out = BufWriter::new(std::io::stdout());
        if let Some(header) = header {
            println!("{header}");
        }

        Self { writer: out }
    }
}

impl TsvWriter<BufWriter<File>> {
    pub fn new_path(
        path: &PathBuf,
        force: bool,
        header: Option<String>,
    ) -> anyhow::Result<Self> {
        if path.exists() && !force {
            return Err(anyhow!(
                "refusing to write over existing file {path:?}"
            ));
        }
        let fh = File::create(path)?;
        let mut buf_writer = BufWriter::new(fh);
        if let Some(header) = header {
            buf_writer.write(format!("{header}\n").as_bytes())?;
        }
        Ok(Self { writer: buf_writer })
    }

    pub fn new_file(
        fp: &str,
        force: bool,
        header: Option<String>,
    ) -> AnyhowResult<Self> {
        let p = Path::new(fp).to_path_buf();
        Self::new_path(&p, force, header)
    }
}

impl TsvWriter<ParCompress<Bgzf>> {
    pub fn new_gzip(
        fp: &str,
        force: bool,
        threads: usize,
        header: Option<String>,
    ) -> anyhow::Result<Self> {
        let fp = Path::new(fp);
        let out_fh = if force {
            File::create(fp)?
        } else {
            File::create_new(fp).context("refusing to overwrite {fp:?}")?
        };
        let mut writer = ParCompressBuilder::<Bgzf>::new()
            .num_threads(threads)
            .unwrap()
            .from_writer(out_fh);
        if let Some(header) = header {
            writer.write(header.as_bytes())?;
            writer.write(&['\n' as u8])?;
        }

        Ok(Self { writer })
    }
}

impl<W: Write> OutWriter<String> for TsvWriter<W> {
    fn write(&mut self, item: String) -> anyhow::Result<u64> {
        self.writer
            .write(item.as_bytes())
            .map(|b| b as u64)
            .map_err(|e| anyhow!("{e}"))
    }
}

impl<'a, W: Write> OutWriter<ModSummary<'a>> for TsvWriter<W> {
    fn write(&mut self, item: ModSummary) -> AnyhowResult<u64> {
        warn!(
            "this output format will not be default in the next version, the \
             table output (set with --table) will become default and this \
             format will require the --tsv option"
        );
        let mut report = String::new();
        let output_bases = summary_output_bases(&item);
        let mod_called_bases = item.mod_bases();
        report.push_str(&format!("mod_bases\t{}\n", mod_called_bases));
        for (dna_base, read_count) in item
            .reads_with_mod_calls
            .iter()
            .sorted_by(|(a, _), (b, _)| a.cmp(b))
        {
            report.push_str(&format!(
                "count_reads_{}\t{}\n",
                dna_base.char(),
                read_count
            ));
        }
        for canonical_base in output_bases {
            let pass_counts = item.mod_call_counts.get(&canonical_base);
            let filtered_counts =
                item.filtered_mod_call_counts.get(&canonical_base);
            let total_calls = pass_counts
                .map(|counts| counts.values().sum::<u64>())
                .unwrap_or(0);
            let total_filtered_calls = filtered_counts
                .map(|counts| counts.values().sum::<u64>())
                .unwrap_or(0);
            let mut states = vec![(
                BaseState::Canonical(canonical_base),
                "unmodified".to_string(),
            )];
            states.extend(
                summary_mod_codes(&item, canonical_base).into_iter().map(
                    |mod_code| {
                        (
                            BaseState::Modified(mod_code),
                            format!("modified_{mod_code}"),
                        )
                    },
                ),
            );

            for (base_state, label) in states {
                let counts = pass_counts
                    .and_then(|counts| counts.get(&base_state))
                    .copied()
                    .unwrap_or(0);
                let filtered = filtered_counts
                    .and_then(|counts| counts.get(&base_state))
                    .copied()
                    .unwrap_or(0);
                let pass_fraction = if total_calls == 0 {
                    0.0
                } else {
                    counts as f64 / total_calls as f64
                };
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
                    pass_fraction
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
                total_calls
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

        self.writer.write(report.as_bytes())?;
        Ok(1)
    }
}

#[derive(new)]
pub struct MultiTableWriter {
    out_dir: PathBuf,
}

#[derive(new)]
pub struct SampledProbs {
    histograms: Option<ProbHistogram>,
    percentiles: HashMap<DnaBase, Vec<PercentileThresholdQual>>,
    prefix: Option<String>,
    primary_base_colors: HashMap<DnaBase, String>,
    mod_base_colors: HashMap<ModCodeRepr, String>,
}

impl SampledProbs {
    fn get_thresholds_filename_prefix(prefix: Option<&String>) -> String {
        if let Some(prefix) = prefix {
            format!("{prefix}_thresholds.tsv")
        } else {
            format!("thresholds.tsv")
        }
    }

    fn get_probabilities_filenames(
        prefix: Option<&String>,
    ) -> (String, String, String) {
        if let Some(prefix) = prefix {
            (
                format!("{prefix}_probabilities.tsv"),
                format!("{prefix}_counts.html"),
                format!("{prefix}_proportion.html"),
            )
        } else {
            (
                "probabilities.tsv".into(),
                "counts.html".into(),
                "proportion.html".into(),
            )
        }
    }

    fn get_thresholds_filename(&self) -> String {
        Self::get_thresholds_filename_prefix(self.prefix.as_ref())
    }

    pub fn check_files(
        p: &PathBuf,
        prefix: Option<&String>,
        force: bool,
        with_histograms: bool,
    ) -> anyhow::Result<()> {
        let filename = Self::get_thresholds_filename_prefix(prefix);
        let fp = p.join(filename);
        if fp.exists() && !force {
            return Err(anyhow!("refusing to overwrite {:?}", fp));
        } else if fp.exists() && force {
            debug!("thresholds file at {:?} will be overwritten", fp);
        }
        if with_histograms {
            let (probs_table_fn, counts_plot_fn, prop_plot_fn) =
                Self::get_probabilities_filenames(prefix);
            let probs_table_fp = p.join(probs_table_fn);
            let counts_plot_fp = p.join(counts_plot_fn);
            let prop_plot_fp = p.join(prop_plot_fn);
            for fp in [probs_table_fp, counts_plot_fp, prop_plot_fp] {
                if fp.exists() && !force {
                    bail!("refusing to overwrite {:?}", fp)
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

    pub fn check_path(&self, p: &PathBuf, force: bool) -> AnyhowResult<()> {
        Self::check_files(
            p,
            self.prefix.as_ref(),
            force,
            self.histograms.is_some(),
        )
    }

    fn thresholds_table(&self) -> Table {
        let mut table = Table::new();
        table.set_format(*prettytable::format::consts::FORMAT_CLEAN);
        table.set_titles(row!["base", "percentile", "threshold", "mod_qual"]);
        for (base, percentiles) in
            self.percentiles.iter().sorted_by(|(a, _), (b, _)| a.cmp(b))
        {
            for x in percentiles.iter() {
                let (q, p, qu) = (x.percentile, x.threshold, x.qual);
                table.add_row(row![base.char(), q, p, qu]);
            }
        }
        table
    }
}

impl ProbHistogram {
    #[inline]
    fn qual_to_bins(q: u8) -> (f32, f32) {
        let q = q as f32;
        (q / 256f32, (q + 1f32) / 256f32)
    }

    fn get_blank_chart(
        name: &str,
        qual_bins: &[u8],
        y_axis_name: &str,
    ) -> Chart {
        let categories = qual_bins
            .iter()
            .map(|x| {
                let (from, to) = Self::qual_to_bins(*x);
                let from = from * 100f32;
                let to = to * 100f32;
                format!("[{from:.2}, {to:.2})")
            })
            .collect();
        Chart::new()
            .data_zoom(DataZoom::new().type_(DataZoomType::Slider))
            .legend(Legend::new())
            .title(Title::new().text(name))
            .tooltip(Tooltip::new().trigger(Trigger::Axis).axis_pointer(
                AxisPointer::new().type_(AxisPointerType::Shadow),
            ))
            .toolbox(
                Toolbox::new().feature(
                    Feature::new()
                        .data_zoom(ToolboxDataZoom::new().y_axis_index("none"))
                        .restore(Restore::new())
                        .save_as_image(SaveAsImage::new()),
                ),
            )
            .x_axis(
                Axis::new()
                    .type_(AxisType::Category)
                    .data(categories)
                    .name("bin"),
            )
            .y_axis(Axis::new().type_(AxisType::Value).name(y_axis_name))
    }

    fn get_artifacts(
        &self,
        extra_dna_colors: &HashMap<DnaBase, String>,
        extra_mod_colors: &HashMap<ModCodeRepr, String>,
    ) -> (Table, Chart, Chart) {
        info!("preparing plots and tables");
        let mut table = Table::new();
        table.set_titles(row![
            "code",
            "primary_base",
            "range_start",
            "range_end",
            "count",
            "frac",
            "percentile_rank",
        ]);
        let bins = self
            .prob_counts
            .values()
            .flat_map(|x| x.keys())
            .unique()
            .sorted()
            .copied()
            .collect::<Vec<u8>>();
        let mut counts_chart = Self::get_blank_chart("Counts", &bins, "counts");
        let mut prop_chart =
            Self::get_blank_chart("Proportion", &bins, "proportion");
        let mut colors = Vec::new();

        let iter =
            self.prob_counts.iter().sorted_by(|((b, bs), _), ((c, cs), _)| {
                match b.cmp(c) {
                    Ordering::Equal => bs.cmp(cs),
                    o @ _ => o,
                }
            });
        for ((primary_base, base_state), counts) in iter {
            let (label, color) = match base_state {
                BaseState::Modified(x) => (
                    format!("{primary_base}:{x}"),
                    extra_mod_colors.get(x).or(MOD_COLORS.get(x)),
                ),
                BaseState::Canonical(x) => (
                    format!("{primary_base}:-"),
                    extra_dna_colors.get(x).or(DNA_BASE_COLORS.get(x)),
                ),
            };
            // dbg!(label, color);
            let color = if let Some(c) = color {
                c.to_string()
            } else {
                let mut gen = RandomColor::new();
                gen.seed(label.as_str());
                gen.to_rgb_string()
            };
            // dbg!(label, color);
            colors.push(color);
            let total = counts.values().sum::<u64>() as f64;
            // todo could this be a .scan?
            let (stats, _) = counts.iter().fold(
                (BTreeMap::new(), 0f64),
                |(mut acc, cum_sum), (b, c)| {
                    let n = *c as f64;
                    let f = n / total;
                    let cum_sum = cum_sum + n;
                    let percentile_rank =
                        ((cum_sum - (0.5f64 * n)) / total) * 100f64;
                    acc.insert(*b, (*c, f, percentile_rank));
                    (acc, cum_sum)
                },
            );

            let dat_counts = bins
                .iter()
                .map(|b| *counts.get(b).unwrap_or(&0) as i64)
                .collect::<Vec<i64>>();
            let tot = dat_counts.iter().sum::<i64>();
            let dat_prop = dat_counts
                .iter()
                .map(|x| *x as f32 / tot as f32)
                .collect::<Vec<f32>>();
            counts_chart =
                counts_chart.series(Bar::new().name(&label).data(dat_counts));
            prop_chart =
                prop_chart.series(Bar::new().name(&label).data(dat_prop));

            for (b, (count, frac, rank)) in stats {
                let (range_start, range_end) = Self::qual_to_bins(b);
                table.add_row(row![
                    base_state,
                    primary_base,
                    range_start,
                    range_end,
                    count,
                    frac,
                    rank
                ]);
            }
        }
        counts_chart = counts_chart.color(
            colors.iter().map(|c| Color::Value(c.to_string())).collect(),
        );
        prop_chart = prop_chart.color(
            colors.iter().map(|c| Color::Value(c.to_string())).collect(),
        );

        (table, counts_chart, prop_chart)
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
            let (probs_table_fn, counts_plot_fn, prop_plot_fn) =
                SampledProbs::get_probabilities_filenames(item.prefix.as_ref());
            let probs_table_fh =
                File::create(self.out_dir.join(probs_table_fn))?;
            let mut counts_plot_fh = BufWriter::new(File::create(
                self.out_dir.join(counts_plot_fn),
            )?);
            let mut prop_plot_fh =
                BufWriter::new(File::create(self.out_dir.join(prop_plot_fn))?);

            let csv_writer = csv::WriterBuilder::new()
                .has_headers(true)
                .delimiter('\t' as u8)
                .from_writer(probs_table_fh);

            let (tab, counts_chart, prop_chart) = histograms.get_artifacts(
                &item.primary_base_colors,
                &item.mod_base_colors,
            );
            tab.to_csv_writer(csv_writer)?;
            match HtmlRenderer::new("Counts", 800, 800).render(&counts_chart) {
                Ok(blob) => {
                    counts_plot_fh.write(blob.as_bytes()).map(|_x| ())?
                }
                Err(e) => debug!("failed to render counts plot, {e:?}"),
            }
            match HtmlRenderer::new("Proportions", 800, 800).render(&prop_chart)
            {
                Ok(blob) => prop_plot_fh.write(blob.as_bytes()).map(|_x| ())?,
                Err(e) => debug!("failed to render proportions plot, {e:?}"),
            }
        }

        Ok(rows_written)
    }
}

impl OutWriter<SampledProbs> for TsvWriter<BufWriter<Stdout>> {
    fn write(&mut self, item: SampledProbs) -> AnyhowResult<u64> {
        let mut rows_written = 0u64;
        let thresholds_table = item.thresholds_table();
        let n_written = thresholds_table.print(&mut self.writer)?;
        rows_written += n_written as u64;
        Ok(rows_written)
    }
}

#[inline]
fn format_feature_counts2(
    chrom_name: &String,
    buff: &mut Cursor<Vec<u8>>,
    feature_count: &PileupFeatureCounts2,
    bedrmod_spec: bool,
) -> anyhow::Result<()> {
    let pos = feature_count.position;
    let pp1 = pos + 1;
    let fraction_modified = feature_count.n_modified as f32
        / feature_count.filtered_coverage as f32;
    write!(buff, "{}{TAB}", chrom_name)?;
    write!(buff, "{pos}{TAB}{pp1}{TAB}")?;
    write!(buff, "{}{TAB}", feature_count.mod_code)?;
    write!(buff, "{}{TAB}", feature_count.filtered_coverage)?;
    write!(buff, "{}{TAB}", feature_count.raw_strand)?;
    write!(buff, "{pos}{TAB}{pp1}{TAB}")?;
    write!(buff, "255,0,0{TAB}")?;
    let cov_val = if bedrmod_spec {
        feature_count.filtered_coverage + feature_count.n_filtered
    } else {
        feature_count.filtered_coverage
    };
    write!(buff, "{}{TAB}", cov_val)?;
    write!(buff, "{:.2}{TAB}", fraction_modified * 100f32)?;
    write!(buff, "{}{TAB}", feature_count.n_modified)?;
    write!(buff, "{}{TAB}", feature_count.n_canonical)?;
    write!(buff, "{}{TAB}", feature_count.n_other_modified)?;
    write!(buff, "{}{TAB}", feature_count.n_delete)?;
    write!(buff, "{}{TAB}", feature_count.n_filtered)?;
    write!(buff, "{}{TAB}", feature_count.n_diff)?;
    write!(buff, "{}\n", feature_count.n_nocall)?;
    Ok(())
}

pub(crate) struct RecordingWriter<T: Write> {
    inner: T,
    pb: ProgressBar,
}

impl RecordingWriter<BufWriter<std::io::Stdout>> {
    pub(crate) fn new_stdout(pb: ProgressBar) -> Self {
        Self { inner: BufWriter::with_capacity(1 << 20, std::io::stdout()), pb }
    }
}

impl RecordingWriter<BufWriter<File>> {
    pub(crate) fn new_file(file: File, pb: ProgressBar) -> Self {
        Self { inner: BufWriter::with_capacity(1 << 20, file), pb }
    }
}

impl RecordingWriter<ParCompress<Bgzf>> {
    pub(crate) fn new_bgzf(
        file: File,
        compress_threads: usize,
        pb: ProgressBar,
    ) -> Self {
        let inner = ParCompressBuilder::<Bgzf>::new()
            .num_threads(compress_threads)
            .unwrap()
            .from_writer(file);
        Self { inner, pb }
    }
}
impl<T: Write> RecordingWriter<T> {
    pub(crate) fn write(&mut self, bulk: &[u8]) -> anyhow::Result<()> {
        let n = bulk.len();
        self.inner.write_all(bulk)?;
        self.pb.inc(n as u64);
        self.flush()?;
        Ok(())
    }

    pub(crate) fn flush(&mut self) -> anyhow::Result<()> {
        self.inner.flush()?;
        Ok(())
    }
}

impl<T: Write> Drop for RecordingWriter<T> {
    fn drop(&mut self) {
        self.pb.finish_and_clear();
        let _ = self.inner.flush();
    }
}

pub struct PhasedBedMethylWriter<T: Write> {
    hp1_writer: RecordingWriter<T>,
    hp2_writer: RecordingWriter<T>,
    combined_writer: RecordingWriter<T>,
    return_mem: Sender<ModBasePileup2>,
}

impl PhasedBedMethylWriter<BufWriter<File>> {
    fn make_writer(
        out_dir: &PathBuf,
        prefix: &str,
        hp_label: &str,
        force: bool,
        multi_progress: MultiProgress,
    ) -> anyhow::Result<RecordingWriter<BufWriter<File>>> {
        let filename = format!("{prefix}{hp_label}.bedmethyl");
        let path = out_dir.join(Path::new(&filename));
        let file =
            if force { File::create(path) } else { File::create_new(path) }?;
        let write_pb = multi_progress.add(get_ticker_with_rate());
        write_pb.set_message(format!("B written: {hp_label}"));
        write_pb.set_position(0);
        let writer = RecordingWriter::new_file(file, write_pb);
        Ok(writer)
    }

    pub fn new_file(
        out_dir: &PathBuf,
        prefix: Option<&String>,
        force: bool,
        multi_progress: MultiProgress,
        return_mem: Sender<ModBasePileup2>,
    ) -> anyhow::Result<Self> {
        let prefix = prefix.map(|x| format!("{x}_")).unwrap_or("".to_string());
        let hp1_writer = Self::make_writer(
            &out_dir,
            &prefix,
            "hp1",
            force,
            multi_progress.clone(),
        )?;
        let hp2_writer = Self::make_writer(
            &out_dir,
            &prefix,
            "hp2",
            force,
            multi_progress.clone(),
        )?;
        let combined_writer = Self::make_writer(
            &out_dir,
            &prefix,
            "combined",
            force,
            multi_progress.clone(),
        )?;
        Ok(Self { hp1_writer, hp2_writer, combined_writer, return_mem })
    }
}

impl PhasedBedMethylWriter<ParCompress<Bgzf>> {
    fn make_writer(
        out_dir: &PathBuf,
        prefix: &str,
        hp_label: &str,
        force: bool,
        compression_threads: usize,
        multi_progress: MultiProgress,
    ) -> anyhow::Result<RecordingWriter<ParCompress<Bgzf>>> {
        let filename = format!("{prefix}{hp_label}.bed.gz");
        let path = out_dir.join(Path::new(&filename));
        let file =
            if force { File::create(path) } else { File::create_new(path) }?;
        let write_pb = multi_progress.add(get_ticker_with_rate());
        write_pb.set_message(format!("B written: {hp_label}"));
        write_pb.set_position(0);
        let writer =
            RecordingWriter::new_bgzf(file, compression_threads, write_pb);
        Ok(writer)
    }

    pub fn new_bgzf(
        out_dir: &PathBuf,
        prefix: Option<&String>,
        force: bool,
        multi_progress: MultiProgress,
        return_mem: Sender<ModBasePileup2>,
        compression_threads: usize,
    ) -> anyhow::Result<Self> {
        let prefix = prefix.map(|x| format!("{x}_")).unwrap_or("".to_string());
        let hp1_writer = Self::make_writer(
            &out_dir,
            &prefix,
            "hp1",
            force,
            compression_threads,
            multi_progress.clone(),
        )?;
        let hp2_writer = Self::make_writer(
            &out_dir,
            &prefix,
            "hp2",
            force,
            compression_threads,
            multi_progress.clone(),
        )?;
        let combined_writer = Self::make_writer(
            &out_dir,
            &prefix,
            "combined",
            force,
            compression_threads,
            multi_progress.clone(),
        )?;
        Ok(Self { hp1_writer, hp2_writer, combined_writer, return_mem })
    }
}

impl<T> PileupWriter<ModBasePileup2> for PhasedBedMethylWriter<T>
where
    T: Write + Send,
{
    fn write(
        &mut self,
        item: ModBasePileup2,
        _motif_labels: &[String],
    ) -> anyhow::Result<u64> {
        let chrom_name = &item.chrom_name;
        let combined_counts = &item.position_feature_counts;
        let [hp1, hp2] = &item.phased_feature_counts;
        let total_rows = combined_counts.len() + hp1.len() + hp2.len();

        // TODO: make the "buff"s part of the object.
        std::thread::scope(|scope| {
            let hp1_handle = scope.spawn(|| {
                let mut buff = Cursor::new(vec![0u8; 1 << 20]);
                for pfc in hp1.iter().filter(|x| x.is_valid()) {
                    format_feature_counts2(chrom_name, &mut buff, pfc, false)
                        .unwrap();
                    let pos = buff.position() as usize;
                    if pos >= 1 << 20 {
                        self.hp1_writer.write(&buff.get_ref()[..pos]).unwrap();
                        buff.set_position(0);
                    }
                }
                let pos = buff.position() as usize;
                self.hp1_writer.write(&buff.get_ref()[..pos]).unwrap();
            });
            let hp2_handle = scope.spawn(|| {
                let mut buff = Cursor::new(vec![0u8; 1 << 20]);
                for pfc in hp2.iter().filter(|x| x.is_valid()) {
                    format_feature_counts2(chrom_name, &mut buff, pfc, false)
                        .unwrap();
                    let pos = buff.position() as usize;
                    if pos >= 1 << 20 {
                        self.hp2_writer.write(&buff.get_ref()[..pos]).unwrap();
                        buff.set_position(0);
                    }
                }
                let pos = buff.position() as usize;
                self.hp2_writer.write(&buff.get_ref()[..pos]).unwrap();
            });
            let combined_handle = scope.spawn(|| {
                let mut buff = Cursor::new(vec![0u8; 1 << 20]);
                for pfc in combined_counts.iter().filter(|x| x.is_valid()) {
                    format_feature_counts2(&chrom_name, &mut buff, pfc, false)
                        .unwrap();
                    let pos = buff.position() as usize;
                    if pos >= 1 << 20 {
                        self.combined_writer
                            .write(&buff.get_ref()[..pos])
                            .unwrap();
                        buff.set_position(0);
                    }
                }
                let pos = buff.position() as usize;
                self.combined_writer.write(&buff.get_ref()[..pos]).unwrap();
            });
            let _ = hp1_handle.join().unwrap();
            let _ = hp2_handle.join().unwrap();
            let _ = combined_handle.join().unwrap();
        });
        let _ = self.return_mem.send(item);

        Ok(total_rows as u64)
    }
}

#[cfg(test)]
mod tests {
    use std::collections::{HashMap, HashSet};
    use std::io::{BufWriter, Write};

    use super::{OutWriter, TableWriter, TsvWriter};
    use crate::mod_base_code::{BaseState, DnaBase, ModCodeRepr};
    use crate::summarize::ModSummary;

    const EXPECTED_TSV: &str = concat!(
        "mod_bases\tA,C,G,T\n",
        "count_reads_A\t1\n",
        "count_reads_C\t2\n",
        "count_reads_G\t3\n",
        "count_reads_T\t4\n",
        "A_pass_calls_unmodified\t8\n",
        "A_pass_frac_unmodified\t1\n",
        "A_fail_calls_unmodified\t1\n",
        "A_pass_calls_modified_a\t0\n",
        "A_pass_frac_modified_a\t0\n",
        "A_fail_calls_modified_a\t0\n",
        "A_total_mod_calls\t8\n",
        "A_total_fail_mod_calls\t1\n",
        "C_pass_calls_unmodified\t6\n",
        "C_pass_frac_unmodified\t0.75\n",
        "C_fail_calls_unmodified\t1\n",
        "C_pass_calls_modified_f\t0\n",
        "C_pass_frac_modified_f\t0\n",
        "C_fail_calls_modified_f\t0\n",
        "C_pass_calls_modified_h\t0\n",
        "C_pass_frac_modified_h\t0\n",
        "C_fail_calls_modified_h\t0\n",
        "C_pass_calls_modified_m\t2\n",
        "C_pass_frac_modified_m\t0.25\n",
        "C_fail_calls_modified_m\t1\n",
        "C_total_mod_calls\t8\n",
        "C_total_fail_mod_calls\t2\n",
        "G_pass_calls_unmodified\t5\n",
        "G_pass_frac_unmodified\t1\n",
        "G_fail_calls_unmodified\t0\n",
        "G_total_mod_calls\t5\n",
        "G_total_fail_mod_calls\t0\n",
        "T_pass_calls_unmodified\t4\n",
        "T_pass_frac_unmodified\t1\n",
        "T_fail_calls_unmodified\t1\n",
        "T_total_mod_calls\t4\n",
        "T_total_fail_mod_calls\t1\n",
        "total_reads_used\t10\n",
    );

    const EXPECTED_TABLE: &str = concat!(
        "# bases                     A,C,G,T \n",
        "# total_reads_used          10 \n",
        "# count_reads_A             1 \n",
        "# count_reads_C             2 \n",
        "# count_reads_G             3 \n",
        "# count_reads_T             4 \n",
        "# pass_threshold_A          0.1 \n",
        "# pass_threshold_C          0.2 \n",
        "# pass_threshold_G          0.3 \n",
        "# pass_threshold_T          0.4 \n",
        "# modification_codes_for_A  a \n",
        "# modification_codes_for_C  f,h,m \n",
        " base  code  pass_count  pass_frac  all_count  all_frac \n",
        " A     -     8           1          9          1 \n",
        " A     a     0           0          0          0 \n",
        " C     -     6           0.75       7          0.7 \n",
        " C     f     0           0          0          0 \n",
        " C     h     0           0          0          0 \n",
        " C     m     2           0.25       3          0.3 \n",
        " G     -     5           1          5          1 \n",
        " T     -     4           1          5          1 \n",
    );

    const EXPECTED_MIXED_CODE_TABLE: &str = concat!(
        "# bases                     C \n",
        "# total_reads_used          1 \n",
        "# count_reads_C             1 \n",
        "# pass_threshold_C          0.5 \n",
        "# modification_codes_for_C  m,123 \n",
        " base  code  pass_count  pass_frac  all_count  all_frac \n",
        " C     -     1           0.5        1          0.5 \n",
        " C     m     0           0          0          0 \n",
        " C     123   1           0.5        1          0.5 \n",
    );

    const EXPECTED_MIXED_CODE_TSV: &str = concat!(
        "mod_bases\tC\n",
        "count_reads_C\t1\n",
        "C_pass_calls_unmodified\t1\n",
        "C_pass_frac_unmodified\t0.5\n",
        "C_fail_calls_unmodified\t0\n",
        "C_pass_calls_modified_m\t0\n",
        "C_pass_frac_modified_m\t0\n",
        "C_fail_calls_modified_m\t0\n",
        "C_pass_calls_modified_123\t1\n",
        "C_pass_frac_modified_123\t0.5\n",
        "C_fail_calls_modified_123\t0\n",
        "C_total_mod_calls\t2\n",
        "C_total_fail_mod_calls\t0\n",
        "total_reads_used\t1\n",
    );

    const EXPECTED_FILTERED_ONLY_TSV: &str = concat!(
        "mod_bases\tC\n",
        "count_reads_C\t1\n",
        "C_pass_calls_unmodified\t0\n",
        "C_pass_frac_unmodified\t0\n",
        "C_fail_calls_unmodified\t1\n",
        "C_pass_calls_modified_m\t0\n",
        "C_pass_frac_modified_m\t0\n",
        "C_fail_calls_modified_m\t3\n",
        "C_pass_calls_modified_123\t0\n",
        "C_pass_frac_modified_123\t0\n",
        "C_fail_calls_modified_123\t2\n",
        "C_total_mod_calls\t0\n",
        "C_total_fail_mod_calls\t6\n",
        "total_reads_used\t1\n",
    );

    const EXPECTED_FILTERED_ONLY_TABLE: &str = concat!(
        "# bases                     C \n",
        "# total_reads_used          1 \n",
        "# count_reads_C             1 \n",
        "# pass_threshold_C          0.5 \n",
        "# modification_codes_for_C  m,123 \n",
        " base  code  pass_count  pass_frac  all_count  all_frac \n",
        " C     -     0           0          1          0.16666667 \n",
        " C     m     0           0          3          0.5 \n",
        " C     123   0           0          2          0.33333334 \n",
    );

    fn summary_with_insertion_order(reverse: bool) -> ModSummary<'static> {
        let base_order = if reverse {
            [DnaBase::T, DnaBase::G, DnaBase::C, DnaBase::A]
        } else {
            [DnaBase::A, DnaBase::C, DnaBase::G, DnaBase::T]
        };

        let mut reads_with_mod_calls = HashMap::new();
        let mut mod_call_counts = HashMap::new();
        let mut filtered_mod_call_counts = HashMap::new();
        let mut per_base_thresholds = HashMap::new();
        let mut per_base_mod_codes = HashMap::new();
        for base in base_order {
            reads_with_mod_calls.insert(
                base,
                match base {
                    DnaBase::A => 1,
                    DnaBase::C => 2,
                    DnaBase::G => 3,
                    DnaBase::T => 4,
                },
            );
            per_base_thresholds.insert(
                base,
                match base {
                    DnaBase::A => 0.1,
                    DnaBase::C => 0.2,
                    DnaBase::G => 0.3,
                    DnaBase::T => 0.4,
                },
            );

            let canonical = BaseState::Canonical(base);
            let mut pass = HashMap::new();
            let mut filtered = HashMap::new();
            if base == DnaBase::C && reverse {
                pass.insert(BaseState::Modified('m'.into()), 2);
                pass.insert(canonical, 6);
                filtered.insert(BaseState::Modified('m'.into()), 1);
                filtered.insert(canonical, 1);
            } else {
                pass.insert(
                    canonical,
                    match base {
                        DnaBase::A => 8,
                        DnaBase::C => 6,
                        DnaBase::G => 5,
                        DnaBase::T => 4,
                    },
                );
                if base == DnaBase::C {
                    pass.insert(BaseState::Modified('m'.into()), 2);
                }
                filtered.insert(
                    canonical,
                    match base {
                        DnaBase::A | DnaBase::C | DnaBase::T => 1,
                        DnaBase::G => 0,
                    },
                );
                if base == DnaBase::C {
                    filtered.insert(BaseState::Modified('m'.into()), 1);
                }
            }
            mod_call_counts.insert(base, pass);
            filtered_mod_call_counts.insert(base, filtered);

            if matches!(base, DnaBase::A | DnaBase::C) {
                let codes: Vec<ModCodeRepr> = match (base, reverse) {
                    (DnaBase::A, false) => vec!['a'.into()],
                    (DnaBase::A, true) => vec!['a'.into()],
                    (DnaBase::C, false) => {
                        vec!['f'.into(), 'h'.into(), 'm'.into()]
                    }
                    (DnaBase::C, true) => {
                        vec!['m'.into(), 'h'.into(), 'f'.into()]
                    }
                    _ => unreachable!(),
                };
                per_base_mod_codes.insert(base, HashSet::from_iter(codes));
            }
        }

        ModSummary::new(
            reads_with_mod_calls,
            mod_call_counts,
            filtered_mod_call_counts,
            10,
            per_base_thresholds,
            None,
            per_base_mod_codes,
        )
    }

    fn render_tsv(summary: ModSummary<'static>) -> String {
        let mut writer = TsvWriter { writer: Vec::new() };
        OutWriter::write(&mut writer, summary).unwrap();
        String::from_utf8(writer.writer).unwrap()
    }

    fn render_table(summary: ModSummary<'static>) -> String {
        let mut writer =
            TableWriter { writer: BufWriter::new(Vec::<u8>::new()) };
        OutWriter::write(&mut writer, summary).unwrap();
        writer.writer.flush().unwrap();
        String::from_utf8(writer.writer.into_inner().unwrap()).unwrap()
    }

    fn mixed_code_summary(reverse: bool) -> ModSummary<'static> {
        let canonical = BaseState::Canonical(DnaBase::C);
        let chebi = BaseState::Modified(ModCodeRepr::ChEbi(123));
        let pass_counts = if reverse {
            HashMap::from([(chebi, 1), (canonical, 1)])
        } else {
            HashMap::from([(canonical, 1), (chebi, 1)])
        };
        let mod_codes = if reverse {
            HashSet::from([ModCodeRepr::ChEbi(123), 'm'.into()])
        } else {
            HashSet::from(['m'.into(), ModCodeRepr::ChEbi(123)])
        };
        ModSummary::new(
            HashMap::from([(DnaBase::C, 1)]),
            HashMap::from([(DnaBase::C, pass_counts)]),
            HashMap::new(),
            1,
            HashMap::from([(DnaBase::C, 0.5)]),
            None,
            HashMap::from([(DnaBase::C, mod_codes)]),
        )
    }

    fn filtered_only_summary() -> ModSummary<'static> {
        let canonical = BaseState::Canonical(DnaBase::C);
        let methyl = BaseState::Modified('m'.into());
        let chebi = BaseState::Modified(ModCodeRepr::ChEbi(123));
        ModSummary::new(
            HashMap::from([(DnaBase::C, 1)]),
            HashMap::new(),
            HashMap::from([(
                DnaBase::C,
                HashMap::from([(canonical, 1), (methyl, 3), (chebi, 2)]),
            )]),
            1,
            HashMap::from([(DnaBase::C, 0.5)]),
            None,
            HashMap::from([(
                DnaBase::C,
                HashSet::from(['m'.into(), ModCodeRepr::ChEbi(123)]),
            )]),
        )
    }

    #[test]
    fn summary_tsv_is_deterministic_across_insertion_orders() {
        let forward = render_tsv(summary_with_insertion_order(false));
        let reverse = render_tsv(summary_with_insertion_order(true));

        assert_eq!(forward, reverse);
        assert_eq!(forward, EXPECTED_TSV);
    }

    #[test]
    fn summary_table_is_deterministic_across_insertion_orders() {
        let forward = render_table(summary_with_insertion_order(false));
        let reverse = render_table(summary_with_insertion_order(true));

        assert_eq!(forward, reverse);
        assert_eq!(forward, EXPECTED_TABLE);
    }

    #[test]
    fn summary_table_uses_one_code_and_chebi_order() {
        let forward = render_table(mixed_code_summary(false));
        let reverse = render_table(mixed_code_summary(true));
        assert_eq!(forward, reverse);
        assert_eq!(forward, EXPECTED_MIXED_CODE_TABLE);
    }

    #[test]
    fn summary_tsv_uses_one_code_and_chebi_order() {
        let forward = render_tsv(mixed_code_summary(false));
        let reverse = render_tsv(mixed_code_summary(true));
        assert_eq!(forward, reverse);
        assert_eq!(forward, EXPECTED_MIXED_CODE_TSV);
    }

    #[test]
    fn summary_outputs_filtered_only_states_exactly() {
        assert_eq!(
            render_tsv(filtered_only_summary()),
            EXPECTED_FILTERED_ONLY_TSV
        );
        assert_eq!(
            render_table(filtered_only_summary()),
            EXPECTED_FILTERED_ONLY_TABLE
        );
    }
}
