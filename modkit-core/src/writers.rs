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

        let iter = item.per_base_mod_codes.into_iter().map(
            |(primary_base, mod_codes)| {
                let pass_counts = item.mod_call_counts.get(&primary_base);
                let filtered_counts =
                    item.filtered_mod_call_counts.get(&primary_base);
                (primary_base, pass_counts, filtered_counts, mod_codes)
            },
        );
        for (
            canonical_base,
            pass_mod_to_counts,
            filtered_counts,
            mut mod_codes,
        ) in iter
        {
            let total_pass_calls = pass_mod_to_counts
                .map(|counts| counts.values().sum::<u64>())
                .unwrap_or(0);
            let total_filtered_calls = filtered_counts
                .map(|counts| counts.values().sum::<u64>())
                .unwrap_or(0);
            let total_calls = total_filtered_calls + total_pass_calls;

            let mut seen_canonical = false;
            if let Some(pass_counts) = pass_mod_to_counts {
                for (base_state, pass_counts) in
                    pass_counts.iter().sorted_by(|(a, _), (b, _)| a.cmp(b))
                {
                    let label = match base_state {
                        BaseState::Canonical(_) => {
                            seen_canonical = true;
                            format!("-") // could be a const..
                        }
                        BaseState::Modified(repr) => {
                            mod_codes.remove(repr);
                            format!("{repr}")
                        }
                    };
                    let filtered = *item
                        .filtered_mod_call_counts
                        .get(&canonical_base)
                        .and_then(|filtered_counts| {
                            filtered_counts.get(&base_state)
                        })
                        .unwrap_or(&0);
                    let all_counts = *pass_counts + filtered;
                    let all_frac = all_counts as f32 / total_calls as f32;
                    let pass_frac =
                        *pass_counts as f32 / total_pass_calls as f32;
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

            if !seen_canonical {
                report_table.add_row(row![
                    canonical_base.char(),
                    format!("-"),
                    0u64,
                    0f32,
                    0u64,
                    0f32
                ]);
            }
            for mod_code in mod_codes {
                report_table.add_row(row![
                    canonical_base.char(),
                    format!("{mod_code}"),
                    0u64,
                    0f32,
                    0u64,
                    0f32
                ]);
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
            for (base_state, counts) in mod_counts {
                let label = match base_state {
                    BaseState::Canonical(_) => format!("unmodified"),
                    BaseState::Modified(repr) => format!("modified_{repr}"),
                };
                let filtered = *item
                    .filtered_mod_call_counts
                    .get(&canonical_base)
                    .and_then(|filtered_counts| {
                        filtered_counts.get(&base_state)
                    })
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

struct RecordingWriter<T: Write> {
    inner: T,
    pb: ProgressBar,
}

impl RecordingWriter<BufWriter<std::io::Stdout>> {
    fn new_stdout(pb: ProgressBar) -> Self {
        Self { inner: BufWriter::with_capacity(1 << 20, std::io::stdout()), pb }
    }
}

impl RecordingWriter<BufWriter<File>> {
    fn new_file(file: File, pb: ProgressBar) -> Self {
        Self { inner: BufWriter::with_capacity(1 << 20, file), pb }
    }
}

impl RecordingWriter<ParCompress<Bgzf>> {
    fn new_bgzf(file: File, compress_threads: usize, pb: ProgressBar) -> Self {
        let inner = ParCompressBuilder::<Bgzf>::new()
            .num_threads(compress_threads)
            .unwrap()
            .from_writer(file);
        Self { inner, pb }
    }
}
impl<T: Write> RecordingWriter<T> {
    fn write(&mut self, bulk: &[u8]) -> anyhow::Result<()> {
        let n = bulk.len();
        self.inner.write_all(bulk)?;
        self.pb.inc(n as u64);
        self.flush()?;
        Ok(())
    }

    fn flush(&mut self) -> anyhow::Result<()> {
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
