use crate::errs::MkError;
use crate::mod_bam::{
    CollapseMethod, EdgeFilter, MmTagInfo, ModBaseInfo, RawModTags, SkipMode,
    ML_TAGS, MM_TAGS,
};
use crate::mod_base_code::{DnaBase, ModCodeRepr};
use crate::monoid::Moniod;
use crate::position_filter::StrandedPositionFilter;
use crate::reads_sampler::record_sampler::{Indicator, RecordSampler};
use crate::record_processor::{RecordProcessor, WithRecords};
use crate::util::{
    get_forward_sequence, get_human_readable_table, get_ticker,
    record_is_not_primary, Strand,
};
use anyhow::bail;
use derive_new::new;
use itertools::Itertools;
use log::{error, info};
use prettytable::row;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{self, Read, Records};
use rustc_hash::{FxHashMap, FxHashSet};
use std::cmp::Ordering;
use std::fs::File;
use std::ops::ControlFlow;
use std::path::PathBuf;

#[allow(non_upper_case_globals)]
pub(super) mod output_filenames {
    pub const error_counts: &str = "error_counts.tsv";
    pub const record_counts: &str = "record_counts.tsv";
    pub const valid_headers: &str = "valid_mm_headers.tsv";
    pub const invalid_headers: &str = "invalid_mm_headers.tsv";
    pub const modified_bases: &str = "modified_bases.tsv";
    pub const filenames: [&str; 5] = [
        error_counts,
        record_counts,
        valid_headers,
        invalid_headers,
        modified_bases,
    ];
}

#[derive(Default, new, Debug)]
pub(crate) struct ModTagViews {
    valid_mm_tag_headers: FxHashMap<String, usize>,
    invalid_mm_tag_headers: FxHashMap<String, usize>,
    modified_bases: FxHashMap<
        Strand,
        FxHashMap<DnaBase, FxHashSet<(ModCodeRepr, SkipMode)>>,
    >,
    error_counts: FxHashMap<String, usize>,
    ok_record_counts: u64,
    num_records: u64,
}

impl ModTagViews {
    fn add_tag_state(&mut self, tag_state: TagState) {
        match tag_state {
            TagState::Error(e) => self.add_error(e),
            TagState::ValidTagsInvalidInfo {
                mmtag_infos,
                modbase_info_err,
            } => {
                mmtag_infos
                    .iter()
                    .filter(|mmi| mmi.has_positions())
                    .map(|mmi| mmi.header())
                    .unique()
                    .for_each(|header| {
                        self.invalid_mm_tag_headers
                            .entry(header)
                            .and_modify(|x| *x = x.saturating_add(1))
                            .or_insert(1);
                    });
                self.add_error(modbase_info_err);
            }
            TagState::Valid { mmtag_infos, mod_base_info } => {
                mmtag_infos
                    .iter()
                    .filter(|mmi| mmi.has_positions())
                    .map(|mmi| mmi.header())
                    .unique()
                    .for_each(|header| {
                        self.valid_mm_tag_headers
                            .entry(header)
                            .and_modify(|x| *x = x.saturating_add(1))
                            .or_insert(1);
                    });
                for (base, strand, probs) in
                    mod_base_info.iter_seq_base_mod_probs()
                {
                    let agg = self
                        .modified_bases
                        .entry(strand)
                        .or_insert_with(FxHashMap::default)
                        .entry(base)
                        .or_insert_with(FxHashSet::default);
                    let mode = probs.skip_mode;
                    probs
                        .pos_to_base_mod_probs
                        .values()
                        .flat_map(|bmp| bmp.iter_probs().map(|(&x, _)| x))
                        .unique()
                        .map(|code| (code, mode))
                        .for_each(|t| {
                            agg.insert(t);
                        })
                }
                self.ok_record_counts = self.ok_record_counts.saturating_add(1);
            }
        }
        self.num_records = self.num_records.saturating_add(1);
    }

    fn add_error(&mut self, e: MkError) {
        self.error_counts
            .entry(e.to_string())
            .and_modify(|x| *x = x.saturating_add(1))
            .or_insert(1);
    }

    fn make_header_table(
        counts: &FxHashMap<String, usize>,
    ) -> Option<prettytable::Table> {
        if counts.is_empty() {
            None
        } else {
            let mut table = get_human_readable_table();
            table.set_titles(row!["tag_header", "count"]);
            for (header, c) in
                counts.iter().sorted_by(|(x, a), (y, b)| match b.cmp(a) {
                    Ordering::Equal => x.cmp(y),
                    o @ _ => o,
                })
            {
                table.add_row(row![header, c]);
            }
            Some(table)
        }
    }

    fn maybe_write_table(
        out_dir: Option<&PathBuf>,
        force: bool,
        file_name: &str,
        prefix: Option<&String>,
        table: &prettytable::Table,
    ) -> anyhow::Result<()> {
        if let Some(out_d) = out_dir {
            let out_fn = if let Some(p) = prefix {
                out_d.join(format!("{p}_{file_name}"))
            } else {
                out_d.join(file_name)
            };
            let fh = if force {
                File::create(out_fn)?
            } else {
                File::create_new(out_fn)?
            };
            let writer =
                csv::WriterBuilder::new().delimiter('\t' as u8).from_writer(fh);
            table.to_csv_writer(writer)?;
        }

        Ok(())
    }

    pub(super) fn report(
        &self,
        out_dir: Option<&PathBuf>,
        prefix: Option<&String>,
        force: bool,
        permissive: bool,
    ) -> anyhow::Result<()> {
        let total_error = self.error_counts.values().sum::<usize>();
        let pass_rate = self.ok_record_counts as f32 / self.num_records as f32;
        let fail_rate = total_error as f32 / self.num_records as f32;
        let msg = format!(
            "input modBAM contains {total_error} ({:.2}%) failed records",
            fail_rate * 100f32
        );
        if total_error > 0 {
            error!("{msg}");
        } else {
            info!("no errors");
        }
        info!(
            "num PASS records: {} ({:.2}%)",
            self.ok_record_counts,
            pass_rate * 100f32
        );
        info!("num records: {}", self.num_records);

        let mut err_table = get_human_readable_table();
        err_table.set_titles(row!["error", "count", "pct"]);
        let total_error = self.error_counts.values().sum::<usize>();
        for (er, c) in
            self.error_counts.iter().sorted_by(|(_, a), (_, b)| b.cmp(a))
        {
            let pct = ((*c as f32) / total_error as f32) * 100f32;
            let pct = format!("{pct:.2}");
            err_table.add_row(row![er, c, pct]);
        }
        err_table.add_row(row!["total", total_error, 100f32]);

        if !self.error_counts.is_empty() {
            info!("errors:\n{err_table}\n");
            Self::maybe_write_table(
                out_dir,
                force,
                output_filenames::error_counts,
                prefix,
                &err_table,
            )?;
        }

        let valid_headers_table =
            Self::make_header_table(&self.valid_mm_tag_headers);
        if let Some(tab) = valid_headers_table {
            info!("valid record tag headers:\n{tab}\n");
            Self::maybe_write_table(
                out_dir,
                force,
                output_filenames::valid_headers,
                prefix,
                &tab,
            )?;
        }
        let invalid_headers_table =
            Self::make_header_table(&self.invalid_mm_tag_headers);
        if let Some(tab) = invalid_headers_table {
            info!("invalid record tag headers:\n{tab}\n");
            Self::maybe_write_table(
                out_dir,
                force,
                output_filenames::invalid_headers,
                prefix,
                &tab,
            )?;
        }

        let mut mods_table = get_human_readable_table();
        mods_table.set_titles(row![
            "strand",
            "primary_base",
            "mod_code",
            "mode"
        ]);
        if let Some(pos_strand_mods) =
            self.modified_bases.get(&Strand::Positive)
        {
            for (base, codes) in
                pos_strand_mods.iter().sorted_by(|(a, _), (b, _)| b.cmp(a))
            {
                for (code, mode) in codes.iter().sorted() {
                    mods_table.add_row(row![
                        Strand::Positive.to_char(),
                        base.char(),
                        code.to_string(),
                        mode.to_string()
                    ]);
                }
            }
        }
        if let Some(neg_strand_mods) =
            self.modified_bases.get(&Strand::Negative)
        {
            for (base, codes) in
                neg_strand_mods.iter().sorted_by(|(a, _), (b, _)| a.cmp(b))
            {
                for (code, mode) in codes.iter().sorted() {
                    mods_table.add_row(row![
                        Strand::Negative.to_char(),
                        base.char(),
                        code.to_string(),
                        mode.to_string()
                    ]);
                }
            }
        }
        info!("modified bases:\n{mods_table}");
        Self::maybe_write_table(
            out_dir,
            force,
            output_filenames::modified_bases,
            prefix,
            &mods_table,
        )?;

        if total_error == 0 || permissive {
            Ok(())
        } else {
            bail!(msg)
        }
    }
}

enum TagState {
    // failed to parse the tags at all
    Error(MkError),
    // valid MM tags, but information is invalid
    ValidTagsInvalidInfo {
        mmtag_infos: Vec<MmTagInfo>,
        modbase_info_err: MkError,
    },
    Valid {
        mmtag_infos: Vec<MmTagInfo>,
        mod_base_info: ModBaseInfo,
    },
}

#[inline]
fn extract_mm_tag_info(record: &bam::Record) -> TagState {
    let raw_mod_tags = match RawModTags::new_from_record(record) {
        Ok(raw_mod_tags) => raw_mod_tags,
        Err(e) => return TagState::Error(e),
    };

    let n_mmml = record
        .aux_iter()
        .filter_ok(|(tag, _)| {
            tag == &MM_TAGS[0].as_bytes()
                || tag == &MM_TAGS[1].as_bytes()
                || tag == &ML_TAGS[0].as_bytes()
                || tag == &ML_TAGS[1].as_bytes()
        })
        .count();
    if n_mmml != 2 {
        return TagState::Error(MkError::MultipleTagInstances);
    }

    match MmTagInfo::parse_mm_tag(&raw_mod_tags.raw_mm) {
        Ok(mm_tags) => {
            let forward_sequence = get_forward_sequence(record);
            match ModBaseInfo::new(&mm_tags, &raw_mod_tags, &forward_sequence) {
                Ok(mod_base_info) => {
                    TagState::Valid { mmtag_infos: mm_tags, mod_base_info }
                }
                Err(e) => TagState::ValidTagsInvalidInfo {
                    mmtag_infos: mm_tags,
                    modbase_info_err: e,
                },
            }
        }
        Err(e) => TagState::Error(e),
    }
}

impl RecordProcessor for ModTagViews {
    type Output = Self;

    fn process_records<T: Read>(
        records: Records<T>,
        with_progress: bool,
        mut record_sampler: RecordSampler,
        _collapse_method: Option<&CollapseMethod>,
        _edge_filter: Option<&EdgeFilter>,
        _position_filter: Option<&StrandedPositionFilter<()>>,
        only_mapped: bool,
        allow_non_primary: bool,
        prev_end: Option<u32>,
        _kmer_size: Option<usize>,
    ) -> anyhow::Result<Self::Output> {
        let pb = if with_progress { Some(get_ticker()) } else { None };
        let tag_views = records
            // .progress_with(pb)
            .map(|res| res.map_err(|e| MkError::HtsLibError(e)))
            .filter_ok(|record| {
                prev_end
                    .map(|cut| record.reference_start() >= cut as i64)
                    .unwrap_or(true)
            })
            .filter_ok(
                |record| if only_mapped { !record.is_unmapped() } else { true },
            )
            .filter_ok(|record| {
                if allow_non_primary {
                    true
                } else {
                    !record_is_not_primary(record)
                }
            })
            .map(|res| res.map(|record| extract_mm_tag_info(&record)))
            .try_fold(Self::default(), |mut agg, res| {
                match record_sampler.ask() {
                    Indicator::Use(tok) => {
                        match res {
                            Ok(tag_state) => {
                                agg.add_tag_state(tag_state);
                            }
                            Err(e) => {
                                agg.add_error(e);
                                agg.num_records =
                                    agg.num_records.saturating_add(1);
                            }
                        }
                        record_sampler.used(tok);
                        if let Some(pb) = pb.as_ref() {
                            pb.inc(1)
                        };
                        ControlFlow::Continue(agg)
                    }
                    Indicator::Skip => ControlFlow::Continue(agg),
                    Indicator::Done => ControlFlow::Break(agg),
                }
            });
        let tag_views = match tag_views {
            ControlFlow::Continue(x) => x,
            ControlFlow::Break(x) => x,
        };

        Ok(tag_views)
    }
}

impl WithRecords for ModTagViews {
    fn size(&self) -> u64 {
        self.num_records
    }

    fn num_reads(&self) -> usize {
        self.num_records as usize
    }
}

impl Moniod for ModTagViews {
    fn zero() -> Self {
        Self::default()
    }

    fn op(self, other: Self) -> Self {
        let valid_mm_tag_headers =
            self.valid_mm_tag_headers.op(other.valid_mm_tag_headers);
        let invalid_mm_tag_headers =
            self.invalid_mm_tag_headers.op(other.invalid_mm_tag_headers);
        let modified_bases = self.modified_bases.op(other.modified_bases);
        let error_counts = self.error_counts.op(other.error_counts);
        let num_records = self.num_records.saturating_add(other.num_records);
        let num_ok =
            self.ok_record_counts.saturating_add(other.ok_record_counts);

        Self {
            valid_mm_tag_headers,
            invalid_mm_tag_headers,
            modified_bases,
            error_counts,
            num_records,
            ok_record_counts: num_ok,
        }
    }

    fn op_mut(&mut self, other: Self) {
        self.invalid_mm_tag_headers.op_mut(other.invalid_mm_tag_headers);
        self.valid_mm_tag_headers.op_mut(other.valid_mm_tag_headers);
        self.modified_bases.op_mut(other.modified_bases);
        self.error_counts.op_mut(other.error_counts);
        self.num_records = self.num_records.saturating_add(other.num_records);
        self.ok_record_counts =
            self.ok_record_counts.saturating_add(other.ok_record_counts);
    }

    fn len(&self) -> usize {
        self.num_records as usize
    }
}
