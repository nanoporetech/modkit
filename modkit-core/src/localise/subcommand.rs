use std::fs::File;
use std::io::{stdout, BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::sync::Arc;

use anyhow::{anyhow, bail, Context};
use clap::Args;
use indicatif::{MultiProgress, ParallelProgressIterator, ProgressIterator};
use log::{debug, info};
use rayon::prelude::*;
use rustc_hash::FxHashMap;

use modkit_logging::init_logging;

use crate::dmr::bedmethyl::BedMethylLine;
use crate::localise::util::{LocalizedModCounts, StrandedFeatures};
use crate::monoid::Moniod;
use crate::tabix::HtsTabixHandler;
use crate::util::{
    get_master_progress_bar, get_ticker, read_sequence_lengths_file,
    GenomeRegion, StrandRule,
};

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct EntryLocalize {
    /// Input bedMethyl table. Should be bgzip-compressed and have an
    /// associated Tabix index. The tabix index will be assumed to be
    /// $this_file.tbi
    in_bedmethyl: PathBuf,
    /// BED file of regions to calculate enrichment around. These BED records
    /// serve as the points from which the `--window` number of bases is
    /// centered.
    #[arg(long)]
    regions: PathBuf,
    /// Create plots showing %-modification vs. offset. Argument should be a
    /// path to a file.
    #[clap(help_heading = "Output Options")]
    #[arg(long = "chart")]
    chart_filepath: Option<PathBuf>,
    /// Give the HTML document and chart a name.
    #[clap(help_heading = "Output Options")]
    #[arg(long = "name", requires = "chart_filepath")]
    chart_name: Option<String>,
    /// Number of base pairs to search around, for example if your BED region
    /// records are single positions, a window of 500 will look 500 base
    /// pairs upstream and downstream of that position. If your region BED
    /// records are larger regions, this will expand from the midpoint of
    /// that region.
    #[arg(short = 'w', long = "window", default_value_t = 2000)]
    expand_window: u64,
    // todo
    // /// Expand the aggregation window `expand-window` base pairs from the
    // start /// and end of each region instead of the midpoint.
    // #[arg(long, alias = "widen")]
    // expand_from_edge: bool,
    // /// Only aggregate for a single modification code, default is to
    // aggregate /// for all modification codes at once.
    // #[arg(long)]
    // mod_code: Option<String>,
    /// Whether to only keep bedMethyl records on the "same" strand or
    /// "opposite" strand.
    #[arg(short = 's', long)]
    stranded: Option<StrandedFeatures>,
    /// Force use bedMethyl records from a particular strand, default is to use
    /// the strand as given in the BED file (will use BOTH for BED3).
    #[arg(long)]
    stranded_features: Option<StrandRule>,
    /// Minimum valid coverage to use a bedMethyl record
    #[arg(long, default_value_t = 3)]
    min_coverage: u64,
    /// TSV of genome sizes, should be <chrom>\t<size_in_bp>
    #[arg(long, short = 'r')]
    genome_sizes: PathBuf,
    /// Optionally specify a file to write output to, default is stdout.
    #[arg(long, short = 'o')]
    out_file: Option<PathBuf>,
    /// Specify a file to write debug logs to.
    #[clap(help_heading = "Logging Options")]
    #[arg(long, alias = "log")]
    log_filepath: Option<PathBuf>,
    /// Number of threads to use.
    #[clap(help_heading = "Compute Options")]
    #[arg(short = 't', long, default_value_t = 4)]
    threads: usize,
    /// Number of tabix/bgzf IO threads to use.
    #[clap(help_heading = "Compute Options")]
    #[arg(long, default_value_t = 2)]
    io_threads: usize,
    /// Force overwrite of existing output file.
    #[clap(help_heading = "Output Options")]
    #[arg(long, short = 'f', default_value_t = false)]
    force: bool,
    #[clap(help_heading = "Compute Options")]
    #[arg(
        long = "batch-size",
        hide_short_help = true,
        default_value_t = 500_000
    )]
    batch_size_bp: u64,
}

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
enum BedDelimiter {
    Tabs,
    Whitespace,
}

impl BedDelimiter {
    fn from_first_data_line(line: &str) -> Self {
        if line.contains('\t') {
            Self::Tabs
        } else {
            Self::Whitespace
        }
    }

    fn matches(self, line: &str) -> bool {
        match self {
            Self::Tabs => line.contains('\t'),
            Self::Whitespace => !line.contains('\t'),
        }
    }

    fn split<'a>(self, line: &'a str) -> Vec<&'a str> {
        match self {
            Self::Tabs => line.split('\t').collect(),
            Self::Whitespace => line.split_whitespace().collect(),
        }
    }
}

fn parse_region_core_fields(fields: &[&str]) -> anyhow::Result<GenomeRegion> {
    if fields.len() < 3 {
        bail!("expected at least 3 BED core columns, found {}", fields.len())
    }
    let chrom = fields[0];
    if chrom.is_empty() || chrom.chars().any(char::is_whitespace) {
        bail!("invalid chromosome name {chrom:?}")
    }
    let start = fields[1].parse::<u64>().map_err(|_| {
        anyhow!(
            "invalid start coordinate {:?}, expected a complete u64",
            fields[1]
        )
    })?;
    let end = fields[2].parse::<u64>().map_err(|_| {
        anyhow!(
            "invalid end coordinate {:?}, expected a complete u64",
            fields[2]
        )
    })?;
    let name = fields
        .get(3)
        .map(|name| {
            if name.is_empty() {
                bail!("BED name must not be empty")
            }
            Ok((*name).to_string())
        })
        .transpose()?;
    if let Some(score) = fields.get(4) {
        if *score != "." {
            let score_value = score.parse::<f64>().map_err(|_| {
                anyhow!(
                    "invalid BED score {score:?}, expected '.' or a complete \
                     finite number"
                )
            })?;
            if !score_value.is_finite() {
                bail!("expected a finite BED score, found {score:?}")
            }
        }
    }
    let strand = match fields.get(5) {
        None => StrandRule::Both,
        Some(&"+") => StrandRule::Positive,
        Some(&"-") => StrandRule::Negative,
        Some(&".") => StrandRule::Both,
        Some(strand) => bail!(
            "invalid strand {strand:?}, expected exactly '+', '-', or '.'"
        ),
    };
    if start > end {
        bail!("start {start} is greater than end {end}")
    }

    Ok(GenomeRegion { chrom: chrom.to_string(), start, end, strand, name })
}

impl EntryLocalize {
    fn load_focus_regions(
        &self,
        sequence_lengths: &FxHashMap<String, u64>,
        index: &HtsTabixHandler<BedMethylLine>,
        multi_progress: &MultiProgress,
    ) -> anyhow::Result<Vec<GenomeRegion>> {
        let pb = multi_progress.add(get_ticker());
        pb.set_message("regions parsed");
        let mut data_lines = BufReader::new(File::open(&self.regions)?)
            .lines()
            .enumerate()
            .progress_with(pb)
            .filter_map(|(line_idx, result)| {
                let line_number = line_idx + 1;
                match result {
                    Ok(mut line) => {
                        if line.trim().is_empty()
                            || line.trim_start().starts_with('#')
                        {
                            None
                        } else {
                            let trimmed_len = line
                                .trim_end_matches(|c: char| {
                                    c.is_ascii_whitespace()
                                })
                                .len();
                            line.truncate(trimmed_len);
                            Some(Ok((line_number, line)))
                        }
                    }
                    Err(e) => Some(Err(anyhow!(
                        "failed to read regions BED file {} at line \
                         {line_number}: {e}",
                        self.regions.display()
                    ))),
                }
            });

        let Some(first_line) = data_lines.next().transpose()? else {
            bail!(
                "failed to inspect regions BED {}: no data lines",
                self.regions.display()
            )
        };
        let first_line_number = first_line.0;
        let delimiter = BedDelimiter::from_first_data_line(&first_line.1);
        let expected_fields = delimiter.split(&first_line.1).len();
        let regions = std::iter::once(Ok(first_line))
            .chain(data_lines)
            .map(|result| {
                let (line_number, line) = result?;
                if !delimiter.matches(&line) {
                    bail!(
                        "failed to parse regions BED file {} at line \
                         {line_number}: delimiter mode changed from \
                         {delimiter:?}, selected at line {first_line_number}",
                        self.regions.display()
                    )
                }
                let fields = delimiter.split(&line);
                let actual_fields = fields.len();
                if actual_fields != expected_fields {
                    bail!(
                        "failed to parse regions BED file {} at line \
                         {line_number}: expected {expected_fields} BED fields \
                         to match line {first_line_number}, found \
                         {actual_fields}",
                        self.regions.display()
                    )
                }
                parse_region_core_fields(&fields)
                    .map(|region| (line_number, region))
                    .map_err(|e| {
                        anyhow!(
                            "failed to parse regions BED file {} at line \
                             {line_number}: {e}",
                            self.regions.display()
                        )
                    })
            })
            .collect::<anyhow::Result<Vec<_>>>()?;

        let (regions, missing_from_sizes, missing_from_index) =
            regions.into_iter().try_fold(
                (Vec::new(), 0usize, 0usize),
                |(mut acc, no_size, no_index), (line_number, mut next)| {
                    let contig_length =
                        sequence_lengths.get(&next.chrom).copied();
                    if let Some(contig_length) = contig_length {
                        if next.start > contig_length
                            || next.end > contig_length
                        {
                            return Err(anyhow!(
                                "invalid region in {} at line {line_number}: \
                                 {}:{}-{} exceeds contig length \
                                 {contig_length}",
                                self.regions.display(),
                                next.chrom,
                                next.start,
                                next.end
                            ));
                        }
                    }
                    if contig_length.is_none() {
                        Ok((acc, no_size + 1, no_index))
                    } else if !index.has_contig(&next.chrom) {
                        Ok((acc, no_size, no_index + 1))
                    } else {
                        let contig_length = contig_length.unwrap();
                        let mp = next.midpoint();
                        let start = mp
                            .checked_sub(self.expand_window + 1)
                            .unwrap_or(0u64);
                        // safe because of conditional above
                        let end = std::cmp::min(
                            mp.saturating_add(self.expand_window),
                            contig_length,
                        );
                        next.start = start;
                        next.end = end;
                        acc.push(next);
                        Ok((acc, no_size, no_index))
                    }
                },
            )?;

        if missing_from_sizes > 0usize {
            debug!(
                "{missing_from_sizes} contigs in regions that are missing \
                 from contig sizes file"
            );
        }
        if missing_from_index > 0usize {
            debug!(
                "{missing_from_index} contigs in regions that are missing \
                 from tabix header"
            );
        }
        if regions.is_empty() {
            bail!("failed to find any valid regions")
        } else {
            debug!("{} valid regions parsed", regions.len());
        }

        info!("loaded {} regions", regions.len());
        Ok(regions)
    }

    pub fn run(&self) -> anyhow::Result<()> {
        let _ = init_logging(self.log_filepath.as_ref());

        // log out some of the settings...
        let min_cov = self.min_coverage;
        debug!("requiring min valid coverage of {min_cov}");
        if let Some(feature_strand) = self.stranded_features.as_ref() {
            info!("only using features from {feature_strand}");
        }
        if let Some(stranded) = self.stranded {
            info!(
                "only using features that are on {stranded:?} strand as the \
                 region of interest"
            );
        }

        let multi_progress = MultiProgress::new();
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()?;
        info!("loading sequence lengths from {:?}", &self.genome_sizes);

        let sequence_lengths = read_sequence_lengths_file(&self.genome_sizes)?
            .into_iter()
            .collect::<FxHashMap<String, u64>>();
        let tabix_index = HtsTabixHandler::from_path(&self.in_bedmethyl)
            .map(|x| Arc::new(x))?;
        let genome_regions = self.load_focus_regions(
            &sequence_lengths,
            &tabix_index,
            &multi_progress,
        )?;
        info!("loaded {} regions", genome_regions.len());

        let successes =
            multi_progress.add(get_master_progress_bar(genome_regions.len()));

        let stranded_features = self.stranded;
        let counts = pool.install(|| {
            genome_regions
                .into_par_iter()
                .progress_with(successes)
                .map(|gr| {
                    let region =
                        format!("{}:{}-{}", gr.chrom, gr.start, gr.end);
                    gr.into_localized_mod_counts(
                        &tabix_index,
                        self.stranded_features,
                        stranded_features,
                        self.io_threads,
                    )
                    .with_context(|| {
                        format!("failed to localize region {region}")
                    })
                })
                .try_fold(
                    || LocalizedModCounts::zero(),
                    |counts, next| -> anyhow::Result<_> {
                        Ok(counts.op(next?))
                    },
                )
                .try_reduce(|| LocalizedModCounts::zero(), |a, b| Ok(a.op(b)))
        })?;

        let table = counts.get_table(&multi_progress);
        let writer: Box<dyn Write> =
            if let Some(out_fp) = self.out_file.as_ref() {
                if self.force {
                    Box::new(BufWriter::new(File::create(out_fp)?))
                } else {
                    Box::new(BufWriter::new(File::create_new(out_fp)?))
                }
            } else {
                Box::new(BufWriter::new(stdout()))
            };
        let writer =
            csv::WriterBuilder::new().delimiter('\t' as u8).from_writer(writer);
        table.to_csv_writer(writer)?;

        if let Some(p) = self.chart_filepath.as_ref() {
            let fh = if self.force {
                File::create(p)?
            } else {
                File::create_new(p)?
            };
            let mut writer = BufWriter::new(fh);
            let blob = counts.get_plot(self.chart_name.as_ref())?;
            writer.write(blob.as_bytes())?;
        }

        multi_progress.clear()?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use std::fs::write;
    use std::path::PathBuf;

    use indicatif::MultiProgress;
    use rustc_hash::FxHashMap;

    use super::EntryLocalize;
    use crate::dmr::bedmethyl::BedMethylLine;
    use crate::tabix::HtsTabixHandler;
    use crate::util::{GenomeRegion, StrandRule};

    fn test_bedmethyl() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(
            "../tests/resources/\
            lung_00733-m_adjacent-normal_5mc-5hmc_chr20_cpg_pileup.bed.gz",
        )
    }

    fn load_regions_with_sequence_lengths(
        regions: PathBuf,
        sequence_lengths: FxHashMap<String, u64>,
    ) -> anyhow::Result<Vec<GenomeRegion>> {
        let entry = EntryLocalize {
            in_bedmethyl: test_bedmethyl(),
            regions,
            chart_filepath: None,
            chart_name: None,
            expand_window: 10,
            stranded: None,
            stranded_features: None,
            min_coverage: 1,
            genome_sizes: PathBuf::new(),
            out_file: None,
            log_filepath: None,
            threads: 1,
            io_threads: 1,
            force: false,
            batch_size_bp: 1,
        };
        let index =
            HtsTabixHandler::<BedMethylLine>::from_path(&entry.in_bedmethyl)?;

        entry.load_focus_regions(
            &sequence_lengths,
            &index,
            &MultiProgress::new(),
        )
    }

    fn load_regions(regions: PathBuf) -> anyhow::Result<Vec<GenomeRegion>> {
        load_regions_with_sequence_lengths(
            regions,
            FxHashMap::from_iter([("chr20".to_string(), 100_000_000)]),
        )
    }

    #[test]
    fn load_focus_regions_rejects_malformed_row_after_valid_row() {
        let temp_dir = tempfile::tempdir().unwrap();
        let regions = temp_dir.path().join("mixed-regions.bed");
        write(
            &regions,
            "# test regions\n\
             chr20\t9681998\t9681999\n\
             not-a-valid-bed-row\n\
             chr20\t9682013\t9682014\n",
        )
        .unwrap();

        let error = load_regions(regions.clone()).unwrap_err();
        let message = error.to_string();

        assert!(
            message.contains(regions.to_string_lossy().as_ref()),
            "missing region filepath in error: {message}"
        );
        assert!(message.contains("line 3"), "missing line number: {message}");
        assert!(
            message.contains("delimiter mode changed"),
            "missing parse reason: {message}"
        );
    }

    #[test]
    fn load_focus_regions_accepts_valid_rows() {
        let temp_dir = tempfile::tempdir().unwrap();
        let regions = temp_dir.path().join("valid-regions.bed");
        write(
            &regions,
            "# test regions\n\
             chr20\t9681998\t9681999\n\
             chr20\t9682013\t9682014\n",
        )
        .unwrap();

        assert_eq!(load_regions(regions).unwrap().len(), 2);
    }

    #[test]
    fn load_focus_regions_preserves_bed4_names_with_spaces() {
        let temp_dir = tempfile::tempdir().unwrap();
        let regions = temp_dir.path().join("bed4-regions.bed");
        write(&regions, "chr20\t9681998\t9681999\tregion name with spaces\n")
            .unwrap();

        let regions = load_regions(regions).unwrap();

        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].name.as_deref(), Some("region name with spaces"));
        assert_eq!(regions[0].strand, StrandRule::Both);
    }

    #[test]
    fn load_focus_regions_accepts_valid_bed5_score_and_opaque_extensions() {
        for (name, row, expected_strand) in [
            (
                "bed5",
                "chr20\t9681998\t9681999\tregion name\t500.5",
                StrandRule::Both,
            ),
            (
                "bed8",
                "chr20\t9681998\t9681999\tregion name\t.\t-\t\topaque",
                StrandRule::Negative,
            ),
        ] {
            let temp_dir = tempfile::tempdir().unwrap();
            let regions = temp_dir.path().join(format!("{name}-regions.bed"));
            write(&regions, format!("{row}\n")).unwrap();

            let parsed = load_regions(regions).unwrap();

            assert_eq!(parsed.len(), 1);
            assert_eq!(parsed[0].name.as_deref(), Some("region name"));
            assert_eq!(parsed[0].strand, expected_strand);
        }
    }

    #[test]
    fn load_focus_regions_ignores_spacing_and_preserves_exact_strands() {
        let temp_dir = tempfile::tempdir().unwrap();
        let regions = temp_dir.path().join("stranded-regions.bed");
        write(
            &regions,
            "  \n\
               # leading comment\n\
             chr20\t9681998\t9681999\tplus\t0\t+\n\
             \t\n\
               # interleaved comment\n\
             chr20\t9682013\t9682014\tminus\t0\t-\n\
             chr20\t9682030\t9682030\tzero-length\t0\t.\n",
        )
        .unwrap();

        let regions = load_regions(regions).unwrap();

        assert_eq!(regions.len(), 3);
        assert_eq!(regions[0].strand, StrandRule::Positive);
        assert_eq!(regions[1].strand, StrandRule::Negative);
        assert_eq!(regions[2].strand, StrandRule::Both);
    }

    #[test]
    fn load_focus_regions_rejects_mixed_bed_field_counts() {
        let temp_dir = tempfile::tempdir().unwrap();
        let regions = temp_dir.path().join("mixed-schema-regions.bed");
        write(
            &regions,
            "chr20\t9681998\t9681999\n\
             chr20\t9682013\t9682014\tminus\t0\t-\n",
        )
        .unwrap();

        let error = load_regions(regions.clone()).unwrap_err();
        let message = error.to_string();

        assert!(
            message.contains(regions.to_string_lossy().as_ref()),
            "missing region filepath in error: {message}"
        );
        assert!(message.contains("line 2"), "missing line number: {message}");
        assert!(
            message.contains("expected 3 BED fields")
                && message.contains("found 6"),
            "missing field-count reason: {message}"
        );
    }

    #[test]
    fn load_focus_regions_accepts_legacy_space_delimited_bed6() {
        let temp_dir = tempfile::tempdir().unwrap();
        let regions = temp_dir.path().join("space-delimited-regions.bed");
        write(&regions, "chr20 9681998 9681999 name 0 +\n").unwrap();

        let parsed = load_regions(regions).unwrap();

        assert_eq!(parsed.len(), 1);
        assert_eq!(parsed[0].name.as_deref(), Some("name"));
        assert_eq!(parsed[0].strand, StrandRule::Positive);
    }

    #[test]
    fn load_focus_regions_rejects_delimiter_mode_changes() {
        for (name, rows) in [
            (
                "space-then-tabs",
                "chr20 9681998 9681999 name 0 +\n\
                 chr20\t9682013\t9682014\tname\t0\t-\n",
            ),
            (
                "tabs-then-space",
                "chr20\t9681998\t9681999\tname\t0\t+\n\
                 chr20 9682013 9682014 name 0 -\n",
            ),
        ] {
            let temp_dir = tempfile::tempdir().unwrap();
            let regions = temp_dir.path().join(format!("{name}.bed"));
            write(&regions, rows).unwrap();

            let error = load_regions(regions.clone()).unwrap_err();
            let message = error.to_string();

            assert!(
                message.contains(regions.to_string_lossy().as_ref()),
                "missing region filepath in error: {message}"
            );
            assert!(message.contains("line 2"), "{message}");
            assert!(message.contains("delimiter mode"), "{message}");
        }
    }

    #[test]
    fn load_focus_regions_rejects_invalid_consumed_bed_fields() {
        for (name, row, reason) in [
            (
                "partial-end",
                "chr20\t9681998\t9681999x",
                "invalid end coordinate",
            ),
            (
                "invalid-score",
                "chr20\t9681998\t9681999\tname\t500x",
                "invalid BED score",
            ),
            (
                "non-finite-score",
                "chr20\t9681998\t9681999\tname\tNaN",
                "finite BED score",
            ),
            (
                "partial-strand",
                "chr20\t9681998\t9681999\tname\t0\t+garbage",
                "invalid strand",
            ),
        ] {
            let temp_dir = tempfile::tempdir().unwrap();
            let regions = temp_dir.path().join(format!("{name}-regions.bed"));
            write(&regions, format!("{row}\n")).unwrap();

            let error = load_regions(regions.clone()).unwrap_err();
            let message = error.to_string();

            assert!(
                message.contains(regions.to_string_lossy().as_ref()),
                "missing region filepath in error: {message}"
            );
            assert!(message.contains("line 1"), "{message}");
            assert!(message.contains(reason), "{message}");
        }
    }

    #[test]
    fn load_focus_regions_rejects_interior_empty_bed_name() {
        let temp_dir = tempfile::tempdir().unwrap();
        let regions = temp_dir.path().join("empty-name-regions.bed");
        write(&regions, "chr20\t9681998\t9681999\t\t500\n").unwrap();

        let error = load_regions(regions.clone()).unwrap_err();
        let message = error.to_string();

        assert!(
            message.contains(regions.to_string_lossy().as_ref()),
            "missing region filepath in error: {message}"
        );
        assert!(message.contains("line 1"), "{message}");
        assert!(message.contains("BED name must not be empty"), "{message}");
    }

    #[test]
    fn load_focus_regions_preserves_trailing_delimiter_compatibility() {
        for (name, row) in [
            ("tab", "chr20\t9681998\t9681999\t\t  "),
            ("space", "chr20 9681998 9681999   "),
        ] {
            let temp_dir = tempfile::tempdir().unwrap();
            let regions = temp_dir.path().join(format!("{name}-trailing.bed"));
            write(&regions, format!("{row}\n")).unwrap();

            let parsed = load_regions(regions).unwrap();

            assert_eq!(parsed.len(), 1);
            assert_eq!(parsed[0].name, None);
            assert_eq!(parsed[0].strand, StrandRule::Both);
        }
    }

    #[test]
    fn load_focus_regions_rejects_missing_strand_after_bed6_schema() {
        let temp_dir = tempfile::tempdir().unwrap();
        let regions = temp_dir.path().join("missing-strand-regions.bed");
        write(
            &regions,
            "chr20\t9681998\t9681999\tplus\t0\t+\n\
             chr20\t9682013\t9682014\tmissing\t0\t\n",
        )
        .unwrap();

        let error = load_regions(regions.clone()).unwrap_err();
        let message = error.to_string();

        assert!(
            message.contains(regions.to_string_lossy().as_ref()),
            "missing region filepath in error: {message}"
        );
        assert!(message.contains("line 2"), "{message}");
        assert!(
            message.contains("expected 6 BED fields")
                && message.contains("found 5"),
            "{message}"
        );
    }

    #[test]
    fn load_focus_regions_rejects_headers_and_browser_directives() {
        for (name, first_line) in [
            ("header", "chrom\tstart\tend"),
            ("track", "track name=test"),
            ("browser", "browser position chr20:9681998-9682014"),
        ] {
            let temp_dir = tempfile::tempdir().unwrap();
            let regions = temp_dir.path().join(format!("{name}-regions.bed"));
            write(
                &regions,
                format!(
                    "{first_line}\n\
                     chr20\t9681998\t9681999\n"
                ),
            )
            .unwrap();

            let error = load_regions(regions.clone()).unwrap_err();
            let message = error.to_string();

            assert!(
                message.contains(regions.to_string_lossy().as_ref()),
                "missing region filepath in error: {message}"
            );
            assert!(message.contains("line 1"), "{message}");
            assert!(message.contains("failed to parse"), "{message}");
        }
    }

    #[test]
    fn load_focus_regions_rejects_invalid_coordinate_ranges() {
        for (name, row, reason) in [
            (
                "reversed",
                "chr20\t9682014\t9682013",
                "start 9682014 is greater than end 9682013",
            ),
            (
                "past-contig-end",
                "chr20\t100000001\t100000001",
                "exceeds contig length 100000000",
            ),
        ] {
            let temp_dir = tempfile::tempdir().unwrap();
            let regions = temp_dir.path().join(format!("{name}-regions.bed"));
            write(&regions, format!("{row}\n")).unwrap();

            let error = load_regions(regions.clone()).unwrap_err();
            let message = error.to_string();

            assert!(
                message.contains(regions.to_string_lossy().as_ref()),
                "missing region filepath in error: {message}"
            );
            assert!(message.contains("line 1"), "{message}");
            assert!(message.contains(reason), "{message}");
        }
    }

    #[test]
    fn load_focus_regions_accepts_zero_length_at_contig_end() {
        let temp_dir = tempfile::tempdir().unwrap();
        let regions = temp_dir.path().join("boundary-region.bed");
        write(&regions, "chr20\t100000000\t100000000\n").unwrap();

        assert_eq!(load_regions(regions).unwrap().len(), 1);
    }

    #[test]
    fn load_focus_regions_validates_coordinates_before_missing_contigs() {
        for (name, row, sequence_lengths, reason) in [
            (
                "reversed-unknown-contig",
                "unknown\t10\t5",
                FxHashMap::default(),
                "start 10 is greater than end 5",
            ),
            (
                "past-known-unindexed-contig",
                "known-unindexed\t101\t101",
                FxHashMap::from_iter([("known-unindexed".to_string(), 100)]),
                "exceeds contig length 100",
            ),
        ] {
            let temp_dir = tempfile::tempdir().unwrap();
            let regions = temp_dir.path().join(format!("{name}.bed"));
            write(&regions, format!("{row}\n")).unwrap();

            let error = load_regions_with_sequence_lengths(
                regions.clone(),
                sequence_lengths,
            )
            .unwrap_err();
            let message = error.to_string();

            assert!(
                message.contains(regions.to_string_lossy().as_ref()),
                "missing region filepath in error: {message}"
            );
            assert!(message.contains("line 1"), "{message}");
            assert!(message.contains(reason), "{message}");
        }
    }
}
