use std::cmp::Ordering;
use std::fs::{self, File};
use std::io::{Read as IoRead, Seek, SeekFrom};
#[cfg(unix)]
use std::os::unix::fs::{MetadataExt, PermissionsExt};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{anyhow, bail, Context};
use clap::Args;
use crossbeam_channel::unbounded;
use derive_new::new;
use indicatif::{MultiProgress, ProgressBar};
use log::{debug, info, warn};
use memchr::memmem::Finder;
use rust_htslib::bam::record::{Aux, AuxArray};
use rust_htslib::bam::{self, Read};
use rustc_hash::FxHashMap;
use tempfile::{Builder as TempFileBuilder, NamedTempFile, TempDir};

use modkit_logging::init_logging;

use crate::mod_bam::{
    format_mm_ml_tag, BaseModProbs, DeltaListConverter, ModBaseInfo,
    SeqPosBaseModProbs, ML_TAGS, MM_TAGS, MN_TAG,
};
use crate::ordered_scheduler::{run_ordered_scheduler, OrderedWorker};
use crate::util::{
    get_forward_sequence_str, get_query_name_string, get_ticker,
    record_is_not_primary,
};

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct RepairTags {
    /// Donor modBAM with original MM/ML tags. Must be sorted by read name.
    #[arg(long, short = 'd', alias = "donor")]
    donor_bam: PathBuf,
    /// Acceptor modBAM with reads to have MM/ML base modification data
    /// projected on to. Must be sorted by read name.
    #[arg(long, short = 'a', alias = "acceptor")]
    acceptor_bam: PathBuf,
    /// output modBAM location.
    #[arg(long, short = 'o', alias = "output")]
    output_bam: PathBuf,
    /// File to write logs to, it is recommended to use this option as some
    /// reads may be rejected and logged here.
    #[arg(long)]
    log_filepath: Option<PathBuf>,
    /// The number of threads to use.
    #[arg(long, short = 't', default_value_t = 4)]
    threads: usize,
}

const REPAIR_TEMP_PREFIX: &str = ".modkit-repair-";
const BGZF_EOF_BLOCK: [u8; 28] = [
    31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, 66, 67, 2, 0, 27, 0, 3, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
];

trait RepairSink {
    fn write(&mut self, record: &bam::Record) -> anyhow::Result<()>;
    fn finish(self, expected_records: usize) -> anyhow::Result<()>;
}

trait StagedBamValidator {
    fn validate(
        &self,
        file: &mut File,
        path: &Path,
        expected_identity: StagedFileIdentity,
        expected_records: usize,
    ) -> anyhow::Result<()>;
}

struct CompleteBamValidator;

impl StagedBamValidator for CompleteBamValidator {
    fn validate(
        &self,
        file: &mut File,
        path: &Path,
        expected_identity: StagedFileIdentity,
        expected_records: usize,
    ) -> anyhow::Result<()> {
        validate_complete_bam_file(
            file,
            path,
            expected_identity,
            expected_records,
        )
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct StagedFileIdentity {
    #[cfg(unix)]
    device: u64,
    #[cfg(unix)]
    inode: u64,
}

impl StagedFileIdentity {
    fn from_metadata(metadata: &fs::Metadata) -> Self {
        #[cfg(unix)]
        {
            Self { device: metadata.dev(), inode: metadata.ino() }
        }
        #[cfg(not(unix))]
        {
            let _ = metadata;
            Self {}
        }
    }
}

fn require_staged_path_identity(
    file: &File,
    path: &Path,
    expected: StagedFileIdentity,
) -> anyhow::Result<()> {
    let retained_metadata = file.metadata().with_context(|| {
        format!(
            "failed to inspect retained staged repair BAM descriptor for {}",
            path.display()
        )
    })?;
    if !retained_metadata.file_type().is_file() {
        bail!(
            "retained staged repair BAM descriptor for {} is not a regular file",
            path.display()
        );
    }
    let path_metadata = fs::symlink_metadata(path).with_context(|| {
        format!(
            "staged repair BAM pathname {} no longer identifies the retained file",
            path.display()
        )
    })?;
    if !path_metadata.file_type().is_file() {
        bail!(
            "staged repair BAM pathname {} no longer identifies the retained regular file",
            path.display()
        );
    }
    #[cfg(unix)]
    {
        let retained_identity =
            StagedFileIdentity::from_metadata(&retained_metadata);
        let path_identity = StagedFileIdentity::from_metadata(&path_metadata);
        if retained_identity != expected || path_identity != expected {
            bail!(
                "staged repair BAM pathname {} no longer identifies the retained file \
                 (expected dev {} ino {}, retained dev {} ino {}, path dev {} ino {})",
                path.display(),
                expected.device,
                expected.inode,
                retained_identity.device,
                retained_identity.inode,
                path_identity.device,
                path_identity.inode
            );
        }
    }
    #[cfg(not(unix))]
    let _ = expected;
    Ok(())
}

struct StagedBamOutput<V> {
    // Keep the writer and temporary file before the directory so early-drop
    // closes htslib, removes the staged file, then removes the private dir.
    writer: Option<bam::Writer>,
    temp_file: NamedTempFile,
    _staging_dir: TempDir,
    destination: PathBuf,
    destination_permissions: Option<fs::Permissions>,
    staged_identity: StagedFileIdentity,
    validator: V,
}

fn existing_destination_permissions(
    destination: &Path,
) -> anyhow::Result<Option<fs::Permissions>> {
    match fs::symlink_metadata(destination) {
        Ok(metadata) if metadata.file_type().is_symlink() => bail!(
            "repair output destination {} must not be a symbolic link",
            destination.display()
        ),
        Ok(metadata) if !metadata.file_type().is_file() => bail!(
            "repair output destination {} must be a regular file when it already exists",
            destination.display()
        ),
        Ok(metadata) => Ok(Some(metadata.permissions())),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(None),
        Err(error) => Err(error).with_context(|| {
            format!(
                "failed to inspect repair output destination {}",
                destination.display()
            )
        }),
    }
}

impl StagedBamOutput<CompleteBamValidator> {
    fn new(destination: &Path, header: &bam::Header) -> anyhow::Result<Self> {
        Self::with_validator(destination, header, CompleteBamValidator)
    }
}

impl<V> StagedBamOutput<V> {
    fn with_validator(
        destination: &Path,
        header: &bam::Header,
        validator: V,
    ) -> anyhow::Result<Self> {
        let destination_permissions =
            existing_destination_permissions(destination)?;
        let destination_dir = destination
            .parent()
            .filter(|parent| !parent.as_os_str().is_empty())
            .unwrap_or_else(|| Path::new("."));

        let mut staging_dir_builder = TempFileBuilder::new();
        staging_dir_builder.prefix(REPAIR_TEMP_PREFIX);
        #[cfg(unix)]
        staging_dir_builder.permissions(fs::Permissions::from_mode(0o700));
        let staging_dir = staging_dir_builder
            .tempdir_in(destination_dir)
            .with_context(|| {
                format!(
                    "failed to create private repair staging directory in {}",
                    destination_dir.display()
                )
            })?;
        #[cfg(unix)]
        {
            // Builder applies the process umask. A restrictive umask may
            // remove owner bits, but group/other access must never be added.
            let staging_mode = fs::symlink_metadata(staging_dir.path())
                .with_context(|| {
                    format!(
                        "failed to inspect repair staging directory {}",
                        staging_dir.path().display()
                    )
                })?
                .permissions()
                .mode()
                & 0o7777;
            if staging_mode & 0o077 != 0 {
                bail!(
                    "repair staging directory {} has non-owner mode bits {staging_mode:o}",
                    staging_dir.path().display()
                );
            }
        }

        let mut temp_file_builder = TempFileBuilder::new();
        temp_file_builder.prefix("staged-").suffix(".bam");
        #[cfg(unix)]
        {
            // Match an ordinary File::create rather than tempfile's 0600
            // default; Builder applies the process umask to this mode.
            temp_file_builder.permissions(fs::Permissions::from_mode(0o666));
        }
        let temp_file = temp_file_builder
            .tempfile_in(staging_dir.path())
            .with_context(|| {
                format!(
                    "failed to create staged repair BAM in private directory {}",
                    staging_dir.path().display()
                )
            })?;
        let staged_identity = StagedFileIdentity::from_metadata(
            &temp_file.as_file().metadata().with_context(|| {
                format!(
                    "failed to inspect retained staged repair BAM {}",
                    temp_file.path().display()
                )
            })?,
        );
        let writer =
            bam::Writer::from_path(temp_file.path(), header, bam::Format::Bam)
                .with_context(|| {
                    format!(
                        "failed to open staged repair BAM {}",
                        temp_file.path().display()
                    )
                })?;
        require_staged_path_identity(
            temp_file.as_file(),
            temp_file.path(),
            staged_identity,
        )
        .context("staged repair BAM changed while opening its htslib writer")?;
        Ok(Self {
            writer: Some(writer),
            temp_file,
            _staging_dir: staging_dir,
            destination: destination.to_path_buf(),
            destination_permissions,
            staged_identity,
            validator,
        })
    }
}

impl<V: StagedBamValidator> RepairSink for StagedBamOutput<V> {
    fn write(&mut self, record: &bam::Record) -> anyhow::Result<()> {
        self.writer
            .as_mut()
            .ok_or_else(|| {
                anyhow!("staged repair BAM writer is already closed")
            })?
            .write(record)
            .context("failed to write staged repair BAM record")
    }

    fn finish(mut self, expected_records: usize) -> anyhow::Result<()> {
        let writer = self.writer.take().ok_or_else(|| {
            anyhow!("staged repair BAM writer is already closed")
        })?;
        // rust-htslib 0.46 does not expose hts_close's return status. Drop the
        // writer, then independently require a readable BAM, canonical EOF,
        // and the expected record count before making it visible.
        drop(writer);
        let staged_path = self.temp_file.path().to_path_buf();
        require_staged_path_identity(
            self.temp_file.as_file(),
            &staged_path,
            self.staged_identity,
        )
        .context("staged repair BAM changed after closing its htslib writer")?;
        self.validator
            .validate(
                self.temp_file.as_file_mut(),
                &staged_path,
                self.staged_identity,
                expected_records,
            )
            .context("failed to validate staged repair BAM")?;
        require_staged_path_identity(
            self.temp_file.as_file(),
            &staged_path,
            self.staged_identity,
        )
        .context("staged repair BAM changed during validation")?;
        if let Some(permissions) = self.destination_permissions.take() {
            // Restore after all writes because a write may clear special mode
            // bits. Any failure still precedes the atomic replacement.
            self.temp_file
                .as_file()
                .set_permissions(permissions)
                .with_context(|| {
                    format!(
                        "failed to preserve permissions for existing repair output {}",
                        self.destination.display()
                    )
                })?;
        }
        self.temp_file.as_file().sync_all().with_context(|| {
            format!(
                "failed to sync staged repair BAM {}",
                staged_path.display()
            )
        })?;
        // NamedTempFile::persist still names the source path. These checks
        // reject observed substitutions, although a hostile writer of a
        // non-sticky parent directory can still race the check and rename.
        require_staged_path_identity(
            self.temp_file.as_file(),
            &staged_path,
            self.staged_identity,
        )
        .context("staged repair BAM changed before atomic replacement")?;
        self.temp_file.persist(&self.destination).map_err(|error| {
            anyhow!(
                "failed to atomically replace repair output {}: {error}",
                self.destination.display()
            )
        })?;
        Ok(())
    }
}

#[cfg(test)]
fn validate_complete_bam(
    path: &Path,
    expected_records: usize,
) -> anyhow::Result<()> {
    let mut file = File::open(path).with_context(|| {
        format!("failed to open staged repair BAM {}", path.display())
    })?;
    let expected_identity = StagedFileIdentity::from_metadata(
        &file.metadata().with_context(|| {
            format!("failed to inspect staged repair BAM {}", path.display())
        })?,
    );
    validate_complete_bam_file(
        &mut file,
        path,
        expected_identity,
        expected_records,
    )
}

fn validate_complete_bam_file(
    file: &mut File,
    path: &Path,
    expected_identity: StagedFileIdentity,
    expected_records: usize,
) -> anyhow::Result<()> {
    let file_len = file
        .metadata()
        .with_context(|| {
            format!("failed to inspect staged repair BAM {}", path.display())
        })?
        .len();
    if file_len < BGZF_EOF_BLOCK.len() as u64 {
        bail!(
            "staged repair BAM {} is too short to contain a BGZF EOF block",
            path.display()
        )
    }
    file.seek(SeekFrom::End(-(BGZF_EOF_BLOCK.len() as i64))).with_context(
        || format!("failed to seek staged repair BAM {}", path.display()),
    )?;
    let mut eof = [0u8; BGZF_EOF_BLOCK.len()];
    file.read_exact(&mut eof).with_context(|| {
        format!("failed to read staged repair BAM EOF at {}", path.display())
    })?;
    if eof != BGZF_EOF_BLOCK {
        bail!(
            "staged repair BAM {} is missing the canonical BGZF EOF block",
            path.display()
        )
    }

    require_staged_path_identity(file, path, expected_identity)
        .context("staged repair BAM changed before full BAM decoding")?;
    let mut reader = bam::Reader::from_path(path).with_context(|| {
        format!("failed to reopen staged repair BAM {}", path.display())
    })?;
    let mut observed_records = 0usize;
    for (index, record) in reader.records().enumerate() {
        record.with_context(|| {
            format!(
                "failed to decode staged repair BAM record {} at {}",
                index + 1,
                path.display()
            )
        })?;
        observed_records += 1;
    }
    drop(reader);
    require_staged_path_identity(file, path, expected_identity)
        .context("staged repair BAM changed during full BAM decoding")?;
    if observed_records != expected_records {
        bail!(
            "staged repair BAM {} contains {observed_records} records, \
             expected {expected_records}",
            path.display()
        )
    }
    Ok(())
}

fn open_bam_reader(
    path: &Path,
    input_label: &str,
    threads: usize,
) -> anyhow::Result<bam::Reader> {
    let mut reader = bam::Reader::from_path(path).with_context(|| {
        format!("failed to open {input_label} BAM {}", path.display())
    })?;
    reader.set_threads(threads).with_context(|| {
        format!(
            "failed to configure {input_label} BAM reader {}",
            path.display()
        )
    })?;
    Ok(reader)
}

impl RepairTags {
    pub fn run(&self) -> anyhow::Result<()> {
        let _handle = init_logging(self.log_filepath.as_ref());
        self.run_with_output_factory(|header| {
            StagedBamOutput::new(&self.output_bam, header)
        })
    }

    fn run_with_output_factory<S, F>(
        &self,
        make_output: F,
    ) -> anyhow::Result<()>
    where
        S: RepairSink,
        F: FnOnce(&bam::Header) -> anyhow::Result<S>,
    {
        let reader_threads = {
            let half = self.threads / 2;
            std::cmp::min(half, 16)
        };
        let threads_per_reader = std::cmp::max(reader_threads / 2, 1);
        let pool_threads = self.threads.saturating_sub(reader_threads).max(1);
        debug!(
            "assigning {threads_per_reader} to each reader and using \
             {pool_threads} to process records"
        );

        let donor_records =
            open_bam_reader(&self.donor_bam, "donor", threads_per_reader)?;
        let acceptor_records = open_bam_reader(
            &self.acceptor_bam,
            "acceptor",
            threads_per_reader,
        )?;
        let donor_order = query_name_order(donor_records.header(), "donor")?;
        let acceptor_order =
            query_name_order(acceptor_records.header(), "acceptor")?;
        let merge_order = match donor_order.common_merge_order(acceptor_order) {
            Some(order) => order,
            None => bail!(
                "donor and acceptor BAMs use incompatible query-name sorting: \
                 donor is {}, acceptor is {}",
                donor_order.label(),
                acceptor_order.label()
            ),
        };
        let header = bam::Header::from_template(acceptor_records.header());
        let mut output = make_output(&header)?;
        info!(
            "repairing records in {} with base modification information in {}",
            &self.acceptor_bam.to_str().unwrap_or_else(|| "??"),
            &self.donor_bam.to_str().unwrap_or_else(|| "??")
        );

        // pb stuff
        let master_progress = MultiProgress::new();
        let donor_ticker = master_progress.add(get_ticker());
        donor_ticker.set_message("~donor records processed");
        let acceptor_ticker = master_progress.add(get_ticker());
        acceptor_ticker.set_message("~acceptor records processed");
        let repaired_ticker = master_progress.add(get_ticker());
        repaired_ticker.set_message("~records repaired");
        let written_ticker =
            master_progress.add(get_ticker()).with_message("~records written");

        let pair_iter = ZipRecordsIter::new(
            donor_records,
            acceptor_records,
            donor_ticker,
            acceptor_ticker,
            merge_order,
        );
        let workers =
            (0..pool_threads).map(|_| RepairWorker).collect::<Vec<_>>();
        let mut n_repaired = 0usize;
        let mut n_failed = 0usize;
        run_repair_scheduler(pair_iter, workers, 1000, |res| {
            repaired_ticker.inc(1);
            match res {
                Ok(record) => {
                    output.write(&record).with_context(|| {
                        format!(
                            "failed to write repaired record {}",
                            String::from_utf8_lossy(record.qname())
                        )
                    })?;
                    written_ticker.inc(1);
                    n_repaired += 1;
                }
                Err(e) => {
                    debug!("record failed to be repaired: {}", e.to_string());
                    n_failed += 1;
                }
            }
            Ok(())
        })?;

        output.finish(n_repaired)?;
        info!("finished, repaired {n_repaired} records, {n_failed} failed.");
        Ok(())
    }
}

struct RepairSlot {
    outcome: Option<anyhow::Result<bam::Record>>,
}

enum RepairJob {
    Matched(RecordPair),
    MissingDonor { read_name: Vec<u8> },
}

struct RepairWorker;

impl OrderedWorker<RepairJob, RepairSlot> for RepairWorker {
    fn process(
        &mut self,
        job: RepairJob,
        mut slot: RepairSlot,
    ) -> anyhow::Result<RepairSlot> {
        debug_assert!(slot.outcome.is_none());
        slot.outcome = Some(match job {
            RepairJob::Matched(record_pair) => repair_record_pair(record_pair),
            RepairJob::MissingDonor { read_name } => Err(anyhow!(
                "record {} failed, no primary donor record",
                String::from_utf8_lossy(&read_name)
            )),
        });
        Ok(slot)
    }
}

fn run_repair_scheduler<I, W, C>(
    feeder: I,
    workers: Vec<W>,
    output_queue_size: usize,
    mut consume: C,
) -> anyhow::Result<()>
where
    I: Iterator<Item = anyhow::Result<RepairJob>> + Send + 'static,
    W: OrderedWorker<RepairJob, RepairSlot> + 'static,
    C: FnMut(anyhow::Result<bam::Record>) -> anyhow::Result<()>,
{
    let (empty_sender, empty_buffers) = unbounded();
    for _ in 0..(workers.len() * 2) {
        empty_sender
            .send(RepairSlot { outcome: None })
            .expect("unbounded repair buffer channel should be connected");
    }
    let recycle = empty_sender.clone();
    run_ordered_scheduler(
        "repair",
        feeder,
        workers,
        empty_buffers,
        output_queue_size,
        move |mut slot| {
            let outcome = slot.outcome.take().ok_or_else(|| {
                anyhow!("repair worker returned an empty slot")
            })?;
            let consumed = consume(outcome);
            // The source drops its receiver after it has fed the last job, so
            // recycling the final ordered slots may legitimately disconnect.
            let _ = recycle.send(slot);
            consumed
        },
    )
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum QueryNameOrder {
    /// `SS:queryname:natural`, matching samtools' numeric-run comparator.
    Natural,
    /// SO-only bundled inputs were sorted by an older natural comparator that
    /// orders numerically tied runs by leading-zero count.
    LegacyNatural,
    Lexicographical,
}

impl QueryNameOrder {
    fn common_merge_order(self, other: Self) -> Option<Self> {
        match (self, other) {
            (Self::Natural, Self::Natural) => Some(Self::Natural),
            (Self::LegacyNatural, Self::LegacyNatural) => {
                Some(Self::LegacyNatural)
            }
            (Self::Lexicographical, Self::Lexicographical) => {
                Some(Self::Lexicographical)
            }
            _ => None,
        }
    }

    fn compare(self, left: &[u8], right: &[u8]) -> Ordering {
        match self {
            Self::Natural => natural_query_name_cmp(left, right),
            Self::LegacyNatural => legacy_natural_query_name_cmp(left, right),
            Self::Lexicographical => left.cmp(right),
        }
    }

    fn label(self) -> &'static str {
        match self {
            Self::Natural => "natural",
            Self::LegacyNatural => "natural (SO-only legacy fallback)",
            Self::Lexicographical => "lexicographical",
        }
    }
}

fn natural_query_name_cmp(left: &[u8], right: &[u8]) -> Ordering {
    natural_query_name_cmp_impl(left, right, false)
}

fn legacy_natural_query_name_cmp(left: &[u8], right: &[u8]) -> Ordering {
    natural_query_name_cmp_impl(left, right, true)
}

fn natural_query_name_cmp_impl(
    left: &[u8],
    right: &[u8],
    tie_break_leading_zeros: bool,
) -> Ordering {
    let mut left_idx = 0usize;
    let mut right_idx = 0usize;
    while left_idx < left.len() && right_idx < right.len() {
        let left_byte = left[left_idx];
        let right_byte = right[right_idx];
        if left_byte.is_ascii_digit() && right_byte.is_ascii_digit() {
            let left_end = left[left_idx..]
                .iter()
                .position(|byte| !byte.is_ascii_digit())
                .map_or(left.len(), |offset| left_idx + offset);
            let right_end = right[right_idx..]
                .iter()
                .position(|byte| !byte.is_ascii_digit())
                .map_or(right.len(), |offset| right_idx + offset);
            let left_digits = &left[left_idx..left_end];
            let right_digits = &right[right_idx..right_end];
            let left_significant = left_digits
                .iter()
                .position(|byte| *byte != b'0')
                .map_or(&left_digits[left_digits.len()..], |idx| {
                    &left_digits[idx..]
                });
            let right_significant = right_digits
                .iter()
                .position(|byte| *byte != b'0')
                .map_or(&right_digits[right_digits.len()..], |idx| {
                    &right_digits[idx..]
                });
            match left_significant.len().cmp(&right_significant.len()) {
                Ordering::Equal => {}
                ordering => return ordering,
            }
            match left_significant.cmp(right_significant) {
                Ordering::Equal => {}
                ordering => return ordering,
            }
            if tie_break_leading_zeros {
                let left_zero_count =
                    left_digits.len() - left_significant.len();
                let right_zero_count =
                    right_digits.len() - right_significant.len();
                match right_zero_count.cmp(&left_zero_count) {
                    Ordering::Equal => {}
                    // Older samtools natural sorting placed more leading
                    // zeros before fewer leading zeros.
                    ordering => return ordering,
                }
            }
            // Keep numerically equal runs comparator-equivalent. The merge
            // groups these names and still matches their exact QNAME bytes.
            left_idx = left_end;
            right_idx = right_end;
        } else {
            match left_byte.cmp(&right_byte) {
                Ordering::Equal => {
                    left_idx += 1;
                    right_idx += 1;
                }
                ordering => return ordering,
            }
        }
    }
    left.len()
        .saturating_sub(left_idx)
        .cmp(&right.len().saturating_sub(right_idx))
}

fn header_tag<'a>(line: &'a [u8], tag: &[u8; 2]) -> Option<&'a [u8]> {
    line.split(|byte| *byte == b'\t').skip(1).find_map(|field| {
        field.strip_prefix(tag).and_then(|value| value.strip_prefix(b":"))
    })
}

fn query_name_order(
    header: &bam::HeaderView,
    input_label: &str,
) -> anyhow::Result<QueryNameOrder> {
    let hd = header
        .as_bytes()
        .split(|byte| *byte == b'\n')
        .find(|line| line.starts_with(b"@HD\t"))
        .ok_or_else(|| anyhow!("{input_label} BAM has no @HD header record"))?;
    let sort_order = header_tag(hd, b"SO")
        .ok_or_else(|| anyhow!("{input_label} BAM @HD record has no SO tag"))?;
    if sort_order != b"queryname" {
        bail!(
            "{input_label} BAM must be query-name sorted, found SO:{}",
            String::from_utf8_lossy(sort_order)
        )
    }
    let Some(sub_sort) = header_tag(hd, b"SS") else {
        warn!(
            "{input_label} BAM declares SO:queryname without SS; assuming \
             natural query-name order with the legacy leading-zero tie-break"
        );
        return Ok(QueryNameOrder::LegacyNatural);
    };
    let mut parts = sub_sort.split(|byte| *byte == b':');
    match (parts.next(), parts.next(), parts.next()) {
        (Some(b"queryname"), Some(b"natural"), None) => {
            Ok(QueryNameOrder::Natural)
        }
        (Some(b"queryname"), Some(b"lexicographical"), None) => {
            Ok(QueryNameOrder::Lexicographical)
        }
        _ => bail!(
            "{input_label} BAM uses unsupported query-name sub-sort SS:{}",
            String::from_utf8_lossy(sub_sort)
        ),
    }
}

#[derive(new)]
struct RecordPair {
    donor: Arc<bam::Record>,
    acceptor: bam::Record,
}

struct DonorGroup {
    representative: Vec<u8>,
    by_exact_name: FxHashMap<Vec<u8>, Arc<bam::Record>>,
}

struct ZipRecordsIter<D: Read, A: Read> {
    donor_records: D,
    acceptor_records: A,
    pending_donor: Option<bam::Record>,
    donor_group: Option<DonorGroup>,
    cur_acceptor_record: Option<bam::Record>,
    last_donor_name: Option<Vec<u8>>,
    last_acceptor_name: Option<Vec<u8>>,
    donor_records_read: usize,
    acceptor_records_read: usize,
    donor_ticker: ProgressBar,
    acceptor_ticker: ProgressBar,
    order: QueryNameOrder,
    terminated: bool,
}

impl<D: Read, A: Read> ZipRecordsIter<D, A> {
    fn new(
        donor: D,
        acceptor: A,
        donor_ticker: ProgressBar,
        acceptor_ticker: ProgressBar,
        order: QueryNameOrder,
    ) -> Self {
        Self {
            donor_records: donor,
            acceptor_records: acceptor,
            pending_donor: None,
            donor_group: None,
            cur_acceptor_record: None,
            last_donor_name: None,
            last_acceptor_name: None,
            donor_records_read: 0,
            acceptor_records_read: 0,
            donor_ticker,
            acceptor_ticker,
            order,
            terminated: false,
        }
    }

    fn validate_monotonic(
        order: QueryNameOrder,
        previous: Option<&[u8]>,
        current: &[u8],
        input_label: &str,
    ) -> anyhow::Result<()> {
        if let Some(previous) = previous {
            if order.compare(previous, current) == Ordering::Greater {
                bail!(
                    "{input_label} BAM is not {} query-name sorted: {} \
                     appears after {}",
                    order.label(),
                    String::from_utf8_lossy(current),
                    String::from_utf8_lossy(previous)
                )
            }
        }
        Ok(())
    }

    fn read_primary_donor(&mut self) -> anyhow::Result<Option<bam::Record>> {
        loop {
            let Some(record) = get_next_record(
                &mut self.donor_records,
                "donor",
                &mut self.donor_records_read,
            )?
            else {
                return Ok(None);
            };
            let name = record.qname().to_vec();
            Self::validate_monotonic(
                self.order,
                self.last_donor_name.as_deref(),
                &name,
                "donor",
            )?;
            self.last_donor_name = Some(name);
            self.donor_ticker.inc(1);
            if !record_is_not_primary(&record) {
                return Ok(Some(record));
            }
        }
    }

    fn load_donor_group(&mut self) -> anyhow::Result<Option<DonorGroup>> {
        let first = match self.pending_donor.take() {
            Some(record) => record,
            None => match self.read_primary_donor()? {
                Some(record) => record,
                None => return Ok(None),
            },
        };
        let representative = first.qname().to_vec();
        let mut by_exact_name = FxHashMap::default();
        by_exact_name.insert(representative.clone(), Arc::new(first));
        while let Some(record) = self.read_primary_donor()? {
            let name = record.qname().to_vec();
            match self.order.compare(&representative, &name) {
                Ordering::Equal => {
                    if by_exact_name
                        .insert(name.clone(), Arc::new(record))
                        .is_some()
                    {
                        bail!(
                            "donor BAM has multiple primary records for \
                             query name {}",
                            String::from_utf8_lossy(&name)
                        )
                    }
                }
                Ordering::Less => {
                    self.pending_donor = Some(record);
                    break;
                }
                Ordering::Greater => unreachable!(
                    "donor monotonicity is checked before grouping"
                ),
            }
        }
        Ok(Some(DonorGroup { representative, by_exact_name }))
    }

    fn load_acceptor(&mut self) -> anyhow::Result<bool> {
        if self.cur_acceptor_record.is_some() {
            return Ok(true);
        }
        let Some(record) = get_next_record(
            &mut self.acceptor_records,
            "acceptor",
            &mut self.acceptor_records_read,
        )?
        else {
            return Ok(false);
        };
        let name = record.qname().to_vec();
        Self::validate_monotonic(
            self.order,
            self.last_acceptor_name.as_deref(),
            &name,
            "acceptor",
        )?;
        self.last_acceptor_name = Some(name);
        self.cur_acceptor_record = Some(record);
        Ok(true)
    }

    fn consume_acceptor(&mut self) -> bam::Record {
        self.acceptor_ticker.inc(1);
        self.cur_acceptor_record
            .take()
            .expect("acceptor is loaded before it is consumed")
    }

    fn next_job(&mut self) -> anyhow::Result<Option<RepairJob>> {
        if !self.load_acceptor()? {
            debug!("exhausted acceptor BAM reader, finished.");
            return Ok(None);
        }
        let acceptor_name = self
            .cur_acceptor_record
            .as_ref()
            .expect("acceptor was loaded")
            .qname()
            .to_vec();
        loop {
            if self.donor_group.is_none() {
                self.donor_group = self.load_donor_group()?;
            }
            let Some(group) = self.donor_group.as_ref() else {
                self.consume_acceptor();
                return Ok(Some(RepairJob::MissingDonor {
                    read_name: acceptor_name,
                }));
            };
            match self.order.compare(&group.representative, &acceptor_name) {
                Ordering::Less => {
                    self.donor_group = None;
                }
                Ordering::Greater => {
                    self.consume_acceptor();
                    return Ok(Some(RepairJob::MissingDonor {
                        read_name: acceptor_name,
                    }));
                }
                Ordering::Equal => {
                    let donor =
                        group.by_exact_name.get(&acceptor_name).cloned();
                    let acceptor = self.consume_acceptor();
                    return Ok(Some(match donor {
                        Some(donor) => {
                            RepairJob::Matched(RecordPair::new(donor, acceptor))
                        }
                        None => {
                            RepairJob::MissingDonor { read_name: acceptor_name }
                        }
                    }));
                }
            }
        }
    }
}

fn get_next_record<T: Read>(
    records: &mut T,
    input_label: &str,
    records_read: &mut usize,
) -> anyhow::Result<Option<bam::Record>> {
    let mut record = bam::Record::new();
    match records.read(&mut record) {
        Some(Ok(())) => {
            *records_read += 1;
            Ok(Some(record))
        }
        Some(Err(error)) => bail!(
            "failed to decode {input_label} BAM record {}: {error}",
            *records_read + 1
        ),
        None => Ok(None),
    }
}

impl<D: Read, A: Read> Iterator for ZipRecordsIter<D, A> {
    type Item = anyhow::Result<RepairJob>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.terminated {
            return None;
        }
        match self.next_job() {
            Ok(Some(job)) => Some(Ok(job)),
            Ok(None) => {
                self.terminated = true;
                None
            }
            Err(error) => {
                self.terminated = true;
                Some(Err(error))
            }
        }
    }
}

fn repair_record_pair(record_pair: RecordPair) -> anyhow::Result<bam::Record> {
    let read_name =
        get_query_name_string(&record_pair.donor).unwrap_or_else(|e| {
            format!("failed to parse query name, {}", e.to_string())
        });
    let modbase_info = ModBaseInfo::new_from_record(&record_pair.donor)
        .map_err(|e| anyhow!("record {read_name} failed, {}", e.to_string()))?;

    let donor_seq =
        get_forward_sequence_str(&record_pair.donor).map_err(|e| {
            anyhow!(
                "donor sequence for record {read_name} failed, {}",
                e.to_string()
            )
        })?;
    let acceptor_seq = get_forward_sequence_str(&record_pair.acceptor)
        .map_err(|e| {
            anyhow!(
                "acceptor sequence for record {read_name} failed, {}",
                e.to_string()
            )
        })?;

    if donor_seq.len() < acceptor_seq.len() {
        bail!(
            "donor sequence for {read_name} is shorter than acceptor sequence"
        )
    }
    let finder = Finder::new(acceptor_seq.as_bytes());
    let Some(start) = finder.find(donor_seq.as_bytes()) else {
        bail!("acceptor sequence is not a substring of the donor sequence")
    };
    if finder.find(&donor_seq.as_bytes()[start + 1..]).is_some() {
        bail!("multiple potential corrections found for {read_name}")
    } else {
        let acceptor_seq_len = acceptor_seq.len();
        let end = start + acceptor_seq_len;

        let mm_style = modbase_info.mm_style;
        let ml_style = modbase_info.ml_style;

        let mut mm_agg = String::new();
        let mut ml_agg = Vec::new();

        let (_, base_mod_probs_iter) = modbase_info.into_iter_base_mod_probs();
        for (dna_base, strand, seq_pos_base_mod_probs) in base_mod_probs_iter {
            let converter =
                DeltaListConverter::new_base(acceptor_seq.as_bytes(), dna_base);
            let skip_mode = seq_pos_base_mod_probs.get_skip_mode();
            let adjusted = seq_pos_base_mod_probs
                .pos_to_base_mod_probs
                .into_iter()
                .filter_map(|(pos, base_mod_probs)| {
                    if pos >= start && pos < end {
                        Some((pos - start, base_mod_probs))
                    } else {
                        None
                    }
                })
                .collect::<FxHashMap<usize, BaseModProbs>>();
            let repaired_seq_pos_base_mod_probs =
                SeqPosBaseModProbs::new(skip_mode, adjusted);
            let (mm, mut ml) = format_mm_ml_tag(
                repaired_seq_pos_base_mod_probs,
                dna_base,
                &converter.cumulative_counts,
                strand,
            );
            mm_agg.push_str(&mm);
            ml_agg.extend_from_slice(&mut ml);
        }

        let mn = Aux::U32(acceptor_seq_len as u32);
        let mm = Aux::String(&mm_agg);
        let ml_arr: AuxArray<u8> = {
            let sl = &ml_agg;
            sl.into()
        };
        let ml = Aux::ArrayU8(ml_arr);

        let mut repaired_record = record_pair.acceptor;
        for tag in MM_TAGS.iter().chain(ML_TAGS.iter()) {
            let _ = repaired_record.remove_aux(tag.as_bytes());
        }
        let _ = repaired_record.remove_aux(MN_TAG.as_bytes());

        repaired_record.push_aux(mm_style.as_bytes(), mm)?;
        repaired_record.push_aux(ml_style.as_bytes(), ml)?;
        repaired_record.push_aux(MN_TAG.as_bytes(), mn)?;

        Ok(repaired_record)
    }
}

#[cfg(test)]
mod tests {
    use std::env;
    use std::fs::{self, File, OpenOptions};
    #[cfg(unix)]
    use std::os::unix::fs::{symlink, FileTypeExt, PermissionsExt};
    use std::path::{Path, PathBuf};
    use std::process::{Command, Stdio};
    use std::sync::{Arc, Mutex};
    use std::thread;
    use std::time::{Duration, Instant};

    use crossbeam_channel::{bounded, Receiver, Sender};
    use rust_htslib::bam::header::HeaderRecord;
    use rust_htslib::bam::{self, record::Aux};

    use crate::ordered_scheduler::OrderedWorker;

    use super::{
        query_name_order, repair_record_pair, run_repair_scheduler,
        validate_complete_bam, CompleteBamValidator, QueryNameOrder,
        RecordPair, RepairJob, RepairSink, RepairSlot, RepairTags,
        StagedBamOutput, StagedBamValidator, StagedFileIdentity,
        REPAIR_TEMP_PREFIX,
    };

    fn make_record(name: &[u8], sequence: &str) -> bam::Record {
        let mut record = bam::Record::new();
        record.set(name, None, sequence.as_bytes(), &vec![255; sequence.len()]);
        record
    }

    fn make_donor(name: &[u8], sequence: &str) -> bam::Record {
        let mut record = make_record(name, sequence);
        record.push_aux(b"MM", Aux::String("A+a.,0;")).unwrap();
        record.push_aux(b"ML", Aux::ArrayU8((&[200u8][..]).into())).unwrap();
        record
    }

    fn repair(donor: &str, acceptor: &str) -> anyhow::Result<bam::Record> {
        repair_record_pair(RecordPair::new(
            Arc::new(make_donor(b"read", donor)),
            make_record(b"read", acceptor),
        ))
    }

    fn repair_reverse(
        donor_stored_sequence: &str,
        acceptor_stored_sequence: &str,
    ) -> anyhow::Result<bam::Record> {
        let mut donor = make_donor(b"reverse", donor_stored_sequence);
        donor.set_reverse();
        let mut acceptor = make_record(b"reverse", acceptor_stored_sequence);
        acceptor.set_reverse();
        repair_record_pair(RecordPair::new(Arc::new(donor), acceptor))
    }

    #[test]
    fn repair_rejects_overlapping_ambiguous_placements() {
        let error = repair("ACACAC", "ACAC").unwrap_err();
        assert!(
            error.to_string().contains("multiple potential corrections"),
            "unexpected error: {error:#}"
        );

        assert!(repair("ACACGGACAC", "ACAC")
            .unwrap_err()
            .to_string()
            .contains("multiple potential corrections"));
    }

    #[test]
    fn repair_preserves_unique_and_absent_placement_behavior() {
        let repaired = repair("TACACG", "ACAC").unwrap();
        assert!(matches!(repaired.aux(b"MM").unwrap(), Aux::String("A+a.,0;")));
        match repaired.aux(b"ML").unwrap() {
            Aux::ArrayU8(qualities) => {
                assert_eq!(qualities.iter().collect::<Vec<_>>(), vec![200]);
            }
            other => panic!("unexpected ML tag: {other:?}"),
        }
        assert!(matches!(repaired.aux(b"MN").unwrap(), Aux::U32(4)));

        assert!(repair("ACACAC", "CGCG")
            .unwrap_err()
            .to_string()
            .contains("not a substring"));
    }

    #[test]
    fn repair_searches_forward_sequences_for_reverse_records() {
        // Stored placement starts at 3, while the forward-sequence placement
        // TACGAAA -> ACG starts at 1.
        let repaired = repair_reverse("TTTCGTA", "CGT").unwrap();
        assert!(repaired.is_reverse());
        assert!(matches!(repaired.aux(b"MM").unwrap(), Aux::String("A+a.,0;")));
        match repaired.aux(b"ML").unwrap() {
            Aux::ArrayU8(qualities) => {
                assert_eq!(qualities.iter().collect::<Vec<_>>(), vec![200]);
            }
            other => panic!("unexpected ML tag: {other:?}"),
        }
        assert!(matches!(repaired.aux(b"MN").unwrap(), Aux::U32(3)));

        assert!(repair_reverse("GTGTGT", "GTGT")
            .unwrap_err()
            .to_string()
            .contains("multiple potential corrections"));
    }

    #[test]
    fn repair_returns_errors_for_empty_sequences() {
        for (donor, acceptor, expected) in [
            (
                "ACAC",
                "",
                "acceptor sequence for record read failed, empty-read-sequence",
            ),
            ("", "ACAC", "record read failed, invalid-MM-tag"),
            ("", "", "record read failed, invalid-MM-tag"),
        ] {
            let error = repair(donor, acceptor).unwrap_err();
            assert_eq!(error.to_string(), expected);
        }
    }

    struct FailNthWriteSink<S> {
        inner: S,
        fail_at: usize,
        writes: usize,
    }

    impl<S: RepairSink> RepairSink for FailNthWriteSink<S> {
        fn write(&mut self, record: &bam::Record) -> anyhow::Result<()> {
            self.writes += 1;
            if self.writes == self.fail_at {
                anyhow::bail!(
                    "injected repair sink write failure at record {}",
                    self.fail_at
                )
            }
            self.inner.write(record)
        }

        fn finish(self, expected_records: usize) -> anyhow::Result<()> {
            self.inner.finish(expected_records)
        }
    }

    struct CorruptBeforeValidation;

    impl StagedBamValidator for CorruptBeforeValidation {
        fn validate(
            &self,
            file: &mut File,
            path: &Path,
            expected_identity: StagedFileIdentity,
            expected_records: usize,
        ) -> anyhow::Result<()> {
            let file_len = file.metadata()?.len();
            file.set_len(file_len.checked_sub(1).expect("BAM is non-empty"))?;
            CompleteBamValidator.validate(
                file,
                path,
                expected_identity,
                expected_records,
            )
        }
    }

    #[cfg(unix)]
    struct SubstituteBeforeDecode;

    #[cfg(unix)]
    impl StagedBamValidator for SubstituteBeforeDecode {
        fn validate(
            &self,
            file: &mut File,
            path: &Path,
            expected_identity: StagedFileIdentity,
            expected_records: usize,
        ) -> anyhow::Result<()> {
            let replacement = path.with_file_name("replacement.bam");
            let replacement_writer = bam::Writer::from_path(
                &replacement,
                &bam::Header::new(),
                bam::Format::Bam,
            )?;
            drop(replacement_writer);
            fs::rename(&replacement, path)?;
            CompleteBamValidator.validate(
                file,
                path,
                expected_identity,
                expected_records,
            )
        }
    }

    const REPAIR_CHILD_CASE_ENV: &str = "MODKIT_REPAIR_CHILD_CASE";
    const REPAIR_CHILD_DIR_ENV: &str = "MODKIT_REPAIR_CHILD_DIR";
    const REPAIR_CHILD_TEST: &str =
        "repair_tags::tests::repair_lifecycle_child";

    fn write_lifecycle_bam(path: &Path, records: usize, donor: bool) {
        let mut header = bam::Header::new();
        let mut hd = HeaderRecord::new(b"HD");
        hd.push_tag(b"VN", "1.6")
            .push_tag(b"SO", "queryname")
            .push_tag(b"SS", "queryname:natural");
        header.push_record(&hd);
        let mut writer =
            bam::Writer::from_path(path, &header, bam::Format::Bam).unwrap();
        for index in 0..records {
            let name = format!("read{index:06}");
            let sequence = if donor { "TACGT" } else { "ACG" };
            let mut record = make_record(name.as_bytes(), sequence);
            if donor {
                record.push_aux(b"MM", Aux::String("A+a.,0;")).unwrap();
                let probability = [((index % 200) + 1) as u8];
                record
                    .push_aux(b"ML", Aux::ArrayU8((&probability[..]).into()))
                    .unwrap();
            }
            writer.write(&record).unwrap();
        }
        drop(writer);
        validate_complete_bam(path, records).unwrap();
    }

    fn truncate_bam_data(path: &Path) {
        let file = OpenOptions::new().write(true).open(path).unwrap();
        let file_len = file.metadata().unwrap().len();
        assert!(file_len > 40);
        // Remove the canonical 28-byte EOF plus part of the final data block.
        file.set_len(file_len - 40).unwrap();
    }

    fn lifecycle_paths(root: &Path) -> (PathBuf, PathBuf, PathBuf) {
        (
            root.join("donor.bam"),
            root.join("acceptor.bam"),
            root.join("output.bam"),
        )
    }

    fn run_lifecycle_child(case: &str, root: &Path) {
        let mut child = Command::new(env::current_exe().unwrap())
            .args(["--exact", REPAIR_CHILD_TEST, "--nocapture"])
            .env(REPAIR_CHILD_CASE_ENV, case)
            .env(REPAIR_CHILD_DIR_ENV, root)
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .unwrap();
        let deadline = Instant::now() + Duration::from_secs(10);
        loop {
            if child.try_wait().unwrap().is_some() {
                break;
            }
            if Instant::now() >= deadline {
                let _ = child.kill();
                let output = child.wait_with_output().unwrap();
                panic!(
                    "repair lifecycle child {case} exceeded watchdog; \
                     stdout={} stderr={}",
                    String::from_utf8_lossy(&output.stdout),
                    String::from_utf8_lossy(&output.stderr)
                );
            }
            thread::sleep(Duration::from_millis(10));
        }
        let output = child.wait_with_output().unwrap();
        assert!(
            output.status.success(),
            "repair lifecycle child {case} failed; stdout={} stderr={}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
    }

    fn staged_repair_entries(root: &Path) -> Vec<PathBuf> {
        fs::read_dir(root)
            .unwrap()
            .map(|entry| entry.unwrap().path())
            .filter(|path| {
                path.file_name()
                    .unwrap()
                    .to_string_lossy()
                    .starts_with(REPAIR_TEMP_PREFIX)
            })
            .collect()
    }

    #[test]
    fn repair_lifecycle_child() {
        let Ok(case) = env::var(REPAIR_CHILD_CASE_ENV) else {
            return;
        };
        let root = PathBuf::from(env::var_os(REPAIR_CHILD_DIR_ENV).unwrap());
        let (donor_bam, acceptor_bam, output_bam) = lifecycle_paths(&root);
        let repair = RepairTags {
            donor_bam,
            acceptor_bam,
            output_bam,
            log_filepath: None,
            threads: 2,
        };
        let result = match case.as_str() {
            "nth-write" => repair.run_with_output_factory(|header| {
                Ok(FailNthWriteSink {
                    inner: StagedBamOutput::new(&repair.output_bam, header)?,
                    fail_at: 2,
                    writes: 0,
                })
            }),
            "validation" => repair.run_with_output_factory(|header| {
                StagedBamOutput::with_validator(
                    &repair.output_bam,
                    header,
                    CorruptBeforeValidation,
                )
            }),
            _ => repair.run_with_output_factory(|header| {
                StagedBamOutput::new(&repair.output_bam, header)
            }),
        };
        match case.as_str() {
            "success" => result.expect("successful repair should finish"),
            "truncated-donor" => {
                let error = result.expect_err("truncated donor must fail");
                assert!(
                    format!("{error:#}")
                        .contains("failed to decode donor BAM record"),
                    "unexpected truncated-donor error: {error:#}"
                );
            }
            "truncated-acceptor" => {
                let error = result.expect_err("truncated acceptor must fail");
                assert!(
                    format!("{error:#}")
                        .contains("failed to decode acceptor BAM record"),
                    "unexpected truncated-acceptor error: {error:#}"
                );
            }
            "nth-write" => {
                let error = result.expect_err("injected write must fail");
                assert!(
                    format!("{error:#}").contains(
                        "injected repair sink write failure at record 2"
                    ),
                    "unexpected injected-write error: {error:#}"
                );
            }
            "validation" => {
                let error = result.expect_err("corrupt staged BAM must fail");
                assert!(
                    format!("{error:#}")
                        .contains("missing the canonical BGZF EOF block"),
                    "unexpected validation error: {error:#}"
                );
            }
            other => panic!("unknown repair lifecycle child case {other}"),
        }
    }

    #[test]
    fn repair_lifecycle_failures_preserve_outputs_and_clean_staging() {
        const SENTINEL: &[u8] = b"existing repair output sentinel";
        for case in [
            "truncated-donor",
            "truncated-acceptor",
            "nth-write",
            "validation",
            "success",
        ] {
            for output_exists in [false, true] {
                let temp_dir = tempfile::tempdir().unwrap();
                let (donor, acceptor, output) =
                    lifecycle_paths(temp_dir.path());
                let records =
                    if case.starts_with("truncated-") { 2_000 } else { 3 };
                write_lifecycle_bam(&donor, records, true);
                write_lifecycle_bam(&acceptor, records, false);
                match case {
                    "truncated-donor" => truncate_bam_data(&donor),
                    "truncated-acceptor" => truncate_bam_data(&acceptor),
                    _ => {}
                }
                if output_exists {
                    fs::write(&output, SENTINEL).unwrap();
                }

                run_lifecycle_child(case, temp_dir.path());

                if case == "success" {
                    validate_complete_bam(&output, records).unwrap();
                    assert_ne!(fs::read(&output).unwrap(), SENTINEL);
                } else if output_exists {
                    assert_eq!(fs::read(&output).unwrap(), SENTINEL);
                } else {
                    assert!(!output.exists());
                }
                let staged = staged_repair_entries(temp_dir.path());
                assert!(
                    staged.is_empty(),
                    "staged repair entries remain: {staged:?}"
                );
            }
        }
    }

    #[cfg(unix)]
    #[test]
    fn staged_output_matches_file_create_mode_and_preserves_existing_mode() {
        let temp_dir = tempfile::tempdir().unwrap();
        // Model a shared, non-sticky destination directory. The staged BAM
        // must still live under an owner-only child directory.
        fs::set_permissions(temp_dir.path(), fs::Permissions::from_mode(0o777))
            .unwrap();
        let header = bam::Header::new();

        let control = temp_dir.path().join("ordinary-file-create");
        drop(File::create(&control).unwrap());
        let ordinary_create_mode =
            fs::metadata(&control).unwrap().permissions().mode() & 0o7777;

        let new_output = temp_dir.path().join("new-output.bam");
        let staged_output = StagedBamOutput::new(&new_output, &header).unwrap();
        let staging_dir = staged_output._staging_dir.path().to_path_buf();
        assert_eq!(staging_dir.parent(), Some(temp_dir.path()));
        assert_eq!(
            fs::metadata(&staging_dir).unwrap().permissions().mode() & 0o077,
            0
        );
        assert_eq!(
            staged_output.temp_file.path().parent(),
            Some(staging_dir.as_path())
        );
        assert_eq!(
            staged_output
                .temp_file
                .as_file()
                .metadata()
                .unwrap()
                .permissions()
                .mode()
                & 0o7777,
            ordinary_create_mode
        );
        staged_output.finish(0).unwrap();
        assert!(!staging_dir.exists());
        let new_output_mode =
            fs::metadata(&new_output).unwrap().permissions().mode() & 0o7777;
        assert_eq!(new_output_mode, ordinary_create_mode);

        let existing_output = temp_dir.path().join("existing-output.bam");
        fs::write(&existing_output, b"sentinel").unwrap();
        let deliberate_mode = 0o404;
        fs::set_permissions(
            &existing_output,
            fs::Permissions::from_mode(deliberate_mode),
        )
        .unwrap();
        StagedBamOutput::new(&existing_output, &header)
            .unwrap()
            .finish(0)
            .unwrap();
        let replacement_mode =
            fs::metadata(&existing_output).unwrap().permissions().mode()
                & 0o7777;
        assert_eq!(replacement_mode, deliberate_mode);
        assert!(staged_repair_entries(temp_dir.path()).is_empty());
    }

    #[cfg(unix)]
    #[test]
    fn staged_inode_substitution_is_fatal_and_cleans_private_directory() {
        const SENTINEL: &[u8] = b"existing repair output sentinel";
        for output_exists in [false, true] {
            let temp_dir = tempfile::tempdir().unwrap();
            let output = temp_dir.path().join("output.bam");
            if output_exists {
                fs::write(&output, SENTINEL).unwrap();
            }
            let staged_output = StagedBamOutput::with_validator(
                &output,
                &bam::Header::new(),
                SubstituteBeforeDecode,
            )
            .unwrap();
            let staging_dir = staged_output._staging_dir.path().to_path_buf();

            let error = staged_output
                .finish(0)
                .expect_err("replacement inode must not be persisted");
            assert!(
                format!("{error:#}")
                    .contains("no longer identifies the retained file"),
                "unexpected substitution error: {error:#}"
            );
            if output_exists {
                assert_eq!(fs::read(&output).unwrap(), SENTINEL);
            } else {
                assert!(!output.exists());
            }
            assert!(!staging_dir.exists());
            assert!(staged_repair_entries(temp_dir.path()).is_empty());
        }
    }

    #[cfg(unix)]
    #[test]
    fn staged_output_rejects_symlink_and_nonregular_destinations() {
        const SENTINEL: &[u8] = b"symlink target sentinel";
        let temp_dir = tempfile::tempdir().unwrap();
        let header = bam::Header::new();

        let target = temp_dir.path().join("target.bam");
        fs::write(&target, SENTINEL).unwrap();
        let symlink_output = temp_dir.path().join("symlink-output.bam");
        symlink(&target, &symlink_output).unwrap();
        let symlink_error = match StagedBamOutput::new(&symlink_output, &header)
        {
            Ok(_) => panic!("symlink destination should be rejected"),
            Err(error) => error,
        };
        assert!(
            format!("{symlink_error:#}")
                .contains("must not be a symbolic link"),
            "unexpected symlink error: {symlink_error:#}"
        );
        assert_eq!(fs::read(&target).unwrap(), SENTINEL);
        assert!(fs::symlink_metadata(&symlink_output)
            .unwrap()
            .file_type()
            .is_symlink());

        let dangling_output = temp_dir.path().join("dangling-output.bam");
        symlink(temp_dir.path().join("missing-target.bam"), &dangling_output)
            .unwrap();
        let dangling_error =
            match StagedBamOutput::new(&dangling_output, &header) {
                Ok(_) => {
                    panic!("dangling symlink destination should be rejected")
                }
                Err(error) => error,
            };
        assert!(
            format!("{dangling_error:#}")
                .contains("must not be a symbolic link"),
            "unexpected dangling-symlink error: {dangling_error:#}"
        );
        assert!(fs::symlink_metadata(&dangling_output)
            .unwrap()
            .file_type()
            .is_symlink());

        let fifo_output = temp_dir.path().join("fifo-output.bam");
        let fifo_status =
            Command::new("mkfifo").arg(&fifo_output).status().unwrap();
        assert!(fifo_status.success(), "failed to create FIFO fixture");
        let fifo_error = match StagedBamOutput::new(&fifo_output, &header) {
            Ok(_) => panic!("FIFO destination should be rejected"),
            Err(error) => error,
        };
        assert!(
            format!("{fifo_error:#}")
                .contains("must be a regular file when it already exists"),
            "unexpected FIFO error: {fifo_error:#}"
        );
        assert!(fs::symlink_metadata(&fifo_output)
            .unwrap()
            .file_type()
            .is_fifo());

        let directory_output = temp_dir.path().join("directory-output.bam");
        fs::create_dir(&directory_output).unwrap();
        let directory_error =
            match StagedBamOutput::new(&directory_output, &header) {
                Ok(_) => panic!("directory destination should be rejected"),
                Err(error) => error,
            };
        assert!(
            format!("{directory_error:#}")
                .contains("must be a regular file when it already exists"),
            "unexpected directory error: {directory_error:#}"
        );
        assert!(directory_output.is_dir());
        assert!(staged_repair_entries(temp_dir.path()).is_empty());
    }

    struct GatedRepairWorker {
        job_one_completed: Sender<()>,
        release_job_zero: Receiver<()>,
        completions: Arc<Mutex<Vec<usize>>>,
    }

    impl OrderedWorker<RepairJob, RepairSlot> for GatedRepairWorker {
        fn process(
            &mut self,
            job: RepairJob,
            mut slot: RepairSlot,
        ) -> anyhow::Result<RepairSlot> {
            let RepairJob::Matched(record_pair) = job else {
                panic!("ordering fixture only uses matched jobs")
            };
            let job = match record_pair.acceptor.qname() {
                b"job0" => 0,
                b"job1" => 1,
                name => panic!("unexpected job name: {:?}", name),
            };
            if job == 0 {
                self.release_job_zero
                    .recv_timeout(Duration::from_secs(2))
                    .expect("job 0 was not released");
            }
            self.completions.lock().unwrap().push(job);
            if job == 1 {
                let _ = self.job_one_completed.send(());
            }
            slot.outcome = Some(Ok(record_pair.acceptor));
            Ok(slot)
        }
    }

    fn matched_job(name: &[u8]) -> anyhow::Result<RepairJob> {
        Ok(RepairJob::Matched(RecordPair::new(
            Arc::new(make_record(name, "ACG")),
            make_record(name, "ACG"),
        )))
    }

    #[test]
    fn repair_scheduler_preserves_acceptor_order_after_reordered_completion() {
        let completions = Arc::new(Mutex::new(Vec::new()));
        let consumed = Arc::new(Mutex::new(Vec::new()));
        let (job_one_completed, wait_for_job_one) = bounded(1);
        let (release_job_zero, wait_for_job_zero_release) = bounded(1);
        let workers = (0..2)
            .map(|_| GatedRepairWorker {
                job_one_completed: job_one_completed.clone(),
                release_job_zero: wait_for_job_zero_release.clone(),
                completions: completions.clone(),
            })
            .collect();
        let consumed_probe = consumed.clone();
        let (finished, wait_for_finish) = bounded(1);
        let handle = thread::spawn(move || {
            let result = run_repair_scheduler(
                vec![matched_job(b"job0"), matched_job(b"job1")].into_iter(),
                workers,
                0,
                move |record| {
                    consumed_probe
                        .lock()
                        .unwrap()
                        .push(record?.qname().to_vec());
                    Ok(())
                },
            );
            let _ = finished.send(result);
        });

        wait_for_job_one
            .recv_timeout(Duration::from_secs(2))
            .expect("job 1 did not complete while job 0 was blocked");
        assert_eq!(*completions.lock().unwrap(), vec![1]);
        assert!(
            consumed.lock().unwrap().is_empty(),
            "later repair was consumed before sequence 0"
        );
        release_job_zero.send(()).unwrap();
        wait_for_finish
            .recv_timeout(Duration::from_secs(2))
            .expect("repair scheduler watchdog expired")
            .unwrap();
        handle.join().expect("repair scheduler thread panicked");

        assert_eq!(*completions.lock().unwrap(), vec![1, 0]);
        assert_eq!(
            *consumed.lock().unwrap(),
            vec![b"job0".to_vec(), b"job1".to_vec()]
        );
    }

    #[test]
    fn query_name_comparators_cover_numeric_and_equivalent_names() {
        assert_eq!(
            QueryNameOrder::Natural
                .common_merge_order(QueryNameOrder::LegacyNatural),
            None
        );
        assert_eq!(
            QueryNameOrder::LegacyNatural
                .common_merge_order(QueryNameOrder::Natural),
            None
        );
        assert_eq!(
            QueryNameOrder::LegacyNatural
                .common_merge_order(QueryNameOrder::Lexicographical),
            None
        );
        assert_eq!(
            QueryNameOrder::Lexicographical
                .common_merge_order(QueryNameOrder::LegacyNatural),
            None
        );
        assert_eq!(
            QueryNameOrder::Natural
                .common_merge_order(QueryNameOrder::Lexicographical),
            None
        );
        assert_eq!(
            QueryNameOrder::Lexicographical
                .common_merge_order(QueryNameOrder::Natural),
            None
        );
        assert_eq!(
            QueryNameOrder::Natural.compare(b"read1", b"read2"),
            std::cmp::Ordering::Less
        );
        assert_eq!(
            QueryNameOrder::Natural.compare(b"read2", b"read10"),
            std::cmp::Ordering::Less
        );
        assert_eq!(
            QueryNameOrder::Natural.compare(b"read1", b"read01"),
            std::cmp::Ordering::Equal
        );
        assert_eq!(
            QueryNameOrder::LegacyNatural.compare(b"read01", b"read1"),
            std::cmp::Ordering::Less
        );
        assert_eq!(
            QueryNameOrder::Natural.compare(
                b"read999999999999999999999999",
                b"read1000000000000000000000000",
            ),
            std::cmp::Ordering::Less
        );
        assert_eq!(
            QueryNameOrder::Lexicographical.compare(b"read10", b"read2"),
            std::cmp::Ordering::Less
        );
    }

    #[test]
    fn query_name_sort_is_read_from_header_and_must_be_supported() {
        let natural = bam::HeaderView::from_bytes(
            b"@HD\tVN:1.6\tSO:queryname\tSS:queryname:natural\n",
        );
        let lexical = bam::HeaderView::from_bytes(
            b"@HD\tVN:1.6\tSO:queryname\tSS:queryname:lexicographical\n",
        );
        let fallback =
            bam::HeaderView::from_bytes(b"@HD\tVN:1.5\tSO:queryname\n");
        let coordinate =
            bam::HeaderView::from_bytes(b"@HD\tVN:1.6\tSO:coordinate\n");
        let malformed = bam::HeaderView::from_bytes(
            b"@HD\tVN:1.6\tSO:queryname\tSS:queryname:natural:extra\n",
        );

        assert_eq!(
            query_name_order(&natural, "test").unwrap(),
            QueryNameOrder::Natural
        );
        assert_eq!(
            query_name_order(&lexical, "test").unwrap(),
            QueryNameOrder::Lexicographical
        );
        assert_eq!(
            query_name_order(&fallback, "test").unwrap(),
            QueryNameOrder::LegacyNatural
        );
        assert!(query_name_order(&coordinate, "test")
            .unwrap_err()
            .to_string()
            .contains("must be query-name sorted"));
        assert!(query_name_order(&malformed, "test")
            .unwrap_err()
            .to_string()
            .contains("unsupported query-name sub-sort"));
    }
}
