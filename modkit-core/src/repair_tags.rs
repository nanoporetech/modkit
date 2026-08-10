use std::cmp::Ordering;
use std::path::PathBuf;
use std::sync::Arc;

use anyhow::{anyhow, bail};
use clap::Args;
use crossbeam_channel::unbounded;
use derive_new::new;
use indicatif::{MultiProgress, ProgressBar};
use log::{debug, error, info, warn};
use rust_htslib::bam::record::{Aux, AuxArray};
use rust_htslib::bam::{self, Read};
use rustc_hash::FxHashMap;

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

impl RepairTags {
    pub fn run(&self) -> anyhow::Result<()> {
        let _handle = init_logging(self.log_filepath.as_ref());

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

        let mut donor_records = bam::Reader::from_path(&self.donor_bam)?;
        donor_records.set_threads(threads_per_reader)?;
        let mut acceptor_records = bam::Reader::from_path(&self.acceptor_bam)?;
        acceptor_records.set_threads(threads_per_reader)?;
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
        let mut writer = bam::Writer::from_path(
            &self.output_bam,
            &header,
            bam::Format::Bam,
        )?;
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
                    if let Err(e) = writer.write(&record) {
                        error!("failed to write record {}", e.to_string());
                        n_failed += 1;
                    } else {
                        written_ticker.inc(1);
                        n_repaired += 1;
                    }
                }
                Err(e) => {
                    debug!("record failed to be repaired: {}", e.to_string());
                    n_failed += 1;
                }
            }
            Ok(())
        })?;

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
            let Some(record) =
                get_next_record(&mut self.donor_records, "donor")
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
        let Some(record) =
            get_next_record(&mut self.acceptor_records, "acceptor")
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
) -> Option<bam::Record> {
    loop {
        let mut record = bam::Record::new();
        match records.read(&mut record) {
            Some(Ok(())) => return Some(record),
            Some(Err(error)) => {
                // Commit 2 makes malformed input fatal and contextual. Keep
                // the current skip-and-warn policy in this ordering-only unit.
                warn!("failed to parse record from {input_label} BAM, {error}");
            }
            None => return None,
        }
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
        bail!("donor sequence for {read_name} is longer than acceptor sequence")
    }
    let matches = donor_seq.match_indices(&acceptor_seq);

    let starts =
        matches.into_iter().map(|(start, _)| start).collect::<Vec<usize>>();
    if starts.len() > 1 {
        bail!("multiple potential corrections found for {read_name}")
    } else if starts.is_empty() {
        bail!("acceptor sequence is not a substring of the donor sequence")
    } else {
        let acceptor_seq_len = acceptor_seq.len();
        let start = *starts.get(0).unwrap();
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
    use std::sync::{Arc, Mutex};
    use std::thread;
    use std::time::Duration;

    use crossbeam_channel::{bounded, Receiver, Sender};
    use rust_htslib::bam;

    use crate::ordered_scheduler::OrderedWorker;

    use super::{
        query_name_order, run_repair_scheduler, QueryNameOrder, RecordPair,
        RepairJob, RepairSlot,
    };

    fn make_record(name: &[u8], sequence: &str) -> bam::Record {
        let mut record = bam::Record::new();
        record.set(name, None, sequence.as_bytes(), &vec![255; sequence.len()]);
        record
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
