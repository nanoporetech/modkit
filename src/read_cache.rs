use std::collections::HashSet;

use log::{debug, error};
use rust_htslib::bam;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::errs::RunError;
use crate::mod_bam::{
    collapse_mod_probs, BaseModCall, CollapseMethod, DuplexModCall, EdgeFilter,
    ModBaseInfo, SeqPosBaseModProbs, SkipMode,
};
use crate::mod_base_code::{DnaBase, ModCode};
use crate::motif_bed::{MotifLocations, MultipleMotifLocations};
use crate::threshold_mod_caller::MultipleThresholdModCaller;
use crate::util::{self, Strand};

/// Mapping of _reference position_ to base mod calls as determined by the aligned pairs for the
/// read
type RefPosBaseModCalls = FxHashMap<u64, BaseModCall>;

pub(crate) struct ReadCache<'a> {
    /// Mapping of read_id to reference position <> base mod calls for that read
    /// organized by the canonical base (the 'char') todo: should use DnaBase here
    pos_reads:
        FxHashMap<String, FxHashMap<char, (RefPosBaseModCalls, SkipMode)>>,
    neg_reads:
        FxHashMap<String, FxHashMap<char, (RefPosBaseModCalls, SkipMode)>>,
    /// these reads don't have mod tags or should be skipped for some other reason
    skip_set: HashSet<String>,
    /// mapping of read_id (query_name) to the mod codes contained in that read
    pos_mod_codes: FxHashMap<String, FxHashSet<ModCode>>,
    neg_mod_codes: FxHashMap<String, FxHashSet<ModCode>>,
    /// collapse method
    method: Option<&'a CollapseMethod>,
    /// Force allowing of implicit canonical
    force_allow: bool,
    caller: &'a MultipleThresholdModCaller,
    /// Edge filter to remove base mod calls at the ends of reads
    edge_filter: Option<&'a EdgeFilter>,
}

impl<'a> ReadCache<'a> {
    pub(crate) fn new(
        method: Option<&'a CollapseMethod>,
        caller: &'a MultipleThresholdModCaller,
        edge_filter: Option<&'a EdgeFilter>,
        force_allow: bool,
    ) -> Self {
        Self {
            pos_reads: FxHashMap::default(),
            neg_reads: FxHashMap::default(),
            skip_set: HashSet::new(),
            pos_mod_codes: FxHashMap::default(),
            neg_mod_codes: FxHashMap::default(),
            method,
            force_allow,
            caller,
            edge_filter,
        }
    }

    /// Subroutine that adds read's mod base calls to the cache (or error),
    /// in the case of an error the caller could remove this read from
    /// future consideration
    #[inline]
    fn add_modbase_probs_for_record_and_canonical_base(
        &mut self,
        record_name: &str,
        record: &bam::Record,
        seq_pos_base_mod_probs: SeqPosBaseModProbs,
        mod_strand: Strand,
        canonical_base: DnaBase,
        threshold_base: DnaBase,
    ) -> Result<(), RunError> {
        let aligned_pairs = util::get_aligned_pairs_forward(&record)
            .filter_map(|ap| ap.ok())
            .collect::<FxHashMap<usize, u64>>();

        let mut warned = false;
        let ref_pos_base_mod_calls = seq_pos_base_mod_probs
            .pos_to_base_mod_probs
            .into_iter() // par iter?
            // here the q_pos is the forward-oriented position
            .flat_map(|(q_pos, bmp)| {
                if let Some(r_pos) = aligned_pairs.get(&q_pos) {
                    // filtering happens here.
                    match self.caller.call(&threshold_base, &bmp) {
                        Ok(base_mod_call) => Some((*r_pos, base_mod_call)),
                        Err(e) => {
                            if !warned {
                                let warning = e.to_string();
                                let e = e.context(format!(
                                    "{} has out of spec mod code, {warning}",
                                    record_name
                                ));
                                debug!("{}, update tags with modkit adjust-mods --convert",
                                e.to_string());
                                warned = true
                            }
                            None
                        }
                    }
                } else {
                    None
                }
            })
            .collect::<FxHashMap<u64, BaseModCall>>();

        let read_table = match mod_strand {
            Strand::Positive => &mut self.pos_reads,
            Strand::Negative => &mut self.neg_reads,
        };
        read_table
            .entry(record_name.to_owned())
            .or_insert(FxHashMap::default())
            .insert(
                canonical_base.char(),
                (ref_pos_base_mod_calls, seq_pos_base_mod_probs.skip_mode),
            );
        Ok(())
    }

    /// Add a record to the cache.
    fn add_record(&mut self, record: &bam::Record) -> Result<(), RunError> {
        let record_name = util::get_query_name_string(record)
            .map_err(|e| RunError::new_input_error(e.to_string()))?;

        let mod_base_info = ModBaseInfo::new_from_record(record)?;
        if mod_base_info.is_empty() {
            let msg = format!("record {} has no mod calls", &record_name);
            return Err(RunError::Skipped(msg));
        }

        // todo(ar) shouln't have to perform this sweep, should be able to keep
        //  a temporary container and update the cache after all seq_pos_probs
        //  are checked to be OK.
        for (_base, _strand, seq_pos_probs) in
            mod_base_info.iter_seq_base_mod_probs()
        {
            if seq_pos_probs.skip_mode == SkipMode::ImplicitProbModified
                && !self.force_allow
            {
                let msg = format!(
                    "record {} has un-allowed mode ({:?}), use \
                '--force-allow-implicit' or 'modkit update-tags --mode ambiguous'",
                    &record_name, seq_pos_probs.skip_mode
                );
                return Err(RunError::Skipped(msg));
            }
        }

        // need a flag here, change to true iff we add probs for at least a single canonical
        // base if they are all filtered out (due to edge filter), return an Err so that we
        // don't re-process this read.
        let mut added_base_mod_probs = false;
        let (_, mod_prob_iter) = mod_base_info.into_iter_base_mod_probs();
        for (base, mod_strand, seq_base_mod_probs) in mod_prob_iter {
            match DnaBase::parse(base) {
                Ok(dna_base) => {
                    let threshold_base = match mod_strand {
                        Strand::Positive => dna_base,
                        Strand::Negative => dna_base.complement(),
                    };
                    let seq_base_mod_probs =
                        if let Some(edge_filter) = &self.edge_filter {
                            seq_base_mod_probs.edge_filter_positions(
                                edge_filter,
                                record.seq_len(),
                            )
                        } else {
                            Some(seq_base_mod_probs)
                        };
                    // not idiomatic, but fights rightward drift..?
                    if seq_base_mod_probs.is_none() {
                        debug!("all base mod positions were removed by edge filter \
                                for {record_name} and base {base}");
                        continue;
                    }
                    let mut seq_base_mod_probs = seq_base_mod_probs.unwrap();
                    if let Some(method) = &self.method {
                        seq_base_mod_probs =
                            collapse_mod_probs(seq_base_mod_probs, method);
                    }
                    // could move this into it's own routine..?
                    let mod_codes = seq_base_mod_probs
                        .pos_to_base_mod_probs
                        .values()
                        .flat_map(|base_mod_probs| {
                            base_mod_probs.iter_probs().filter_map(
                                |(&raw_mod_code, _)| {
                                    ModCode::parse_raw_mod_code(raw_mod_code)
                                        .ok()
                                },
                            )
                        })
                        .collect::<FxHashSet<ModCode>>();

                    let mod_codes_for_read =
                        match (mod_strand, record.is_reverse()) {
                            // C+C positive stranded
                            (Strand::Positive, false) => {
                                &mut self.pos_mod_codes
                            }
                            (Strand::Positive, true) => &mut self.neg_mod_codes,
                            // G-C negative stranded
                            (Strand::Negative, false) => {
                                &mut self.neg_mod_codes
                            }
                            (Strand::Negative, true) => &mut self.pos_mod_codes,
                        };
                    mod_codes_for_read
                        .entry(record_name.to_owned())
                        .or_insert(FxHashSet::default())
                        .extend(mod_codes);

                    self.add_modbase_probs_for_record_and_canonical_base(
                        &record_name,
                        record,
                        seq_base_mod_probs,
                        mod_strand,
                        dna_base,
                        threshold_base,
                    )?;
                    added_base_mod_probs = true
                }
                Err(e) => {
                    let message = format!(
                        "record {record_name} has unallowed DNA base {base}, {}", e.to_string()
                    );
                    debug!("{}", &message);
                    // short circuit
                    return Err(RunError::new_failed(message));
                }
            }
        }
        if added_base_mod_probs {
            Ok(())
        } else {
            Err(RunError::Skipped(format!(
                "all base mod positions removed in filtering"
            )))
        }
    }

    #[inline]
    fn get_mod_call_from_mapping(
        strand_calls: &FxHashMap<char, (RefPosBaseModCalls, SkipMode)>,
        canonical_base: char,
        position: u32,
    ) -> Option<BaseModCall> {
        strand_calls.get(&canonical_base).and_then(
            |(ref_pos_mod_calls, skip_mode)| {
                let mod_base_call =
                    ref_pos_mod_calls.get(&(position as u64)).map(|bmc| *bmc);
                // .map(|base_mod_call| {
                //     Self::filter_base_mod_call(base_mod_call, threshold)
                // });
                match skip_mode {
                    SkipMode::Ambiguous => mod_base_call,
                    SkipMode::ImplicitProbModified | SkipMode::ProbModified => {
                        Some(
                            mod_base_call
                                .unwrap_or(BaseModCall::Canonical(1.0f32)),
                        )
                    }
                }
            },
        )
    }

    /// Get the mod call for a reference position from a read. If this read is
    /// in the cache, look it up, if not parse the tags, add it to the cache
    /// and return the mod call (if present). (+ strand, - strand) calls
    /// returned, in this case the strand is the strand of the read. So most
    /// of the time this is positive, because reads are single stranded DNA.
    /// In the case of duplex or any situation where there is a reporting of a
    /// modification on the opposite strand (or both) you can get negative
    /// strand calls or positive and negative strand calls.
    pub(crate) fn get_mod_call(
        &mut self,
        record: &bam::Record,
        position: u32,
        canonical_base: char, // todo make this DnaBase
    ) -> (Option<BaseModCall>, Option<BaseModCall>) {
        let read_id = String::from_utf8(record.qname().to_vec()).unwrap();
        if self.skip_set.contains(&read_id) {
            (None, None)
        } else {
            match (self.pos_reads.get(&read_id), self.neg_reads.get(&read_id)) {
                (Some(pos_base_mod_calls), Some(neg_base_mod_calls)) => (
                    Self::get_mod_call_from_mapping(
                        pos_base_mod_calls,
                        canonical_base,
                        position,
                    ),
                    Self::get_mod_call_from_mapping(
                        neg_base_mod_calls,
                        canonical_base,
                        position,
                    ),
                ),
                (Some(pos_base_mod_calls), None) => (
                    Self::get_mod_call_from_mapping(
                        pos_base_mod_calls,
                        canonical_base,
                        position,
                    ),
                    None,
                ),
                (None, Some(neg_base_mod_calls)) => (
                    None,
                    Self::get_mod_call_from_mapping(
                        neg_base_mod_calls,
                        canonical_base,
                        position,
                    ),
                ),
                (None, None) => {
                    match self.add_record(record) {
                        Ok(_) => {}
                        Err(run_error) => {
                            debug!(
                                "read {read_id} failed to get mod tags {}",
                                run_error.to_string()
                            );
                            self.skip_set.insert(read_id.clone());
                        }
                    }
                    if !(self.skip_set.contains(&read_id)
                        || self.pos_reads.contains_key(&read_id)
                        || self.neg_reads.contains_key(&read_id))
                    {
                        error!("didn't add failed read id to skip sets, likely a bug");
                    }
                    assert!(
                        self.skip_set.contains(&read_id)
                            || self.pos_reads.contains_key(&read_id)
                            || self.neg_reads.contains_key(&read_id),
                    );
                    self.get_mod_call(record, position, canonical_base)
                }
            }
        }
    }

    pub(crate) fn add_mod_codes_for_record(
        &mut self,
        record: &bam::Record,
        pos_strand_mod_codes: &mut HashSet<ModCode>,
        neg_strand_mod_codes: &mut HashSet<ModCode>,
    ) {
        // optimize, could use a better implementation here - pass the read_id
        // from the calling function perhaps
        let read_id = String::from_utf8(record.qname().to_vec()).unwrap();
        if self.skip_set.contains(&read_id) {
            return;
        } else {
            match (
                self.pos_mod_codes.get(&read_id),
                self.neg_mod_codes.get(&read_id),
            ) {
                (Some(pos_codes), Some(neg_codes)) => {
                    pos_strand_mod_codes.extend(pos_codes.iter().map(|mc| *mc));
                    neg_strand_mod_codes.extend(neg_codes.iter().map(|mc| *mc));
                }
                (Some(pos_codes), None) => {
                    pos_strand_mod_codes.extend(pos_codes.iter().map(|mc| *mc));
                }
                (None, Some(neg_codes)) => {
                    neg_strand_mod_codes.extend(neg_codes.iter().map(|mc| *mc));
                }
                (None, None) => {
                    match self.add_record(record) {
                        Ok(_) => {}
                        Err(run_error) => {
                            debug!("read {read_id}, {}", run_error.to_string());
                            self.skip_set.insert(read_id.clone());
                        }
                    }
                    if !(self.skip_set.contains(&read_id)
                        || self.pos_reads.contains_key(&read_id)
                        || self.neg_reads.contains_key(&read_id))
                    {
                        error!("didn't add failed read id to skip sets, likely a bug");
                    }
                    assert!(
                        self.skip_set.contains(&read_id)
                            || self.pos_mod_codes.contains_key(&read_id)
                            || self.neg_mod_codes.contains_key(&read_id)
                    );
                    self.add_mod_codes_for_record(
                        record,
                        pos_strand_mod_codes,
                        neg_strand_mod_codes,
                    );
                }
            }
        }
    }

    pub(crate) fn get_records_used_and_skipped(&self) -> (usize, usize) {
        let used = self
            .pos_reads
            .keys()
            .chain(self.neg_reads.keys())
            .collect::<HashSet<&String>>();
        let n_skipped = self.skip_set.len();
        (used.len(), n_skipped)
    }
}

pub(crate) struct DuplexReadCache<'a> {
    read_cache: ReadCache<'a>,
}

impl<'a> DuplexReadCache<'a> {
    pub(crate) fn new(
        method: Option<&'a CollapseMethod>,
        caller: &'a MultipleThresholdModCaller,
        edge_filter: Option<&'a EdgeFilter>,
        force_allow: bool,
    ) -> Self {
        let read_cache =
            ReadCache::new(method, caller, edge_filter, force_allow);

        Self { read_cache }
    }

    fn get_pos_strand_base_mod_call(
        &mut self,
        record: &bam::Record,
        position: u32,
        read_base: char,
    ) -> Option<BaseModCall> {
        if record.is_reverse() {
            match self.read_cache.get_mod_call(&record, position, read_base) {
                (_, Some(base_mod_call)) => Some(base_mod_call),
                _ => None,
            }
        } else {
            match self.read_cache.get_mod_call(&record, position, read_base) {
                (Some(base_mod_call), _) => Some(base_mod_call),
                _ => None,
            }
        }
    }

    fn get_neg_strand_base_mod_call(
        &mut self,
        record: &bam::Record,
        position: u32,
        read_base: char,
    ) -> Option<BaseModCall> {
        if record.is_reverse() {
            match self.read_cache.get_mod_call(&record, position, read_base) {
                (Some(base_mod_call), _) => Some(base_mod_call),
                _ => None,
            }
        } else {
            match self.read_cache.get_mod_call(&record, position, read_base) {
                (_, Some(base_mod_call)) => Some(base_mod_call),
                _ => None,
            }
        }
    }

    pub(crate) fn get_duplex_mod_call(
        &mut self,
        record: &bam::Record,
        position: u32,
        read_base: DnaBase,
        motif: &MotifLocations,
    ) -> Option<DuplexModCall> {
        let read_id = util::get_query_name_string(&record).ok()?;
        if self.read_cache.skip_set.contains(&read_id) {
            return None;
        }
        let (pos_base, neg_base) = if record.is_reverse() {
            (read_base.complement(), read_base)
        } else {
            (read_base, read_base.complement())
        };
        debug!("pos {} neg {}", pos_base.char(), neg_base.char());

        let pos_base_mod_call = self.get_pos_strand_base_mod_call(
            record,
            position,
            pos_base.char(),
        );
        if pos_base_mod_call.is_none() {
            debug!("no pos strand call");
        };
        let negative_position =
            motif.motif().negative_strand_position(position);

        if negative_position.is_none() {
            debug!("no call");
            return Some(DuplexModCall::NoCall {
                primary_base: read_base.char(),
            });
        }
        let negative_position = negative_position.unwrap();
        let neg_strand_base_mod_call = self.get_neg_strand_base_mod_call(
            record,
            negative_position,
            neg_base.char(), // double check this should be complemented
        );
        if neg_strand_base_mod_call.is_none() {
            debug!("negative stand mod call missing");
        }
        Some(DuplexModCall::from_base_mod_calls(
            pos_base_mod_call?,
            neg_strand_base_mod_call?,
            read_base.char(),
        ))
    }

    pub(crate) fn get_records_used_and_skipped(&self) -> (usize, usize) {
        self.read_cache.get_records_used_and_skipped()
    }
}

#[cfg(test)]
mod read_cache_tests {
    use std::collections::HashMap;

    use rust_htslib::bam::{self, FetchDefinition, Read, Reader as BamReader};

    use crate::mod_bam::{
        base_mod_probs_from_record, DeltaListConverter, SkipMode,
    };
    use crate::read_cache::ReadCache;
    use crate::test_utils::dna_complement;
    use crate::threshold_mod_caller::MultipleThresholdModCaller;
    use crate::util;

    fn tests_record(record: &bam::Record) {
        let query_name = String::from_utf8(record.qname().to_vec()).unwrap();
        let x = record.tid();
        assert!(x >= 0);

        // mapping of _reference position_ to forward read position
        let forward_aligned_pairs = util::get_aligned_pairs_forward(&record)
            .filter_map(|r| r.ok())
            .map(|(forward_q_pos, r_pos)| (r_pos, forward_q_pos))
            .collect::<HashMap<u64, usize>>();

        let forward_sequence = util::get_forward_sequence(&record)
            .map(|seq| seq.chars().collect::<Vec<char>>())
            .unwrap();

        let caller = MultipleThresholdModCaller::new_passthrough();
        let mut cache = ReadCache::new(None, &caller, None, false);
        cache.add_record(&record).unwrap();
        let converter =
            DeltaListConverter::new_from_record(&record, 'C').unwrap();
        let base_mod_probs =
            base_mod_probs_from_record(&record, &converter).unwrap();

        let (read_base_mod_probs, skip_mode) = cache
            .pos_reads
            .get(&query_name)
            .and_then(|base_to_calls| base_to_calls.get(&'C'))
            .unwrap();

        assert_eq!(skip_mode, &SkipMode::Ambiguous);

        // let read_base_mod_probs = match record.is_reverse() {
        //     true => cache
        //         .neg_reads
        //         .get(&query_name)
        //         .and_then(|base_to_calls| base_to_calls.get(&'C'))
        //         .unwrap(),
        //     false => cache
        //         .pos_reads
        //         .get(&query_name)
        //         .and_then(|base_to_calls| base_to_calls.get(&'C'))
        //         .unwrap(),
        // };

        // let read_base_mod_probs = cache
        //     .pos_reads
        //     .get(&query_name)
        //     .and_then(|base_to_calls| base_to_calls.get(&'C'))
        //     .unwrap();

        assert_eq!(
            base_mod_probs.pos_to_base_mod_probs.len(),
            read_base_mod_probs.len()
        );
        for (ref_pos, _probs) in read_base_mod_probs.iter() {
            assert!(ref_pos >= &0);
            let forward_read_pos = forward_aligned_pairs.get(ref_pos).unwrap();
            let read_base = forward_sequence[*forward_read_pos];
            assert_eq!(read_base, 'C');
        }
    }

    #[test]
    fn test_read_cache_aligned_pairs() {
        let mut reader =
            BamReader::from_path("tests/resources/fwd_rev_modbase_records.bam")
                .unwrap();
        for r in reader.records() {
            let record = r.unwrap();
            tests_record(&record);
        }
    }

    #[test]
    fn test_read_cache_stack_overflow_empty_tags() {
        let mut reader =
            BamReader::from_path("tests/resources/empty-tags.sorted.bam")
                .unwrap();

        let caller = MultipleThresholdModCaller::new_passthrough();
        let mut cache = ReadCache::new(None, &caller, None, false);
        for r in reader.records() {
            let record = r.unwrap();
            assert!(cache.add_record(&record).is_err());
        }
    }

    #[test]
    #[ignore = "verbose, used for development"]
    fn test_read_cache_get_mod_calls() {
        let mut reader = bam::IndexedReader::from_path(
            "tests/resources/fwd_rev_modbase_records.sorted.bam",
        )
        .unwrap();
        let header = reader.header().to_owned();
        let tid = 0;
        let target_name =
            String::from_utf8(header.tid2name(tid).to_vec()).unwrap();
        assert_eq!(target_name, "oligo_1512_adapters");
        let target_length = header.target_len(tid).unwrap();

        reader
            .fetch(FetchDefinition::Region(tid as i32, 0, target_length as i64))
            .unwrap();

        let caller = MultipleThresholdModCaller::new_passthrough();
        let mut read_cache = ReadCache::new(None, &caller, None, false);
        for p in reader.pileup() {
            let pileup = p.unwrap();
            for alignment in pileup.alignments() {
                if alignment.is_del() || alignment.is_refskip() {
                    continue;
                }
                let record = alignment.record();
                let read_base = if record.is_reverse() {
                    dna_complement(
                        record.seq()[alignment.qpos().unwrap()] as char,
                    )
                    .unwrap()
                } else {
                    record.seq()[alignment.qpos().unwrap()] as char
                };
                let mod_base_call =
                    read_cache.get_mod_call(&record, pileup.pos(), read_base);
                let read_id = String::from_utf8_lossy(record.qname());
                println!("{}\t{}\t{:?}", read_id, pileup.pos(), mod_base_call);
            }
        }
    }
}
