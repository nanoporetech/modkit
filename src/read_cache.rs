use std::collections::{HashMap, HashSet};

use log::debug;
use rust_htslib::bam;

use crate::errs::RunError;
use crate::mod_bam::{
    collapse_mod_probs, BaseModCall, CollapseMethod, ModBaseInfo,
    SeqPosBaseModProbs,
};
use crate::mod_base_code::ModCode;
use crate::util;

/// Mapping of _reference position_ to base mod calls as determined by the aligned pairs for the
/// read
type RefPosBaseModCalls = HashMap<u64, BaseModCall>; // todo use FxHasher

// todo last position (for gc)
pub(crate) struct ReadCache<'a> {
    /// Mapping of read_id to reference position <> base mod calls for that read
    reads: HashMap<String, HashMap<char, RefPosBaseModCalls>>,
    /// these reads don't have mod tags or should be skipped for some other reason
    skip_set: HashSet<String>,
    /// mapping of read_id (query_name) to the mod codes contained in that read
    mod_codes: HashMap<String, HashSet<ModCode>>,
    /// collapse method
    method: Option<&'a CollapseMethod>,
}

impl<'a> ReadCache<'a> {
    // todo garbage collect freq
    pub(crate) fn new(method: Option<&'a CollapseMethod>) -> Self {
        Self {
            reads: HashMap::new(),
            skip_set: HashSet::new(),
            mod_codes: HashMap::new(),
            method,
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
        canonical_base: char,
    ) -> Result<(), RunError> {
        let aligned_pairs = util::get_aligned_pairs_forward(&record)
            .collect::<HashMap<usize, u64>>();
        let ref_pos_base_mod_calls = seq_pos_base_mod_probs
            .into_iter()
            // here the q_pos is the forward-oriented position
            .flat_map(|(q_pos, bmp)| {
                if let Some(r_pos) = aligned_pairs.get(&q_pos) {
                    Some((*r_pos, bmp.base_mod_call()))
                } else {
                    None
                }
            })
            .collect::<HashMap<u64, BaseModCall>>();
        let bases_to_mod_calls = self
            .reads
            .entry(record_name.to_owned())
            .or_insert(HashMap::new());
        bases_to_mod_calls.insert(canonical_base, ref_pos_base_mod_calls);
        Ok(())
    }

    /// Add a record to the cache.
    fn add_record(&mut self, record: &bam::Record) -> Result<(), RunError> {
        let record_name = util::get_query_name_string(record)
            .map_err(|e| RunError::new_input_error(e.to_string()))?;
        let mod_base_info = ModBaseInfo::new_from_record(record)?;

        if mod_base_info.is_empty() {
            let msg = format!("record {} has no mod calls", &record_name);
            debug!("{}", &msg);
            return Err(RunError::Skipped(msg));
        }

        for (converter, mut seq_base_mod_probs) in
            mod_base_info.into_iter_base_mod_probs()
        {
            if let Some(method) = &self.method {
                seq_base_mod_probs =
                    collapse_mod_probs(seq_base_mod_probs, method);
            }

            let mod_codes = seq_base_mod_probs
                .values()
                .flat_map(|base_mod_probs| {
                    base_mod_probs.mod_codes.iter().filter_map(|raw_mod_code| {
                        ModCode::parse_raw_mod_code(*raw_mod_code).ok()
                    })
                })
                .collect::<HashSet<ModCode>>();

            let record_mod_codes = self
                .mod_codes
                .entry(record_name.to_owned())
                .or_insert(HashSet::new());
            record_mod_codes.extend(mod_codes);

            self.add_modbase_probs_for_record_and_canonical_base(
                &record_name,
                record,
                seq_base_mod_probs,
                converter.canonical_base,
            )?;
        }
        Ok(())
    }

    /// Get the mod call for a reference position from a read. If this read is
    /// in the cache, look it up, if not parse the tags, add it to the cache
    /// and return the mod call (if present).
    pub(crate) fn get_mod_call(
        &mut self,
        record: &bam::Record,
        position: u32,
        canonical_base: char, // todo make this DnaBase
        threshold: f32,
    ) -> Option<BaseModCall> {
        let read_id = String::from_utf8(record.qname().to_vec()).unwrap();
        if self.skip_set.contains(&read_id) {
            None
        } else {
            if let Some(canonical_base_to_calls) = self.reads.get(&read_id) {
                // todo(arand) this is ugly, make it easier to follow or at least comment up the
                //  logic also maybe remove the need for the read cache all together?
                canonical_base_to_calls.get(&canonical_base).and_then(
                    |ref_pos_mod_base_calls| {
                        ref_pos_mod_base_calls.get(&(position as u64)).map(
                            |base_mod_call| match base_mod_call {
                                BaseModCall::Canonical(p)
                                | BaseModCall::Modified(p, _) => {
                                    if *p > threshold {
                                        *base_mod_call
                                    } else {
                                        BaseModCall::Filtered
                                    }
                                }
                                BaseModCall::Filtered => *base_mod_call,
                            },
                        )
                    },
                )
            } else {
                match self.add_record(record) {
                    Ok(_) => {}
                    Err(err) => {
                        debug!(
                            "read {read_id} failed to get mod tags {}",
                            err.to_string()
                        );
                        self.skip_set.insert(read_id);
                    }
                }
                self.get_mod_call(record, position, canonical_base, threshold)
            }
        }
    }

    pub(crate) fn get_mod_codes_for_record(
        &mut self,
        record: &bam::Record,
    ) -> HashSet<ModCode> {
        let read_id = String::from_utf8(record.qname().to_vec()).unwrap();
        if self.skip_set.contains(&read_id) {
            HashSet::new()
        } else {
            if let Some(mod_codes) = self.mod_codes.get(&read_id) {
                mod_codes.iter().map(|mc| *mc).collect()
            } else {
                match self.add_record(record) {
                    Ok(_) => {}
                    Err(err) => {
                        debug!(
                            "read {read_id} failed to get mod tags {}",
                            err.to_string()
                        );
                        self.skip_set.insert(read_id.clone());
                    }
                }
                assert!(
                    self.skip_set.contains(&read_id)
                        || self.mod_codes.contains_key(&read_id),
                    "should have been added..."
                );
                self.get_mod_codes_for_record(record)
            }
        }
    }
}

#[cfg(test)]
mod read_cache_tests {
    use std::collections::HashMap;

    use rust_htslib::bam::{self, FetchDefinition, Read, Reader as BamReader};

    use crate::mod_bam::{base_mod_probs_from_record, DeltaListConverter};
    use crate::read_cache::ReadCache;
    use crate::test_utils::dna_complement;
    use crate::util;

    fn tests_record(record: &bam::Record) {
        let query_name = String::from_utf8(record.qname().to_vec()).unwrap();
        let x = record.tid();
        assert!(x >= 0);

        // mapping of _reference position_ to forward read position
        let forward_aligned_pairs = util::get_aligned_pairs_forward(&record)
            .map(|(forward_q_pos, r_pos)| (r_pos, forward_q_pos))
            .collect::<HashMap<u64, usize>>();

        let forward_sequence = util::get_forward_sequence(&record)
            .map(|seq| seq.chars().collect::<Vec<char>>())
            .unwrap();

        let mut cache = ReadCache::new(None);
        cache.add_record(&record).unwrap();
        let converter =
            DeltaListConverter::new_from_record(&record, 'C').unwrap();
        let base_mod_probs =
            base_mod_probs_from_record(&record, &converter).unwrap();
        let read_base_mod_probs = cache
            .reads
            .get(&query_name)
            .and_then(|base_to_calls| base_to_calls.get(&'C'))
            .unwrap();

        assert_eq!(base_mod_probs.len(), read_base_mod_probs.len());
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
            BamReader::from_path("tests/resources/empty-tags.bam").unwrap();

        let mut cache = ReadCache::new(None);
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

        let mut read_cache = ReadCache::new(None);
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
                let mod_base_call = read_cache.get_mod_call(
                    &record,
                    pileup.pos(),
                    read_base,
                    0f32,
                );
                let read_id = String::from_utf8_lossy(record.qname());
                println!("{}\t{}\t{:?}", read_id, pileup.pos(), mod_base_call);
            }
        }
    }
}
