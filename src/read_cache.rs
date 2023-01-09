use crate::errs::RunError;
use crate::mod_bam::{
    base_mod_probs_from_record, extract_mod_probs, get_mm_tag_from_record, parse_raw_mod_tags,
    BaseModCall, BaseModPositions, DeltaListConverter, SeqPosBaseModProbs,
};
use crate::util;
use rust_htslib::bam;
use std::collections::{HashMap, HashSet};

/// Mapping of _reference position_ to base mod calls as determined by the aligned pairs for the
/// read
type ReadModBaseLookup = HashMap<u64, BaseModCall>; // todo use FxHasher

pub(crate) struct ReadCache {
    // todo last position (for gc)
    reads: HashMap<String, HashMap<char, ReadModBaseLookup>>,
    // these reads don't have mod tags
    skip_set: HashSet<String>,
}

impl ReadCache {
    // todo garbage collect freq
    pub(crate) fn new() -> Self {
        Self {
            reads: HashMap::new(),
            skip_set: HashSet::new(),
        }
    }

    fn add_base_mod_probs_for_base(
        &mut self,
        record_name: &str,
        record: &bam::Record,
        seq_pos_base_mod_probs: SeqPosBaseModProbs,
        canonical_base: char,
    ) -> Result<(), RunError> {
        let aligned_pairs =
            util::get_aligned_pairs_forward(&record).collect::<HashMap<usize, u64>>();
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

    fn add_record(&mut self, record: &bam::Record) -> Result<(), RunError> {
        let record_name = String::from_utf8(record.qname().to_vec())
            .map_err(|e| RunError::new_input_error(e.to_string()))?;
        match parse_raw_mod_tags(record) {
            Some(Ok((mm, ml))) => {
                let bases_with_mod_calls = mm
                    .split(';')
                    .filter_map(|raw_mm| {
                        if raw_mm.is_empty() {
                            None
                        } else {
                            Some(BaseModPositions::parse(raw_mm).map(|pos| pos.canonical_base))
                        }
                    })
                    .collect::<Result<Vec<char>, _>>()?;
                for canonical_base in bases_with_mod_calls {
                    let converter = DeltaListConverter::new_from_record(record, canonical_base)?;
                    let seq_pos_base_mod_probs =
                        extract_mod_probs(&mm, &ml, canonical_base, &converter)?;
                    self.add_base_mod_probs_for_base(
                        &record_name,
                        record,
                        seq_pos_base_mod_probs,
                        canonical_base,
                    )?;
                }
            }
            Some(Err(run_error)) => {
                return Err(run_error);
            }
            None => {
                // no mod tags, make a sentinel so we don't check again
                self.skip_set.insert(record_name);
            }
        }

        Ok(())
    }

    pub(crate) fn get_mod_call(
        &mut self,
        record: &bam::Record,
        position: u32,
        canonical_base: char,
        threshold: f32,
    ) -> Option<BaseModCall> {
        let read_id = String::from_utf8(record.qname().to_vec()).unwrap();
        if self.skip_set.contains(&read_id) {
            None
        } else {
            if let Some(canonical_base_to_calls) = self.reads.get(&read_id) {
                canonical_base_to_calls
                    .get(&canonical_base)
                    .and_then(|ref_pos_mod_base_calls| {
                        ref_pos_mod_base_calls
                            .get(&(position as u64))
                            .map(|base_mod_call| match base_mod_call {
                                BaseModCall::Canonical(p) | BaseModCall::Modified(p, _) => {
                                    if *p > threshold {
                                        *base_mod_call
                                    } else {
                                        BaseModCall::Filtered
                                    }
                                }
                                BaseModCall::Filtered => *base_mod_call,
                            })
                    })
            } else {
                match self.add_record(record) {
                    Ok(_) => {}
                    Err(_) => {
                        self.skip_set.insert(read_id);
                    }
                }
                self.get_mod_call(record, position, canonical_base, threshold)
            }
        }
    }
}

#[cfg(test)]
mod read_cache_tests {
    use crate::mod_bam::{base_mod_probs_from_record, BaseModCall, DeltaListConverter};
    use crate::read_cache::ReadCache;
    use crate::util;
    use crate::util::dna_complement;
    use rust_htslib::bam::{self, FetchDefinition, Read, Reader as BamReader};
    use std::collections::HashMap;

    fn tests_record(record: &bam::Record, header: &bam::HeaderView) {
        let query_name = String::from_utf8(record.qname().to_vec()).unwrap();
        let tid = {
            let x = record.tid();
            assert!(x >= 0);
            x as u32
        };

        // mapping of _reference position_ to forward read position
        let forward_aligned_pairs = util::get_aligned_pairs_forward(&record)
            .map(|(forward_q_pos, r_pos)| (r_pos, forward_q_pos))
            .collect::<HashMap<u64, usize>>();

        let forward_sequence = util::get_forward_sequence(&record)
            .map(|seq| seq.chars().collect::<Vec<char>>())
            .unwrap();

        // let reference_name = String::from_utf8_lossy(header.tid2name(tid));
        // let ref_len = header.target_len(tid).unwrap();
        // let seq = fasta
        //     .fetch_seq_string(reference_name, 0, ref_len as usize)
        //     .unwrap()
        //     .chars()
        //     .collect::<Vec<_>>();

        let mut cache = ReadCache::new();
        cache.add_record(&record).unwrap();
        let converter = DeltaListConverter::new_from_record(&record, 'C').unwrap();
        let base_mod_probs = base_mod_probs_from_record(&record, &converter, 'C').unwrap();
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
            BamReader::from_path("tests/resources/fwd_rev_modbase_records.bam").unwrap();
        let header = reader.header().to_owned();
        for r in reader.records() {
            let record = r.unwrap();
            tests_record(&record, &header);
        }
    }

    #[test]
    fn test_read_cache_get_mod_calls() {
        let mut reader =
            bam::IndexedReader::from_path("tests/resources/fwd_rev_modbase_records.sorted.bam")
                .unwrap();
        let header = reader.header().to_owned();
        let tid = 0;
        let target_name = String::from_utf8(header.tid2name(tid).to_vec()).unwrap();
        assert_eq!(target_name, "oligo_1512_adapters");
        let target_length = header.target_len(tid).unwrap();

        reader
            .fetch(FetchDefinition::Region(tid as i32, 0, target_length as i64))
            .unwrap();

        let mut read_cache = ReadCache::new();
        for p in reader.pileup() {
            let pileup = p.unwrap();
            for alignment in pileup.alignments() {
                if alignment.is_del() || alignment.is_refskip() {
                    continue;
                }
                let record = alignment.record();
                let read_base = if record.is_reverse() {
                    dna_complement(record.seq()[alignment.qpos().unwrap()] as char).unwrap()
                } else {
                    record.seq()[alignment.qpos().unwrap()] as char
                };
                let mod_base_call = read_cache.get_mod_call(&record, pileup.pos(), read_base, 0f32);
                let read_id = String::from_utf8_lossy(record.qname());
                println!("{}\t{}\t{:?}", read_id, pileup.pos(), mod_base_call);
            }
        }
    }
}
