use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap, VecDeque};
use std::marker::PhantomData;
use std::path::Path;

use rust_htslib::bam::{
    self, FetchDefinition, HeaderView, IndexedReader as IndexedBamReader, Read,
};
use rust_htslib::faidx::{self, Reader as FastaReader};

use crate::errs::{InputError, RunError};
use crate::mod_bam::{
    base_mod_probs_from_record, extract_mod_probs, BaseModCall, BaseModProbs, DeltaListConverter,
};
use crate::mod_base_code::{ModBaseCode, ModificationMotif};
use crate::util;

/// number of features
/// 0. N_canonical
/// 1. N_modified
/// 2. N_delete
/// 3. N_substitution
/// 4. N_filtered
/// 5. N_no_call
const NUM_FEATURES: usize = 6usize;

#[derive(Debug)]
enum Feature {
    Canonical,
    Modified(char),
    Delete,
    Substitution,
    Filtered,
    NoCall,
}

#[derive(Debug)]
pub struct PileupFeatureCounts {
    strand: util::Strand,
    position: u32,
    n_canonical: u16,
    n_modified: HashMap<char, u16>,
    n_delete: u16,
    n_substitution: u16,
    n_filtered: u16,
    n_nocall: u16,
}

impl PileupFeatureCounts {
    fn calc_total_modified_count(&self) -> u16 {
        self.n_modified.values().sum::<u16>()
    }

    pub(crate) fn coverage(&self) -> u16 {
        let n_modified = self.calc_total_modified_count();
        n_modified + self.n_canonical + self.n_filtered + self.n_substitution + self.n_delete
    }

    pub(crate) fn percent_modified(&self, mod_code: char) -> f32 {
        let n_modified = self.calc_n_modified(mod_code);
        let n_modified = n_modified as f32;
        let denom = self.n_canonical as f32 + self.calc_total_modified_count() as f32;
        let percent_mod = 100f32 * (n_modified / denom);
        percent_mod
    }

    pub(crate) fn calc_n_modified(&self, mod_code: char) -> u16 {
        *self.n_modified.get(&mod_code).unwrap_or(&0u16)
    }

    pub(crate) fn score(&self, mod_code: char) -> u32 {
        let n_modified = self.calc_n_modified(mod_code);
        let base_score_numerator = n_modified + self.n_canonical;
        let base_score_denominator =
            n_modified + self.n_canonical + self.n_substitution + self.n_delete + self.n_filtered;
        let score = base_score_numerator as f32 / base_score_denominator as f32;
        if score.is_nan() {
            0u32
        } else {
            (1000f32 * score).floor() as u32
        }
    }
}

struct FeatureVector {
    matrix: Vec<u16>,
    positions: Vec<Option<(u32, util::Strand)>>,
    feature_size: usize,
}

impl FeatureVector {
    fn new(chunk_size: usize, mod_base_code: &dyn ModBaseCode) -> Self {
        let feature_size = NUM_FEATURES + mod_base_code.num_mods() - 1;
        let matrix = vec![0u16; chunk_size * feature_size];
        let positions = vec![None; chunk_size];

        Self {
            matrix,
            positions,
            feature_size,
        }
    }

    #[inline(always)]
    fn feature_index(&self, feature: &Feature, mod_base_code: &dyn ModBaseCode) -> usize {
        match feature {
            Feature::Canonical => 0,
            Feature::Delete => 1,
            Feature::Substitution => 2,
            Feature::Filtered => 3,
            Feature::NoCall => 4,
            Feature::Modified(raw_code) => 5 + mod_base_code.idx_for_mod_code(*raw_code).unwrap(),
        }
    }

    fn init_position(&mut self, idx: usize, strand: util::Strand, position: u32) {
        self.positions[idx] = Some((position, strand));
    }

    fn add_feature_to_position(
        &mut self,
        idx: usize,
        position: u32,
        strand: util::Strand,
        feature: Feature,
        mod_base_code: &dyn ModBaseCode,
    ) {
        assert!(idx < self.positions.len());
        if self.positions[idx].is_none() {
            self.positions[idx] = Some((position, strand));
        } else {
            assert_eq!(self.positions[idx], Some((position, strand)));
        }

        let offset = idx * self.feature_size + self.feature_index(&feature, mod_base_code);
        assert!(offset < self.matrix.len());
        self.matrix[offset] += 1;
    }

    pub fn iter_counts<'a>(
        &'a self,
        mod_base_code: &'a dyn ModBaseCode,
    ) -> impl Iterator<Item = PileupFeatureCounts> + 'a {
        self.matrix
            .chunks(self.feature_size)
            .zip(self.positions.iter())
            .filter_map(|(counts, position_strand)| {
                position_strand.map(|(pos, strand)| {
                    let n_canonical = counts[0];
                    let n_delete = counts[1];
                    let n_substitution = counts[2];
                    let n_filtered = counts[3];
                    let n_nocall = counts[4];
                    let mut n_modified = HashMap::<char, u16>::new();
                    for raw_code in mod_base_code.raw_mod_codes() {
                        let idx = self.feature_index(&Feature::Modified(*raw_code), mod_base_code);
                        assert!(idx < counts.len());
                        let count_modified = counts[idx];
                        n_modified.insert(*raw_code, count_modified);
                    }
                    PileupFeatureCounts {
                        strand,
                        position: pos,
                        n_canonical,
                        n_modified,
                        n_delete,
                        n_substitution,
                        n_filtered,
                        n_nocall,
                    }
                })
            })
    }
}

enum MotifPileup {
    CoveredMotif(bam::pileup::Pileup),
    NoCoverage(u32),
}

struct MotifPileupIter<'a> {
    pileups: bam::pileup::Pileups<'a, bam::IndexedReader>,
    motif_positions: VecDeque<(usize, util::Strand)>,
    start_pos: u32,
    end_pos: u32,
    next_pos: Option<(u32, util::Strand)>,
}

impl<'a> MotifPileupIter<'a> {
    fn new(
        pileups: bam::pileup::Pileups<'a, bam::IndexedReader>,
        start_pos: u32,
        end_pos: u32,
        mut motif_positions: VecDeque<(usize, util::Strand)>,
    ) -> Self {
        let next_pos = if motif_positions.is_empty() {
            None
        } else {
            motif_positions
                .pop_front()
                .map(|(pos, strand)| (pos as u32, strand))
        };
        Self {
            pileups,
            start_pos,
            end_pos,
            motif_positions,
            next_pos,
        }
    }

    fn update_next_position(&mut self) {
        self.next_pos = self
            .motif_positions
            .pop_front()
            .map(|(pos, strand)| (pos as u32, strand));
    }
}

impl<'a> Iterator for MotifPileupIter<'a> {
    type Item = (util::Strand, MotifPileup);

    fn next(&mut self) -> Option<Self::Item> {
        if let Some((next_pos, strand)) = self.next_pos {
            // dbg!(next_pos, strand);
            let mut pileup: Option<Self::Item> = None;
            while let Some(plp) = self.pileups.next().map(|res| res.unwrap()) {
                // todo: document how end_pos works, [start_pos, end_pos)
                let off_end = plp.pos() > self.end_pos;
                if off_end {
                    return None;
                }
                match plp.pos().cmp(&next_pos) {
                    Ordering::Equal => {
                        pileup = Some((strand, MotifPileup::CoveredMotif(plp)));
                        self.update_next_position();
                        break;
                    }
                    Ordering::Less => continue,
                    Ordering::Greater => {
                        println!(">greater {next_pos}");
                        pileup = Some((strand, MotifPileup::NoCoverage(next_pos)));
                        self.update_next_position();
                        break;
                    }
                }
            }
            return pileup;
        } else {
            None
        }
    }
}

pub struct ModBasePileup {
    chrom_name: String,
    start_pos: u32,
    features: FeatureVector,
}

impl ModBasePileup {
    pub fn decode(self, mod_base_code: &dyn ModBaseCode, sep: char) -> Option<(String, u64)> {
        let mut decoded = String::new();
        let mut n = 0u64;
        for counts in self.features.iter_counts(mod_base_code) {
            for raw_code in mod_base_code.raw_mod_codes() {
                let row = format!(
                    "\
                {}{sep}{}{sep}{}{sep}{}{sep}{}{sep}{}{sep}{}{sep}{}\
                {sep}{}{sep}{}{sep}{}{sep}{}{sep}{}{sep}{}{sep}{}\
                {sep}{}{sep}{}\n",
                    self.chrom_name,                                      // 1
                    counts.position,                                      // 2
                    counts.position + 1,                                  // 3
                    raw_code,                                             // 4
                    counts.score(*raw_code),                              // 5
                    counts.strand.to_char(),                              // 6
                    counts.position,                                      // 7
                    counts.position + 1,                                  // 8
                    "0,0,0",                                              // 9
                    counts.coverage(),                                    // 10
                    format!("{:.2}", counts.percent_modified(*raw_code)), // 11
                    counts.n_canonical,                                   // 12
                    counts.calc_n_modified(*raw_code),                    // 13
                    counts.n_filtered,                                    // 14
                    counts.n_substitution,                                // 15
                    counts.n_delete,                                      // 16
                    counts.n_nocall,                                      // 17
                );
                decoded.push_str(&row);
                n += 1;
            }
        }
        if decoded.is_empty() {
            assert_eq!(n, 0);
            None
        } else {
            Some((decoded, n))
        }
    }
}

pub struct ModBasePileupProcessor<T: AsRef<Path>> {
    bam_fp: T,
    chrom_tid: u32,
    start_pos: u32,
    end_pos: u32,
    read_cache: ReadCache,
}

impl<T: AsRef<Path>> ModBasePileupProcessor<T> {
    pub fn new(bam_fp: T, chrom_tid: u32, start_pos: u32, end_pos: u32) -> Result<Self, String> {
        if end_pos <= start_pos {
            Err("end_pos needs to be larger (after) start_pos".to_owned())
        } else {
            let read_cache = ReadCache::new();
            Ok(Self {
                bam_fp,
                chrom_tid,
                start_pos,
                end_pos,
                read_cache,
            })
        }
    }

    pub fn process_region(
        &mut self,
        mod_base_code: &dyn ModBaseCode,
        modification_motif: &dyn ModificationMotif,
        ref_name: &str,
        ref_seq: &str,
    ) -> ModBasePileup {
        let mut bam_reader = IndexedBamReader::from_path(&self.bam_fp).unwrap();
        let mut motif_positions = modification_motif
            .find_matches(ref_seq)
            .into_iter()
            .map(|(pos, strand)| (pos + (self.start_pos as usize), strand))
            .collect::<VecDeque<(usize, util::Strand)>>();

        let reference_seq = ref_seq.chars().collect::<Vec<char>>();
        bam_reader
            .fetch(FetchDefinition::Region(
                self.chrom_tid as i32,
                self.start_pos as i64,
                self.end_pos as i64,
            ))
            .unwrap();

        let mut features = FeatureVector::new(motif_positions.len(), mod_base_code);
        // todo log out where the motifs are
        // if ref_name == "oligo_1512_adapters" {
        //     eprintln!(
        //         "> ref name {ref_name}, start={}, end={}",
        //         self.start_pos, self.end_pos
        //     );
        //     for pos in motif_positions.iter() {
        //         dbg!(pos);
        //     }
        // }
        let pileup_iter = MotifPileupIter::new(
            bam_reader.pileup(),
            self.start_pos,
            self.end_pos,
            motif_positions,
        );
        for (idx, (strand, motif_pileup)) in pileup_iter.enumerate() {
            match motif_pileup {
                MotifPileup::CoveredMotif(pileup) => {
                    features.init_position(idx, strand, pileup.pos());
                    // println!(">pos {}", pileup.pos());
                    let adjusted_pos = (pileup.pos() - self.start_pos) as usize;
                    assert_eq!(pileup.tid(), self.chrom_tid);
                    assert!(
                        pileup.pos() >= self.start_pos,
                        "pileup pos too small {}, {}",
                        pileup.pos(),
                        self.start_pos
                    );
                    assert!(
                        pileup.pos() < self.end_pos + 1,
                        "pileup pos too large {}, {}",
                        pileup.pos(),
                        self.end_pos
                    );
                    for alignment in pileup.alignments() {
                        // this is a hack because we only allow calling on the "top" a.k.a. the
                        // read strand, in duplex mode we'd have to remove this logic and account
                        // for the call in the get_mod_call method below
                        match (strand, alignment.record().is_reverse()) {
                            (util::Strand::Positive, true) => {
                                continue;
                            }
                            (util::Strand::Negative, false) => {
                                continue;
                            }
                            _ => {}
                        }
                        let feature = if alignment.is_del() || alignment.is_refskip() {
                            Feature::Delete
                        } else {
                            let record = alignment.record();
                            // n.b. this can be the reverse sequence
                            let read_base = record.seq()[alignment.qpos().unwrap()] as char;
                            let ref_base = reference_seq[adjusted_pos];
                            if read_base != ref_base {
                                // println!(
                                //     "> sub ref={} > read={} ({})",
                                //     ref_base,
                                //     read_base,
                                //     alignment.qpos().unwrap()
                                // );
                                Feature::Substitution
                            } else {
                                self.read_cache
                                    .get_mod_call(
                                        // todo stand should be included here for duplex
                                        &record,
                                        pileup.pos(),
                                        modification_motif.canonical_base(),
                                        0.0f32,
                                    )
                                    // todo this should really be a switch based on the mod-calling mode,
                                    //  for the typical `?` mode, when there is no modification call, we
                                    //  do not default to canonical
                                    .unwrap_or(Feature::NoCall)
                            }
                        };
                        features.add_feature_to_position(
                            idx,
                            pileup.pos(),
                            strand,
                            feature,
                            mod_base_code,
                        );
                    }
                }
                MotifPileup::NoCoverage(pos) => {
                    // println!("> no coverage {pos}");
                    features.add_feature_to_position(
                        idx,
                        pos,
                        strand,
                        Feature::NoCall,
                        mod_base_code,
                    );
                }
            }
        }

        ModBasePileup {
            chrom_name: ref_name.to_owned(),
            start_pos: self.start_pos,
            features,
        }
    }
}

pub fn estimate_quantile_threshold() {
    unimplemented!()
}

/// Mapping of _reference position_ to base mod calls as determined by the aligned pairs for the
/// read
type ReadModBaseLookup = HashMap<u64, BaseModCall>;

struct ReadCache {
    // todo last position (for gc)
    reads: HashMap<String, ReadModBaseLookup>,
}

impl ReadCache {
    // todo garbage collect freq
    fn new() -> Self {
        Self {
            reads: HashMap::new(),
        }
    }

    fn add_read(&mut self, record: &bam::Record, canonical_base: char) -> Result<(), RunError> {
        let converter = DeltaListConverter::new_from_record(record, canonical_base)?;
        // todo need collapse option in here
        let base_mod_probs = base_mod_probs_from_record(&record, &converter, canonical_base)?;
        let aligned_pairs =
            util::get_aligned_pairs_forward(&record).collect::<HashMap<usize, u64>>();
        let base_mod_probs = base_mod_probs
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
        let record_name = String::from_utf8(record.qname().to_vec())
            .map_err(|e| RunError::new_input_error(e.to_string()))?;
        self.reads.insert(record_name, base_mod_probs);

        Ok(())
    }

    fn get_mod_call(
        &mut self,
        record: &bam::Record,
        position: u32,
        canonical_base: char,
        threshold: f32,
    ) -> Option<Feature> {
        let read_id = String::from_utf8(record.qname().to_vec()).unwrap();
        if let Some(read_mod_lookup) = self.reads.get(&read_id) {
            read_mod_lookup
                .get(&(position as u64))
                .map(|base_mod_call| match base_mod_call {
                    BaseModCall::Canonical(p) => {
                        if *p > threshold {
                            Feature::Canonical
                        } else {
                            Feature::Filtered
                        }
                    }
                    BaseModCall::Modified(p, mod_code) => {
                        if *p > threshold {
                            Feature::Modified(*mod_code)
                        } else {
                            Feature::Filtered
                        }
                    }
                })
        } else {
            self.add_read(record, canonical_base).unwrap();
            self.get_mod_call(record, position, canonical_base, threshold)
        }
    }
}

#[cfg(test)]
mod mod_pileup_tests {
    use std::collections::HashMap;

    use crate::interval_chunks::IntervalChunks;
    use rust_htslib::bam::{self, HeaderView, Read, Reader as BamReader};
    use rust_htslib::faidx::{self, Reader as FastaReader};

    use crate::mod_bam::{base_mod_probs_from_record, DeltaListConverter};
    use crate::mod_base_code::{
        CHHModificationMotif, CpGModificationMotif, CytosineModBaseCode, HydroxyMethylCytosineCode,
        ModBaseCode, ModificationMotif,
    };
    use crate::mod_pileup::{
        Feature, FeatureVector, ModBasePileupProcessor, PileupFeatureCounts, ReadCache,
        NUM_FEATURES,
    };
    use crate::util;
    use crate::util::Strand;

    fn load_test_sequence(name: &str) -> String {
        let fasta_fp = "tests/resources/CGI_ladder_3.6kb_ref.fa";
        let mut fasta_reader = faidx::Reader::from_path(fasta_fp).unwrap();
        let dna = fasta_reader.fetch_seq_string(name, 0, 156).unwrap();
        dna
    }

    fn tests_record(record: &bam::Record, fasta: &mut FastaReader, header: &HeaderView) {
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

        let reference_name = String::from_utf8_lossy(header.tid2name(tid));
        let ref_len = header.target_len(tid).unwrap();
        let seq = fasta
            .fetch_seq_string(reference_name, 0, ref_len as usize)
            .unwrap()
            .chars()
            .collect::<Vec<_>>();

        let mut cache = ReadCache::new();
        cache.add_read(&record, 'C').unwrap();
        let converter = DeltaListConverter::new_from_record(&record, 'C').unwrap();
        let base_mod_probs = base_mod_probs_from_record(&record, &converter, 'C').unwrap();
        let read_base_mod_probs = cache.reads.get(&query_name).unwrap();

        assert_eq!(base_mod_probs.len(), read_base_mod_probs.len());

        for (ref_pos, _probs) in read_base_mod_probs.iter() {
            assert!(ref_pos >= &0);
            let forward_read_pos = forward_aligned_pairs.get(ref_pos).unwrap();
            let read_base = forward_sequence[*forward_read_pos];
            assert_eq!(read_base, 'C');
        }
    }

    #[test]
    fn test_read_cache() {
        let mut fasta_reader =
            FastaReader::from_path("tests/resources/CGI_ladder_3.6kb_ref.fa").unwrap();
        let mut reader =
            BamReader::from_path("tests/resources/fwd_rev_modbase_records.bam").unwrap();
        let header = reader.header().to_owned();
        for r in reader.records() {
            let record = r.unwrap();
            tests_record(&record, &mut fasta_reader, &header);
        }
    }

    #[test]
    fn test_feature_vector() {
        let code = if 1 > 10 {
            // never true, just to show dynamism
            &CytosineModBaseCode as &dyn ModBaseCode
        } else {
            &HydroxyMethylCytosineCode as &dyn ModBaseCode
        };

        let chunk_size = 4;
        let num_features = NUM_FEATURES + code.num_mods() - 1;
        let mut feature_vector = FeatureVector::new(chunk_size, code);
        assert_eq!(feature_vector.feature_size, num_features);
        assert_eq!(feature_vector.matrix.len(), chunk_size * num_features);

        let seq = "ACGCGTT";
        let motif = &CpGModificationMotif;
        let mut motif_positions = motif.find_matches(seq);
        assert_eq!(motif_positions.len(), chunk_size); // tested in module

        let (pos, strand) = motif_positions.pop_front().unwrap();
        assert_eq!(pos, 1);
        assert_eq!(strand, Strand::Positive);
        feature_vector.add_feature_to_position(0, pos as u32, strand, Feature::Modified('h'), code);
        feature_vector.add_feature_to_position(0, pos as u32, strand, Feature::Canonical, code);
        feature_vector.add_feature_to_position(0, pos as u32, strand, Feature::Filtered, code);

        let counts = feature_vector.iter_counts(code).collect::<Vec<_>>();
        assert_eq!(counts.len(), 1);
        let counts = counts.iter().next().unwrap();
        assert_eq!(counts.n_canonical, 1);
        assert_eq!(counts.n_filtered, 1);
        assert_eq!(counts.position, pos as u32);
        assert_eq!(counts.n_delete, 0);
        assert_eq!(counts.n_modified.get(&'h').unwrap(), &1u16);
    }

    #[test]
    fn test_check_sequence_slicing_is_same_as_fetch() {
        let fasta_fp = "tests/resources/CGI_ladder_3.6kb_ref.fa";
        let mut fasta_reader = faidx::Reader::from_path(fasta_fp).unwrap();
        let name = "oligo_1512_adapters";
        let dna = load_test_sequence(name);
        let start = 49;
        let end = 99;
        let slice_a = util::slice_dna_sequence(&dna, start, end);
        let slice_b = fasta_reader.fetch_seq_string(name, start, end).unwrap();
        assert_eq!(slice_a, slice_b);
    }

    #[test]
    fn test_mod_pileup_processor() {
        let name = "oligo_1512_adapters";
        let dna = load_test_sequence(name);
        let code = &HydroxyMethylCytosineCode;
        let start_pos = 70;
        let end_pos = 95;
        let mut processor = ModBasePileupProcessor::new(
            "tests/resources/fwd_rev_modbase_records.sorted.bam",
            0,
            start_pos,
            end_pos,
        )
        .unwrap();
        let ref_seq = util::slice_dna_sequence(&dna, start_pos as usize, end_pos as usize);
        let pileup = processor.process_region(code, &CpGModificationMotif, name, &ref_seq);
        let counts = pileup
            .features
            .iter_counts(code)
            .map(|pileup| {
                let key = (pileup.position, pileup.strand);
                (key, pileup)
            })
            .collect::<HashMap<_, _>>();

        assert_eq!(counts.len(), 6);
        let pileup = counts.get(&(72, Strand::Positive)).unwrap();
        assert_eq!(pileup.n_modified.get(&'m').unwrap(), &1u16);
        assert_eq!(pileup.n_modified.get(&'h').unwrap(), &0u16);
        assert_eq!(pileup.n_canonical, 0u16);
        let pileup = counts.get(&(73, Strand::Negative)).unwrap();
        assert_eq!(pileup.n_modified.get(&'m').unwrap(), &1u16);
        assert_eq!(pileup.n_modified.get(&'h').unwrap(), &0u16);
        assert_eq!(pileup.n_canonical, 0u16);
        let pileup = counts.get(&(90, Strand::Positive)).unwrap();
        assert_eq!(pileup.n_modified.get(&'m').unwrap(), &1u16);
        assert_eq!(pileup.n_modified.get(&'h').unwrap(), &0u16);
        assert_eq!(pileup.n_canonical, 0u16);
        let pileup = counts.get(&(91, Strand::Negative)).unwrap();
        assert_eq!(pileup.n_modified.get(&'m').unwrap(), &0u16);
        assert_eq!(pileup.n_modified.get(&'h').unwrap(), &0u16);
        assert_eq!(pileup.n_canonical, 0u16);
        assert_eq!(pileup.n_nocall, 1u16);
        let pileup = counts.get(&(93, Strand::Positive)).unwrap();
        assert_eq!(pileup.n_modified.get(&'m').unwrap(), &1u16);
        assert_eq!(pileup.n_modified.get(&'h').unwrap(), &0u16);
        assert_eq!(pileup.n_canonical, 0u16);
        let pileup = counts.get(&(94, Strand::Negative)).unwrap();
        assert_eq!(pileup.n_modified.get(&'m').unwrap(), &0u16);
        assert_eq!(pileup.n_modified.get(&'h').unwrap(), &0u16);
        assert_eq!(pileup.n_canonical, 0u16);
        assert_eq!(pileup.n_nocall, 0u16);
        assert_eq!(pileup.n_delete, 1u16);
    }
}
