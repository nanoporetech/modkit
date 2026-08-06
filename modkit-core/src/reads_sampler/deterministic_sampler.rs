use anyhow::{bail, Result};
use rust_htslib::bam;

/// Domain separation and compatibility version for the stable alignment
/// sampler. Changing the identity fields, byte encoding, digest, or threshold
/// rule requires a new version.
const SAMPLER_V1_DOMAIN: &[u8] = b"modkit-alignment-fraction-sampler-v1\0";

/// Alignment-identity flag bits: reverse, read1, read2, secondary, and
/// supplementary. Mutable annotation flags such as duplicate, QC failure, and
/// proper pair are intentionally excluded.
const ALIGNMENT_IDENTITY_FLAG_MASK: u16 = 0x09d0;
const DRAW_SPACE: f64 = (1u64 << 53) as f64;
const FNV1A_OFFSET_BASIS: u64 = 0xcbf2_9ce4_8422_2325;
const FNV1A_PRIME: u64 = 0x0000_0100_0000_01b3;

/// The v1 portable hash is FNV-1a-64 followed by the SplitMix64 finalizer.
/// These operations and constants are part of the sampler compatibility
/// contract, rather than an implementation-selected `Hash` algorithm.
struct StableV1Hasher(u64);

impl StableV1Hasher {
    fn new(master_seed: u64) -> Self {
        let mut hasher = Self(FNV1A_OFFSET_BASIS);
        hasher.update(SAMPLER_V1_DOMAIN);
        hasher.update(&master_seed.to_le_bytes());
        hasher
    }

    #[inline]
    fn update(&mut self, bytes: &[u8]) {
        for byte in bytes {
            self.0 ^= *byte as u64;
            self.0 = self.0.wrapping_mul(FNV1A_PRIME);
        }
    }

    fn finish(mut self) -> u64 {
        self.0 ^= self.0 >> 30;
        self.0 = self.0.wrapping_mul(0xbf58_476d_1ce4_e5b9);
        self.0 ^= self.0 >> 27;
        self.0 = self.0.wrapping_mul(0x94d0_49bb_1331_11eb);
        self.0 ^ (self.0 >> 31)
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum InclusionThreshold {
    None,
    Below(u64),
    All,
}

/// A versioned, stateless Bernoulli sampler for BAM alignment records.
///
/// The v1 sampling unit is one alignment record, not one molecule or QNAME.
/// Its stable identity consists of the raw QNAME, numeric reference ID,
/// alignment start, canonical CIGAR operations, and the masked flag bits in
/// [`ALIGNMENT_IDENTITY_FLAG_MASK`]. This distinguishes ordinary primary,
/// secondary, supplementary, strand, and read-end alignments sharing a QNAME,
/// while ensuring an alignment refetched in multiple genomic intervals gets
/// one decision. Byte-identical duplicate alignment records intentionally
/// share a decision because BAM records do not expose a portable stable
/// ordinal. Numeric reference IDs make this contract stable for refetches of
/// the same BAM/CRAM header; reproducibility after reference-header reordering
/// is not promised.
///
/// MAPQ, mate fields, SEQ, QUAL, and auxiliary tags (including MM/ML) are not
/// identity. They are mutable annotations or measured content and must not
/// influence whether the observation is sampled.
#[derive(Clone, Copy, Debug)]
pub(crate) struct DeterministicFractionSampler {
    master_seed: u64,
    threshold: InclusionThreshold,
}

impl DeterministicFractionSampler {
    pub(crate) fn new(master_seed: u64, fraction: f64) -> Result<Self> {
        if !fraction.is_finite() || !(0.0..=1.0).contains(&fraction) {
            bail!(
                "sampling fraction must be a finite number in the inclusive \
                 range [0, 1]; got '{fraction}'"
            );
        }
        let threshold = if fraction == 0.0 {
            InclusionThreshold::None
        } else if fraction == 1.0 {
            InclusionThreshold::All
        } else {
            InclusionThreshold::Below((fraction * DRAW_SPACE).floor() as u64)
        };
        Ok(Self { master_seed, threshold })
    }

    #[inline]
    pub(crate) fn include(&self, record: &bam::Record) -> bool {
        match self.threshold {
            InclusionThreshold::None => false,
            InclusionThreshold::All => true,
            InclusionThreshold::Below(threshold) => {
                (self.alignment_score(record) >> 11) < threshold
            }
        }
    }

    /// Return the v1 score used by the inclusion threshold. Integer fields use
    /// fixed-width little-endian encoding (signed fields use their two's
    /// complement bytes), and variable-length fields have an explicit u64
    /// length. Inclusion compares the score's top 53 bits with
    /// `floor(fraction * 2^53)`.
    pub(crate) fn alignment_score(&self, record: &bam::Record) -> u64 {
        let mut hash = StableV1Hasher::new(self.master_seed);

        let qname = record.qname();
        hash.update(&(qname.len() as u64).to_le_bytes());
        hash.update(qname);
        hash.update(&record.tid().to_le_bytes());
        hash.update(&record.pos().to_le_bytes());
        hash.update(
            &(record.flags() & ALIGNMENT_IDENTITY_FLAG_MASK).to_le_bytes(),
        );

        let cigar = record.cigar();
        hash.update(&(cigar.len() as u64).to_le_bytes());
        for op in cigar.iter() {
            hash.update(&op.len().to_le_bytes());
            hash.update(&[op.char() as u8]);
        }
        hash.finish()
    }
}

/// Resolve one seed for the whole indexed sampling job. The resolved value
/// must be shared by every worker; resolving per worker would reintroduce
/// scheduling-dependent decisions. An omitted seed remains nondeterministic as
/// documented by the CLI.
pub(crate) fn resolve_master_seed(seed: Option<u64>) -> u64 {
    seed.unwrap_or_else(rand::random)
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use rust_htslib::bam::{
        self,
        record::{Aux, Cigar, CigarString},
    };

    use super::DeterministicFractionSampler;

    fn record(name: &str, pos: i64, flags: u16) -> bam::Record {
        let mut record = bam::Record::new();
        let cigar = CigarString(vec![Cigar::Match(8)]);
        record.set(name.as_bytes(), Some(&cigar), b"ACGTACGT", &[20; 8]);
        record.set_tid(3);
        record.set_pos(pos);
        record.set_flags(flags);
        record.push_aux(b"MM", Aux::String("C+m?,0;")).unwrap();
        record.push_aux(b"ML", Aux::ArrayU8((&[128][..]).into())).unwrap();
        record
    }

    #[test]
    fn boundary_fractions_are_exact() {
        let record = record("read", 10, 0);
        for seed in [0, 7, u64::MAX] {
            assert!(!DeterministicFractionSampler::new(seed, 0.0)
                .unwrap()
                .include(&record));
            assert!(DeterministicFractionSampler::new(seed, 1.0)
                .unwrap()
                .include(&record));
        }
    }

    #[test]
    fn rejects_non_finite_or_out_of_range_fractions() {
        for fraction in [f64::NEG_INFINITY, -0.1, 1.1, f64::INFINITY, f64::NAN]
        {
            assert!(DeterministicFractionSampler::new(7, fraction).is_err());
        }
    }

    #[test]
    fn identity_distinguishes_alignment_records_sharing_qname() {
        let sampler = DeterministicFractionSampler::new(7, 0.5).unwrap();
        let mut different_cigar = record("shared", 10, 0);
        different_cigar.set(
            b"shared",
            Some(&CigarString(vec![Cigar::Equal(8)])),
            b"ACGTACGT",
            &[20; 8],
        );
        let records = [
            record("shared", 10, 0),
            record("shared", 11, 0),
            different_cigar,
            record("shared", 10, 0x10),
            record("shared", 10, 0x40),
            record("shared", 10, 0x80),
            record("shared", 10, 0x100),
            record("shared", 10, 0x800),
        ];
        let scores = records
            .iter()
            .map(|record| sampler.alignment_score(record))
            .collect::<Vec<_>>();

        assert_eq!(scores.iter().copied().collect::<HashSet<_>>().len(), 8);
    }

    #[test]
    fn mutable_annotations_and_measured_content_are_not_identity() {
        let sampler = DeterministicFractionSampler::new(7, 0.5).unwrap();
        let original = record("read", 10, 0);
        let expected = sampler.alignment_score(&original);

        let mut changed = record("read", 10, 0);
        changed.set_mapq(42);
        changed.set_mtid(8);
        changed.set_mpos(900);
        changed.set_insert_size(1234);
        changed.set(
            b"read",
            Some(&CigarString(vec![Cigar::Match(8)])),
            b"TTTTTTTT",
            &[40; 8],
        );
        changed.remove_aux(b"MM").unwrap();
        changed.remove_aux(b"ML").unwrap();
        changed.push_aux(b"MM", Aux::String("T+m?,0;")).unwrap();
        changed.push_aux(b"ML", Aux::ArrayU8((&[250][..]).into())).unwrap();
        changed.set_flags(0x2 | 0x200 | 0x400);

        assert_eq!(sampler.alignment_score(&changed), expected);
    }

    #[test]
    fn fixed_v1_score_locks_hash_compatibility() {
        let sampler = DeterministicFractionSampler::new(7, 0.5).unwrap();
        assert_eq!(
            sampler.alignment_score(&record("read", 10, 0x10 | 0x800)),
            6_384_807_673_264_054_871
        );
    }

    #[test]
    fn master_seed_and_integer_threshold_are_part_of_v1() {
        let record = record("read", 10, 0x10 | 0x800);
        let sampler = DeterministicFractionSampler::new(7, 0.5).unwrap();
        assert_ne!(
            sampler.alignment_score(&record),
            DeterministicFractionSampler::new(8, 0.5)
                .unwrap()
                .alignment_score(&record)
        );

        let draw = sampler.alignment_score(&record) >> 11;
        let at_draw = draw as f64 / super::DRAW_SPACE;
        let just_above = (draw + 1) as f64 / super::DRAW_SPACE;
        assert!(!DeterministicFractionSampler::new(7, at_draw)
            .unwrap()
            .include(&record));
        assert!(DeterministicFractionSampler::new(7, just_above)
            .unwrap()
            .include(&record));
    }

    #[test]
    fn fixed_population_is_monotone_and_statistically_representative() {
        const POPULATION_SIZE: usize = 10_000;
        const INTERIOR_FRACTION: f64 = 0.37;
        let low = DeterministicFractionSampler::new(7, 0.1).unwrap();
        let interior =
            DeterministicFractionSampler::new(7, INTERIOR_FRACTION).unwrap();
        let high = DeterministicFractionSampler::new(7, 0.9).unwrap();
        let mut included = 0usize;

        for i in 0..POPULATION_SIZE {
            let record = record(&format!("population-{i}"), i as i64, 0);
            let at_low = low.include(&record);
            let at_interior = interior.include(&record);
            let at_high = high.include(&record);
            assert!(!at_low || at_interior);
            assert!(!at_interior || at_high);
            included += at_interior as usize;
        }

        let expected = POPULATION_SIZE as f64 * INTERIOR_FRACTION;
        let standard_deviation = (expected * (1.0 - INTERIOR_FRACTION)).sqrt();
        let conservative_tolerance = 8.0 * standard_deviation;
        assert!(
            (included as f64 - expected).abs() <= conservative_tolerance,
            "fixed population included {included} records; expected about \
             {expected} within {conservative_tolerance}"
        );
    }
}
