use std::collections::BTreeMap;

use bitvec::{bitvec, order::Lsb0};
use rust_htslib::bam::{
    self,
    record::{Aux, Cigar, CigarString},
};

use super::{
    AlignedBaseAndModArgmaxProbs, AlignedBaseArgmaxProbs,
    BaseAndModArgmaxProbs, BaseArgmaxProbs, ExtractsMleProbs, ProbsExtractor,
    QualHist,
};
use crate::{
    interval_chunks::{ChromCoordinates, FocusPositions2},
    mod_base_code::DnaBase,
    reads_sampler::deterministic_sampler::DeterministicFractionSampler,
};

fn make_record(name: &str, pos: i64) -> bam::Record {
    let mut record = bam::Record::new();
    let cigar = CigarString(vec![Cigar::Match(1)]);
    record.set(name.as_bytes(), Some(&cigar), b"C", &[255]);
    record.set_tid(0);
    record.set_pos(pos);
    record.push_aux(b"MM", Aux::String("C+m?,0;")).unwrap();
    record.push_aux(b"ML", Aux::ArrayU8((&[128][..]).into())).unwrap();
    record
}

fn make_long_record(name: &str) -> bam::Record {
    const READ_LENGTH: usize = 12;
    let mut record = bam::Record::new();
    let cigar = CigarString(vec![Cigar::Match(READ_LENGTH as u32)]);
    record.set(
        name.as_bytes(),
        Some(&cigar),
        &vec![b'C'; READ_LENGTH],
        &vec![255; READ_LENGTH],
    );
    record.set_tid(0);
    record.set_pos(0);
    let mm = format!("C+m?,{};", vec!["0"; READ_LENGTH].join(","));
    record.push_aux(b"MM", Aux::String(&mm)).unwrap();
    record
        .push_aux(b"ML", Aux::ArrayU8((&vec![128; READ_LENGTH]).into()))
        .unwrap();
    record
}

fn process_record<T>(
    extractor: &mut ProbsExtractor,
    record: &bam::Record,
    interval_start: u32,
) -> bool
where
    ProbsExtractor: ExtractsMleProbs<T>,
{
    let mut mask = bitvec![usize, Lsb0; 0; 2];
    mask.set(0, true);
    let coords = ChromCoordinates::new(
        0,
        interval_start,
        interval_start + 1,
        FocusPositions2::MaskedPositions { mask },
        true,
    );
    let mut hist = QualHist::default();
    <ProbsExtractor as ExtractsMleProbs<T>>::process_record(
        extractor,
        record,
        &coords,
        &mut hist.explicit_canonical_probs,
        &mut hist.hist,
        &mut hist.base_totals,
        &mut hist.mods_hists,
        &mut hist.num_records_with_base_mods,
    )
    .unwrap()
}

fn sampling_decisions<T>(
    records: &[bam::Record],
    assignments: &[usize],
    worker_count: usize,
) -> BTreeMap<Vec<u8>, bool>
where
    ProbsExtractor: ExtractsMleProbs<T>,
{
    let mut workers = (0..worker_count)
        .map(|_| {
            <ProbsExtractor as ExtractsMleProbs<T>>::new(
                7,
                true,
                0.5,
                [DnaBase::A; 4],
                None,
            )
            .unwrap()
        })
        .collect::<Vec<_>>();

    records
        .iter()
        .zip(assignments)
        .map(|(record, worker)| {
            (
                record.qname().to_vec(),
                process_record::<T>(
                    &mut workers[*worker],
                    record,
                    record.pos() as u32,
                ),
            )
        })
        .collect()
}

#[test]
fn seeded_fraction_sampling_is_independent_of_worker_assignment() {
    let records = (0..12)
        .map(|i| make_record(&format!("read-{i}"), i))
        .collect::<Vec<_>>();
    fn assert_stable<T>(records: &[bam::Record])
    where
        ProbsExtractor: ExtractsMleProbs<T>,
    {
        let serial = sampling_decisions::<T>(records, &[0; 12], 1);
        let scheduled = sampling_decisions::<T>(
            records,
            &[0, 1, 2, 0, 2, 1, 1, 0, 2, 2, 0, 1],
            3,
        );
        assert_eq!(serial, scheduled);
    }

    assert_stable::<BaseArgmaxProbs>(&records);
    assert_stable::<BaseAndModArgmaxProbs>(&records);
    assert_stable::<AlignedBaseArgmaxProbs>(&records);
    assert_stable::<AlignedBaseAndModArgmaxProbs>(&records);
}

#[test]
fn seeded_fraction_sampling_is_consistent_across_interval_refetches() {
    let sampler = DeterministicFractionSampler::new(7, 0.5).unwrap();
    for should_include in [false, true] {
        let record = (0..100)
            .map(|i| make_long_record(&format!("long-read-{i}")))
            .find(|record| sampler.include(record) == should_include)
            .unwrap();
        let mut extractor = <ProbsExtractor as ExtractsMleProbs<
            AlignedBaseArgmaxProbs,
        >>::new(7, true, 0.5, [DnaBase::A; 4], None)
        .unwrap();
        let mut hist = QualHist::default();
        let mut counted = Vec::new();

        for (start, final_interval) in [(0, false), (4, false), (8, true)] {
            let mut mask = bitvec![usize, Lsb0; 0; 8];
            for position in (0..8).step_by(2) {
                mask.set(position, true);
            }
            let coords = ChromCoordinates::new(
                0,
                start,
                start + 4,
                FocusPositions2::MaskedPositions { mask },
                final_interval,
            );
            counted.push(
                <ProbsExtractor as ExtractsMleProbs<
                    AlignedBaseArgmaxProbs,
                >>::process_record(
                    &mut extractor,
                    &record,
                    &coords,
                    &mut hist.explicit_canonical_probs,
                    &mut hist.hist,
                    &mut hist.base_totals,
                    &mut hist.mods_hists,
                    &mut hist.num_records_with_base_mods,
                )
                .unwrap(),
            );
        }

        if should_include {
            assert_eq!(counted, [false, false, true]);
            assert_eq!(hist.base_totals[DnaBase::C as usize], 12);
            assert_eq!(hist.num_records_with_base_mods[DnaBase::C as usize], 1);
        } else {
            assert_eq!(counted, [false, false, false]);
            assert_eq!(hist.base_totals, [0; 4]);
            assert_eq!(hist.num_records_with_base_mods, [0; 4]);
        }
    }
}
