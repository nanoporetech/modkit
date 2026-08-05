use std::collections::{BTreeMap, HashMap};

use rust_htslib::bam::record::{Aux, Cigar, CigarString};
use rust_htslib::bam::Record;

use super::{
    machine_parseable_table, process_bam_record, BaseStatus,
    ChromStrandPositionNames, StatusProbs, TidToChrom,
};
use crate::mod_bam::EdgeFilter;
use crate::mod_base_code::{DnaBase, ModCodeRepr};
use crate::util::Strand;

const M: BaseStatus = BaseStatus::Modified(ModCodeRepr::Code('m'));

fn make_record(
    name: &str,
    sequence: &str,
    cigar: Vec<Cigar>,
    mm: &str,
    ml: &[u8],
) -> Record {
    let mut record = Record::new();
    let cigar = CigarString(cigar);
    record.set(
        name.as_bytes(),
        Some(&cigar),
        sequence.as_bytes(),
        &vec![30; sequence.len()],
    );
    record.set_tid(0);
    record.set_pos(0);
    record.set_mapq(60);
    record.push_aux(b"MM", Aux::String(mm)).unwrap();
    record.push_aux(b"ML", Aux::ArrayU8(ml.into())).unwrap();
    record.push_aux(b"MN", Aux::U32(sequence.len() as u32)).unwrap();
    record
}

fn truth(
    strands: impl IntoIterator<Item = Strand>,
    status: BaseStatus,
) -> ChromStrandPositionNames {
    let strand_positions = strands
        .into_iter()
        .map(|strand| (strand, BTreeMap::from([(1, status)])))
        .collect();
    HashMap::from([("chr1".to_string(), strand_positions)])
}

fn truth_positions(
    strand: Strand,
    positions: impl IntoIterator<Item = i64>,
    status: BaseStatus,
) -> ChromStrandPositionNames {
    HashMap::from([(
        "chr1".to_string(),
        HashMap::from([(
            strand,
            positions.into_iter().map(|pos| (pos, status)).collect(),
        )]),
    )])
}

fn classify(
    records: impl IntoIterator<Item = Record>,
    truth: &ChromStrandPositionNames,
) -> StatusProbs {
    classify_with_edge_filter(records, truth, None)
}

fn classify_with_edge_filter(
    records: impl IntoIterator<Item = Record>,
    truth: &ChromStrandPositionNames,
    edge_filter: Option<&EdgeFilter>,
) -> StatusProbs {
    let tid_to_chrom: TidToChrom = HashMap::from([(0, "chr1".to_string())]);
    let mut combined = StatusProbs::new();
    for record in records {
        let observed = process_bam_record(
            &record,
            truth,
            &tid_to_chrom,
            DnaBase::C,
            None,
            edge_filter,
        )
        .unwrap();
        for (status, probs) in observed {
            combined.entry(status).or_default().extend(probs);
        }
    }
    combined
}

fn called_record(name: &str, mod_strand: char) -> Record {
    let (sequence, fundamental_base) = match mod_strand {
        '+' => ("CC", 'C'),
        '-' => ("GG", 'G'),
        _ => panic!("invalid modification strand"),
    };
    make_record(
        name,
        sequence,
        vec![Cigar::Match(2)],
        &format!("{fundamental_base}{mod_strand}m?,1;"),
        &[255],
    )
}

fn uncalled_record(name: &str, sequence: &str, mod_strand: char) -> Record {
    let fundamental_base = match mod_strand {
        '+' => 'C',
        '-' => 'G',
        _ => panic!("invalid modification strand"),
    };
    make_record(
        name,
        sequence,
        vec![Cigar::Match(2)],
        &format!("{fundamental_base}{mod_strand}m?,0;"),
        &[255],
    )
}

fn deletion_record(name: &str) -> Record {
    make_record(
        name,
        "C",
        vec![Cigar::Match(1), Cigar::Del(1)],
        "C+m?,0;",
        &[255],
    )
}

fn empty_descriptor_record(
    name: &str,
    sequence: &str,
    mod_strand: char,
) -> Record {
    let fundamental_base = match mod_strand {
        '+' => 'C',
        '-' => 'G',
        _ => panic!("invalid modification strand"),
    };
    make_record(
        name,
        sequence,
        vec![Cigar::Match(sequence.len() as u32)],
        &format!("{fundamental_base}{mod_strand}m?;"),
        &[],
    )
}

fn assert_exact_counts(
    observed: &StatusProbs,
    expected: &[((BaseStatus, BaseStatus), usize)],
) {
    let observed = observed
        .iter()
        .map(|(statuses, probs)| (*statuses, probs.len()))
        .collect::<BTreeMap<_, _>>();
    let expected = expected.iter().copied().collect::<BTreeMap<_, _>>();
    assert_eq!(observed, expected);
}

#[test]
fn explicit_mode_counts_no_calls_without_truth_overlapping_seed_call() {
    let records =
        (0..3).map(|i| called_record(&format!("called-{i}"), '+')).chain(
            (0..6)
                .map(|i| uncalled_record(&format!("uncalled-{i}"), "CC", '+')),
        );
    let observed = classify(records, &truth([Strand::Positive], M));

    assert_exact_counts(
        &observed,
        &[((M, M), 3), ((M, BaseStatus::NoCall), 6)],
    );
    assert_eq!(
        machine_parseable_table(DnaBase::C, &observed),
        "[[\"ground_truth_label\",\"m\",\"No Call\"],[\"m\",3,6]]"
    );
}

#[test]
fn explicit_mode_counts_unseeded_no_calls_mismatches_and_deletions() {
    let records = (0..2)
        .map(|i| called_record(&format!("called-{i}"), '+'))
        .chain(
            (0..6)
                .map(|i| uncalled_record(&format!("uncalled-{i}"), "CC", '+')),
        )
        .chain(
            (0..2)
                .map(|i| uncalled_record(&format!("mismatch-{i}"), "CA", '+')),
        )
        .chain((0..2).map(|i| deletion_record(&format!("deletion-{i}"))));
    let observed = classify(records, &truth([Strand::Positive], M));

    assert_exact_counts(
        &observed,
        &[
            ((M, M), 2),
            ((M, BaseStatus::NoCall), 6),
            ((M, BaseStatus::Mismatch(DnaBase::A)), 2),
            ((M, BaseStatus::Deletion), 2),
        ],
    );
    assert_eq!(
        machine_parseable_table(DnaBase::C, &observed),
        "[[\"ground_truth_label\",\"m\",\"No Call\",\"A\",\"Deletion\"],[\"m\",2,6,2,2]]"
    );
}

#[test]
fn explicit_mode_counts_only_the_descriptor_reference_strand() {
    let records = std::iter::once(called_record("positive-called", '+'))
        .chain((0..5).map(|i| {
            uncalled_record(&format!("positive-uncalled-{i}"), "CC", '+')
        }))
        .chain(std::iter::once(called_record("negative-called", '-')))
        .chain((0..5).map(|i| {
            uncalled_record(&format!("negative-uncalled-{i}"), "GG", '-')
        }));
    let observed =
        classify(records, &truth([Strand::Positive, Strand::Negative], M));

    assert_exact_counts(
        &observed,
        &[((M, M), 2), ((M, BaseStatus::NoCall), 10)],
    );
    assert_eq!(
        machine_parseable_table(DnaBase::C, &observed),
        "[[\"ground_truth_label\",\"m\",\"No Call\"],[\"m\",2,10]]"
    );
}

#[test]
fn implicit_mode_canonical_calls_are_not_reclassified_as_no_calls() {
    let records = (0..3).map(|i| {
        make_record(
            &format!("implicit-{i}"),
            "CC",
            vec![Cigar::Match(2)],
            "C+m.,0;",
            &[255],
        )
    });
    let observed =
        classify(records, &truth([Strand::Positive], BaseStatus::Canonical));

    assert_exact_counts(
        &observed,
        &[((BaseStatus::Canonical, BaseStatus::Canonical), 3)],
    );
    assert_eq!(
        machine_parseable_table(DnaBase::C, &observed),
        "[[\"ground_truth_label\",\"C\"],[\"C\",3]]"
    );
}

#[test]
fn explicit_empty_positive_descriptor_counts_no_call() {
    let observed = classify(
        [empty_descriptor_record("empty-positive", "CC", '+')],
        &truth([Strand::Positive], M),
    );

    assert_exact_counts(&observed, &[((M, BaseStatus::NoCall), 1)]);
}

#[test]
fn explicit_empty_negative_descriptor_counts_no_call() {
    let observed = classify(
        [empty_descriptor_record("empty-negative", "GG", '-')],
        &truth([Strand::Negative], M),
    );

    assert_exact_counts(&observed, &[((M, BaseStatus::NoCall), 1)]);
}

#[test]
fn explicit_empty_n_descriptor_expands_to_target_base() {
    let record =
        make_record("empty-n", "CC", vec![Cigar::Match(2)], "N+m?;", &[]);
    let observed = classify([record], &truth([Strand::Positive], M));

    assert_exact_counts(&observed, &[((M, BaseStatus::NoCall), 1)]);
}

#[test]
fn unrelated_canonical_base_descriptor_does_not_create_observations() {
    let unrelated = make_record(
        "adenine-only",
        "AA",
        vec![Cigar::Match(2)],
        "A+a?,1;",
        &[255],
    );
    let observed =
        classify([unrelated], &truth([Strand::Positive, Strand::Negative], M));

    assert!(observed.is_empty());
}

#[test]
fn minus_strand_call_on_reverse_alignment_uses_modified_primary_base() {
    let mut record = make_record(
        "reverse-g-minus",
        "CC",
        vec![Cigar::Match(2)],
        "G-m?,0;",
        &[255],
    );
    record.set_reverse();

    let observed = classify([record], &truth([Strand::Positive], M));

    assert_exact_counts(&observed, &[((M, M), 1)]);
}

#[test]
fn unrelated_empty_descriptor_does_not_create_observations() {
    let unrelated = make_record(
        "empty-adenine-only",
        "AA",
        vec![Cigar::Match(2)],
        "A+a?;",
        &[],
    );
    let observed =
        classify([unrelated], &truth([Strand::Positive, Strand::Negative], M));

    assert!(observed.is_empty());
}

#[test]
fn explicit_unseeded_sites_compose_reference_skips_and_iupac_no_calls() {
    let record = make_record(
        "unseeded-refskip-iupac",
        "CNRAC",
        vec![Cigar::Match(1), Cigar::RefSkip(3), Cigar::Match(4)],
        "C+m?,0;",
        &[255],
    );
    let observed =
        classify([record], &truth_positions(Strand::Positive, 1..8, M));

    assert_exact_counts(
        &observed,
        &[
            ((M, BaseStatus::NoCall), 3),
            ((M, BaseStatus::Mismatch(DnaBase::A)), 1),
        ],
    );
}

#[test]
fn reverse_minus_fallback_preserves_no_call_mismatch_and_deletion_orientation()
{
    let mut record = make_record(
        "reverse-minus-fallback",
        "CAC",
        vec![Cigar::Match(2), Cigar::Del(1), Cigar::Match(1)],
        "G-m?;",
        &[],
    );
    record.set_reverse();
    let observed =
        classify([record], &truth_positions(Strand::Positive, 0..3, M));

    assert_exact_counts(
        &observed,
        &[
            ((M, BaseStatus::NoCall), 1),
            ((M, BaseStatus::Mismatch(DnaBase::A)), 1),
            ((M, BaseStatus::Deletion), 1),
        ],
    );
}

fn positive_truth_range(end: i64) -> ChromStrandPositionNames {
    HashMap::from([(
        "chr1".to_string(),
        HashMap::from([(
            Strand::Positive,
            (0..end).map(|position| (position, M)).collect(),
        )]),
    )])
}

fn forward_edge_filter_record() -> Record {
    make_record(
        "forward-edge-filter",
        "CCCACCC",
        vec![Cigar::Match(7)],
        "C+m?,0,0,1,1;",
        &[255; 4],
    )
}

#[test]
fn ordinary_edge_filter_excludes_filtered_fallback_sites() {
    let edge_filter = EdgeFilter::new(1, 1, false);
    let observed = classify_with_edge_filter(
        [forward_edge_filter_record()],
        &positive_truth_range(7),
        Some(&edge_filter),
    );

    assert_exact_counts(
        &observed,
        &[
            ((M, M), 2),
            ((M, BaseStatus::NoCall), 2),
            ((M, BaseStatus::Mismatch(DnaBase::A)), 1),
        ],
    );
}

#[test]
fn inverted_edge_filter_excludes_filtered_fallback_sites() {
    let edge_filter = EdgeFilter::new(1, 1, true);
    let observed = classify_with_edge_filter(
        [forward_edge_filter_record()],
        &positive_truth_range(7),
        Some(&edge_filter),
    );

    assert_exact_counts(&observed, &[((M, M), 2)]);
}

#[test]
fn reverse_alignment_edge_filter_uses_forward_query_coordinates() {
    let mut record = make_record(
        "reverse-edge-filter",
        "CACNCC",
        vec![Cigar::Match(6)],
        "G-m?;",
        &[],
    );
    record.set_reverse();
    let truth = positive_truth_range(6);

    let ordinary = EdgeFilter::new(1, 2, false);
    let observed =
        classify_with_edge_filter([record.clone()], &truth, Some(&ordinary));
    assert_exact_counts(&observed, &[((M, BaseStatus::NoCall), 3)]);

    let inverted = EdgeFilter::new(1, 2, true);
    let observed = classify_with_edge_filter([record], &truth, Some(&inverted));
    assert_exact_counts(
        &observed,
        &[
            ((M, BaseStatus::NoCall), 2),
            ((M, BaseStatus::Mismatch(DnaBase::A)), 1),
        ],
    );
}

#[test]
fn edge_filter_does_not_reclassify_deletion_without_query_position() {
    let record = make_record(
        "edge-filter-deletion",
        "CCC",
        vec![Cigar::Match(1), Cigar::Del(1), Cigar::Match(2)],
        "C+m?,0;",
        &[255],
    );
    let truth = HashMap::from([(
        "chr1".to_string(),
        HashMap::from([(Strand::Positive, BTreeMap::from([(1, M)]))]),
    )]);
    let edge_filter = EdgeFilter::new(1, 1, false);

    let observed =
        classify_with_edge_filter([record], &truth, Some(&edge_filter));

    assert_exact_counts(&observed, &[((M, BaseStatus::Deletion), 1)]);
}
