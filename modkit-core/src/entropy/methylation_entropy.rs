use derive_new::new;
use itertools::Itertools;
use log::debug;
use rustc_hash::FxHashSet;
use std::collections::{BTreeSet, HashMap};

#[derive(Debug, Copy, Clone, Eq, PartialEq, Hash, Ord, PartialOrd)]
#[repr(transparent)]
pub(super) struct EntropySymbol(u32);

impl EntropySymbol {
    pub(super) const CANONICAL: Self = Self(0);
    pub(super) const FILTERED: Self = Self(u32::MAX);

    pub(super) fn called(id: u32) -> Self {
        assert!(
            id > Self::CANONICAL.0 && id < Self::FILTERED.0,
            "modified state ID must not collide with canonical or wildcard"
        );
        Self(id)
    }

    fn is_filtered(self) -> bool {
        self == Self::FILTERED
    }
}

pub(super) type EntropyPattern = Box<[EntropySymbol]>;

#[derive(new)]
struct AlphabetInfo {
    columns: Box<[Box<[EntropySymbol]>]>,
}

impl AlphabetInfo {
    fn from_sequences(
        sequences: &[EntropyPattern],
        window_size: usize,
    ) -> Self {
        assert_eq!(
            sequences.iter().map(|x| x.len()).unique().count(),
            1,
            "all sequences should be the same length {sequences:?}"
        );
        let mut columns = vec![BTreeSet::new(); window_size];
        for sequence in sequences {
            for (position, symbol) in sequence.iter().enumerate() {
                if !symbol.is_filtered() {
                    columns[position].insert(*symbol);
                }
            }
        }
        let columns = columns
            .into_iter()
            .map(|elements| {
                debug_assert!(
                    !elements.is_empty(),
                    "column with zero coverage in {sequences:?}, \
                     {elements:?}"
                );
                elements.into_iter().collect::<Vec<_>>().into_boxed_slice()
            })
            .collect::<Vec<_>>()
            .into_boxed_slice();

        Self { columns }
    }

    fn get_column(
        &self,
        idx: usize,
    ) -> impl Iterator<Item = EntropySymbol> + '_ {
        self.columns[idx].iter().copied()
    }
}

fn sequence_matches_pattern(
    sequence: &[EntropySymbol],
    pattern: &[EntropySymbol],
) -> bool {
    debug_assert_eq!(sequence.len(), pattern.len());
    debug_assert!(pattern.iter().all(|symbol| !symbol.is_filtered()));
    sequence
        .iter()
        .zip(pattern)
        .all(|(observed, called)| observed.is_filtered() || observed == called)
}

fn extend_pattern(
    pattern: &[EntropySymbol],
    symbol: EntropySymbol,
) -> EntropyPattern {
    let mut extended = Vec::with_capacity(pattern.len() + 1);
    extended.extend_from_slice(pattern);
    extended.push(symbol);
    extended.into_boxed_slice()
}

fn all_patterns_dp(
    sequences: &[EntropyPattern],
    window_size: usize,
    alphabet_info: &AlphabetInfo,
) -> Vec<EntropyPattern> {
    let sequences = sequences.iter().collect::<BTreeSet<&EntropyPattern>>();
    debug_assert!(
        sequences.iter().all(|x| x.len() == window_size),
        "all sequences should be the same length, {sequences:?}"
    );
    // easy case
    if !sequences.iter().any(|x| x.contains(&EntropySymbol::FILTERED)) {
        return sequences.into_iter().cloned().collect();
    }

    let basecase = alphabet_info
        .get_column(0)
        .map(|symbol| vec![symbol].into_boxed_slice())
        .collect::<FxHashSet<EntropyPattern>>();
    debug_assert!(basecase.len() >= 1, "first column has zero valid coverage");
    #[cfg(debug_assertions)]
    {
        let basecase_check = sequences
            .iter()
            .filter_map(|sequence| {
                let symbol = sequence[0];
                (!symbol.is_filtered()).then(|| vec![symbol].into_boxed_slice())
            })
            .collect::<FxHashSet<EntropyPattern>>();
        assert_eq!(basecase, basecase_check);
    }

    let all_combs = (1..window_size).fold(basecase, |acc, idx| {
        let mut acc_patterns = FxHashSet::default();
        for sequence in &sequences {
            let prefix = &sequence[..idx];
            for pattern in acc.iter() {
                if sequence_matches_pattern(prefix, pattern) {
                    let observed = sequence[idx];
                    if observed.is_filtered() {
                        for symbol in alphabet_info.get_column(idx) {
                            acc_patterns
                                .insert(extend_pattern(pattern, symbol));
                        }
                    } else {
                        acc_patterns.insert(extend_pattern(pattern, observed));
                    }
                }
            }
        }

        acc_patterns
    });

    all_combs.into_iter().sorted().collect::<Vec<EntropyPattern>>()
}

fn calc_entropy(sequences: &[EntropyPattern], window_size: usize) -> f32 {
    let alphabet_info = AlphabetInfo::from_sequences(sequences, window_size);
    let patterns = all_patterns_dp(sequences, window_size, &alphabet_info);

    let counts = sequences.iter().fold(HashMap::new(), |mut acc, seq| {
        let matches = patterns
            .iter()
            .filter(|pattern| sequence_matches_pattern(seq, pattern))
            .collect::<Vec<&EntropyPattern>>();
        assert!(matches.len() > 0, "no matches for {seq:?} in {patterns:?}");
        let factor = 1f32 / matches.len() as f32;
        for pattern in matches {
            *acc.entry(pattern).or_insert(0f32) += factor;
        }
        acc
    });

    let total = counts.values().sum::<f32>();
    if total - sequences.len() as f32 > 0.1f32 {
        let n_seqs = sequences.len() as f32;
        if total > n_seqs {
            debug!(
                "encountered discordant total value calculation, too high \
                 total={total} n_seqs={n_seqs}"
            );
        } else {
            debug!(
                "encountered discordant total value calculation, too low \
                 total={total} n_seqs={n_seqs}"
            );
        }
    }
    debug_assert!((total - sequences.len() as f32) < 1f32);
    counts
        .values()
        .map(|&x| {
            let p = x / total;
            p * (p.log2())
        })
        .sum::<f32>()
        * -1f32
}

pub(super) fn calc_me_entropy(
    sequences: &[EntropyPattern],
    window_size: usize,
    constant: f32,
) -> f32 {
    let shannons = calc_entropy(sequences, window_size);
    let me_entropy = constant * shannons;
    if me_entropy == -0f32 {
        0f32
    } else {
        me_entropy
    }
}

#[cfg(test)]
mod methylation_entropy_tests {
    use super::{
        all_patterns_dp, calc_entropy, calc_me_entropy, AlphabetInfo,
        EntropyPattern, EntropySymbol,
    };
    use assert_approx_eq::assert_approx_eq;
    use std::mem::{size_of, size_of_val};

    // Keep the original 0/1/* fixtures readable while ensuring that the
    // production entropy path has no text representation.
    fn legacy_pattern(raw: &str) -> EntropyPattern {
        raw.chars()
            .map(|symbol| match symbol {
                '0' => EntropySymbol::CANONICAL,
                '1' => EntropySymbol::called(1),
                '*' => EntropySymbol::FILTERED,
                _ => panic!("unsupported legacy test symbol {symbol}"),
            })
            .collect::<Vec<_>>()
            .into_boxed_slice()
    }

    fn legacy_patterns(raw: &[&str]) -> Vec<EntropyPattern> {
        raw.iter().map(|pattern| legacy_pattern(pattern)).collect()
    }

    fn pattern(states: &[Option<usize>]) -> EntropyPattern {
        states
            .iter()
            .map(|state| match state {
                Some(0) => EntropySymbol::CANONICAL,
                Some(id) => EntropySymbol::called((*id).try_into().unwrap()),
                None => EntropySymbol::FILTERED,
            })
            .collect::<Vec<_>>()
            .into_boxed_slice()
    }

    fn called(ids: &[usize]) -> EntropyPattern {
        ids.iter()
            .map(|id| match id {
                0 => EntropySymbol::CANONICAL,
                _ => EntropySymbol::called((*id).try_into().unwrap()),
            })
            .collect::<Vec<_>>()
            .into_boxed_slice()
    }

    fn p_q_w_oracle(code_count: usize) -> Vec<EntropyPattern> {
        let p = called(&(1..=code_count).collect::<Vec<_>>());
        let mut q = p.to_vec();
        q[0] = EntropySymbol::CANONICAL;
        let mut wildcard = p.to_vec();
        wildcard[0] = EntropySymbol::FILTERED;
        vec![p, q.into_boxed_slice(), wildcard.into_boxed_slice()]
    }

    fn assert_per_read_mass_and_counts(
        sequences: &[EntropyPattern],
        concrete_patterns: &[EntropyPattern],
    ) {
        let mut counts = vec![0.0f32; concrete_patterns.len()];
        for sequence in sequences {
            let matching_indices = concrete_patterns
                .iter()
                .enumerate()
                .filter_map(|(idx, concrete)| {
                    super::sequence_matches_pattern(sequence, concrete)
                        .then_some(idx)
                })
                .collect::<Vec<_>>();
            assert!(!matching_indices.is_empty());
            let factor = 1.0 / matching_indices.len() as f32;
            let read_mass = factor * matching_indices.len() as f32;
            assert_eq!(read_mass, 1.0);
            for idx in matching_indices {
                counts[idx] += factor;
            }
        }
        assert_eq!(counts, vec![1.5, 1.5]);
        assert_eq!(counts.iter().sum::<f32>(), sequences.len() as f32);
    }

    #[test]
    fn entropy_symbols_are_compact_and_patterns_have_no_capacity_field() {
        assert_eq!(size_of::<EntropySymbol>(), size_of::<u32>());
        assert_eq!(
            size_of::<EntropyPattern>(),
            size_of::<Box<[EntropySymbol]>>()
        );
        let sequence = called(&[0, 1, 10, 100]);
        assert_eq!(
            size_of_val(sequence.as_ref()),
            sequence.len() * size_of::<u32>()
        );
        assert!(std::panic::catch_unwind(|| EntropySymbol::called(0)).is_err());
        assert!(std::panic::catch_unwind(|| EntropySymbol::called(u32::MAX))
            .is_err());
    }

    #[test]
    fn test_calc_entropy() {
        let sequences = legacy_patterns(&["0000", "0000", "0000", "0000"]);
        assert_eq!(calc_me_entropy(&sequences, 4, 0.25), 0.0);
        let sequences = legacy_patterns(&["1111", "1111", "1111", "1111"]);
        assert_eq!(calc_me_entropy(&sequences, 4, 0.25), 0.0);
        let sequences = legacy_patterns(&["0010", "0010", "0010", "0010"]);
        assert_eq!(calc_me_entropy(&sequences, 4, 0.25), 0.0);
        let sequences = legacy_patterns(&[
            "1111", "1111", "1111", "1111", "0000", "0000", "0000", "0000",
        ]);
        assert_eq!(calc_me_entropy(&sequences, 4, 0.25), 0.25);
        let sequences = legacy_patterns(&[
            "1111", "1111", "0011", "0011", "1100", "1100", "0000", "0000",
        ]);
        assert_eq!(calc_me_entropy(&sequences, 4, 0.25), 0.50);
        let sequences = legacy_patterns(&[
            "0000", "1111", "0101", "0111", "0111", "0111", "0000", "1111",
        ]);
        assert_eq!(calc_me_entropy(&sequences, 4, 0.25), 0.47640976);
    }

    #[test]
    fn test_calc_entropy_wildcards() {
        let sequences = legacy_patterns(&["1*01", "1111", "1011", "1111"]);

        let alphabet_info = AlphabetInfo::from_sequences(&sequences, 4);
        let patterns = all_patterns_dp(&sequences, 4, &alphabet_info);
        assert_eq!(
            patterns,
            legacy_patterns(&["1001", "1011", "1101", "1111"])
        );
        let entropy = calc_entropy(&sequences, 4);
        assert_eq!(entropy, 1.75);

        let sequences = legacy_patterns(&["1*11", "1111", "1011", "1111"]);

        let alphabet_info = AlphabetInfo::from_sequences(&sequences, 4);
        let patterns = all_patterns_dp(&sequences, 4, &alphabet_info);
        assert_eq!(patterns, legacy_patterns(&["1011", "1111"]));
        let entropy = calc_entropy(&sequences, 4);
        assert_eq!(entropy, 0.95443404);

        let sequences = legacy_patterns(&["1*01", "1101", "1011", "1111"]);

        let alphabet_info = AlphabetInfo::from_sequences(&sequences, 4);
        let patterns = all_patterns_dp(&sequences, 4, &alphabet_info);
        assert_eq!(
            patterns,
            legacy_patterns(&["1001", "1011", "1101", "1111"])
        );
        let entropy = calc_entropy(&sequences, 4);
        assert_approx_eq!(entropy, 1.9, 0.01);

        let sequences = legacy_patterns(&["*010", "1010", "0010"]);

        let alphabet_info = AlphabetInfo::from_sequences(&sequences, 4);
        let patterns = all_patterns_dp(&sequences, 4, &alphabet_info);
        assert_eq!(patterns, legacy_patterns(&["0010", "1010"]));
        let entropy = calc_entropy(&sequences, 4);
        assert_eq!(entropy, 1.0f32);

        let sequences = legacy_patterns(&["1010", "1010", "1010", "1010"]);

        let _alphabet_info = AlphabetInfo::from_sequences(&sequences, 4);
        let entropy = calc_entropy(&sequences, 4);
        assert_eq!(entropy, 0f32);
    }

    #[test]
    fn repeated_one_code_control_has_zero_entropy() {
        let sequences = vec![called(&[12]); 8];
        assert_eq!(calc_entropy(&sequences, 1), 0.0);
    }

    #[test]
    fn p_q_w_oracles_cover_large_state_sets_and_read_order_permutations() {
        for code_count in [9, 12, 17] {
            let sequences = p_q_w_oracle(code_count);
            let alphabet = AlphabetInfo::from_sequences(&sequences, code_count);
            let concrete = all_patterns_dp(&sequences, code_count, &alphabet);
            assert_eq!(concrete.len(), 2);
            assert_per_read_mass_and_counts(&sequences, &concrete);

            let expected = 1.0 / code_count as f32;
            assert_eq!(
                calc_me_entropy(
                    &sequences,
                    code_count,
                    1.0 / code_count as f32,
                ),
                expected
            );

            let mut reversed = sequences.clone();
            reversed.reverse();
            assert_eq!(
                calc_me_entropy(&reversed, code_count, 1.0 / code_count as f32,),
                expected
            );
        }
    }

    #[test]
    fn twelve_code_wildcard_matches_uniform_hand_calculation() {
        let mut sequences =
            (1..=12).map(|id| called(&[id])).collect::<Vec<EntropyPattern>>();
        sequences.push(pattern(&[None]));

        let entropy = calc_entropy(&sequences, 1);
        assert_approx_eq!(entropy, 12f32.log2(), 0.000_01);
    }

    #[test]
    #[should_panic]
    fn test_alphabet_info() {
        // test that assert fires when columns are all *
        let sequences = legacy_patterns(&["*111", "*111", "*111", "*111"]);
        AlphabetInfo::from_sequences(&sequences, 4);
    }
}
