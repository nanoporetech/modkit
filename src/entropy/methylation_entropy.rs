use derive_new::new;
use itertools::Itertools;
use log_once::debug_once;
use regex::Regex;
use rustc_hash::{FxHashMap, FxHashSet};
use std::collections::{BTreeSet, HashMap};
use std::str::Chars;
use substring::Substring;

#[derive(new)]
struct AlphabetInfo {
    columns: FxHashMap<usize, String>,
}

impl AlphabetInfo {
    fn from_sequences(sequences: &[String], window_size: usize) -> Self {
        assert_eq!(
            sequences.iter().map(|x| x.len()).unique().count(),
            1,
            "all sequences should be the same length {sequences:?}"
        );
        let columns =
            (0..window_size).fold(FxHashMap::default(), |mut acc, idx| {
                acc.insert(idx, BTreeSet::new());
                acc
            });
        let columns = sequences
            .iter()
            .fold(columns, |mut acc, seq| {
                for (i, c) in seq.chars().enumerate().filter(|(_, c)| *c != '*')
                {
                    acc.entry(i).or_insert_with(BTreeSet::new).insert(c);
                }
                acc
            })
            .into_iter()
            .inspect(|(_, elems)| {
                debug_assert!(
                    elems.len() >= 1,
                    "column with zero coverage in {sequences:?}, {elems:?}"
                )
            })
            .map(|(pos, elements)| {
                (pos, elements.into_iter().collect::<String>())
            })
            .collect::<FxHashMap<usize, String>>();

        Self { columns }
    }

    fn seq_to_regex(&self, seq: &str) -> Regex {
        let pattern = seq
            .chars()
            .enumerate()
            .map(|(pos, c)| match c {
                '*' => {
                    format!("[{}]", self.columns.get(&pos).unwrap())
                }
                _ => {
                    format!("{c}")
                }
            })
            .collect::<String>();
        Regex::new(&pattern).unwrap()
    }

    fn get_column(&self, idx: usize) -> Chars {
        self.columns.get(&idx).unwrap().chars()
    }
}

fn all_patterns_dp(
    sequences: &[String],
    window_size: usize,
    alphabet_info: &AlphabetInfo,
) -> Vec<String> {
    let sequences = sequences.iter().collect::<BTreeSet<&String>>();
    debug_assert!(
        sequences.iter().all(|x| x.len() == window_size),
        "all sequences should be the same length, {sequences:?}"
    );
    // easy case
    if !sequences.iter().any(|x| x.contains('*')) {
        return sequences.iter().map(|x| x.to_string()).collect();
    }

    let basecase = alphabet_info
        .get_column(0)
        .map(|c| format!("{c}"))
        .collect::<FxHashSet<String>>();
    debug_assert!(basecase.len() >= 1, "first column has zero valid coverage");
    #[cfg(debug_assertions)]
    {
        let basecase_check = sequences
            .iter()
            .map(|x| x.substring(0, 1).to_string())
            .filter(|x| x != "*")
            .collect::<FxHashSet<String>>();
        assert_eq!(basecase, basecase_check);
    }

    let mut cache = FxHashMap::default();
    let all_combs = (1..window_size).fold(basecase, |acc, idx| {
        let mut acc_patterns = FxHashSet::default();
        for seq in sequences.iter() {
            let prefix = seq.substring(0, idx);
            let re = if let Some(re) = cache.get(prefix) {
                re
            } else {
                let re = alphabet_info.seq_to_regex(prefix);
                cache.insert(prefix, re);
                cache.get(prefix).unwrap()
            };
            for pattern in acc.iter() {
                if re.is_match(&pattern) {
                    let last_letter = seq
                        .chars()
                        .nth(idx)
                        .expect(&format!("should get last letter at {idx}"));
                    match last_letter {
                        '*' => {
                            for x in alphabet_info.get_column(idx) {
                                let new_pattern = format!("{pattern}{x}");
                                acc_patterns.insert(new_pattern);
                            }
                        }
                        _ => {
                            let new_pattern = format!("{pattern}{last_letter}");
                            acc_patterns.insert(new_pattern);
                        }
                    }
                }
            }
        }

        acc_patterns
    });

    all_combs.into_iter().sorted().collect::<Vec<String>>()
}

fn calc_entropy(sequences: &[String], window_size: usize) -> f32 {
    let mut alphabet_info =
        AlphabetInfo::from_sequences(sequences, window_size);
    let patterns = all_patterns_dp(sequences, window_size, &mut alphabet_info);

    let mut cache = FxHashMap::default();
    let counts = sequences.iter().fold(HashMap::new(), |mut acc, seq| {
        // let re = seq_to_regex(seq, &alphabet_info.wildcard_regex);
        let re = if let Some(re) = cache.get(seq) {
            re
        } else {
            let re = alphabet_info.seq_to_regex(seq);
            cache.insert(seq, re);
            cache.get(seq).unwrap()
        };
        let matches = patterns
            .iter()
            .filter(|p| re.is_match(p))
            .collect::<Vec<&String>>();
        assert!(matches.len() > 0, "no matches for {seq} in {patterns:?}");
        let factor = 1f32 / matches.len() as f32;
        for pattern in matches {
            *acc.entry(pattern).or_insert(0f32) += factor;
        }
        acc
    });

    let total = counts.values().sum::<f32>();
    if total - sequences.len() as f32 > 1e-3 {
        if total > sequences.len() as f32 {
            debug_once!(
                "encountered discordant total value calculation, too high"
            );
        } else {
            debug_once!(
                "encountered discordant total value calculation, too low"
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
    sequences: &[String],
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
    use crate::entropy::methylation_entropy::{
        all_patterns_dp, calc_entropy, calc_me_entropy, AlphabetInfo,
    };
    use assert_approx_eq::assert_approx_eq;

    #[test]
    fn test_calc_entropy() {
        let sequences = vec![
            "0000".to_string(),
            "0000".to_string(),
            "0000".to_string(),
            "0000".to_string(),
        ];
        assert_eq!(calc_me_entropy(&sequences, 4, 0.25), 0.0);
        let sequences = vec![
            "1111".to_string(),
            "1111".to_string(),
            "1111".to_string(),
            "1111".to_string(),
        ];
        assert_eq!(calc_me_entropy(&sequences, 4, 0.25), 0.0);
        let sequences = vec![
            "0010".to_string(),
            "0010".to_string(),
            "0010".to_string(),
            "0010".to_string(),
        ];
        assert_eq!(calc_me_entropy(&sequences, 4, 0.25), 0.0);
        let sequences = vec![
            "1111".to_string(),
            "1111".to_string(),
            "1111".to_string(),
            "1111".to_string(),
            "0000".to_string(),
            "0000".to_string(),
            "0000".to_string(),
            "0000".to_string(),
        ];
        assert_eq!(calc_me_entropy(&sequences, 4, 0.25), 0.25);
        let sequences = vec![
            "1111".to_string(),
            "1111".to_string(),
            "0011".to_string(),
            "0011".to_string(),
            "1100".to_string(),
            "1100".to_string(),
            "0000".to_string(),
            "0000".to_string(),
        ];
        assert_eq!(calc_me_entropy(&sequences, 4, 0.25), 0.50);
        let sequences = vec![
            "0000".to_string(),
            "1111".to_string(),
            "0101".to_string(),
            "0111".to_string(),
            "0111".to_string(),
            "0111".to_string(),
            "0000".to_string(),
            "1111".to_string(),
        ];
        assert_eq!(calc_me_entropy(&sequences, 4, 0.25), 0.47640976);
    }

    #[test]
    fn test_calc_entropy_wildcards() {
        let sequences = vec!["1*01", "1111", "1011", "1111"]
            .into_iter()
            .map(|x| x.to_string())
            .collect::<Vec<String>>();

        let alphabet_info = AlphabetInfo::from_sequences(&sequences, 4);
        let patterns = all_patterns_dp(&sequences, 4, &alphabet_info);
        assert_eq!(
            patterns,
            vec![
                "1001".to_string(),
                "1011".to_string(),
                "1101".to_string(),
                "1111".to_string(),
            ]
        );
        let entropy = calc_entropy(&sequences, 4);
        assert_eq!(entropy, 1.75);

        let sequences = vec!["1*11", "1111", "1011", "1111"]
            .into_iter()
            .map(|x| x.to_string())
            .collect::<Vec<String>>();

        let alphabet_info = AlphabetInfo::from_sequences(&sequences, 4);
        let patterns = all_patterns_dp(&sequences, 4, &alphabet_info);
        assert_eq!(patterns, vec!["1011".to_string(), "1111".to_string(),]);
        let entropy = calc_entropy(&sequences, 4);
        assert_eq!(entropy, 0.95443404);

        let sequences = vec!["1*01", "1101", "1011", "1111"]
            .into_iter()
            .map(|x| x.to_string())
            .collect::<Vec<String>>();

        let alphabet_info = AlphabetInfo::from_sequences(&sequences, 4);
        let patterns = all_patterns_dp(&sequences, 4, &alphabet_info);
        assert_eq!(
            patterns,
            vec![
                "1001".to_string(),
                "1011".to_string(),
                "1101".to_string(),
                "1111".to_string(),
            ]
        );
        let entropy = calc_entropy(&sequences, 4);
        assert_approx_eq!(entropy, 1.9, 0.01);

        let sequences = vec!["*010", "1010", "0010"]
            .into_iter()
            .map(|x| x.to_string())
            .collect::<Vec<String>>();

        let alphabet_info = AlphabetInfo::from_sequences(&sequences, 4);
        let patterns = all_patterns_dp(&sequences, 4, &alphabet_info);
        assert_eq!(patterns, vec!["0010".to_string(), "1010".to_string(),]);
        let entropy = calc_entropy(&sequences, 4);
        assert_eq!(entropy, 1.0f32);

        let sequences = vec!["1010", "1010", "1010", "1010"]
            .into_iter()
            .map(|x| x.to_string())
            .collect::<Vec<String>>();

        let _alphabet_info = AlphabetInfo::from_sequences(&sequences, 4);
        let entropy = calc_entropy(&sequences, 4);
        assert_eq!(entropy, 0f32);
    }

    #[test]
    #[should_panic]
    fn test_alphabet_info() {
        // test that assert fires when columns are all *
        let sequences = vec![
            "*111".to_string(),
            "*111".to_string(),
            "*111".to_string(),
            "*111".to_string(),
        ];
        AlphabetInfo::from_sequences(&sequences, 4);
    }
}
