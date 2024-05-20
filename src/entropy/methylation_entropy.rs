use derive_new::new;
use std::collections::{BTreeSet, HashMap, HashSet};

use itertools::Itertools;
use log_once::debug_once;
use regex::Regex;
use rustc_hash::FxHashSet;
use substring::Substring;

#[inline]
fn seq_to_regex(seq: &str, alphabet: &str) -> Regex {
    let pattern = seq
        .chars()
        .map(|x| match x {
            '*' => format!("{alphabet}"),
            _ => format!("{x}"),
        })
        .collect::<String>();
    Regex::new(&pattern).unwrap()
}

#[derive(new)]
struct AlphabetInfo {
    alphabet: BTreeSet<char>,
    wildcard_regex: String,
}

impl AlphabetInfo {
    fn from_sequences(sequences: &[String]) -> Self {
        let alphabet = sequences
            .iter()
            .flat_map(|seq| {
                seq.chars().filter(|&x| x != '*').collect::<FxHashSet<char>>()
            })
            .collect::<BTreeSet<char>>();
        let wildcard_regex = {
            let alpha = alphabet.iter().sorted().join("");
            format!("[{alpha}]")
        };
        Self::new(alphabet, wildcard_regex)
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
        "all sequences should be the same length"
    );
    // easy case
    if !sequences.iter().any(|x| x.contains('*')) {
        return sequences.iter().map(|x| x.to_string()).collect();
    }

    let basecase = sequences
        .iter()
        .map(|x| x.substring(0, 1).to_string())
        .filter(|x| x != "*")
        .collect::<HashSet<String>>();

    let all_combs = (1..window_size).fold(basecase, |acc, idx| {
        let mut acc_patterns = HashSet::new();
        for seq in sequences.iter() {
            let prefix = seq.substring(0, idx);
            let re = seq_to_regex(prefix, &alphabet_info.wildcard_regex);
            for pattern in acc.iter() {
                if re.is_match(&pattern) {
                    let last_letter = seq
                        .chars()
                        .nth(idx)
                        .expect(&format!("should get last letter at {idx}"));
                    match last_letter {
                        '*' => {
                            for x in &alphabet_info.alphabet {
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
    let alphabet_info = AlphabetInfo::from_sequences(sequences);
    let patterns = all_patterns_dp(sequences, window_size, &alphabet_info);

    let counts = sequences.iter().fold(HashMap::new(), |mut acc, seq| {
        let re = seq_to_regex(seq, &alphabet_info.wildcard_regex);
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
    use std::collections::VecDeque;

    fn all_sequences(k: usize, alphabet: &[char]) -> Vec<String> {
        all_sequences_recur(k, vec![], alphabet)
    }

    fn all_sequences_recur(
        k: usize,
        subsequences: Vec<String>,
        alphabet: &[char],
    ) -> Vec<String> {
        if k < 1 {
            vec![]
        } else if k == 1 {
            alphabet.iter().map(|c| format!("{c}")).collect()
        } else {
            let subsequences =
                all_sequences_recur(k - 1, subsequences, alphabet);
            subsequences
                .into_iter()
                .flat_map(|x| {
                    alphabet
                        .iter()
                        .map(|c| format!("{x}{c}"))
                        .collect::<Vec<String>>()
                })
                .collect::<Vec<String>>()
        }
    }

    struct WildcardSequence<'a> {
        base_sequence: &'a str,
        patterns: VecDeque<Vec<char>>,
        wc_positions: Vec<usize>,
    }

    impl<'a> WildcardSequence<'a> {
        fn new(base_sequence: &'a str, alphabet: &[char]) -> Self {
            let wc_positions = base_sequence
                .char_indices()
                .filter_map(|(pos, c)| if c == '*' { Some(pos) } else { None })
                .collect::<Vec<usize>>();
            let patterns = all_sequences(wc_positions.len(), alphabet)
                .into_iter()
                .map(|x| x.chars().collect::<Vec<char>>())
                .collect::<VecDeque<Vec<char>>>();
            Self { base_sequence, patterns, wc_positions }
        }
    }

    impl<'a> Iterator for WildcardSequence<'a> {
        type Item = String;

        fn next(&mut self) -> Option<Self::Item> {
            if let Some(pattern) = self.patterns.pop_front() {
                assert_eq!(pattern.len(), self.wc_positions.len());
                let mut seq = self.base_sequence.chars().collect::<Vec<char>>();
                for (i, &pos) in self.wc_positions.iter().enumerate() {
                    seq[pos] = pattern[i];
                }
                Some(seq.into_iter().collect::<String>())
            } else {
                None
            }
        }
    }
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
        // let sequences = vec![
        //     "1101".to_string(),
        //     "1101".to_string(),
        //     "0*11".to_string(),
        //     "01*1".to_string(),
        // ];
        // let patterns = all_patterns_dp(&sequences, 4);
        // assert_eq!(patterns, vec!["0111".to_string(), "1101".to_string()]);
        // assert_eq!(calc_entropy(&sequences, 4), 1.0);

        let sequences = vec![
            "2222", "2222", "2221", "221*", "2*21", "1212", "1221", "*221",
            "*221", "2221",
        ]
        .into_iter()
        .map(|x| x.to_string())
        .collect::<Vec<String>>();

        let alphabet_info = AlphabetInfo::from_sequences(&sequences);
        let patterns = all_patterns_dp(&sequences, 4, &alphabet_info);

        let entropy = calc_entropy(&sequences, 4);
        // assert_eq!(entropy, 2.439354);
        assert_approx_eq!(entropy, 2.439354, 0.001)
    }

    #[test]
    fn test_patterns_bf() {
        let alphabet = ['0', '1'];

        let seq = "0*11*1";
        let wc_seq = WildcardSequence::new(seq, &alphabet);
        let expanded = wc_seq.into_iter().collect::<Vec<String>>();

        assert_eq!(expanded, vec!["001101", "001111", "011101", "011111",])
    }
}
