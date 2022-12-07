use crate::errs::InputError;
use indexmap::{indexset, IndexSet};
use std::collections::HashMap;

#[derive(Debug, PartialEq, Eq)]
enum Strand {
    Positive,
    Negative,
}

impl Strand {
    fn parse(raw: &str) -> Option<Self> {
        match raw {
            "+" => Some(Self::Positive),
            "-" => Some(Self::Negative),
            _ => None,
        }
    }

    fn parse_char(x: char) -> Result<Self, InputError> {
        match x {
            '+' => Ok(Self::Positive),
            '-' => Ok(Self::Negative),
            _ => Err(format!("failed to parse strand {}", x).into()),
        }
    }
}

#[derive(Debug, Eq, PartialEq)]
enum SkipMode {
    Ambiguous,
    ProbModified,
}

impl SkipMode {
    fn parse(raw_mode: char) -> Result<Self, InputError> {
        match raw_mode {
            '?' => Ok(Self::Ambiguous),
            '.' => Ok(Self::ProbModified),
            _ => Err(InputError::new(&format!("unknown mode {}", raw_mode))),
        }
    }
}

#[derive(Debug)]
pub struct BaseModProbs {
    mod_codes: IndexSet<char>,
    probs: Vec<f32>,
}

impl BaseModProbs {
    fn new(mod_code: char, prob: f32) -> Self {
        let mod_codes = indexset! {mod_code};
        let probs = vec![prob];
        Self { mod_codes, probs }
    }

    fn insert_base_mod_prob(&mut self, mod_code: char, prob: f32) {
        if let Some(idx) = self.mod_codes.get_index_of(&mod_code) {
            self.probs[idx] += prob;
        } else {
            self.mod_codes.insert(mod_code);
            self.probs.push(prob);
        }
    }

    fn collapse(self, mod_to_collapse: char) -> BaseModProbs {
        let canonical_prob = 1f32 - self.probs.iter().sum::<f32>();
        let marginal_collapsed_prob = self
            .mod_codes
            .iter()
            .zip(self.probs.iter())
            .filter(|(mod_code, prob)| **mod_code != mod_to_collapse)
            .collect::<Vec<(&char, &f32)>>();
        let total_marginal_collapsed_prob =
            marginal_collapsed_prob.iter().map(|(_, p)| *p).sum::<f32>() + canonical_prob;

        let mut mod_codes = IndexSet::new();
        let mut probs = Vec::new();
        for (mod_code, mod_prob) in marginal_collapsed_prob {
            let collapsed_prob = mod_prob / total_marginal_collapsed_prob;
            mod_codes.insert(*mod_code);
            probs.push(collapsed_prob)
        }

        Self { mod_codes, probs }
    }

    fn combine(&mut self, other: Self) {
        for (mod_code, prob) in other.mod_codes.into_iter().zip(other.probs.into_iter()) {
            if let Some(idx) = self.mod_codes.get_index_of(&mod_code) {
                self.probs[idx] += prob
            } else {
                self.mod_codes.insert(mod_code);
                self.probs.push(prob);
            }
        }
    }
}

pub struct DeltaListConverter {
    cumulative_counts: Vec<u32>,
}

impl DeltaListConverter {
    pub fn new(read_sequence: &str, base: char) -> Self {
        let cumulative_counts = read_sequence
            .chars()
            .scan(0, |count, nt| {
                if nt == base {
                    *count = *count + 1;
                }
                Some(*count)
            })
            .collect::<Vec<u32>>();

        assert_eq!(cumulative_counts.len(), read_sequence.len());
        Self { cumulative_counts }
    }

    pub fn to_positions(&self, delta_list: &[u32]) -> Vec<usize> {
        let mut finger = 0usize;
        let mut n_skips = 0u32;
        let mut positions = Vec::with_capacity(delta_list.len());
        for d in delta_list {
            while self.cumulative_counts[finger] <= (*d + n_skips) {
                finger += 1;
                assert!(
                    finger < self.cumulative_counts.len(),
                    "{:?} >= {:?},\ndelta_list: {:?}\ncumulative counts: {:?}",
                    finger,
                    self.cumulative_counts.len(),
                    delta_list,
                    self.cumulative_counts
                );
            }
            positions.push(finger);
            n_skips += d + 1;
        }
        positions
    }

    pub fn to_delta_list(&self, positions: &[usize]) -> Vec<u32> {
        let mut last = 0;
        let mut delta_list = Vec::new();
        for pos in positions {
            let cumulative_count = self.cumulative_counts[*pos];
            let d = cumulative_count - last - 1;
            delta_list.push(d);
            last = cumulative_count
        }
        delta_list
    }
}

#[inline]
fn qual_to_prob(qual: u16) -> f32 {
    let q = qual as f32;
    (q + 0.5f32) / 256f32
}

#[inline]
fn prob_to_qual(prob: f32) -> u16 {
    if prob == 1.0f32 {
        255u16
    } else {
        let p = prob * 256f32;
        let q = p.floor() as u16;
        assert!(q <= 255);
        q
    }
}

fn quals_to_probs(quals: &mut [f32]) {
    let arch = pulp::Arch::new();
    arch.dispatch(|| {
        for q in quals {
            let qual = *q as f32;
            *q = (qual + 0.5f32) / 256f32;
        }
    });
}

#[derive(Debug, Eq, PartialEq)]
struct BaseModPositions {
    canonical_base: char,
    mode: SkipMode,
    strand: Strand,
    mod_base_codes: Vec<char>,
    delta_list: Vec<u32>,
}

impl BaseModPositions {
    fn parse(mod_positions: &str) -> Result<Self, InputError> {
        let mut parts = mod_positions.split(',');
        let mut header = parts
            .nth(0)
            .ok_or(InputError::new(
                "failed to get leader for base mod position line",
            ))?
            .chars();

        let canonical_base = header
            .nth(0)
            .ok_or(InputError::new("failed to get canonical base"))?;

        let raw_stand = header
            .nth(0)
            .ok_or(InputError::new("failed to get strand"))?;

        // TODO handle duplex
        let strand = Strand::parse_char(raw_stand)?;

        let mut mod_base_codes = Vec::new();
        let mut mode: Option<SkipMode> = None;

        while let Some(c) = header.next() {
            match c {
                '?' | '.' => {
                    mode = Some(SkipMode::parse(c).unwrap());
                }
                _ => mod_base_codes.push(c),
            }
        }
        // default to the "old version"
        let mode = mode.unwrap_or(SkipMode::ProbModified);

        // taking the liberty to think that a read wouldn't be larger
        // than 2**32 - 1 bases long
        let delta_list = parts
            .into_iter()
            .map(|raw_pos| raw_pos.replace(";", "").parse::<u32>())
            .collect::<Result<Vec<u32>, _>>()
            .map_err(|e| {
                InputError::new(&format!("failed to parse position list, {}", e.to_string()))
            })?;

        Ok(Self {
            canonical_base,
            mod_base_codes,
            mode,
            strand,
            delta_list,
        })
    }
}

fn combine_positions_to_probs(
    agg: &mut HashMap<usize, BaseModProbs>,
    to_add: HashMap<usize, BaseModProbs>,
) {
    for (position, base_mod_probs) in to_add.into_iter() {
        if let Some(probs) = agg.get_mut(&position) {
            probs.combine(base_mod_probs);
        } else {
            agg.insert(position, base_mod_probs);
        }
    }
}

pub fn extract_mod_probs(
    raw_mm: &str,
    mod_quals: &[u16],
    canonical_base: char,
    converter: &DeltaListConverter,
) -> Result<HashMap<usize, BaseModProbs>, InputError> {
    let splited = raw_mm.split(";");
    let mut positions_to_probs = HashMap::new();
    let mut pointer = 0usize;
    for mod_positions in splited {
        if mod_positions.len() == 0 {
            continue;
        }
        let base_mod_positions = BaseModPositions::parse(mod_positions)?;
        if base_mod_positions.canonical_base == canonical_base {
            let base_mod_probs =
                get_base_mod_probs(&base_mod_positions, &mod_quals, pointer, converter).unwrap();
            combine_positions_to_probs(&mut positions_to_probs, base_mod_probs);
        }
        pointer += base_mod_positions.delta_list.len();
    }

    Ok(positions_to_probs)
}

fn get_base_mod_probs(
    base_mod_positions: &BaseModPositions,
    mod_quals: &[u16],
    pointer: usize,
    converter: &DeltaListConverter,
) -> Result<HashMap<usize, BaseModProbs>, InputError> {
    let positions = converter.to_positions(&base_mod_positions.delta_list);
    let probs = {
        let mut probs = mod_quals[pointer..pointer + base_mod_positions.delta_list.len()]
            .iter()
            .map(|qual| *qual as f32)
            .collect::<Vec<f32>>();
        quals_to_probs(&mut probs);
        probs
    };
    assert_eq!(probs.len(), positions.len());

    let mut positions_to_probs = HashMap::<usize, BaseModProbs>::new();
    for mod_base_code in &base_mod_positions.mod_base_codes {
        for (position, prob) in positions.iter().zip(&probs) {
            if let Some(base_mod_probs) = positions_to_probs.get_mut(position) {
                base_mod_probs.insert_base_mod_prob(*mod_base_code, *prob);
            } else {
                positions_to_probs.insert(*position, BaseModProbs::new(*mod_base_code, *prob));
            }
        }
    }

    Ok(positions_to_probs)
}

pub fn get_mod_probs_for_query_positions(
    mm: &str,
    canonical_base: char,
    mod_quals: &[u16],
    converter: &DeltaListConverter,
) -> Result<HashMap<usize, BaseModProbs>, InputError> {
    // todo move this outside this function. should handle the case where mods for another base
    //  come first and the offset of mod_quals is already handled
    let filtered_mod_positions = mm
        .split(';')
        .filter(|positions| positions.starts_with(canonical_base))
        .collect::<Vec<&str>>();

    let mut probs_for_positions = HashMap::<usize, BaseModProbs>::new();
    let mut prob_array_idx = 0usize;
    for mod_positions in filtered_mod_positions {
        let mut parts = mod_positions.split(',');
        let mut header = parts
            .nth(0)
            .ok_or(InputError::new(
                "failed to get leader for base mod position line",
            ))?
            .chars();

        let raw_stand = header
            .nth(1)
            .ok_or(InputError::new("failed to get strand"))?;

        // TODO handle duplex
        let _strand = Strand::parse_char(raw_stand);

        let mut mod_base_codes = Vec::new();
        let mut mode: Option<char> = None;

        while let Some(c) = header.next() {
            match c {
                '?' | '.' => {
                    mode = Some(c);
                }
                _ => mod_base_codes.push(c),
            }
        }

        // taking the liberty to think that a read wouldn't be larger
        // than 2**32 - 1 bases long
        let delta_list = parts
            .into_iter()
            .map(|raw_pos| raw_pos.parse::<u32>())
            .collect::<Result<Vec<u32>, _>>()
            .map_err(|e| {
                InputError::new(&format!("failed to parse position list, {}", e.to_string()))
            })?;

        let positions = converter.to_positions(&delta_list);
        for mod_base in mod_base_codes {
            for pos in &positions {
                let qual = mod_quals[prob_array_idx];
                let prob = qual_to_prob(qual);
                if let Some(base_mod_probs) = probs_for_positions.get_mut(pos) {
                    base_mod_probs.insert_base_mod_prob(mod_base, prob);
                } else {
                    probs_for_positions.insert(*pos, BaseModProbs::new(mod_base, prob));
                }
                // consume from the ML array
                prob_array_idx += 1;
            }
        }
    }

    Ok(probs_for_positions)
}

pub fn collapse_mod_probs(
    positions_to_probs: HashMap<usize, BaseModProbs>,
    mod_base_to_remove: char,
) -> HashMap<usize, BaseModProbs> {
    // let mut collapsed_probs = HashMap::new();
    positions_to_probs
        .into_iter()
        .map(|(pos, mod_base_probs)| (pos, mod_base_probs.collapse(mod_base_to_remove)))
        .collect()
}

pub fn format_mm_ml_tag(
    positions_to_probs: HashMap<usize, BaseModProbs>,
    canonical_base: char,
    converter: &DeltaListConverter,
) -> (String, Vec<u16>) {
    let mut mod_code_to_position = HashMap::new();
    for (position, mod_base_probs) in positions_to_probs {
        for (mod_base_code, mod_base_prob) in mod_base_probs
            .mod_codes
            .iter()
            .zip(mod_base_probs.probs.iter())
        {
            let entry = mod_code_to_position
                .entry(*mod_base_code)
                .or_insert(Vec::new());
            entry.push((position, *mod_base_prob));
        }
    }

    let mut mm_tag = String::new();
    let mut ml_tag = Vec::new();

    for (mod_code, mut positions_and_probs) in mod_code_to_position.into_iter() {
        positions_and_probs.sort_by(|(x_pos, _), (y_pos, _)| x_pos.cmp(&y_pos));
        let header = format!("{}+{}?,", canonical_base, mod_code);
        let positions = positions_and_probs
            .iter()
            .map(|(pos, _prob)| *pos)
            .collect::<Vec<usize>>();
        let delta_list = converter.to_delta_list(&positions);
        let delta_list = delta_list
            .into_iter()
            .map(|d| d.to_string())
            .collect::<Vec<String>>()
            .join(",");
        mm_tag.push_str(&header);
        mm_tag.push_str(&delta_list);
        mm_tag.push(';');
        let quals = positions_and_probs
            .iter()
            .map(|(_pos, prob)| prob_to_qual(*prob))
            .collect::<Vec<u16>>();
        ml_tag.extend(quals.into_iter());
    }

    (mm_tag, ml_tag)
}

// fn parse_mod_probs(record: &bam::Record) {
//     let mm = record.aux("MM".as_bytes()).unwrap();
//     let ml = record.aux("ML".as_bytes()).unwrap();
//     dbg!(mm);
//     dbg!(ml);
// }

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam;
    use rust_htslib::bam::{Read, Reader};

    #[test]
    fn test_delta_list_to_positions() {
        let canonical_base = 'C';
        let read_sequence = "ACCGCCGTCGTCG";
        let converter = DeltaListConverter::new(read_sequence, canonical_base);

        let ds = [1, 1, 0];
        let expected = [2, 5, 8];
        let obs = converter.to_positions(&ds);
        assert_eq!(&obs, &expected);
        let obs = converter.to_delta_list(&expected);
        assert_eq!(obs, &ds);

        let ds = [3, 0, 0];
        let expected = [5, 8, 11];
        let obs = converter.to_positions(&ds);
        assert_eq!(&obs, &expected);
        let obs = converter.to_delta_list(&expected);
        assert_eq!(obs, &ds);

        let ds = [3, 1];
        let expected = [5, 11];
        let obs = converter.to_positions(&ds);
        assert_eq!(&obs, &expected);
        let obs = converter.to_delta_list(&expected);
        assert_eq!(obs, &ds);
    }

    #[test]
    fn test_mod_prob_collapse() {
        let mod_base_probs = BaseModProbs {
            mod_codes: indexset! {'h', 'm'},
            probs: vec![0.05273438, 0.03320312],
        };
        let collapsed = mod_base_probs.collapse('h');
        assert_eq!(collapsed.probs, vec![0.035051543,]);
        assert_eq!(collapsed.mod_codes, indexset! {'m'});
    }

    #[test]
    fn test_mod_prob_combine() {
        let mut a = BaseModProbs {
            mod_codes: indexset! {'h', 'm'},
            probs: vec![0.05273438, 0.03320312],
        };
        let b = BaseModProbs {
            mod_codes: indexset! {'m'},
            probs: vec![0.03320312],
        };
        a.combine(b);
        assert_eq!(&a.probs, &[0.05273438, 0.06640624]);

        let mut a = BaseModProbs {
            mod_codes: indexset! {'m'},
            probs: vec![0.03320312],
        };

        let b = BaseModProbs {
            mod_codes: indexset! {'h'},
            probs: vec![0.05273438],
        };
        a.combine(b);
        assert_eq!(&a.probs, &[0.03320312, 0.05273438]);
        assert_eq!(&a.mod_codes, &indexset! {'m', 'h'});
    }

    #[test]
    fn test_parse_mm_tag() {
        let tag = "C+h?,5,2,1,3,1,2,3,1,2,1,11,5;C+m?,5,2,1,3,1,2,3,1,2,1,11,5;";
        let dna = "ATGTGCCTGCTGGACATGTTTATGCTCGTCTACTTCGTTCAGTTACGTATTGCTCCAG\
            CGCTCGAACTGTAGCCGCTGCTGCTGGGTGAAGTTGTGGCGGTACACGAGCTCCGCCGGCTGCAGCAGCTTC\
            TCCCCATCCTGGCGCTTCTCCCCGAGCAATTGGTG";
        let mod_quals = vec![
            197, 13, 156, 1, 3, 5, 9, 26, 8, 1, 0, 13, 10, 67, 1, 0, 1, 0, 5, 5, 5, 0, 0, 8,
        ];

        let converter = DeltaListConverter::new(dna, 'C');
        let positions_to_probs =
            get_mod_probs_for_query_positions(tag, 'C', &mod_quals, &converter).unwrap();
        assert_eq!(positions_to_probs.len(), 12);
    }

    #[test]
    fn test_mod_probs_to_tags() {
        let canonical_base = 'C';
        let read_sequence = "ACCGCCGTCGTCG";
        let converter = DeltaListConverter::new(read_sequence, canonical_base);

        let positions_and_probs = vec![
            (5, BaseModProbs::new('m', 0.9)),
            (2, BaseModProbs::new('m', 0.1)),
            (8, BaseModProbs::new('m', 0.2)),
        ]
        .into_iter()
        .collect::<HashMap<usize, BaseModProbs>>();

        let (mm, ml) = format_mm_ml_tag(positions_and_probs, 'C', &converter);
        assert_eq!(mm, "C+m?,1,1,0;");
        assert_eq!(ml, vec![25, 230, 51,]);
    }

    #[test]
    fn test_mod_base_positions() {
        let raw_positions = "C+h?,5,2,1,3,1,2,3,1,2,1,11,5;";
        let base_mod_positions = BaseModPositions::parse(raw_positions).unwrap();
        let expected = BaseModPositions {
            canonical_base: 'C',
            mode: SkipMode::Ambiguous,
            strand: Strand::Positive,
            mod_base_codes: vec!['h'],
            delta_list: vec![5, 2, 1, 3, 1, 2, 3, 1, 2, 1, 11, 5],
        };

        assert_eq!(base_mod_positions, expected);
    }

    #[test]
    fn test_quals_to_probs() {
        let mut quals = vec![100f32; 10_000];
        let start_time = std::time::Instant::now();
        let probs2 = quals_to_probs(&mut quals);
        let elapsed = start_time.elapsed();
        dbg!(elapsed);

        let quals = vec![100u16; 10_000];
        let start_time = std::time::Instant::now();
        let probs1 = quals.iter().map(|q| qual_to_prob(*q)).collect::<Vec<_>>();
        let elapsed = start_time.elapsed();
        dbg!(elapsed);
    }

    #[test]
    fn test_extract_positions_to_probs() {
        let dna = "GATCGACTACGTCGA";
        let tag = "C+h?,0,1,0;A+a?,0,1,0;C+m?,0,1,0;";
        let quals = vec![1, 1, 1, 200, 200, 200, 1, 1, 1];
        let canonical_base = 'C';
        let converter = DeltaListConverter::new(dna, canonical_base);

        let positions_to_probs =
            extract_mod_probs(tag, &quals, canonical_base, &converter).unwrap();
        assert_eq!(positions_to_probs.len(), 3);
        for (_pos, base_mod_probs) in positions_to_probs.iter() {
            assert_eq!(&base_mod_probs.probs, &[0.005859375, 0.005859375]);
            assert_eq!(&base_mod_probs.mod_codes, &indexset! { 'h', 'm' });
        }

        let tag = "C+hm?,0,1,0;A+a?,0,1,0;";
        let quals = vec![1, 1, 1, 1, 1, 1, 200, 200, 200];
        let positions_to_probs_comb =
            extract_mod_probs(tag, &quals, canonical_base, &converter).unwrap();
        assert_eq!(positions_to_probs_comb.len(), positions_to_probs.len());
        for (position, base_mod_probs) in positions_to_probs_comb.iter() {
            let other = positions_to_probs.get(position).unwrap();
            assert_eq!(base_mod_probs.probs, other.probs);
            assert_eq!(base_mod_probs.mod_codes, other.mod_codes);
        }
    }

    // #[test]
    // fn test_iter_records() {
    //     let fp = "data/hac_v4_C_20220804_1429_X3_FAT62446_e21856a4.bam";
    //     // let fp = "data/tmp_old_tags.bam";
    //     // let fp = "data/tmp.sorted.bam";
    //     let mut reader = Reader::from_path(fp).unwrap();
    //
    //     let alpha = "Chm";
    //     let n = alpha.len();
    //
    //     for (i, result) in reader.records().enumerate() {
    //         let record = result.unwrap();
    //         let seq = if record.is_reverse() {
    //             bio::alphabets::dna::revcomp(record.seq().as_bytes())
    //         } else {
    //             record.seq().as_bytes()
    //         };
    //         let dna = String::from_utf8(seq).unwrap();
    //
    //
    //
    //         dbg!(record.is_reverse());
    //         dbg!(dna);
    //         // parse_mod_probs(&record);
    //         break;
    //         // if i > 1 {
    //         //     break;
    //         // }
    //     }
    //
    //     println!("done");
    // }
}
