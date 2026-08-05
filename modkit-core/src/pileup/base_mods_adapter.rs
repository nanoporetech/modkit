use anyhow::{anyhow, bail};
use bitvec::{order::Lsb0, view::BitView};
use derive_new::new;
use memchr::{memchr, memchr_iter};
use rust_htslib::bam::{
    self,
    ext::{BamRecordExtensions, IterAlignedPairs},
    record::AuxArray,
};
use smallvec::SmallVec;

use crate::{
    errs::MkResult,
    mod_base_code::{DnaBase, ModCodeRepr},
    util::qual_to_prob,
};

#[derive(Debug, Copy, Clone, new)]
pub(crate) struct ModState {
    pub mod_position: usize,
    pub modified: bool,
    pub filtered: bool,
    pub mod_code: ModCodeRepr,
    pub primary_base: DnaBase,
    pub inferred: bool,
    pub mod_qual: u8,
}

#[derive(Debug, Copy, Clone)]
struct ModCodeState {
    mm_pos: usize,
    ml_pos: usize,
    ml_stride: usize,
    mod_code: ModCodeRepr,
    mm_next: Option<u32>,
    num_explicit_positions: u32,
    canonical_base: u8,
    implicit: bool,
    explicit_at_position: bool,
}

#[derive(Debug)]
pub(crate) struct BaseModsAdapter<'a, const SIZE: usize = 16> {
    mm: &'a [u8],
    seq: bam::record::Seq<'a>,
    ml: AuxArray<'a, u8>,
    code_states: SmallVec<[ModCodeState; SIZE]>,
    reverse: bool,
    left_to_right_seq_pos: usize,
}

impl<'a, const SIZE: usize> BaseModsAdapter<'a, SIZE> {
    pub fn new(record: &'a bam::Record) -> anyhow::Result<Self> {
        let seq = record.seq();
        let mm = match record
            .aux("MM".as_bytes())
            .or(record.aux("Mm".as_bytes()))?
        {
            bam::record::Aux::String(x) => x.as_bytes(),
            _ => bail!("MM tag must be a string"),
        };
        if mm.is_empty() {
            bail!("no mods")
        }
        let ml = match record
            .aux("ML".as_bytes())
            .or(record.aux("Ml".as_bytes()))?
        {
            bam::record::Aux::ArrayU8(aux_array) => aux_array,
            _ => bail!("ML tag must be an array"),
        };

        let reverse = record.is_reverse();
        let freq = if reverse {
            // [A, C, G, T]
            //  0  1  2  3
            (0..record.seq_len()).map(|idx| seq[idx]).fold(
                [0u32; 4],
                |mut agg, next| {
                    match next {
                        b'A' | b'a' => agg[3] += 1u32,
                        b'C' | b'c' => agg[2] += 1u32,
                        b'G' | b'g' => agg[1] += 1u32,
                        b'T' | b't' => agg[0] += 1u32,
                        _ => unreachable!(),
                    }
                    agg
                },
            )
        } else {
            [0u32; 4]
        };

        let mut i = 0;
        let mut ml_start = 0usize;
        let mut code_states = SmallVec::<[ModCodeState; SIZE]>::new();
        while i < mm.len() {
            let base = mm[i];
            i += 1;
            let strand = mm[i];
            if strand != b'+' {
                bail!("duplex data not currently supported")
            }
            i += 1;
            let (parsed_mod_codes, offset) = parse_mod_code(&mm[i..]);
            let mods_in_rec = parsed_mod_codes.len();
            i += offset;
            assert!(i < mm.len());
            let implicit_mode = match mm[i] {
                // no mode provided
                b';' | b',' => true,
                b'?' => {
                    i += 1;
                    false
                }
                b'.' => {
                    i += 1;
                    true
                }
                _ => bail!("malformed mm"),
            };
            // now we've moved into the list of deltas (if there are any)
            let record_end =
                memchr(b';', &mm[i..]).ok_or(anyhow!("missing end"))?;

            let (delta, mm_idx, n_deltas, ml_idx) = if record_end > 0 {
                if !reverse {
                    let n_deltas =
                        memchr_iter(b',', &mm[i..(record_end + i)]).count();
                    assert!(n_deltas > 0);
                    let mut j = i;
                    // this is unnecessary..
                    if mm[j] == b',' {
                        j += 1;
                    }
                    let (delta, offset) = parse_int::<b',', b';'>(&mm[j..]);
                    (Some(delta), offset + j, n_deltas as u32, ml_start)
                } else {
                    // consume the mm string forwards and calculate the cum-sum
                    // of accounted bases
                    let mut j = i;
                    let mut count = 0u32;
                    let mut n_deltas = 0u32;
                    while j < (record_end + i) {
                        let (d, o) = parse_int::<b',', b';'>(&mm[j..]);
                        count += d + 1;
                        j += o;
                        n_deltas += 1;
                    }
                    assert!(n_deltas > 0);
                    let total = match base {
                        b'A' => freq[0],
                        b'C' => freq[1],
                        b'G' => freq[2],
                        b'T' => freq[3],
                        _ => unreachable!(),
                    };
                    let remainder = total.checked_sub(count).expect(&format!(
                        "count should be less than or equal total, MN tag \
                         incorrect?"
                    ));
                    let ml_idx =
                        ml_start + ((n_deltas - 1) as usize * mods_in_rec);
                    (Some(remainder), record_end + i, n_deltas, ml_idx)
                }
            } else {
                (None, 0usize, 0u32, 0usize)
            };

            for j in 0..mods_in_rec {
                code_states.push(ModCodeState {
                    mm_pos: mm_idx,
                    ml_pos: ml_idx + j,
                    ml_stride: mods_in_rec,
                    mod_code: parsed_mod_codes.get(j),
                    mm_next: delta,
                    num_explicit_positions: n_deltas,
                    canonical_base: base,
                    implicit: implicit_mode,
                    explicit_at_position: false,
                });
            }
            ml_start += n_deltas as usize * mods_in_rec;
            i += record_end + 1;
        }

        Ok(Self { mm, ml, seq, code_states, reverse, left_to_right_seq_pos: 0 })
    }

    pub fn next_modified_position_no_thresh(
        &mut self,
    ) -> MkResult<Option<ModState>> {
        self.next_modified_position([0f32; 4], &Vec::new())
    }

    pub fn next_modified_position(
        &mut self,
        filter_thresholds: [f32; 4],
        mod_thresholds: &[(ModCodeRepr, f32)],
    ) -> MkResult<Option<ModState>> {
        let mut mod_pos = None;
        let mut pos = self.left_to_right_seq_pos;
        let mut done = false;
        let mut inferred = true;
        // A different MM group can count down to zero while this position is
        // scanned. Only code states already at zero have an explicit call here.
        // Store readiness with each dynamically sized state so this scan does
        // not allocate, including when the inline capacity is exceeded.
        self.code_states
            .iter_mut()
            .for_each(|state| state.explicit_at_position = false);

        while !done && pos < self.seq.len() {
            let base = if self.reverse {
                base_complement(self.seq[pos])
            } else {
                self.seq[pos]
            };
            for code_state in &mut self.code_states {
                match code_state.mm_next {
                    Some(0) => {
                        // there is at least one explicit modification here
                        if code_state.canonical_base == base {
                            done = true;
                            mod_pos = Some(pos);
                            inferred = false;
                            code_state.explicit_at_position = true;
                        }
                    }
                    Some(skip_count) => {
                        if code_state.canonical_base == base {
                            if code_state.implicit {
                                done = true;
                                mod_pos = Some(pos);
                            } else {
                                code_state.mm_next = Some(
                                    skip_count
                                        .checked_sub(1u32)
                                        .expect("should not go off the end"),
                                );
                            }
                        }
                    }
                    None => {
                        if code_state.canonical_base == base {
                            if code_state.implicit {
                                done = true;
                                mod_pos = Some(pos);
                            }
                        }
                    }
                }
            }
            pos += 1;
        }
        let mod_state = if let Some(mod_pos) = mod_pos {
            let mut mod_qual = 0u8;
            let mut total_mod_qual = 0u8;
            let base = if self.reverse {
                base_complement(self.seq[mod_pos])
            } else {
                self.seq[mod_pos]
            };
            let mut mod_code = ModCodeRepr::Code(base as char);
            for code_state in &self.code_states {
                if code_state.canonical_base == base
                    && code_state.explicit_at_position
                    && code_state.mm_next.map(|x| x == 0).unwrap_or(false)
                {
                    let q = self.ml.get(code_state.ml_pos).unwrap();
                    if q > mod_qual {
                        mod_code = code_state.mod_code;
                        mod_qual = q;
                    }
                    total_mod_qual = total_mod_qual.saturating_add(q);
                }
            }
            let canonical_qual = 255u8.checked_sub(total_mod_qual).unwrap();
            let primary_base = DnaBase::parse(base as char).unwrap();
            let threshold = filter_thresholds[primary_base as usize];
            let mod_state = if canonical_qual > mod_qual {
                Some(ModState::new(
                    mod_pos,
                    false,
                    qual_to_prob(canonical_qual as i32) < threshold,
                    ModCodeRepr::Code(base as char),
                    primary_base,
                    inferred,
                    canonical_qual,
                ))
            } else {
                let mod_threshold = mod_thresholds
                    .iter()
                    .find_map(|(code, p)| {
                        if code == &mod_code || code.matches_any(&primary_base)
                        {
                            Some(*p)
                        } else {
                            None
                        }
                    })
                    .unwrap_or(threshold);
                Some(ModState::new(
                    mod_pos,
                    true,
                    qual_to_prob(mod_qual as i32) < mod_threshold,
                    mod_code,
                    primary_base,
                    inferred,
                    mod_qual,
                ))
            };

            self.move_forward(mod_pos, base);

            mod_state
        } else {
            None
        };

        Ok(mod_state)
    }

    pub fn primary_bases_in_record(&self) -> u8 {
        let mut bs = 0u8;
        for code_state in &self.code_states {
            match code_state.canonical_base {
                b'A' => bs.view_bits_mut::<Lsb0>().set(0, true),
                b'C' => bs.view_bits_mut::<Lsb0>().set(1, true),
                b'G' => bs.view_bits_mut::<Lsb0>().set(2, true),
                b'T' => bs.view_bits_mut::<Lsb0>().set(3, true),
                _ => unreachable!(),
            }
        }
        bs
    }

    #[inline]
    fn move_forward(&mut self, last_pos: usize, base: u8) {
        self.left_to_right_seq_pos = last_pos.saturating_add(1);
        for code_state in self
            .code_states
            .iter_mut()
            .filter(|state| state.canonical_base == base)
        {
            match code_state.mm_next {
                Some(0) if code_state.explicit_at_position => {
                    let num_explicit_positions =
                        code_state.num_explicit_positions.saturating_sub(1);
                    if num_explicit_positions == 0 {
                        // done
                        code_state.mm_next = None;
                    } else {
                        let mm_pos = code_state.mm_pos;
                        if self.reverse {
                            let mut p = mm_pos.saturating_sub(1);
                            while self.mm[p] != b',' {
                                p = p.saturating_sub(1);
                            }
                            let mut val = 0;
                            for idx in (p + 1)..mm_pos {
                                let c = self.mm[idx];
                                val = val * 10 + (c - b'0') as u32;
                            }
                            code_state.mm_next = Some(val);
                            code_state.mm_pos = p;
                            code_state.ml_pos -= code_state.ml_stride;
                            // if num_explicit_positions > 1 {
                            // }
                        } else {
                            let (delta, offset) =
                                parse_int::<b',', b';'>(&self.mm[mm_pos..]);
                            code_state.mm_next = Some(delta);
                            code_state.mm_pos = offset + mm_pos;
                            code_state.ml_pos += code_state.ml_stride;
                            // if num_explicit_positions > 1 {
                            // }
                        }
                        code_state.num_explicit_positions =
                            num_explicit_positions;
                    }
                }
                Some(skip_count) if code_state.implicit => {
                    code_state.mm_next = Some(
                        skip_count
                            .checked_sub(1u32)
                            .expect("should not go less than zero"),
                    );
                }
                _ => {}
            }
        }
    }
}

enum ParsedModCodes<'a> {
    Alphabetic(&'a [u8]),
    ChEbi(u32),
}

impl ParsedModCodes<'_> {
    fn len(&self) -> usize {
        match self {
            Self::Alphabetic(codes) => codes.len(),
            Self::ChEbi(_) => 1,
        }
    }

    fn get(&self, index: usize) -> ModCodeRepr {
        match self {
            Self::Alphabetic(codes) => ModCodeRepr::Code(codes[index] as char),
            Self::ChEbi(code) => {
                debug_assert_eq!(index, 0);
                ModCodeRepr::ChEbi(*code)
            }
        }
    }
}

#[inline(always)]
fn parse_mod_code(bytes: &[u8]) -> (ParsedModCodes<'_>, usize) {
    let mut idx = 0;
    let mut val: u32 = 0;

    let mut all_digits = true;
    while idx < bytes.len()
        && (bytes[idx] != b'.'
            && bytes[idx] != b'?'
            && bytes[idx] != b','
            && bytes[idx] != b';')
    {
        assert!(idx < bytes.len());
        let c = bytes[idx];
        if (b'0'..=b'9').contains(&c) {
            val = val * 10 + (c - b'0') as u32;
        } else {
            all_digits = false;
            break;
        }
        idx += 1;
    }

    if !all_digits {
        let start = idx;
        while bytes[idx] != b'.'
            && bytes[idx] != b'?'
            && bytes[idx] != b','
            && bytes[idx] != b';'
        {
            idx += 1;
        }

        (ParsedModCodes::Alphabetic(&bytes[start..idx]), idx)
    } else {
        (ParsedModCodes::ChEbi(val), idx)
    }
}

fn parse_int<const DELIM: u8, const END: u8>(bs: &[u8]) -> (u32, usize) {
    let mut val = 0;
    let mut idx = 0;
    if bs[idx] == DELIM {
        idx += 1;
    }
    loop {
        let c = bs[idx];
        if c == DELIM || c == END {
            break;
        }
        val = val * 10 + (c - b'0') as u32;
        idx += 1;
    }
    (val, idx)
}

fn base_complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => panic!("not allowed base"),
    }
}

pub(crate) struct AlignedBaseModsIterator<'a> {
    aligned_pairs_iter: IterAlignedPairs,
    modbase_iter: BaseModsAdapter<'a>,
    mod_state: Option<ModState>,
}

impl<'a> AlignedBaseModsIterator<'a> {
    pub(crate) fn new(record: &'a bam::Record) -> anyhow::Result<Self> {
        let aligned_pairs_iter = record.aligned_pairs();
        let mut modbase_iter = BaseModsAdapter::new(record)?;
        let mod_state = modbase_iter.next_modified_position_no_thresh()?;
        Ok(Self { aligned_pairs_iter, modbase_iter, mod_state })
    }

    fn next_res(&mut self) -> MkResult<Option<(u32, u32, ModState)>> {
        if self.mod_state.is_none() {
            return Ok(None);
        }
        'scan: loop {
            let Some((qpos, r_pos)) =
                self.aligned_pairs_iter.next().and_then(|x| {
                    if x[0] < 0i64 || x[1] < 0i64 {
                        None
                    } else {
                        Some((x[0] as u32, x[1] as u32))
                    }
                })
            else {
                return Ok(None);
            };
            let mod_pos = self.mod_state.unwrap().mod_position as u32;
            if mod_pos == qpos {
                return Ok(Some((qpos, r_pos, self.mod_state.unwrap())));
            } else if qpos < mod_pos {
                continue 'scan;
            } else {
                assert!(qpos > mod_pos);
                'advance_mods: loop {
                    self.mod_state =
                        self.modbase_iter.next_modified_position_no_thresh()?;
                    match self.mod_state.as_ref() {
                        Some(ms) => {
                            if ms.mod_position < qpos as usize {
                                continue 'advance_mods;
                            } else if ms.mod_position == qpos as usize {
                                let ret = Ok(Some((qpos, r_pos, *ms)));
                                self.mod_state = self
                                    .modbase_iter
                                    .next_modified_position_no_thresh()?;
                                return ret;
                            } else {
                                continue 'scan;
                            }
                        }
                        None => return Ok(None),
                    }
                }
            }
        }
    }

    pub fn primary_bases_in_record(&self) -> u8 {
        self.modbase_iter.primary_bases_in_record()
    }
}

impl<'a> Iterator for AlignedBaseModsIterator<'a> {
    type Item = MkResult<(u32, u32, ModState)>;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_res().transpose()
    }
}

#[cfg(test)]
mod base_mods_adapter_tests {
    use rust_htslib::bam::{self, record::Aux};

    use crate::mod_base_code::{
        ModCodeRepr, ANY_ADENINE, ANY_CYTOSINE, HYDROXY_METHYL_CYTOSINE,
        METHYL_CYTOSINE, SIX_METHYL_ADENINE,
    };

    use super::BaseModsAdapter;
    fn make_record(
        mm: &str,
        ml: &[u8],
        seq: &str,
        name: Option<&str>,
        reverse: bool,
    ) -> bam::Record {
        let qname = name.map(|x| x.as_bytes()).unwrap_or("r1".as_bytes());
        let mut record = bam::Record::new();
        record.set(qname, None, seq.as_bytes(), &vec![255u8; seq.len()]);
        record.push_aux("MM".as_bytes(), Aux::String(mm)).unwrap();
        record.push_aux("ML".as_bytes(), Aux::ArrayU8((&ml).into())).unwrap();
        if reverse {
            record.set_reverse();
        }
        record
    }

    fn make_forward_record_with_mn(
        mm: &str,
        ml: &[u8],
        seq: &str,
    ) -> bam::Record {
        let mut record = make_record(mm, ml, seq, None, false);
        record.push_aux(b"MN", Aux::I32(seq.len() as i32)).unwrap();
        record
    }

    fn collect_states<const SIZE: usize>(
        scanner: &mut BaseModsAdapter<SIZE>,
    ) -> Vec<(usize, ModCodeRepr, u8, bool, bool)> {
        std::iter::from_fn(|| {
            scanner.next_modified_position_no_thresh().unwrap()
        })
        .map(|state| {
            (
                state.mod_position,
                state.mod_code,
                state.mod_qual,
                state.modified,
                state.inferred,
            )
        })
        .collect()
    }

    fn alphabetic_codes(count: usize) -> String {
        (b'a'..).take(count).map(char::from).collect()
    }

    #[test]
    fn test_reverse() {
        let mm = "C+m.,0,0;C+h.,0,0;";
        let seq = "CGGATGGAGTC";
        let ml = vec![100, 25, 101, 24];
        let record = make_record(mm, &ml, seq, None, true);

        let mut scanner = BaseModsAdapter::<2>::new(&record).unwrap();
        let mod_state = std::iter::from_fn(|| {
            scanner.next_modified_position([0f32; 4], &[]).unwrap()
        })
        .find_map(
            |mod_state| {
                if !mod_state.inferred {
                    Some(mod_state)
                } else {
                    None
                }
            },
        )
        .unwrap();
        assert_eq!(mod_state.mod_position, 6);
        assert_eq!(mod_state.mod_qual, 255 - (25 + 24));
        assert!(!mod_state.modified);
        assert!(!mod_state.inferred);

        let mod_state = std::iter::from_fn(|| {
            scanner.next_modified_position([0f32; 4], &[]).unwrap()
        })
        .find_map(
            |mod_state| {
                if !mod_state.inferred {
                    Some(mod_state)
                } else {
                    None
                }
            },
        )
        .unwrap();
        assert_eq!(mod_state.mod_position, 8);
        assert!(!mod_state.inferred);
        assert!(mod_state.modified);
        assert_eq!(mod_state.mod_code, ModCodeRepr::Code('h'));
        let mod_state = scanner.next_modified_position([0f32; 4], &[]).unwrap();
        assert!(mod_state.is_none());
    }

    #[test]
    fn test_forward() {
        let seq = "ATCATCATTCCTACCGCTATAGCCT";
        let mm = "C+mh,2,0,1;";
        let ml = vec![200, 10, 50, 170, 160, 20];
        let record = make_record(mm, &ml, seq, None, false);
        let thresholds = [0f32; 4];
        let mut scanner = BaseModsAdapter::<2>::new(&record).unwrap();
        let mod_state = std::iter::from_fn(|| {
            scanner.next_modified_position(thresholds, &[]).unwrap()
        })
        .find_map(
            |mod_state| {
                if !mod_state.inferred {
                    Some(mod_state)
                } else {
                    None
                }
            },
        )
        .unwrap();
        assert_eq!(mod_state.mod_position, 9);
        assert_eq!(mod_state.mod_qual, 200);
        assert!(mod_state.modified);
        assert!(!mod_state.inferred);
        assert_eq!(mod_state.mod_code, METHYL_CYTOSINE);

        let mod_state = std::iter::from_fn(|| {
            scanner.next_modified_position(thresholds, &[]).unwrap()
        })
        .find_map(
            |mod_state| {
                if !mod_state.inferred {
                    Some(mod_state)
                } else {
                    None
                }
            },
        )
        .unwrap();
        assert_eq!(mod_state.mod_position, 10);
        assert_eq!(mod_state.mod_qual, 170);
        assert!(mod_state.modified);
        assert!(!mod_state.inferred);
        assert_eq!(mod_state.mod_code, HYDROXY_METHYL_CYTOSINE);

        let mod_state = std::iter::from_fn(|| {
            scanner.next_modified_position(thresholds, &[]).unwrap()
        })
        .find_map(
            |mod_state| {
                if !mod_state.inferred {
                    Some(mod_state)
                } else {
                    None
                }
            },
        )
        .unwrap();
        assert_eq!(mod_state.mod_position, 14);
        assert_eq!(mod_state.mod_qual, 160);
        assert!(mod_state.modified);
        assert!(!mod_state.inferred);
        assert_eq!(mod_state.mod_code, METHYL_CYTOSINE);
        let mod_state =
            scanner.next_modified_position(thresholds, &[]).unwrap().unwrap();
        assert_eq!(mod_state.mod_position, 16);
        assert_eq!(mod_state.mod_qual, 255);
        assert!(!mod_state.modified);
        assert!(mod_state.inferred);
    }

    #[test]
    fn test_independent_same_base_groups_keep_their_delta_progression() {
        let mut record =
            make_record("C+m?,0;C+h?,1;", &[200, 250], "CC", None, false);
        record.push_aux(b"MN", Aux::I32(2)).unwrap();

        let mut scanner = BaseModsAdapter::<2>::new(&record).unwrap();
        let mod_states = std::iter::from_fn(|| {
            scanner.next_modified_position_no_thresh().unwrap()
        })
        .map(|state| {
            (state.mod_position, state.mod_code, state.mod_qual, state.modified)
        })
        .collect::<Vec<_>>();

        assert_eq!(
            mod_states,
            vec![
                (0, METHYL_CYTOSINE, 200, true),
                (1, HYDROXY_METHYL_CYTOSINE, 250, true),
            ]
        );
    }

    #[test]
    fn test_implicit_calls() {
        let seq = "ATCATCATTCCTACCGCTATAGCCT";
        let mm = "C+mh.,2,0,1;";
        let ml = vec![200, 10, 50, 170, 160, 20];
        let record = make_record(mm, &ml, seq, None, false);
        let thresholds = [0f32; 4];
        let mut scanner = BaseModsAdapter::<2>::new(&record).unwrap();
        let mod_state =
            scanner.next_modified_position(thresholds, &[]).unwrap().unwrap();
        assert_eq!(mod_state.mod_position, 2);
        assert_eq!(mod_state.mod_qual, 255);
        assert!(!mod_state.modified);
        assert!(mod_state.inferred);
        let mod_state =
            scanner.next_modified_position(thresholds, &[]).unwrap().unwrap();
        assert_eq!(mod_state.mod_position, 5);
        assert_eq!(mod_state.mod_qual, 255);
        assert!(!mod_state.modified);
        assert!(mod_state.inferred);
        let mod_state =
            scanner.next_modified_position(thresholds, &[]).unwrap().unwrap();
        assert_eq!(mod_state.mod_position, 9);
        assert_eq!(mod_state.mod_qual, 200);
        assert!(mod_state.modified);
        assert!(!mod_state.inferred);
        let mod_state =
            scanner.next_modified_position(thresholds, &[]).unwrap().unwrap();
        assert_eq!(mod_state.mod_position, 10);
        assert_eq!(mod_state.mod_qual, 170);
        assert_eq!(mod_state.mod_code, HYDROXY_METHYL_CYTOSINE);
        assert!(mod_state.modified);
        assert!(!mod_state.inferred);
        let mod_state =
            scanner.next_modified_position(thresholds, &[]).unwrap().unwrap();
        assert_eq!(mod_state.mod_position, 13);
        assert_eq!(mod_state.mod_qual, 255);
        assert!(!mod_state.modified);
        assert!(mod_state.inferred);
        let mod_state =
            scanner.next_modified_position(thresholds, &[]).unwrap().unwrap();
        assert_eq!(mod_state.mod_position, 14);
        assert_eq!(mod_state.mod_qual, 160);
        assert_eq!(mod_state.mod_code, METHYL_CYTOSINE);
        assert!(mod_state.modified);
        assert!(!mod_state.inferred);
        let mod_state =
            scanner.next_modified_position(thresholds, &[]).unwrap().unwrap();
        assert_eq!(mod_state.mod_position, 16);
        assert_eq!(mod_state.mod_qual, 255);
        assert!(!mod_state.modified);
        assert!(mod_state.inferred);
        let mod_state =
            scanner.next_modified_position(thresholds, &[]).unwrap().unwrap();
        assert_eq!(mod_state.mod_position, 22);
        assert_eq!(mod_state.mod_qual, 255);
        assert!(!mod_state.modified);
        assert!(mod_state.inferred);
        let mod_state =
            scanner.next_modified_position(thresholds, &[]).unwrap().unwrap();
        assert_eq!(mod_state.mod_position, 23);
        assert_eq!(mod_state.mod_qual, 255);
        assert!(!mod_state.modified);
        assert!(mod_state.inferred);
    }

    #[test]
    fn test_implicit_reverse() {
        let mm = "C+mh.,1,0;";
        let seq = "GAGACGGA";
        let ml = vec![100, 25, 101, 24];

        let record = make_record(mm, &ml, seq, None, true);
        let thresholds = [0f32; 4];
        let mut scanner = BaseModsAdapter::<2>::new(&record).unwrap();
        assert_eq!(
            scanner
                .code_states
                .iter()
                .map(|state| state.ml_pos)
                .collect::<Vec<_>>(),
            vec![2, 3]
        );
        let mod_state =
            scanner.next_modified_position(thresholds, &[]).unwrap().unwrap();
        assert_eq!(mod_state.mod_position, 0);
        assert!(!mod_state.modified);
        assert!(mod_state.inferred);
        let mod_state =
            scanner.next_modified_position(thresholds, &[]).unwrap().unwrap();
        assert_eq!(mod_state.mod_position, 2);
        assert!(!mod_state.modified);
        assert!(!mod_state.inferred);
        assert_eq!(mod_state.mod_qual, 255u8 - (100 + 25));
        let mod_state =
            scanner.next_modified_position(thresholds, &[]).unwrap().unwrap();
        assert_eq!(mod_state.mod_position, 5);
        assert!(!mod_state.modified);
        assert!(!mod_state.inferred);
        assert_eq!(mod_state.mod_qual, 255u8 - (101 + 24));
        let mod_state =
            scanner.next_modified_position(thresholds, &[]).unwrap().unwrap();
        assert_eq!(mod_state.mod_position, 6);
        assert!(!mod_state.modified);
        assert!(mod_state.inferred);
    }

    #[test]
    fn test_empty() {
        let mm = "C+mh?;";
        let seq = "GAGACGGA";
        let ml = vec![];
        let thresholds = [0f32; 4];

        let record = make_record(mm, &ml, seq, None, true);
        let mut scanner = BaseModsAdapter::<2>::new(&record).unwrap();
        assert!(scanner
            .next_modified_position(thresholds, &[])
            .unwrap()
            .is_none());
        let mut scanner = BaseModsAdapter::<2>::new(&record).unwrap();
        assert!(scanner
            .next_modified_position(thresholds, &[])
            .unwrap()
            .is_none());

        let mm = "C+mh.;";
        let seq = "GAGACGGA";
        let ml = vec![];
        let record = make_record(mm, &ml, seq, None, false);
        let mut scanner = BaseModsAdapter::<2>::new(&record).unwrap();
        let mod_state =
            scanner.next_modified_position(thresholds, &[]).unwrap().unwrap();
        assert_eq!(mod_state.mod_position, 4);
        assert!(!mod_state.modified);
        assert!(mod_state.inferred);

        let record = make_record(mm, &ml, seq, None, true);
        let mut scanner = BaseModsAdapter::<2>::new(&record).unwrap();
        let mod_positions = std::iter::from_fn(|| {
            scanner.next_modified_position(thresholds, &[]).unwrap()
        })
        .map(|x| x.mod_position)
        .collect::<Vec<_>>();
        assert_eq!(&mod_positions, &[0, 2, 5, 6]);
    }

    #[test]
    fn test_multiple() {
        let mm = "C+mh.,2;A+a.,1,0;";
        let ml: Vec<u8> = vec![200, 10, 20, 250];
        let seq = "CAGCATCGAT";
        let record = make_record(mm, &ml, seq, None, false);
        let mut scanner = BaseModsAdapter::<3>::new(&record).unwrap();
        let thresholds = [0f32; 4];
        let mod_state =
            scanner.next_modified_position(thresholds, &[]).unwrap().unwrap();
        assert_eq!(mod_state.mod_position, 0);
        assert!(!mod_state.modified);
        assert!(mod_state.inferred);
        assert_eq!(mod_state.mod_code, ANY_CYTOSINE);
        let mod_state =
            scanner.next_modified_position(thresholds, &[]).unwrap().unwrap();
        assert_eq!(mod_state.mod_position, 1);
        assert!(!mod_state.modified);
        assert!(mod_state.inferred);
        assert_eq!(mod_state.mod_code, ANY_ADENINE);
        let mod_state =
            scanner.next_modified_position(thresholds, &[]).unwrap().unwrap();
        assert_eq!(mod_state.mod_position, 3);
        assert!(!mod_state.modified);
        assert!(mod_state.inferred);
        assert_eq!(mod_state.mod_code, ANY_CYTOSINE);
        let mod_state =
            scanner.next_modified_position(thresholds, &[]).unwrap().unwrap();
        assert_eq!(mod_state.mod_position, 4);
        assert!(!mod_state.modified);
        assert!(!mod_state.inferred);
        assert_eq!(mod_state.mod_code, ANY_ADENINE);
        assert_eq!(mod_state.mod_qual, 255u8 - 20);
        let mod_state =
            scanner.next_modified_position(thresholds, &[]).unwrap().unwrap();
        assert_eq!(mod_state.mod_position, 6);
        assert!(mod_state.modified);
        assert!(!mod_state.inferred);
        assert_eq!(mod_state.mod_code, METHYL_CYTOSINE);
        assert_eq!(mod_state.mod_qual, 200);
        let mod_state =
            scanner.next_modified_position(thresholds, &[]).unwrap().unwrap();
        assert_eq!(mod_state.mod_position, 8);
        assert!(mod_state.modified);
        assert!(!mod_state.inferred);
        assert_eq!(mod_state.mod_code, SIX_METHYL_ADENINE);
        assert_eq!(mod_state.mod_qual, 250);
    }

    #[test]
    fn test_multiple_reverse() {
        let mm = "C+mh.,2;A+a.,1,0;";
        let ml: Vec<u8> = vec![200, 10, 20, 250];
        let seq = "ATCGATGCTG";
        let record = make_record(mm, &ml, seq, None, true);
        let mut scanner = BaseModsAdapter::<3>::new(&record).unwrap();
        let thresholds = [0f32; 4];
        let mod_state =
            scanner.next_modified_position(thresholds, &[]).unwrap().unwrap();
        assert_eq!(mod_state.mod_position, 1);
        assert!(mod_state.modified);
        assert!(!mod_state.inferred);
        assert_eq!(mod_state.mod_code, SIX_METHYL_ADENINE);
        assert_eq!(mod_state.mod_qual, 250);
        let mod_state =
            scanner.next_modified_position(thresholds, &[]).unwrap().unwrap();
        assert_eq!(mod_state.mod_position, 3);
        assert!(mod_state.modified);
        assert!(!mod_state.inferred);
        assert_eq!(mod_state.mod_code, METHYL_CYTOSINE);
        assert_eq!(mod_state.mod_qual, 200);
        let mod_state =
            scanner.next_modified_position(thresholds, &[]).unwrap().unwrap();
        assert_eq!(mod_state.mod_position, 5);
        assert!(!mod_state.modified);
        assert!(!mod_state.inferred);
        assert_eq!(mod_state.mod_code, ANY_ADENINE);
        assert_eq!(mod_state.mod_qual, 235);
        let mod_state =
            scanner.next_modified_position(thresholds, &[]).unwrap().unwrap();
        assert_eq!(mod_state.mod_position, 6);
        assert!(!mod_state.modified);
        assert!(mod_state.inferred);
        assert_eq!(mod_state.mod_code, ANY_CYTOSINE);
        let mod_state =
            scanner.next_modified_position(thresholds, &[]).unwrap().unwrap();
        assert_eq!(mod_state.mod_position, 8);
        assert!(!mod_state.modified);
        assert!(mod_state.inferred);
        assert_eq!(mod_state.mod_code, ANY_ADENINE);
        let mod_state =
            scanner.next_modified_position(thresholds, &[]).unwrap().unwrap();
        assert_eq!(mod_state.mod_position, 9);
        assert!(!mod_state.modified);
        assert!(mod_state.inferred);
        assert_eq!(mod_state.mod_code, ANY_CYTOSINE);
    }

    #[test]
    fn test_single_group_grows_beyond_inline_code_capacity() {
        for count in [16, 17] {
            let codes = alphabetic_codes(count);
            let mm = format!("C+{codes}?,0;");
            let mut ml = vec![0; count];
            ml[count - 1] = 200;
            let record = make_forward_record_with_mn(&mm, &ml, "C");

            let mut scanner = BaseModsAdapter::<16>::new(&record).unwrap();
            assert_eq!(
                scanner
                    .code_states
                    .iter()
                    .map(|state| state.mod_code)
                    .collect::<Vec<_>>(),
                codes
                    .bytes()
                    .map(|code| ModCodeRepr::Code(code as char))
                    .collect::<Vec<_>>()
            );
            assert!(scanner
                .code_states
                .iter()
                .all(|state| state.ml_stride == count));
            assert_eq!(scanner.code_states.spilled(), count > 16);
            let expected_code =
                ModCodeRepr::Code(codes.as_bytes()[count - 1] as char);
            assert_eq!(
                collect_states(&mut scanner),
                vec![(0, expected_code, 200, true, false)]
            );
        }
    }

    #[test]
    fn test_spilled_single_group_reverse_cursor_and_ml_stride() {
        let codes = alphabetic_codes(17);
        let mm = format!("C+{codes}?,0,0;");
        let mut ml = vec![0; 34];
        ml[0] = 200;
        ml[33] = 220;
        let mut record = make_record(&mm, &ml, "GG", None, true);
        record.push_aux(b"MN", Aux::I32(2)).unwrap();

        let mut scanner = BaseModsAdapter::<16>::new(&record).unwrap();
        assert!(scanner.code_states.spilled());
        assert_eq!(scanner.code_states.len(), 17);
        assert!(scanner.code_states.iter().all(|state| state.ml_stride == 17));
        assert_eq!(
            scanner
                .code_states
                .iter()
                .map(|state| state.ml_pos)
                .collect::<Vec<_>>(),
            (17..34).collect::<Vec<_>>()
        );
        let initial_mm_pos = scanner.code_states[0].mm_pos;

        let state =
            scanner.next_modified_position_no_thresh().unwrap().unwrap();
        assert_eq!(
            (
                state.mod_position,
                state.mod_code,
                state.mod_qual,
                state.modified,
                state.inferred,
            ),
            (0, ModCodeRepr::Code('q'), 220, true, false)
        );
        assert!(scanner.code_states.iter().all(|state| {
            state.mm_pos == initial_mm_pos - 2
                && state.mm_next == Some(0)
                && state.num_explicit_positions == 1
        }));
        assert_eq!(
            scanner
                .code_states
                .iter()
                .map(|state| state.ml_pos)
                .collect::<Vec<_>>(),
            (0..17).collect::<Vec<_>>()
        );

        let state =
            scanner.next_modified_position_no_thresh().unwrap().unwrap();
        assert_eq!(
            (
                state.mod_position,
                state.mod_code,
                state.mod_qual,
                state.modified,
                state.inferred,
            ),
            (1, ModCodeRepr::Code('a'), 200, true, false)
        );
        assert!(scanner.next_modified_position_no_thresh().unwrap().is_none());
    }

    #[test]
    fn test_one_code_groups_grow_beyond_inline_capacity() {
        for count in [16, 17] {
            // Separate MM groups for the same canonical base must target
            // different positions to avoid conflicting calls.
            let codes =
                (0..count).map(|idx| 10_000 + idx as u32).collect::<Vec<_>>();
            let mm = codes
                .iter()
                .enumerate()
                .map(|(idx, code)| format!("C+{code}?,{idx};"))
                .collect::<String>();
            let ml = vec![200; count];
            let seq = "C".repeat(count);
            let record = make_forward_record_with_mn(&mm, &ml, &seq);

            let mut scanner = BaseModsAdapter::<16>::new(&record).unwrap();
            assert_eq!(
                scanner
                    .code_states
                    .iter()
                    .map(|state| state.mod_code)
                    .collect::<Vec<_>>(),
                codes
                    .iter()
                    .copied()
                    .map(ModCodeRepr::ChEbi)
                    .collect::<Vec<_>>()
            );
            assert_eq!(
                scanner
                    .code_states
                    .iter()
                    .map(|state| state.mm_next)
                    .collect::<Vec<_>>(),
                (0..count).map(|idx| Some(idx as u32)).collect::<Vec<_>>()
            );
            assert_eq!(
                scanner
                    .code_states
                    .iter()
                    .map(|state| state.ml_pos)
                    .collect::<Vec<_>>(),
                (0..count).collect::<Vec<_>>()
            );
            assert!(scanner
                .code_states
                .iter()
                .all(|state| state.ml_stride == 1));
            assert_eq!(scanner.code_states.spilled(), count > 16);
            let expected = codes
                .iter()
                .enumerate()
                .map(|(position, code)| {
                    (position, ModCodeRepr::ChEbi(*code), 200, true, false)
                })
                .collect::<Vec<_>>();
            assert_eq!(collect_states(&mut scanner), expected);
        }
    }

    #[test]
    fn test_explicit_position_count_is_not_code_capacity() {
        let seq = "C".repeat(17);
        let mm = format!("C+m?,{};", vec!["0"; 17].join(","));
        let ml = (200..217).collect::<Vec<u8>>();
        let record = make_forward_record_with_mn(&mm, &ml, &seq);

        let mut scanner = BaseModsAdapter::<16>::new(&record).unwrap();
        assert_eq!(scanner.code_states.len(), 1);
        assert_eq!(scanner.code_states[0].num_explicit_positions, 17);
        assert!(!scanner.code_states.spilled());
        let expected = ml
            .iter()
            .enumerate()
            .map(|(position, qual)| {
                (position, METHYL_CYTOSINE, *qual, true, false)
            })
            .collect::<Vec<_>>();
        assert_eq!(collect_states(&mut scanner), expected);
    }

    #[test]
    fn test_chebi() {
        let mm = "C+m.,0,0;C+76792.,0,0;";
        let seq = "CGGATGGAGTC";
        let ml = vec![100, 25, 101, 24];
        let record = make_record(mm, &ml, seq, None, true);
        let scanner = BaseModsAdapter::<2>::new(&record).unwrap();
        assert_eq!(
            scanner
                .code_states
                .iter()
                .map(|state| state.mod_code)
                .collect::<Vec<_>>(),
            vec![METHYL_CYTOSINE, ModCodeRepr::ChEbi(76792)]
        )
    }

    // TODO: Add tests for:
    // - [x] Forard
    // - [x] Implicit
    // - [x] Implicit reverse
    // - [ ] Explicit
    // - [x] Empty
    // - [x] Multi-base
    // - [ ] Remaining edge cases from htslib
}
