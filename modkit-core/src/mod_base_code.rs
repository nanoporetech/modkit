use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap};
use std::fmt::{Display, Formatter};
use std::str::FromStr;

use crate::errs::{MkError, MkResult};
use crate::motifs::iupac::nt_bytes;
use anyhow::{anyhow, bail};
use clap::ValueEnum;
use common_macros::hash_map;
use derive_new::new;
use lazy_static::lazy_static;
use rustc_hash::FxHashMap;

pub trait ParseChar {
    fn parse_char(c: char) -> MkResult<Self>
    where
        Self: Sized;
    fn char(&self) -> char;
}

// Cytosine mods
pub const METHYL_CYTOSINE: ModCodeRepr = ModCodeRepr::Code('m');
pub const HYDROXY_METHYL_CYTOSINE: ModCodeRepr = ModCodeRepr::Code('h');
pub const FORMYL_CYTOSINE: ModCodeRepr = ModCodeRepr::Code('f');
pub const CARBOXY_CYTOSINE: ModCodeRepr = ModCodeRepr::Code('c');
pub const FOUR_METHYL_CYTOSINE: ModCodeRepr = ModCodeRepr::ChEbi(21839);
pub const TWO_OME_CYTOSINE: ModCodeRepr = ModCodeRepr::ChEbi(19228);
pub const ANY_CYTOSINE: ModCodeRepr = ModCodeRepr::Code('C');

// Adenine mods
pub const SIX_METHYL_ADENINE: ModCodeRepr = ModCodeRepr::Code('a');
pub const ANY_ADENINE: ModCodeRepr = ModCodeRepr::Code('A');
pub const INOSINE: ModCodeRepr = ModCodeRepr::ChEbi(17596);
pub const TWO_OME_ADENINE: ModCodeRepr = ModCodeRepr::ChEbi(69426);

// Thymine(/Uracil) mods
pub const HYDROXY_METHYL_URACIL: ModCodeRepr = ModCodeRepr::Code('g');
pub const FORMYL_URACIL: ModCodeRepr = ModCodeRepr::Code('e');
pub const CARBOXY_URACIL: ModCodeRepr = ModCodeRepr::Code('b');
pub const ANY_THYMINE: ModCodeRepr = ModCodeRepr::Code('T');
pub const PSEUDOURIDINE: ModCodeRepr = ModCodeRepr::ChEbi(17802);
pub const DEOXY_URACIL: ModCodeRepr = ModCodeRepr::ChEbi(16450);
pub const TWO_OME_URACIL: ModCodeRepr = ModCodeRepr::ChEbi(19227);

// Guanine mods
pub const OXO_GUANINE: ModCodeRepr = ModCodeRepr::Code('o');
pub const TWO_OME_GUANINE: ModCodeRepr = ModCodeRepr::ChEbi(19229);
pub const ANY_GUANINE: ModCodeRepr = ModCodeRepr::Code('G');

pub const ANY_MOD_CODES: [ModCodeRepr; 4] =
    [ANY_ADENINE, ANY_CYTOSINE, ANY_GUANINE, ANY_THYMINE];
pub const SUPPORTED_CODES: [ModCodeRepr; 21] = [
    METHYL_CYTOSINE,
    HYDROXY_METHYL_CYTOSINE,
    FORMYL_CYTOSINE,
    CARBOXY_CYTOSINE,
    FOUR_METHYL_CYTOSINE,
    ANY_CYTOSINE,
    TWO_OME_CYTOSINE,
    SIX_METHYL_ADENINE,
    ANY_ADENINE,
    INOSINE,
    TWO_OME_ADENINE,
    HYDROXY_METHYL_URACIL,
    FORMYL_URACIL,
    CARBOXY_URACIL,
    TWO_OME_URACIL,
    ANY_THYMINE,
    PSEUDOURIDINE,
    OXO_GUANINE,
    TWO_OME_GUANINE,
    ANY_GUANINE,
    DEOXY_URACIL,
];

lazy_static! {
    pub static ref MOD_CODE_TO_DNA_BASE: FxHashMap<ModCodeRepr, DnaBase> = {
        let hm = hash_map! {
            METHYL_CYTOSINE => DnaBase::C,
            HYDROXY_METHYL_CYTOSINE => DnaBase::C,
            FORMYL_CYTOSINE => DnaBase::C,
            CARBOXY_CYTOSINE => DnaBase::C,
            FOUR_METHYL_CYTOSINE => DnaBase::C,
            TWO_OME_CYTOSINE => DnaBase::C,
            ANY_CYTOSINE => DnaBase::C,
            SIX_METHYL_ADENINE => DnaBase::A,
            TWO_OME_ADENINE => DnaBase::A,
            ANY_ADENINE => DnaBase::A,
            INOSINE => DnaBase::A,
            HYDROXY_METHYL_URACIL => DnaBase::T,
            FORMYL_URACIL => DnaBase::T,
            CARBOXY_URACIL => DnaBase::T,
            PSEUDOURIDINE => DnaBase::T,
            TWO_OME_URACIL => DnaBase::T,
            DEOXY_URACIL => DnaBase::T,
            ANY_THYMINE => DnaBase::T,
            OXO_GUANINE => DnaBase::G,
            TWO_OME_GUANINE => DnaBase::G,
            ANY_GUANINE => DnaBase::G,
        };
        hm.into_iter().collect()
    };
}

lazy_static! {
    pub static ref MOD_COLORS: HashMap<ModCodeRepr, String> = hash_map! {
            METHYL_CYTOSINE => "#FF0000".to_string(),
            HYDROXY_METHYL_CYTOSINE => "#FF00FF".to_string(),
            SIX_METHYL_ADENINE => "#0084A9".to_string(),
            FOUR_METHYL_CYTOSINE => "#FFA100".to_string()
    };
    pub static ref DNA_BASE_COLORS: HashMap<DnaBase, String> = hash_map! {
            DnaBase::C => "#0000FF".to_string(),
            DnaBase::A => "#009600".to_string(),
    };
}

#[derive(Debug, Copy, Clone, Eq, PartialEq, Ord, Hash)]
pub enum ModCodeRepr {
    Code(char),
    ChEbi(u32),
}

impl ModCodeRepr {
    pub fn parse(raw: &str) -> anyhow::Result<Self> {
        if let Ok(code) = raw.parse::<char>() {
            Ok(Self::Code(code))
        } else {
            if let Ok(chebi) = raw.parse::<u32>() {
                Ok(Self::ChEbi(chebi))
            } else {
                Err(anyhow!("failed to parse mod code {raw}"))
            }
        }
    }

    pub fn check_base(&self, dna_base: DnaBase) -> bool {
        if let Some(self_base) = MOD_CODE_TO_DNA_BASE.get(self) {
            *self_base == dna_base
        } else {
            false
        }
    }

    pub fn is_any(&self) -> bool {
        ANY_MOD_CODES.contains(self)
    }

    pub(crate) fn any_mod_code(dna_base: &DnaBase) -> Self {
        Self::Code(dna_base.char())
    }

    pub(crate) fn matches_any(&self, dna_base: &DnaBase) -> bool {
        match self {
            ModCodeRepr::Code(c) => *c == dna_base.char(),
            ModCodeRepr::ChEbi(_) => false,
        }
    }
}

impl PartialOrd for ModCodeRepr {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Display for ModCodeRepr {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ChEbi(x) => write!(f, "{}", x),
            Self::Code(x) => write!(f, "{}", x),
        }
    }
}

impl From<char> for ModCodeRepr {
    fn from(value: char) -> Self {
        Self::Code(value)
    }
}

impl From<u32> for ModCodeRepr {
    fn from(value: u32) -> Self {
        Self::ChEbi(value)
    }
}

impl From<i32> for ModCodeRepr {
    #[inline]
    fn from(value: i32) -> Self {
        if value < 0 {
            let pos_value = value.abs() as u32;
            Self::ChEbi(pos_value)
        } else {
            let ch: Result<u8, _> = value.try_into();
            match ch {
                Ok(c) => Self::Code(c as char),
                Err(_) => Self::ChEbi(value as u32),
            }
        }
    }
}

#[derive(
    Debug, Copy, Clone, Eq, PartialEq, Hash, PartialOrd, Ord, ValueEnum,
)]
#[repr(usize)]
pub enum DnaBase {
    #[clap(name = "A")]
    A = 0,
    #[clap(name = "C")]
    C = 1,
    #[clap(name = "G")]
    G = 2,
    #[clap(name = "T")]
    T = 3,
}

impl DnaBase {
    pub fn parse(nt: char) -> MkResult<Self> {
        match nt {
            'A' => Ok(Self::A),
            'C' => Ok(Self::C),
            'G' => Ok(Self::G),
            'T' => Ok(Self::T),
            _ => Err(MkError::InvalidDnaBase),
        }
    }

    pub fn complement(self) -> Self {
        match self {
            Self::A => Self::T,
            Self::C => Self::G,
            Self::G => Self::C,
            Self::T => Self::A,
        }
    }

    pub fn char(&self) -> char {
        match self {
            Self::A => 'A',
            Self::C => 'C',
            Self::G => 'G',
            Self::T => 'T',
        }
    }

    pub(crate) fn as_byte(&self) -> u8 {
        self.char() as u8
    }
}

impl TryFrom<u8> for DnaBase {
    type Error = MkError;

    #[inline]
    fn try_from(value: u8) -> Result<Self, Self::Error> {
        match value {
            nt_bytes::A => Ok(Self::A),
            nt_bytes::C => Ok(Self::C),
            nt_bytes::G => Ok(Self::G),
            nt_bytes::T => Ok(Self::T),
            _ => Err(MkError::InvalidDnaBase),
        }
    }
}

impl TryFrom<i32> for DnaBase {
    type Error = MkError;

    #[inline]
    fn try_from(value: i32) -> Result<Self, Self::Error> {
        let value: u8 =
            value.try_into().map_err(|_| MkError::InvalidDnaBase)?;
        match value {
            nt_bytes::A => Ok(Self::A),
            nt_bytes::C => Ok(Self::C),
            nt_bytes::G => Ok(Self::G),
            nt_bytes::T => Ok(Self::T),
            _ => Err(MkError::InvalidDnaBase),
        }
    }
}

impl TryFrom<usize> for DnaBase {
    type Error = MkError;

    #[inline]
    fn try_from(value: usize) -> Result<Self, Self::Error> {
        match value {
            0usize => Ok(Self::A),
            1usize => Ok(Self::C),
            2usize => Ok(Self::G),
            3usize => Ok(Self::T),
            _ => Err(MkError::InvalidDnaBase),
        }
    }
}

impl TryFrom<&str> for DnaBase {
    type Error = MkError;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        match value.to_ascii_uppercase().as_str() {
            "A" => Ok(DnaBase::A),
            "C" => Ok(DnaBase::C),
            "G" => Ok(DnaBase::G),
            "T" => Ok(DnaBase::T),
            _ => Err(MkError::InvalidDnaBase),
        }
    }
}

impl ParseChar for DnaBase {
    fn parse_char(c: char) -> MkResult<Self> {
        DnaBase::parse(c)
    }
    fn char(&self) -> char {
        self.char()
    }
}

impl Display for DnaBase {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.char())
    }
}

// TODO this little enum is ripe for a refactor, try to make it just { DnaBase,
//  Modified(code) | Canonical }
#[derive(Debug, Copy, Clone, Eq, PartialEq, Hash, PartialOrd, Ord)]
pub enum BaseState {
    Canonical(DnaBase),
    Modified(ModCodeRepr),
}

pub type BaseAndState = (DnaBase, BaseState);

impl Display for BaseState {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Canonical(_dna_base) => write!(f, "-"),
            Self::Modified(mod_code) => write!(f, "{mod_code}"),
        }
    }
}

#[derive(new)]
pub struct ProbHistogram {
    pub prob_counts: HashMap<BaseAndState, BTreeMap<u8, u64>>,
}

const LN_FIVE_HYDROXY_METHYL_CYTOSINE: &str = "5hmC";
const LN_FIVE_METHYL_CYTOSINE: &str = "5mC";
const LN_FOUR_METHYL_CYTOSINE: &str = "4mC";
const LN_SIX_METHYL_ADENINE: &str = "6mA";
const LN_ANY_ADENINE: &str = "A";
const LN_ANY_CYTOSINE: &str = "C";
const LN_ANY_THYMINE: &str = "T";
const LN_ANY_GUANINE: &str = "G";

const LN_RNA_SIX_METHYL_ADENINE: &str = "m6A";
const LN_RNA_FIVE_METHYL_CYTOSINE: &str = "m5C";
const LN_PSEUDOURIDINE: &str = "pseU";
const LN_INOSINE: &str = "inosine";
const LN_TWO_OME_GUANINE: &str = "2OmeG";
const LN_TWO_OME_CYTOSINE: &str = "2OmeC";
const LN_TWO_OME_ADENINE: &str = "2OmeA";
const LN_TWO_OME_URACIL: &str = "2OmeU";

pub static LONG_NAME_TO_CODE: std::sync::LazyLock<
    HashMap<&'static str, ModCodeRepr>,
> = std::sync::LazyLock::new(|| {
    hash_map! {
        LN_FIVE_HYDROXY_METHYL_CYTOSINE => HYDROXY_METHYL_CYTOSINE,
        LN_FIVE_METHYL_CYTOSINE => METHYL_CYTOSINE,
        LN_FOUR_METHYL_CYTOSINE => FOUR_METHYL_CYTOSINE,
        LN_SIX_METHYL_ADENINE => SIX_METHYL_ADENINE,
        LN_RNA_SIX_METHYL_ADENINE => SIX_METHYL_ADENINE,
        LN_RNA_FIVE_METHYL_CYTOSINE => METHYL_CYTOSINE,
        LN_PSEUDOURIDINE => PSEUDOURIDINE,
        LN_INOSINE => INOSINE,
        LN_TWO_OME_GUANINE => TWO_OME_GUANINE,
        LN_TWO_OME_CYTOSINE => TWO_OME_CYTOSINE,
        LN_TWO_OME_ADENINE => TWO_OME_ADENINE,
        LN_TWO_OME_URACIL => TWO_OME_URACIL,
        LN_ANY_ADENINE => ANY_ADENINE,
        LN_ANY_CYTOSINE => ANY_CYTOSINE,
        LN_ANY_THYMINE => ANY_THYMINE,
        LN_ANY_GUANINE => ANY_GUANINE,
    }
});

lazy_static! {
    pub(crate) static ref RNA_CODES_TO_MODOMICS_NAMES: FxHashMap<ModifiedBasesOptions, &'static str> = {
        let hm = hash_map! {
            // LN_TWO_OME_ADENINE => "Am",
            ModifiedBasesOptions::new(TWO_OME_ADENINE, DnaBase::A) => "Am",
            // LN_TWO_OME_CYTOSINE => "Cm",
            ModifiedBasesOptions::new(TWO_OME_CYTOSINE, DnaBase::C) => "Cm",
            // LN_TWO_OME_GUANINE => "Gm",
            ModifiedBasesOptions::new(TWO_OME_GUANINE, DnaBase::G) => "Gm",
            // LN_TWO_OME_URACIL => "Um",
            ModifiedBasesOptions::new(TWO_OME_URACIL, DnaBase::T) => "Um",
            // LN_INOSINE => "I",
            ModifiedBasesOptions::new(INOSINE, DnaBase::A) => "I",
            // LN_RNA_FIVE_METHYL_CYTOSINE => "m5C",
            ModifiedBasesOptions::new(METHYL_CYTOSINE, DnaBase::C) => "m5C",
            // LN_RNA_SIX_METHYL_ADENINE => "m6A",
            ModifiedBasesOptions::new(SIX_METHYL_ADENINE, DnaBase::A) => "m6A",
            // LN_PSEUDOURIDINE => "Y",
            ModifiedBasesOptions::new(PSEUDOURIDINE, DnaBase::T) => "Y",
        };
        hm.into_iter().collect()
    };
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Hash, new)]
pub(crate) struct ModifiedBasesOptions {
    pub mod_code: ModCodeRepr,
    pub primary_base: DnaBase,
}

impl FromStr for ModifiedBasesOptions {
    type Err = anyhow::Error;

    fn from_str(raw: &str) -> Result<Self, Self::Err> {
        if raw.contains(":") {
            let parts = raw.split(":").collect::<Vec<&str>>();
            if parts.len() != 2 {
                bail!(
                    "invalid mod specification {raw}, should be \
                     <primary_base>:<mod_code>, e.g. C:m"
                )
            }
            let primary_base = parts[0]
                .parse::<char>()
                .map_err(|e| anyhow!("invalid DNA base {}, {e}", parts[0]))
                .and_then(|b| DnaBase::parse(b).map_err(|e| e.into()))?;
            let mod_code = ModCodeRepr::parse(parts[1])?;

            Ok(Self { mod_code, primary_base })
        } else {
            let (mod_code, primary_base) = LONG_NAME_TO_CODE
                .get(raw)
                .and_then(|mod_code| {
                    MOD_CODE_TO_DNA_BASE.get(mod_code).map(|b| (*mod_code, *b))
                })
                .ok_or(anyhow!("unknown long-name base modification {raw}"))?;
            Ok(Self { mod_code, primary_base })
        }
    }
}

#[cfg(test)]
mod tests {
    use std::cmp::Ordering;

    use super::ModCodeRepr;

    const ORDERED_REPRESENTATIVES: [ModCodeRepr; 6] = [
        ModCodeRepr::Code('A'),
        ModCodeRepr::Code('a'),
        ModCodeRepr::Code('m'),
        ModCodeRepr::ChEbi(0),
        ModCodeRepr::ChEbi(1),
        ModCodeRepr::ChEbi(u32::MAX),
    ];

    #[test]
    fn mod_code_ordering_obeys_total_order_laws() {
        for left in ORDERED_REPRESENTATIVES {
            for right in ORDERED_REPRESENTATIVES {
                let ordering = left.cmp(&right);
                assert_eq!(
                    ordering == Ordering::Equal,
                    left == right,
                    "comparison equality differs for {left} and {right}"
                );
                assert_eq!(
                    ordering,
                    right.cmp(&left).reverse(),
                    "comparison is not antisymmetric for {left} and {right}"
                );
                assert_eq!(
                    left.partial_cmp(&right),
                    Some(ordering),
                    "Ord and PartialOrd differ for {left} and {right}"
                );

                for last in ORDERED_REPRESENTATIVES {
                    let left_to_right = left.cmp(&right);
                    let right_to_last = right.cmp(&last);
                    if left_to_right != Ordering::Greater
                        && right_to_last != Ordering::Greater
                    {
                        assert_ne!(
                            left.cmp(&last),
                            Ordering::Greater,
                            "ordering is not transitive for {left}, {right}, \
                             and {last}"
                        );
                    }

                    let partial_left_to_right =
                        left.partial_cmp(&right).unwrap();
                    let partial_right_to_last =
                        right.partial_cmp(&last).unwrap();
                    if partial_left_to_right != Ordering::Greater
                        && partial_right_to_last != Ordering::Greater
                    {
                        assert_ne!(
                            left.partial_cmp(&last).unwrap(),
                            Ordering::Greater,
                            "partial ordering is not transitive for {left}, \
                             {right}, and {last}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn mod_code_sorted_output_places_letter_codes_before_chebi_ids() {
        let mut codes = vec![
            ModCodeRepr::ChEbi(2),
            ModCodeRepr::Code('m'),
            ModCodeRepr::ChEbi(1),
            ModCodeRepr::Code('a'),
        ];
        codes.sort();

        let rendered =
            codes.into_iter().map(|code| code.to_string()).collect::<Vec<_>>();
        assert_eq!(rendered.join(","), "a,m,1,2");
    }
}
