use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap};
use std::fmt::{Display, Formatter};

use crate::errs::{MkError, MkResult};
use crate::motifs::iupac::nt_bytes;
use anyhow::anyhow;
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

pub const RNA_TWO_OME_ADENINE_NAME: &'static str = "2OmeA";
pub const RNA_TWO_OME_CYTOSINE_NAME: &'static str = "2OmeC";
pub const RNA_TWO_OME_GUANINE_NAME: &'static str = "2OmeG";
pub const RNA_TWO_OME_URACIL_NAME: &'static str = "2OmeU";
pub const RNA_INOSINE_NAME: &'static str = "inosine";
pub const RNA_METHYL_CYTOSINE_NAME: &'static str = "m5C";
pub const RNA_SIX_METHYL_ADENINE_NAME: &'static str = "m6A";
pub const RNA_PSEUDOURIDINE_NAME: &'static str = "pseU";

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
    pub static ref RNA_ONT_NAMES_TO_MODOMICS_NAMES: FxHashMap<&'static str, &'static str> = {
        let hm = hash_map! {
            RNA_TWO_OME_ADENINE_NAME => "Am",
            RNA_TWO_OME_CYTOSINE_NAME => "Cm",
            RNA_TWO_OME_GUANINE_NAME => "Gm",
            RNA_TWO_OME_URACIL_NAME => "Um",
            RNA_INOSINE_NAME => "I",
            RNA_METHYL_CYTOSINE_NAME => "m5C",
            RNA_SIX_METHYL_ADENINE_NAME => "m6A",
            RNA_PSEUDOURIDINE_NAME => "Y",
        };
        hm.into_iter().collect()
    };
}

lazy_static! {
    pub static ref RNA_LONG_NAME_TO_CODE: FxHashMap<&'static str, ModCodeRepr> = {
        let hm = hash_map! {
            RNA_TWO_OME_ADENINE_NAME => TWO_OME_ADENINE,
            RNA_TWO_OME_CYTOSINE_NAME => TWO_OME_CYTOSINE,
            RNA_TWO_OME_GUANINE_NAME => TWO_OME_GUANINE,
            RNA_TWO_OME_URACIL_NAME => TWO_OME_URACIL,
            RNA_INOSINE_NAME => INOSINE,
            RNA_METHYL_CYTOSINE_NAME => METHYL_CYTOSINE,
            RNA_SIX_METHYL_ADENINE_NAME => SIX_METHYL_ADENINE,
            RNA_PSEUDOURIDINE_NAME => PSEUDOURIDINE,
        };
        hm.into_iter().collect()
    };
}

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
}

impl PartialOrd for ModCodeRepr {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        match (self, other) {
            (Self::Code(x), Self::Code(y)) => x.partial_cmp(y),
            (Self::Code(_), Self::ChEbi(_)) => Some(Ordering::Greater),
            (Self::ChEbi(x), Self::ChEbi(y)) => x.partial_cmp(y),
            (Self::ChEbi(_), Self::Code(_)) => Some(Ordering::Less),
        }
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

#[derive(
    Debug, Copy, Clone, Eq, PartialEq, Hash, PartialOrd, Ord, ValueEnum,
)]
pub enum DnaBase {
    #[clap(name = "A")]
    A,
    #[clap(name = "C")]
    C,
    #[clap(name = "G")]
    G,
    #[clap(name = "T")]
    T,
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
    pub prob_counts: HashMap<BaseAndState, BTreeMap<u8, usize>>,
}
