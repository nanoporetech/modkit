use anyhow::{anyhow, Result as AnyhowResult};
use std::cmp::Ordering;
use std::fmt::{Display, Formatter};

pub trait ParseChar {
    fn parse_char(c: char) -> AnyhowResult<Self>
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
pub const ANY_CYTOSINE: ModCodeRepr = ModCodeRepr::Code('C');

// Adenine mods
pub const SIX_METHYL_ADENINE: ModCodeRepr = ModCodeRepr::Code('a');
pub const ANY_ADENINE: ModCodeRepr = ModCodeRepr::Code('A');

// Thymine(/Uracil) mods
pub const HYDROXY_METHYL_URACIL: ModCodeRepr = ModCodeRepr::Code('g');
pub const FORMYL_URACIL: ModCodeRepr = ModCodeRepr::Code('e');
pub const CARBOXY_URACIL: ModCodeRepr = ModCodeRepr::Code('b');
pub const ANY_THYMINE: ModCodeRepr = ModCodeRepr::Code('T');

// Guanine mods
pub const OXO_GUANINE: ModCodeRepr = ModCodeRepr::Code('o');
pub const ANY_GUANINE: ModCodeRepr = ModCodeRepr::Code('G');

pub const ANY_MOD_CODES: [ModCodeRepr; 4] =
    [ANY_ADENINE, ANY_CYTOSINE, ANY_GUANINE, ANY_THYMINE];
pub const SUPPORTED_CODES: [ModCodeRepr; 14] = [
    METHYL_CYTOSINE,
    HYDROXY_METHYL_CYTOSINE,
    FOUR_METHYL_CYTOSINE,
    CARBOXY_CYTOSINE,
    FOUR_METHYL_CYTOSINE,
    ANY_CYTOSINE,
    SIX_METHYL_ADENINE,
    ANY_ADENINE,
    HYDROXY_METHYL_URACIL,
    FORMYL_URACIL,
    CARBOXY_URACIL,
    ANY_THYMINE,
    OXO_GUANINE,
    ANY_GUANINE,
];

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
        match self {
            &METHYL_CYTOSINE
            | &HYDROXY_METHYL_CYTOSINE
            | &FORMYL_CYTOSINE
            | &CARBOXY_CYTOSINE
            | &FOUR_METHYL_CYTOSINE
            | &ANY_CYTOSINE => dna_base == DnaBase::C,
            &SIX_METHYL_ADENINE | &ANY_ADENINE => dna_base == DnaBase::A,
            &HYDROXY_METHYL_URACIL
            | &FORMYL_URACIL
            | &CARBOXY_URACIL
            | &ANY_THYMINE => dna_base == DnaBase::T,
            &OXO_GUANINE | &ANY_GUANINE => dna_base == DnaBase::G,
            _ => false,
        }
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

impl ModCodeRepr {
    pub(crate) fn any_mod_code(dna_base: &DnaBase) -> Self {
        Self::Code(dna_base.char())
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

#[derive(Debug, Copy, Clone, Eq, PartialEq, Hash, PartialOrd, Ord)]
pub enum DnaBase {
    A,
    C,
    G,
    T,
}

impl DnaBase {
    pub(crate) fn parse(nt: char) -> AnyhowResult<Self> {
        match nt {
            'A' => Ok(Self::A),
            'C' => Ok(Self::C),
            'G' => Ok(Self::G),
            'T' => Ok(Self::T),
            _ => Err(anyhow!("unknown? {nt}".to_string())),
        }
    }

    pub(crate) fn complement(self) -> Self {
        match self {
            Self::A => Self::T,
            Self::C => Self::G,
            Self::G => Self::C,
            Self::T => Self::A,
        }
    }

    pub(crate) fn char(&self) -> char {
        match self {
            Self::A => 'A',
            Self::C => 'C',
            Self::G => 'G',
            Self::T => 'T',
        }
    }
}

impl ParseChar for DnaBase {
    fn parse_char(c: char) -> AnyhowResult<Self> {
        DnaBase::parse(c)
    }
    fn char(&self) -> char {
        self.char()
    }
}

#[derive(Debug, Copy, Clone, Eq, PartialEq, Hash, PartialOrd, Ord)]
pub enum BaseState {
    Canonical(DnaBase),
    Modified(ModCodeRepr),
}

impl Display for BaseState {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Canonical(dna_base) => write!(f, "{}", dna_base.char()),
            Self::Modified(mod_code) => write!(f, "{}", mod_code),
        }
    }
}
