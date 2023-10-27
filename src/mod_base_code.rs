use anyhow::{anyhow, Result as AnyhowResult};
use std::cmp::Ordering;
use std::fmt::{Display, Formatter};

pub trait ParseChar {
    fn parse_char(c: char) -> AnyhowResult<Self>
    where
        Self: Sized;
    fn char(&self) -> char;
}

pub const HYDROXY_METHYL_CYTOSINE: ModCodeRepr = ModCodeRepr::Code('h');
pub const METHYL_CYTOSINE: ModCodeRepr = ModCodeRepr::Code('m');
pub const SIX_METHYL_ADENINE: ModCodeRepr = ModCodeRepr::Code('a');
pub const ANY_ADENINE: ModCodeRepr = ModCodeRepr::Code('A');
pub const ANY_CYTOSINE: ModCodeRepr = ModCodeRepr::Code('C');

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

    // pub(crate) fn into_base_state(self) -> AnyhowResult<ModCode> {
    //     match self {
    //         Self::A => Ok(ModCode::A),
    //         Self::C => Ok(ModCode::C),
    //         Self::G => {
    //             Err(anyhow!("no mod code for canonical base {}", self.char()))
    //         }
    //         Self::T => {
    //             Err(anyhow!("no mod code for canonical base {}", self.char()))
    //         }
    //     }
    // }
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
