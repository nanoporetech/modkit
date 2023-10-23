use anyhow::{anyhow, Result as AnyhowResult};
use std::cmp::Ordering;
use std::fmt::{Display, Formatter};

pub trait ParseChar {
    fn parse_char(c: char) -> AnyhowResult<Self>
    where
        Self: Sized;
    fn char(&self) -> char;
}

// #[allow(non_camel_case_types)]
// #[derive(Debug, Copy, Clone, Eq, PartialEq, Hash)]
// pub enum ModCode {
//     /// canonical A
//     A,
//     /// canonical C
//     C,
//     a,
//     h,
//     m,
//     G,
//     T,
//     /// Any C mod
//     anyC,
//     /// Any A mod
//     anyA,
// }

// // todo consider renaming..
// #[derive(Debug, Copy, Clone, Eq, PartialEq, Hash)]
// pub struct ModCode {
//     repr: ModCodeRepr,
//     primary_base: DnaBase,
//     // true when this probability is canonical probability
//     canonical: bool,
// }
//
// impl ModCode {
//     // pub(crate) fn parse_raw_mod_code(raw_mod_code: char) -> AnyhowResult<Self> {
//     //     match raw_mod_code {
//     //         'a' => Ok(Self::a),
//     //         'h' => Ok(Self::h),
//     //         'm' => Ok(Self::m),
//     //         'C' => Ok(Self::anyC),
//     //         'A' => Ok(Self::anyA),
//     //         _ => Err(anyhow!("no mod code for {raw_mod_code}")),
//     //     }
//     // }
//
//     // pub fn char(&self) -> char {
//     //     unimplemented!()
//     // }
//
//     pub(crate) fn is_canonical(&self) -> bool {
//         match self {
//             Self::Canonical => true,
//             Self::Modified => false,
//         }
//     }
//
//     pub(crate) fn primary_base(&self) -> DnaBase {
//         match self {
//             ModCode::Canonical { primary_base } | ModCode::Modified { repr: _, primary_base} => {
//                 *primary_base
//             }
//         }
//
//     }
//
//     pub(crate) fn repr(&self) -> ModCodeRepr {
//
//
//         unimplemented!()
//     }
// }

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
        unimplemented!()
    }
}

impl PartialOrd for ModCodeRepr {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        todo!()
    }
}

impl ModCodeRepr {
    pub(crate) fn any_mod_code(dna_base: &DnaBase) -> Self {
        unimplemented!()
    }
}

impl Display for ModCodeRepr {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        todo!()
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
        todo!()
    }
}
