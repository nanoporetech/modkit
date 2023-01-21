use std::collections::HashSet;

#[allow(non_camel_case_types)]
#[derive(Debug, Copy, Clone, Eq, PartialEq, Hash)]
pub enum ModCode {
    A,
    C,
    a,
    h,
    m,
    G,
    T,
}

impl ModCode {
    pub(crate) fn parse_raw_mod_code(
        raw_mod_code: char,
    ) -> Result<Self, String> {
        match raw_mod_code {
            'a' => Ok(Self::a),
            'h' => Ok(Self::h),
            'm' => Ok(Self::m),
            _ => Err("no mod code for {raw_mod_code}".to_string()),
        }
    }

    pub fn char(&self) -> char {
        match self {
            Self::A => 'A',
            Self::C => 'C',
            Self::a => 'a',
            Self::h => 'h',
            Self::m => 'm',
            Self::G => 'G',
            Self::T => 'T',
        }
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
    pub(crate) fn parse(nt: char) -> Result<Self, String> {
        match nt {
            'A' => Ok(Self::A),
            'C' => Ok(Self::C),
            'G' => Ok(Self::G),
            'T' => Ok(Self::T),
            _ => Err("unknown? {nt}".to_string()),
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

    pub(crate) fn canonical_mod_code(self) -> Result<ModCode, String> {
        match self {
            Self::A => Ok(ModCode::A),
            Self::C => Ok(ModCode::C),
            Self::G => {
                Err(format!("no mod code for canonical base {}", self.char()))
            }
            Self::T => {
                Err(format!("no mod code for canonical base {}", self.char()))
            }
        }
    }

    pub(crate) fn get_mod_codes(&self) -> HashSet<ModCode> {
        match self {
            Self::A => HashSet::from([ModCode::a]),
            Self::C => HashSet::from([ModCode::m, ModCode::h]),
            _ => HashSet::new(),
        }
    }
}
