use log::debug;

#[allow(non_camel_case_types)]
#[derive(Debug, Copy, Clone, Eq, PartialEq, Hash)]
pub enum ModCode {
    /// canonical A
    A,
    /// canonical C
    C,
    a,
    h,
    m,
    G,
    T,
    /// Any C mod
    anyC,
    /// Any A mod
    anyA,
}

impl ModCode {
    pub(crate) fn parse_raw_mod_code(
        raw_mod_code: char,
    ) -> Result<Self, String> {
        match raw_mod_code {
            'a' => Ok(Self::a),
            'h' => Ok(Self::h),
            'm' => Ok(Self::m),
            'C' => Ok(Self::anyC),
            'A' => Ok(Self::anyA),
            _ => Err(format!("no mod code for {raw_mod_code}")),
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
            Self::anyA => 'A',
            Self::anyC => 'C',
        }
    }

    pub fn canonical_base(&self) -> DnaBase {
        match self {
            Self::A => DnaBase::A,
            Self::C => DnaBase::C,
            Self::a => DnaBase::A,
            Self::h => DnaBase::C,
            Self::m => DnaBase::C,
            Self::G => DnaBase::G,
            Self::T => DnaBase::T,
            Self::anyC => DnaBase::C,
            Self::anyA => DnaBase::A,
        }
    }

    pub(crate) fn is_canonical(&self) -> bool {
        match self {
            Self::A | Self::C | Self::G | Self::T => true,
            _ => false,
        }
    }

    pub fn get_any_mod_code(raw_dna_base: char) -> Option<ModCode> {
        match raw_dna_base {
            'C' => Some(ModCode::anyC),
            'A' => Some(ModCode::anyA),
            _ => {
                debug!("raw dna base {raw_dna_base} does not have a 'any mod' code");
                None
            }
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
}
