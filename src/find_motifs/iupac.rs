use anyhow::bail;
use std::fmt::{Display, Formatter};

pub(crate) mod nt_bytes {
    pub const A: u8 = 65;
    pub const C: u8 = 67;
    pub const G: u8 = 71;
    pub const T: u8 = 84;
    pub const BASES: [u8; 4] = [A, C, G, T];
}

pub(super) mod iupac_offsets {
    const A: usize = 0;
    const C: usize = 1;
    const G: usize = 2;
    const T: usize = 3;

    pub const ALPHABET_SIZE: usize = 4;
    pub const A_IDX: &'static [usize] = &[A];
    pub const C_IDX: &'static [usize] = &[C];
    pub const G_IDX: &'static [usize] = &[G];
    pub const T_IDX: &'static [usize] = &[T];
    pub const R: &'static [usize] = &[A, G];
    pub const Y: &'static [usize] = &[C, T];
    pub const S: &'static [usize] = &[C, G];
    pub const W: &'static [usize] = &[A, T];
    pub const K: &'static [usize] = &[G, T];
    pub const M: &'static [usize] = &[A, C];
    pub const B: &'static [usize] = &[C, G, T];
    pub const D: &'static [usize] = &[A, G, T];
    pub const H: &'static [usize] = &[A, C, T];
    pub const V: &'static [usize] = &[A, C, G];
    pub const N_MOTIF: &'static [usize] = &[];
}

#[derive(Debug, Copy, Clone, Eq, PartialEq, Hash)]
pub(super) enum IupacBase {
    A,
    C,
    G,
    T,
    R,
    Y,
    S,
    W,
    K,
    M,
    B,
    D,
    H,
    V,
    N,
    // matches nothing
    Hole,
}

impl Display for IupacBase {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            IupacBase::A => write!(f, "A"),
            IupacBase::C => write!(f, "C"),
            IupacBase::G => write!(f, "G"),
            IupacBase::T => write!(f, "T"),
            IupacBase::R => write!(f, "R"),
            IupacBase::Y => write!(f, "Y"),
            IupacBase::S => write!(f, "S"),
            IupacBase::W => write!(f, "W"),
            IupacBase::K => write!(f, "K"),
            IupacBase::M => write!(f, "M"),
            IupacBase::B => write!(f, "B"),
            IupacBase::D => write!(f, "D"),
            IupacBase::H => write!(f, "H"),
            IupacBase::V => write!(f, "V"),
            IupacBase::N => write!(f, "N"),
            IupacBase::Hole => write!(f, "-"),
        }
    }
}

impl IupacBase {
    pub(super) fn parse_char(b: char) -> anyhow::Result<Self> {
        match b {
            'A' => Ok(Self::A),
            'C' => Ok(Self::C),
            'G' => Ok(Self::G),
            'T' => Ok(Self::T),
            'R' => Ok(Self::R),
            'Y' => Ok(Self::Y),
            'S' => Ok(Self::S),
            'W' => Ok(Self::W),
            'K' => Ok(Self::K),
            'M' => Ok(Self::M),
            'B' => Ok(Self::B),
            'D' => Ok(Self::D),
            'H' => Ok(Self::H),
            'V' => Ok(Self::V),
            'N' => Ok(Self::N),
            _ => bail!("unrecognized IUPAC code {b}"),
        }
    }

    pub(super) fn from_base(nt: u8) -> anyhow::Result<Self> {
        match nt {
            nt_bytes::A => Ok(Self::A),
            nt_bytes::C => Ok(Self::C),
            nt_bytes::G => Ok(Self::G),
            nt_bytes::T => Ok(Self::T),
            _ => bail!("should not be creating a non-DNA base"),
        }
    }

    pub(super) fn from_base_unchecked(nt: u8) -> Self {
        Self::from_base(nt).unwrap()
    }

    pub(super) fn intersect(self, other: Self) -> Self {
        match self {
            IupacBase::A => match other {
                IupacBase::A => self,
                _ => IupacBase::Hole,
            },
            IupacBase::C => match other {
                IupacBase::C => self,
                _ => IupacBase::Hole,
            },
            IupacBase::G => match other {
                IupacBase::G => self,
                _ => IupacBase::Hole,
            },
            IupacBase::T => match other {
                IupacBase::T => self,
                _ => IupacBase::Hole,
            },
            // {A, G}
            IupacBase::R => match other {
                IupacBase::A => Self::A,
                IupacBase::G => Self::G,
                IupacBase::R => self,
                IupacBase::S => Self::G,
                IupacBase::W => Self::A,
                IupacBase::K => Self::G,
                IupacBase::M => Self::A,
                IupacBase::B => Self::G,
                IupacBase::D => self,
                IupacBase::H => Self::A,
                IupacBase::V => self,
                IupacBase::N => self,
                _ => IupacBase::Hole,
            },
            // {C, T}
            IupacBase::Y => match other {
                IupacBase::C => Self::C,
                IupacBase::T => Self::T,
                IupacBase::Y => self,
                IupacBase::S => Self::C,
                IupacBase::W => Self::T,
                IupacBase::K => Self::T,
                IupacBase::M => Self::C,
                IupacBase::B => self,
                IupacBase::D => Self::T,
                IupacBase::H => self,
                IupacBase::V => Self::C,
                IupacBase::N => self,
                _ => IupacBase::Hole,
            },
            // {C, G}
            IupacBase::S => match other {
                IupacBase::C => Self::C,
                IupacBase::G => Self::G,
                IupacBase::R => Self::G,
                IupacBase::Y => Self::C,
                IupacBase::S => self,
                IupacBase::K => Self::G,
                IupacBase::M => Self::C,
                IupacBase::B => self,
                IupacBase::D => Self::G,
                IupacBase::H => Self::C,
                IupacBase::V => self,
                IupacBase::N => self,
                _ => IupacBase::Hole,
            },
            // {A, T}
            IupacBase::W => match other {
                IupacBase::A => Self::A,
                IupacBase::T => Self::T,
                IupacBase::R => Self::A,
                IupacBase::Y => Self::T,
                IupacBase::W => self,
                IupacBase::K => Self::T,
                IupacBase::M => Self::A,
                IupacBase::B => Self::T,
                IupacBase::D => Self::A,
                IupacBase::H => self,
                IupacBase::V => Self::A,
                IupacBase::N => self,
                _ => IupacBase::Hole,
            },
            // {G, T}
            IupacBase::K => match other {
                IupacBase::G => Self::G,
                IupacBase::T => Self::T,
                IupacBase::R => Self::G,
                IupacBase::Y => Self::T,
                IupacBase::S => Self::G,
                IupacBase::W => Self::T,
                IupacBase::K => self,
                IupacBase::B => self,
                IupacBase::D => self,
                IupacBase::H => Self::T,
                IupacBase::V => Self::G,
                IupacBase::N => self,
                _ => IupacBase::Hole,
            },
            // {A, C}
            IupacBase::M => match other {
                IupacBase::A => Self::A,
                IupacBase::C => Self::C,
                IupacBase::R => Self::A,
                IupacBase::Y => Self::C,
                IupacBase::S => Self::C,
                IupacBase::W => Self::A,
                IupacBase::M => Self::A,
                IupacBase::B => Self::C,
                IupacBase::D => Self::A,
                IupacBase::H => self,
                IupacBase::V => self,
                IupacBase::N => self,
                _ => IupacBase::Hole,
            },
            // {C, G, T}
            IupacBase::B => match other {
                IupacBase::A => IupacBase::Hole,
                IupacBase::C => Self::C,
                IupacBase::G => Self::G,
                IupacBase::T => Self::T,
                IupacBase::R => Self::G,
                IupacBase::Y => Self::Y,
                IupacBase::S => Self::S,
                IupacBase::W => Self::T,
                IupacBase::K => Self::K,
                IupacBase::M => Self::C,
                IupacBase::B => self,
                IupacBase::D => Self::K,
                IupacBase::H => Self::Y,
                IupacBase::V => Self::S,
                IupacBase::N => self,
                IupacBase::Hole => IupacBase::Hole,
            },
            // {A, G, T}
            IupacBase::D => match other {
                IupacBase::A => Self::A,
                IupacBase::C => IupacBase::Hole,
                IupacBase::G => Self::G,
                IupacBase::T => Self::T,
                IupacBase::R => Self::R,
                IupacBase::Y => Self::T,
                IupacBase::S => Self::G,
                IupacBase::W => Self::W,
                IupacBase::K => Self::K,
                IupacBase::M => Self::A,
                IupacBase::B => Self::K,
                IupacBase::D => self,
                IupacBase::H => Self::W,
                IupacBase::V => Self::R,
                IupacBase::N => self,
                IupacBase::Hole => IupacBase::Hole,
            },
            // {A, C, T}
            IupacBase::H => match other {
                IupacBase::A => Self::A,
                IupacBase::C => Self::C,
                IupacBase::G => IupacBase::Hole,
                IupacBase::T => Self::T,
                IupacBase::R => Self::A,
                IupacBase::Y => Self::Y,
                IupacBase::S => Self::C,
                IupacBase::W => Self::W,
                IupacBase::K => Self::T,
                IupacBase::M => Self::M,
                IupacBase::B => Self::Y,
                IupacBase::D => Self::W,
                IupacBase::H => self,
                IupacBase::V => Self::M,
                IupacBase::N => self,
                IupacBase::Hole => IupacBase::Hole,
            },
            // {A, C, G}
            IupacBase::V => match other {
                IupacBase::A => Self::A,
                IupacBase::C => Self::C,
                IupacBase::G => Self::G,
                IupacBase::T => IupacBase::Hole,
                IupacBase::R => Self::R,
                IupacBase::Y => Self::C,
                IupacBase::S => Self::S,
                IupacBase::W => Self::A,
                IupacBase::K => Self::G,
                IupacBase::M => Self::M,
                IupacBase::B => Self::S,
                IupacBase::D => Self::R,
                IupacBase::H => Self::M,
                IupacBase::V => self,
                IupacBase::N => self,
                IupacBase::Hole => IupacBase::Hole,
            },
            IupacBase::N => self,
            IupacBase::Hole => IupacBase::Hole,
        }
    }

    pub(super) fn matches(&self, nt: u8) -> bool {
        // https://www.bioinformatics.org/sms/iupac.html
        use nt_bytes as nts;
        match self {
            IupacBase::A => nt == nts::A,
            IupacBase::C => nt == nts::C,
            IupacBase::G => nt == nts::G,
            IupacBase::T => nt == nts::T,
            IupacBase::R => nt == nts::A || nt == nts::G,
            IupacBase::Y => nt == nts::C || nt == nts::T,
            IupacBase::S => nt == nts::C || nt == nts::G,
            IupacBase::W => nt == nts::A || nt == nts::T,
            IupacBase::K => nt == nts::G || nt == nts::T,
            IupacBase::M => nt == nts::A || nt == nts::C,
            // these could be just !=
            IupacBase::B => nt == nts::C || nt == nts::G || nt == nts::T,
            IupacBase::D => nt == nts::A || nt == nts::G || nt == nts::T,
            IupacBase::H => nt == nts::A || nt == nts::C || nt == nts::T,
            IupacBase::V => nt == nts::A || nt == nts::C || nt == nts::G,
            IupacBase::N => {
                nt == nts::A || nt == nts::C || nt == nts::G || nt == nts::T
            }
            IupacBase::Hole => false,
        }
    }

    pub(super) fn is_fixed(&self) -> bool {
        match self {
            IupacBase::A | IupacBase::C | IupacBase::G | IupacBase::T => true,
            _ => false,
        }
    }

    pub(super) fn complete_match(&self, nt: u8) -> bool {
        // maybe could ne neater
        if self.matches(nt) {
            match self {
                Self::A | Self::C | Self::G | Self::T => true,
                _ => false,
            }
        } else {
            false
        }
    }

    pub(super) fn is_n(&self) -> bool {
        match self {
            IupacBase::N => true,
            _ => false,
        }
    }

    pub(super) fn remove_to_n(self, nt: u8) -> Self {
        if self == IupacBase::Hole {
            todo!("remove from Hole?")
        } else {
            match self.remove(nt) {
                IupacBase::Hole => IupacBase::N,
                b @ _ => b,
            }
        }
    }

    pub(super) fn remove(self, nt: u8) -> Self {
        use nt_bytes as nts;
        if self.complete_match(nt) {
            Self::Hole
        } else if self.matches(nt) {
            match self {
                IupacBase::R => match nt {
                    nts::A => Self::G,
                    nts::G => Self::A,
                    _ => unreachable!(),
                },
                IupacBase::Y => match nt {
                    nts::C => Self::T,
                    nts::T => Self::C,
                    _ => unreachable!(),
                },
                IupacBase::S => match nt {
                    nts::C => Self::G,
                    nts::G => Self::C,
                    _ => unreachable!(),
                },
                IupacBase::W => match nt {
                    nts::A => Self::T,
                    nts::T => Self::A,
                    _ => unreachable!(),
                },
                IupacBase::K => match nt {
                    nts::G => Self::T,
                    nts::T => Self::G,
                    _ => unreachable!(),
                },
                IupacBase::M => match nt {
                    nts::A => Self::C,
                    nts::C => Self::A,
                    _ => unreachable!(),
                },
                IupacBase::B => match nt {
                    nts::C => Self::K,
                    nts::G => Self::Y,
                    nts::T => Self::S,
                    _ => unreachable!(),
                },
                IupacBase::D => match nt {
                    nts::A => Self::K,
                    nts::G => Self::W,
                    nts::T => Self::R,
                    _ => unreachable!(),
                },
                IupacBase::H => match nt {
                    nts::A => Self::Y,
                    nts::C => Self::W,
                    nts::T => Self::M,
                    _ => unreachable!(),
                },
                IupacBase::V => match nt {
                    nts::A => Self::S,
                    nts::C => Self::R,
                    nts::G => Self::M,
                    _ => unreachable!(),
                },
                IupacBase::N => match nt {
                    nts::A => Self::B,
                    nts::C => Self::D,
                    nts::G => Self::H,
                    nts::T => Self::V,
                    _ => unreachable!(),
                },
                _ => unreachable!(),
            }
        } else {
            // do nothing, maybe log here? should this be allowed?
            self
        }
    }

    pub(super) fn add_mut(&mut self, nt: u8) {
        let this = self.add(nt);
        *self = this
    }

    pub(super) fn add(self, nt: u8) -> Self {
        use nt_bytes as nts;
        match self {
            IupacBase::A => match nt {
                nts::A => self,
                nts::C => Self::M,
                nts::G => Self::R,
                nts::T => Self::W,
                _ => unreachable!(),
            },
            IupacBase::C => match nt {
                nts::A => Self::M,
                nts::C => self,
                nts::G => Self::S,
                nts::T => Self::Y,
                _ => unreachable!(),
            },
            IupacBase::G => match nt {
                nts::A => Self::R,
                nts::C => Self::S,
                nts::G => self,
                nts::T => Self::K,
                _ => unreachable!(),
            },
            IupacBase::T => match nt {
                nts::A => Self::W,
                nts::C => Self::Y,
                nts::G => Self::K,
                nts::T => self,
                _ => unreachable!(),
            },
            IupacBase::R => match nt {
                nts::A | nts::G => self,
                nts::C => Self::V,
                nts::T => Self::D,
                _ => unreachable!(),
            },
            IupacBase::Y => match nt {
                nts::C | nts::T => self,
                nts::A => Self::H,
                nts::G => Self::B,
                _ => unreachable!(),
            },
            IupacBase::S => match nt {
                nts::C | nts::G => self,
                nts::A => Self::V,
                nts::T => Self::B,
                _ => unreachable!(),
            },
            IupacBase::W => match nt {
                nts::A | nts::T => self,
                nts::G => Self::D,
                nts::C => Self::H,
                _ => unreachable!(),
            },
            IupacBase::K => match nt {
                nts::G | nts::T => self,
                nts::A => Self::D,
                nts::C => Self::B,
                _ => unreachable!(),
            },
            IupacBase::M => match nt {
                nts::A | nts::C => self,
                nts::G => Self::V,
                nts::T => Self::H,
                _ => unreachable!(),
            },
            IupacBase::B => match nt {
                nts::A => Self::N,
                _ => self,
            },
            IupacBase::D => match nt {
                nts::C => Self::N,
                _ => self,
            },
            IupacBase::H => match nt {
                nts::G => Self::N,
                _ => self,
            },
            IupacBase::V => match nt {
                nts::T => Self::N,
                _ => self,
            },
            IupacBase::N => self,
            IupacBase::Hole => match nt {
                nts::A => Self::A,
                nts::C => Self::C,
                nts::G => Self::G,
                nts::T => Self::T,
                _ => unreachable!(),
            },
        }
    }

    pub(super) fn to_char(&self) -> char {
        match self {
            IupacBase::A => 'A',
            IupacBase::C => 'C',
            IupacBase::G => 'G',
            IupacBase::T => 'T',
            IupacBase::R => 'R',
            IupacBase::Y => 'Y',
            IupacBase::S => 'S',
            IupacBase::W => 'W',
            IupacBase::K => 'K',
            IupacBase::M => 'M',
            IupacBase::B => 'B',
            IupacBase::D => 'D',
            IupacBase::H => 'H',
            IupacBase::V => 'V',
            IupacBase::N => 'N',
            IupacBase::Hole => '-',
        }
    }

    pub(super) fn union(self, other: Self) -> Self {
        match self {
            IupacBase::A => match other {
                IupacBase::A => self,
                IupacBase::C => Self::M,
                IupacBase::G => Self::R,
                IupacBase::T => Self::W,
                IupacBase::R => other,
                IupacBase::Y => Self::H,
                IupacBase::S => Self::V,
                IupacBase::W => other,
                IupacBase::K => Self::D,
                IupacBase::M => other,
                IupacBase::B => Self::N,
                IupacBase::D => other,
                IupacBase::H => other,
                IupacBase::V => other,
                IupacBase::N => other,
                IupacBase::Hole => self,
            },
            IupacBase::C => match other {
                IupacBase::A => Self::M,
                IupacBase::C => self,
                IupacBase::G => Self::S,
                IupacBase::T => Self::Y,
                IupacBase::R => Self::V,
                IupacBase::Y => other,
                IupacBase::S => other,
                IupacBase::W => Self::H,
                IupacBase::K => Self::B,
                IupacBase::M => other,
                IupacBase::B => other,
                IupacBase::D => Self::N,
                IupacBase::H => other,
                IupacBase::V => other,
                IupacBase::N => other,
                IupacBase::Hole => self,
            },
            IupacBase::G => match other {
                IupacBase::A => Self::R,
                IupacBase::C => Self::S,
                IupacBase::G => self,
                IupacBase::T => Self::K,
                IupacBase::R => other,
                IupacBase::Y => Self::B,
                IupacBase::S => other,
                IupacBase::W => Self::D,
                IupacBase::K => other,
                IupacBase::M => Self::V,
                IupacBase::B => other,
                IupacBase::D => other,
                IupacBase::H => Self::N,
                IupacBase::V => other,
                IupacBase::N => other,
                IupacBase::Hole => self,
            },
            IupacBase::T => match other {
                IupacBase::A => Self::W,
                IupacBase::C => Self::Y,
                IupacBase::G => Self::K,
                IupacBase::T => self,
                IupacBase::R => Self::D,
                IupacBase::Y => other,
                IupacBase::S => Self::B,
                IupacBase::W => other,
                IupacBase::K => other,
                IupacBase::M => Self::H,
                IupacBase::B => other,
                IupacBase::D => other,
                IupacBase::H => other,
                IupacBase::V => Self::N,
                IupacBase::N => Self::N,
                IupacBase::Hole => self,
            },
            IupacBase::R => match other {
                IupacBase::A => self,
                IupacBase::C => Self::V,
                IupacBase::G => self,
                IupacBase::T => Self::D,
                IupacBase::R => self,
                IupacBase::Y => Self::N,
                IupacBase::S => Self::V,
                IupacBase::W => Self::D,
                IupacBase::K => Self::D,
                IupacBase::M => Self::V,
                IupacBase::B => Self::N,
                IupacBase::D => other,
                IupacBase::H => Self::N,
                IupacBase::V => other,
                IupacBase::N => other,
                IupacBase::Hole => self,
            },
            IupacBase::Y => match other {
                IupacBase::A => Self::H,
                IupacBase::C => self,
                IupacBase::G => Self::B,
                IupacBase::T => self,
                IupacBase::R => Self::N,
                IupacBase::Y => self,
                IupacBase::S => Self::B,
                IupacBase::W => Self::H,
                IupacBase::K => Self::B,
                IupacBase::M => Self::H,
                IupacBase::B => other,
                IupacBase::D => Self::N,
                IupacBase::H => other,
                IupacBase::V => Self::N,
                IupacBase::N => other,
                IupacBase::Hole => self,
            },
            IupacBase::S => match other {
                IupacBase::A => Self::V,
                IupacBase::C => self,
                IupacBase::G => self,
                IupacBase::T => Self::B,
                IupacBase::R => Self::V,
                IupacBase::Y => Self::B,
                IupacBase::S => self,
                IupacBase::W => Self::N,
                IupacBase::K => Self::B,
                IupacBase::M => Self::V,
                IupacBase::B => other,
                IupacBase::D => Self::N,
                IupacBase::H => Self::N,
                IupacBase::V => other,
                IupacBase::N => other,
                IupacBase::Hole => self,
            },
            IupacBase::W => match other {
                IupacBase::A => self,
                IupacBase::C => Self::H,
                IupacBase::G => Self::D,
                IupacBase::T => self,
                IupacBase::R => Self::D,
                IupacBase::Y => Self::H,
                IupacBase::S => Self::N,
                IupacBase::W => self,
                IupacBase::K => Self::D,
                IupacBase::M => Self::H,
                IupacBase::B => Self::N,
                IupacBase::D => other,
                IupacBase::H => other,
                IupacBase::V => Self::N,
                IupacBase::N => other,
                IupacBase::Hole => self,
            },
            IupacBase::K => match other {
                IupacBase::A => Self::D,
                IupacBase::C => Self::B,
                IupacBase::G => self,
                IupacBase::T => self,
                IupacBase::R => Self::D,
                IupacBase::Y => Self::B,
                IupacBase::S => Self::B,
                IupacBase::W => Self::D,
                IupacBase::K => self,
                IupacBase::M => Self::N,
                IupacBase::B => other,
                IupacBase::D => other,
                IupacBase::H => Self::N,
                IupacBase::V => Self::N,
                IupacBase::N => other,
                IupacBase::Hole => self,
            },
            IupacBase::M => match other {
                IupacBase::A => self,
                IupacBase::C => self,
                IupacBase::G => Self::V,
                IupacBase::T => Self::H,
                IupacBase::R => Self::V,
                IupacBase::Y => Self::H,
                IupacBase::S => Self::V,
                IupacBase::W => Self::H,
                IupacBase::K => Self::N,
                IupacBase::M => self,
                IupacBase::B => Self::N,
                IupacBase::D => Self::N,
                IupacBase::H => other,
                IupacBase::V => other,
                IupacBase::N => other,
                IupacBase::Hole => self,
            },
            IupacBase::B => {
                if other.matches(nt_bytes::A) {
                    Self::N
                } else {
                    self
                }
            }
            IupacBase::D => {
                if other.matches(nt_bytes::C) {
                    Self::N
                } else {
                    self
                }
            }
            IupacBase::H => {
                if other.matches(nt_bytes::G) {
                    Self::N
                } else {
                    self
                }
            }
            IupacBase::V => {
                if other.matches(nt_bytes::T) {
                    Self::N
                } else {
                    self
                }
            }
            IupacBase::N => other,
            IupacBase::Hole => self,
        }
    }

    pub(super) fn is_superset(&self, other: &Self) -> bool {
        match self {
            IupacBase::A | IupacBase::C | IupacBase::G | IupacBase::T => {
                self == other || other == &IupacBase::Hole
            }
            IupacBase::R => match other {
                IupacBase::A
                | IupacBase::G
                | IupacBase::R
                | IupacBase::Hole => true,
                _ => false,
            },
            IupacBase::Y => match other {
                IupacBase::C
                | IupacBase::T
                | IupacBase::Y
                | IupacBase::Hole => true,
                _ => false,
            },
            IupacBase::S => match other {
                IupacBase::C
                | IupacBase::G
                | IupacBase::S
                | IupacBase::Hole => true,
                _ => false,
            },
            IupacBase::W => match other {
                IupacBase::A
                | IupacBase::T
                | IupacBase::W
                | IupacBase::Hole => true,
                _ => false,
            },
            IupacBase::K => match other {
                IupacBase::G
                | IupacBase::T
                | IupacBase::K
                | IupacBase::Hole => true,
                _ => false,
            },
            IupacBase::M => match other {
                IupacBase::A
                | IupacBase::C
                | IupacBase::M
                | IupacBase::Hole => true,
                _ => false,
            },
            IupacBase::B => match other {
                IupacBase::N => false,
                IupacBase::Hole => true,
                _ => !other.matches(nt_bytes::A),
            },
            IupacBase::D => match other {
                IupacBase::N => false,
                IupacBase::Hole => true,
                _ => !other.matches(nt_bytes::C),
            },
            IupacBase::H => match other {
                IupacBase::N => false,
                IupacBase::Hole => true,
                _ => !other.matches(nt_bytes::G),
            },
            IupacBase::V => match other {
                IupacBase::N => false,
                IupacBase::Hole => true,
                _ => !other.matches(nt_bytes::T),
            },
            IupacBase::N => true,
            IupacBase::Hole => false,
        }
    }

    pub(super) fn to_offsets(&self) -> &'static [usize] {
        use iupac_offsets::*;
        match self {
            Self::A => A_IDX,
            Self::C => C_IDX,
            Self::G => G_IDX,
            Self::T => T_IDX,
            Self::R => R,
            Self::Y => Y,
            Self::S => S,
            Self::W => W,
            Self::K => K,
            Self::M => M,
            Self::B => B,
            Self::D => D,
            Self::H => H,
            Self::V => V,
            Self::N | Self::Hole => N_MOTIF,
        }
    }
}

#[cfg(test)]
mod iupac_tests {
    use crate::find_motifs::iupac::IupacBase;

    #[test]
    fn test_superset() {
        assert!(IupacBase::B.is_superset(&IupacBase::S));
        assert!(IupacBase::B.is_superset(&IupacBase::Y));
        assert!(IupacBase::B.is_superset(&IupacBase::K));

        assert!(IupacBase::H.is_superset(&IupacBase::W));
        assert!(IupacBase::H.is_superset(&IupacBase::M));
        assert!(!IupacBase::H.is_superset(&IupacBase::K));
    }
}
