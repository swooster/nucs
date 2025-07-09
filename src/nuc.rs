use std::cmp::Ordering;
use std::fmt::{Display, Formatter};
use std::hash::Hash;
use std::ops::{BitOr, BitOrAssign};
use std::str::FromStr;

use crate::error::ParseSymbolError;
use crate::translation::GeneticCode;
use crate::{AmbiAmino, Amino, Symbol};

/// Concrete nucleotide
///
/// [`Nuc`]s are ordered by the value of their ASCII representation.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum Nuc {
    /// [`A`](Self::A)denine
    A = 0b0001,
    /// [`C`](Self::C)ytosine
    C = 0b0010,
    /// [`G`](Self::G)uanine
    G = 0b0100,
    /// [`T`](Self::T)hymine
    T = 0b1000,
}

impl Nuc {
    /// All [`Nuc`]s sorted in ascending order
    pub const ALL: [Self; 4] = Self::lit(b"ACGT");

    /// Swap [`A`](Self::A) and [`T`](Self::T), as well as [`C`](Self::C) and [`G`](Self::G).
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::Nuc;
    ///
    /// assert_eq!(Nuc::A.complement(), Nuc::T);
    /// assert_eq!(Nuc::C.complement(), Nuc::G);
    /// assert_eq!(Nuc::G.complement(), Nuc::C);
    /// assert_eq!(Nuc::T.complement(), Nuc::A);
    /// ```
    #[must_use]
    pub const fn complement(self) -> Self {
        match self {
            Self::A => Self::T,
            Self::C => Self::G,
            Self::G => Self::C,
            Self::T => Self::A,
        }
    }

    /// Return uppercase string representation
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::Nuc;
    ///
    /// assert_eq!(Nuc::A.to_str(), "A");
    /// ```
    #[must_use]
    pub const fn to_str(self) -> &'static str {
        match self {
            Self::A => "A",
            Self::C => "C",
            Self::G => "G",
            Self::T => "T",
        }
    }

    /// Construct from ASCII representation
    ///
    /// # Errors
    ///
    /// Returns [`ParseSymbolError`] if the given byte isn't `A`, `C`, `G` or `T`
    /// (case-insensitive).
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::Nuc;
    ///
    /// assert_eq!(Nuc::from_ascii(b'A'), Ok(Nuc::A));
    /// assert!(Nuc::from_ascii(b'B').is_err());
    /// ```
    pub const fn from_ascii(ascii: u8) -> Result<Self, ParseSymbolError> {
        Ok(match ascii {
            b'a' | b'A' => Self::A,
            b'c' | b'C' => Self::C,
            b'g' | b'G' => Self::G,
            b't' | b'T' => Self::T,
            _ => return Err(ParseSymbolError::new::<Self>()),
        })
    }

    /// Return uppercase ASCII representation
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::Nuc;
    ///
    /// assert_eq!(Nuc::A.to_ascii(), b'A');
    /// ```
    #[must_use]
    pub const fn to_ascii(self) -> u8 {
        self.to_str().as_bytes()[0]
    }

    /// Construct [`Nuc`] array from literal without allocating.
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::Nuc;
    ///
    /// let dna1 = Nuc::lit(b"TACT");
    /// // ...is shorthand for...
    /// use Nuc::{T, A, C};
    /// let dna2 = [T, A, C, T];
    ///
    /// assert_eq!(dna1, dna2);
    /// ```
    ///
    /// # Panics
    ///
    /// This panics if the supplied literal isn't valid. **Whitespace is NOT allowed**
    /// because the returned array must have the same length.
    #[must_use]
    #[track_caller]
    pub const fn lit<const N: usize>(literal: &[u8; N]) -> [Nuc; N] {
        let mut nucs = [Self::A; N];
        let mut i = 0;
        while i < literal.len() {
            let Ok(nuc) = Self::from_ascii(literal[i]) else {
                panic!("Invalid Nuc in literal");
            };
            nucs[i] = nuc;
            i += 1;
        }
        nucs
    }
}

impl Display for Nuc {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        self.to_str().fmt(f)
    }
}

impl FromStr for Nuc {
    type Err = ParseSymbolError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.as_bytes() {
            [b] => Self::from_ascii(*b),
            _ => Err(ParseSymbolError::new::<Self>()),
        }
    }
}

impl TryFrom<AmbiNuc> for Nuc {
    type Error = AmbiNuc;

    fn try_from(nuc: AmbiNuc) -> Result<Self, Self::Error> {
        Ok(match nuc {
            AmbiNuc::A => Self::A,
            AmbiNuc::C => Self::C,
            AmbiNuc::G => Self::G,
            AmbiNuc::T => Self::T,
            other => return Err(other),
        })
    }
}

impl AsRef<Nuc> for Nuc {
    fn as_ref(&self) -> &Nuc {
        self
    }
}

impl AsMut<Nuc> for Nuc {
    fn as_mut(&mut self) -> &mut Nuc {
        self
    }
}

/// Ambiguous nucleotide
///
/// Note that [`AmbiNuc`] always represents *at least one* possible nucleotide;
/// there is no such thing as an empty/null [`AmbiNuc`].
///
/// For details, see: <https://en.wikipedia.org/wiki/FASTA_format#Sequence_representation>
///
/// <div class="warning">
///
/// Beware that [`AmbiNuc`]'s order isn't alphabetic.
///
/// </div>
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum AmbiNuc {
    /// [`A`](Self::A)denine
    A = Nuc::A as u8,
    /// [`C`](Self::C)ytosine
    C = Nuc::C as u8,
    /// [`G`](Self::G)uanine
    G = Nuc::G as u8,
    /// [`T`](Self::T)hymine
    T = Nuc::T as u8,

    /// [`A`](Self::A)/[`G`](Self::G): pu[`R`](Self::R)ine
    R = Nuc::A as u8 | Nuc::G as u8,
    /// [`C`](Self::C)/[`T`](Self::T): p[`Y`](Self::Y)rimidines
    Y = Nuc::C as u8 | Nuc::T as u8,

    /// [`A`](Self::A)/[`C`](Self::C): a[`M`](Self::M)ino groups
    M = Nuc::A as u8 | Nuc::C as u8,
    /// [`G`](Self::G)/[`T`](Self::T): [`K`](Self::K)etones
    K = Nuc::G as u8 | Nuc::T as u8,

    /// [`A`](Self::A)/[`T`](Self::T): [`W`](Self::W)eak interaction
    W = Nuc::A as u8 | Nuc::T as u8,
    /// [`C`](Self::C)/[`G`](Self::G): [`S`](Self::S)trong interaction
    S = Nuc::C as u8 | Nuc::G as u8,

    /// [`C`](Self::C)/[`G`](Self::G)/[`T`](Self::T): Not [`A`](Self::A) ([`B`](Self::B) comes after [`A`](Self::A))
    B = Nuc::C as u8 | Nuc::G as u8 | Nuc::T as u8,
    /// [`A`](Self::A)/[`G`](Self::G)/[`T`](Self::T): Not [`C`](Self::C) ([`D`](Self::D) comes after [`C`](Self::C))
    D = Nuc::A as u8 | Nuc::G as u8 | Nuc::T as u8,
    /// [`A`](Self::A)/[`C`](Self::C)/[`T`](Self::T): Not [`G`](Self::G) ([`H`](Self::H) comes after [`G`](Self::G))
    H = Nuc::A as u8 | Nuc::C as u8 | Nuc::T as u8,
    /// [`A`](Self::A)/[`C`](Self::C)/[`G`](Self::G): Not [`T`](Self::T) ([`V`](Self::V) comes after U/[`T`](Self::T))
    V = Nuc::A as u8 | Nuc::C as u8 | Nuc::G as u8,

    /// [`A`](Self::A)/[`C`](Self::C)/[`G`](Self::G)/[`T`](Self::T): [`N`](Self::N)ucleic acid
    N = Nuc::A as u8 | Nuc::C as u8 | Nuc::G as u8 | Nuc::T as u8,
}

impl AmbiNuc {
    /// All [`AmbiNuc`]s sorted in ascending order
    pub const ALL: [Self; 15] = Self::lit(b"ACMGRSVTWYHKDBN");

    /// Return ambiguous nucleotide containing complements for each potential value.
    ///
    /// ```
    /// use nucs::AmbiNuc;
    ///
    /// assert_eq!(AmbiNuc::A.complement(), AmbiNuc::T);
    /// assert_eq!(AmbiNuc::C.complement(), AmbiNuc::G);
    /// assert_eq!(AmbiNuc::G.complement(), AmbiNuc::C);
    /// assert_eq!(AmbiNuc::T.complement(), AmbiNuc::A);
    ///
    /// // D -> A|G|T -> T|C|A -> H
    /// assert_eq!(AmbiNuc::D.complement(), AmbiNuc::H);
    /// ```
    #[must_use]
    pub const fn complement(self) -> Self {
        match self {
            Self::A => Self::T,
            Self::C => Self::G,
            Self::G => Self::C,
            Self::T => Self::A,

            Self::R => Self::Y,
            Self::Y => Self::R,
            Self::M => Self::K,
            Self::K => Self::M,
            Self::W => Self::W,
            Self::S => Self::S,

            Self::B => Self::V,
            Self::D => Self::H,
            Self::H => Self::D,
            Self::V => Self::B,

            Self::N => Self::N,
        }
    }

    /// Return uppercase string representation
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::AmbiNuc;
    ///
    /// assert_eq!(AmbiNuc::A.to_str(), "A");
    /// assert_eq!(AmbiNuc::N.to_str(), "N");
    /// ```
    #[must_use]
    pub const fn to_str(self) -> &'static str {
        match self {
            Self::A => "A",
            Self::C => "C",
            Self::G => "G",
            Self::T => "T",

            Self::R => "R",
            Self::Y => "Y",
            Self::M => "M",
            Self::K => "K",
            Self::W => "W",
            Self::S => "S",

            Self::B => "B",
            Self::D => "D",
            Self::H => "H",
            Self::V => "V",

            Self::N => "N",
        }
    }

    /// Construct from ASCII representation
    ///
    /// # Errors
    ///
    /// Returns [`ParseSymbolError`] if the given byte isn't `A`, `B`, `C`, `D`, `G`, `H`,
    /// `K`, `M`, `N`, `R`, `S`, `T`, `V`, `W` or `Y` (case-insensitve).
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::AmbiNuc;
    ///
    /// assert_eq!(AmbiNuc::from_ascii(b'A'), Ok(AmbiNuc::A));
    /// assert!(AmbiNuc::from_ascii(b'E').is_err());
    /// ```
    pub const fn from_ascii(ascii: u8) -> Result<Self, ParseSymbolError> {
        Ok(match ascii {
            b'a' | b'A' => Self::A,
            b'c' | b'C' => Self::C,
            b'g' | b'G' => Self::G,
            b't' | b'T' => Self::T,

            b'r' | b'R' => Self::R,
            b'y' | b'Y' => Self::Y,
            b'm' | b'M' => Self::M,
            b'k' | b'K' => Self::K,
            b'w' | b'W' => Self::W,
            b's' | b'S' => Self::S,

            b'b' | b'B' => Self::B,
            b'd' | b'D' => Self::D,
            b'h' | b'H' => Self::H,
            b'v' | b'V' => Self::V,

            b'n' | b'N' => Self::N,

            _ => return Err(ParseSymbolError::new::<Self>()),
        })
    }

    /// Return uppercase ASCII representation
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::AmbiNuc;
    ///
    /// assert_eq!(AmbiNuc::A.to_ascii(), b'A');
    /// assert_eq!(AmbiNuc::N.to_ascii(), b'N');
    /// ```
    #[must_use]
    pub const fn to_ascii(self) -> u8 {
        self.to_str().as_bytes()[0]
    }

    /// Return iterator of [`Nuc`]s that this ambiguous nucleotide could be.
    ///
    /// The iterator is guaranteed to return things in sorted order and without duplicates,
    /// and its contents are guaranteed to recompose into this [`AmbiNuc`] via [`BitOr`].
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{AmbiNuc, Nuc};
    ///
    /// assert!(AmbiNuc::B.iter().eq(Nuc::lit(b"CGT")));
    /// ```
    pub fn iter(self) -> std::iter::Copied<std::slice::Iter<'static, Nuc>> {
        match self {
            Self::A => &[Nuc::A] as &[_],
            Self::C => &[Nuc::C],
            Self::G => &[Nuc::G],
            Self::T => &[Nuc::T],

            Self::R => &[Nuc::A, Nuc::G],
            Self::Y => &[Nuc::C, Nuc::T],
            Self::M => &[Nuc::A, Nuc::C],
            Self::K => &[Nuc::G, Nuc::T],
            Self::W => &[Nuc::A, Nuc::T],
            Self::S => &[Nuc::C, Nuc::G],

            Self::B => &[Nuc::C, Nuc::G, Nuc::T],
            Self::D => &[Nuc::A, Nuc::G, Nuc::T],
            Self::H => &[Nuc::A, Nuc::C, Nuc::T],
            Self::V => &[Nuc::A, Nuc::C, Nuc::G],

            Self::N => &[Nuc::A, Nuc::C, Nuc::G, Nuc::T],
        }
        .iter()
        .copied()
    }

    /// Construct [`AmbiNuc`] array from literal without allocating.
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::AmbiNuc;
    ///
    /// let dna1 = AmbiNuc::lit(b"NASTYGRAM");
    /// // ...is shorthand for...
    /// use AmbiNuc::{N, A, S, T, Y, G, R, M};
    /// let dna2 = [N, A, S, T, Y, G, R, A, M];
    ///
    /// assert_eq!(dna1, dna2);
    /// ```
    ///
    /// # Panics
    ///
    /// This panics if the supplied literal isn't valid. **Whitespace is NOT allowed**
    /// because the returned array must have the same length.
    #[must_use]
    #[track_caller]
    pub const fn lit<const N: usize>(literal: &[u8; N]) -> [AmbiNuc; N] {
        let mut nucs = [Self::A; N];
        let mut i = 0;
        while i < literal.len() {
            let Ok(nuc) = Self::from_ascii(literal[i]) else {
                panic!("Invalid AmbiNuc in literal");
            };
            nucs[i] = nuc;
            i += 1;
        }
        nucs
    }

    fn from_u8(byte: u8) -> Option<Self> {
        // (hopefully) efficiently match u8 against repr(u8) of enum
        // e.g. from_u8!(byte, MyEnum { A B C })
        macro_rules! from_u8 {
            ($byte:expr, $($variant:ident)+) => {{
                $(const $variant: u8 = AmbiNuc::$variant as u8;)+
                match $byte {
                    $($variant => Some(AmbiNuc::$variant),)+
                    _ => None,
                }
            }};
        }
        from_u8!(byte, A C M G R S V T W Y H K D B N)
    }
}

impl BitOr for AmbiNuc {
    type Output = AmbiNuc;

    fn bitor(self, rhs: AmbiNuc) -> Self::Output {
        // TODO: check performance of using .unwrap_unchecked() when unsafe is enabled?
        Self::from_u8(self as u8 | rhs as u8).expect("BUG: invalid nucleotide encountered")
    }
}

impl BitOr for &AmbiNuc {
    type Output = AmbiNuc;

    fn bitor(self, rhs: &AmbiNuc) -> Self::Output {
        *self | *rhs
    }
}

impl BitOr<Nuc> for AmbiNuc {
    type Output = AmbiNuc;

    fn bitor(self, rhs: Nuc) -> Self::Output {
        self | AmbiNuc::from(rhs)
    }
}

impl BitOr<&Nuc> for &AmbiNuc {
    type Output = AmbiNuc;

    fn bitor(self, rhs: &Nuc) -> Self::Output {
        *self | *rhs
    }
}

impl BitOr<AmbiNuc> for Nuc {
    type Output = AmbiNuc;

    fn bitor(self, rhs: AmbiNuc) -> Self::Output {
        AmbiNuc::from(self) | rhs
    }
}

impl BitOr<&AmbiNuc> for &Nuc {
    type Output = AmbiNuc;

    fn bitor(self, rhs: &AmbiNuc) -> Self::Output {
        *self | *rhs
    }
}

impl BitOr for Nuc {
    type Output = AmbiNuc;

    fn bitor(self, rhs: Nuc) -> Self::Output {
        AmbiNuc::from(self) | rhs
    }
}

impl BitOr for &Nuc {
    type Output = AmbiNuc;

    fn bitor(self, rhs: &Nuc) -> Self::Output {
        *self | *rhs
    }
}

impl BitOrAssign for AmbiNuc {
    fn bitor_assign(&mut self, rhs: AmbiNuc) {
        *self = *self | rhs;
    }
}

impl BitOrAssign<&AmbiNuc> for AmbiNuc {
    fn bitor_assign(&mut self, rhs: &AmbiNuc) {
        *self |= *rhs;
    }
}

impl BitOrAssign<Nuc> for AmbiNuc {
    fn bitor_assign(&mut self, rhs: Nuc) {
        *self = *self | rhs;
    }
}

impl BitOrAssign<&Nuc> for AmbiNuc {
    fn bitor_assign(&mut self, rhs: &Nuc) {
        *self |= *rhs;
    }
}

impl PartialEq<Nuc> for AmbiNuc {
    fn eq(&self, other: &Nuc) -> bool {
        *self == AmbiNuc::from(*other)
    }
}

impl PartialEq<AmbiNuc> for Nuc {
    fn eq(&self, other: &AmbiNuc) -> bool {
        other == self
    }
}

impl PartialOrd<Nuc> for AmbiNuc {
    fn partial_cmp(&self, other: &Nuc) -> Option<Ordering> {
        Some(self.cmp(&AmbiNuc::from(*other)))
    }
}

impl PartialOrd<AmbiNuc> for Nuc {
    fn partial_cmp(&self, other: &AmbiNuc) -> Option<Ordering> {
        other.partial_cmp(self).map(Ordering::reverse)
    }
}

impl Display for AmbiNuc {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        self.to_str().fmt(f)
    }
}

impl FromStr for AmbiNuc {
    type Err = ParseSymbolError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.as_bytes() {
            [b] => Self::from_ascii(*b),
            _ => Err(ParseSymbolError::new::<Self>()),
        }
    }
}

impl From<Nuc> for AmbiNuc {
    fn from(nuc: Nuc) -> Self {
        match nuc {
            Nuc::A => Self::A,
            Nuc::C => Self::C,
            Nuc::G => Self::G,
            Nuc::T => Self::T,
        }
    }
}

impl AsRef<AmbiNuc> for AmbiNuc {
    fn as_ref(&self) -> &AmbiNuc {
        self
    }
}

impl AsMut<AmbiNuc> for AmbiNuc {
    fn as_mut(&mut self) -> &mut AmbiNuc {
        self
    }
}

impl IntoIterator for AmbiNuc {
    type IntoIter = std::iter::Copied<std::slice::Iter<'static, Nuc>>;
    type Item = Nuc;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

impl IntoIterator for &AmbiNuc {
    type IntoIter = std::iter::Copied<std::slice::Iter<'static, Nuc>>;
    type Item = Nuc;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

/// A nucleotide; either [`Nuc`] or [`AmbiNuc`].
pub trait Nucleotide: Symbol {
    /// Amino acid; either [`Amino`] or [`AmbiAmino`]
    type Amino: Symbol;

    /// All possible nucleotides in ascending order
    const ALL: &[Self];

    /// Return this nucleotide's complement.
    ///
    /// See [`Nuc::complement`] or [`AmbiNuc::complement`] for details.
    #[must_use]
    fn complement(self) -> Self;

    /// Translate a codon into an amino acid
    fn translate<G: GeneticCode>(codon: [Self; 3], genetic_code: G) -> Self::Amino;
}

impl Nucleotide for Nuc {
    type Amino = Amino;

    const ALL: &[Self] = &Nuc::ALL;

    fn complement(self) -> Self {
        Nuc::complement(self)
    }

    fn translate<G: GeneticCode>(codon: [Self; 3], genetic_code: G) -> Self::Amino {
        genetic_code.map_codon(codon)
    }
}

impl Nucleotide for AmbiNuc {
    type Amino = AmbiAmino;

    const ALL: &[Self] = &AmbiNuc::ALL;

    fn complement(self) -> Self {
        AmbiNuc::complement(self)
    }

    fn translate<G: GeneticCode>(codon: [Self; 3], genetic_code: G) -> Self::Amino {
        genetic_code.map_ambi_codon(codon)
    }
}

impl Symbol for Nuc {
    type Concrete = Nuc;
    type Ambiguous = AmbiNuc;

    fn to_str(self) -> &'static str {
        Self::to_str(self)
    }

    fn from_ascii(ascii: u8) -> Result<Self, ParseSymbolError> {
        Self::from_ascii(ascii)
    }

    fn to_ascii(self) -> u8 {
        Self::to_ascii(self)
    }

    fn lit<const N: usize>(literal: &[u8; N]) -> [Self; N] {
        Self::lit(literal)
    }
}

impl Symbol for AmbiNuc {
    type Concrete = Nuc;
    type Ambiguous = AmbiNuc;

    fn to_str(self) -> &'static str {
        Self::to_str(self)
    }

    fn from_ascii(ascii: u8) -> Result<Self, ParseSymbolError> {
        Self::from_ascii(ascii)
    }

    fn to_ascii(self) -> u8 {
        Self::to_ascii(self)
    }

    fn lit<const N: usize>(literal: &[u8; N]) -> [Self; N] {
        Self::lit(literal)
    }
}

impl crate::symbol::sealed::Sealed for Nuc {
    const NAME: &str = "nucleotide";
    const EXPECTED: &str = "one of A/C/G/T (case-insensitive)";

    #[cfg(feature = "unsafe")]
    fn as_ambi_nucs(nucs: &[Self]) -> &[AmbiNuc] {
        crate::casts::nucs_as_ambi(nucs)
    }

    #[cfg(feature = "unsafe")]
    fn to_nucs(nucs: &[Self]) -> Option<&[Nuc]> {
        Some(nucs)
    }

    #[cfg(feature = "unsafe")]
    fn to_nucs_mut(nucs: &mut [Self]) -> Option<&mut [Nuc]> {
        Some(nucs)
    }
}

impl crate::symbol::sealed::Sealed for AmbiNuc {
    const NAME: &str = "ambiguous nucleotide";
    const EXPECTED: &str = "one of A/C/G/T or B/D/H/K/M/N/R/S/V/W/Y (case-insensitive)";

    #[cfg(feature = "unsafe")]
    fn as_ambi_nucs(nucs: &[Self]) -> &[AmbiNuc] {
        nucs
    }

    #[cfg(feature = "unsafe")]
    fn to_nucs(nucs: &[Self]) -> Option<&[Nuc]> {
        crate::casts::ambi_to_nucs(nucs)
    }

    #[cfg(feature = "unsafe")]
    fn to_nucs_mut(nucs: &mut [Self]) -> Option<&mut [Nuc]> {
        crate::casts::ambi_to_nucs_mut(nucs)
    }
}

#[cfg(test)]
mod tests {
    use crate::DnaSlice;

    use super::*;

    #[test]
    fn complementation_is_involution() {
        for nuc in Nuc::ALL {
            assert_eq!(nuc.complement().complement(), nuc);
        }

        for nuc in AmbiNuc::ALL {
            assert_eq!(nuc.complement().complement(), nuc);
        }
    }

    #[test]
    fn complementation_commutes_with_possibilities() {
        // Essentially checks that `AmbiNuc::complement` is correct by
        // defining it in terms of `Nuc::complement` and `AmbiNuc::possibilities`
        for nuc in AmbiNuc::ALL {
            let mut expected: Vec<_> = nuc.iter().collect();
            expected.complement();
            expected.sort();
            assert!(nuc.complement().iter().eq(expected));
        }
    }

    #[test]
    fn unambiguous_nuc_is_its_only_possiblity() {
        for nuc in Nuc::ALL {
            assert!(AmbiNuc::from(nuc).iter().eq([nuc]));
        }
    }

    #[test]
    fn unambiguous_pair_of_nucs_are_their_only_possiblities() {
        for nuc1 in Nuc::ALL {
            for nuc2 in Nuc::ALL {
                if nuc1 < nuc2 {
                    assert!((nuc1 | nuc2).iter().eq([nuc1, nuc2]));
                }
            }
        }
    }

    #[test]
    fn nuc_possiblities_are_sorted_and_unique() {
        for nuc in AmbiNuc::ALL {
            assert!(nuc.iter().zip(nuc.iter().skip(1)).all(|(a, b)| a < b));
        }
    }

    #[test]
    fn ambi_nuc_bitor_is_itempotent() {
        for nuc in AmbiNuc::ALL {
            assert_eq!(nuc | nuc, nuc);
        }
    }

    #[test]
    fn possibilities_can_be_recomposed_into_ambi_nucs() {
        for nuc in AmbiNuc::ALL {
            assert_eq!(
                nuc.iter().map(AmbiNuc::from).reduce(|a, b| a | b),
                Some(nuc)
            );
        }
    }

    #[test]
    fn ambi_nuc_bitor_is_consistent_with_possibilities() {
        for ambi1 in AmbiNuc::ALL {
            for ambi2 in AmbiNuc::ALL {
                let mut expected: Vec<_> = ambi1.iter().chain(ambi2).collect();
                expected.sort();
                expected.dedup();
                assert!((ambi1 | ambi2).iter().eq(expected));
            }
        }
    }

    #[test]
    fn ambi_nuc_bitor_with_n_is_n() {
        for nuc in AmbiNuc::ALL {
            assert_eq!(nuc | AmbiNuc::N, AmbiNuc::N);
            assert_eq!(AmbiNuc::N | nuc, AmbiNuc::N);
        }
    }

    #[test]
    fn str_roundtrips() {
        for nuc in Nuc::ALL {
            assert_eq!(Nuc::from_str(nuc.to_str()), Ok(nuc));
        }
        for nuc in AmbiNuc::ALL {
            assert_eq!(AmbiNuc::from_str(nuc.to_str()), Ok(nuc));
        }
    }

    #[test]
    fn ascii_roundtrips() {
        for nuc in Nuc::ALL {
            assert_eq!(Nuc::from_ascii(nuc.to_ascii()), Ok(nuc));
        }
        for nuc in AmbiNuc::ALL {
            assert_eq!(AmbiNuc::from_ascii(nuc.to_ascii()), Ok(nuc));
        }
    }

    #[test]
    fn all_is_sorted() {
        let mut sorted = Nuc::ALL;
        sorted.sort();
        assert_eq!(Nuc::ALL, sorted);
        let mut sorted = AmbiNuc::ALL;
        sorted.sort();
        assert_eq!(AmbiNuc::ALL, sorted);
    }
}
