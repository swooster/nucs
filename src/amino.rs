use std::cmp::Ordering;
use std::fmt::{Display, Formatter};
use std::hash::Hash;
use std::num::NonZeroU32;
use std::ops::{BitOr, BitOrAssign};
use std::str::FromStr;

use crate::Symbol;
use crate::error::ParseSymbolError;
use crate::iter::AmbiAminoIter;

/// Amino acid
///
/// For details, see: <https://en.wikipedia.org/wiki/FASTA_format#Sequence_representation>
///
/// [`Amino`]s are ordered by the value of their ASCII representation.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum Amino {
    /// Translation stop
    Stop = b'*',
    /// Alanine
    A = b'A',
    /// Cysteine
    C = b'C',
    /// Aspartic acid
    D = b'D',
    /// Glutamic acid
    E = b'E',
    /// Phenylalanine
    F = b'F',
    /// Glycine
    G = b'G',
    /// Histidine
    H = b'H',
    /// Isoleucine
    I = b'I',
    /// Lysine
    K = b'K',
    /// Leucine
    L = b'L',
    /// Methionine
    M = b'M',
    /// Asparagine
    N = b'N',
    /// Pyrrolysine
    O = b'O',
    /// Proline,
    P = b'P',
    /// Glutamine
    Q = b'Q',
    /// Arginine
    R = b'R',
    /// Serine
    S = b'S',
    /// Threonine
    T = b'T',
    /// Selenocysteine
    U = b'U',
    /// Valine
    V = b'V',
    /// Tryptophan
    W = b'W',
    /// Tyrosine
    Y = b'Y',
}

impl Amino {
    /// All [`Amino`]s sorted in ascending order
    pub const ALL: [Self; 23] = Self::lit(b"*ACDEFGHIKLMNOPQRSTUVWY");

    /// Return uppercase string representation
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::Amino;
    ///
    /// assert_eq!(Amino::A.to_str(), "A");
    /// assert_eq!(Amino::Stop.to_str(), "*");
    /// ```
    #[must_use]
    pub fn to_str(self) -> &'static str {
        match self {
            Self::Stop => "*",

            Self::A => "A",
            Self::C => "C",
            Self::D => "D",
            Self::E => "E",
            Self::F => "F",
            Self::G => "G",
            Self::H => "H",
            Self::I => "I",
            Self::K => "K",
            Self::L => "L",
            Self::M => "M",
            Self::N => "N",
            Self::O => "O",
            Self::P => "P",
            Self::Q => "Q",
            Self::R => "R",
            Self::S => "S",
            Self::T => "T",
            Self::U => "U",
            Self::V => "V",
            Self::W => "W",
            Self::Y => "Y",
        }
    }

    /// Construct from ASCII representation
    ///
    /// # Errors
    ///
    /// Returns [`ParseSymbolError`] if the given byte isn't `A`, `C`, `D`, `E`, `F`,
    /// `G`, `H`, `I`, `K`, `L`, `M`, `N`, `O`, `P`, `Q`, `R`, `S`, `T`, `U`, `V`, `W`, `Y`
    /// or `*` (case-insensitive).
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::Amino;
    ///
    /// assert_eq!(Amino::from_ascii(b'A'), Ok(Amino::A));
    /// assert_eq!(Amino::from_ascii(b'*'), Ok(Amino::Stop));
    /// assert!(Amino::from_ascii(b'B').is_err());
    /// ```
    pub const fn from_ascii(ascii: u8) -> Result<Self, ParseSymbolError> {
        Ok(match ascii {
            b'*' => Self::Stop,

            b'A' | b'a' => Self::A,
            b'C' | b'c' => Self::C,
            b'D' | b'd' => Self::D,
            b'E' | b'e' => Self::E,
            b'F' | b'f' => Self::F,
            b'G' | b'g' => Self::G,
            b'H' | b'h' => Self::H,
            b'I' | b'i' => Self::I,
            b'K' | b'k' => Self::K,
            b'L' | b'l' => Self::L,
            b'M' | b'm' => Self::M,
            b'N' | b'n' => Self::N,
            b'O' | b'o' => Self::O,
            b'P' | b'p' => Self::P,
            b'Q' | b'q' => Self::Q,
            b'R' | b'r' => Self::R,
            b'S' | b's' => Self::S,
            b'T' | b't' => Self::T,
            b'U' | b'u' => Self::U,
            b'V' | b'v' => Self::V,
            b'W' | b'w' => Self::W,
            b'Y' | b'y' => Self::Y,

            _ => return Err(ParseSymbolError::new::<Self>()),
        })
    }

    /// Return uppercase ASCII representation
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::Amino;
    ///
    /// assert_eq!(Amino::A.to_ascii(), b'A');
    /// assert_eq!(Amino::Stop.to_ascii(), b'*');
    /// ```
    #[must_use]
    pub fn to_ascii(self) -> u8 {
        self.to_str().as_bytes()[0]
    }

    /// Construct [`Amino`] array from literal without allocating.
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::Amino;
    ///
    /// let aas1 = Amino::lit(b"ACME*WIDGETS");
    /// // ...is shorthand for...
    /// use Amino::{A, C, M, E, Stop, W, I, D, G, T, S};
    /// let aas2 = [A, C, M, E, Stop, W, I, D, G, E, T, S];
    ///
    /// assert_eq!(aas1, aas2);
    /// ```
    ///
    /// # Panics
    ///
    /// This panics if the supplied literal isn't valid. **Whitespace is NOT allowed**
    /// because the returned array must have the same length.
    #[must_use]
    #[track_caller]
    pub const fn lit<const N: usize>(literal: &[u8; N]) -> [Amino; N] {
        let mut aas = [Self::A; N];
        let mut i = 0;
        while i < literal.len() {
            let Ok(aa) = Self::from_ascii(literal[i]) else {
                panic!("Invalid Amino in literal");
            };
            aas[i] = aa;
            i += 1;
        }
        aas
    }
}

impl Display for Amino {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        self.to_str().fmt(f)
    }
}

impl FromStr for Amino {
    type Err = ParseSymbolError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.as_bytes() {
            [b] => Self::from_ascii(*b),
            _ => Err(ParseSymbolError::new::<Self>()),
        }
    }
}

impl AsRef<Amino> for Amino {
    fn as_ref(&self) -> &Amino {
        self
    }
}

impl AsMut<Amino> for Amino {
    fn as_mut(&mut self) -> &mut Amino {
        self
    }
}

/// Ambiguous amino acid
///
/// Note that [`AmbiAmino`] always represents *at least one* possible amino acid;
/// there is no such thing as an empty/null [`AmbiAmino`]. Also beware that
/// `ambi_amino.to_string().parse()` is **lossy** because there are too many [`AmbiAmino`]s
/// to represent in a single character; anything without a standard character will be
/// displayed as `"X"` but `"X"` is parsed as a maximally ambiguous [`AmbiAmino`].
/// However, the [`Debug`](std::fmt::Debug) implementation will exactly represent
/// the [`AmbiAmino`] at the cost of being 1-25 characters long.
///
/// For details, see: <https://en.wikipedia.org/wiki/FASTA_format#Sequence_representation>
///
/// <div class="warning">
///
/// Beware that [`AmbiAmino`]'s order isn't alphabetic.
///
/// </div>
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AmbiAmino(NonZeroU32);

// Implementation notes:
// Bit 0 is always on, allowing niche optimizations.
// Other bits are determined by Self::bit_offset.
// (the low 5 bits of the IUPAC ASCII representation are the bit offset)
// That unfortunately means * falls on bit 10, so the bits aren't in quite in alphabetical order.
// (we fix that up if anybody actually asks us to iterate)

macro_rules! ambi_amino_consts {
    ( $( #[doc = $doc:literal] $ambi_amino:ident = $($amino:ident)|+  ),+ $(,)? ) => {
        $(
            #[doc = $doc]
            pub const $ambi_amino: Self = const {
                let bits = 1 | $(Self::bit_mask(Amino::$amino))|+;
                Self(NonZeroU32::new(bits).expect("BUG: somehow 1 | X == 0"))
            };
        )+
    }
}

impl AmbiAmino {
    /// Translation stop (alias for consistency with [`Amino::Stop`])
    #[allow(non_upper_case_globals)]
    pub const Stop: Self = Self::STOP;

    ambi_amino_consts!(
        /// Translation stop
        STOP = Stop,
        /// Alanine
        A = A,
        /// Aspartic acid ([D](Self::D)) or asparagine ([N](Self::N))
        B = D | N,
        /// Cysteine
        C = C,
        /// Aspartic acid
        D = D,
        /// Glutamic acid
        E = E,
        /// Phenylalanine
        F = F,
        /// Glycine
        G = G,
        /// Histidine
        H = H,
        /// Isoleucine
        I = I,
        /// Leucine ([L](Self::L)) or isoleucine ([I](Self::I))
        J = L | I,
        /// Lysine
        K = K,
        /// Leucine
        L = L,
        /// Methionine
        M = M,
        /// Asparagine
        N = N,
        /// Pyrrolysine
        O = O,
        /// Proline,
        P = P,
        /// Glutamine
        Q = Q,
        /// Arginine
        R = R,
        /// Serine
        S = S,
        /// Threonine
        T = T,
        /// Selenocysteine
        U = U,
        /// Valine
        V = V,
        /// Tryptophan
        W = W,
        /// Tyrosine
        Y = Y,
        /// Glutamic acid ([E](Self::E)) or glutamine ([Q](Self::Q))
        Z = E | Q,
    );

    /// Any amino acid
    ///
    /// Note that this doesn't perfectly round-trip. Any ambiguous amino that can't be represented
    /// with a single-letter code will be displayed as "X" but "X" will be parsed as the
    /// [`X`](Self::X) constant, which refers to a maximally ambiguous [`AmbiAmino`].
    pub const X: Self = const {
        let mut i = 0;
        let mut bits = 1;
        while i < Amino::ALL.len() {
            bits |= Self::bit_mask(Amino::ALL[i]);
            i += 1;
        }
        Self(NonZeroU32::new(bits).expect("BUG: somehow 1 | X == 0"))
    };

    /// Cast [`AmbiAmino`] to [`NonZeroU32`].
    ///
    /// <div class="warning">
    ///
    /// No particular representation is promised and the exact values are not part of
    /// stability guarantees.
    ///
    /// </div>
    #[must_use]
    pub fn to_bits(self) -> NonZeroU32 {
        self.0
    }

    /// Return uppercase string representation
    ///
    /// Note that anything without a single-character representation will result in `"X"`.
    /// Use [`Debug`](std::fmt::Debug) formatting if you want a lossless representation.
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::AmbiAmino;
    ///
    /// assert_eq!(AmbiAmino::A.to_str(), "A");
    /// assert_eq!(AmbiAmino::STOP.to_str(), "*");
    /// assert_eq!((AmbiAmino::D | AmbiAmino::N).to_str(), "B");
    /// assert_eq!((AmbiAmino::D | AmbiAmino::E).to_str(), "X");
    /// ```
    #[must_use]
    pub fn to_str(self) -> &'static str {
        match self {
            Self::STOP => "*",
            Self::A => "A",
            Self::B => "B",
            Self::C => "C",
            Self::D => "D",
            Self::E => "E",
            Self::F => "F",
            Self::G => "G",
            Self::H => "H",
            Self::I => "I",
            Self::J => "J",
            Self::K => "K",
            Self::L => "L",
            Self::M => "M",
            Self::N => "N",
            Self::O => "O",
            Self::P => "P",
            Self::Q => "Q",
            Self::R => "R",
            Self::S => "S",
            Self::T => "T",
            Self::U => "U",
            Self::V => "V",
            Self::W => "W",
            Self::Y => "Y",
            Self::Z => "Z",
            _ => "X",
        }
    }

    /// Construct from ASCII representation
    ///
    /// # Errors
    ///
    /// Returns [`ParseSymbolError`] if the given byte isn't a letter or `*`.
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::AmbiAmino;
    ///
    /// assert_eq!(AmbiAmino::from_ascii(b'A'), Ok(AmbiAmino::A));
    /// assert_eq!(AmbiAmino::from_ascii(b'*'), Ok(AmbiAmino::STOP));
    /// assert_eq!(AmbiAmino::from_ascii(b'X'), Ok(AmbiAmino::X));
    /// assert!(AmbiAmino::from_ascii(b'.').is_err());
    /// ```
    pub const fn from_ascii(ascii: u8) -> Result<Self, ParseSymbolError> {
        Ok(match ascii {
            b'*' => Self::STOP,

            b'A' | b'a' => Self::A,
            b'B' | b'b' => Self::B,
            b'C' | b'c' => Self::C,
            b'D' | b'd' => Self::D,
            b'E' | b'e' => Self::E,
            b'F' | b'f' => Self::F,
            b'G' | b'g' => Self::G,
            b'H' | b'h' => Self::H,
            b'I' | b'i' => Self::I,
            b'J' | b'j' => Self::J,
            b'K' | b'k' => Self::K,
            b'L' | b'l' => Self::L,
            b'M' | b'm' => Self::M,
            b'N' | b'n' => Self::N,
            b'O' | b'o' => Self::O,
            b'P' | b'p' => Self::P,
            b'Q' | b'q' => Self::Q,
            b'R' | b'r' => Self::R,
            b'S' | b's' => Self::S,
            b'T' | b't' => Self::T,
            b'U' | b'u' => Self::U,
            b'V' | b'v' => Self::V,
            b'W' | b'w' => Self::W,
            b'X' | b'x' => Self::X,
            b'Y' | b'y' => Self::Y,
            b'Z' | b'z' => Self::Z,

            _ => return Err(ParseSymbolError::new::<Self>()),
        })
    }

    /// Return uppercase ASCII representation
    ///
    /// See [`to_str`](Self::to_str) for caveats.
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::AmbiAmino;
    ///
    /// assert_eq!(AmbiAmino::A.to_ascii(), b'A');
    /// assert_eq!(AmbiAmino::STOP.to_ascii(), b'*');
    /// assert_eq!((AmbiAmino::D | AmbiAmino::N).to_ascii(), b'B');
    /// assert_eq!((AmbiAmino::D | AmbiAmino::E).to_ascii(), b'X');
    /// ```
    #[must_use]
    pub fn to_ascii(self) -> u8 {
        self.to_str().as_bytes()[0]
    }

    /// Construct [`AmbiAmino`] array from literal without allocating.
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::AmbiAmino;
    ///
    /// let aas1 = AmbiAmino::lit(b"TOO*LONG");
    /// // ...is shorthand for...
    /// let aas2 = [
    ///     AmbiAmino::T, AmbiAmino::O, AmbiAmino::O, AmbiAmino::STOP,
    ///     AmbiAmino::L, AmbiAmino::O, AmbiAmino::N, AmbiAmino::G
    /// ];
    ///
    /// assert_eq!(aas1, aas2);
    /// ```
    ///
    /// # Panics
    ///
    /// This panics if the supplied literal isn't valid. **Whitespace is NOT allowed**
    /// because the returned array must have the same length.
    #[must_use]
    #[track_caller]
    pub const fn lit<const N: usize>(literal: &[u8; N]) -> [AmbiAmino; N] {
        let mut aas = [Self::A; N];
        let mut i = 0;
        while i < literal.len() {
            let Ok(aa) = Self::from_ascii(literal[i]) else {
                panic!("Invalid Amino in literal");
            };
            aas[i] = aa;
            i += 1;
        }
        aas
    }

    /// Return iterator of [`Amino`]s that this ambiguous amino acid could be.
    ///
    /// The iterator is guaranteed to return things in sorted order and without duplicates,
    /// and its contents are guaranteed to recompose into this [`AmbiAmino`] via [`BitOr`].
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{AmbiAmino, Amino};
    ///
    /// let aa = AmbiAmino::A | AmbiAmino::B | AmbiAmino::C;
    /// assert!(aa.iter().eq(Amino::lit(b"ACDN")));
    /// ```
    pub fn iter(self) -> AmbiAminoIter {
        AmbiAminoIter::new(self)
    }

    pub(crate) const fn bit_offset(amino: Amino) -> u8 {
        (amino as u8) % 32
    }

    pub(crate) const fn bit_mask(amino: Amino) -> u32 {
        1 << Self::bit_offset(amino)
    }
}

impl BitOr for AmbiAmino {
    type Output = AmbiAmino;

    fn bitor(self, rhs: AmbiAmino) -> Self::Output {
        Self(self.0 | rhs.0)
    }
}

impl BitOr for &AmbiAmino {
    type Output = AmbiAmino;

    fn bitor(self, rhs: &AmbiAmino) -> Self::Output {
        *self | *rhs
    }
}

impl BitOr<Amino> for AmbiAmino {
    type Output = AmbiAmino;

    fn bitor(self, rhs: Amino) -> Self::Output {
        self | AmbiAmino::from(rhs)
    }
}

impl BitOr<&Amino> for &AmbiAmino {
    type Output = AmbiAmino;

    fn bitor(self, rhs: &Amino) -> Self::Output {
        *self | *rhs
    }
}

impl BitOr<AmbiAmino> for Amino {
    type Output = AmbiAmino;

    fn bitor(self, rhs: AmbiAmino) -> Self::Output {
        AmbiAmino::from(self) | rhs
    }
}

impl BitOr<&AmbiAmino> for &Amino {
    type Output = AmbiAmino;

    fn bitor(self, rhs: &AmbiAmino) -> Self::Output {
        *self | *rhs
    }
}

impl BitOr for Amino {
    type Output = AmbiAmino;

    fn bitor(self, rhs: Amino) -> Self::Output {
        AmbiAmino::from(self) | rhs
    }
}

impl BitOr for &Amino {
    type Output = AmbiAmino;

    fn bitor(self, rhs: &Amino) -> Self::Output {
        *self | *rhs
    }
}

impl BitOrAssign for AmbiAmino {
    fn bitor_assign(&mut self, rhs: AmbiAmino) {
        *self = *self | rhs;
    }
}

impl BitOrAssign<&AmbiAmino> for AmbiAmino {
    fn bitor_assign(&mut self, rhs: &AmbiAmino) {
        *self |= *rhs;
    }
}

impl BitOrAssign<Amino> for AmbiAmino {
    fn bitor_assign(&mut self, rhs: Amino) {
        *self = *self | rhs;
    }
}

impl BitOrAssign<&Amino> for AmbiAmino {
    fn bitor_assign(&mut self, rhs: &Amino) {
        *self |= *rhs;
    }
}

impl PartialEq<Amino> for AmbiAmino {
    fn eq(&self, other: &Amino) -> bool {
        *self == AmbiAmino::from(*other)
    }
}

impl PartialEq<AmbiAmino> for Amino {
    fn eq(&self, other: &AmbiAmino) -> bool {
        other == self
    }
}

impl PartialOrd<Amino> for AmbiAmino {
    fn partial_cmp(&self, other: &Amino) -> Option<Ordering> {
        Some(self.cmp(&AmbiAmino::from(*other)))
    }
}

impl PartialOrd<AmbiAmino> for Amino {
    fn partial_cmp(&self, other: &AmbiAmino) -> Option<Ordering> {
        other.partial_cmp(self).map(Ordering::reverse)
    }
}

impl From<Amino> for AmbiAmino {
    fn from(amino: Amino) -> Self {
        Self(NonZeroU32::MIN | Self::bit_mask(amino))
    }
}

impl TryFrom<AmbiAmino> for Amino {
    type Error = AmbiAmino;

    fn try_from(amino: AmbiAmino) -> Result<Self, Self::Error> {
        Ok(match amino {
            AmbiAmino::STOP => Self::Stop,
            AmbiAmino::A => Amino::A,
            AmbiAmino::C => Amino::C,
            AmbiAmino::D => Amino::D,
            AmbiAmino::E => Amino::E,
            AmbiAmino::F => Amino::F,
            AmbiAmino::G => Amino::G,
            AmbiAmino::H => Amino::H,
            AmbiAmino::I => Amino::I,
            AmbiAmino::K => Amino::K,
            AmbiAmino::L => Amino::L,
            AmbiAmino::M => Amino::M,
            AmbiAmino::N => Amino::N,
            AmbiAmino::O => Amino::O,
            AmbiAmino::P => Amino::P,
            AmbiAmino::Q => Amino::Q,
            AmbiAmino::R => Amino::R,
            AmbiAmino::S => Amino::S,
            AmbiAmino::T => Amino::T,
            AmbiAmino::U => Amino::U,
            AmbiAmino::V => Amino::V,
            AmbiAmino::W => Amino::W,
            AmbiAmino::Y => Amino::Y,
            _ => return Err(amino),
        })
    }
}

#[cfg(any(feature = "proptest", feature = "rand", test))]
impl AmbiAmino {
    pub(crate) const BITS_RANGE: std::ops::Range<u32> = 1..(1 << Amino::ALL.len());

    /// Note: `bits` must be in [`Self::BITS_RANGE`].
    pub(crate) fn from_bits(bits: u32) -> AmbiAmino {
        // The bit layout of all valid AmbiAminos:
        // 000000x0 xxxxxxxx xxxxxxxx xxxxx0x1
        //          <--------21 bits------>
        // Note: Big-endian, and x can be any bit (with the contraint at least one x is true)

        // Given a random number, consisting of 23 (low-order) bits, we need to spread them
        // to match the bit pattern above... Easiest approach is move the second bit up.
        Self(NonZeroU32::MIN | ((bits & 2) << 24) | ((bits & !2) << 1))
    }
}

/// Displays a single-character code, falling back to `X`.
///
/// # Examples
///
/// ```
/// use nucs::AmbiAmino;
///
/// assert_eq!(AmbiAmino::S.to_string(), "S");
/// assert_eq!((AmbiAmino::L | AmbiAmino::I).to_string(), "J");
/// assert_eq!((AmbiAmino::A | AmbiAmino::B).to_string(), "X");
/// assert_eq!(AmbiAmino::X.to_string(), "X");
/// // Note that despite the same representation they're not actually the same:
/// assert_ne!(AmbiAmino::A | AmbiAmino::B, AmbiAmino::X);
/// ```
impl Display for AmbiAmino {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        self.to_str().fmt(f)
    }
}

/// Displays single-character code when possible, otherwise all possibilities.
///
/// Note that when falling back to bracket-enclosed possibilities, this doesn't use
/// ambiguity codes like [`B`](Self::B), [`J`](Self::J) or [`Z`](Self::Z).
///
// # Examples
///
/// ```
/// use nucs::AmbiAmino;
///
/// assert_eq!(format!("{:?}", AmbiAmino::S), "S");
/// assert_eq!(format!("{:?}", AmbiAmino::L | AmbiAmino::I), "J");
/// assert_eq!(format!("{:?}", AmbiAmino::A | AmbiAmino::B), "[ADN]");
/// assert_eq!(format!("{:?}", AmbiAmino::X), "X");
/// ```
impl std::fmt::Debug for AmbiAmino {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        let chr = self.to_str();
        if *self != Self::X && chr == "X" {
            f.write_str("[")?;
            for amino in self {
                f.write_str(amino.to_str())?;
            }
            f.write_str("]")
        } else {
            f.write_str(chr)
        }
    }
}

impl FromStr for AmbiAmino {
    type Err = ParseSymbolError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.as_bytes() {
            [b] => Self::from_ascii(*b),
            _ => Err(ParseSymbolError::new::<Self>()),
        }
    }
}

impl AsRef<AmbiAmino> for AmbiAmino {
    fn as_ref(&self) -> &AmbiAmino {
        self
    }
}

impl AsMut<AmbiAmino> for AmbiAmino {
    fn as_mut(&mut self) -> &mut AmbiAmino {
        self
    }
}

impl IntoIterator for AmbiAmino {
    type IntoIter = AmbiAminoIter;
    type Item = Amino;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

impl IntoIterator for &AmbiAmino {
    type IntoIter = AmbiAminoIter;
    type Item = Amino;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

impl Symbol for Amino {
    type Concrete = Amino;
    type Ambiguous = AmbiAmino;

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

impl Symbol for AmbiAmino {
    type Concrete = Amino;
    type Ambiguous = AmbiAmino;

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

impl crate::symbol::sealed::Sealed for Amino {
    const NAME: &str = "amino acid";
    const EXPECTED: &str =
        "one of A/C/D/E/F/G/H/I/K/L/M/N/O/P/Q/R/S/T/U/V/W/Y/* (case-insensitive)";
}

impl crate::symbol::sealed::Sealed for AmbiAmino {
    const NAME: &str = "ambiguous amino acid";
    const EXPECTED: &str = "a letter or the character *";
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unambiguous_amino_is_its_only_possiblity() {
        for aa in Amino::ALL {
            assert!(AmbiAmino::from(aa).iter().eq([aa]));
        }
    }

    #[test]
    fn unambiguous_pair_of_aminos_are_their_only_possiblities() {
        for aa1 in Amino::ALL {
            for aa2 in Amino::ALL {
                if aa1 < aa2 {
                    assert!((aa1 | aa2).iter().eq([aa1, aa2]));
                }
            }
        }
    }

    fn all_ambi_aminos() -> impl Iterator<Item = AmbiAmino> {
        (1..(1 << Amino::ALL.len()) - 1).filter_map(|amino_i| {
            (0..Amino::ALL.len())
                .filter(|bit| amino_i & (1 << bit) != 0)
                .map(|bit| AmbiAmino::from(Amino::ALL[bit]))
                .reduce(|a, b| a | b)
        })
    }

    #[cfg_attr(debug_assertions, ignore = "slow outside of release mode")]
    #[test]
    fn amino_possiblities_are_sorted_and_unique() {
        for aa in all_ambi_aminos() {
            assert!(aa.iter().zip(aa.iter().skip(1)).all(|(a, b)| a < b));
        }
    }

    #[cfg_attr(debug_assertions, ignore = "slow outside of release mode")]
    #[test]
    fn ambi_amino_bitor_is_itempotent() {
        for aa in all_ambi_aminos() {
            assert_eq!(aa | aa, aa);
        }
    }

    #[cfg_attr(debug_assertions, ignore = "slow outside of release mode")]
    #[test]
    fn amino_possibilities_can_be_recomposed_into_ambi_nucs() {
        for aa in all_ambi_aminos() {
            let possibilities = aa.iter().map(AmbiAmino::from);
            assert_eq!(possibilities.reduce(|a, b| a | b), Some(aa));
        }
    }

    // Definitely way to many AmbiAminos to exhaustively test OR-ing them together...
    #[test]
    fn smoke_test_ambi_amino_bitor() {
        assert!(
            (AmbiAmino::A | AmbiAmino::B)
                .iter()
                .eq([Amino::A, Amino::D, Amino::N])
        );
    }

    // // WAY too many possibilities (70T) to exhaustively test. :/
    // #[ignore]
    // #[test]
    // fn ambi_amino_bitor_is_consistent_with_possibilities() {
    //     for ambi1 in all_ambi_aminos() {
    //         for ambi2 in all_ambi_aminos() {
    //             let mut expected: Vec<_> = ambi1.iter().chain(ambi2).collect();
    //             expected.sort();
    //             expected.dedup();
    //             assert!((ambi1 | ambi2).iter().eq(expected));
    //         }
    //     }
    // }

    #[test]
    fn amino_bitor_with_x_is_x() {
        for aa in Amino::ALL {
            assert_eq!(aa | AmbiAmino::X, AmbiAmino::X);
            assert_eq!(AmbiAmino::X | aa, AmbiAmino::X);
        }
    }

    #[cfg_attr(debug_assertions, ignore = "slow outside of release mode")]
    #[test]
    fn ambi_amino_bitor_with_x_is_x() {
        for aa in all_ambi_aminos() {
            assert_eq!(aa | AmbiAmino::X, AmbiAmino::X);
            assert_eq!(AmbiAmino::X | aa, AmbiAmino::X);
        }
    }

    #[test]
    fn str_roundtrips() {
        for aa in Amino::ALL {
            assert_eq!(Amino::from_str(aa.to_str()), Ok(aa));
        }
        // Notably, this test is NOT applicable to `AmbiAmino`
    }

    #[test]
    fn ascii_roundtrips() {
        for aa in Amino::ALL {
            assert_eq!(Amino::from_ascii(aa.to_ascii()), Ok(aa));
        }
        // Notably, this test is NOT applicable to `AmbiAmino`
    }

    #[test]
    fn all_is_sorted() {
        let mut sorted = Amino::ALL;
        sorted.sort();
        assert_eq!(Amino::ALL, sorted);
    }

    #[test]
    fn sanity_check_from_bits_values() {
        // AmbiAmino::X has all bits enabled, and it's generated at compile-time via a slow
        // but fairly fool-proof method, so it makes a good reference value.
        // The last value in BITS_RANGE should have all bits turned on, so hopefully passing
        // that to AmbiAmino::from_bits should result in AmbiAmino::X.
        let all_bits = AmbiAmino::BITS_RANGE.last().unwrap();
        assert_eq!(AmbiAmino::from_bits(all_bits), AmbiAmino::X);

        // Make sure turning each bit on individually results in the same set of values
        // as the concrete Aminos.
        let mut concrete_ambi_aminos = std::array::from_fn(|i| AmbiAmino::from_bits(1u32 << i));
        let mut expected = Amino::ALL.map(AmbiAmino::from);
        concrete_ambi_aminos.sort();
        expected.sort();
        assert_eq!(concrete_ambi_aminos, expected);
    }
}
