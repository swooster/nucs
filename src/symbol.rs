use std::fmt::Display;
use std::hash::Hash;
use std::str::FromStr;

use crate::error::ParseSymbolError;

/// A sequence element; either [`Nuc`](crate::Nuc), [`AmbiNuc`](crate::AmbiNuc),
/// [`Amino`](crate::Amino) or [`AmbiAmino`](crate::AmbiAmino).
pub trait Symbol:
    'static
    + Clone
    + Copy
    + Default
    + Display
    + std::fmt::Debug
    + FromStr<Err = ParseSymbolError>
    + AsRef<Self>
    + AsMut<Self>
    + PartialEq
    + Eq
    + PartialOrd
    + Ord
    + PartialEq<Self::Concrete>
    + PartialOrd<Self::Concrete>
    + PartialEq<Self::Ambiguous>
    + PartialOrd<Self::Ambiguous>
    + Into<Self::Ambiguous>
    + TryInto<Self::Concrete>
    + Hash
    + sealed::Sealed
{
    /// Concrete symbols, i.e. [`Nuc`](crate::Nuc) or [`Amino`](crate::Amino)
    type Concrete: Symbol;
    /// Ambiguous symbols, i.e. [`AmbiNuc`](crate::AmbiNuc) or [`AmbiAmino`](crate::AmbiAmino)
    type Ambiguous: Symbol;

    /// Return uppercase string representation
    #[must_use]
    fn to_str(self) -> &'static str;

    /// Construct from (case-insensitive) ASCII representation
    ///
    /// # Errors
    ///
    /// Returns [`ParseSymbolError`] if the given byte isn't a valid (case-insensitive).
    fn from_ascii(ascii: u8) -> Result<Self, ParseSymbolError>;

    /// Return uppercase ASCII representation
    #[must_use]
    fn to_ascii(self) -> u8;

    /// Construct array from literal without allocating.
    ///
    /// This is just intended to be a convenience method for testing.
    #[must_use]
    #[track_caller]
    fn lit<const N: usize>(literal: &[u8; N]) -> [Self; N];
}

pub(crate) mod sealed {

    pub trait Sealed: Sized {
        const NAME: &str;
        const EXPECTED: &str;

        // I admit it's super icky that these exist yet aren't applicable to `Amino`s...
        // I'm doing it this way so I don't need to have an additional `Sealed` trait for
        // nucleotides show up in docs. These are an internal detail so that `DnaSlice`
        // can have a blanket impl on `[impl Nucleotide]`, which makes trait inference
        // nicer for users.

        #[cfg(feature = "unsafe")]
        fn as_ambi_nucs(_nucs: &[Self]) -> &[crate::AmbiNuc] {
            unimplemented!()
        }

        #[cfg(feature = "unsafe")]
        fn to_nucs(_nucs: &[Self]) -> Option<&[crate::Nuc]> {
            unimplemented!()
        }

        #[cfg(feature = "unsafe")]
        fn to_nucs_mut(_nucs: &mut [Self]) -> Option<&mut [crate::Nuc]> {
            unimplemented!()
        }
    }
}
