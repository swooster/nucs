//! Error types

use std::error::Error;
use std::fmt::{Debug, Display, Formatter};

use crate::Symbol;

/// Error parsing [`Nuc`](crate::Nuc), [`AmbiNuc`](crate::AmbiNuc),
/// [`Amino`](crate::Amino) or [`AmbiAmino`](crate::AmbiAmino)
#[derive(Clone, PartialEq, Eq)]
pub struct ParseSymbolError {
    pub(crate) kind: &'static str,
    pub(crate) expected: &'static str,
}

impl ParseSymbolError {
    pub(crate) const fn new<S: Symbol>() -> Self {
        ParseSymbolError {
            kind: S::NAME,
            expected: S::EXPECTED,
        }
    }
}

impl Display for ParseSymbolError {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        let Self { kind, expected } = self;
        write!(f, "invalid {kind}: expected {expected}")
    }
}

impl Debug for ParseSymbolError {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        f.debug_struct("ParseSymbolError").finish_non_exhaustive()
    }
}

impl Error for ParseSymbolError {}

/// Error parsing sequence of [`Nuc`](crate::Nuc), [`AmbiNuc`](crate::AmbiNuc),
/// [`Amino`](crate::Amino) or [`AmbiAmino`](crate::AmbiAmino)
#[derive(Clone, PartialEq, Eq)]
pub struct ParseSeqError {
    pub(crate) kind: &'static str,
    pub(crate) expected: &'static str,
    pub(crate) chr: char,
    pub(crate) pos: usize,
}

impl ParseSeqError {
    /// The position (measured in codepoints, not bytes) where parsing failed.
    #[must_use]
    pub fn position(&self) -> usize {
        self.pos
    }

    /// The unexpected character that cause parsing to fail.
    #[must_use]
    pub fn char(&self) -> char {
        self.chr
    }
}

impl Display for ParseSeqError {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        let Self {
            kind,
            expected,
            chr,
            pos,
        } = self;
        write!(
            f,
            "invalid {kind} {chr:?} at position {pos}: expected {expected}"
        )
    }
}

impl Debug for ParseSeqError {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        f.debug_struct("ParseSeqError")
            .field("chr", &self.chr)
            .field("pos", &self.pos)
            .finish_non_exhaustive()
    }
}

impl Error for ParseSeqError {}
