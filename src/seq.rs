use std::fmt::{Debug, Display, Formatter};
use std::hash::Hash;
use std::ops::{Deref, DerefMut};
use std::str::FromStr;

use crate::Symbol;
use crate::error::ParseSeqError;

/// Provides DNA/peptide ergonomics for collections.
///
/// While containers like [`Vec<Nuc>`] are great for interop, their [`Debug`] representations
/// are a bit verbose for DNA or peptides and they don't support parsing.
/// [`Seq`] can wrap sufficiently [`Vec`]-like containers to provide such features.
///
/// # Examples
///
/// ```
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// use std::collections::VecDeque;
/// use nucs::{Amino, Nuc, Seq};
///
/// // `Seq` supports parsing strings. Whitespace is ignored:
/// let mut dna: Seq<Vec<Nuc>> = "
///     A
///     CAT
///     ATAC
/// ".parse()?;
/// // `Seq` implements `Display`:
/// assert_eq!(dna.to_string(), "ACATATAC");
/// // If alternate formatting is used, the sequence is line-wrapped according to the width:
/// assert_eq!(format!("{dna:#5}"), "ACATA\nTAC");
///
/// // You can still work with the underlying collection:
/// dna[2..].fill(Nuc::G);
/// assert_eq!(dna.to_string(), "ACGGGGGG");
///
/// // Other collections than `Vec` are supported:
/// let mut dna = Seq(VecDeque::from_iter(dna));
/// dna.push_front(Nuc::T);
/// assert_eq!(dna.to_string(), "TACGGGGGG");
///
/// // `AmbiNuc` and `Amino` are supported as well:
/// let peptide: Seq<Vec<Amino>> = "INTEROP".parse()?;
/// assert_eq!(peptide.to_string(), "INTEROP");
/// # Ok(())
/// # }
/// ```
///
/// # Requirements and limitations:
///
/// [`Seq<T>`] reads from its collection via [`&T::into_iter`](IntoIterator::into_iter), and
/// expects yielded items to be [`&impl Symbol`](Symbol). This means it works with [`Vec`],
/// [`VecDeque`](std::collections::VecDeque), [`[T; N]`](array), and custom [`Vec`]-like
/// containers such as `SmallVec`, `TinyVec`, `ArrayVec`, etc. **Unfortunately, it doesn't
/// work with [`&[T]`](slice) or [`Arc<[T]>`](std::sync::Arc)** because neither
/// [`&&[T]`](slice) nor [`&Arc<[T]>`](std::sync::Arc) are directly iterable without autoderef.
///
/// Note that parsing requires [`&T::into_iter`](IntoIterator::into_iter) to determine
/// what type of [`Symbol`] (e.g. [`Nuc`](crate::Nuc) vs [`Amino`](crate::Amino)) to expect;
/// [`FromIterator`] cannot be relied on, because a single collection may support multiple
/// [`FromIterator`] implementations.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Seq<T>(pub T);

impl<T> Deref for Seq<T> {
    type Target = T;

    fn deref(&self) -> &T {
        &self.0
    }
}

impl<T> DerefMut for Seq<T> {
    fn deref_mut(&mut self) -> &mut T {
        &mut self.0
    }
}

impl<T: FromIterator<A>, A> FromIterator<A> for Seq<T> {
    fn from_iter<U>(iter: U) -> Self
    where
        U: IntoIterator<Item = A>,
    {
        Self(T::from_iter(iter))
    }
}

impl<T: IntoIterator> IntoIterator for Seq<T> {
    type IntoIter = T::IntoIter;
    type Item = T::Item;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl<T> Display for Seq<T>
where
    for<'a> &'a T: IntoIterator<Item: Display>,
{
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        Display::fmt(&crate::iter::Display::new(&self.0), f)
    }
}

impl<T> Debug for Seq<T>
where
    for<'a> &'a T: IntoIterator<Item: Display>,
{
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        f.debug_tuple("Seq")
            .field(&crate::iter::Display::new(&self.0))
            .finish()
    }
}

impl<T, U> FromStr for Seq<T>
where
    // Need IntoIterator bound to infer type of contained element
    for<'a> &'a T: IntoIterator<Item = &'a U>,
    T: FromIterator<U>,
    U: Symbol,
{
    type Err = ParseSeqError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        iter_symbols(s).collect::<Result<_, _>>().map(Self)
    }
}

fn iter_symbols<S: Symbol>(s: &str) -> impl Iterator<Item = Result<S, ParseSeqError>> {
    s.split("")
        .filter(|c| !c.is_empty())
        .enumerate()
        .filter(|(_, c)| !c.trim().is_empty()) // ignore whitespace; maybe I shouldn't?
        .map(|(pos, chr)| {
            chr.parse().map_err(|_| ParseSeqError {
                kind: S::NAME,
                expected: S::EXPECTED,
                chr: chr.chars().next().expect("BUG: chr was impossibly empty"),
                pos,
            })
        })
}
