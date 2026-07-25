use std::borrow::{Borrow, BorrowMut};
use std::fmt::{Debug, Display, Formatter};
use std::hash::Hash;
use std::ops::{Deref, DerefMut, Index, IndexMut};
use std::str::FromStr;

use ref_cast::{RefCastCustom, ref_cast_custom};

use crate::error::ParseSeqError;
use crate::iter::{Codons, Translated};
use crate::slice::Translation;
use crate::translation::GeneticCode;
use crate::{DnaIterExt, DnaSliceExt, Nucleotide, Symbol};

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
/// assert_eq!(dna, "ACATATAC");
/// // If alternate formatting is used, the sequence is line-wrapped according to the width:
/// assert_eq!(format!("{dna:#5}"), "ACATA\nTAC");
///
/// // For testing convenience, you can directly compare `Seq<T>` to strings:
/// // (whitespace and case are ignored)
/// assert_eq!(dna, "    A cat ATAC");
///
/// // You can still work with the underlying collection:
/// dna[2..].fill(Nuc::G);
/// assert_eq!(dna, "ACGGGGGG");
///
/// // Other collections than `Vec` are supported:
/// let mut dna = Seq(VecDeque::from_iter(dna));
/// dna.push_front(Nuc::T);
/// assert_eq!(dna, "TACGGGGGG");
///
/// // `AmbiNuc` and `Amino` are supported as well:
/// let peptide: Seq<Vec<Amino>> = "INTEROP".parse()?;
/// assert_eq!(peptide, "INTEROP");
/// # Ok(())
/// # }
/// ```
///
/// # Requirements and limitations
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
///
/// # Features
///
/// * **`serde`:** When enabled, [`Seq<T>`] is serializable (as a string) whenever it has a
///   [`Display`] impl, and deserializable whenever it has a [`FromStr`] impl.
#[derive(Clone, Copy, Default, PartialEq, Eq, PartialOrd, Ord, Hash, RefCastCustom)]
#[repr(transparent)]
pub struct Seq<T: ?Sized>(pub T);

impl<T: ?Sized> Seq<T> {
    /// Wrap a reference in a [`Seq`].
    ///
    /// Most of `Seq`'s features (e.g. [`Display`] impl or string comparisons) work for
    /// `&Seq<[T]>` but not `Seq<&[T]>`, so `Seq::wrap(slice)` should be prefered
    /// over `Seq(slice)`. For convenience, there's also [`DnaSliceExt::as_seq`].
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{Nuc, Seq};
    ///
    /// let dna = Nuc::arr(b"ACTGACTG");
    /// let partial_dna = Seq::wrap(&dna[3..6]);
    /// assert_eq!(partial_dna, "GAC");
    /// ```
    #[ref_cast_custom]
    pub fn wrap(slice: &T) -> &Self;

    /// Wrap a mutable reference in a [`Seq`].
    ///
    /// Most of `Seq`'s features (e.g. [`Display`] impl or string comparisons) work for
    /// `&mut Seq<[T]>` but not `Seq<&mut [T]>`, so `Seq::wrap_mut(slice)` should be prefered
    /// over `Seq(slice)`. For convenience, there's also [`DnaSliceExt::as_seq_mut`].
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{Nuc, Seq};
    ///
    /// let mut dna = Nuc::seq(b"ACTGACTG");
    /// let partial_dna = &mut dna[3..6];
    /// assert_eq!(partial_dna, "GAC");
    /// partial_dna[1] = Nuc::C;
    /// assert_eq!(dna, "ACTGCCTG");
    /// ```
    #[ref_cast_custom]
    pub fn wrap_mut(slice: &mut T) -> &mut Self;

    /// Return translation of DNA via [`GeneticCode`].
    ///
    /// This returns a builder for performing translations. See [`Translation`] for details.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::collections::VecDeque;
    ///
    /// use nucs::{AmbiNuc, NCBI1, Nuc, Seq};
    ///
    /// let dna = Nuc::seq(b"TATGCGAGAAAC");
    /// assert_eq!(dna.translated_by(NCBI1).to_seq(), "YARN");
    ///
    /// let rc_dna = AmbiNuc::seq(b"NGCACCGCTAGGTACTGGCGAA");
    /// let peptide = rc_dna.translated_by(NCBI1).reverse_complemented().to_seq();
    /// assert_eq!(peptide, "FAST*RC");
    ///
    /// let peptide: Seq<VecDeque<_>> = dna.translated_by(NCBI1).into_iter().collect();
    /// assert_eq!(peptide, "YARN");
    /// ```
    pub fn translated_by<'a, S: 'a, G>(
        &'a self,
        genetic_code: G,
    ) -> Translation<'a, <[S] as DnaSliceExt>::Nuc, G>
    where
        T: AsRef<[S]>,
        [S]: DnaSliceExt,
        G: GeneticCode,
    {
        self.0.as_ref().translated_by(genetic_code)
    }

    /// Iterate over translated codons.
    ///
    /// This should work with a wider variety of collections than
    /// [`translated_by`](Self::translated_by), but may be slower as a result.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::collections::VecDeque;
    /// use nucs::{NCBI1, Nuc, Seq};
    ///
    /// let dna = Nuc::seq(b"TATGCGAGAAACA");
    /// let peptide: Seq<VecDeque<_>> = dna.iter_translated_by(NCBI1).collect();
    /// assert_eq!(peptide, "YARN");
    /// ```
    pub fn iter_translated_by<S, G>(
        &self,
        genetic_code: G,
    ) -> Translated<G, Codons<S, <&T as IntoIterator>::IntoIter>>
    where
        for<'a> &'a T: IntoIterator<Item = &'a S>,
        S: Nucleotide,
        G: GeneticCode,
    {
        self.0.into_iter().translated_by(genetic_code)
    }
}

impl<T: ?Sized> Deref for Seq<T> {
    type Target = T;

    fn deref(&self) -> &T {
        &self.0
    }
}

impl<T: ?Sized> DerefMut for Seq<T> {
    fn deref_mut(&mut self) -> &mut T {
        &mut self.0
    }
}

impl<T, I> Index<I> for Seq<T>
where
    T: Index<I, Output: SeqWrap + 'static>,
{
    type Output = <T::Output as SeqWrap>::Wrapped;

    fn index(&self, index: I) -> &Self::Output {
        self.0.index(index).wrap_slice()
    }
}

impl<T, I> IndexMut<I> for Seq<T>
where
    T: IndexMut<I, Output: SeqWrap + 'static>,
{
    fn index_mut(&mut self, index: I) -> &mut Self::Output {
        self.0.index_mut(index).wrap_mut_slice()
    }
}

impl<T: Borrow<[S]>, S> Borrow<Seq<[S]>> for Seq<T> {
    fn borrow(&self) -> &Seq<[S]> {
        Seq::wrap(self.0.borrow())
    }
}

impl<T: BorrowMut<[S]>, S> BorrowMut<Seq<[S]>> for Seq<T> {
    fn borrow_mut(&mut self) -> &mut Seq<[S]> {
        Seq::wrap_mut(self.0.borrow_mut())
    }
}

impl<S: Clone> ToOwned for Seq<[S]> {
    type Owned = Seq<Vec<S>>;

    fn to_owned(&self) -> Self::Owned {
        Seq(self.0.to_owned())
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

impl<T: ?Sized, S> PartialEq<&str> for Seq<T>
where
    for<'a> &'a T: IntoIterator<Item = &'a S>,
    S: Symbol,
{
    fn eq(&self, rhs: &&str) -> bool {
        self == *rhs
    }
}

impl<T: ?Sized, S> PartialEq<str> for Seq<T>
where
    for<'a> &'a T: IntoIterator<Item = &'a S>,
    S: Symbol,
{
    fn eq(&self, rhs: &str) -> bool {
        self.into_iter().copied().map(Ok).eq(iter_symbols(rhs))
    }
}

impl<T: ?Sized> Display for Seq<T>
where
    for<'a> &'a T: IntoIterator<Item: Display>,
{
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        Display::fmt(&crate::iter::Display::new(&self.0), f)
    }
}

impl<T: ?Sized> Debug for Seq<T>
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

/// Utility trait to wrap slices in [`Seq`] while leaving other types alone.
///
/// Used in [`Seq`]'s implementation of [`Index`].
pub trait SeqWrap {
    /// `Seq<Self>` for `[T]`, otherwise `Self`.
    type Wrapped: ?Sized;

    /// Wrap slices in [`Seq`], but returns other types unchanged.
    fn wrap_slice(&self) -> &Self::Wrapped;
    /// Wrap mutable slices in [`Seq`], but returns other types unchanged.
    fn wrap_mut_slice(&mut self) -> &mut Self::Wrapped;
}

impl<S> SeqWrap for [S] {
    type Wrapped = Seq<Self>;

    fn wrap_slice(&self) -> &Self::Wrapped {
        Seq::wrap(self)
    }

    fn wrap_mut_slice(&mut self) -> &mut Self::Wrapped {
        Seq::wrap_mut(self)
    }
}

impl<T> SeqWrap for T {
    type Wrapped = Self;

    fn wrap_slice(&self) -> &Self::Wrapped {
        self
    }

    fn wrap_mut_slice(&mut self) -> &mut Self::Wrapped {
        self
    }
}

#[cfg(feature = "serde")]
mod serde_impls {
    use std::fmt::{Display, Formatter};
    use std::marker::PhantomData;

    use serde::{Deserialize, Deserializer, Serialize, Serializer, de::Visitor};

    use super::{Seq, Symbol};

    impl<T> Serialize for Seq<T>
    where
        Self: Display,
    {
        fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
            serializer.serialize_str(&self.to_string())
        }
    }

    impl<'de, T, U> Deserialize<'de> for Seq<T>
    where
        // Need IntoIterator bound to infer type of contained element
        for<'a> &'a T: IntoIterator<Item = &'a U>,
        T: FromIterator<U>,
        U: Symbol,
    {
        fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
            deserializer.deserialize_str(SeqVisitor(PhantomData))
        }
    }

    struct SeqVisitor<T>(PhantomData<T>);

    impl<T, U> Visitor<'_> for SeqVisitor<T>
    where
        // Need IntoIterator bound to infer type of contained element
        for<'a> &'a T: IntoIterator<Item = &'a U>,
        T: FromIterator<U>,
        U: Symbol,
    {
        type Value = Seq<T>;

        fn expecting(&self, f: &mut Formatter) -> std::fmt::Result {
            write!(f, "a string of {}s", U::NAME)
        }

        fn visit_str<E: serde::de::Error>(self, v: &str) -> Result<Self::Value, E> {
            v.parse().map_err(E::custom)
        }

        // TODO: Maybe also accept sequences of symbols?
    }
}

#[cfg(all(test, feature = "serde"))]
mod serde_tests {
    use crate::Dna;

    #[test]
    fn dna_roundtrip() {
        let original_dna: Dna = "CATTAG".parse().unwrap();
        let json = serde_json::to_string(&original_dna).unwrap();
        assert_eq!(json, "\"CATTAG\"");
        let deserialized_dna: Dna = serde_json::from_str(&json).unwrap();
        assert_eq!(deserialized_dna, original_dna);
    }

    #[test]
    fn invalid_dna() {
        let err = serde_json::from_str::<Dna>("[]").unwrap_err();
        assert!(err.to_string().contains("expected a string of nucleotides"));
        let err = serde_json::from_str::<Dna>("\"CATXTAG\"").unwrap_err();
        assert!(
            err.to_string()
                .contains("invalid nucleotide 'X' at position 3")
        );
    }
}

#[cfg(test)]
mod tests {
    use std::collections::VecDeque;

    use crate::{Dna, NCBI1, Nuc, Seq};

    #[test]
    fn sanity_check_that_seq_works_with_arrays() {
        let mut dna = Nuc::seq(b"ACGT");
        assert_eq!(dna, "ACGT");
        let peptide = dna.translated_by(NCBI1).to_seq();
        assert_eq!(peptide, "T");
        dna[2] = Nuc::A;
        assert_eq!(dna[1..], "CAT");
        let owned: Dna = dna[2..].to_owned();
        assert_eq!(owned, "AT");
    }

    #[test]
    fn sanity_check_that_seq_works_with_vecs() {
        let mut dna = Seq(Nuc::arr(b"ACGT").to_vec());
        assert_eq!(dna, "ACGT");
        let peptide = dna.translated_by(NCBI1).to_seq();
        assert_eq!(peptide, "T");
        dna[2] = Nuc::A;
        assert_eq!(dna[1..], "CAT");
        let owned: Dna = dna[2..].to_owned();
        assert_eq!(owned, "AT");
    }

    #[test]
    fn sanity_check_that_seq_works_with_slices() {
        let dna = Seq::wrap(const { &Nuc::arr(b"ACGT") });
        assert_eq!(dna, "ACGT");
        let peptide = dna.translated_by(NCBI1).to_seq();
        assert_eq!(peptide, "T");
        assert_eq!(dna[1..], "CGT");
        let owned: Dna = dna[2..].to_owned();
        assert_eq!(owned, "GT");
    }

    #[test]
    fn sanity_check_that_seq_works_with_vecdeques() {
        type VdSeq<T> = Seq<VecDeque<T>>;
        let mut dna = VdSeq::from_iter(Nuc::arr(b"ACGT"));
        assert_eq!(dna, "ACGT");
        let peptide: VdSeq<_> = dna.iter_translated_by(NCBI1).collect();
        assert_eq!(peptide, "T");
        dna[2] = Nuc::A;
        assert_eq!(dna, "ACAT");
    }
}
