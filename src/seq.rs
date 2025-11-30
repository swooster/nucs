use std::fmt::{Debug, Display, Formatter};
use std::hash::Hash;
use std::ops::{Deref, DerefMut};
use std::str::FromStr;

use ref_cast::{RefCastCustom, ref_cast_custom};

use crate::error::ParseSeqError;
use crate::translation::GeneticCode;
use crate::{DnaIter, DnaSlice, Nucleotide, Symbol};

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
    /// over `Seq(slice)`. For convenience, there's also [`DnaSlice::as_seq`].
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{Nuc, Seq};
    ///
    /// let dna = Nuc::lit(b"ACTGACTG");
    /// let partial_dna = Seq::wrap(&dna[3..6]);
    /// assert_eq!(partial_dna, "GAC");
    /// ```
    #[ref_cast_custom]
    pub fn wrap(slice: &T) -> &Self;

    /// Wrap a mutable reference in a [`Seq`].
    ///
    /// Most of `Seq`'s features (e.g. [`Display`] impl or string comparisons) work for
    /// `&mut Seq<[T]>` but not `Seq<&mut [T]>`, so `Seq::wrap_mut(slice)` should be prefered
    /// over `Seq(slice)`. For convenience, there's also [`DnaSlice::as_seq_mut`].
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{Nuc, Seq};
    ///
    /// let mut dna = Nuc::seq(b"ACTGACTG");
    /// let partial_dna = Seq::wrap_mut(&mut dna[3..6]);
    /// assert_eq!(partial_dna, "GAC");
    /// partial_dna[1] = Nuc::C;
    /// assert_eq!(dna, "ACTGCCTG");
    /// ```
    #[ref_cast_custom]
    pub fn wrap_mut(slice: &mut T) -> &mut Self;

    /// Translate codons into [`Seq`]-wrapped peptide.
    ///
    /// This should work with a wider variety of collections, but may be slower as a result.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::collections::VecDeque;
    /// use nucs::{NCBI1, Nuc, Seq};
    ///
    /// let dna = Nuc::seq(b"TATGCGAGAAACA");
    /// let peptide: Seq<VecDeque<_>> = dna.translated_by(NCBI1);
    /// assert_eq!(peptide, "YARN");
    /// ```
    #[allow(
        clippy::needless_pass_by_value,
        reason = "consistency with other methods, plus G is likely Copy"
    )]
    pub fn translated_by<S, G, U>(&self, genetic_code: G) -> Seq<U>
    where
        for<'a> &'a T: IntoIterator<Item = &'a S>,
        S: Nucleotide,
        G: GeneticCode,
        U: FromIterator<S::Amino>,
    {
        Seq(self
            .0
            .into_iter()
            .codons()
            .map(|codon| genetic_code.translate(codon))
            .collect())
    }

    /// Translate codons into fixed-length [`Seq`]-wrapped peptide.
    ///
    /// # Panics
    ///
    /// Panics if the number of codons to be translated is different from the returned array.
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{NCBI1, Nuc, Seq};
    ///
    /// let dna = Nuc::seq(b"TATGCGAGAAACA");
    /// let peptide: Seq<[_; 4]> = dna.translated_to_array_by(NCBI1);
    /// assert_eq!(peptide, "YARN");
    /// ```
    pub fn translated_to_array_by<S, G, const N: usize>(
        &self,
        genetic_code: G,
    ) -> Seq<[<<[S] as DnaSlice>::Nuc as Nucleotide>::Amino; N]>
    where
        T: AsRef<[S]>,
        [S]: DnaSlice,
        G: GeneticCode,
    {
        Seq(self.0.as_ref().translated_to_array_by(genetic_code))
    }

    /// Translate codons into fixed-length [`Seq`]-wrapped  peptide in reverse order.
    ///
    /// This translates codons starting at the end. If the DNA can be converted to codons without
    /// excess nucleotides, then this produces the exact reverse of the output of
    /// [`translated_to_array_by`](Self::translated_to_array_by). It's intended to be used with
    /// [`FastTranslator::reverse_complement`](crate::translation::FastTranslator::reverse_complement)
    /// as that speeds up translation by folding the complementation into the translator.
    ///
    /// <div class="warning">
    ///
    /// **BEWARE:** This translates the *codons* in reverse order, *not* the nucleotides.
    /// The `SDRAWKCAB`/`BACKWARDS` example below demonstrates this.
    ///
    /// </div>
    ///
    /// # Panics
    ///
    /// Panics if the number of codons to be translated is different from the returned array.
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{AmbiNuc, NCBI1, NCBI1_RC, Seq};
    ///
    /// let dna = AmbiNuc::seq(b"AGCGATAGGGCCTGGAAATGTGCCRAY");
    /// let peptide: Seq<[_; 9]> = dna.translated_to_array_by(NCBI1);
    /// assert_eq!(peptide, "SDRAWKCAB");
    /// let peptide: Seq<[_; 9]> = dna.rev_translated_to_array_by(NCBI1);
    /// assert_eq!(peptide, "BACKWARDS");
    ///
    /// // The proper way to use this for RC translation is with a reverse-complemented
    /// // translation table like `NCBI1_RC`.
    /// let dna = AmbiNuc::seq(b"NGCACCGCTAGGTACTGGCGAA");
    /// let peptide: Seq<[_; 7]> = dna.rev_translated_to_array_by(NCBI1_RC);
    /// assert_eq!(peptide, "FAST*RC");
    /// ```
    pub fn rev_translated_to_array_by<S, G, const N: usize>(
        &self,
        genetic_code: G,
    ) -> Seq<[<<[S] as DnaSlice>::Nuc as Nucleotide>::Amino; N]>
    where
        T: AsRef<[S]>,
        [S]: DnaSlice,
        G: GeneticCode,
    {
        Seq(self.0.as_ref().rev_translated_to_array_by(genetic_code))
    }

    /// Translate codons into [`Seq`]-wrapped peptide [`Vec`].
    ///
    /// For large sequences, this is usually much faster than populating directly from an iterator.
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{NCBI1, Nuc};
    ///
    /// let dna = Nuc::seq(b"TATGCGAGAAACA");
    /// let peptide = dna.translated_to_vec_by(NCBI1);
    /// assert_eq!(peptide, "YARN");
    /// ```
    pub fn translated_to_vec_by<S, G>(
        &self,
        genetic_code: G,
    ) -> Seq<Vec<<<[S] as DnaSlice>::Nuc as Nucleotide>::Amino>>
    where
        T: AsRef<[S]>,
        [S]: DnaSlice,
        G: GeneticCode,
    {
        Seq(self.0.as_ref().translated_to_vec_by(genetic_code))
    }

    /// Translate codons into [`Seq`]-wrapped peptide [`Vec`] in reverse order.
    ///
    /// This translates codons starting at the end. If the DNA can be converted to codons without
    /// excess nucleotides, then this produces the exact reverse of the output of
    /// [`translated_to_vec_by`](Self::translated_to_vec_by). It's intended to be used with
    /// [`FastTranslator::reverse_complement`](crate::translation::FastTranslator::reverse_complement)
    /// as that speeds up translation by folding the complementation into the translator.
    ///
    /// <div class="warning">
    ///
    /// **BEWARE:** This translates the *codons* in reverse order, *not* the nucleotides.
    /// The `SDRAWKCAB`/`BACKWARDS` example below demonstrates this.
    ///
    /// </div>
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{AmbiNuc, NCBI1, NCBI1_RC};
    ///
    /// let dna = AmbiNuc::seq(b"AGCGATAGGGCCTGGAAATGTGCCRAY");
    /// let peptide = dna.translated_to_vec_by(NCBI1);
    /// assert_eq!(peptide, "SDRAWKCAB");
    /// let peptide = dna.rev_translated_to_vec_by(NCBI1);
    /// assert_eq!(peptide, "BACKWARDS");
    ///
    /// // The proper way to use this for RC translation is with a reverse-complemented
    /// // translation table like `NCBI1_RC`.
    /// let dna = AmbiNuc::seq(b"NGCACCGCTAGGTACTGGCGAA");
    /// let peptide = dna.rev_translated_to_vec_by(NCBI1_RC);
    /// assert_eq!(peptide, "FAST*RC");
    /// ```
    pub fn rev_translated_to_vec_by<S, G>(
        &self,
        genetic_code: G,
    ) -> Seq<Vec<<<[S] as DnaSlice>::Nuc as Nucleotide>::Amino>>
    where
        T: AsRef<[S]>,
        [S]: DnaSlice,
        G: GeneticCode,
    {
        Seq(self.0.as_ref().rev_translated_to_vec_by(genetic_code))
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

    use crate::{NCBI1, Nuc, Seq};

    #[test]
    fn sanity_check_that_seq_works_with_arrays() {
        let dna = Nuc::seq(b"ACGT");
        assert_eq!(dna, "ACGT");
        let peptide = dna.translated_to_vec_by(NCBI1);
        assert_eq!(peptide, "T");
    }

    #[test]
    fn sanity_check_that_seq_works_with_vecs() {
        let dna = Seq(Nuc::lit(b"ACGT").to_vec());
        assert_eq!(dna, "ACGT");
        let peptide = dna.translated_to_vec_by(NCBI1);
        assert_eq!(peptide, "T");
    }

    #[test]
    fn sanity_check_that_seq_works_with_slices() {
        let dna = Seq::wrap(const { &Nuc::lit(b"ACGT") });
        assert_eq!(dna, "ACGT");
        let peptide = dna.translated_to_vec_by(NCBI1);
        assert_eq!(peptide, "T");
    }

    #[test]
    fn sanity_check_that_seq_works_with_vecdeques() {
        type VdSeq<T> = Seq<VecDeque<T>>;
        let dna = VdSeq::from_iter(Nuc::lit(b"ACGT"));
        assert_eq!(dna, "ACGT");
        let peptide: VdSeq<_> = dna.translated_by(NCBI1);
        assert_eq!(peptide, "T");
    }
}
