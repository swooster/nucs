//! Types related to packed ambiguous DNA.

use std::cell::Cell;
use std::cmp::Ordering;
use std::fmt::Formatter;
use std::ops::{Deref, DerefMut};

use crate::symbol::iter_symbols;
use crate::{AmbiNuc, Seq, iter::display};

use super::dna::{PackedArrayDna, PackedDna};
use super::packable_array::{ArrayDefault, Sealed as ArrayDivide};
use super::{PackableArray, RefCmp, UnpackingIter};

// Note on storage: AmbiNucs are packed big-endian so naive lexical sorting of bytes is correct.
//
// The last byte of a non-empty `PackedAmbiDna` must have 1-2 elements. If it has 1 element,
// the least-significant 4 bits are b0000, which isn't a valid representation.

/// Like [`Vec<AmbiNuc>`], but takes 50% less space for long DNA sequences.
///
/// # Examples
///
/// ```
/// use nucs::{AmbiNuc, Packed};
///
/// let dna = AmbiNuc::arr(b"ACATATTACK");
/// let packed: Packed<[AmbiNuc]> = dna.as_slice().into();
/// let unpacked: Vec<AmbiNuc> = packed.into();
/// assert_eq!(unpacked, dna);
/// ```
#[derive(Clone, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct PackedAmbiDna(Vec<u8>);

impl PackedAmbiDna {
    /// Returns the number of [`AmbiNuc`]s in the packed DNA.
    #[must_use]
    pub fn len(&self) -> usize {
        // len <= `usize::MAX` is an invariant of this type, so no need to worry about overflow.
        match &*self.0 {
            [] => 0,
            bulk @ [.., tail] => 2 * bulk.len() - (tail.trailing_zeros() / 4) as usize,
        }
    }

    /// Returns `true` if the packed DNA contains no [`AmbiNuc`]s.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Appends an [`AmbiNuc`] to the back of the DNA.
    ///
    /// # Panics
    ///
    /// Panics if the new length exceeds [`isize::MAX`].
    ///
    /// # Examples
    ///
    /// ```
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// use nucs::{AmbiDna, AmbiNuc};
    ///
    /// let dna: AmbiDna = "BARN".parse()?;
    /// let mut packed = dna.pack();
    /// packed.push(AmbiNuc::C);
    /// packed.push(AmbiNuc::A);
    /// packed.push(AmbiNuc::T);
    /// assert_eq!(packed, "BARNCAT");
    /// # Ok(())
    /// # }
    /// ```
    pub fn push(&mut self, ambi_nuc: AmbiNuc) {
        assert_ne!(self.len(), isize::MAX as usize);
        match &mut *self.0 {
            [.., last] if last.trailing_zeros() >= 4 => *last |= ambi_nuc.compress(),
            _ => self.0.push(ambi_nuc.compress() << 4),
        }
    }

    /// Removes the last [`AmbiNuc`] from the DNA and returns it, or [`None`] if it is empty.
    ///
    /// # Examples
    ///
    /// ```
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// use nucs::{AmbiDna, AmbiNuc};
    ///
    /// let dna: AmbiDna = "BAT".parse()?;
    /// let mut packed = dna.pack();
    /// assert_eq!(packed.pop(), Some(AmbiNuc::T));
    /// assert_eq!(packed.pop(), Some(AmbiNuc::A));
    /// assert_eq!(packed.pop(), Some(AmbiNuc::B));
    /// assert_eq!(packed.pop(), None);
    /// # Ok(())
    /// # }
    /// ```
    pub fn pop(&mut self) -> Option<AmbiNuc> {
        match &mut *self.0 {
            [] => None,
            [.., last] if last.trailing_zeros() >= 4 => {
                let ambi_nuc = AmbiNuc::decompress(*last >> 4);
                self.0.pop();
                Some(ambi_nuc)
            }
            [.., last] => {
                let ambi_nuc = AmbiNuc::decompress(*last);
                *last &= !0b1111;
                Some(ambi_nuc)
            }
        }
    }

    /// Returns an iterator over the packed DNA.
    #[must_use]
    pub fn iter(&self) -> PackedAmbiDnaIter<'_> {
        self.into_iter()
    }

    /// Returns a mutable iterator over the packed DNA.
    ///
    /// # Examples
    ///
    /// ```
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// use nucs::AmbiDna;
    ///
    /// let dna: AmbiDna = "VTYNGTA".parse()?;
    /// let mut packed = dna.pack();
    /// for mut nuc in &mut packed {
    ///     *nuc = nuc.complement();
    /// }
    /// assert_eq!(packed, "BARNCAT");
    /// # Ok(())
    /// # }
    /// ```
    #[must_use]
    pub fn iter_mut(&mut self) -> PackedAmbiDnaMutIter<'_> {
        self.into_iter()
    }
}

impl From<Seq<Vec<AmbiNuc>>> for PackedAmbiDna {
    fn from(ambi_dna: Seq<Vec<AmbiNuc>>) -> PackedAmbiDna {
        ambi_dna.0.into()
    }
}

impl From<Vec<AmbiNuc>> for PackedAmbiDna {
    fn from(ambi_dna: Vec<AmbiNuc>) -> PackedAmbiDna {
        (&ambi_dna).into()
    }
}

impl<T: AsRef<[AmbiNuc]> + ?Sized> From<&T> for PackedAmbiDna {
    fn from(ambi_dna: &T) -> PackedAmbiDna {
        let ambi_dna = ambi_dna.as_ref();
        let packed_len = ambi_dna.len().div_ceil(2);
        let mut packed = vec![0; packed_len];
        pack(&mut packed, ambi_dna);
        Self(packed)
    }
}

impl From<PackedAmbiDna> for Seq<Vec<AmbiNuc>> {
    fn from(packed_ambi_dna: PackedAmbiDna) -> Seq<Vec<AmbiNuc>> {
        Seq(packed_ambi_dna.into())
    }
}

impl From<&PackedAmbiDna> for Seq<Vec<AmbiNuc>> {
    fn from(packed_ambi_dna: &PackedAmbiDna) -> Seq<Vec<AmbiNuc>> {
        Seq(packed_ambi_dna.into())
    }
}

impl From<PackedAmbiDna> for Vec<AmbiNuc> {
    fn from(packed_ambi_dna: PackedAmbiDna) -> Vec<AmbiNuc> {
        (&packed_ambi_dna).into()
    }
}

impl From<&PackedAmbiDna> for Vec<AmbiNuc> {
    fn from(packed_ambi_dna: &PackedAmbiDna) -> Vec<AmbiNuc> {
        let mut ambi_dna = vec![AmbiNuc::default(); packed_ambi_dna.len()];
        unpack(&mut ambi_dna, &packed_ambi_dna.0);
        ambi_dna
    }
}

impl<'a> IntoIterator for &'a PackedAmbiDna {
    type Item = AmbiNuc;
    type IntoIter = PackedAmbiDnaIter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        PackedAmbiDnaIter(UnpackingIter::new(0..self.len(), &self.0))
    }
}

impl<'a> IntoIterator for &'a mut PackedAmbiDna {
    type Item = PackedAmbiNuc<'a>;
    type IntoIter = PackedAmbiDnaMutIter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        let len = self.len();
        let backing = Cell::from_mut(self.0.as_mut_slice()).as_slice_of_cells();
        PackedAmbiDnaMutIter(UnpackingIter::new(0..len, backing))
    }
}

impl IntoIterator for PackedAmbiDna {
    type Item = AmbiNuc;
    type IntoIter = PackedAmbiDnaIntoIter;

    fn into_iter(self) -> Self::IntoIter {
        PackedAmbiDnaIntoIter(UnpackingIter::new(0..self.len(), self.0))
    }
}

impl<const N: usize> PartialEq<PackedArrayAmbiDna<N>> for PackedAmbiDna
where
    [(); N]: PackableArray,
{
    fn eq(&self, other: &PackedArrayAmbiDna<N>) -> bool {
        self.0.eq(other.0.as_ref())
    }
}

impl<const N: usize> PartialOrd<PackedArrayAmbiDna<N>> for PackedAmbiDna
where
    [(); N]: PackableArray,
{
    fn partial_cmp(&self, other: &PackedArrayAmbiDna<N>) -> Option<Ordering> {
        self.0.as_slice().partial_cmp(other.0.as_ref())
    }
}

impl PartialEq<PackedDna> for PackedAmbiDna {
    fn eq(&self, other: &PackedDna) -> bool {
        other == self
    }
}

impl PartialOrd<PackedDna> for PackedAmbiDna {
    fn partial_cmp(&self, other: &PackedDna) -> Option<Ordering> {
        other.partial_cmp(self).map(Ordering::reverse)
    }
}

impl<const N: usize> PartialEq<PackedArrayDna<N>> for PackedAmbiDna
where
    [(); N]: PackableArray,
{
    fn eq(&self, other: &PackedArrayDna<N>) -> bool {
        other == self
    }
}

impl<const N: usize> PartialOrd<PackedArrayDna<N>> for PackedAmbiDna
where
    [(); N]: PackableArray,
{
    fn partial_cmp(&self, other: &PackedArrayDna<N>) -> Option<Ordering> {
        other.partial_cmp(self).map(Ordering::reverse)
    }
}

impl<T: PartialEq<AmbiNuc>, const M: usize> PartialEq<[T; M]> for PackedAmbiDna {
    fn eq(&self, other: &[T; M]) -> bool {
        self.eq(&other.as_slice())
    }
}

impl<T: PartialOrd<AmbiNuc>, const M: usize> PartialOrd<[T; M]> for PackedAmbiDna {
    fn partial_cmp(&self, other: &[T; M]) -> Option<Ordering> {
        self.partial_cmp(&other.as_slice())
    }
}

impl<T: PartialEq<AmbiNuc>, const M: usize> PartialEq<&mut [T; M]> for PackedAmbiDna {
    fn eq(&self, other: &&mut [T; M]) -> bool {
        self.eq(&other.as_slice())
    }
}

impl<T: PartialOrd<AmbiNuc>, const M: usize> PartialOrd<&mut [T; M]> for PackedAmbiDna {
    fn partial_cmp(&self, other: &&mut [T; M]) -> Option<Ordering> {
        self.partial_cmp(&other.as_slice())
    }
}

impl<T: PartialEq<AmbiNuc>, const M: usize> PartialEq<&[T; M]> for PackedAmbiDna {
    fn eq(&self, other: &&[T; M]) -> bool {
        self.eq(&other.as_slice())
    }
}

impl<T: PartialOrd<AmbiNuc>, const M: usize> PartialOrd<&[T; M]> for PackedAmbiDna {
    fn partial_cmp(&self, other: &&[T; M]) -> Option<Ordering> {
        self.partial_cmp(&other.as_slice())
    }
}

impl<T: PartialEq<AmbiNuc>> PartialEq<Vec<T>> for PackedAmbiDna {
    fn eq(&self, other: &Vec<T>) -> bool {
        self.eq(&&**other)
    }
}

impl<T: PartialOrd<AmbiNuc>> PartialOrd<Vec<T>> for PackedAmbiDna {
    fn partial_cmp(&self, other: &Vec<T>) -> Option<Ordering> {
        self.partial_cmp(&&**other)
    }
}

impl<T: PartialEq<AmbiNuc>> PartialEq<&mut [T]> for PackedAmbiDna {
    fn eq(&self, other: &&mut [T]) -> bool {
        self.eq(&&**other)
    }
}

impl<T: PartialOrd<AmbiNuc>> PartialOrd<&mut [T]> for PackedAmbiDna {
    fn partial_cmp(&self, other: &&mut [T]) -> Option<Ordering> {
        self.partial_cmp(&&**other)
    }
}

impl<T: PartialEq<AmbiNuc>> PartialEq<&[T]> for PackedAmbiDna {
    fn eq(&self, other: &&[T]) -> bool {
        // Without `TrustedLen` we don't benefit from stdlib's length-check specializations.
        self.len() == other.len() && self.iter().map(RefCmp).eq(*other)
    }
}

impl<T: PartialOrd<AmbiNuc>> PartialOrd<&[T]> for PackedAmbiDna {
    fn partial_cmp(&self, other: &&[T]) -> Option<Ordering> {
        self.iter().map(RefCmp).partial_cmp(*other)
    }
}

impl PartialEq<&str> for PackedAmbiDna {
    fn eq(&self, rhs: &&str) -> bool {
        self == *rhs
    }
}

impl PartialEq<str> for PackedAmbiDna {
    fn eq(&self, rhs: &str) -> bool {
        self.iter().map(Ok).eq(iter_symbols(rhs))
    }
}

impl std::fmt::Display for PackedAmbiDna {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        self.iter().fmt(f)
    }
}

impl std::fmt::Debug for PackedAmbiDna {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        f.debug_tuple("PackedAmbiDna")
            .field(&display(self.iter()))
            .finish()
    }
}

/// Like [`[AmbiNuc; N]`](array), but takes 50% less space.
///
/// # Examples
///
/// ```
/// use nucs::{AmbiNuc, Packed};
///
/// let dna = AmbiNuc::arr(b"ACATATTACK");
/// assert_eq!(size_of_val(&dna), 10);
/// let packed: Packed<[AmbiNuc; 10]> = dna.into();
/// assert_eq!(size_of_val(&packed), 5);
/// let unpacked: [AmbiNuc; 10] = packed.into();
/// assert_eq!(unpacked, dna);
/// ```
#[derive(Clone, Copy, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct PackedArrayAmbiDna<const N: usize>(PackedBuf<N>)
where
    [(); N]: PackableArray;

type PackedBuf<const N: usize> = <[(); N] as ArrayDivide>::By2<u8>;

impl<const N: usize> PackedArrayAmbiDna<N>
where
    [(); N]: PackableArray,
{
    /// Returns an iterator over the packed DNA.
    #[must_use]
    pub fn iter(&self) -> PackedAmbiDnaIter<'_> {
        self.into_iter()
    }

    /// Returns a mutable iterator over the packed DNA.
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::AmbiNuc;
    ///
    /// let dna = AmbiNuc::seq(b"VTYNGTA");
    /// let mut packed = dna.pack();
    /// for mut nuc in &mut packed {
    ///     *nuc = nuc.complement();
    /// }
    /// assert_eq!(packed, "BARNCAT");
    /// ```
    #[must_use]
    pub fn iter_mut(&mut self) -> PackedAmbiDnaMutIter<'_> {
        self.into_iter()
    }
}

impl<const N: usize> Default for PackedArrayAmbiDna<N>
where
    [(); N]: PackableArray,
{
    fn default() -> Self {
        Self(ArrayDefault::array_default())
    }
}

impl<const N: usize> From<Seq<[AmbiNuc; N]>> for PackedArrayAmbiDna<N>
where
    [(); N]: PackableArray,
{
    fn from(ambi_dna: Seq<[AmbiNuc; N]>) -> PackedArrayAmbiDna<N> {
        ambi_dna.0.into()
    }
}

impl<const N: usize> From<&Seq<[AmbiNuc; N]>> for PackedArrayAmbiDna<N>
where
    [(); N]: PackableArray,
{
    fn from(ambi_dna: &Seq<[AmbiNuc; N]>) -> PackedArrayAmbiDna<N> {
        ambi_dna.0.into()
    }
}

impl<const N: usize> From<[AmbiNuc; N]> for PackedArrayAmbiDna<N>
where
    [(); N]: PackableArray,
{
    fn from(ambi_dna: [AmbiNuc; N]) -> PackedArrayAmbiDna<N> {
        (&ambi_dna).into()
    }
}

impl<const N: usize> From<&[AmbiNuc; N]> for PackedArrayAmbiDna<N>
where
    [(); N]: PackableArray,
{
    fn from(ambi_dna: &[AmbiNuc; N]) -> PackedArrayAmbiDna<N> {
        let mut this = Self(ArrayDefault::array_default());
        pack(this.0.as_mut(), ambi_dna);
        this
    }
}

impl<const N: usize> From<PackedArrayAmbiDna<N>> for Seq<[AmbiNuc; N]>
where
    [(); N]: PackableArray,
{
    fn from(packed_ambi_dna: PackedArrayAmbiDna<N>) -> Seq<[AmbiNuc; N]> {
        Seq(packed_ambi_dna.into())
    }
}

impl<const N: usize> From<&PackedArrayAmbiDna<N>> for Seq<[AmbiNuc; N]>
where
    [(); N]: PackableArray,
{
    fn from(packed_ambi_dna: &PackedArrayAmbiDna<N>) -> Seq<[AmbiNuc; N]> {
        Seq(packed_ambi_dna.into())
    }
}

impl<const N: usize> From<PackedArrayAmbiDna<N>> for [AmbiNuc; N]
where
    [(); N]: PackableArray,
{
    fn from(packed_ambi_dna: PackedArrayAmbiDna<N>) -> [AmbiNuc; N] {
        (&packed_ambi_dna).into()
    }
}

impl<const N: usize> From<&PackedArrayAmbiDna<N>> for [AmbiNuc; N]
where
    [(); N]: PackableArray,
{
    fn from(packed_ambi_dna: &PackedArrayAmbiDna<N>) -> [AmbiNuc; N] {
        let mut ambi_dna = [AmbiNuc::default(); N];
        unpack(&mut ambi_dna, packed_ambi_dna.0.as_ref());
        ambi_dna
    }
}

impl<'a, const N: usize> IntoIterator for &'a PackedArrayAmbiDna<N>
where
    [(); N]: PackableArray,
{
    type Item = AmbiNuc;
    type IntoIter = PackedAmbiDnaIter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        PackedAmbiDnaIter(UnpackingIter::new(0..N, self.0.as_ref()))
    }
}

impl<'a, const N: usize> IntoIterator for &'a mut PackedArrayAmbiDna<N>
where
    [(); N]: PackableArray,
{
    type Item = PackedAmbiNuc<'a>;
    type IntoIter = PackedAmbiDnaMutIter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        let backing = Cell::from_mut(self.0.as_mut()).as_slice_of_cells();
        PackedAmbiDnaMutIter(UnpackingIter::new(0..N, backing))
    }
}

impl<const N: usize> IntoIterator for PackedArrayAmbiDna<N>
where
    [(); N]: PackableArray,
{
    type Item = AmbiNuc;
    type IntoIter = PackedArrayAmbiDnaIntoIter<N>;

    fn into_iter(self) -> Self::IntoIter {
        PackedArrayAmbiDnaIntoIter(UnpackingIter::new(0..N, self.0))
    }
}

impl<const N: usize> PartialEq<PackedAmbiDna> for PackedArrayAmbiDna<N>
where
    [(); N]: PackableArray,
{
    fn eq(&self, other: &PackedAmbiDna) -> bool {
        other == self
    }
}

impl<const N: usize> PartialOrd<PackedAmbiDna> for PackedArrayAmbiDna<N>
where
    [(); N]: PackableArray,
{
    fn partial_cmp(&self, other: &PackedAmbiDna) -> Option<Ordering> {
        other.partial_cmp(self).map(Ordering::reverse)
    }
}

impl<const N: usize> PartialEq<PackedDna> for PackedArrayAmbiDna<N>
where
    [(); N]: PackableArray,
{
    fn eq(&self, other: &PackedDna) -> bool {
        other == self
    }
}

impl<const N: usize> PartialOrd<PackedDna> for PackedArrayAmbiDna<N>
where
    [(); N]: PackableArray,
{
    fn partial_cmp(&self, other: &PackedDna) -> Option<Ordering> {
        other.partial_cmp(self).map(Ordering::reverse)
    }
}

impl<const N: usize> PartialEq<PackedArrayDna<N>> for PackedArrayAmbiDna<N>
where
    [(); N]: PackableArray,
{
    fn eq(&self, other: &PackedArrayDna<N>) -> bool {
        other == self
    }
}

impl<const N: usize> PartialOrd<PackedArrayDna<N>> for PackedArrayAmbiDna<N>
where
    [(); N]: PackableArray,
{
    fn partial_cmp(&self, other: &PackedArrayDna<N>) -> Option<Ordering> {
        other.partial_cmp(self).map(Ordering::reverse)
    }
}

impl<T: PartialEq<AmbiNuc>, const N: usize, const M: usize> PartialEq<[T; M]>
    for PackedArrayAmbiDna<N>
where
    [(); N]: PackableArray,
{
    fn eq(&self, other: &[T; M]) -> bool {
        self.eq(&other.as_slice())
    }
}

impl<T: PartialOrd<AmbiNuc>, const N: usize, const M: usize> PartialOrd<[T; M]>
    for PackedArrayAmbiDna<N>
where
    [(); N]: PackableArray,
{
    fn partial_cmp(&self, other: &[T; M]) -> Option<Ordering> {
        self.partial_cmp(&other.as_slice())
    }
}

impl<T: PartialEq<AmbiNuc>, const N: usize, const M: usize> PartialEq<&mut [T; M]>
    for PackedArrayAmbiDna<N>
where
    [(); N]: PackableArray,
{
    fn eq(&self, other: &&mut [T; M]) -> bool {
        self.eq(&other.as_slice())
    }
}

impl<T: PartialOrd<AmbiNuc>, const N: usize, const M: usize> PartialOrd<&mut [T; M]>
    for PackedArrayAmbiDna<N>
where
    [(); N]: PackableArray,
{
    fn partial_cmp(&self, other: &&mut [T; M]) -> Option<Ordering> {
        self.partial_cmp(&other.as_slice())
    }
}

impl<T: PartialEq<AmbiNuc>, const N: usize, const M: usize> PartialEq<&[T; M]>
    for PackedArrayAmbiDna<N>
where
    [(); N]: PackableArray,
{
    fn eq(&self, other: &&[T; M]) -> bool {
        self.eq(&other.as_slice())
    }
}

impl<T: PartialOrd<AmbiNuc>, const N: usize, const M: usize> PartialOrd<&[T; M]>
    for PackedArrayAmbiDna<N>
where
    [(); N]: PackableArray,
{
    fn partial_cmp(&self, other: &&[T; M]) -> Option<Ordering> {
        self.partial_cmp(&other.as_slice())
    }
}

impl<T: PartialEq<AmbiNuc>, const N: usize> PartialEq<Vec<T>> for PackedArrayAmbiDna<N>
where
    [(); N]: PackableArray,
{
    fn eq(&self, other: &Vec<T>) -> bool {
        self.eq(&&**other)
    }
}

impl<T: PartialOrd<AmbiNuc>, const N: usize> PartialOrd<Vec<T>> for PackedArrayAmbiDna<N>
where
    [(); N]: PackableArray,
{
    fn partial_cmp(&self, other: &Vec<T>) -> Option<Ordering> {
        self.partial_cmp(&&**other)
    }
}

impl<T: PartialEq<AmbiNuc>, const N: usize> PartialEq<&mut [T]> for PackedArrayAmbiDna<N>
where
    [(); N]: PackableArray,
{
    fn eq(&self, other: &&mut [T]) -> bool {
        self.eq(&&**other)
    }
}

impl<T: PartialOrd<AmbiNuc>, const N: usize> PartialOrd<&mut [T]> for PackedArrayAmbiDna<N>
where
    [(); N]: PackableArray,
{
    fn partial_cmp(&self, other: &&mut [T]) -> Option<Ordering> {
        self.partial_cmp(&&**other)
    }
}

impl<T: PartialEq<AmbiNuc>, const N: usize> PartialEq<&[T]> for PackedArrayAmbiDna<N>
where
    [(); N]: PackableArray,
{
    fn eq(&self, other: &&[T]) -> bool {
        // Without `TrustedLen` we don't benefit from stdlib's length-check specializations.
        N == other.len() && self.iter().map(RefCmp).eq(*other)
    }
}

impl<T: PartialOrd<AmbiNuc>, const N: usize> PartialOrd<&[T]> for PackedArrayAmbiDna<N>
where
    [(); N]: PackableArray,
{
    fn partial_cmp(&self, other: &&[T]) -> Option<Ordering> {
        self.iter().map(RefCmp).partial_cmp(*other)
    }
}

impl<const N: usize> PartialEq<&str> for PackedArrayAmbiDna<N>
where
    [(); N]: PackableArray,
{
    fn eq(&self, rhs: &&str) -> bool {
        self == *rhs
    }
}

impl<const N: usize> PartialEq<str> for PackedArrayAmbiDna<N>
where
    [(); N]: PackableArray,
{
    fn eq(&self, rhs: &str) -> bool {
        self.iter().map(Ok).eq(iter_symbols(rhs))
    }
}

impl<const N: usize> std::fmt::Display for PackedArrayAmbiDna<N>
where
    [(); N]: PackableArray,
{
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        self.iter().fmt(f)
    }
}

impl<const N: usize> std::fmt::Debug for PackedArrayAmbiDna<N>
where
    [(); N]: PackableArray,
{
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        f.debug_tuple("PackedArrayAmbiDna")
            .field(&display(self.iter()))
            .finish()
    }
}

fn pack(packed: &mut [u8], ambi_dna: &[AmbiNuc]) {
    let (pairs, remainder) = ambi_dna.as_chunks();
    for (byte, &[n1, n2]) in packed.iter_mut().zip(pairs) {
        *byte = n2.compress() | (n1.compress() << 4);
    }
    if let Some(byte) = packed.last_mut()
        && let [ambi_nuc] = remainder
    {
        *byte = ambi_nuc.compress() << 4;
    }
}

fn unpack(ambi_dna: &mut [AmbiNuc], packed: &[u8]) {
    let (pairs, remainder) = ambi_dna.as_chunks_mut();
    for ([n1, n2], &byte) in pairs.iter_mut().zip(packed) {
        *n1 = AmbiNuc::decompress(byte >> 4);
        *n2 = AmbiNuc::decompress(byte);
    }
    if let Some(byte) = packed.last()
        && let [ambi_nuc] = remainder
    {
        *ambi_nuc = AmbiNuc::decompress(byte >> 4);
    }
}

/// Owned [`PackedAmbiDna`] iterator.
#[derive(Clone)]
pub struct PackedAmbiDnaIntoIter(UnpackingIter<4, 8, std::vec::IntoIter<u8>>);

impl PackedAmbiDnaIntoIter {
    fn as_ref(&self) -> PackedAmbiDnaIter<'_> {
        PackedAmbiDnaIter(self.0.as_ref())
    }
}

impl Iterator for PackedAmbiDnaIntoIter {
    type Item = AmbiNuc;

    fn next(&mut self) -> Option<AmbiNuc> {
        self.0
            .next()
            .map(|(shift, byte)| AmbiNuc::decompress(byte >> shift))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.0.size_hint()
    }
}

impl DoubleEndedIterator for PackedAmbiDnaIntoIter {
    fn next_back(&mut self) -> Option<AmbiNuc> {
        self.0
            .next_back()
            .map(|(shift, byte)| AmbiNuc::decompress(byte >> shift))
    }
}

impl ExactSizeIterator for PackedAmbiDnaIntoIter {
    fn len(&self) -> usize {
        self.0.len()
    }
}

impl std::fmt::Display for PackedAmbiDnaIntoIter {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        self.as_ref().fmt(f)
    }
}

impl std::fmt::Debug for PackedAmbiDnaIntoIter {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        f.debug_tuple("PackedAmbiDnaIntoIter")
            .field(&display(self.as_ref()))
            .finish()
    }
}

/// Owned [`PackedArrayAmbiDna`] iterator.
#[derive(Clone)]
pub struct PackedArrayAmbiDnaIntoIter<const N: usize>(
    UnpackingIter<4, 8, <PackedBuf<N> as IntoIterator>::IntoIter>,
)
where
    [(); N]: PackableArray;

impl<const N: usize> PackedArrayAmbiDnaIntoIter<N>
where
    [(); N]: PackableArray,
{
    fn as_ref(&self) -> PackedAmbiDnaIter<'_> {
        PackedAmbiDnaIter(self.0.as_ref())
    }
}

impl<const N: usize> Iterator for PackedArrayAmbiDnaIntoIter<N>
where
    [(); N]: PackableArray,
{
    type Item = AmbiNuc;

    fn next(&mut self) -> Option<AmbiNuc> {
        self.0
            .next()
            .map(|(shift, byte)| AmbiNuc::decompress(byte >> shift))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.0.size_hint()
    }
}

impl<const N: usize> DoubleEndedIterator for PackedArrayAmbiDnaIntoIter<N>
where
    [(); N]: PackableArray,
{
    fn next_back(&mut self) -> Option<AmbiNuc> {
        self.0
            .next_back()
            .map(|(shift, byte)| AmbiNuc::decompress(byte >> shift))
    }
}

impl<const N: usize> ExactSizeIterator for PackedArrayAmbiDnaIntoIter<N>
where
    [(); N]: PackableArray,
{
    fn len(&self) -> usize {
        self.0.len()
    }
}

impl<const N: usize> std::fmt::Display for PackedArrayAmbiDnaIntoIter<N>
where
    [(); N]: PackableArray,
{
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        self.as_ref().fmt(f)
    }
}

impl<const N: usize> std::fmt::Debug for PackedArrayAmbiDnaIntoIter<N>
where
    [(); N]: PackableArray,
{
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        f.debug_tuple("PackedArrayAmbiDnaIntoIter")
            .field(&display(self.as_ref()))
            .finish()
    }
}

/// Borrowed packed ambiguous DNA iterator.
#[derive(Clone)]
pub struct PackedAmbiDnaIter<'a>(UnpackingIter<4, 8, std::slice::Iter<'a, u8>>);

impl Iterator for PackedAmbiDnaIter<'_> {
    type Item = AmbiNuc;

    fn next(&mut self) -> Option<AmbiNuc> {
        self.0
            .next()
            .map(|(shift, byte)| AmbiNuc::decompress(byte >> shift))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.0.size_hint()
    }
}

impl DoubleEndedIterator for PackedAmbiDnaIter<'_> {
    fn next_back(&mut self) -> Option<AmbiNuc> {
        self.0
            .next_back()
            .map(|(shift, byte)| AmbiNuc::decompress(byte >> shift))
    }
}

impl ExactSizeIterator for PackedAmbiDnaIter<'_> {
    fn len(&self) -> usize {
        self.0.len()
    }
}

impl std::fmt::Display for PackedAmbiDnaIter<'_> {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        display(self.clone()).fmt(f)
    }
}

impl std::fmt::Debug for PackedAmbiDnaIter<'_> {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        f.debug_tuple("PackedAmbiDnaIter")
            .field(&display(self.clone()))
            .finish()
    }
}

/// Mutably borrowed packed ambiguous DNA iterator.
///
/// Created via [`PackedAmbiDna::iter_mut`] and [`PackedArrayAmbiDna::iter_mut`].
pub struct PackedAmbiDnaMutIter<'a>(UnpackingIter<4, 8, std::slice::Iter<'a, Cell<u8>>>);

impl PackedAmbiDnaMutIter<'_> {
    // Create temporary iterator over remaining values without consuming this iter.
    fn read_values(&self) -> impl Iterator<Item = AmbiNuc> + Clone {
        self.0
            .clone()
            .map(|(shift, byte)| AmbiNuc::decompress(byte.get() >> shift))
    }
}

impl<'a> Iterator for PackedAmbiDnaMutIter<'a> {
    type Item = PackedAmbiNuc<'a>;

    fn next(&mut self) -> Option<PackedAmbiNuc<'a>> {
        self.0.next().map(PackedAmbiNuc::new)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.0.size_hint()
    }
}

impl<'a> DoubleEndedIterator for PackedAmbiDnaMutIter<'a> {
    fn next_back(&mut self) -> Option<PackedAmbiNuc<'a>> {
        self.0.next_back().map(PackedAmbiNuc::new)
    }
}

impl ExactSizeIterator for PackedAmbiDnaMutIter<'_> {
    fn len(&self) -> usize {
        self.0.len()
    }
}

impl std::fmt::Display for PackedAmbiDnaMutIter<'_> {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        display(self.read_values()).fmt(f)
    }
}

impl std::fmt::Debug for PackedAmbiDnaMutIter<'_> {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        f.debug_tuple("PackedAmbiDnaMutIter")
            .field(&display(self.read_values()))
            .finish()
    }
}

/// Mutable reference to packed [`AmbiNuc`].
///
/// Yielded by [`PackedAmbiDnaMutIter`].
///
/// # Leaking
///
/// If the [`PackedAmbiNuc`] goes out of scope without being dropped (due to [`std::mem::forget`],
/// for example), changes to the [`AmbiNuc`] may fail to be persisted.
pub struct PackedAmbiNuc<'a> {
    packed: &'a Cell<u8>,
    shift: u8,
    unpacked: AmbiNuc,
}

impl<'a> PackedAmbiNuc<'a> {
    fn new((shift, packed): (u8, &'a Cell<u8>)) -> Self {
        let unpacked = AmbiNuc::decompress(packed.get() >> shift);
        Self {
            packed,
            shift,
            unpacked,
        }
    }
}

impl Drop for PackedAmbiNuc<'_> {
    fn drop(&mut self) {
        self.packed
            .update(|p| p & !(0b1111 << self.shift) | self.unpacked.compress() << self.shift);
    }
}

impl Deref for PackedAmbiNuc<'_> {
    type Target = AmbiNuc;

    fn deref(&self) -> &AmbiNuc {
        &self.unpacked
    }
}

impl DerefMut for PackedAmbiNuc<'_> {
    fn deref_mut(&mut self) -> &mut AmbiNuc {
        &mut self.unpacked
    }
}

impl AsRef<AmbiNuc> for PackedAmbiNuc<'_> {
    fn as_ref(&self) -> &AmbiNuc {
        self
    }
}

impl AsMut<AmbiNuc> for PackedAmbiNuc<'_> {
    fn as_mut(&mut self) -> &mut AmbiNuc {
        self
    }
}

impl std::fmt::Display for PackedAmbiNuc<'_> {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        self.unpacked.fmt(f)
    }
}

impl std::fmt::Debug for PackedAmbiNuc<'_> {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        f.debug_tuple("PackedAmbiNuc")
            .field(&self.unpacked)
            .finish()
    }
}

#[cfg(test)]
mod tests {
    use proptest::{arbitrary::any, proptest};

    use super::super::tests::{assert_both_roundtrips, assert_roundtrip};
    use crate::proptest::{any_ambi_dna, any_dna};

    use super::*;

    #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
    #[test]
    fn all_short_roundtrips() {
        // Yes, this is hideously ugly, but arrays need a size known at compile-time,
        // and this is the simplest way to ensure I get them all.
        assert_both_roundtrips(&[] as &[AmbiNuc; 0]);
        for n1 in AmbiNuc::ALL {
            assert_both_roundtrips(&[n1]);
            for n2 in AmbiNuc::ALL {
                assert_both_roundtrips(&[n1, n2]);
                for n3 in AmbiNuc::ALL {
                    assert_both_roundtrips(&[n1, n2, n3]);
                    for n4 in AmbiNuc::ALL {
                        assert_both_roundtrips(&[n1, n2, n3, n4]);
                    }
                }
            }
        }
    }

    // The bulk of the fiddly logic is handled by UnpackingIter, which has more tests.
    // This is just sanity-checking that things were hooked up to it correctly.
    #[test]
    fn smoke_test_iters() {
        let ambi_dna = AmbiNuc::seq(b"AMBACGT")[..].pack();
        assert_eq!(Seq(Vec::from_iter(&ambi_dna)), "AMBACGT");
        assert_eq!(Seq(Vec::from_iter(ambi_dna)), "AMBACGT");
        let ambi_dna = AmbiNuc::seq(b"AMBACGT").pack();
        assert_eq!(Seq(Vec::from_iter(&ambi_dna)), "AMBACGT");
        assert_eq!(Seq(Vec::from_iter(ambi_dna)), "AMBACGT");
    }

    #[test]
    fn packed_ambi_dna_eq_str() {
        let dna = AmbiNuc::seq(b"AMBACGT")[..].pack();
        assert_eq!(dna, "AMBACGT");
        let dna = AmbiNuc::seq(b"AMBACGT").pack();
        assert_eq!(dna, "AMBACGT");
    }

    #[test]
    fn display() {
        let ambi_dna = AmbiNuc::seq(b"AMBACGT")[..].pack();
        assert_eq!(ambi_dna.to_string(), "AMBACGT");
        assert_eq!(ambi_dna.iter().to_string(), "AMBACGT");
        assert_eq!(ambi_dna.into_iter().to_string(), "AMBACGT");
        let ambi_dna = AmbiNuc::seq(b"AMBACGT").pack();
        assert_eq!(ambi_dna.to_string(), "AMBACGT");
        assert_eq!(ambi_dna.iter().to_string(), "AMBACGT");
        assert_eq!(ambi_dna.into_iter().to_string(), "AMBACGT");
    }

    proptest! {
        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn packed_ambi_dna_roundtrip(
            ambi_dna in any_ambi_dna(5..25) // 0..=4 covered above
        ) {
            assert_roundtrip(&*ambi_dna);
        }

        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn packed_array_ambi_dna5_roundtrip(
            ambi_dna in any::<[AmbiNuc; 5]>()
        ) {
            assert_roundtrip(&ambi_dna);
        }

        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn packed_array_ambi_dna6_roundtrip(
            ambi_dna in any::<[AmbiNuc; 6]>()
        ) {
            assert_roundtrip(&ambi_dna);
        }

        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn packed_ambi_dna_lexical_ordering(
            dna1 in any_ambi_dna(0..50),
            dna2 in any_ambi_dna(0..50),
        ) {
            let packed1 = PackedAmbiDna::from(dna1.as_slice());
            let packed2 = PackedAmbiDna::from(dna2.as_slice());
            assert_eq!(packed1.0.cmp(&packed2.0), dna1.cmp(&dna2));
        }

        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn packed_ambi_dna_length(
            ambi_dna in any_ambi_dna(0..50)
        ) {
            let packed = PackedAmbiDna::from(&ambi_dna);
            assert_eq!(packed.len(), ambi_dna.len());
            assert_eq!(packed.is_empty(), ambi_dna.is_empty());
        }

        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn packed_ambi_dna_ord_vs_slice(
            ambi_dna1 in any_ambi_dna(0..50),
            ambi_dna2 in any_ambi_dna(0..50),
        ) {
            let packed1 = PackedAmbiDna::from(&ambi_dna1);
            assert_eq!(packed1.partial_cmp(&ambi_dna2), ambi_dna1.partial_cmp(&ambi_dna2));
            assert_eq!(packed1.eq(&ambi_dna2), ambi_dna1.eq(&ambi_dna2));
        }

        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn packed_ambi_dna_ord_vs_concrete_slice(
            ambi_dna1 in any_ambi_dna(0..50),
            dna2 in any_dna(0..50),
        ) {
            let packed1 = PackedAmbiDna::from(&ambi_dna1);
            assert_eq!(packed1.partial_cmp(&dna2), ambi_dna1.iter().partial_cmp(&dna2));
            assert_eq!(packed1.eq(&dna2), ambi_dna1.eq(&dna2));
        }

        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn packed_ambi_dna5_ord(
            ambi_dna1 in any::<[AmbiNuc; 5]>(),
            ambi_dna2 in any_ambi_dna(0..10),
        ) {
            let packed1 = PackedArrayAmbiDna::from(&ambi_dna1);
            assert_eq!(packed1.partial_cmp(&ambi_dna2), ambi_dna1.as_slice().partial_cmp(&ambi_dna2));
            assert_eq!(packed1.eq(&ambi_dna2), ambi_dna1.as_slice().eq(&ambi_dna2));
        }

        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn push_and_pop(
            mut ambi_dna in any_ambi_dna(0..25),
            ops in proptest::collection::vec(any::<Option<AmbiNuc>>(), 0..50),
        ) {
            let mut packed = PackedAmbiDna::from(&ambi_dna);
            for op in ops {
                if let Some(ambi_nuc) = op {
                    packed.push(ambi_nuc);
                    ambi_dna.push(ambi_nuc);
                } else {
                    assert_eq!(packed.pop(), ambi_dna.pop());
                }
                assert_eq!(packed, ambi_dna);
            }
        }
    }
}
