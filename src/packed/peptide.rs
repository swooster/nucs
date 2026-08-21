//! Types related to packed ambiguous peptides.

use std::cmp::Ordering;
use std::fmt::Formatter;

use crate::{Amino, Seq, iter::display};

use super::packable_array::{ArrayDefault, Sealed as ArrayDivide};
use super::{PackableArray, RefCmp, UnpackingIter};

// Note on storage: Aminos are packed big-endian so naive lexical sorting of bytes is correct.
// The [u8; 2] themselves are big-endian u16s, for the same reason.
//
// The last word of a non-empty `PackedPeptide` must have 1-3 elements. If it has 2 elements,
// the 5 least-significant bits are 0. If it has 1 element, the 10 least-significant bits are 0.

/// Like [`Vec<Amino>`], but takes 33% less space for long peptides.
///
/// # Examples
///
/// ```
/// use nucs::{Amino, Packed};
///
/// let peptide = Amino::arr(b"KITTYTUMSOFDANGER");
/// let packed: Packed<[Amino]> = peptide.as_slice().into();
/// let unpacked: Vec<Amino> = packed.into();
/// assert_eq!(unpacked, peptide);
/// ```
#[derive(Clone, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct PackedPeptide(Vec<[u8; 2]>);

impl PackedPeptide {
    /// Returns the number of [`Amino`]s in the packed peptide.
    #[must_use]
    pub fn len(&self) -> usize {
        // len <= `usize::MAX` is an invariant of this type, so no need to worry about overflow.
        match &*self.0 {
            [] => 0,
            bulk @ [.., tail] => {
                3 * bulk.len() - (u16::from_be_bytes(*tail).trailing_zeros() / 5) as usize
            }
        }
    }

    /// Returns `true` if the packed peptide contains no [`Amino`]s.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns an iterator over the packed peptide.
    #[must_use]
    pub fn iter(&self) -> PackedPeptideIter<'_> {
        self.into_iter()
    }
}

impl From<Seq<Vec<Amino>>> for PackedPeptide {
    fn from(peptide: Seq<Vec<Amino>>) -> PackedPeptide {
        peptide.0.into()
    }
}

impl From<Vec<Amino>> for PackedPeptide {
    fn from(peptide: Vec<Amino>) -> PackedPeptide {
        (&peptide).into()
    }
}

impl<T: AsRef<[Amino]> + ?Sized> From<&T> for PackedPeptide {
    fn from(peptide: &T) -> PackedPeptide {
        let peptide = peptide.as_ref();
        let packed_len = peptide.len().div_ceil(3);
        let mut packed = vec![[0, 0]; packed_len];
        pack(&mut packed, peptide);
        Self(packed)
    }
}

impl From<PackedPeptide> for Seq<Vec<Amino>> {
    fn from(packed_peptide: PackedPeptide) -> Seq<Vec<Amino>> {
        Seq(packed_peptide.into())
    }
}

impl From<&PackedPeptide> for Seq<Vec<Amino>> {
    fn from(packed_peptide: &PackedPeptide) -> Seq<Vec<Amino>> {
        Seq(packed_peptide.into())
    }
}

impl From<PackedPeptide> for Vec<Amino> {
    fn from(packed_peptide: PackedPeptide) -> Vec<Amino> {
        (&packed_peptide).into()
    }
}

impl From<&PackedPeptide> for Vec<Amino> {
    fn from(packed_peptide: &PackedPeptide) -> Vec<Amino> {
        let mut peptide = vec![Amino::default(); packed_peptide.len()];
        unpack(&mut peptide, &packed_peptide.0);
        peptide
    }
}

impl<'a> IntoIterator for &'a PackedPeptide {
    type Item = Amino;
    type IntoIter = PackedPeptideIter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        PackedPeptideIter(UnpackingIter::new(0..self.len(), &self.0))
    }
}

impl IntoIterator for PackedPeptide {
    type Item = Amino;
    type IntoIter = PackedPeptideIntoIter;

    fn into_iter(self) -> Self::IntoIter {
        PackedPeptideIntoIter(UnpackingIter::new(0..self.len(), self.0))
    }
}

impl<T: PartialEq<Amino>, const M: usize> PartialEq<[T; M]> for PackedPeptide {
    fn eq(&self, other: &[T; M]) -> bool {
        self.eq(&other.as_slice())
    }
}

impl<T: PartialOrd<Amino>, const M: usize> PartialOrd<[T; M]> for PackedPeptide {
    fn partial_cmp(&self, other: &[T; M]) -> Option<Ordering> {
        self.partial_cmp(&other.as_slice())
    }
}

impl<T: PartialEq<Amino>, const M: usize> PartialEq<&mut [T; M]> for PackedPeptide {
    fn eq(&self, other: &&mut [T; M]) -> bool {
        self.eq(&other.as_slice())
    }
}

impl<T: PartialOrd<Amino>, const M: usize> PartialOrd<&mut [T; M]> for PackedPeptide {
    fn partial_cmp(&self, other: &&mut [T; M]) -> Option<Ordering> {
        self.partial_cmp(&other.as_slice())
    }
}

impl<T: PartialEq<Amino>, const M: usize> PartialEq<&[T; M]> for PackedPeptide {
    fn eq(&self, other: &&[T; M]) -> bool {
        self.eq(&other.as_slice())
    }
}

impl<T: PartialOrd<Amino>, const M: usize> PartialOrd<&[T; M]> for PackedPeptide {
    fn partial_cmp(&self, other: &&[T; M]) -> Option<Ordering> {
        self.partial_cmp(&other.as_slice())
    }
}

impl<T: PartialEq<Amino>> PartialEq<Vec<T>> for PackedPeptide {
    fn eq(&self, other: &Vec<T>) -> bool {
        self.eq(&&**other)
    }
}

impl<T: PartialOrd<Amino>> PartialOrd<Vec<T>> for PackedPeptide {
    fn partial_cmp(&self, other: &Vec<T>) -> Option<Ordering> {
        self.partial_cmp(&&**other)
    }
}

impl<T: PartialEq<Amino>> PartialEq<&mut [T]> for PackedPeptide {
    fn eq(&self, other: &&mut [T]) -> bool {
        self.eq(&&**other)
    }
}

impl<T: PartialOrd<Amino>> PartialOrd<&mut [T]> for PackedPeptide {
    fn partial_cmp(&self, other: &&mut [T]) -> Option<Ordering> {
        self.partial_cmp(&&**other)
    }
}

impl<T: PartialEq<Amino>> PartialEq<&[T]> for PackedPeptide {
    fn eq(&self, other: &&[T]) -> bool {
        // Without `TrustedLen` we don't benefit from stdlib's length-check specializations.
        self.len() == other.len() && self.iter().map(RefCmp).eq(*other)
    }
}

impl<T: PartialOrd<Amino>> PartialOrd<&[T]> for PackedPeptide {
    fn partial_cmp(&self, other: &&[T]) -> Option<Ordering> {
        self.iter().map(RefCmp).partial_cmp(*other)
    }
}

impl std::fmt::Display for PackedPeptide {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        self.iter().fmt(f)
    }
}

impl std::fmt::Debug for PackedPeptide {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        f.debug_tuple("PackedPeptide")
            .field(&display(self.iter()))
            .finish()
    }
}

/// Like [`[Amino; N]`](array), but takes 33% less space.
///
/// # Examples
///
/// ```
/// use nucs::{Amino, Packed};
///
/// let peptide = Amino::arr(b"KITTYTUMSOFDANGER");
/// assert_eq!(size_of_val(&peptide), 17);
/// let packed: Packed<[Amino; 17]> = peptide.into();
/// assert_eq!(size_of_val(&packed), 12);
/// let unpacked: [Amino; 17] = packed.into();
/// assert_eq!(unpacked, peptide);
/// ```
#[derive(Clone, Copy, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct PackedArrayPeptide<const N: usize>(PackedBuf<N>)
where
    [(); N]: PackableArray;

type PackedBuf<const N: usize> = <[(); N] as ArrayDivide>::By3<[u8; 2]>;

impl<const N: usize> PackedArrayPeptide<N>
where
    [(); N]: PackableArray,
{
    /// Returns an iterator over the packed peptide.
    #[must_use]
    pub fn iter(&self) -> PackedPeptideIter<'_> {
        self.into_iter()
    }
}

impl<const N: usize> Default for PackedArrayPeptide<N>
where
    [(); N]: PackableArray,
{
    fn default() -> Self {
        Self(ArrayDefault::array_default())
    }
}

impl<const N: usize> From<Seq<[Amino; N]>> for PackedArrayPeptide<N>
where
    [(); N]: PackableArray,
{
    fn from(peptide: Seq<[Amino; N]>) -> PackedArrayPeptide<N> {
        peptide.0.into()
    }
}

impl<const N: usize> From<&Seq<[Amino; N]>> for PackedArrayPeptide<N>
where
    [(); N]: PackableArray,
{
    fn from(peptide: &Seq<[Amino; N]>) -> PackedArrayPeptide<N> {
        peptide.0.into()
    }
}

impl<const N: usize> From<[Amino; N]> for PackedArrayPeptide<N>
where
    [(); N]: PackableArray,
{
    fn from(peptide: [Amino; N]) -> PackedArrayPeptide<N> {
        (&peptide).into()
    }
}

impl<const N: usize> From<&[Amino; N]> for PackedArrayPeptide<N>
where
    [(); N]: PackableArray,
{
    fn from(peptide: &[Amino; N]) -> PackedArrayPeptide<N> {
        let mut this = Self(ArrayDefault::array_default());
        pack(this.0.as_mut(), peptide);
        this
    }
}

impl<const N: usize> From<PackedArrayPeptide<N>> for Seq<[Amino; N]>
where
    [(); N]: PackableArray,
{
    fn from(packed_peptide: PackedArrayPeptide<N>) -> Seq<[Amino; N]> {
        Seq(packed_peptide.into())
    }
}

impl<const N: usize> From<&PackedArrayPeptide<N>> for Seq<[Amino; N]>
where
    [(); N]: PackableArray,
{
    fn from(packed_peptide: &PackedArrayPeptide<N>) -> Seq<[Amino; N]> {
        Seq(packed_peptide.into())
    }
}

impl<const N: usize> From<PackedArrayPeptide<N>> for [Amino; N]
where
    [(); N]: PackableArray,
{
    fn from(packed_peptide: PackedArrayPeptide<N>) -> [Amino; N] {
        (&packed_peptide).into()
    }
}

impl<const N: usize> From<&PackedArrayPeptide<N>> for [Amino; N]
where
    [(); N]: PackableArray,
{
    fn from(packed_peptide: &PackedArrayPeptide<N>) -> [Amino; N] {
        let mut peptide = [Amino::default(); N];
        unpack(&mut peptide, packed_peptide.0.as_ref());
        peptide
    }
}

impl<'a, const N: usize> IntoIterator for &'a PackedArrayPeptide<N>
where
    [(); N]: PackableArray,
{
    type Item = Amino;
    type IntoIter = PackedPeptideIter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        PackedPeptideIter(UnpackingIter::new(0..N, self.0.as_ref()))
    }
}

impl<const N: usize> IntoIterator for PackedArrayPeptide<N>
where
    [(); N]: PackableArray,
{
    type Item = Amino;
    type IntoIter = PackedArrayPeptideIntoIter<N>;

    fn into_iter(self) -> Self::IntoIter {
        PackedArrayPeptideIntoIter(UnpackingIter::new(0..N, self.0))
    }
}

impl<T: PartialEq<Amino>, const N: usize, const M: usize> PartialEq<[T; M]>
    for PackedArrayPeptide<N>
where
    [(); N]: PackableArray,
{
    fn eq(&self, other: &[T; M]) -> bool {
        self.eq(&other.as_slice())
    }
}

impl<T: PartialOrd<Amino>, const N: usize, const M: usize> PartialOrd<[T; M]>
    for PackedArrayPeptide<N>
where
    [(); N]: PackableArray,
{
    fn partial_cmp(&self, other: &[T; M]) -> Option<Ordering> {
        self.partial_cmp(&other.as_slice())
    }
}

impl<T: PartialEq<Amino>, const N: usize, const M: usize> PartialEq<&mut [T; M]>
    for PackedArrayPeptide<N>
where
    [(); N]: PackableArray,
{
    fn eq(&self, other: &&mut [T; M]) -> bool {
        self.eq(&other.as_slice())
    }
}

impl<T: PartialOrd<Amino>, const N: usize, const M: usize> PartialOrd<&mut [T; M]>
    for PackedArrayPeptide<N>
where
    [(); N]: PackableArray,
{
    fn partial_cmp(&self, other: &&mut [T; M]) -> Option<Ordering> {
        self.partial_cmp(&other.as_slice())
    }
}

impl<T: PartialEq<Amino>, const N: usize, const M: usize> PartialEq<&[T; M]>
    for PackedArrayPeptide<N>
where
    [(); N]: PackableArray,
{
    fn eq(&self, other: &&[T; M]) -> bool {
        self.eq(&other.as_slice())
    }
}

impl<T: PartialOrd<Amino>, const N: usize, const M: usize> PartialOrd<&[T; M]>
    for PackedArrayPeptide<N>
where
    [(); N]: PackableArray,
{
    fn partial_cmp(&self, other: &&[T; M]) -> Option<Ordering> {
        self.partial_cmp(&other.as_slice())
    }
}

impl<T: PartialEq<Amino>, const N: usize> PartialEq<Vec<T>> for PackedArrayPeptide<N>
where
    [(); N]: PackableArray,
{
    fn eq(&self, other: &Vec<T>) -> bool {
        self.eq(&&**other)
    }
}

impl<T: PartialOrd<Amino>, const N: usize> PartialOrd<Vec<T>> for PackedArrayPeptide<N>
where
    [(); N]: PackableArray,
{
    fn partial_cmp(&self, other: &Vec<T>) -> Option<Ordering> {
        self.partial_cmp(&&**other)
    }
}

impl<T: PartialEq<Amino>, const N: usize> PartialEq<&mut [T]> for PackedArrayPeptide<N>
where
    [(); N]: PackableArray,
{
    fn eq(&self, other: &&mut [T]) -> bool {
        self.eq(&&**other)
    }
}

impl<T: PartialOrd<Amino>, const N: usize> PartialOrd<&mut [T]> for PackedArrayPeptide<N>
where
    [(); N]: PackableArray,
{
    fn partial_cmp(&self, other: &&mut [T]) -> Option<Ordering> {
        self.partial_cmp(&&**other)
    }
}

impl<T: PartialEq<Amino>, const N: usize> PartialEq<&[T]> for PackedArrayPeptide<N>
where
    [(); N]: PackableArray,
{
    fn eq(&self, other: &&[T]) -> bool {
        // Without `TrustedLen` we don't benefit from stdlib's length-check specializations.
        N == other.len() && self.iter().map(RefCmp).eq(*other)
    }
}

impl<T: PartialOrd<Amino>, const N: usize> PartialOrd<&[T]> for PackedArrayPeptide<N>
where
    [(); N]: PackableArray,
{
    fn partial_cmp(&self, other: &&[T]) -> Option<Ordering> {
        self.iter().map(RefCmp).partial_cmp(*other)
    }
}

impl<const N: usize> std::fmt::Display for PackedArrayPeptide<N>
where
    [(); N]: PackableArray,
{
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        self.iter().fmt(f)
    }
}

impl<const N: usize> std::fmt::Debug for PackedArrayPeptide<N>
where
    [(); N]: PackableArray,
{
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        f.debug_tuple("PackedArrayPeptide")
            .field(&display(self.iter()))
            .finish()
    }
}

fn pack(packed: &mut [[u8; 2]], peptide: &[Amino]) {
    let (triplets, remainder) = peptide.as_chunks();
    for (pair, &[a1, a2, a3]) in packed.iter_mut().zip(triplets) {
        *pair = (a3.compress() | (a2.compress() << 5) | (a1.compress() << 10)).to_be_bytes();
    }
    match (packed.last_mut(), remainder) {
        (_, []) => {}
        (Some(pair), [a1]) => *pair = (a1.compress() << 10).to_be_bytes(),
        (Some(pair), [a1, a2]) => {
            *pair = ((a1.compress() << 10) | (a2.compress() << 5)).to_be_bytes();
        }
        _ => panic!(),
    }
}

fn unpack(peptide: &mut [Amino], packed: &[[u8; 2]]) {
    let (triplets, remainder) = peptide.as_chunks_mut();
    for ([a1, a2, a3], &pair) in triplets.iter_mut().zip(packed) {
        let val = u16::from_be_bytes(pair);
        *a1 = Amino::decompress(val >> 10);
        *a2 = Amino::decompress(val >> 5);
        *a3 = Amino::decompress(val);
    }
    match (remainder, packed.last()) {
        ([], _) => {}
        ([a1], Some(pair)) => *a1 = Amino::decompress(u16::from_be_bytes(*pair) >> 10),
        ([a1, a2], Some(pair)) => {
            let val = u16::from_be_bytes(*pair);
            *a1 = Amino::decompress(val >> 10);
            *a2 = Amino::decompress(val >> 5);
        }
        _ => panic!(),
    }
}

/// Owned [`PackedPeptide`] iterator.
#[derive(Clone)]
pub struct PackedPeptideIntoIter(UnpackingIter<5, 15, std::vec::IntoIter<[u8; 2]>>);

impl PackedPeptideIntoIter {
    fn as_ref(&self) -> PackedPeptideIter<'_> {
        PackedPeptideIter(self.0.as_ref())
    }
}

impl Iterator for PackedPeptideIntoIter {
    type Item = Amino;

    fn next(&mut self) -> Option<Amino> {
        self.0
            .next()
            .map(|(shift, bytes)| Amino::decompress(u16::from_be_bytes(bytes) >> shift))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.0.size_hint()
    }
}

impl DoubleEndedIterator for PackedPeptideIntoIter {
    fn next_back(&mut self) -> Option<Amino> {
        self.0
            .next_back()
            .map(|(shift, bytes)| Amino::decompress(u16::from_be_bytes(bytes) >> shift))
    }
}

impl ExactSizeIterator for PackedPeptideIntoIter {
    fn len(&self) -> usize {
        self.0.len()
    }
}

impl std::fmt::Display for PackedPeptideIntoIter {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        self.as_ref().fmt(f)
    }
}

impl std::fmt::Debug for PackedPeptideIntoIter {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        f.debug_tuple("PackedPeptideIntoIter")
            .field(&display(self.as_ref()))
            .finish()
    }
}

/// Owned [`PackedArrayPeptide`] iterator.
#[derive(Clone)]
pub struct PackedArrayPeptideIntoIter<const N: usize>(
    UnpackingIter<5, 15, <PackedBuf<N> as IntoIterator>::IntoIter>,
)
where
    [(); N]: PackableArray;

impl<const N: usize> PackedArrayPeptideIntoIter<N>
where
    [(); N]: PackableArray,
{
    fn as_ref(&self) -> PackedPeptideIter<'_> {
        PackedPeptideIter(self.0.as_ref())
    }
}

impl<const N: usize> Iterator for PackedArrayPeptideIntoIter<N>
where
    [(); N]: PackableArray,
{
    type Item = Amino;

    fn next(&mut self) -> Option<Amino> {
        self.0
            .next()
            .map(|(shift, bytes)| Amino::decompress(u16::from_be_bytes(bytes) >> shift))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.0.size_hint()
    }
}

impl<const N: usize> DoubleEndedIterator for PackedArrayPeptideIntoIter<N>
where
    [(); N]: PackableArray,
{
    fn next_back(&mut self) -> Option<Amino> {
        self.0
            .next_back()
            .map(|(shift, bytes)| Amino::decompress(u16::from_be_bytes(bytes) >> shift))
    }
}

impl<const N: usize> ExactSizeIterator for PackedArrayPeptideIntoIter<N>
where
    [(); N]: PackableArray,
{
    fn len(&self) -> usize {
        self.0.len()
    }
}

impl<const N: usize> std::fmt::Display for PackedArrayPeptideIntoIter<N>
where
    [(); N]: PackableArray,
{
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        self.as_ref().fmt(f)
    }
}

impl<const N: usize> std::fmt::Debug for PackedArrayPeptideIntoIter<N>
where
    [(); N]: PackableArray,
{
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        f.debug_tuple("PackedArrayPeptideIntoIter")
            .field(&display(self.as_ref()))
            .finish()
    }
}

/// Borrowed packed peptide iterator.
#[derive(Clone)]
pub struct PackedPeptideIter<'a>(UnpackingIter<5, 15, std::slice::Iter<'a, [u8; 2]>>);

impl Iterator for PackedPeptideIter<'_> {
    type Item = Amino;

    fn next(&mut self) -> Option<Amino> {
        self.0
            .next()
            .map(|(shift, bytes)| Amino::decompress(u16::from_be_bytes(*bytes) >> shift))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.0.size_hint()
    }
}

impl DoubleEndedIterator for PackedPeptideIter<'_> {
    fn next_back(&mut self) -> Option<Amino> {
        self.0
            .next_back()
            .map(|(shift, bytes)| Amino::decompress(u16::from_be_bytes(*bytes) >> shift))
    }
}

impl ExactSizeIterator for PackedPeptideIter<'_> {
    fn len(&self) -> usize {
        self.0.len()
    }
}

impl std::fmt::Display for PackedPeptideIter<'_> {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        display(self.clone()).fmt(f)
    }
}

impl std::fmt::Debug for PackedPeptideIter<'_> {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        f.debug_tuple("PackedPeptideIter")
            .field(&display(self.clone()))
            .finish()
    }
}

#[cfg(test)]
mod tests {
    use proptest::{arbitrary::any, proptest};

    use super::super::tests::{assert_both_roundtrips, assert_roundtrip};
    use crate::proptest::{any_ambi_peptide, any_peptide};

    use super::*;

    #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
    #[test]
    fn all_short_roundtrips() {
        // Yes, this is hideously ugly, but arrays need a size known at compile-time,
        // and this is the simplest way to ensure I get them all.
        assert_both_roundtrips(&[] as &[Amino; 0]);
        for a1 in Amino::ALL {
            assert_both_roundtrips(&[a1]);
            for a2 in Amino::ALL {
                assert_both_roundtrips(&[a1, a2]);
                for a3 in Amino::ALL {
                    assert_both_roundtrips(&[a1, a2, a3]);
                    for a4 in Amino::ALL {
                        assert_both_roundtrips(&[a1, a2, a3, a4]);
                    }
                }
            }
        }
    }

    // The bulk of the fiddly logic is handled by UnpackingIter, which has more tests.
    // This is just sanity-checking that things were hooked up to it correctly.
    #[test]
    fn smoke_test_iters() {
        let peptide = Amino::seq(b"PEPTIDE")[..].pack();
        assert_eq!(Seq(Vec::from_iter(&peptide)), "PEPTIDE");
        assert_eq!(Seq(Vec::from_iter(peptide)), "PEPTIDE");
        let peptide = Amino::seq(b"PEPTIDE").pack();
        assert_eq!(Seq(Vec::from_iter(&peptide)), "PEPTIDE");
        assert_eq!(Seq(Vec::from_iter(peptide)), "PEPTIDE");
    }

    #[test]
    fn display() {
        let peptide = Amino::seq(b"PEPTIDE")[..].pack();
        assert_eq!(peptide.to_string(), "PEPTIDE");
        assert_eq!(peptide.iter().to_string(), "PEPTIDE");
        assert_eq!(peptide.into_iter().to_string(), "PEPTIDE");
        let peptide = Amino::seq(b"PEPTIDE").pack();
        assert_eq!(peptide.to_string(), "PEPTIDE");
        assert_eq!(peptide.iter().to_string(), "PEPTIDE");
        assert_eq!(peptide.into_iter().to_string(), "PEPTIDE");
    }

    proptest! {
        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn packed_peptide_roundtrip(
            peptide in any_peptide(5..25) // 0..=4 covered above
        ) {
            assert_roundtrip(&*peptide);
        }

        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn packed_array_peptide5_roundtrip(
            peptide in any::<[Amino; 5]>()
        ) {
            assert_roundtrip(&peptide);
        }

        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn packed_array_peptide6_roundtrip(
            peptide in any::<[Amino; 6]>()
        ) {
            assert_roundtrip(&peptide);
        }

        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn packed_array_peptide7_roundtrip(
            peptide in any::<[Amino; 7]>()
        ) {
            assert_roundtrip(&peptide);
        }

        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn packed_peptide_lexical_ordering(
            peptide1 in any_peptide(0..50),
            peptide2 in any_peptide(0..50),
        ) {
            let packed1 = PackedPeptide::from(peptide1.as_slice());
            let packed2 = PackedPeptide::from(peptide2.as_slice());
            assert_eq!(packed1.0.cmp(&packed2.0), peptide1.cmp(&peptide2));
        }

        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn packed_peptide_length(
            peptide in any_peptide(0..50)
        ) {
            let packed = PackedPeptide::from(&peptide);
            assert_eq!(packed.len(), peptide.len());
            assert_eq!(packed.is_empty(), peptide.is_empty());
        }

        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn packed_peptide_ord_vs_slice(
            peptide1 in any_peptide(0..50),
            peptide2 in any_peptide(0..50),
        ) {
            let packed1 = PackedPeptide::from(&peptide1);
            assert_eq!(packed1.partial_cmp(&peptide2), peptide1.partial_cmp(&peptide2));
            assert_eq!(packed1.eq(&peptide2), peptide1.eq(&peptide2));
        }

        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn packed_peptide_ord_vs_ambi_slice(
            peptide1 in any_peptide(0..50),
            ambi_peptide2 in any_ambi_peptide(0..50),
        ) {
            let packed1 = PackedPeptide::from(&peptide1);
            assert_eq!(packed1.partial_cmp(&ambi_peptide2), peptide1.iter().partial_cmp(&ambi_peptide2));
            assert_eq!(packed1.eq(&ambi_peptide2), peptide1.eq(&ambi_peptide2));
        }

        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn packed_peptide5_ord(
            peptide1 in any::<[Amino; 5]>(),
            peptide2 in any_peptide(0..10),
        ) {
            let packed1 = PackedArrayPeptide::from(&peptide1);
            assert_eq!(packed1.partial_cmp(&peptide2), peptide1.as_slice().partial_cmp(&peptide2));
            assert_eq!(packed1.eq(&peptide2), peptide1.as_slice().eq(&peptide2));
        }
    }
}
