//! Types related to packed ambiguous DNA.

use std::cmp::Ordering;
use std::fmt::Formatter;

use crate::{AmbiNuc, Seq, iter::display};

use super::packable_array::{ArrayDefault, Sealed as ArrayDivide};
use super::{PackableArray, UnpackingIter};

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
#[derive(Clone)]
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

    /// Returns an iterator over the packed DNA.
    #[must_use]
    pub fn iter(&self) -> PackedAmbiDnaIter<'_> {
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

impl IntoIterator for PackedAmbiDna {
    type Item = AmbiNuc;
    type IntoIter = PackedAmbiDnaIntoIter;

    fn into_iter(self) -> Self::IntoIter {
        PackedAmbiDnaIntoIter(UnpackingIter::new(0..self.len(), self.0))
    }
}

impl PartialEq for PackedAmbiDna {
    fn eq(&self, other: &Self) -> bool {
        self.0.eq(&other.0)
    }
}

impl Eq for PackedAmbiDna {}

impl PartialOrd for PackedAmbiDna {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for PackedAmbiDna {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.cmp(&other.0)
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
#[derive(Clone, Copy)]
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

impl<const N: usize> PartialEq for PackedArrayAmbiDna<N>
where
    [(); N]: PackableArray,
{
    fn eq(&self, other: &Self) -> bool {
        self.0.as_ref().eq(other.0.as_ref())
    }
}

impl<const N: usize> Eq for PackedArrayAmbiDna<N> where [(); N]: PackableArray {}

impl<const N: usize> PartialOrd for PackedArrayAmbiDna<N>
where
    [(); N]: PackableArray,
{
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<const N: usize> Ord for PackedArrayAmbiDna<N>
where
    [(); N]: PackableArray,
{
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.as_ref().cmp(other.0.as_ref())
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

#[cfg(test)]
mod tests {
    use proptest::{arbitrary::any, proptest};

    use super::super::tests::{assert_both_roundtrips, assert_roundtrip};
    use crate::proptest::any_ambi_dna;

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
        fn packed_ambi_dna_ord(
            ambi_dna1 in any_ambi_dna(0..50),
            ambi_dna2 in any_ambi_dna(0..50),
        ) {
            let packed1 = PackedAmbiDna::from(&ambi_dna1);
            let packed2 = PackedAmbiDna::from(&ambi_dna2);
            assert_eq!(packed1.cmp(&packed2), ambi_dna1.cmp(&ambi_dna2));
            assert_eq!(packed1.eq(&packed2), ambi_dna1.eq(&ambi_dna2));
        }

        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn packed_array_dna1_ord(
            ambi_dna1 in any::<[AmbiNuc; 1]>(),
            ambi_dna2 in any::<[AmbiNuc; 1]>(),
        ) {
            let packed1 = PackedArrayAmbiDna::from(&ambi_dna1);
            let packed2 = PackedArrayAmbiDna::from(&ambi_dna2);
            assert_eq!(packed1.cmp(&packed2), ambi_dna1.cmp(&ambi_dna2));
            assert_eq!(packed1.eq(&packed2), ambi_dna1.eq(&ambi_dna2));
        }

        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn packed_array_dna2_ord(
            ambi_dna1 in any::<[AmbiNuc; 2]>(),
            ambi_dna2 in any::<[AmbiNuc; 2]>(),
        ) {
            let packed1 = PackedArrayAmbiDna::from(&ambi_dna1);
            let packed2 = PackedArrayAmbiDna::from(&ambi_dna2);
            assert_eq!(packed1.cmp(&packed2), ambi_dna1.cmp(&ambi_dna2));
            assert_eq!(packed1.eq(&packed2), ambi_dna1.eq(&ambi_dna2));
        }

        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn packed_array_dna3_ord(
            ambi_dna1 in any::<[AmbiNuc; 3]>(),
            ambi_dna2 in any::<[AmbiNuc; 3]>(),
        ) {
            let packed1 = PackedArrayAmbiDna::from(&ambi_dna1);
            let packed2 = PackedArrayAmbiDna::from(&ambi_dna2);
            assert_eq!(packed1.cmp(&packed2), ambi_dna1.cmp(&ambi_dna2));
            assert_eq!(packed1.eq(&packed2), ambi_dna1.eq(&ambi_dna2));
        }
    }
}
