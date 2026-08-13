//! Types related to packed DNA.

use std::cmp::Ordering;
use std::fmt::Formatter;

use crate::{Nuc, Seq, iter::display};

use super::packable_array::{ArrayDefault, Sealed as ArrayDivide};
use super::{PackableArray, UnpackingIter};

// Note on storage: Nucs are packed big-endian. Sadly, this encoding DOES NOT maintain lexical
// ordering, due to the suffix (see below).
//
// The last byte of a non-empty `PackedDna` must have 0-3 elements. Any missing element is encoded
// 0b00, and the the least significant 2 bits record the number of elements.
//
// We do this rather than store the length externally, both to minimize the memory usage
// (with inline lengths, you don't need any extra memory 75% of the time, whereas with an
// external u8, the Vec overhead would increase from 24 to 32 bytes) and to ensure the
// length can be correctly inferred from serialized data.

/// Like [`Vec<Nuc>`], but takes 75% less space for long DNA sequences.
///
/// # Examples
///
/// ```
/// use nucs::{Nuc, Packed};
///
/// let dna = Nuc::arr(b"ACATATTAC");
/// let packed: Packed<[Nuc]> = dna.as_slice().into();
/// let unpacked: Vec<Nuc> = packed.into();
/// assert_eq!(unpacked, dna);
/// ```
#[derive(Clone)]
pub struct PackedDna(Vec<u8>);

impl PackedDna {
    /// Returns the number of [`Nuc`]s in the packed DNA.
    #[must_use]
    pub fn len(&self) -> usize {
        // len <= `usize::MAX` is an invariant of this type, so no need to worry about overflow.
        match &*self.0 {
            [] => 0,
            [bulk @ .., tail] => 4 * bulk.len() + (tail & 0b11) as usize,
        }
    }

    /// Returns `true` if the packed DNA contains no [`Nuc`]s.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        matches!(&*self.0, [] | [0])
    }

    /// Returns an iterator over the packed DNA.
    #[must_use]
    pub fn iter(&self) -> PackedDnaIter<'_> {
        self.into_iter()
    }
}

impl From<Seq<Vec<Nuc>>> for PackedDna {
    fn from(dna: Seq<Vec<Nuc>>) -> PackedDna {
        dna.0.into()
    }
}

impl From<Vec<Nuc>> for PackedDna {
    fn from(dna: Vec<Nuc>) -> PackedDna {
        (&dna).into()
    }
}

impl<T: AsRef<[Nuc]> + ?Sized> From<&T> for PackedDna {
    fn from(dna: &T) -> PackedDna {
        let dna = dna.as_ref();
        let packed_len = match dna.len() {
            0 => 0, // special-case to reduce allocs
            l => l / 4 + 1,
        };
        let mut packed = vec![0; packed_len];
        pack(&mut packed, dna);
        if let Some(byte) = packed.last_mut() {
            *byte |=
                u8::try_from(dna.len() % 4).unwrap_or_else(|_| unreachable!("x % 4 fits in a u8"));
        }
        Self(packed)
    }
}

impl From<PackedDna> for Seq<Vec<Nuc>> {
    fn from(packed_dna: PackedDna) -> Seq<Vec<Nuc>> {
        Seq(packed_dna.into())
    }
}

impl From<&PackedDna> for Seq<Vec<Nuc>> {
    fn from(packed_dna: &PackedDna) -> Seq<Vec<Nuc>> {
        Seq(packed_dna.into())
    }
}

impl From<PackedDna> for Vec<Nuc> {
    fn from(packed_dna: PackedDna) -> Vec<Nuc> {
        (&packed_dna).into()
    }
}

impl From<&PackedDna> for Vec<Nuc> {
    fn from(packed_dna: &PackedDna) -> Vec<Nuc> {
        let mut dna = vec![Nuc::default(); packed_dna.len()];
        unpack(&mut dna, &packed_dna.0);
        dna
    }
}

impl<'a> IntoIterator for &'a PackedDna {
    type Item = Nuc;
    type IntoIter = PackedDnaIter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        PackedDnaIter(UnpackingIter::new(0..self.len(), &self.0))
    }
}

impl IntoIterator for PackedDna {
    type Item = Nuc;
    type IntoIter = PackedDnaIntoIter;

    fn into_iter(self) -> Self::IntoIter {
        PackedDnaIntoIter(UnpackingIter::new(0..self.len(), self.0))
    }
}

impl PartialEq for PackedDna {
    fn eq(&self, other: &Self) -> bool {
        match (&*self.0, &*other.0) {
            ([], [0]) | ([0], []) => true,
            (x, y) => x == y,
        }
    }
}

impl Eq for PackedDna {}

impl PartialOrd for PackedDna {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for PackedDna {
    fn cmp(&self, other: &Self) -> Ordering {
        match (&*self.0, &*other.0) {
            ([] | [0], [] | [0]) => Ordering::Equal,
            ([] | [0], _) => Ordering::Less,
            (_, [] | [0]) => Ordering::Greater,
            ([x_bulk @ .., x_last], [y_bulk @ .., y_last]) => {
                let common_prefix_len = x_bulk.len().min(y_bulk.len());
                let (x_prefix, x_remainder) = x_bulk.split_at(common_prefix_len);
                let (y_prefix, y_remainder) = y_bulk.split_at(common_prefix_len);
                x_prefix
                    .cmp(y_prefix)
                    .then_with(|| match (x_remainder, y_remainder) {
                        ([], []) => x_last.cmp(y_last),
                        ([], [y_next, ..]) => (x_last & !0b11).cmp(y_next).then(Ordering::Less),
                        ([x_next, ..], []) => x_next.cmp(&(y_last & !0b11)).then(Ordering::Greater),
                        _ => unreachable!(),
                    })
            }
        }
    }
}

impl std::fmt::Display for PackedDna {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        self.iter().fmt(f)
    }
}

impl std::fmt::Debug for PackedDna {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        f.debug_tuple("PackedDna")
            .field(&display(self.iter()))
            .finish()
    }
}

/// Like [`[Nuc; N]`](array), but takes 75% less space.
///
/// # Examples
///
/// ```
/// use nucs::{Nuc, Packed};
///
/// let dna = Nuc::arr(b"ACATATTAC");
/// assert_eq!(size_of_val(&dna), 9);
/// let packed: Packed<[Nuc; 9]> = dna.into();
/// assert_eq!(size_of_val(&packed), 3);
/// let unpacked: [Nuc; 9] = packed.into();
/// assert_eq!(unpacked, dna);
/// ```
#[derive(Clone, Copy)]
pub struct PackedArrayDna<const N: usize>(PackedBuf<N>)
where
    [(); N]: PackableArray;

type PackedBuf<const N: usize> = <[(); N] as ArrayDivide>::By4<u8>;

impl<const N: usize> PackedArrayDna<N>
where
    [(); N]: PackableArray,
{
    /// Returns an iterator over the packed DNA.
    #[must_use]
    pub fn iter(&self) -> PackedDnaIter<'_> {
        self.into_iter()
    }
}

impl<const N: usize> From<Seq<[Nuc; N]>> for PackedArrayDna<N>
where
    [(); N]: PackableArray,
{
    fn from(dna: Seq<[Nuc; N]>) -> PackedArrayDna<N> {
        dna.0.into()
    }
}

impl<const N: usize> From<&Seq<[Nuc; N]>> for PackedArrayDna<N>
where
    [(); N]: PackableArray,
{
    fn from(dna: &Seq<[Nuc; N]>) -> PackedArrayDna<N> {
        dna.0.into()
    }
}

impl<const N: usize> From<[Nuc; N]> for PackedArrayDna<N>
where
    [(); N]: PackableArray,
{
    fn from(dna: [Nuc; N]) -> PackedArrayDna<N> {
        (&dna).into()
    }
}

impl<const N: usize> From<&[Nuc; N]> for PackedArrayDna<N>
where
    [(); N]: PackableArray,
{
    fn from(dna: &[Nuc; N]) -> PackedArrayDna<N> {
        let mut this = Self(ArrayDefault::array_default());
        pack(this.0.as_mut(), dna);
        this
    }
}

impl<const N: usize> From<PackedArrayDna<N>> for Seq<[Nuc; N]>
where
    [(); N]: PackableArray,
{
    fn from(packed_dna: PackedArrayDna<N>) -> Seq<[Nuc; N]> {
        Seq(packed_dna.into())
    }
}

impl<const N: usize> From<&PackedArrayDna<N>> for Seq<[Nuc; N]>
where
    [(); N]: PackableArray,
{
    fn from(packed_dna: &PackedArrayDna<N>) -> Seq<[Nuc; N]> {
        Seq(packed_dna.into())
    }
}

impl<const N: usize> From<PackedArrayDna<N>> for [Nuc; N]
where
    [(); N]: PackableArray,
{
    fn from(packed_dna: PackedArrayDna<N>) -> [Nuc; N] {
        (&packed_dna).into()
    }
}

impl<const N: usize> From<&PackedArrayDna<N>> for [Nuc; N]
where
    [(); N]: PackableArray,
{
    fn from(packed_dna: &PackedArrayDna<N>) -> [Nuc; N] {
        let mut dna = [Nuc::default(); N];
        unpack(&mut dna, packed_dna.0.as_ref());
        dna
    }
}

impl<'a, const N: usize> IntoIterator for &'a PackedArrayDna<N>
where
    [(); N]: PackableArray,
{
    type Item = Nuc;
    type IntoIter = PackedDnaIter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        PackedDnaIter(UnpackingIter::new(0..N, self.0.as_ref()))
    }
}

impl<const N: usize> IntoIterator for PackedArrayDna<N>
where
    [(); N]: PackableArray,
{
    type Item = Nuc;
    type IntoIter = PackedArrayDnaIntoIter<N>;

    fn into_iter(self) -> Self::IntoIter {
        PackedArrayDnaIntoIter(UnpackingIter::new(0..N, self.0))
    }
}

impl<const N: usize> PartialEq for PackedArrayDna<N>
where
    [(); N]: PackableArray,
{
    fn eq(&self, other: &Self) -> bool {
        self.0.as_ref().eq(other.0.as_ref())
    }
}

impl<const N: usize> Eq for PackedArrayDna<N> where [(); N]: PackableArray {}

impl<const N: usize> PartialOrd for PackedArrayDna<N>
where
    [(); N]: PackableArray,
{
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<const N: usize> Ord for PackedArrayDna<N>
where
    [(); N]: PackableArray,
{
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.as_ref().cmp(other.0.as_ref())
    }
}

impl<const N: usize> std::fmt::Display for PackedArrayDna<N>
where
    [(); N]: PackableArray,
{
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        self.iter().fmt(f)
    }
}

impl<const N: usize> std::fmt::Debug for PackedArrayDna<N>
where
    [(); N]: PackableArray,
{
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        f.debug_tuple("PackedArrayDna")
            .field(&display(self.iter()))
            .finish()
    }
}

fn pack(packed: &mut [u8], dna: &[Nuc]) {
    let (quads, remainder) = dna.as_chunks();
    for (byte, &[n1, n2, n3, n4]) in packed.iter_mut().zip(quads) {
        *byte = n4.compress() | (n3.compress() << 2) | (n2.compress() << 4) | (n1.compress() << 6);
    }
    if let Some(byte) = packed.last_mut()
        && !remainder.is_empty()
    {
        *byte = 0;
        for (offset, &nuc) in [6, 4, 2].into_iter().zip(remainder) {
            *byte |= nuc.compress() << offset;
        }
    }
}

fn unpack(dna: &mut [Nuc], packed: &[u8]) {
    let (quads, remainder) = dna.as_chunks_mut();
    for ([n1, n2, n3, n4], &byte) in quads.iter_mut().zip(packed) {
        *n4 = Nuc::decompress(byte);
        *n3 = Nuc::decompress(byte >> 2);
        *n2 = Nuc::decompress(byte >> 4);
        *n1 = Nuc::decompress(byte >> 6);
    }
    if let Some(byte) = packed.last()
        && !remainder.is_empty()
    {
        for (offset, nuc) in [6, 4, 2].into_iter().zip(remainder) {
            *nuc = Nuc::decompress(byte >> offset);
        }
    }
}

/// Owned [`PackedDna`] iterator.
#[derive(Clone)]
pub struct PackedDnaIntoIter(UnpackingIter<2, 8, std::vec::IntoIter<u8>>);

impl PackedDnaIntoIter {
    fn as_ref(&self) -> PackedDnaIter<'_> {
        PackedDnaIter(self.0.as_ref())
    }
}

impl Iterator for PackedDnaIntoIter {
    type Item = Nuc;

    fn next(&mut self) -> Option<Nuc> {
        self.0
            .next()
            .map(|(shift, byte)| Nuc::decompress(byte >> shift))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.0.size_hint()
    }
}

impl DoubleEndedIterator for PackedDnaIntoIter {
    fn next_back(&mut self) -> Option<Nuc> {
        self.0
            .next_back()
            .map(|(shift, byte)| Nuc::decompress(byte >> shift))
    }
}

impl ExactSizeIterator for PackedDnaIntoIter {
    fn len(&self) -> usize {
        self.0.len()
    }
}

impl std::fmt::Display for PackedDnaIntoIter {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        self.as_ref().fmt(f)
    }
}

impl std::fmt::Debug for PackedDnaIntoIter {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        f.debug_tuple("PackedDnaIntoIter")
            .field(&display(self.as_ref()))
            .finish()
    }
}

/// Owned [`PackedArrayDna`] iterator.
#[derive(Clone)]
pub struct PackedArrayDnaIntoIter<const N: usize>(
    UnpackingIter<2, 8, <PackedBuf<N> as IntoIterator>::IntoIter>,
)
where
    [(); N]: PackableArray;

impl<const N: usize> PackedArrayDnaIntoIter<N>
where
    [(); N]: PackableArray,
{
    fn as_ref(&self) -> PackedDnaIter<'_> {
        PackedDnaIter(self.0.as_ref())
    }
}

impl<const N: usize> Iterator for PackedArrayDnaIntoIter<N>
where
    [(); N]: PackableArray,
{
    type Item = Nuc;

    fn next(&mut self) -> Option<Nuc> {
        self.0
            .next()
            .map(|(shift, byte)| Nuc::decompress(byte >> shift))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.0.size_hint()
    }
}

impl<const N: usize> DoubleEndedIterator for PackedArrayDnaIntoIter<N>
where
    [(); N]: PackableArray,
{
    fn next_back(&mut self) -> Option<Nuc> {
        self.0
            .next_back()
            .map(|(shift, byte)| Nuc::decompress(byte >> shift))
    }
}

impl<const N: usize> ExactSizeIterator for PackedArrayDnaIntoIter<N>
where
    [(); N]: PackableArray,
{
    fn len(&self) -> usize {
        self.0.len()
    }
}

impl<const N: usize> std::fmt::Display for PackedArrayDnaIntoIter<N>
where
    [(); N]: PackableArray,
{
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        self.as_ref().fmt(f)
    }
}

impl<const N: usize> std::fmt::Debug for PackedArrayDnaIntoIter<N>
where
    [(); N]: PackableArray,
{
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        f.debug_tuple("PackedArrayDnaIntoIter")
            .field(&display(self.as_ref()))
            .finish()
    }
}

/// Borrowed packed DNA iterator.
#[derive(Clone)]
pub struct PackedDnaIter<'a>(UnpackingIter<2, 8, std::slice::Iter<'a, u8>>);

impl Iterator for PackedDnaIter<'_> {
    type Item = Nuc;

    fn next(&mut self) -> Option<Nuc> {
        self.0
            .next()
            .map(|(shift, byte)| Nuc::decompress(byte >> shift))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.0.size_hint()
    }
}

impl DoubleEndedIterator for PackedDnaIter<'_> {
    fn next_back(&mut self) -> Option<Nuc> {
        self.0
            .next_back()
            .map(|(shift, byte)| Nuc::decompress(byte >> shift))
    }
}

impl ExactSizeIterator for PackedDnaIter<'_> {
    fn len(&self) -> usize {
        self.0.len()
    }
}

impl std::fmt::Display for PackedDnaIter<'_> {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        display(self.clone()).fmt(f)
    }
}

impl std::fmt::Debug for PackedDnaIter<'_> {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        f.debug_tuple("PackedDnaIter")
            .field(&display(self.clone()))
            .finish()
    }
}

#[cfg(test)]
mod tests {
    use proptest::{arbitrary::any, proptest};

    use super::super::tests::{assert_both_roundtrips, assert_roundtrip};
    use crate::proptest::any_dna;

    use super::*;

    #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
    #[test]
    fn all_short_roundtrips() {
        // Yes, this is hideously ugly, but arrays need a size known at compile-time,
        // and this is the simplest way to ensure I get them all.
        assert_both_roundtrips(&[] as &[Nuc; 0]);
        for n1 in Nuc::ALL {
            assert_both_roundtrips(&[n1]);
            for n2 in Nuc::ALL {
                assert_both_roundtrips(&[n1, n2]);
                for n3 in Nuc::ALL {
                    assert_both_roundtrips(&[n1, n2, n3]);
                    for n4 in Nuc::ALL {
                        assert_both_roundtrips(&[n1, n2, n3, n4]);
                        for n5 in Nuc::ALL {
                            assert_both_roundtrips(&[n1, n2, n3, n4, n5]);
                            for n6 in Nuc::ALL {
                                assert_both_roundtrips(&[n1, n2, n3, n4, n5, n6]);
                                for n7 in Nuc::ALL {
                                    assert_both_roundtrips(&[n1, n2, n3, n4, n5, n6, n7]);
                                    for n8 in Nuc::ALL {
                                        assert_both_roundtrips(&[n1, n2, n3, n4, n5, n6, n7, n8]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // The bulk of the fiddly logic is handled by UnpackingIter, which has more tests.
    // This is just sanity-checking that things were hooked up to it correctly.
    #[test]
    fn smoke_test_iters() {
        let dna = Nuc::seq(b"ACGT")[..].pack();
        assert_eq!(Seq(Vec::from_iter(&dna)), "ACGT");
        assert_eq!(Seq(Vec::from_iter(dna)), "ACGT");
        let dna = Nuc::seq(b"ACGT").pack();
        assert_eq!(Seq(Vec::from_iter(&dna)), "ACGT");
        assert_eq!(Seq(Vec::from_iter(dna)), "ACGT");
    }

    #[test]
    fn display() {
        let dna = Nuc::seq(b"ACGT")[..].pack();
        assert_eq!(dna.to_string(), "ACGT");
        assert_eq!(dna.iter().to_string(), "ACGT");
        assert_eq!(dna.into_iter().to_string(), "ACGT");
        let dna = Nuc::seq(b"ACGT").pack();
        assert_eq!(dna.to_string(), "ACGT");
        assert_eq!(dna.iter().to_string(), "ACGT");
        assert_eq!(dna.into_iter().to_string(), "ACGT");
    }

    proptest! {
        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn packed_dna_roundtrip(
            dna in any_dna(9..50) // 0..=8 covered above
        ) {
            assert_roundtrip(&*dna);
        }

        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn packed_array_dna9_roundtrip(
            dna in any::<[Nuc; 9]>()
        ) {
            assert_roundtrip(&dna);
        }

        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn packed_array_dna10_roundtrip(
            dna in any::<[Nuc; 10]>()
        ) {
            assert_roundtrip(&dna);
        }

        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn packed_array_dna11_roundtrip(
            dna in any::<[Nuc; 11]>()
        ) {
            assert_roundtrip(&dna);
        }

        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn packed_array_dna12_roundtrip(
            dna in any::<[Nuc; 12]>()
        ) {
            assert_roundtrip(&dna);
        }

        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn packed_dna_length(
            dna in any_dna(0..50)
        ) {
            let packed = PackedDna::from(&dna);
            assert_eq!(packed.len(), dna.len());
            assert_eq!(packed.is_empty(), dna.is_empty());
        }

        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn packed_dna_ord(
            dna1 in any_dna(0..50),
            dna2 in any_dna(0..50),
        ) {
            let packed1 = PackedDna::from(&dna1);
            let packed2 = PackedDna::from(&dna2);
            assert_eq!(packed1.cmp(&packed2), dna1.cmp(&dna2));
            assert_eq!(packed1.eq(&packed2), dna1.eq(&dna2));
        }

        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn packed_array_dna1_ord(
            dna1 in any::<[Nuc; 1]>(),
            dna2 in any::<[Nuc; 1]>(),
        ) {
            let packed1 = PackedArrayDna::from(&dna1);
            let packed2 = PackedArrayDna::from(&dna2);
            assert_eq!(packed1.cmp(&packed2), dna1.cmp(&dna2));
            assert_eq!(packed1.eq(&packed2), dna1.eq(&dna2));
        }

        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn packed_array_dna2_ord(
            dna1 in any::<[Nuc; 2]>(),
            dna2 in any::<[Nuc; 2]>(),
        ) {
            let packed1 = PackedArrayDna::from(&dna1);
            let packed2 = PackedArrayDna::from(&dna2);
            assert_eq!(packed1.cmp(&packed2), dna1.cmp(&dna2));
            assert_eq!(packed1.eq(&packed2), dna1.eq(&dna2));
        }

        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn packed_array_dna3_ord(
            dna1 in any::<[Nuc; 3]>(),
            dna2 in any::<[Nuc; 3]>(),
        ) {
            let packed1 = PackedArrayDna::from(&dna1);
            let packed2 = PackedArrayDna::from(&dna2);
            assert_eq!(packed1.cmp(&packed2), dna1.cmp(&dna2));
            assert_eq!(packed1.eq(&packed2), dna1.eq(&dna2));
        }

        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn packed_array_dna4_ord(
            dna1 in any::<[Nuc; 4]>(),
            dna2 in any::<[Nuc; 4]>(),
        ) {
            let packed1 = PackedArrayDna::from(&dna1);
            let packed2 = PackedArrayDna::from(&dna2);
            assert_eq!(packed1.cmp(&packed2), dna1.cmp(&dna2));
            assert_eq!(packed1.eq(&packed2), dna1.eq(&dna2));
        }

        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn packed_array_dna5_ord(
            dna1 in any::<[Nuc; 5]>(),
            dna2 in any::<[Nuc; 5]>(),
        ) {
            let packed1 = PackedArrayDna::from(&dna1);
            let packed2 = PackedArrayDna::from(&dna2);
            assert_eq!(packed1.cmp(&packed2), dna1.cmp(&dna2));
            assert_eq!(packed1.eq(&packed2), dna1.eq(&dna2));
        }
    }
}
