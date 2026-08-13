//! Types related to packed ambiguous peptides.

use std::fmt::Formatter;

use crate::{AmbiAmino, Seq, iter::display, packed::packable_array::ArrayDefault};

// Note on storage: There's not much room for packing `AmbiAmino`s... we can shave off one byte,
// but that's about it. Again, we store thing big-endian for lexical sorting reasons.
// The good news is that this means a lot of the logic is pretty simple.

/// Like [`Vec<AmbiAmino>`], but takes 25% less space for long peptides.
///
/// # Examples
///
/// ```
/// use nucs::{AmbiAmino, Packed};
///
/// let peptide = AmbiAmino::arr(b"KITTYJUMPSOFDANGER");
/// let packed: Packed<[AmbiAmino]> = peptide.as_slice().into();
/// let unpacked: Vec<AmbiAmino> = packed.into();
/// assert_eq!(unpacked, peptide);
/// ```
#[derive(Clone, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct PackedAmbiPeptide(Vec<[u8; 3]>);

impl PackedAmbiPeptide {
    /// Returns the number of [`AmbiAmino`]s in the packed peptide.
    #[must_use]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Returns `true` if the packed peptide contains no [`AmbiAmino`]s.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns an iterator over the packed peptide.
    #[must_use]
    pub fn iter(&self) -> PackedAmbiPeptideIter<'_> {
        self.into_iter()
    }
}

impl From<Seq<Vec<AmbiAmino>>> for PackedAmbiPeptide {
    fn from(ambi_peptide: Seq<Vec<AmbiAmino>>) -> PackedAmbiPeptide {
        ambi_peptide.0.into()
    }
}

impl From<Vec<AmbiAmino>> for PackedAmbiPeptide {
    fn from(ambi_peptide: Vec<AmbiAmino>) -> PackedAmbiPeptide {
        (&ambi_peptide).into()
    }
}

impl<T: AsRef<[AmbiAmino]> + ?Sized> From<&T> for PackedAmbiPeptide {
    fn from(ambi_peptide: &T) -> PackedAmbiPeptide {
        let ambi_peptide = ambi_peptide.as_ref();
        Self(ambi_peptide.iter().map(|a| a.compress()).collect())
    }
}

impl From<PackedAmbiPeptide> for Seq<Vec<AmbiAmino>> {
    fn from(packed_ambi_peptide: PackedAmbiPeptide) -> Seq<Vec<AmbiAmino>> {
        Seq(packed_ambi_peptide.into())
    }
}

impl From<&PackedAmbiPeptide> for Seq<Vec<AmbiAmino>> {
    fn from(packed_ambi_peptide: &PackedAmbiPeptide) -> Seq<Vec<AmbiAmino>> {
        Seq(packed_ambi_peptide.into())
    }
}

impl From<PackedAmbiPeptide> for Vec<AmbiAmino> {
    fn from(packed_ambi_peptide: PackedAmbiPeptide) -> Vec<AmbiAmino> {
        (&packed_ambi_peptide).into()
    }
}

impl From<&PackedAmbiPeptide> for Vec<AmbiAmino> {
    fn from(packed_ambi_peptide: &PackedAmbiPeptide) -> Vec<AmbiAmino> {
        let packed = packed_ambi_peptide.0.iter();
        packed.map(|&bits| AmbiAmino::decompress(bits)).collect()
    }
}

impl<'a> IntoIterator for &'a PackedAmbiPeptide {
    type Item = AmbiAmino;
    type IntoIter = PackedAmbiPeptideIter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        PackedAmbiPeptideIter(self.0.iter())
    }
}

impl IntoIterator for PackedAmbiPeptide {
    type Item = AmbiAmino;
    type IntoIter = PackedAmbiPeptideIntoIter;

    fn into_iter(self) -> Self::IntoIter {
        PackedAmbiPeptideIntoIter(self.0.into_iter())
    }
}

impl std::fmt::Display for PackedAmbiPeptide {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        self.iter().fmt(f)
    }
}

impl std::fmt::Debug for PackedAmbiPeptide {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        f.debug_tuple("PackedAmbiPeptide")
            .field(&display(self.iter()))
            .finish()
    }
}

/// Like [`[AmbiAmino; N]`](array), but takes 25% less space.
///
/// # Examples
///
/// ```
/// use nucs::{AmbiAmino, Packed};
///
/// let peptide = AmbiAmino::arr(b"KITTYJUMPSOFDANGER");
/// assert_eq!(size_of_val(&peptide), 72);
/// let packed: Packed<[AmbiAmino; 18]> = peptide.into();
/// assert_eq!(size_of_val(&packed), 54);
/// let unpacked: [AmbiAmino; 18] = packed.into();
/// assert_eq!(unpacked, peptide);
/// ```
#[derive(Clone, Copy, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct PackedArrayAmbiPeptide<const N: usize>([[u8; 3]; N]);

impl<const N: usize> PackedArrayAmbiPeptide<N> {
    /// Returns an iterator over the packed peptide.
    #[must_use]
    pub fn iter(&self) -> PackedAmbiPeptideIter<'_> {
        self.into_iter()
    }
}

impl<const N: usize> Default for PackedArrayAmbiPeptide<N> {
    fn default() -> Self {
        Self(ArrayDefault::array_default())
    }
}

impl<const N: usize> From<Seq<[AmbiAmino; N]>> for PackedArrayAmbiPeptide<N> {
    fn from(ambi_peptide: Seq<[AmbiAmino; N]>) -> PackedArrayAmbiPeptide<N> {
        ambi_peptide.0.into()
    }
}

impl<const N: usize> From<&Seq<[AmbiAmino; N]>> for PackedArrayAmbiPeptide<N> {
    fn from(ambi_peptide: &Seq<[AmbiAmino; N]>) -> PackedArrayAmbiPeptide<N> {
        ambi_peptide.0.into()
    }
}

impl<const N: usize> From<[AmbiAmino; N]> for PackedArrayAmbiPeptide<N> {
    fn from(ambi_peptide: [AmbiAmino; N]) -> PackedArrayAmbiPeptide<N> {
        (&ambi_peptide).into()
    }
}

impl<const N: usize> From<&[AmbiAmino; N]> for PackedArrayAmbiPeptide<N> {
    fn from(ambi_peptide: &[AmbiAmino; N]) -> PackedArrayAmbiPeptide<N> {
        Self(ambi_peptide.map(AmbiAmino::compress))
    }
}

impl<const N: usize> From<PackedArrayAmbiPeptide<N>> for Seq<[AmbiAmino; N]> {
    fn from(packed_ambi_peptide: PackedArrayAmbiPeptide<N>) -> Seq<[AmbiAmino; N]> {
        Seq(packed_ambi_peptide.into())
    }
}

impl<const N: usize> From<&PackedArrayAmbiPeptide<N>> for Seq<[AmbiAmino; N]> {
    fn from(packed_ambi_peptide: &PackedArrayAmbiPeptide<N>) -> Seq<[AmbiAmino; N]> {
        Seq(packed_ambi_peptide.into())
    }
}

impl<const N: usize> From<PackedArrayAmbiPeptide<N>> for [AmbiAmino; N] {
    fn from(packed_ambi_peptide: PackedArrayAmbiPeptide<N>) -> [AmbiAmino; N] {
        (&packed_ambi_peptide).into()
    }
}

impl<const N: usize> From<&PackedArrayAmbiPeptide<N>> for [AmbiAmino; N] {
    fn from(packed_ambi_peptide: &PackedArrayAmbiPeptide<N>) -> [AmbiAmino; N] {
        packed_ambi_peptide.0.map(AmbiAmino::decompress)
    }
}

impl<const N: usize> std::fmt::Display for PackedArrayAmbiPeptide<N> {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        self.iter().fmt(f)
    }
}

impl<const N: usize> std::fmt::Debug for PackedArrayAmbiPeptide<N> {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        f.debug_tuple("PackedArrayAmbiPeptide")
            .field(&display(self.iter()))
            .finish()
    }
}

impl<'a, const N: usize> IntoIterator for &'a PackedArrayAmbiPeptide<N> {
    type Item = AmbiAmino;
    type IntoIter = PackedAmbiPeptideIter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        PackedAmbiPeptideIter(self.0.iter())
    }
}

impl<const N: usize> IntoIterator for PackedArrayAmbiPeptide<N> {
    type Item = AmbiAmino;
    type IntoIter = PackedArrayAmbiPeptideIntoIter<N>;

    fn into_iter(self) -> Self::IntoIter {
        PackedArrayAmbiPeptideIntoIter(self.0.into_iter())
    }
}

/// Owned [`PackedAmbiPeptide`] iterator.
#[derive(Clone)]
pub struct PackedAmbiPeptideIntoIter(std::vec::IntoIter<[u8; 3]>);

impl PackedAmbiPeptideIntoIter {
    fn as_ref(&self) -> PackedAmbiPeptideIter<'_> {
        PackedAmbiPeptideIter(self.0.as_slice().iter())
    }
}

impl Iterator for PackedAmbiPeptideIntoIter {
    type Item = AmbiAmino;

    fn next(&mut self) -> Option<AmbiAmino> {
        self.0.next().map(AmbiAmino::decompress)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.0.size_hint()
    }
}

impl DoubleEndedIterator for PackedAmbiPeptideIntoIter {
    fn next_back(&mut self) -> Option<AmbiAmino> {
        self.0.next_back().map(AmbiAmino::decompress)
    }
}

impl ExactSizeIterator for PackedAmbiPeptideIntoIter {
    fn len(&self) -> usize {
        self.0.len()
    }
}

impl std::fmt::Display for PackedAmbiPeptideIntoIter {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        self.as_ref().fmt(f)
    }
}

impl std::fmt::Debug for PackedAmbiPeptideIntoIter {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        f.debug_tuple("PackedAmbiPeptideIntoIter")
            .field(&display(self.as_ref()))
            .finish()
    }
}

/// Owned [`PackedArrayAmbiPeptide`] iterator.
#[derive(Clone)]
pub struct PackedArrayAmbiPeptideIntoIter<const N: usize>(<[[u8; 3]; N] as IntoIterator>::IntoIter);

impl<const N: usize> PackedArrayAmbiPeptideIntoIter<N> {
    fn as_ref(&self) -> PackedAmbiPeptideIter<'_> {
        PackedAmbiPeptideIter(self.0.as_slice().iter())
    }
}

impl<const N: usize> Iterator for PackedArrayAmbiPeptideIntoIter<N> {
    type Item = AmbiAmino;

    fn next(&mut self) -> Option<AmbiAmino> {
        self.0.next().map(AmbiAmino::decompress)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.0.size_hint()
    }
}

impl<const N: usize> DoubleEndedIterator for PackedArrayAmbiPeptideIntoIter<N> {
    fn next_back(&mut self) -> Option<AmbiAmino> {
        self.0.next_back().map(AmbiAmino::decompress)
    }
}

impl<const N: usize> ExactSizeIterator for PackedArrayAmbiPeptideIntoIter<N> {
    fn len(&self) -> usize {
        self.0.len()
    }
}

impl<const N: usize> std::fmt::Display for PackedArrayAmbiPeptideIntoIter<N> {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        self.as_ref().fmt(f)
    }
}

impl<const N: usize> std::fmt::Debug for PackedArrayAmbiPeptideIntoIter<N> {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        f.debug_tuple("PackedArrayAmbiPeptideIntoIter")
            .field(&display(self.as_ref()))
            .finish()
    }
}

/// Borrowed packed ambiguous peptide iterator.
#[derive(Clone)]
pub struct PackedAmbiPeptideIter<'a>(std::slice::Iter<'a, [u8; 3]>);

impl Iterator for PackedAmbiPeptideIter<'_> {
    type Item = AmbiAmino;

    fn next(&mut self) -> Option<AmbiAmino> {
        self.0.next().map(|&bytes| AmbiAmino::decompress(bytes))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.0.size_hint()
    }
}

impl DoubleEndedIterator for PackedAmbiPeptideIter<'_> {
    fn next_back(&mut self) -> Option<AmbiAmino> {
        self.0
            .next_back()
            .map(|&bytes| AmbiAmino::decompress(bytes))
    }
}

impl ExactSizeIterator for PackedAmbiPeptideIter<'_> {
    fn len(&self) -> usize {
        self.0.len()
    }
}

impl std::fmt::Display for PackedAmbiPeptideIter<'_> {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        display(self.clone()).fmt(f)
    }
}

impl std::fmt::Debug for PackedAmbiPeptideIter<'_> {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        f.debug_tuple("PackedAmbiPeptideIter")
            .field(&display(self.clone()))
            .finish()
    }
}

#[cfg(test)]
mod tests {
    use proptest::proptest;

    use super::super::tests::{assert_both_roundtrips, assert_roundtrip};
    use crate::proptest::any_ambi_peptide;

    use super::*;

    #[test]
    fn sanity_check_array_roundtrips() {
        assert_both_roundtrips(&[] as &[AmbiAmino; 0]);
        assert_both_roundtrips(&[AmbiAmino::A]);
        assert_both_roundtrips(&[AmbiAmino::B, AmbiAmino::C]);
        assert_both_roundtrips(&[AmbiAmino::D, AmbiAmino::E, AmbiAmino::F]);
        assert_both_roundtrips(&AmbiAmino::arr(b"*SIGN"));
    }

    // This is just sanity-checking that things were hooked up to it correctly.
    #[test]
    fn smoke_test_iters() {
        let ambi_peptide = AmbiAmino::seq(b"AMBI*PEPTIDE")[..].pack();
        assert_eq!(Seq(Vec::from_iter(&ambi_peptide)), "AMBI*PEPTIDE");
        assert_eq!(Seq(Vec::from_iter(ambi_peptide)), "AMBI*PEPTIDE");
        let ambi_peptide = AmbiAmino::seq(b"AMBI*PEPTIDE").pack();
        assert_eq!(Seq(Vec::from_iter(&ambi_peptide)), "AMBI*PEPTIDE");
        assert_eq!(Seq(Vec::from_iter(ambi_peptide)), "AMBI*PEPTIDE");
    }

    #[test]
    fn display() {
        let ambi_peptide = AmbiAmino::seq(b"AMBI*PEPTIDE")[..].pack();
        assert_eq!(ambi_peptide.to_string(), "AMBI*PEPTIDE");
        assert_eq!(ambi_peptide.iter().to_string(), "AMBI*PEPTIDE");
        assert_eq!(ambi_peptide.into_iter().to_string(), "AMBI*PEPTIDE");
        let ambi_peptide = AmbiAmino::seq(b"AMBI*PEPTIDE").pack();
        assert_eq!(ambi_peptide.to_string(), "AMBI*PEPTIDE");
        assert_eq!(ambi_peptide.iter().to_string(), "AMBI*PEPTIDE");
        assert_eq!(ambi_peptide.into_iter().to_string(), "AMBI*PEPTIDE");
    }

    proptest! {
        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn packed_ambi_peptide_roundtrip(
            ambi_peptide in any_ambi_peptide(1..25) // 0 covered above
        ) {
            assert_roundtrip(&*ambi_peptide);
        }

        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn packed_ambi_peptide_lexical_ordering(
            peptide1 in any_ambi_peptide(0..50),
            peptide2 in any_ambi_peptide(0..50),
        ) {
            let packed1 = PackedAmbiPeptide::from(peptide1.as_slice());
            let packed2 = PackedAmbiPeptide::from(peptide2.as_slice());
            assert_eq!(packed1.0.cmp(&packed2.0), peptide1.cmp(&peptide2));
        }

        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn packed_ambi_peptide_length(
            ambi_peptide in any_ambi_peptide(0..50)
        ) {
            let packed = PackedAmbiPeptide::from(&ambi_peptide);
            assert_eq!(packed.len(), ambi_peptide.len());
            assert_eq!(packed.is_empty(), ambi_peptide.is_empty());
        }
    }
}
