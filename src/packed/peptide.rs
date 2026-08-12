use crate::{Amino, Seq};

use super::PackableArray;
use super::packable_array::{ArrayDefault, Sealed as ArrayDivide};

// Note on storage: Aminos are packed big-endian so naive lexical sorting of bytes is correct.
// The [u8; 2] themselves are big-endian u16s, for the same reason.
//
// The last word of a non-empty `PackedPeptide` must have 1-3 elements. If it has 2 elements,
// the 5 least-significant bits are 0. If it has 1 element, the 10 least-significant bits would
// be 0, so only the upper byte is stored; this means the amino is at bits 2..7, not 0..5.

// NOTE to future self: when implementing iterators, it'll make the most sense to cast to
// &[[u8; 2]] then chain its iter to an Option<Amino>.

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
#[derive(Clone)]
pub struct PackedPeptide(Vec<u8>);

impl PackedPeptide {
    /// Returns the number of [`Amino`]s in the packed peptide.
    #[must_use]
    pub fn len(&self) -> usize {
        // len <= `usize::MAX` is an invariant of this type, so no need to worry about overflow.
        match self.0.as_chunks() {
            ([], []) => 0,
            (bulk @ [.., last], []) => {
                3 * bulk.len() - (u16::from_be_bytes(*last).trailing_zeros() / 5) as usize
            }
            (bulk, [_]) => 3 * bulk.len() + 1,
            (_, [_, _, ..]) => unreachable!(),
        }
    }

    /// Returns `true` if the packed peptide contains no [`Amino`]s.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
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
        let packed_len = (2 * peptide.len()).div_ceil(3);
        let mut packed = vec![0; packed_len];
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
#[derive(Clone, Copy)]
pub struct PackedArrayPeptide<const N: usize>(<[(); N] as ArrayDivide>::By3_2<u8>)
where
    [(); N]: PackableArray;

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

fn pack(packed: &mut [u8], peptide: &[Amino]) {
    let (pairs, packed_remainder) = packed.as_chunks_mut();
    let (triplets, peptide_remainder) = peptide.as_chunks();
    for (pair, &[a1, a2, a3]) in pairs.iter_mut().zip(triplets) {
        *pair = (a3.compress() | (a2.compress() << 5) | (a1.compress() << 10)).to_be_bytes();
    }
    match (packed_remainder, peptide_remainder) {
        ([], []) => {}
        ([b1], &[a1]) => [*b1, _] = (a1.compress() << 10).to_be_bytes(),
        ([], &[a1, a2]) => {
            *pairs.last_mut().unwrap() =
                ((a2.compress() << 5) | (a1.compress() << 10)).to_be_bytes();
        }
        _ => panic!(),
    }
}

fn unpack(peptide: &mut [Amino], packed: &[u8]) {
    let (pairs, packed_remainder) = packed.as_chunks();
    let (triplets, peptide_remainder) = peptide.as_chunks_mut();
    for ([a1, a2, a3], &pair) in triplets.iter_mut().zip(pairs) {
        let val = u16::from_be_bytes(pair);
        *a1 = Amino::decompress(val >> 10);
        *a2 = Amino::decompress(val >> 5);
        *a3 = Amino::decompress(val);
    }
    match (peptide_remainder, packed_remainder) {
        ([], []) => {}
        ([a1], &[b1]) => {
            *a1 = Amino::decompress(u16::from_be_bytes([b1, 0]) >> 10);
        }
        ([a1, a2], []) => {
            let val = u16::from_be_bytes(*pairs.last().unwrap());
            *a1 = Amino::decompress(val >> 10);
            *a2 = Amino::decompress(val >> 5);
        }
        _ => panic!(),
    }
}

#[cfg(test)]
mod tests {
    use proptest::{arbitrary::any, proptest};

    use super::super::tests::{assert_both_roundtrips, assert_roundtrip};
    use crate::proptest::any_peptide;

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
    }
}
