use crate::{AmbiNuc, Seq};

use super::PackableArray;
use super::packable_array::{ArrayDefault, Sealed as ArrayDivide};

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
        if let Some(byte) = packed_ambi_dna.0.last() {
            let len = 2 * packed_ambi_dna.0.len() - usize::from(byte.trailing_zeros() >= 4);
            let mut ambi_dna = vec![AmbiNuc::default(); len];
            unpack(&mut ambi_dna, &packed_ambi_dna.0);
            ambi_dna
        } else {
            vec![]
        }
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
pub struct PackedArrayAmbiDna<const N: usize>(<[(); N] as ArrayDivide>::By2<u8>)
where
    [(); N]: PackableArray;

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
    }
}
