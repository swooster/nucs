use crate::Nuc;

use super::PackableArray;
use super::packable_array::{ArrayDefault, Sealed as ArrayDivide};

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

impl From<&[Nuc]> for PackedDna {
    fn from(dna: &[Nuc]) -> PackedDna {
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

impl From<PackedDna> for Vec<Nuc> {
    fn from(packed_dna: PackedDna) -> Vec<Nuc> {
        (&packed_dna).into()
    }
}

impl From<&PackedDna> for Vec<Nuc> {
    fn from(packed_dna: &PackedDna) -> Vec<Nuc> {
        if let Some(byte) = packed_dna.0.last() {
            let len = 4 * (packed_dna.0.len() - 1) + (byte & 0b11) as usize;
            let mut dna = vec![Nuc::default(); len];
            unpack(&mut dna, &packed_dna.0);
            dna
        } else {
            vec![]
        }
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
pub struct PackedArrayDna<const N: usize>(<[(); N] as ArrayDivide>::By4<u8>)
where
    [(); N]: PackableArray;

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
    }
}
