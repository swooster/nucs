use crate::{AmbiAmino, Seq};

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
#[derive(Clone)]
pub struct PackedAmbiPeptide(Vec<[u8; 3]>);

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
#[derive(Clone, Copy)]
pub struct PackedArrayAmbiPeptide<const N: usize>([[u8; 3]; N]);

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
    }
}
