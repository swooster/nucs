//! Packed sequence types

use crate::{AmbiAmino, AmbiNuc, Amino, Nuc};

mod ambidna;
mod ambipeptide;
mod dna;
mod peptide;

pub use ambidna::{PackedAmbiDna, PackedArrayAmbiDna};
pub use ambipeptide::{PackedAmbiPeptide, PackedArrayAmbiPeptide};
pub use dna::{PackedArrayDna, PackedDna};
pub use peptide::{PackedArrayPeptide, PackedPeptide};

/// Names packed types in terms of what data they're packing.
///
/// * [`Packed<[Nuc]>`](Self) is [`PackedDna`]
/// * [`Packed<[Nuc; N]>`](Self) is [`PackedArrayDna<N>`]
/// * [`Packed<[AmbiNuc]>`](Self) is [`PackedAmbiDna`]
/// * [`Packed<[AmbiNuc; N]>`](Self) is [`PackedArrayAmbiDna<N>`]
/// * [`Packed<[Amino]>`](Self) is [`PackedPeptide`]
/// * [`Packed<[Amino; N]>`](Self) is [`PackedArrayPeptide<N>`]
/// * [`Packed<[AmbiAmino]>`](Self) is [`PackedAmbiPeptide`]
/// * [`Packed<[AmbiAmino; N]>`](Self) is [`PackedArrayAmbiPeptide<N>`]
pub type Packed<T> = <T as Packable>::Packed;

/// Determine packed version of type.
///
/// This should work for all `[impl Symbol]`, and sufficiently short `[impl Symbol; N]`
/// (currently, `N <= 64`).
pub trait Packable: ToOwned {
    /// Packed version of [`Self`].
    type Packed: Into<Self::Owned> + for<'a> From<&'a Self>;
}

impl Packable for [Nuc] {
    type Packed = PackedDna;
}

impl Packable for Vec<Nuc> {
    type Packed = PackedDna;
}

impl<const N: usize> Packable for [Nuc; N]
where
    [(); N]: PackableArray,
{
    type Packed = PackedArrayDna<N>;
}

impl Packable for [AmbiNuc] {
    type Packed = PackedAmbiDna;
}

impl Packable for Vec<AmbiNuc> {
    type Packed = PackedAmbiDna;
}

impl<const N: usize> Packable for [AmbiNuc; N]
where
    [(); N]: PackableArray,
{
    type Packed = PackedArrayAmbiDna<N>;
}

impl Packable for [Amino] {
    type Packed = PackedPeptide;
}

impl Packable for Vec<Amino> {
    type Packed = PackedPeptide;
}

impl<const N: usize> Packable for [Amino; N]
where
    [(); N]: PackableArray,
{
    type Packed = PackedArrayPeptide<N>;
}

impl Packable for [AmbiAmino] {
    type Packed = PackedAmbiPeptide;
}

impl Packable for Vec<AmbiAmino> {
    type Packed = PackedAmbiPeptide;
}

impl<const N: usize> Packable for [AmbiAmino; N] {
    type Packed = PackedArrayAmbiPeptide<N>;
}

/// An array for which packed sizes are known.
pub trait PackableArray: packable_array::Sealed {}

pub(crate) mod packable_array {
    /// Implementation details of [`PackedArray`].
    ///
    /// This uses GATs rather than constants to work around limitations in the current
    /// trait solver: A type/trait can't declare an array with a length obtained from a
    /// trait's associated constant.
    ///
    /// In theory, we could use non-generic associated types by setting them to something
    /// like `[(); N]` and matching a const generic parameter against them, but that ends up
    /// introducing additional generic parameters into types, and more implementation-specific
    /// trait constraints into everything else, complicating the public APIs.
    pub trait Sealed {
        /// `[T; N.div_ceil(2)]`
        type By2<T: Copy + Default>: AsRef<[T]> + AsMut<[T]> + Copy + ArrayDefault;
        /// `[T; N.div_ceil(3)]`
        type By3<T: Copy + Default>: AsRef<[T]> + AsMut<[T]> + Copy + ArrayDefault;
        /// `[T; N.div_ceil(4)]`
        type By4<T: Copy + Default>: AsRef<[T]> + AsMut<[T]> + Copy + ArrayDefault;
    }

    /// Workaround for `[T; N]: Default` being limited to `N <= 32`.
    pub trait ArrayDefault {
        fn array_default() -> Self;
    }

    impl<T: Default, const N: usize> ArrayDefault for [T; N] {
        fn array_default() -> Self {
            std::array::from_fn(|_| T::default())
        }
    }
}

macro_rules! div_table {
    ( $( $d:literal => $q2:literal $q3:literal $q4:literal),+ , ) => {
        $(
            impl PackableArray for [(); $d] {}

            impl packable_array::Sealed for [(); $d] {
                type By2<T: Copy + Default> = [T; $q2];
                type By3<T: Copy + Default> = [T; $q3];
                type By4<T: Copy + Default> = [T; $q4];
            }
        )+
    };
}

div_table!(
    // x => x/2 x/3 x/4
       0 =>   0   0   0,
       1 =>   1   1   1,
       2 =>   1   1   1,
       3 =>   2   1   1,
       4 =>   2   2   1,
       5 =>   3   2   2,
       6 =>   3   2   2,
       7 =>   4   3   2,
       8 =>   4   3   2,
       9 =>   5   3   3,
      10 =>   5   4   3,
      11 =>   6   4   3,
      12 =>   6   4   3,
      13 =>   7   5   4,
      14 =>   7   5   4,
      15 =>   8   5   4,
      16 =>   8   6   4,
      17 =>   9   6   5,
      18 =>   9   6   5,
      19 =>  10   7   5,
      20 =>  10   7   5,
      21 =>  11   7   6,
      22 =>  11   8   6,
      23 =>  12   8   6,
      24 =>  12   8   6,
      25 =>  13   9   7,
      26 =>  13   9   7,
      27 =>  14   9   7,
      28 =>  14  10   7,
      29 =>  15  10   8,
      30 =>  15  10   8,
      31 =>  16  11   8,
      32 =>  16  11   8,
      33 =>  17  11   9,
      34 =>  17  12   9,
      35 =>  18  12   9,
      36 =>  18  12   9,
      37 =>  19  13  10,
      38 =>  19  13  10,
      39 =>  20  13  10,
      40 =>  20  14  10,
      41 =>  21  14  11,
      42 =>  21  14  11,
      43 =>  22  15  11,
      44 =>  22  15  11,
      45 =>  23  15  12,
      46 =>  23  16  12,
      47 =>  24  16  12,
      48 =>  24  16  12,
      49 =>  25  17  13,
      50 =>  25  17  13,
      51 =>  26  17  13,
      52 =>  26  18  13,
      53 =>  27  18  14,
      54 =>  27  18  14,
      55 =>  28  19  14,
      56 =>  28  19  14,
      57 =>  29  19  15,
      58 =>  29  20  15,
      59 =>  30  20  15,
      60 =>  30  20  15,
      61 =>  31  21  16,
      62 =>  31  21  16,
      63 =>  32  21  16,
      64 =>  32  22  16,
    // TODO: Maybe consider up to 256? Probably overkill though.
);

#[cfg(test)]
mod tests {
    use std::fmt::Debug;

    use super::*;

    #[track_caller]
    pub(crate) fn assert_both_roundtrips<T, const N: usize>(arr: &[T; N])
    where
        T: Debug + PartialEq,
        [T]: Packable<Owned = Vec<T>>,
        [T; N]: Packable<Owned = [T; N]>,
    {
        assert_roundtrip(arr.as_slice());
        assert_roundtrip(arr);
    }

    #[track_caller]
    pub(crate) fn assert_roundtrip<T>(packable: &T)
    where
        T: Packable + Debug + ?Sized,
        T::Owned: Debug + PartialEq<T>,
    {
        let roundtripped: T::Owned = <T::Packed>::from(packable).into();
        assert_eq!(
            roundtripped,
            *packable,
            "roundtrip failed for {}",
            std::any::type_name::<T>()
        );
    }
}
