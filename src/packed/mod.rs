//! Packed sequence types

use std::hash::Hash;
use std::ops::Range;

use crate::{AmbiAmino, AmbiNuc, Amino, Nuc, packed::packable_array::ContiguousIterator};

pub mod ambidna;
pub mod ambipeptide;
pub mod dna;
pub mod peptide;

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
    use std::hash::Hash;

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
        type By2<T: Copy + Default + Eq + Ord + Hash>: AsRef<[T]>
            + AsMut<[T]>
            + Copy
            + Eq
            + Ord
            + Hash
            + ArrayDefault
            + IntoIterator<Item = T, IntoIter: ContiguousIterator<SliceItem = T>>;
        /// `[T; N.div_ceil(3)]`
        type By3<T: Copy + Default + Eq + Ord + Hash>: AsRef<[T]>
            + AsMut<[T]>
            + Copy
            + Eq
            + Ord
            + Hash
            + ArrayDefault
            + IntoIterator<Item = T, IntoIter: ContiguousIterator<SliceItem = T>>;
        /// `[T; N.div_ceil(4)]`
        type By4<T: Copy + Default + Eq + Ord + Hash>: AsRef<[T]>
            + AsMut<[T]>
            + Copy
            + Eq
            + Ord
            + Hash
            + ArrayDefault
            + IntoIterator<Item = T, IntoIter: ContiguousIterator<SliceItem = T>>;
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

    pub trait ContiguousIterator:
        DoubleEndedIterator + ExactSizeIterator + Clone + Default
    {
        type SliceItem;
        fn as_slice(&self) -> &[Self::SliceItem];
        fn peek_front(&self) -> Option<Self::Item>;
        fn peek_back(&self) -> Option<Self::Item>;
    }

    impl<T: Copy, const N: usize> ContiguousIterator for std::array::IntoIter<T, N> {
        type SliceItem = T;

        fn as_slice(&self) -> &[Self::SliceItem] {
            self.as_slice()
        }

        fn peek_front(&self) -> Option<Self::Item> {
            self.as_slice().first().copied()
        }

        fn peek_back(&self) -> Option<Self::Item> {
            self.as_slice().last().copied()
        }
    }

    impl<T: Copy> ContiguousIterator for std::vec::IntoIter<T> {
        type SliceItem = T;

        fn as_slice(&self) -> &[Self::SliceItem] {
            self.as_slice()
        }

        fn peek_front(&self) -> Option<Self::Item> {
            self.as_slice().first().copied()
        }

        fn peek_back(&self) -> Option<Self::Item> {
            self.as_slice().last().copied()
        }
    }

    impl<T> ContiguousIterator for std::slice::Iter<'_, T> {
        type SliceItem = T;

        fn as_slice(&self) -> &[Self::SliceItem] {
            self.as_slice()
        }

        fn peek_front(&self) -> Option<Self::Item> {
            self.as_slice().first()
        }

        fn peek_back(&self) -> Option<Self::Item> {
            self.as_slice().last()
        }
    }
}

macro_rules! div_table {
    ( $( $d:literal => $q2:literal $q3:literal $q4:literal),+ , ) => {
        $(
            impl PackableArray for [(); $d] {}

            impl packable_array::Sealed for [(); $d] {
                type By2<T: Copy + Default + Eq + Ord + Hash> = [T; $q2];
                type By3<T: Copy + Default + Eq + Ord + Hash> = [T; $q3];
                type By4<T: Copy + Default + Eq + Ord + Hash> = [T; $q4];
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

/// Foundation of all packed iteration.
///
/// `UnpackingIter<STEP, MAX, I>` behaves like
/// `itertools::iproduct(iter, shifts).map(|(x, s)| (s, x))` where
/// `shifts` is `(0..MAX).step_by(STEP).rev()` and `iter` is an `I`.
///
/// Note that `STEP` must evenly divide `MAX`.
#[derive(Clone, Default, Debug)]
pub(crate) struct UnpackingIter<const STEP: u8, const MAX: u8, I> {
    iter: I,
    front_shift: u8,
    back_shift: u8,
}

impl<const STEP: u8, const MAX: u8, I: ContiguousIterator> UnpackingIter<STEP, MAX, I> {
    // NOTE: caller must check that range is valid
    pub(crate) fn new(range: Range<usize>, iterable: impl IntoIterator<IntoIter = I>) -> Self {
        let Range { start, end } = range;
        if start < end {
            let (front_idx, front_shift) = Self::idx_and_shift(start);
            let (back_idx, back_shift) = Self::idx_and_shift(end - 1);
            let mut iter = iterable.into_iter();
            if iter.len() >= 2 && iter.len() - 2 >= back_idx {
                iter.nth_back(iter.len() - 2 - back_idx);
            }
            if front_idx >= 1 {
                iter.nth(front_idx - 1);
            }
            Self {
                iter,
                front_shift,
                back_shift,
            }
        } else {
            Self::default()
        }
    }

    fn idx_and_shift(packed_idx: usize) -> (usize, u8) {
        let idx = packed_idx / (MAX / STEP) as usize;
        let word_idx = packed_idx % (MAX / STEP) as usize;
        let word_idx = u8::try_from(word_idx).expect("division by u8 should produce u8 remainder");
        let shift = MAX - STEP - STEP * word_idx;
        (idx, shift)
    }

    fn as_ref(&self) -> UnpackingIter<STEP, MAX, std::slice::Iter<'_, I::SliceItem>> {
        UnpackingIter {
            iter: self.iter.as_slice().iter(),
            front_shift: self.front_shift,
            back_shift: self.back_shift,
        }
    }
}

impl<const STEP: u8, const MAX: u8, I> Iterator for UnpackingIter<STEP, MAX, I>
where
    I: packable_array::ContiguousIterator<Item: Copy>,
{
    type Item = (u8, I::Item);

    fn next(&mut self) -> Option<Self::Item> {
        let item = (self.front_shift, self.iter.peek_front()?);
        self.front_shift = (self.front_shift + MAX - STEP) % MAX;
        if self.front_shift == MAX - STEP
            || (self.iter.len() == 1 && self.front_shift < self.back_shift)
        {
            self.iter.next();
        }
        Some(item)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = self.len();
        (len, Some(len))
    }
}

impl<const STEP: u8, const MAX: u8, I> DoubleEndedIterator for UnpackingIter<STEP, MAX, I>
where
    I: packable_array::ContiguousIterator<Item: Copy>,
{
    fn next_back(&mut self) -> Option<Self::Item> {
        let item = (self.back_shift, self.iter.peek_back()?);
        self.back_shift = (self.back_shift + STEP) % MAX;
        if self.back_shift == 0 || (self.iter.len() == 1 && self.front_shift < self.back_shift) {
            self.iter.next_back();
        }
        Some(item)
    }
}

impl<const STEP: u8, const MAX: u8, I> ExactSizeIterator for UnpackingIter<STEP, MAX, I>
where
    I: packable_array::ContiguousIterator<Item: Copy>,
{
    fn len(&self) -> usize {
        let skipped = (MAX - STEP - self.front_shift + self.back_shift) / STEP;
        // We use saturating sub because we don't reset the shifts on exhaustion.
        ((MAX / STEP) as usize * self.iter.len()).saturating_sub(skipped as usize)
    }
}

#[cfg(test)]
mod tests {
    use std::fmt::Debug;

    use proptest::arbitrary::any;
    use proptest::proptest;

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

    #[test]
    fn smoke_test_unpacking_iter() {
        let expected_values = [
            (12, "foobar"),
            (9, "foobar"),
            (6, "foobar"),
            (3, "foobar"),
            (0, "foobar"),
            (12, "XYZZY"),
            (9, "XYZZY"),
            (6, "XYZZY"),
            (3, "XYZZY"),
            (0, "XYZZY"),
            (12, "baz"),
            (9, "baz"),
            (6, "baz"),
            (3, "baz"),
            (0, "baz"),
            (12, "QUX"),
            (9, "QUX"),
            (6, "QUX"),
            (3, "QUX"),
            (0, "QUX"),
        ];
        for start in 0..=expected_values.len() {
            for end in start..=expected_values.len() {
                let data = ["foobar", "XYZZY", "baz", "QUX"];
                let iter = UnpackingIter::<3, 15, _>::new(start..end, data);
                assert_eq!(iter.collect::<Vec<_>>(), expected_values[start..end]);
            }
        }
    }

    #[test]
    fn empty_unpacking_iter() {
        let iter = || UnpackingIter::<2, 4, _>::new(0..0, [0u32; 0]);
        assert_eq!(iter().len(), 0);
        assert_eq!(iter().next(), None);
        assert_eq!(iter().next_back(), None);
    }

    proptest! {
        #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
        #[test]
        fn unpacking_iter_matches_reference(
            data_and_ends in proptest::collection::vec((any::<u32>(), any::<bool>()), 0..20),
            last_end in any::<bool>(),
            i in any::<usize>(),
            j in any::<usize>(),
        ) {
            let (data, ends): (Vec<_>, Vec<_>) = data_and_ends.into_iter().unzip();
            let i = i % (data.len() + 1);
            let j = j % (data.len() + 1);
            let range = i.min(j)..i.max(j);

            let shifts = [6, 4, 2, 0u8];
            let reference: Vec<_> = data
                .iter()
                .copied()
                .flat_map(|v| shifts.map(|s| (s, v)))
                .collect();

            let mut iter = UnpackingIter::<2, 8, _>::new(range.clone(), data);
            let mut reference_iter = reference[range].iter().copied();

            assert_eq!(iter.len(), reference_iter.len(), "initial lengths");
            for end in ends {
                if end {
                    assert_eq!(iter.next(), reference_iter.next(), "items");
                } else {
                    assert_eq!(iter.next_back(), reference_iter.next_back(), "items");
                }
                assert_eq!(iter.len(), reference_iter.len(), "lengths");
            }
            assert_eq!(iter.len(), 0, "empty iter length");
            assert_eq!(reference_iter.len(), 0, "empty refernece iter length");
            if last_end {
                assert_eq!(iter.next(), None, "iter is exhausted");
                assert_eq!(reference_iter.next(), None, "reference is exhausted");
            } else {
                assert_eq!(iter.next_back(), None, "iter is exhausted");
                assert_eq!(reference_iter.next_back(), None, "reference is exhausted");
            }
        }
    }
}
