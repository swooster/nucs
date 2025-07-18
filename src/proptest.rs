//! Helpers for working with [`proptest`]

use proptest::arbitrary::{Arbitrary, arbitrary};
use proptest::collection::SizeRange;
use proptest::sample::Select;
use proptest::strategy::{SBoxedStrategy, Strategy};

use crate::{AmbiAmino, AmbiNuc, Amino, Nuc};

impl Arbitrary for Nuc {
    type Parameters = ();
    type Strategy = Select<Nuc>;

    fn arbitrary_with(_args: ()) -> Self::Strategy {
        proptest::sample::select(&Nuc::ALL)
    }
}

impl Arbitrary for AmbiNuc {
    type Parameters = ();
    type Strategy = Select<AmbiNuc>;

    fn arbitrary_with(_args: ()) -> Self::Strategy {
        proptest::sample::select(&AmbiNuc::ALL)
    }
}

impl Arbitrary for Amino {
    type Parameters = ();
    type Strategy = Select<Amino>;

    fn arbitrary_with(_args: ()) -> Self::Strategy {
        proptest::sample::select(&Amino::ALL)
    }
}

impl Arbitrary for AmbiAmino {
    type Parameters = ();
    type Strategy = SBoxedStrategy<AmbiAmino>;

    fn arbitrary_with(_args: ()) -> Self::Strategy {
        (1u32..1 << Amino::ALL.len())
            .prop_map(|bits| {
                Amino::ALL
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| bits & (1 << i) != 0)
                    .map(|(_, aa)| AmbiAmino::from(*aa))
                    .reduce(|a, b| a | b)
                    .expect("BUG: bits must have been zero despite requesting at least 1")
            })
            .sboxed()
    }
}

/// Helper for generating arbitrary [`Vec<Nuc>`]
pub fn any_dna<S: Into<SizeRange>>(size: S) -> VecStrategy<Nuc> {
    proptest::collection::vec(arbitrary(), size)
}

/// Helper for generating arbitrary [`Vec<AmbiNuc>`]
pub fn any_ambi_dna<S: Into<SizeRange>>(size: S) -> VecStrategy<AmbiNuc> {
    proptest::collection::vec(arbitrary(), size)
}

/// Helper for generating arbitrary [`Vec<Amino>`]
pub fn any_peptide<S: Into<SizeRange>>(size: S) -> VecStrategy<Amino> {
    proptest::collection::vec(arbitrary(), size)
}

/// Helper for generating arbitrary [`Vec<AmbiAmino>`]
pub fn any_ambi_peptide<S: Into<SizeRange>>(size: S) -> VecStrategy<AmbiAmino> {
    proptest::collection::vec(arbitrary(), size)
}

/// [`Strategy`] for generating [`Vec<S>`]
pub type VecStrategy<S> = proptest::collection::VecStrategy<<S as Arbitrary>::Strategy>;
