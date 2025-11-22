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
        AmbiAmino::BITS_RANGE
            .prop_map(AmbiAmino::from_bits)
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
