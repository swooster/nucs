//! `rand` trait implementations

use rand::Rng;
use rand::distr::{Distribution, StandardUniform};
use rand::seq::IndexedRandom;

use crate::{AmbiAmino, AmbiNuc, Amino, Nuc};

impl Distribution<Nuc> for StandardUniform {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Nuc {
        *(Nuc::ALL).choose(rng).expect("BUG: Nuc::ALL was empty")
    }
}

impl Distribution<AmbiNuc> for StandardUniform {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> AmbiNuc {
        *(AmbiNuc::ALL)
            .choose(rng)
            .expect("BUG: AmbiNuc::ALL was empty")
    }
}

impl Distribution<Amino> for StandardUniform {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Amino {
        *(Amino::ALL).choose(rng).expect("BUG: Amino::ALL was empty")
    }
}

impl Distribution<AmbiAmino> for StandardUniform {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> AmbiAmino {
        AmbiAmino::from_bits(rng.random_range(AmbiAmino::BITS_RANGE))
    }
}
