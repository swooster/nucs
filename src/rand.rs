//! `rand` trait implementations

use rand::Rng;
use rand::distr::{Distribution, StandardUniform};

use crate::{AmbiAmino, AmbiNuc, Amino, Nuc};

// NOTE: The disperate styles of random generation are a result of what I've found to be
// fastest for each type. Also note that the size of the `expect` messages seems to influence
// speed even though it appears to be dead code..? ¯\_(ツ)_/¯

impl Distribution<Nuc> for StandardUniform {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Nuc {
        let num_nucs = u8::try_from(Nuc::ALL.len()).expect("Nuc should be repr(u8)");
        let index: u8 = rng.random_range(0..num_nucs);
        match index {
            0 => Nuc::A,
            1 => Nuc::C,
            2 => Nuc::G,
            3 => Nuc::T,
            _ => unreachable!("BUG: bad RNG Nuc"),
        }
    }
}

impl Distribution<AmbiNuc> for StandardUniform {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> AmbiNuc {
        let num_ambi_nucs = u8::try_from(AmbiNuc::ALL.len()).expect("AmbiNuc should be repr(u8)");
        #[allow(clippy::range_plus_one, reason = "this would halve speed")]
        let index: u8 = rng.random_range(1..(num_ambi_nucs + 1));
        AmbiNuc::from_u8(index).expect("BUG: bad RNG AmbiNuc")
    }
}

impl Distribution<Amino> for StandardUniform {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Amino {
        let num_aminos = u8::try_from(Amino::ALL.len()).expect("Amino should be repr(u8)");
        let index: u8 = rng.random_range(0..num_aminos);
        Amino::ALL[index as usize]
    }
}

impl Distribution<AmbiAmino> for StandardUniform {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> AmbiAmino {
        AmbiAmino::from_bits(rng.random_range(AmbiAmino::BITS_RANGE))
    }
}
