//! Types related to translation of codons into amino acids.

use crate::{AmbiAmino, AmbiNuc, Amino, Nuc};

/// Trait representing any type that can be used to translate codons into amino acids.
///
/// The most commonly used [`GeneticCode`] is [`NCBI1`] but there are others in
/// [`nucs::translation`](crate::translation).
pub trait GeneticCode {
    /// Map a concrete codon to an amino acid
    fn map_codon(&self, codon: [Nuc; 3]) -> Amino;

    /// Map an ambiguous codon to an amino acid
    fn map_ambi_codon(&self, codon: [AmbiNuc; 3]) -> AmbiAmino {
        let [ambi_n1, ambi_n2, ambi_n3] = codon;
        ambi_n1
            .iter()
            .flat_map(|n1| ambi_n2.iter().map(move |n2| [n1, n2]))
            .flat_map(|[n1, n2]| ambi_n3.iter().map(move |n3| [n1, n2, n3]))
            .map(|codon| AmbiAmino::from(self.map_codon(codon)))
            .reduce(|a, b| a | b)
            .expect("BUG: null nucleotide encountered")
    }
}

impl<F: Fn([Nuc; 3]) -> Amino> GeneticCode for F {
    fn map_codon(&self, codon: [Nuc; 3]) -> Amino {
        self(codon)
    }
}

/// A table of amino acids can be used for genetic coding.
///
/// The [`Amino`]s must be ordered to correspond with:
/// `AAA`, `AAC`, `AAG`, `AAT`, `ACA`, `ACC`, `ACG`, `ACT`,
/// `AGA`, `AGC`, `AGG`, `AGT`, `ATA`, `ATC`, `ATG`, `ATT`,
/// `CAA`, `CAC`, `CAG`, `CAT`, `CCA`, `CCC`, `CCG`, `CCT`,
/// `CGA`, `CGC`, `CGG`, `CGT`, `CTA`, `CTC`, `CTG`, `CTT`,
/// `GAA`, `GAC`, `GAG`, `GAT`, `GCA`, `GCC`, `GCG`, `GCT`,
/// `GGA`, `GGC`, `GGG`, `GGT`, `GTA`, `GTC`, `GTG`, `GTT`,
/// `TAA`, `TAC`, `TAG`, `TAT`, `TCA`, `TCC`, `TCG`, `TCT`,
/// `TGA`, `TGC`, `TGG`, `TGT`, `TTA`, `TTC`, `TTG`, `TTT`,
impl GeneticCode for &[Amino; 64] {
    fn map_codon(&self, codon: [Nuc; 3]) -> Amino {
        let [n1, n2, n3] = codon;
        let idx = 16 * (n1 as u8).trailing_zeros()
            + 4 * (n2 as u8).trailing_zeros()
            + (n3 as u8).trailing_zeros();
        self[idx as usize]
    }
}

/// [`GeneticCode`] implementation optimized for speed over space.
///
/// Each [`FastTranslator`] comes with about 18KiB of additional lookup tables,
/// dramatically improving translation speed, particularly for ambigous codons.
pub struct FastTranslator {
    // `lookup` is for `Nuc`s and `ambi_lookup` is for `AmbiNuc`s. They have similar layouts;
    // each `Nuc`/`AmbiNuc` can be thought of as a hexidecimal digit, and each codon is thought of
    // as a big-endian 3-digit hexadecimal index into the lookup table. A `Nuc` caps out at 0x8,
    // but an `AmbiNuc` caps out at 0xf, hence the difference in table lengths. Note that this
    // layout still requires entries for indices containing the digit 0, so we need an extra
    // entry at the start for 0x000.
    lookup: [Amino; 1 + 0x888],          // about 2KiB in size
    ambi_lookup: [AmbiAmino; 1 + 0xfff], // 16KiB in size because `AmbiAmino` is 4 bytes.
}

impl FastTranslator {
    /// Build [`FastTranslator`] from a table of [`Amino`]s.
    ///
    /// The [`Amino`]s must correspond to all codons in ascending lexicographical order.
    ///
    /// This requires up-front work, but can be done at compile-time.
    #[must_use]
    pub const fn from_table(table: [Amino; 64]) -> Self {
        let lookup = Self::build_concrete_lookup(table);
        Self {
            ambi_lookup: Self::build_ambi_lookup(&lookup),
            lookup,
        }
    }

    const fn build_concrete_lookup(table: [Amino; 64]) -> [Amino; 1 + 0x888] {
        let mut lookup = [Amino::A; _];
        let mut i = 0;
        while i < table.len() {
            let n1 = 1 << (0b11 & (i >> 4));
            let n2 = 1 << (0b11 & (i >> 2));
            let n3 = 1 << (0b11 & i);
            lookup[(n1 << 8) | (n2 << 4) | n3] = table[i];
            i += 1;
        }
        lookup
    }

    const fn build_ambi_lookup(concrete_lookup: &[Amino; 1 + 0x888]) -> [AmbiAmino; 1 + 0xfff] {
        // This is 20x faster than calculating each lookup entry in isolation because it dedupes
        // the work of building merged table results. It treats the lookup like a lattice of sets:
        // Instead of having to build each set from scratch, it just builds them out of two
        // smaller sets. And counting upwards has the nice property of visiting the sets in
        // topological order so subsets will be populated before the sets that they feed into.
        // When an index contains multiple bits (when _b variables are non-zero), that means it's
        // ambiguous and should composed out of the results of smaller indices. In theory this
        // could be flattened into a single unnested loop that populates the indices in order,
        // but then the iterations couldn't share as much work, so it'd take about 50% longer.
        let mut lookup = [AmbiAmino::X; _];
        let mut n1 = 0x100;
        while n1 <= 0xf00 {
            let (n1_a, n1_b) = Self::split_lowest_bit(n1);
            if n1_b == 0 {
                let mut n2 = 0x010;
                while n2 <= 0x0f0 {
                    let n1_n2 = n1 | n2;
                    let (n2_a, n2_b) = Self::split_lowest_bit(n2);
                    if n2_b == 0 {
                        let mut n3 = 0x001;
                        while n3 <= 0x00f {
                            let n1_n2_n3 = n1_n2 | n3;
                            let (n3_a, n3_b) = Self::split_lowest_bit(n3);
                            if n3_b == 0 {
                                lookup[n1_n2_n3] = AmbiAmino::from_amino(concrete_lookup[n1_n2_n3]);
                            } else {
                                lookup[n1_n2_n3] = lookup[n1_n2 | n3_a].or(lookup[n1_n2 | n3_b]);
                            }
                            n3 += 0x001;
                        }
                    } else {
                        let n1_n2_a = n1 | n2_a;
                        let n1_n2_b = n1 | n2_b;
                        let mut n3 = 0x001;
                        while n3 <= 0x00f {
                            lookup[n1_n2 | n3] = lookup[n1_n2_a | n3].or(lookup[n1_n2_b | n3]);
                            n3 += 1;
                        }
                    }
                    n2 += 0x010;
                }
            } else {
                let mut n2_n3 = 0x011;
                while n2_n3 <= 0x0ff {
                    lookup[n1 | n2_n3] = lookup[n1_a | n2_n3].or(lookup[n1_b | n2_n3]);
                    n2_n3 += 1;
                }
            }
            n1 += 0x100;
        }
        lookup
    }

    // Split a number into its lowest bit and other bits.
    const fn split_lowest_bit(i: usize) -> (usize, usize) {
        let lowest = 1 << i.trailing_zeros();
        (lowest, i & !lowest)
    }
}

impl GeneticCode for &FastTranslator {
    #[inline(always)]
    fn map_codon(&self, codon: [Nuc; 3]) -> Amino {
        let [n1, n2, n3] = codon;
        self.lookup[((n1 as usize) << 8) | ((n2 as usize) << 4) | (n3 as usize)]
    }

    #[inline(always)]
    fn map_ambi_codon(&self, codon: [AmbiNuc; 3]) -> AmbiAmino {
        let [n1, n2, n3] = codon;
        self.ambi_lookup[((n1 as usize) << 8) | ((n2 as usize) << 4) | (n3 as usize)]
    }
}

/// Standard code
pub const NCBI1: &FastTranslator = &ncbi(0);
/// Vertebrate mitochondrial code
pub const NCBI2: &FastTranslator = &ncbi(1);
/// Yeast mitochondrial code
pub const NCBI3: &FastTranslator = &ncbi(2);
/// Mold, protozoan, and coelenterate mitochondrial code as well as mycoplasma and spiroplasma code
pub const NCBI4: &FastTranslator = &ncbi(3);
/// Invertebrate mitochondrial code
pub const NCBI5: &FastTranslator = &ncbi(4);
/// Ciliate dasycladacean and hexamita nuclear code
pub const NCBI6: &FastTranslator = &ncbi(5);
/// Echinoderm and flatworm mitochondrial code
pub const NCBI9: &FastTranslator = &ncbi(6);
/// Euplotid nuclear code
pub const NCBI10: &FastTranslator = &ncbi(7);
/// Bacterial, archaeal and plant plastid code
pub const NCBI11: &FastTranslator = &ncbi(8);
/// Alternative yeast nuclear code
pub const NCBI12: &FastTranslator = &ncbi(9);
/// Ascidian mitochondrial code
pub const NCBI13: &FastTranslator = &ncbi(10);
/// Alternative flatworm mitochondrial code
pub const NCBI14: &FastTranslator = &ncbi(11);
/// Blepharisma nuclear code
pub const NCBI15: &FastTranslator = &ncbi(12);
/// Chlorophycean mitochondrial code
pub const NCBI16: &FastTranslator = &ncbi(13);
/// Trematode mitochondrial code
pub const NCBI21: &FastTranslator = &ncbi(14);
/// Scenedesmus obliquus mitochondrial code
pub const NCBI22: &FastTranslator = &ncbi(15);
/// Thraustochytrium mitochondrial code
pub const NCBI23: &FastTranslator = &ncbi(16);
/// Pterobranchia mitochondrial code
pub const NCBI24: &FastTranslator = &ncbi(17);
/// Candidate division SR1 and gracilibacteria code
pub const NCBI25: &FastTranslator = &ncbi(18);
/// Pachysolen tannophilus nuclear code
pub const NCBI26: &FastTranslator = &ncbi(19);
/// Karyorelict nuclear code
pub const NCBI27: &FastTranslator = &ncbi(20);
/// Condylostoma nuclear code
pub const NCBI28: &FastTranslator = &ncbi(21);
/// Mesodinium nuclear code
pub const NCBI29: &FastTranslator = &ncbi(22);
/// Peritrich nuclear code
pub const NCBI30: &FastTranslator = &ncbi(23);
/// Blastocrithidia nuclear code
pub const NCBI31: &FastTranslator = &ncbi(24);
/// Balanophoraceae plastid code
pub const NCBI32: &FastTranslator = &ncbi(25);
/// Cephalodiscidae mitochondrial code
pub const NCBI33: &FastTranslator = &ncbi(26);
/// Enterosoma code
pub const NCBI34: &FastTranslator = &ncbi(27);
/// Peptacetobacter code
pub const NCBI35: &FastTranslator = &ncbi(28);
/// Anaerococcus and onthovivens code
pub const NCBI36: &FastTranslator = &ncbi(29);
/// Absconditabacterales genetic code
pub const NCBI37: &FastTranslator = &ncbi(30);

const fn ncbi(i: usize) -> FastTranslator {
    FastTranslator::from_table(NCBI_DATA[i])
}

macro_rules! aminos {
    [ $( [$( $c:tt )+ ]),+ $(,)? ] => {{
        [$(
            [$(to_amino!($c)),+]
        ),+]
    }}
}

macro_rules! to_amino {
    (*) => {{ Amino::Stop }};
    ($v:tt) => {{ Amino::$v }};
}

const NCBI_DATA: [[Amino; 64]; 31] = transpose(aminos! [
    // Tables are obtained from: https://en.wikipedia.org/wiki/List_of_genetic_codes
    //         1  2  3  4  5  6  9 10 11 12 13 14 15 16 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37
    /* AAA */ [K  K  K  K  K  K  N  K  K  K  K  N  K  K  N  K  K  K  K  K  K  K  K  K  K  K  K  K  K  K  K],
    /* AAC */ [N  N  N  N  N  N  N  N  N  N  N  N  N  N  N  N  N  N  N  N  N  N  N  N  N  N  N  N  N  N  N],
    /* AAG */ [K  K  K  K  K  K  K  K  K  K  K  K  K  K  K  K  K  K  K  K  K  K  K  K  K  K  K  K  K  K  K],
    /* AAT */ [N  N  N  N  N  N  N  N  N  N  N  N  N  N  N  N  N  N  N  N  N  N  N  N  N  N  N  N  N  N  N],
    /* ACA */ [T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T],
    /* ACC */ [T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T],
    /* ACG */ [T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T],
    /* ACT */ [T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T  T],
    /* AGA */ [R  *  R  R  S  R  S  R  R  R  G  S  R  R  S  R  R  S  R  R  R  R  R  R  R  R  S  R  R  R  R],
    /* AGC */ [S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S],
    /* AGG */ [R  *  R  R  S  R  S  R  R  R  G  S  R  R  S  R  R  K  R  R  R  R  R  R  R  R  K  M  R  R  R],
    /* AGT */ [S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S],
    /* ATA */ [I  M  M  I  M  I  I  I  I  I  M  I  I  I  M  I  I  I  I  I  I  I  I  I  I  I  I  I  I  I  I],
    /* ATC */ [I  I  I  I  I  I  I  I  I  I  I  I  I  I  I  I  I  I  I  I  I  I  I  I  I  I  I  I  I  I  I],
    /* ATG */ [M  M  M  M  M  M  M  M  M  M  M  M  M  M  M  M  M  M  M  M  M  M  M  M  M  M  M  M  M  M  M],
    /* ATT */ [I  I  I  I  I  I  I  I  I  I  I  I  I  I  I  I  I  I  I  I  I  I  I  I  I  I  I  I  I  I  I],
    /* CAA */ [Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q],
    /* CAC */ [H  H  H  H  H  H  H  H  H  H  H  H  H  H  H  H  H  H  H  H  H  H  H  H  H  H  H  H  H  H  H],
    /* CAG */ [Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q  Q],
    /* CAT */ [H  H  H  H  H  H  H  H  H  H  H  H  H  H  H  H  H  H  H  H  H  H  H  H  H  H  H  H  H  H  H],
    /* CCA */ [P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P],
    /* CCC */ [P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P],
    /* CCG */ [P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P],
    /* CCT */ [P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P  P],
    /* CGA */ [R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  W],
    /* CGC */ [R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R],
    /* CGG */ [R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  Q  W  W],
    /* CGT */ [R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R  R],
    /* CTA */ [L  L  T  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L],
    /* CTC */ [L  L  T  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L],
    /* CTG */ [L  L  T  L  L  L  L  L  L  S  L  L  L  L  L  L  L  L  L  A  L  L  L  L  L  L  L  L  L  L  L],
    /* CTT */ [L  L  T  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L],
    /* GAA */ [E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E],
    /* GAC */ [D  D  D  D  D  D  D  D  D  D  D  D  D  D  D  D  D  D  D  D  D  D  D  D  D  D  D  D  D  D  D],
    /* GAG */ [E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E],
    /* GAT */ [D  D  D  D  D  D  D  D  D  D  D  D  D  D  D  D  D  D  D  D  D  D  D  D  D  D  D  D  D  D  D],
    /* GCA */ [A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A],
    /* GCC */ [A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A],
    /* GCG */ [A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A],
    /* GCT */ [A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A  A],
    /* GGA */ [G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G],
    /* GGC */ [G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G],
    /* GGG */ [G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G],
    /* GGT */ [G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G  G],
    /* GTA */ [V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V],
    /* GTC */ [V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V],
    /* GTG */ [V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V],
    /* GTT */ [V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V  V],
    /* TAA */ [*  *  *  *  *  Q  *  *  *  *  *  Y  *  *  *  *  *  *  *  *  Q  Q  Y  E  E  *  Y  *  *  *  *],
    /* TAC */ [Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y],
    /* TAG */ [*  *  *  *  *  Q  *  *  *  *  *  *  Q  L  *  L  *  *  *  *  Q  Q  Y  E  E  W  *  *  *  *  *],
    /* TAT */ [Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y  Y],
    /* TCA */ [S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  *  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S],
    /* TCC */ [S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S],
    /* TCG */ [S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S],
    /* TCT */ [S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S  S],
    /* TGA */ [*  W  W  W  W  *  W  C  *  *  W  W  *  *  W  *  *  W  G  *  W  W  *  *  W  *  W  *  *  *  G],
    /* TGC */ [C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C],
    /* TGG */ [W  W  W  W  W  W  W  W  W  W  W  W  W  W  W  W  W  W  W  W  W  W  W  W  W  W  W  W  W  W  W],
    /* TGT */ [C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C  C],
    /* TTA */ [L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  *  L  L  L  L  L  L  L  L  L  L  L  L  L  L],
    /* TTC */ [F  F  F  F  F  F  F  F  F  F  F  F  F  F  F  F  F  F  F  F  F  F  F  F  F  F  F  F  F  F  F],
    /* TTG */ [L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L  L],
    /* TTT */ [F  F  F  F  F  F  F  F  F  F  F  F  F  F  F  F  F  F  F  F  F  F  F  F  F  F  F  F  F  F  F],
]);

const fn transpose<const N: usize, const M: usize>(table: [[Amino; N]; M]) -> [[Amino; M]; N] {
    let mut output = [[Amino::A; M]; N];
    let mut i = 0;
    while i < N {
        let mut j = 0;
        while j < M {
            output[i][j] = table[j][i];
            j += 1;
        }
        i += 1;
    }
    output
}

#[cfg(test)]
mod tests {
    use crate::Nucleotide;

    use super::*;

    #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
    #[test]
    fn compare_fast_lookup_against_reference_implementation() {
        // In order to be able to build `FastLookup`s for all NCBI tables at compile-time
        // in 1 second, the ambiguous lookup generation code is convoluted. (doing things the
        // obvious way would result in `nucs` taking 30 seconds to compile) To improve confidence,
        // we check that for all codons and NCBI tables, `FastLookup` gives the same results as
        // the simpler implementations provided by `&[Amino; 64]` and `GeneticCode`.
        for raw in &NCBI_DATA {
            let fast = &FastTranslator::from_table(*raw);
            for n1 in Nuc::ALL {
                for n2 in Nuc::ALL {
                    for n3 in Nuc::ALL {
                        let codon = [n1, n2, n3];
                        assert_eq!(Nuc::translate(&fast, codon), Nuc::translate(&raw, codon));
                    }
                }
            }
            for n1 in AmbiNuc::ALL {
                for n2 in AmbiNuc::ALL {
                    for n3 in AmbiNuc::ALL {
                        let codon = [n1, n2, n3];
                        assert_eq!(
                            AmbiNuc::translate(&fast, codon),
                            AmbiNuc::translate(&raw, codon)
                        );
                    }
                }
            }
        }
    }
}
