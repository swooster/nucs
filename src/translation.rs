//! Types related to translation of codons into amino acids.

use crate::{AmbiAmino, AmbiNuc, Amino, DnaSlice, Nuc, Nucleotide};

/// Trait representing any type that can be used to translate codons into amino acids.
///
/// The most commonly used [`GeneticCode`] is [`NCBI1`] but there are others in
/// [`nucs::translation`](crate::translation).
pub trait GeneticCode {
    /// Translate a codon to an amino acid
    ///
    /// Examples
    ///
    /// ```
    /// use nucs::{AmbiAmino, Amino, AmbiNuc, Nuc, translation::GeneticCode};
    ///
    /// let genetic_code = |_codon| Amino::W; // translate everything to tryptophan
    /// assert_eq!(genetic_code.translate([Nuc::A; 3]), Amino::W);
    /// assert_eq!(genetic_code.translate([AmbiNuc::N; 3]), AmbiAmino::W);
    /// ```
    fn translate<N: Nucleotide>(&self, codon: [N; 3]) -> N::Amino {
        N::translate(self, codon)
    }

    /// Translate reverse complement of a codon to an amino acid
    ///
    /// Examples
    ///
    /// ```
    /// use nucs::{AmbiAmino, Amino, AmbiNuc, Nuc, translation::GeneticCode};
    /// use Nuc::{A, G, T};
    ///
    /// let genetic_code = |codon: [Nuc; 3]| match codon {
    ///    [A, T, G] => Amino::W,
    ///    _ => Amino::Stop,
    /// };
    /// assert_eq!(genetic_code.translate_rc(Nuc::lit(b"CAT")), Amino::W);
    /// assert_eq!(genetic_code.translate_rc(AmbiNuc::lit(b"CAT")), AmbiAmino::W);
    /// assert_eq!(genetic_code.translate_rc(AmbiNuc::lit(b"ANT")), AmbiAmino::Stop);
    /// ```
    fn translate_rc<N: Nucleotide>(&self, codon: [N; 3]) -> N::Amino {
        N::translate_rc(self, codon)
    }

    /// Translate a concrete codon to an amino acid
    ///
    /// Consider calling [`GeneticCode::translate`] instead; this primarily exists to be provided
    /// by implementors of [`GeneticCode`].
    fn translate_concrete_codon(&self, codon: [Nuc; 3]) -> Amino;

    /// Translate an ambiguous codon to an amino acid
    ///
    /// Consider calling [`GeneticCode::translate`] instead; this primarily exists to be provided
    /// by implementors of [`GeneticCode`].
    fn translate_ambiguous_codon(&self, codon: [AmbiNuc; 3]) -> AmbiAmino {
        let [ambi_n1, ambi_n2, ambi_n3] = codon;
        ambi_n1
            .iter()
            .flat_map(|n1| ambi_n2.iter().map(move |n2| [n1, n2]))
            .flat_map(|[n1, n2]| ambi_n3.iter().map(move |n3| [n1, n2, n3]))
            .map(|codon| AmbiAmino::from(self.translate_concrete_codon(codon)))
            .reduce(|a, b| a | b)
            .expect("BUG: null nucleotide encountered")
    }

    /// Translate reverse complement of a concrete codon to an amino acid
    ///
    /// Consider calling [`GeneticCode::translate_rc`] instead; this primarily exists to be provided
    /// by implementors of [`GeneticCode`].
    fn translate_rc_concrete_codon(&self, mut codon: [Nuc; 3]) -> Amino {
        codon.reverse();
        codon.complement();
        self.translate_concrete_codon(codon)
    }

    /// Translate reverse complement of an ambiguous codon to an amino acid
    ///
    /// Consider calling [`GeneticCode::translate_rc`] instead; this primarily exists to be provided
    /// by implementors of [`GeneticCode`].
    fn translate_rc_ambiguous_codon(&self, mut codon: [AmbiNuc; 3]) -> AmbiAmino {
        codon.reverse();
        codon.complement();
        self.translate_ambiguous_codon(codon)
    }
}

impl<F: Fn([Nuc; 3]) -> Amino> GeneticCode for F {
    fn translate_concrete_codon(&self, codon: [Nuc; 3]) -> Amino {
        self(codon)
    }
}

/// A table of amino acids can be used for genetic coding.
///
/// The [`Amino`]s must be ordered to correspond with codons in ascending lexicographical order
/// (`AAA`, `AAC`, `AAG`, `AAT`, `ACA`, ... `TTT`).
impl GeneticCode for &[Amino; 64] {
    fn translate_concrete_codon(&self, codon: [Nuc; 3]) -> Amino {
        let [n1, n2, n3] = codon;
        let idx = 16 * (n1 as u8).trailing_zeros()
            + 4 * (n2 as u8).trailing_zeros()
            + (n3 as u8).trailing_zeros();
        self[idx as usize]
    }
}

/// [`GeneticCode`] implementation optimized for speed over space.
///
/// Each [`FullLookup`] comes with about 18KiB of additional lookup tables,
/// dramatically improving translation speed, particularly for ambigous codons.
#[derive(Clone)]
pub struct FullLookup {
    lookup: ConcreteLookup,
    lookup_rc: ConcreteLookup,
    ambi_lookup: AmbiLookup,
    ambi_lookup_rc: AmbiLookup,
}

impl FullLookup {
    /// Build [`FullLookup`] from another [`GeneticCode`].
    ///
    /// This precalculates every possible ambiguous codon's translation.
    #[must_use]
    pub fn from_genetic_code<G: GeneticCode>(genetic_code: &G) -> Self {
        Self::from_table(&std::array::from_fn(|i| {
            let codon = [4, 2, 0].map(|offset| Nuc::ALL[(i >> offset) & 0b11]);
            genetic_code.translate(codon)
        }))
    }

    /// Build [`FullLookup`] from a table of [`Amino`]s.
    ///
    /// The [`Amino`]s must correspond to all codons in ascending lexicographical order.
    ///
    /// This precalculates every possible ambiguous codon's translation,
    /// but can be done at compile-time.
    #[must_use]
    pub const fn from_table(table: &[Amino; 64]) -> Self {
        let lookup = ConcreteLookup::from_table(table);
        let lookup_rc = lookup.reverse_complement();
        let ambi_lookup = lookup.to_ambi_lookup();
        let ambi_lookup_rc = ambi_lookup.reverse_complement();
        Self {
            lookup,
            lookup_rc,
            ambi_lookup,
            ambi_lookup_rc,
        }
    }

    /// Convert to simple amino acid table.
    ///
    /// The returned [`Amino`]s correspond to all codons in ascending lexicographical order.
    /// Although such a table can be used as a [`GeneticCode`], it's intended more for interop
    /// or inspection.
    ///
    /// ```
    /// use nucs::{Amino, NCBI1};
    ///
    /// let table = NCBI1.to_table();
    /// assert_eq!(table[0], Amino::K);  // AAA -> K
    /// assert_eq!(table[1], Amino::N);  // AAC -> N
    /// assert_eq!(table[2], Amino::K);  // AAG -> K
    /// assert_eq!(table[3], Amino::N);  // AAT -> N
    /// assert_eq!(table[4], Amino::T);  // ACA -> T
    /// // ...
    /// assert_eq!(table[16], Amino::Q); // CAA -> Q
    /// // ...
    /// assert_eq!(table[32], Amino::E); // GAA -> E
    /// // ...
    /// assert_eq!(table[62], Amino::L); // TTG -> L
    /// assert_eq!(table[63], Amino::F); // TTT -> F
    /// ```
    #[must_use]
    pub const fn to_table(&self) -> [Amino; 64] {
        self.lookup.to_table()
    }
}

impl GeneticCode for &FullLookup {
    #[inline(always)]
    fn translate_concrete_codon(&self, codon: [Nuc; 3]) -> Amino {
        (&self.lookup).translate_concrete_codon(codon)
    }

    #[inline(always)]
    fn translate_ambiguous_codon(&self, codon: [AmbiNuc; 3]) -> AmbiAmino {
        (&self.ambi_lookup).translate_ambiguous_codon(codon)
    }

    #[inline(always)]
    fn translate_rc_concrete_codon(&self, codon: [Nuc; 3]) -> Amino {
        (&self.lookup_rc).translate_concrete_codon(codon)
    }

    #[inline(always)]
    fn translate_rc_ambiguous_codon(&self, codon: [AmbiNuc; 3]) -> AmbiAmino {
        (&self.ambi_lookup_rc).translate_ambiguous_codon(codon)
    }
}

impl std::fmt::Debug for FullLookup {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.lookup.fmt(f)
    }
}

/// A lookup to efficiently translate concrete codons
///
/// The lookup doesn't allocate, but is about 2KiB in size. It supports very fast translation of
/// concrete codons (`[Nuc; 3]`) but is slower at translating ambiguous codons, or
/// reverse-complement codons.
#[derive(Clone)]
pub struct ConcreteLookup([Amino; 1 + 0x888]);
// The table layout:
// Each `Nuc` can be thought of as a hexidecimal digit, and each codon as a big-endian
// 3-digit hexadecimal index into the lookup table. A `Nuc` caps out at 0x8, so TTT
// corresponds to 0x888, plus an extra spot is needed at the start for 0x000.

impl ConcreteLookup {
    /// Build [`ConcreteLookup`] from a table of [`Amino`]s.
    ///
    /// The [`Amino`]s must correspond to all codons in ascending lexicographical order.
    #[must_use]
    pub const fn from_table(table: &[Amino; 64]) -> Self {
        let mut lookup = [Amino::A; _];
        let mut i = 0;
        while i < table.len() {
            lookup[Self::table_index_to_fast_lookup_index(i)] = table[i];
            i += 1;
        }
        Self(lookup)
    }

    /// Convert to simple amino acid table.
    ///
    /// The returned [`Amino`]s correspond to all codons in ascending lexicographical order.
    /// Although such a table can be used as a [`GeneticCode`], it's intended more for interop
    /// or inspection.
    ///
    /// ```
    /// use nucs::{Amino, NCBI1};
    /// use nucs::translation::ConcreteLookup;
    ///
    /// let table = Amino::lit(b"ACDEFGHIKLMNOPQRSTUVWY*ACDEFGHIKLMNOPQRSTUVWY*ACDEFGHIKLMNOPQRS*");
    /// let lookup = ConcreteLookup::from_table(&table);
    /// assert_eq!(lookup.to_table(), table);
    /// ```
    #[must_use]
    pub const fn to_table(&self) -> [Amino; 64] {
        let mut table = [Amino::A; _];
        let mut i = 0;
        while i < table.len() {
            table[i] = self.0[Self::table_index_to_fast_lookup_index(i)];
            i += 1;
        }
        table
    }

    // (by "table", I mean list of one amino per concrete codon, in lexicographical order)
    const fn table_index_to_fast_lookup_index(idx: usize) -> usize {
        let n1 = 1 << (0b11 & (idx >> 4));
        let n2 = 1 << (0b11 & (idx >> 2));
        let n3 = 1 << (0b11 & idx);
        (n1 << 8) | (n2 << 4) | n3
    }

    /// Return the reverse complement of this translation table.
    ///
    /// It produces a copy of the translation table where the reverse complement of each codon
    /// maps to the same amino acid as before. For example:
    ///
    /// ```
    /// use nucs::{Amino, Nuc};
    /// use nucs::translation::{ConcreteLookup, GeneticCode};
    /// use Nuc::{A, C, G, T};
    ///
    /// let table = Amino::lit(b"ACDEFGHIKLMNOPQRSTUVWY*ACDEFGHIKLMNOPQRSTUVWY*ACDEFGHIKLMNOPQRS*");
    /// let lookup = ConcreteLookup::from_table(&table);
    /// assert_eq!((&lookup).translate([A, T, G]), Amino::Q);
    /// let lookup_rc = lookup.reverse_complement();
    /// // ATG reverse-complemented is CAT, so...
    /// assert_eq!((&lookup_rc).translate([C, A, T]), Amino::Q);
    /// ```
    #[must_use]
    pub const fn reverse_complement(&self) -> Self {
        let mut lookup = [Amino::A; _];
        // As before, this could be flattened into a single unnested loop, but this nested
        // loop results in much sparser work.
        let mut n1 = 0x100;
        while n1 <= 0x800 {
            let mut n2 = 0x010;
            while n2 <= 0x080 {
                let mut n3 = 0x001;
                while n3 <= 0x008 {
                    let n1_n2_n3 = n1 | n2 | n3;
                    lookup[n1_n2_n3] = self.0[reverse_complement_index(n1_n2_n3)];
                    n3 <<= 1;
                }
                n2 <<= 1;
            }
            n1 <<= 1;
        }
        Self(lookup)
    }

    /// Build [`AmbiLookup`] from this [`ConcreteLookup`].
    ///
    /// This precalculates every possible ambiguous codon's translation,
    /// but can be done at compile-time.
    #[must_use]
    pub const fn to_ambi_lookup(&self) -> AmbiLookup {
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
                                lookup[n1_n2_n3] = AmbiAmino::from_amino(self.0[n1_n2_n3]);
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
        AmbiLookup(lookup)
    }

    // Split a number into its lowest bit and other bits.
    const fn split_lowest_bit(i: usize) -> (usize, usize) {
        let lowest = 1 << i.trailing_zeros();
        (lowest, i & !lowest)
    }
}

impl GeneticCode for &ConcreteLookup {
    #[inline(always)]
    fn translate_concrete_codon(&self, codon: [Nuc; 3]) -> Amino {
        let [n1, n2, n3] = codon;
        self.0[((n1 as usize) << 8) | ((n2 as usize) << 4) | (n3 as usize)]
    }
}

impl std::fmt::Debug for ConcreteLookup {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        fmt_genetic_code(&self, f)
    }
}

/// A lookup to efficiently translate ambiguous codons
///
/// The lookup doesn't allocate, but is about 2KiB in size. It supports very fast translation of
/// ambiguous codons (`[AmbiNuc; 3]`) but is slower at translating concrete codons, or
/// reverse-complement codons.
#[derive(Clone)]
pub struct AmbiLookup([AmbiAmino; 1 + 0xfff]);
// The table layout is similar to `ConcreteLookup`:
// Each `AmbiNuc` can be thought of as a hexidecimal digit, and each codon as a big-endian
// 3-digit hexadecimal index into the lookup table. An `AmbiNuc` caps out at 0xf, so NNN
// corresponds to 0xfff, plus an extra spot is needed at the start for 0x000.

impl AmbiLookup {
    /// Return the reverse complement of this translation table.
    ///
    /// It produces a copy of the translation table where the reverse complement of each codon
    /// maps to the same amino acid as before. For example:
    ///
    /// ```
    /// use nucs::{Amino, Nuc, NCBI1};
    /// use nucs::translation::{ConcreteLookup, GeneticCode};
    /// use Nuc::{A, C, G, T};
    ///
    /// let table = Amino::lit(b"ACDEFGHIKLMNOPQRSTUVWY*ACDEFGHIKLMNOPQRSTUVWY*ACDEFGHIKLMNOPQRS*");
    /// let lookup = ConcreteLookup::from_table(&table).to_ambi_lookup();
    /// assert_eq!((&lookup).translate([A, T, G]), Amino::Q);
    /// let lookup_rc = lookup.reverse_complement();
    /// // ATG reverse-complemented is CAT, so...
    /// assert_eq!((&lookup_rc).translate([C, A, T]), Amino::Q);
    /// ```
    #[must_use]
    pub const fn reverse_complement(&self) -> Self {
        let mut ambi_lookup = [AmbiAmino::X; _];
        // As before, this could be flattened into a single unnested loop, but this nested
        // loop results in much sparser work.
        let mut n1 = 0x100;
        while n1 <= 0xf00 {
            let mut n2 = 0x010;
            while n2 <= 0x0f0 {
                let mut n3 = 0x001;
                while n3 <= 0x00f {
                    let n1_n2_n3 = n1 | n2 | n3;
                    ambi_lookup[n1_n2_n3] = self.0[reverse_complement_index(n1_n2_n3)];
                    n3 <<= 1;
                }
                n2 <<= 1;
            }
            n1 <<= 1;
        }
        Self(ambi_lookup)
    }
}

impl GeneticCode for &AmbiLookup {
    #[inline(always)]
    fn translate_concrete_codon(&self, codon: [Nuc; 3]) -> Amino {
        let [n1, n2, n3] = codon;
        let amino = self.translate_ambiguous_codon([n1.into(), n2.into(), n3.into()]);
        amino
            .try_into()
            .expect("BUG: to_ambi_lookup mapped concrete codon to ambi amino")
    }

    #[inline(always)]
    fn translate_ambiguous_codon(&self, codon: [AmbiNuc; 3]) -> AmbiAmino {
        let [n1, n2, n3] = codon;
        self.0[((n1 as usize) << 8) | ((n2 as usize) << 4) | (n3 as usize)]
    }
}

impl std::fmt::Debug for AmbiLookup {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        fmt_genetic_code(&self, f)
    }
}

#[allow(
    clippy::cast_possible_truncation,
    reason = "valid lookup indexes are 12 bits"
)]
const fn reverse_complement_index(idx: usize) -> usize {
    ((idx as u16).reverse_bits() >> 4) as usize
}

fn fmt_genetic_code(
    genetic_code: &impl GeneticCode,
    f: &mut std::fmt::Formatter<'_>,
) -> std::fmt::Result {
    let mut map = f.debug_map();
    for n1 in Nuc::ALL {
        for n2 in Nuc::ALL {
            for n3 in Nuc::ALL {
                let codon = [n1, n2, n3];
                map.entry(&codon.display(), &genetic_code.translate([n1, n2, n3]));
            }
        }
    }
    map.finish()
}

/// Standard code
pub const NCBI1: &FullLookup = &ncbi(0);
/// Vertebrate mitochondrial code
pub const NCBI2: &FullLookup = &ncbi(1);
/// Yeast mitochondrial code
pub const NCBI3: &FullLookup = &ncbi(2);
/// Mold, protozoan, and coelenterate mitochondrial code as well as mycoplasma and spiroplasma code
pub const NCBI4: &FullLookup = &ncbi(3);
/// Invertebrate mitochondrial code
pub const NCBI5: &FullLookup = &ncbi(4);
/// Ciliate dasycladacean and hexamita nuclear code
pub const NCBI6: &FullLookup = &ncbi(5);
/// Echinoderm and flatworm mitochondrial code
pub const NCBI9: &FullLookup = &ncbi(6);
/// Euplotid nuclear code
pub const NCBI10: &FullLookup = &ncbi(7);
/// Bacterial, archaeal and plant plastid code
pub const NCBI11: &FullLookup = &ncbi(8);
/// Alternative yeast nuclear code
pub const NCBI12: &FullLookup = &ncbi(9);
/// Ascidian mitochondrial code
pub const NCBI13: &FullLookup = &ncbi(10);
/// Alternative flatworm mitochondrial code
pub const NCBI14: &FullLookup = &ncbi(11);
/// Blepharisma nuclear code
pub const NCBI15: &FullLookup = &ncbi(12);
/// Chlorophycean mitochondrial code
pub const NCBI16: &FullLookup = &ncbi(13);
/// Trematode mitochondrial code
pub const NCBI21: &FullLookup = &ncbi(14);
/// Scenedesmus obliquus mitochondrial code
pub const NCBI22: &FullLookup = &ncbi(15);
/// Thraustochytrium mitochondrial code
pub const NCBI23: &FullLookup = &ncbi(16);
/// Pterobranchia mitochondrial code
pub const NCBI24: &FullLookup = &ncbi(17);
/// Candidate division SR1 and gracilibacteria code
pub const NCBI25: &FullLookup = &ncbi(18);
/// Pachysolen tannophilus nuclear code
pub const NCBI26: &FullLookup = &ncbi(19);
/// Karyorelict nuclear code
pub const NCBI27: &FullLookup = &ncbi(20);
/// Condylostoma nuclear code
pub const NCBI28: &FullLookup = &ncbi(21);
/// Mesodinium nuclear code
pub const NCBI29: &FullLookup = &ncbi(22);
/// Peritrich nuclear code
pub const NCBI30: &FullLookup = &ncbi(23);
/// Blastocrithidia nuclear code
pub const NCBI31: &FullLookup = &ncbi(24);
/// Balanophoraceae plastid code
pub const NCBI32: &FullLookup = &ncbi(25);
/// Cephalodiscidae mitochondrial code
pub const NCBI33: &FullLookup = &ncbi(26);
/// Enterosoma code
pub const NCBI34: &FullLookup = &ncbi(27);
/// Peptacetobacter code
pub const NCBI35: &FullLookup = &ncbi(28);
/// Anaerococcus and onthovivens code
pub const NCBI36: &FullLookup = &ncbi(29);
/// Absconditabacterales genetic code
pub const NCBI37: &FullLookup = &ncbi(30);

const fn ncbi(i: usize) -> FullLookup {
    FullLookup::from_table(&NCBI_DATA[i])
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
    use super::*;

    #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
    #[test]
    fn compare_fast_lookup_against_reference_implementation() {
        // In order to be able to build `FullLookup`s for all NCBI tables at compile-time
        // in 1 second, the ambiguous lookup generation code is convoluted. (doing things the
        // obvious way would result in `nucs` taking 30 seconds to compile) To improve confidence,
        // we check that for all codons and NCBI tables, `FullLookup` gives the same results as
        // the simpler implementations provided by `&[Amino; 64]` and `GeneticCode`.
        for raw in &NCBI_DATA {
            let fast = &FullLookup::from_table(raw);
            assert_genetic_codes_eq(&fast, &raw);
        }
    }

    #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
    #[test]
    fn fast_lookup_can_be_built_from_other_genetic_code() {
        let new = &FullLookup::from_genetic_code(&NCBI1);
        assert_genetic_codes_eq(&new, &NCBI1);
    }

    #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
    #[test]
    fn fast_lookup_conversion_roundtrips() {
        let table = &NCBI_DATA[0];
        assert_eq!(&FullLookup::from_table(table).to_table(), table);
    }

    #[test]
    fn fast_lookup_reverse_complement() {
        let dna = Nuc::lit(b"ATCTTCGGGGGGAATTAAAAACTAATAAAGTTCAACAATGGTTGGCATCTCTTCCCGGGG");
        let peptide = dna.rc_translated_to_vec_by(NCBI1);
        assert_eq!(peptide, Amino::lit(b"PREEMPTIVELY*FLIPPED"));
    }

    #[test]
    fn exhaustively_test_reverse_complement_concrete_lookups() {
        for n1 in Nuc::ALL {
            for n2 in Nuc::ALL {
                for n3 in Nuc::ALL {
                    let codon = [n1, n2, n3];
                    let mut codon_rc = codon;
                    codon_rc.revcomp();
                    assert_eq!(
                        NCBI1.translate_rc(codon),
                        NCBI1.translate(codon_rc),
                        "Mismatch for {codon:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn exhaustively_test_reverse_complement_ambiguous_lookups() {
        for n1 in AmbiNuc::ALL {
            for n2 in AmbiNuc::ALL {
                for n3 in AmbiNuc::ALL {
                    let codon = [n1, n2, n3];
                    let mut codon_rc = codon;
                    codon_rc.revcomp();
                    assert_eq!(
                        NCBI1.translate_rc(codon),
                        NCBI1.translate(codon_rc),
                        "Mismatch for {codon:?}"
                    );
                }
            }
        }
    }

    fn assert_genetic_codes_eq(g1: &impl GeneticCode, g2: &impl GeneticCode) {
        for n1 in Nuc::ALL {
            for n2 in Nuc::ALL {
                for n3 in Nuc::ALL {
                    let codon = [n1, n2, n3];
                    assert_eq!(g1.translate(codon), g2.translate(codon));
                }
            }
        }
        for n1 in AmbiNuc::ALL {
            for n2 in AmbiNuc::ALL {
                for n3 in AmbiNuc::ALL {
                    let codon = [n1, n2, n3];
                    assert_eq!(g1.translate(codon), g2.translate(codon));
                }
            }
        }
    }
}
