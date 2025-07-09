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

// TODO: type representing a pair of precomputed tables? (one for concrete, one for ambiguous)

/// Standard code
pub const NCBI1: &[Amino; 64] = &NCBI_DATA[0];
/// Vertebrate mitochondrial code
pub const NCBI2: &[Amino; 64] = &NCBI_DATA[1];
/// Yeast mitochondrial code
pub const NCBI3: &[Amino; 64] = &NCBI_DATA[2];
/// Mold, protozoan, and coelenterate mitochondrial code as well as mycoplasma and spiroplasma code
pub const NCBI4: &[Amino; 64] = &NCBI_DATA[3];
/// Invertebrate mitochondrial code
pub const NCBI5: &[Amino; 64] = &NCBI_DATA[4];
/// Ciliate dasycladacean and hexamita nuclear code
pub const NCBI6: &[Amino; 64] = &NCBI_DATA[5];
/// Echinoderm and flatworm mitochondrial code
pub const NCBI9: &[Amino; 64] = &NCBI_DATA[6];
/// Euplotid nuclear code
pub const NCBI10: &[Amino; 64] = &NCBI_DATA[7];
/// Bacterial, archaeal and plant plastid code
pub const NCBI11: &[Amino; 64] = &NCBI_DATA[8];
/// Alternative yeast nuclear code
pub const NCBI12: &[Amino; 64] = &NCBI_DATA[9];
/// Ascidian mitochondrial code
pub const NCBI13: &[Amino; 64] = &NCBI_DATA[10];
/// Alternative flatworm mitochondrial code
pub const NCBI14: &[Amino; 64] = &NCBI_DATA[11];
/// Blepharisma nuclear code
pub const NCBI15: &[Amino; 64] = &NCBI_DATA[12];
/// Chlorophycean mitochondrial code
pub const NCBI16: &[Amino; 64] = &NCBI_DATA[13];
/// Trematode mitochondrial code
pub const NCBI21: &[Amino; 64] = &NCBI_DATA[14];
/// Scenedesmus obliquus mitochondrial code
pub const NCBI22: &[Amino; 64] = &NCBI_DATA[15];
/// Thraustochytrium mitochondrial code
pub const NCBI23: &[Amino; 64] = &NCBI_DATA[16];
/// Pterobranchia mitochondrial code
pub const NCBI24: &[Amino; 64] = &NCBI_DATA[17];
/// Candidate division SR1 and gracilibacteria code
pub const NCBI25: &[Amino; 64] = &NCBI_DATA[18];
/// Pachysolen tannophilus nuclear code
pub const NCBI26: &[Amino; 64] = &NCBI_DATA[19];
/// Karyorelict nuclear code
pub const NCBI27: &[Amino; 64] = &NCBI_DATA[20];
/// Condylostoma nuclear code
pub const NCBI28: &[Amino; 64] = &NCBI_DATA[21];
/// Mesodinium nuclear code
pub const NCBI29: &[Amino; 64] = &NCBI_DATA[22];
/// Peritrich nuclear code
pub const NCBI30: &[Amino; 64] = &NCBI_DATA[23];
/// Blastocrithidia nuclear code
pub const NCBI31: &[Amino; 64] = &NCBI_DATA[24];
/// Balanophoraceae plastid code
pub const NCBI32: &[Amino; 64] = &NCBI_DATA[25];
/// Cephalodiscidae mitochondrial code
pub const NCBI33: &[Amino; 64] = &NCBI_DATA[26];
/// Enterosoma code
pub const NCBI34: &[Amino; 64] = &NCBI_DATA[27];
/// Peptacetobacter code
pub const NCBI35: &[Amino; 64] = &NCBI_DATA[28];
/// Anaerococcus and onthovivens code
pub const NCBI36: &[Amino; 64] = &NCBI_DATA[29];
/// Absconditabacterales genetic code
pub const NCBI37: &[Amino; 64] = &NCBI_DATA[30];

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
