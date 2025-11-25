//! Slice-related types

use crate::iter::{Codons, Translated};
use crate::translation::GeneticCode;
use crate::{DnaIter, Nucleotide};

#[cfg(feature = "unsafe")]
use crate::symbol::sealed::Sealed;
#[cfg(feature = "unsafe")]
use crate::{AmbiNuc, Nuc};

/// Helpers for working with slices of [`Nucleotide`]s.
pub trait DnaSlice {
    /// The type of [`Nucleotide`] this slice is made of.
    type Nuc: Nucleotide;

    /// Cast to slice of nucleotides
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{DnaSlice, Nuc};
    /// use Nuc::{A, C, G, T};
    ///
    /// let codons = [
    ///     [C, A, T],
    ///     [T, A, G],
    ///     [A, C, T],
    /// ];
    /// assert_eq!(codons.as_flat_dna(), Nuc::lit(b"CATTAGACT"));
    /// ```
    fn as_flat_dna(&self) -> &[Self::Nuc] {
        self.as_codons().as_flattened()
    }

    /// Cast to mutable slice of nucleotides
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{DnaSlice, Nuc};
    /// use Nuc::{A, C, G, T};
    ///
    /// let mut codons = [
    ///     [C, A, T],
    ///     [T, A, G],
    ///     [A, C, T],
    /// ];
    /// codons.as_flat_dna_mut()[3] = Nuc::G;
    /// assert_eq!(
    ///     codons,
    ///     [
    ///         [C, A, T],
    ///         [G, A, G],
    ///         [A, C, T],
    ///     ]
    /// );
    /// ```
    fn as_flat_dna_mut(&mut self) -> &mut [Self::Nuc] {
        self.as_codons_mut().as_flattened_mut()
    }

    /// Cast to slice of codons, discarding leftover trailing nucleotides.
    ///
    /// This returns the first reading frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{DnaSlice, Nuc};
    /// use Nuc::{A, C, T};
    ///
    /// let dna = Nuc::lit(b"CATATTAC");
    /// assert_eq!(
    ///     dna.as_codons(),
    ///     [[C, A, T], [A, T, T]]
    /// );
    /// ```
    fn as_codons(&self) -> &[[Self::Nuc; 3]] {
        self.as_flat_dna().as_chunks().0
    }

    /// Cast to mutable slice of codons, discarding leftover trailing nucleotides.
    ///
    /// This returns the first reading frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{DnaSlice, Nuc};
    ///
    /// let mut dna = Nuc::lit(b"CATATTAC");
    /// let codons = dna.as_codons_mut();
    /// // Set the second codon's first nucleotide...
    /// codons[1][0] = Nuc::G;
    /// assert_eq!(dna, Nuc::lit(b"CATGTTAC"));
    /// ```
    fn as_codons_mut(&mut self) -> &mut [[Self::Nuc; 3]] {
        self.as_flat_dna_mut().as_chunks_mut().0
    }

    /// Return all 3 reading frames of codons
    ///
    /// Akin to non-panicking version of:
    /// `[dna[0..].as_codons(), dna[1..].as_codons(), dna[2..].as_codons()]`
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{DnaSlice, Nuc};
    /// use Nuc::{A, C, T};
    ///
    /// let dna = Nuc::lit(b"ACATATTAC");
    /// assert_eq!(
    ///     dna.reading_frames(),
    ///     [
    ///         &[[A, C, A], [T, A, T], [T, A, C]] as &[_],
    ///         &[[C, A, T], [A, T, T]],
    ///         &[[A, T, A], [T, T, A]],
    ///     ]
    /// );
    /// ```
    fn reading_frames(&self) -> [&[[Self::Nuc; 3]]; 3] {
        std::array::from_fn(|i| self.as_flat_dna().get(i..).unwrap_or_default().as_codons())
    }

    /// Return iterator translating codons into amino acids.
    ///
    /// The given [`GeneticCode`] is applied to the first reading frame's codons
    /// (discarding leftover trailing nucleotides).
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{DnaSlice, NCBI1, Nuc, Seq};
    ///
    /// let peptide: Seq<Vec<_>> = Nuc::lit(b"TATGCGAGAAAC").translate(NCBI1).collect();
    /// assert_eq!(peptide.to_string(), "YARN");
    /// ```
    fn translate<G: GeneticCode>(
        &self,
        genetic_code: G,
    ) -> Translated<G, Codons<Self::Nuc, std::slice::Iter<'_, Self::Nuc>>> {
        self.as_flat_dna().iter().translate(genetic_code)
    }

    /// Translate codons into peptide [`Vec`].
    ///
    /// For large sequences, this is usually much faster than populating directly from an iterator.
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{DnaSlice, NCBI1, Nuc, Seq};
    ///
    /// let peptide = Seq(Nuc::lit(b"TATGCGAGAAAC").translate_to_vec(NCBI1));
    /// assert_eq!(peptide.to_string(), "YARN");
    /// ```
    fn translate_to_vec<G: GeneticCode>(
        &self,
        genetic_code: G,
    ) -> Vec<<Self::Nuc as Nucleotide>::Amino> {
        let codons = self.as_codons();
        let mut peptide = vec![Default::default(); codons.len()];
        codons.translate_to_buf(genetic_code, &mut peptide);
        peptide
    }

    /// Translate codons into fixed-length peptide.
    ///
    /// # Panics
    ///
    /// Panics if the number of codons to be translated is different from the returned array.
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{Amino, DnaSlice, NCBI1, Nuc, Seq};
    ///
    /// let peptide: Seq<[Amino; 4]> = Seq(Nuc::lit(b"TATGCGAGAAAC").translate_to_array(NCBI1));
    /// assert_eq!(peptide.to_string(), "YARN");
    /// ```
    fn translate_to_array<G: GeneticCode, const N: usize>(
        &self,
        genetic_code: G,
    ) -> [<Self::Nuc as Nucleotide>::Amino; N] {
        let mut buf = [Default::default(); _];
        self.translate_to_buf(genetic_code, &mut buf);
        buf
    }

    /// Fill a buffer with amino acids built from translated codons.
    ///
    /// For large sequences, this is usually much faster than populating directly from an iterator.
    ///
    /// # Panics
    ///
    /// Panics if the number of codons to be translated is different from the length of `buf`.
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{Amino, DnaSlice, NCBI1, Nuc, Seq};
    ///
    /// let mut peptide: Seq<[Amino; 4]> = Default::default();
    /// Nuc::lit(b"TATGCGAGAAAC").translate_to_buf(NCBI1, &mut *peptide);
    /// assert_eq!(peptide.to_string(), "YARN");
    /// ```
    fn translate_to_buf<G: GeneticCode>(
        &self,
        genetic_code: G,
        buf: &mut [<Self::Nuc as Nucleotide>::Amino],
    ) {
        const CHUNK_LEN: usize = 16;
        let codons = self.as_codons();
        assert_eq!(codons.len(), buf.len());
        let (codon_chunks, codon_remainder) = codons.as_chunks::<CHUNK_LEN>();
        let (amino_chunks, amino_remainder) = buf.as_chunks_mut::<CHUNK_LEN>();
        for (aminos, codons) in std::iter::zip(amino_chunks, codon_chunks) {
            for (amino, codon) in std::iter::zip(aminos, codons) {
                *amino = Self::Nuc::translate(&genetic_code, *codon);
            }
        }
        for (amino, codon) in std::iter::zip(amino_remainder, codon_remainder) {
            *amino = Self::Nuc::translate(&genetic_code, *codon);
        }
    }

    /// Return object that implements [`Display`](std::fmt::Display)
    /// for printing sequence compactly. See [`nucs::iter::Display`](crate::iter::Display)
    /// for more details.
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{DnaSlice, Nuc};
    ///
    /// let dna = Nuc::lit(b"CATATTAC");
    /// assert_eq!(dna.display().to_string(), "CATATTAC");
    /// assert_eq!(format!("{:#4}", dna.display()), "CATA\nTTAC");
    /// ```
    fn display(&self) -> crate::iter::Display<std::slice::Iter<'_, Self::Nuc>> {
        self.as_flat_dna().iter().display()
    }

    /// Cast to slice of [`AmbiNuc`].
    ///
    /// <div class="warning">
    ///
    /// This requires the `unsafe` feature.
    ///
    /// </div>
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{AmbiNuc, DnaSlice, Nuc};
    ///
    /// let dna = Nuc::lit(b"CATATTAC");
    /// assert_eq!(dna.as_ambi_nucs(), AmbiNuc::lit(b"CATATTAC"));
    /// ```
    #[cfg(feature = "unsafe")]
    #[cfg_attr(docsrs, doc(cfg(feature = "unsafe")))]
    fn as_ambi_nucs(&self) -> &[AmbiNuc] {
        Self::Nuc::as_ambi_nucs(self.as_flat_dna())
    }

    /// Attempt to cast to slice of [`Nuc`].
    ///
    /// [`None`] is returned if any nucleotides are degenerate (inexpressible by [`Nuc`]).
    ///
    /// <div class="warning">
    ///
    /// This requires the `unsafe` feature.
    ///
    /// </div>
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{AmbiNuc, DnaSlice, Nuc};
    ///
    /// let dna = AmbiNuc::lit(b"CATATTAC");
    /// assert_eq!(dna.to_nucs().unwrap(), Nuc::lit(b"CATATTAC"));
    ///
    /// let dna = AmbiNuc::lit(b"CATTY");
    /// assert!(dna.to_nucs().is_none());
    /// ```
    #[cfg(feature = "unsafe")]
    #[cfg_attr(docsrs, doc(cfg(feature = "unsafe")))]
    fn to_nucs(&self) -> Option<&[Nuc]> {
        Self::Nuc::to_nucs(self.as_flat_dna())
    }

    /// Attempt to cast to mutable slice of [`Nuc`].
    ///
    /// [`None`] is returned if any nucleotides are degenerate (inexpressible by [`Nuc`]).
    ///
    /// <div class="warning">
    ///
    /// This requires the `unsafe` feature.
    ///
    /// </div>
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{AmbiNuc, DnaSlice, Nuc};
    ///
    /// let mut dna = AmbiNuc::lit(b"CATATTAC");
    /// if let Some(nucs) = dna.to_nucs_mut() {
    ///     nucs[7] = Nuc::G;
    /// }
    /// assert_eq!(dna, AmbiNuc::lit(b"CATATTAG"));
    ///
    /// let mut dna = AmbiNuc::lit(b"CATTY");
    /// assert!(dna.to_nucs_mut().is_none());
    /// ```
    #[cfg(feature = "unsafe")]
    #[cfg_attr(docsrs, doc(cfg(feature = "unsafe")))]
    fn to_nucs_mut(&mut self) -> Option<&mut [Nuc]> {
        Self::Nuc::to_nucs_mut(self.as_flat_dna_mut())
    }

    /// Perform in-place complement.
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{DnaSlice, Nuc};
    ///
    /// let mut dna = Nuc::lit(b"CATATTAC");
    /// dna.complement();
    /// assert_eq!(dna, Nuc::lit(b"GTATAATG"));
    /// ```
    fn complement(&mut self) {
        self.as_flat_dna_mut().iter_mut().complement();
    }

    /// Perform in-place reverse-complement.
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{DnaSlice, Nuc};
    ///
    /// let mut dna = Nuc::lit(b"CATATTAC");
    /// dna.revcomp();
    /// assert_eq!(dna, Nuc::lit(b"GTAATATG"));
    /// ```
    fn revcomp(&mut self) {
        self.as_flat_dna_mut().iter_mut().revcomp();
    }
}

impl<N: Nucleotide> DnaSlice for [N] {
    type Nuc = N;

    fn as_flat_dna(&self) -> &[N] {
        self
    }

    fn as_flat_dna_mut(&mut self) -> &mut [N] {
        self
    }
}

impl<N: Nucleotide> DnaSlice for [[N; 3]] {
    type Nuc = N;

    fn as_codons(&self) -> &[[N; 3]] {
        self
    }

    fn as_codons_mut(&mut self) -> &mut [[N; 3]] {
        self
    }
}

#[cfg(test)]
mod tests {
    use crate::{NCBI1, Nuc, Seq};

    use super::*;

    // My previous attempt at the iter appeared to succeed, but later turned out not to be
    // callable in certain generic situations. In order to guard against that, I'm using these
    // functions to test everything with minimal available type info.

    fn anon_slice(dna: &[impl Nucleotide]) -> &[impl Nucleotide] {
        dna
    }

    fn anon_mut_slice(dna: &mut [impl Nucleotide]) -> &mut [impl Nucleotide] {
        dna
    }

    #[test]
    fn translated_type_inference() {
        let mut dna = Nuc::lit(b"AAAACCCGGT");
        let peptide: Seq<Vec<_>> = anon_slice(&dna).translate(NCBI1).collect();
        assert_eq!(peptide.to_string(), "KTR");
        let peptide: Seq<Vec<_>> = anon_mut_slice(&mut dna).translate(NCBI1).collect();
        assert_eq!(peptide.to_string(), "KTR");
    }

    #[test]
    fn display_type_inference() {
        let mut dna = Nuc::lit(b"AAAACCCGGT");
        assert_eq!(anon_slice(&dna).display().to_string(), "AAAACCCGGT");
        assert_eq!(anon_mut_slice(&mut dna).display().to_string(), "AAAACCCGGT");
    }
}
