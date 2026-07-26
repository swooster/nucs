//! Slice-related types

use crate::translation::{GeneticCode, Translation};
use crate::{DnaIterExt, Nucleotide, Seq};

#[cfg(feature = "unsafe")]
use crate::symbol::sealed::Sealed;
#[cfg(feature = "unsafe")]
use crate::{AmbiNuc, Nuc};

/// Helpers for working with slices of [`Nucleotide`]s.
pub trait DnaSliceExt {
    /// The type of [`Nucleotide`] this slice is made of.
    type Nuc: Nucleotide;

    /// Cast to slice of nucleotides
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{DnaSliceExt, Nuc, Seq};
    /// use Nuc::{A, C, G, T};
    ///
    /// let codons = [
    ///     [C, A, T],
    ///     [T, A, G],
    ///     [A, C, T],
    /// ];
    /// assert_eq!(codons.as_flat_dna(), Nuc::arr(b"CATTAGACT"));
    /// ```
    fn as_flat_dna(&self) -> &[Self::Nuc] {
        self.as_codons().as_flattened()
    }

    /// Cast to mutable slice of nucleotides
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{DnaSliceExt, Nuc};
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

    /// Cast to slice of codons, discarding excess trailing nucleotides.
    ///
    /// This returns the first reading frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{DnaSliceExt, Nuc};
    /// use Nuc::{A, C, T};
    ///
    /// let dna = Nuc::arr(b"CATATTAC");
    /// assert_eq!(
    ///     dna.as_codons(),
    ///     [[C, A, T], [A, T, T]]
    /// );
    /// ```
    fn as_codons(&self) -> &[[Self::Nuc; 3]] {
        self.as_flat_dna().as_chunks().0
    }

    /// Cast to mutable slice of codons, discarding excess trailing nucleotides.
    ///
    /// This returns the first reading frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{DnaSliceExt, Nuc};
    ///
    /// let mut dna = Nuc::arr(b"CATATTAC");
    /// let codons = dna.as_codons_mut();
    /// // Set the second codon's first nucleotide...
    /// codons[1][0] = Nuc::G;
    /// assert_eq!(dna, Nuc::arr(b"CATGTTAC"));
    /// ```
    fn as_codons_mut(&mut self) -> &mut [[Self::Nuc; 3]] {
        self.as_flat_dna_mut().as_chunks_mut().0
    }

    /// Cast to slice of codons, discarding excess leading nucleotides.
    ///
    /// This returns the first reading frame from the end.
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{DnaSliceExt, Nuc};
    /// use Nuc::{A, C, T};
    ///
    /// let dna = Nuc::arr(b"CATATTAC");
    /// assert_eq!(
    ///     dna.as_rcodons(),
    ///     [[T, A, T], [T, A, C]]
    /// );
    /// ```
    fn as_rcodons(&self) -> &[[Self::Nuc; 3]] {
        self.as_flat_dna().as_rchunks().1
    }

    /// Cast to mutable slice of codons, discarding excess leading nucleotides.
    ///
    /// This returns the first reading frame from the end.
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{DnaSliceExt, Nuc};
    ///
    /// let mut dna = Nuc::arr(b"CATATTAC");
    /// let codons = dna.as_rcodons_mut();
    /// // Set the second codon's first nucleotide...
    /// codons[1][0] = Nuc::G;
    /// assert_eq!(dna, Nuc::arr(b"CATATGAC"));
    /// ```
    fn as_rcodons_mut(&mut self) -> &mut [[Self::Nuc; 3]] {
        self.as_flat_dna_mut().as_rchunks_mut().1
    }

    /// Wrap slice in [`Seq`].
    ///
    /// Alias for [`Seq::wrap`].
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{DnaSliceExt, Nuc};
    ///
    /// let dna = Nuc::seq(b"ACTGACTG");
    /// let partial_dna = dna[3..6].as_seq();
    /// assert_eq!(partial_dna, "GAC");
    /// ```
    fn as_seq(&self) -> &Seq<Self> {
        Seq::wrap(self)
    }

    /// Wrap mutable slice in [`Seq`].
    ///
    /// Alias for [`Seq::wrap`].
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{DnaSliceExt, Nuc};
    ///
    /// let mut dna = Nuc::seq(b"ACTGACTG");
    /// let partial_dna = dna[3..6].as_seq_mut();
    /// assert_eq!(partial_dna, "GAC");
    /// partial_dna[1] = Nuc::C;
    /// assert_eq!(dna, "ACTGCCTG");
    /// ```
    fn as_seq_mut(&mut self) -> &mut Seq<Self> {
        Seq::wrap_mut(self)
    }

    /// Return all 3 reading frames of codons
    ///
    /// Akin to non-panicking version of:
    /// `[dna[0..].as_codons(), dna[1..].as_codons(), dna[2..].as_codons()]`
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{DnaSliceExt, Nuc};
    /// use Nuc::{A, C, T};
    ///
    /// let dna = Nuc::arr(b"ACATATTAC");
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

    /// Return translation of DNA via [`GeneticCode`].
    ///
    /// This returns a builder for performing translations. See [`Translation`] for details.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::collections::VecDeque;
    ///
    /// use nucs::{AmbiNuc, DnaSliceExt, NCBI1, Nuc, Seq};
    ///
    /// let dna = Nuc::arr(b"TATGCGAGAAAC");
    /// assert_eq!(dna.translated_by(NCBI1).to_seq(), "YARN");
    ///
    /// let rc_dna = AmbiNuc::arr(b"NGCACCGCTAGGTACTGGCGAA");
    /// let peptide = rc_dna.translated_by(NCBI1).reverse_complemented().to_seq();
    /// assert_eq!(peptide, "FAST*RC");
    ///
    /// let peptide: Seq<VecDeque<_>> = dna.translated_by(NCBI1).into_iter().collect();
    /// assert_eq!(peptide, "YARN");
    /// ```
    fn translated_by<G: GeneticCode>(&self, genetic_code: G) -> Translation<'_, Self::Nuc, G> {
        Translation {
            dna: self.as_flat_dna(),
            genetic_code,
        }
    }

    /// Return object that implements [`Display`](std::fmt::Display)
    /// for printing sequence compactly. See [`nucs::iter::Display`](crate::iter::Display)
    /// for more details.
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{DnaSliceExt, Nuc};
    ///
    /// let dna = Nuc::arr(b"CATATTAC");
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
    /// use nucs::{AmbiNuc, DnaSliceExt, Nuc};
    ///
    /// let dna = Nuc::arr(b"CATATTAC");
    /// assert_eq!(dna.as_ambi_nucs(), AmbiNuc::arr(b"CATATTAC"));
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
    /// use nucs::{AmbiNuc, DnaSliceExt, Nuc};
    ///
    /// let dna = AmbiNuc::arr(b"CATATTAC");
    /// assert_eq!(dna.to_nucs().unwrap(), Nuc::arr(b"CATATTAC"));
    ///
    /// let dna = AmbiNuc::arr(b"CATTY");
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
    /// use nucs::{AmbiNuc, DnaSliceExt, Nuc};
    ///
    /// let mut dna = AmbiNuc::arr(b"CATATTAC");
    /// if let Some(nucs) = dna.to_nucs_mut() {
    ///     nucs[7] = Nuc::G;
    /// }
    /// assert_eq!(dna, AmbiNuc::arr(b"CATATTAG"));
    ///
    /// let mut dna = AmbiNuc::arr(b"CATTY");
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
    /// use nucs::{DnaSliceExt, Nuc};
    ///
    /// let mut dna = Nuc::arr(b"CATATTAC");
    /// dna.complement();
    /// assert_eq!(dna, Nuc::arr(b"GTATAATG"));
    /// ```
    fn complement(&mut self) {
        self.as_flat_dna_mut().iter_mut().complement();
    }

    /// Perform in-place reverse-complement.
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{DnaSliceExt, Nuc};
    ///
    /// let mut dna = Nuc::arr(b"CATATTAC");
    /// dna.reverse_complement();
    /// assert_eq!(dna, Nuc::arr(b"GTAATATG"));
    /// ```
    fn reverse_complement(&mut self) {
        self.as_flat_dna_mut().iter_mut().reverse_complement();
    }
}

impl<N: Nucleotide> DnaSliceExt for [N] {
    type Nuc = N;

    fn as_flat_dna(&self) -> &[N] {
        self
    }

    fn as_flat_dna_mut(&mut self) -> &mut [N] {
        self
    }
}

impl<N: Nucleotide> DnaSliceExt for [[N; 3]] {
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
        let mut dna = Nuc::arr(b"AAAACCCGGT");
        let peptide: Seq<Vec<_>> = anon_slice(&dna).translated_by(NCBI1).into_iter().collect();
        assert_eq!(peptide, "KTR");
        let peptide: Seq<Vec<_>> = anon_mut_slice(&mut dna)
            .translated_by(NCBI1)
            .into_iter()
            .collect();
        assert_eq!(peptide, "KTR");
    }

    #[test]
    fn display_type_inference() {
        let mut dna = Nuc::arr(b"AAAACCCGGT");
        assert_eq!(anon_slice(&dna).display().to_string(), "AAAACCCGGT");
        assert_eq!(anon_mut_slice(&mut dna).display().to_string(), "AAAACCCGGT");
    }
}
