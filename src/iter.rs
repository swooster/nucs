//! Iterator-related types

use std::fmt::Formatter;
use std::marker::PhantomData;

use crate::translation::GeneticCode;
use crate::{AmbiAmino, Amino, Nucleotide};

/// Helpers for working with iterators of [`Nucleotide`]s.
pub trait DnaIter: Iterator {
    /// Return iterator over reverse-complemented nucleotides.
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{DnaIter, Nuc};
    ///
    /// let complement = Nuc::lit(b"GATTACA").into_iter().revcomped();
    /// assert!(complement.eq(Nuc::lit(b"TGTAATC")));
    /// ```
    fn revcomped<N>(self) -> Complemented<N, std::iter::Rev<Self>>
    where
        Self: Sized + DoubleEndedIterator<Item: AsRef<N>>,
        N: Nucleotide,
    {
        self.rev().complemented()
    }

    /// Return iterator over complemented nucleotides.
    ///
    /// This is like mapping the iterator through [`Nucleotide::complement`].
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{DnaIter, Nuc};
    ///
    /// let complement = Nuc::lit(b"GATTACA").into_iter().complemented();
    /// assert!(complement.eq(Nuc::lit(b"CTAATGT")));
    /// ```
    fn complemented<N>(self) -> Complemented<N, Self>
    where
        Self: Sized + Iterator<Item: AsRef<N>>,
        N: Nucleotide,
    {
        Complemented {
            nuc_type: PhantomData,
            iter: self,
        }
    }

    /// Perform in-place reverse-complement, consuming the iterator.
    ///
    /// This swaps the front and back of the iterator and sets them to their
    /// [`Nucleotide::complement`].
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{DnaIter, Nuc};
    ///
    /// let mut dna = Nuc::lit(b"GATTACA");
    /// dna.iter_mut().revcomp();
    /// assert_eq!(dna, Nuc::lit(b"TGTAATC"));
    /// ```
    fn revcomp<'a, N>(mut self)
    where
        Self: Sized + DoubleEndedIterator<Item = &'a mut N>,
        N: Nucleotide,
    {
        while let Some(a) = self.next() {
            *a = a.complement();
            let Some(b) = self.next_back() else { break };
            *a = std::mem::replace(b, *a).complement();
        }
    }

    /// Perform in-place complement, consuming the iterator.
    ///
    /// This is sets each element to its [`Nucleotide::complement`].
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{DnaIter, Nuc};
    ///
    /// let mut dna = Nuc::lit(b"GATTACA");
    /// dna.iter_mut().complement();
    /// assert_eq!(dna, Nuc::lit(b"CTAATGT"));
    /// ```
    fn complement<'a, N>(self)
    where
        Self: Sized + Iterator<Item = &'a mut N>,
        N: Nucleotide,
    {
        for nuc in self {
            *nuc = nuc.complement();
        }
    }

    /// Return iterator over codons of first reading frame.
    ///
    /// This discards any leftover trailing nucleotides.
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{DnaIter, Nuc};
    /// use Nuc::{A, C, G, T};
    ///
    /// let codons = Nuc::lit(b"GATTACA").into_iter().codons();
    /// assert!(codons.eq([
    ///     [G, A, T],
    ///     [T, A, C],
    /// ]));
    /// ```
    fn codons<N>(self) -> Codons<N, Self>
    where
        Self: Sized + Iterator<Item: AsRef<N>>,
        N: Nucleotide,
    {
        Codons {
            nuc_type: PhantomData,
            iter: self,
        }
    }

    /// Discard trailing nucleotides that aren't part of the first reading frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{DnaIter, Nuc};
    /// use Nuc::{A, C, G, T};
    ///
    /// let codons = Nuc::lit(b"GATTACA").into_iter().trimmed_to_codon();
    /// assert!(codons.eq([
    ///     G, A, T,
    ///     T, A, C,
    /// ]));
    /// ```
    #[must_use]
    fn trimmed_to_codon<N>(mut self) -> Self
    where
        Self: Sized + DoubleEndedIterator<Item: AsRef<N>> + ExactSizeIterator,
        N: Nucleotide,
    {
        for _ in 0..(self.len() % 3) {
            self.next_back();
        }
        self
    }

    /// Return iterator translating codons into amino acids.
    ///
    /// The given [`GeneticCode`] is applied to the first reading frame's codons
    /// (discarding leftover trailing nucleotides).
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{DnaIter, NCBI1, Nuc, Seq};
    ///
    /// let peptide: Seq<Vec<_>> =
    ///     Nuc::lit(b"TATGCGAGAAAC").into_iter().translate(NCBI1).collect();
    /// assert_eq!(peptide.to_string(), "YARN");
    /// ```
    fn translate<N, G>(self, genetic_code: G) -> Translated<G, Codons<N, Self>>
    where
        Self: Sized + Iterator<Item: AsRef<N>>,
        N: Nucleotide,
        G: GeneticCode,
    {
        Translated {
            genetic_code,
            iter: self.codons(),
        }
    }

    /// Return object that implements [`Display`] for printing the iterated sequence compactly.
    ///
    /// # Examples
    ///
    /// ```
    /// use nucs::{DnaIter, Nuc};
    ///
    /// let dna = Nuc::lit(b"GATTACA").into_iter().display();
    /// assert_eq!(format!("{dna:#4}"), "GATT\nACA");
    /// ```
    fn display<N>(&self) -> Display<Self>
    where
        Self: Clone + Iterator<Item: AsRef<N>>,
        N: Nucleotide,
    {
        Display(self.clone())
    }
}

impl<I: Iterator> DnaIter for I {}

/// An iterator that complements the nucleotides of an underlying iterator.
///
/// This is created by the [`DnaIter::complemented`] method. See its documentation for details.
#[must_use = "iterators are lazy and do nothing unless consumed"]
#[derive(Clone, Debug, Default)]
pub struct Complemented<N, I> {
    nuc_type: PhantomData<N>,
    iter: I,
}

impl<N, I> Iterator for Complemented<N, I>
where
    N: Nucleotide,
    I: Iterator<Item: AsRef<N>>,
{
    type Item = N;

    fn next(&mut self) -> Option<N> {
        self.iter.next().map(|n| n.as_ref().complement())
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.iter.size_hint()
    }
}

impl<N, I> DoubleEndedIterator for Complemented<N, I>
where
    N: Nucleotide,
    I: DoubleEndedIterator<Item: AsRef<N>>,
{
    fn next_back(&mut self) -> Option<N> {
        self.iter.next_back().map(|n| n.as_ref().complement())
    }
}

impl<N, I> ExactSizeIterator for Complemented<N, I>
where
    N: Nucleotide,
    I: ExactSizeIterator<Item: AsRef<N>>,
{
    fn len(&self) -> usize {
        self.iter.len()
    }
}

/// An iterator that groups nucleotides of an underlying iterator into codons.
///
/// This is created by the [`DnaIter::codons`] method. See its documentation for details.
#[must_use = "iterators are lazy and do nothing unless consumed"]
#[derive(Clone, Debug, Default)]
pub struct Codons<N, I> {
    nuc_type: PhantomData<N>,
    iter: I,
}

impl<N, I> Iterator for Codons<N, I>
where
    N: Nucleotide,
    I: Iterator<Item: AsRef<N>>,
{
    type Item = [N; 3];

    fn next(&mut self) -> Option<Self::Item> {
        let n1 = self.iter.next()?;
        let n2 = self.iter.next()?;
        let n3 = self.iter.next()?;
        Some([n1, n2, n3].map(|n| *n.as_ref()))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let (min, max) = self.iter.size_hint();
        (min / 3, max.map(|m| m / 3))
    }
}

impl<N, I> DoubleEndedIterator for Codons<N, I>
where
    N: Nucleotide,
    I: DoubleEndedIterator<Item: AsRef<N>> + ExactSizeIterator,
{
    fn next_back(&mut self) -> Option<Self::Item> {
        let _ = (&mut self.iter).trimmed_to_codon();
        let n3 = self.iter.next_back()?;
        let n2 = self.iter.next_back()?;
        let n1 = self.iter.next_back()?;
        Some([n1, n2, n3].map(|n| *n.as_ref()))
    }
}

impl<N, I> ExactSizeIterator for Codons<N, I>
where
    N: Nucleotide,
    I: ExactSizeIterator<Item: AsRef<N>>,
{
    fn len(&self) -> usize {
        self.iter.len() / 3
    }
}

/// An iterator that translates codons via a [`GeneticCode`].
///
/// This can be created by the [`DnaIter::translate`] method. See its documentation for details.
#[must_use = "iterators are lazy and do nothing unless consumed"]
#[derive(Clone, Debug, Default)]
pub struct Translated<G, I> {
    genetic_code: G,
    iter: I,
}

impl<I, N, G> Iterator for Translated<G, I>
where
    I: Iterator<Item = [N; 3]>,
    N: Nucleotide,
    G: GeneticCode + Clone,
{
    type Item = N::Amino;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter
            .next()
            .map(|codon| N::translate(codon, self.genetic_code.clone()))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.iter.size_hint()
    }
}

impl<I, N, G> DoubleEndedIterator for Translated<G, I>
where
    I: DoubleEndedIterator<Item = [N; 3]>,
    N: Nucleotide,
    G: GeneticCode + Clone,
{
    fn next_back(&mut self) -> Option<Self::Item> {
        self.iter
            .next_back()
            .map(|codon| N::translate(codon, self.genetic_code.clone()))
    }
}

impl<I, N, G> ExactSizeIterator for Translated<G, I>
where
    I: ExactSizeIterator<Item = [N; 3]>,
    N: Nucleotide,
    G: GeneticCode + Clone,
{
    fn len(&self) -> usize {
        self.iter.len()
    }
}

/// A wrapper that provides [`Display`](std::fmt::Display) and [`Debug`](std::fmt::Debug)
/// impls contatenating the items of its contained cloneable iterables.
///
/// # Examples
///
/// ```
/// use nucs::{Amino, DnaIter, Nuc};
///
/// let dna = Nuc::lit(b"GATTACA").into_iter().display();
///
/// // The default formatting just prints everything on a single line:
/// assert_eq!(dna.to_string(), "GATTACA");
///
/// // Alternate formatting enables line-wrap with a default width of 80 characters,
/// // but said width can be customized:
/// assert_eq!(format!("{dna:#4}"), "GATT\nACA");
/// assert_eq!(format!("{dna:#3}"), "GAT\nTAC\nA");
///
/// // The `Debug` impl is identical to the `Display` impl:
/// assert_eq!(format!("{dna:#4?}"), "GATT\nACA");
///
/// // `Amino`s are supported too, though there's not currently a helper trait for that.
/// let peptide = Amino::lit(b"PEPTIDE");
/// let peptide = nucs::iter::Display::new(peptide);
/// assert_eq!(format!("{peptide:#4}"), "PEPT\nIDE");
/// ```
#[must_use]
#[derive(Clone)]
pub struct Display<I>(I);

impl<I> Display<I> {
    /// Construct a new [`Display`] from a cloneable iterable.
    ///
    /// This is most often accessed via the [`DnaIter::display`] method.
    pub fn new(iterable: I) -> Self
    where
        // make it easy to notice problems earlier
        I: IntoIterator<Item: std::fmt::Display> + Clone,
    {
        Self(iterable)
    }
}

impl<I> std::fmt::Display for Display<I>
where
    I: IntoIterator<Item: std::fmt::Display> + Clone,
{
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        if f.alternate() {
            // Ok, realistically, lots of DNA may be a bit too large to fit on one line...
            // for convenience, let's make it possible to request wrapping via the alternate
            // format, using the width the control wrap length. I'm not 100% sure this is
            // the best way to go. Maybe I should instead make a wrapper type?
            // We'll assume that each element is one character.
            let width = f.width().unwrap_or(80).max(1);
            let mut at_start = true;
            for (i, symbol) in (0..width).cycle().zip(self.0.clone()) {
                if i == 0 && !at_start {
                    f.write_str("\n")?;
                }
                at_start = false;
                write!(f, "{symbol}")?;
            }
        } else {
            for symbol in self.0.clone() {
                write!(f, "{symbol}")?;
            }
        }
        Ok(())
    }
}

impl<I: IntoIterator<Item: std::fmt::Display> + Clone> std::fmt::Debug for Display<I> {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        std::fmt::Display::fmt(self, f)
    }
}

/// Iterator of possible [`Amino`]s an [`AmbiAmino`] could be.
///
/// This is created by [`AmbiAmino::iter`].
#[must_use]
pub struct AmbiAminoIter(u32);

impl AmbiAminoIter {
    pub(crate) fn new(amino: AmbiAmino) -> Self {
        let bits = amino.to_bits().get();
        // Clear bit 0 and move the Amino::Stop bit to it. We do this for two reasons:
        // * so we don't have to worry about the dummy bit
        // * so iteration order matches `Amino` order
        let stop_mask = AmbiAmino::bit_mask(Amino::Stop);
        let shifted_stop = (bits >> AmbiAmino::bit_offset(Amino::Stop)) & 1;
        Self(bits & !1 & !stop_mask | shifted_stop)
    }

    fn amino_from_bit_offset(bit: u32) -> Amino {
        let table = const {
            let mut table = [Amino::Stop; 26];
            let mut i = 0;
            while i < Amino::ALL.len() {
                let amino = Amino::ALL[i];
                table[AmbiAmino::bit_offset(amino) as usize] = amino;
                i += 1;
            }
            table
        };
        table[bit as usize]
    }
}

impl Iterator for AmbiAminoIter {
    type Item = Amino;

    fn next(&mut self) -> Option<Self::Item> {
        let next_bit = self.0.trailing_zeros();
        if next_bit >= 32 {
            return None;
        }
        self.0 &= !(1 << next_bit);
        Some(Self::amino_from_bit_offset(next_bit))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = self.0.count_ones() as usize;
        (len, Some(len))
    }
}

impl DoubleEndedIterator for AmbiAminoIter {
    fn next_back(&mut self) -> Option<Self::Item> {
        let next_bit = self.0.leading_zeros();
        if next_bit >= 32 {
            return None;
        }
        let next_bit = 31 - next_bit;
        self.0 &= !(1 << next_bit);
        Some(Self::amino_from_bit_offset(next_bit))
    }
}

impl ExactSizeIterator for AmbiAminoIter {}

#[cfg(test)]
mod tests {
    use crate::test_util::all_dna;
    use crate::{Dna, DnaSlice, NCBI1, Nuc, Seq};

    use super::*;

    // My previous attempt at this iter appeared to succeed, but later turned out not to be
    // callable in certain generic situations. In order to guard against that, I'm using these
    // functions to test everything with minimal available type info.

    fn anon_vals(iter: impl IntoIterator<Item: Nucleotide>) -> impl IntoIterator<Item: Nucleotide> {
        iter
    }

    fn anon_refs<'a>(
        iter: impl IntoIterator<Item = &'a impl Nucleotide>,
    ) -> impl IntoIterator<Item = &'a impl Nucleotide> {
        iter
    }

    fn anon_muts<'a>(
        iter: impl IntoIterator<Item = &'a mut impl Nucleotide>,
    ) -> impl IntoIterator<Item = &'a mut impl Nucleotide> {
        iter
    }

    #[test]
    fn complemented_type_inference() {
        let mut dna = Nuc::lit(b"AAAACCCGGT");
        let vals: Seq<Vec<_>> = anon_vals(dna).into_iter().complemented().collect();
        assert_eq!(vals.to_string(), "TTTTGGGCCA");
        let refs: Seq<Vec<_>> = anon_refs(&dna).into_iter().complemented().collect();
        assert_eq!(refs.to_string(), "TTTTGGGCCA");
        let muts: Seq<Vec<_>> = anon_muts(&mut dna).into_iter().complemented().collect();
        assert_eq!(muts.to_string(), "TTTTGGGCCA");
    }

    #[test]
    fn complement_type_inference() {
        let mut dna = Nuc::lit(b"AAAACCCGGT");
        anon_muts(&mut dna).into_iter().complement();
        assert_eq!(dna, Nuc::lit(b"TTTTGGGCCA"));
    }

    #[test]
    fn revcomp_type_inference() {
        fn anon_muts<'a>(
            iter: impl IntoIterator<IntoIter: DoubleEndedIterator, Item = &'a mut impl Nucleotide>,
        ) -> impl IntoIterator<IntoIter: DoubleEndedIterator, Item = &'a mut impl Nucleotide>
        {
            iter
        }

        let mut dna = Nuc::lit(b"AAAACCCGGT");
        anon_muts(&mut dna).into_iter().revcomp();
        assert_eq!(dna, Nuc::lit(b"ACCGGGTTTT"));
    }

    #[test]
    fn codons_type_inference() {
        let mut dna = Nuc::lit(b"AAAACCCGGT");
        let vals: Vec<_> = anon_vals(dna)
            .into_iter()
            .codons()
            .map(|c| Seq(c).to_string())
            .collect();
        assert_eq!(vals, ["AAA", "ACC", "CGG"]);
        let refs: Vec<_> = anon_refs(&dna)
            .into_iter()
            .codons()
            .map(|c| Seq(c).to_string())
            .collect();
        assert_eq!(refs, ["AAA", "ACC", "CGG"]);
        let muts: Vec<_> = anon_muts(&mut dna)
            .into_iter()
            .codons()
            .map(|c| Seq(c).to_string())
            .collect();
        assert_eq!(muts, ["AAA", "ACC", "CGG"]);
    }

    #[test]
    fn trimmed_to_codon_type_inference() {
        fn anon_vals(
            iter: impl IntoIterator<IntoIter: DoubleEndedIterator + ExactSizeIterator, Item: Nucleotide>,
        ) -> impl IntoIterator<IntoIter: DoubleEndedIterator + ExactSizeIterator, Item: Nucleotide>
        {
            iter
        }

        fn anon_refs<'a>(
            iter: impl IntoIterator<
                IntoIter: DoubleEndedIterator + ExactSizeIterator,
                Item = &'a impl Nucleotide,
            >,
        ) -> impl IntoIterator<
            IntoIter: DoubleEndedIterator + ExactSizeIterator,
            Item = &'a impl Nucleotide,
        > {
            iter
        }

        fn anon_muts<'a>(
            iter: impl IntoIterator<
                IntoIter: DoubleEndedIterator + ExactSizeIterator,
                Item = &'a mut impl Nucleotide,
            >,
        ) -> impl IntoIterator<
            IntoIter: DoubleEndedIterator + ExactSizeIterator,
            Item = &'a mut impl Nucleotide,
        > {
            iter
        }

        let mut dna = Nuc::lit(b"AAAACCCGGT");
        let vals: Seq<Vec<_>> = anon_vals(dna).into_iter().trimmed_to_codon().collect();
        assert_eq!(vals.to_string(), "AAAACCCGG");
        let refs: Seq<Vec<_>> = anon_refs(&dna).into_iter().trimmed_to_codon().collect();
        assert_eq!(refs.to_string(), "AAAACCCGG");
        let muts: Seq<Vec<_>> = anon_muts(&mut dna).into_iter().trimmed_to_codon().collect();
        assert_eq!(muts.to_string(), "AAAACCCGG");
    }

    #[test]
    fn translate_type_inference() {
        let mut dna = Nuc::lit(b"AAAACCCGGT");
        let vals: Seq<Vec<_>> = anon_vals(dna).into_iter().translate(NCBI1).collect();
        assert_eq!(vals.to_string(), "KTR");
        let refs: Seq<Vec<_>> = anon_refs(&dna).into_iter().translate(NCBI1).collect();
        assert_eq!(refs.to_string(), "KTR");
        let muts: Seq<Vec<_>> = anon_muts(&mut dna).into_iter().translate(NCBI1).collect();
        assert_eq!(muts.to_string(), "KTR");
    }

    #[test]
    fn display_type_inference() {
        fn anon_vals(
            iter: impl IntoIterator<IntoIter: Clone, Item: Nucleotide>,
        ) -> impl IntoIterator<IntoIter: Clone, Item: Nucleotide> {
            iter
        }

        fn anon_refs<'a>(
            iter: impl IntoIterator<IntoIter: Clone, Item = &'a impl Nucleotide>,
        ) -> impl IntoIterator<IntoIter: Clone, Item = &'a impl Nucleotide> {
            iter
        }

        let dna = Nuc::lit(b"AAAACCCGGT");
        let vals = anon_vals(dna).into_iter().display();
        assert_eq!(vals.to_string(), "AAAACCCGGT");
        let refs = anon_refs(&dna).into_iter().display();
        assert_eq!(refs.to_string(), "AAAACCCGGT");
    }

    #[test]
    fn revcomp() {
        let mut dna = Nuc::lit(b"");
        dna.iter_mut().revcomp(); // just sanity check a lack of panics

        let mut dna = Nuc::lit(b"A");
        dna.iter_mut().revcomp();
        assert_eq!(dna, Nuc::lit(b"T"));

        let mut dna = Nuc::lit(b"AC");
        dna.iter_mut().revcomp();
        assert_eq!(dna, Nuc::lit(b"GT"));

        let mut dna = Nuc::lit(b"ACT");
        dna.iter_mut().revcomp();
        assert_eq!(dna, Nuc::lit(b"AGT"));

        let mut dna = Nuc::lit(b"GACT");
        dna.iter_mut().revcomp();
        assert_eq!(dna, Nuc::lit(b"AGTC"));
    }

    #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
    #[test]
    fn revcomp_is_involution() {
        for dna in all_dna::<Nuc>().take(10000) {
            let mut dna2: Vec<_> = dna.iter().revcomped().revcomped().collect();
            assert_eq!(dna, dna2);
            dna2.iter_mut().revcomp();
            dna2.iter_mut().revcomp();
            assert_eq!(dna, dna2);
        }
    }

    #[cfg_attr(miri, ignore = "slow in miri; shouldn't touch unsafe code anyway")]
    #[test]
    fn revcomp_matches_reverse_and_complement() {
        // Iterator revcomp isn't trivial, so compare against something that is
        for mut dna in all_dna::<Nuc>().take(10000) {
            let mut expected = dna.clone();
            expected.reverse();
            expected.complement();
            dna.iter_mut().revcomp();
            assert_eq!(dna, expected);
        }
    }

    #[test]
    fn display_formats() {
        let disp = |s: &str| s.parse::<Dna>().unwrap().into_iter().display();
        assert_eq!(disp("").to_string(), "");
        assert_eq!(format!("{:5}", disp("")), "");
        assert_eq!(disp("C").to_string(), "C");
        assert_eq!(format!("{:5}", disp("C")), "C");
        let display = disp("AAAACCCGGTACGTACGT");
        assert_eq!(display.to_string(), "AAAACCCGGTACGTACGT");
        assert_eq!(format!("{display:?}"), "AAAACCCGGTACGTACGT");
        // Don't allow widths of 0. We treat it as 1 instead of panicking or showing nothing.
        let width = 0;
        assert_eq!(
            format!("{display:#width$}"),
            "A\nA\nA\nA\nC\nC\nC\nG\nG\nT\nA\nC\nG\nT\nA\nC\nG\nT"
        );
        assert_eq!(
            format!("{display:#1}"),
            "A\nA\nA\nA\nC\nC\nC\nG\nG\nT\nA\nC\nG\nT\nA\nC\nG\nT"
        );
        assert_eq!(
            format!("{display:#1}"),
            "A\nA\nA\nA\nC\nC\nC\nG\nG\nT\nA\nC\nG\nT\nA\nC\nG\nT"
        );
        assert_eq!(
            format!("{display:#2}"),
            "AA\nAA\nCC\nCG\nGT\nAC\nGT\nAC\nGT"
        );
        assert_eq!(format!("{display:#3}"), "AAA\nACC\nCGG\nTAC\nGTA\nCGT");
        assert_eq!(format!("{display:#4}"), "AAAA\nCCCG\nGTAC\nGTAC\nGT");
        // Default width is 80
        let display = disp(
            "
            AAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTT
            AAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTTACGT
        ",
        );
        assert_eq!(
            format!("{display:#}"),
            "\
            AAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTT\
            AAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTT\nACGT\
        "
        );
    }
}
