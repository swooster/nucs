#![warn(missing_docs)]
#![warn(clippy::pedantic)]
#![allow(clippy::inline_always)]
#![cfg_attr(not(feature = "unsafe"), forbid(unsafe_code))]
#![cfg_attr(docsrs, feature(doc_cfg))]

//! [`nucs`](crate) is a library for working with nucleotide and amino acid sequences.
//!
//! The goal is to supply useful tools for working with DNA/peptides while attempting to
//! integrate with the rest of Rust.
//!
//! ```
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! // `Nuc` represents concrete nucleotides, and `Dna` holds `Nuc`s
//! use nucs::{Dna, Nuc};
//!
//! // `Dna` can be parsed, modified and displayed.
//! let mut dna: Dna = "CATG".parse()?;
//! dna.extend([Nuc::A, Nuc::G]);
//! assert_eq!(dna.to_string(), "CATGAG");
//!
//! // For convenience, there's a helper to build const literals:
//! const CAT: &[Nuc] = &Nuc::lit(b"CAT");
//! assert!(dna.starts_with(CAT));
//!
//! // `Seq` is a wrapper to add convenience features to `Vec`-like collections
//! use nucs::Seq;
//! // and `Dna` is actually just an alias for `Seq<Vec<Nuc>>`
//! let dna: Seq<Vec<Nuc>> = dna;
//!
//! // ...but `Seq` can wrap any sufficiently `Vec`-like collection:
//! let mut dna = Seq(std::collections::VecDeque::from_iter(dna));
//! dna[3] = Nuc::T;
//! dna.push_front(Nuc::A);
//! // Displayed `Seq`s can be line-wrapped by using alternate formatting:
//! assert_eq!(format!("{dna:#4}"), "ACAT\nTAG");
//!
//! // `DnaSlice` supplies helpers for working with slices:
//! use nucs::DnaSlice;
//! use Nuc::{A, C, G, T};
//! let slice = dna.make_contiguous();
//! assert_eq!(
//!     slice.reading_frames(),
//!     [
//!         &[[A, C, A], [T, T, A]],
//!         &[[C, A, T], [T, A ,G]],
//!         &[[A, T, T]],
//!     ] as [&[_]; 3]
//! );
//! slice.revcomp(); // in-place reverse-complement
//! assert_eq!(dna.to_string(), "CTAATGT");
//!
//! // `DnaIter` supplies helpers for working with DNA iterators non-destructively:
//! use nucs::DnaIter;
//!
//! let iter = dna
//!     .iter()
//!     .trimmed_to_codon()
//!     .revcomped();
//! // (cloneable) DNA iterators can be displayed too:
//! let wrapped = format!("{:#3}", iter.display());
//! assert_eq!(wrapped, "CAT\nTAG");
//!
//! // Ambiguous nucleotides represent non-empty sets of nucleotides.
//! use nucs::AmbiNuc;
//!
//! // `Nuc`s can be composed into `AmbiNuc`s...
//! assert_eq!(C | A | T, AmbiNuc::H);
//! // ...which can be decomposed back into `Nuc`s
//! let dna = AmbiNuc::lit(b"STRAYGYMNAST");
//! assert!(dna[0].iter().eq([C, G]));
//! assert!(dna[1].iter().eq([T]));
//! assert!(dna[8].iter().eq(Nuc::ALL));
//!
//! // Both concrete and ambiguous amino acids are supported as well:
//! use nucs::{Amino, AmbiAmino};
//!
//! let peptide = Seq(Amino::lit(b"KITTY*PAWS"));
//! assert_eq!(format!("{peptide:#5}"), "KITTY\n*PAWS");
//!
//! assert_eq!(Amino::I | Amino::L, AmbiAmino::J);
//! assert!((Amino::C | Amino::A | Amino::T).iter().eq(Amino::lit(b"ACT")));
//!
//! // And it's easy to translate DNA into peptides:
//! use nucs::NCBI1;
//!
//! let dna = Nuc::lit(b"TTTGAGCTCATAAACGAGA");
//! let peptide: Seq<Vec<_>> = dna.translate(NCBI1).collect();
//! assert_eq!(peptide.to_string(), "FELINE");
//!
//! // Even ambiguous DNA can be translated:
//! let dna = AmbiNuc::lit(b"MTTGCGTCTCCCGAGCGC");
//! let peptide: Seq<Vec<_>> = dna.translate(NCBI1).collect();
//! assert_eq!(peptide.to_string(), "JASPER");
//! # Ok(())
//! # }
//! ```
//!
//! # Features
//!
//! * **`proptest`:** Enables [`proptest`](https://crates.io/crates/proptest) integration
//!   and utils, particularly `Arbitrary` generation of [`Nuc`], [`AmbiNuc`], [`Amino`] and
//!   [`AmbiAmino`].
//! * **`rand`:** Enables [`rand`](https://crates.io/crates/rand) integration, particularly
//!   `StandardUniform` generation of [`Nuc`], [`AmbiNuc`], [`Amino`] and [`AmbiAmino`].
//! * **`serde`:** Enables [`serde`](https://crates.io/crates/serde) integration for [`Seq<T>`].
//! * **`unsafe`:** (experimental) This enables casting between [`&[Nuc]`](crate::Nuc)
//!   and [`&[AmbiNuc]`](crate::AmbiNuc).

#[doc = include_str!("../README.md")]
#[cfg(doctest)]
pub struct ReadmeDoctests;

mod amino;
mod nuc;
mod seq;
mod symbol;

#[cfg(feature = "unsafe")]
#[cfg_attr(docsrs, doc(cfg(feature = "unsafe")))]
pub mod casts;
pub mod error;
pub mod iter;
#[cfg(any(feature = "proptest", test))]
#[cfg_attr(docsrs, doc(cfg(feature = "proptest")))]
pub mod proptest;
#[cfg(feature = "rand")]
mod rand;
pub mod slice;
pub mod translation;

pub use amino::{AmbiAmino, Amino};
pub use iter::DnaIter;
pub use nuc::{AmbiNuc, Nuc, Nucleotide};
pub use seq::Seq;
pub use slice::DnaSlice;
pub use symbol::Symbol;
pub use translation::{NCBI1, NCBI1_RC};

/// Common nucleotide sequence type
pub type Dna = Seq<Vec<Nuc>>;
/// Common ambiguous nucleotide sequence type
pub type AmbiDna = Seq<Vec<AmbiNuc>>;
/// Common amino acid sequence type
pub type Peptide = Seq<Vec<Amino>>;
/// Common ambiguous amino acid sequence type
pub type AmbiPeptide = Seq<Vec<AmbiAmino>>;
