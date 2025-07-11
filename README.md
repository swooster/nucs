# `nucs`

[![Crates.io Version](https://img.shields.io/crates/v/nucs.svg)](https://crates.io/crates/nucs)
[![CI](https://github.com/swooster/nucs/actions/workflows/ci.yml/badge.svg?event=push)](https://github.com/swooster/nucs/actions/workflows/ci.yml)
[![Documentation](https://docs.rs/nucs/badge.svg)](https://docs.rs/nucs)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE-MIT)
[![License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE-APACHE)

This is a personal experiment with an API for working with nucleotide and amino acid sequences.
Its design is based off of my experience using and helping maintain
<https://github.com/SecureDNA/quickdna>. My goals were to make an API that...

* ...is solely focused on Rust.
* ...integrates with the Rust `std` library (e.g. representing codons as `[Nuc; 3]` allows the
  `std` library to understand that codons can be cheaply flattened into nucleotides).
* ...is (largely) collection-agnostic.
* ...tries to be consistent with Zipf's law of abbreviation when naming things.

```rust
use nucs::{Dna, DnaSlice, NCBI1, Nuc, Peptide};

let mut dna: Dna = "ACACACATATCTTACGCTTAGGAAATCTGACCCGAA"
    .parse().unwrap();

let codons = dna[4..].as_codons_mut();
// Selects this: v-------------------------------------v
//       ACAC    ACA TAT CTT ACG CTT AGG AAA TCT GAC CCG    AA

codons[3..8].revcomp();
// Reverse complements this: v-----------------v
//       ACAC    ACA TAT CTT ACG CTT AGG AAA TCT GAC CCG    AA
// Changing it to:           AGA TTT CCT AAG CGT

dna.extend(const { Nuc::lit(b"CCAACCATTGATGAG") });

let peptide: Peptide = dna.translate(NCBI1).collect();
assert_eq!(peptide.to_string(), "THIS*IS*A*PEPTIDE");
```

Non-`Vec` containers are supported too, and it's possible to work with DNA non-destructively
via iterators:
```rust
use std::collections::VecDeque;
use nucs::{DnaIter, NCBI1, Nuc, Peptide, Seq};

let mut dna: Seq<VecDeque<Nuc>> =
    "ACTCTATCACCTACTCAGAGCGCTCCACCGCGCGTGT".parse().unwrap();
// Prepend things to the `VecDeque`;
// it's no longer stored contiguously.
for _ in 0..4 {
    dna.push_front(Nuc::C);
}

let immutable_dna = dna;
// Apply reverse compliment and NCBI1 non-destructively.
let peptide: Peptide = immutable_dna
    .iter()
    .revcomped()
    .translate(NCBI1)
    .collect();
assert_eq!(peptide.to_string(), "TRAVERSE*VIEW");
```

Ambiguous nucleotides and amino acids are supported:
```rust
use nucs::{AmbiAmino, AmbiNuc, AmbiPeptide, DnaSlice, NCBI1};
use AmbiNuc::{A, C};

let mut dna = AmbiNuc::lit(b"TTAGCGGACGATTAT");

// Because `dna` contains ambiguous nucleotides,
// translating it produces an ambiguous peptide
let peptide: AmbiPeptide = dna.translate(NCBI1).collect();
assert_eq!(peptide.to_string(), "LADDY");

dna[0] |= A | C;
dna[6] |= A;
dna[9] |= A;
assert_eq!(dna.display().to_string(), "HTAGCGRACRATTAT");
let peptide: AmbiPeptide = dna.translate(NCBI1).collect();
assert_eq!(peptide.to_string(), "JABBY");
```

## Planned functionality

* Packing
* FASTA parsing
* `serde` integration
* Expansion of ambiguous k-mers into concrete k-mers
* Base canonicalization
* Unsafe casts for `Vec` and `Arc`
* Better efficiency

## Incompatibility with `quickdna`

Note that while `nucs` is heavily inspired by <https://github.com/SecureDNA/quickdna>,
there are subtle-yet-important incompatibilities with the order and representation of
nucleotides and amino acids. In particular nucleotides are ordered alphabetically
in `nucs`, to keep the ordering identical to strings as well as (hopefully) making
future bit-packing work easier.

## License

Licensed under either of

 * Apache License, Version 2.0
   ([LICENSE-APACHE](LICENSE-APACHE) or <http://www.apache.org/licenses/LICENSE-2.0>)
 * MIT license
   ([LICENSE-MIT](LICENSE-MIT) or <http://opensource.org/licenses/MIT>)

at your option.

## Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.

