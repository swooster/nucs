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
* ...integrates with the Rust `std` library.
* ...is collection-agnostic where feasible.
* ...tries to be consistent with Zipf's law of abbrevation for type names.

The high-level API tries to declutter basic usage:
```rust
use nucs::{Dna, NCBI1};

let dna: Dna = "TGTGCCACCAATATTCCC".parse().unwrap();
// For testing convience, string comparisons are supported:
assert_eq!(dna, "tgtgccaccaatattccc");
assert_eq!(dna, "TGT GCC ACC AAT ATT CCC");
assert_eq!(dna[5..15], "CACCAATATT");
assert_eq!(dna.translated_by(NCBI1), "CATNIP");
```

...but tools are provided to aid using other collections. For example,
the following code uses arrays to completely avoid allocations:
```rust
use nucs::{Amino, DnaSliceExt, NCBI1, Nuc};

const SNEK: [Nuc; 23] = Nuc::arr(b"CCCGGGAATGGTTGGTCCTAGAG");

let mut dna = SNEK;

// Select this:         v---------------------v
//                CCCG  GGA ATG GTT GGT CCT AGA  G
let codons = dna[4..].as_codons_mut();

// Reverse complement this: v---------v
//                CCCG  GGA ATG GTT GGT CCT AGA  G
// Changing it to:          ACC AAC CAT
codons[1..4].reverse_complement();

// (performs no allocations)
let mut peptide: [Amino; 7] =
    dna[2..].translated_by(NCBI1).try_into().unwrap();

assert_eq!(peptide, Amino::arr(b"REPTILE"));
peptide[0] = Amino::P;
peptide[5] = Amino::D;
assert_eq!(peptide, Amino::arr(b"PEPTIDE"));
```

Ambiguous nucleotides and amino acids are also supported:
```rust
use nucs::{AmbiDna, AmbiNuc, NCBI1};

let mut dna: AmbiDna = "TTAGCGGACGATTAT".parse().unwrap();

// Because `dna` contains ambiguous nucleotides,
// translating it produces an ambiguous peptide
assert_eq!(dna.translated_by(NCBI1), "LADDY");

use AmbiNuc::{A, C};
dna[0] |= A | C;
dna[6] |= A;
dna[9] |= A;

assert_eq!(dna.translated_by(NCBI1), "JABBY");
```

## Planned functionality

* Packing
* FASTA parsing
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

