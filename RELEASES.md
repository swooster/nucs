# Releases

## Unreleased

* Add `proptest` integration for `Nuc`, `AmbiNuc`, `Amino` and `AmbiAmino`,
  as well as DNA and peptides.
* Add `rand` integration for `Nuc`, `AmbiNuc`, `Amino` and `AmbiAmino`.
* Expand `DnaIter::revcomp` to work with iters of any `impl AsMut<Nuc>` rather than
  just `&mut Nuc`.

## Version 0.1.2 (2025-07-12)

* Add `serde` integration for `Seq<T>`.

## Version 0.1.1 (2025-07-10)

* Fix documentation on docs.rs.

## Version 0.1.0 (2025-07-10)

* First release!
