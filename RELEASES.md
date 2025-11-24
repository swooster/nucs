# Releases

## Unreleased

* Add `proptest` integration for `Nuc`, `AmbiNuc`, `Amino` and `AmbiAmino`,
  as well as DNA and peptides. (#11)
* Add `rand` integration for `Nuc`, `AmbiNuc`, `Amino` and `AmbiAmino`. (#14)
* Expand `DnaIter::revcomp` to work with iters of any `impl AsMut<Nuc>` rather than
  just `&mut Nuc`. (#12)
* `Symbol` now requires `Default`. Additionally `Seq` impls `Default`. (#20)
* Add `FastTranslator`, `DnaSlice::translate_to_vec`, `DnaSlice::translate_to_array`
  and `DnaSlice::translate_to_buf` to dramatically improve translation speed. (#18)

## Version 0.1.2 (2025-07-12)

* Add `serde` integration for `Seq<T>`. (#5)

## Version 0.1.1 (2025-07-10)

* Fix documentation on docs.rs.

## Version 0.1.0 (2025-07-10)

* First release!
