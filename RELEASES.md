# Releases

## Unreleased

* Tweak `Nucleotide::translate`, swapping order of arguments and taking reference. This eliminates
  a lot of `Clone` bounds in APIs. (#23)
* Add `GeneticCode::translate` to provide a simpler alternative to `Nucleotide::translate`.
  Give longer names to `GeneticCode::map_codon` and `GeneticCode::map_ambi_codon`, to better
  reflect that they're not intended to be used directly. (#24)
* Add `Clone` and `Debug` impls for `FastTranslator`. (#25)
* Add `FastTranslator::from_genetic_code` and tweak `FastTranslator::from_table` to take
  table reference instead of value. (#26)
* Add `Seq::translated_*` methods. (#30)
* Allow `Seq` to be compared to strings. (#28)
* Add `Symbol::seq` helper. (#29)

## Version 0.2.0 (2025-11-24)

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
