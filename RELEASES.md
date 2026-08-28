# Releases

## Unreleased

* Packed types now support comparison to slices, arrays and `Vec`s, strings
  and each other. (#79, #80, #81)
* Packed types now support mutable iteration. (#82)
* Packed vec-like types now support `push` and `pop`. (#83)

## Version 0.4.1 (2026-08-13)

* `Seq` now has `complement` and `reverse_complement` methods. (#61)
* `Seq` now supports `AsRef` and `AsMut`. (#67)
* Add basic support for packed types! (#62, #66, #69, #70, #73, #74, #75, #76)

## Version 0.4.0 (2026-07-28)

* `Seq` now wraps subslices so e.g. `assert_eq!(&dna[2..], "ATCG");` now works. (#46)
* `Seq` now has limited support for `ToOwned`, `Borrow` and `BorrowMut`. (#47)
* `DnaSlice` and `DnaIter` were renamed to `DnaSliceExt` and `DnaIterExt` to reflect that
  they're extension traits. (#48, #49)
* Rename `revcomp{,ed}` to `reverse_complement{,ed}`. The latter is annoying to type,
  but easier to remember. (#51)
* Rework translation to use more of a builder pattern. (#52, #53, #54)
* Allow comparison between `Symbol`s and strings. (#57)
* Make `DnaSliceExt::reading_frames` return flat DNA and add `Seq::reading_frames`. (#58)
* Upgrade optional `rand` dependency to 0.10.0. (#59)

## Version 0.3.1 (2026-04-20)

* Tweak [`docs.rs`](https://docs.rs/nucs/latest/nucs/) settings to build with all features
  again. (#44)

## Version 0.3.0 (2026-04-10)

* Loads of adjustments to translation. There are now methods on `DnaSlice` and `Seq` to perform
  reverse complement translation almost as fast as ordinary translation, and a number of
  translation-related APIs have changed. (#23, #24, #25, #26, #27, #30, #31, #38, #41)
* Usability improvements to `Seq`:
  * You can now do e.g. `Nuc::seq(b"ATCG")` as a shorthand for `Seq(Nuc::lit(b"ATCG"))`. (#29)
  * `Seq` can now wrap slices, not just owned data. (#32)
  * `Seq` can be compared to strings for easier testing. (#28)
  * `Seq` now has translation methods. (#30)
* Rename `T::lit` to `T::arr` to improve clarity. (#42)

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
