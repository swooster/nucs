//! Unsafe casts
//!
//! <div class="warning">
//!
//! This module (and everything in it) requires the `unsafe` feature.
//!
//! </div>

use crate::{AmbiNuc, Nuc};

/// Asserts that `[Nuc]` and `[AmbiNuc]` have interchangeable memory representations if:
/// * Both have the same length.
/// * The contained nucleotides consist only of A, C, G or T.
/// * Both `Nuc` and `AmbiNuc` are `repr(u8)` (they should be).
const fn assert_nuc_and_concrete_ambi_nuc_buffers_have_same_representation() {
    assert!(size_of::<Nuc>() == size_of::<AmbiNuc>());
    assert!(align_of::<Nuc>() == align_of::<AmbiNuc>());
    assert!(size_of::<Nuc>() == 1);
    // Demonstrate that Nuc has no values other than A/C/G/T:
    let _ = |(Nuc::A | Nuc::C | Nuc::G | Nuc::T)| {};
    assert!(Nuc::A as u8 == AmbiNuc::A as u8);
    assert!(Nuc::C as u8 == AmbiNuc::C as u8);
    assert!(Nuc::G as u8 == AmbiNuc::G as u8);
    assert!(Nuc::T as u8 == AmbiNuc::T as u8);
}

/// Checks if it would be safe to cast an `[AmbiNuc]` to `[Nuc]`
fn is_representable_as_nucs(nucs: &[AmbiNuc]) -> bool {
    nucs.iter().all(|nuc| Nuc::try_from(*nuc).is_ok())
}

/// Cast [`&[Nuc]`](slice) to [`&[AmbiNuc]`](slice).
///
/// You may wish to use [`DnaSlice::as_ambi_nucs`](crate::DnaSlice::as_ambi_nucs) instead.
///
/// # Examples
///
/// ```
/// use nucs::{AmbiNuc, Nuc};
///
/// let dna = Nuc::lit(b"CATATTAC");
/// assert_eq!(nucs::casts::nucs_as_ambi(&dna), AmbiNuc::lit(b"CATATTAC"));
/// ```
#[cfg_attr(docsrs, doc(cfg(feature = "unsafe")))]
#[must_use]
pub fn nucs_as_ambi(nucs: &[Nuc]) -> &[AmbiNuc] {
    const { assert_nuc_and_concrete_ambi_nuc_buffers_have_same_representation() };
    // It's valid to reinterpret this pointer due to the above assertion.
    let data = nucs.as_ptr().cast::<AmbiNuc>();
    let len = nucs.len();
    // SAFETY:
    // * `data` points to a valid contiguous buffer of `AmbiNuc` (see above)
    // * `data`'s lifetime is sufficient because it was obtained from a slice with
    //   the same lifetime as the one we're returning
    unsafe { std::slice::from_raw_parts(data, len) }
}

/// Attempt to cast [`&[AmbiNuc]`](slice) to [`&[Nuc]`](slice).
///
/// [`None`] is returned if any nucleotides are degenerate (inexpressible by [`Nuc`]).
///
/// You may wish to use [`DnaSlice::to_nucs`](crate::DnaSlice::to_nucs) instead.
///
/// # Examples
///
/// ```
/// use nucs::{AmbiNuc, Nuc};
///
/// let dna = AmbiNuc::lit(b"CATATTAC");
/// assert_eq!(nucs::casts::ambi_to_nucs(&dna).unwrap(), Nuc::lit(b"CATATTAC"));
///
/// let dna = AmbiNuc::lit(b"CATTY");
/// assert!(nucs::casts::ambi_to_nucs(&dna).is_none());
/// ```
#[cfg_attr(docsrs, doc(cfg(feature = "unsafe")))]
#[must_use]
pub fn ambi_to_nucs(nucs: &[AmbiNuc]) -> Option<&[Nuc]> {
    if is_representable_as_nucs(nucs) {
        const { assert_nuc_and_concrete_ambi_nuc_buffers_have_same_representation() };
        // It's valid to reinterpret this pointer due to the above check and assertion.
        let data = nucs.as_ptr().cast::<Nuc>();
        let len = nucs.len();
        // SAFETY:
        // * `data` points to a valid contiguous buffer of `Nuc` (see above)
        //   In particular, we've checked that every element of `nucs` is also a valid value
        //   of `Nuc` and asserted `Nuc` and `Nucs` have the same memory representations.
        // * `data`'s lifetime is sufficient because it was obtained from a slice with
        //    the same lifetime as the one we're returning
        Some(unsafe { std::slice::from_raw_parts(data, len) })
    } else {
        None
    }
}

/// Attempt to cast [`&mut [AmbiNuc]`](slice) to [`&mut [Nuc]`](slice).
///
/// [`None`] is returned if any nucleotides are degenerate (inexpressible by [`Nuc`]).
///
/// You may wish to use [`DnaSlice::to_nucs_mut`](crate::DnaSlice::to_nucs_mut) instead.
///
/// # Examples
///
/// ```
/// use nucs::{AmbiNuc, DnaSlice, Nuc};
///
/// let mut dna = AmbiNuc::lit(b"CATATTAC");
/// if let Some(nucs) = nucs::casts::ambi_to_nucs_mut(&mut dna) {
///     nucs[7] = Nuc::G;
/// }
/// assert_eq!(dna, AmbiNuc::lit(b"CATATTAG"));
///
/// let mut dna = AmbiNuc::lit(b"CATTY");
/// assert!(nucs::casts::ambi_to_nucs_mut(&mut dna).is_none());
/// ```
#[cfg_attr(docsrs, doc(cfg(feature = "unsafe")))]
#[must_use]
pub fn ambi_to_nucs_mut(nucs: &mut [AmbiNuc]) -> Option<&mut [Nuc]> {
    if is_representable_as_nucs(nucs) {
        const { assert_nuc_and_concrete_ambi_nuc_buffers_have_same_representation() };
        // It's valid to reinterpret this pointer due to the above check and assertion.
        let data = nucs.as_mut_ptr().cast::<Nuc>();
        let len = nucs.len();
        // SAFETY:
        // * `data` points to a valid contiguous buffer of `Nuc` (see above)
        //   In particular, we've checked that every element of `nucs` is also a valid value
        //   of `Nuc` and asserted `Nuc` and `Nucs` have the same memory representations.
        // * `data`'s lifetime is sufficient because it was obtained from a slice with
        //    the same lifetime as the one we're returning
        // * `data` will not be accessed through any other pointer while this slice exists
        //   because we obtained it from an exclusive reference.
        // * No matter how the callar mutates the returned slice, valid values of
        //   `Nuc` are also valid values of `AmbiNuc`, so the mutations should be valid changes
        //   to the memory representation of the original slice.
        Some(unsafe { std::slice::from_raw_parts_mut(data, len) })
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use std::ops::Range;

    use proptest::prelude::{Just, ProptestConfig, Strategy, any};
    use proptest::{prop_compose, prop_oneof, proptest};

    use crate::proptest::{any_ambi_dna, any_dna};

    use super::*;

    #[test]
    fn nucs_as_ambi_smoke_test() {
        assert_eq!(nucs_as_ambi(&[] as &[Nuc]), &[] as &[AmbiNuc]);

        assert_eq!(
            nucs_as_ambi(&Nuc::lit(b"GATTACA")),
            AmbiNuc::lit(b"GATTACA")
        );
    }

    #[test]
    fn ambi_to_nucs_smoke_test() {
        assert_eq!(ambi_to_nucs(&[] as &[AmbiNuc]).unwrap(), &[] as &[Nuc]);

        assert_eq!(
            ambi_to_nucs(&AmbiNuc::lit(b"CATTAG")).unwrap(),
            Nuc::lit(b"CATTAG")
        );

        assert_eq!(ambi_to_nucs(&AmbiNuc::lit(b"CATSTAG")), None);
    }

    #[test]
    fn ambi_to_nucs_mut_smoke_test() {
        assert_eq!(
            ambi_to_nucs_mut(&mut [] as &mut [AmbiNuc]).unwrap(),
            &[] as &[Nuc]
        );

        let mut buffer = AmbiNuc::lit(b"CATTAG");
        let nuc_buf = ambi_to_nucs_mut(&mut buffer).unwrap();
        nuc_buf[0] = Nuc::A;
        nuc_buf[1] = Nuc::C;
        nuc_buf[5] = Nuc::T;
        assert_eq!(nuc_buf, Nuc::lit(b"ACTTAT"));
        assert_eq!(buffer, AmbiNuc::lit(b"ACTTAT"));

        assert_eq!(ambi_to_nucs_mut(&mut AmbiNuc::lit(b"CATSTAG")), None);
    }

    // Like any_ambi_dna(...) except convertible to `&[Nuc]` more than half the time.
    fn often_convertible_dna(size: Range<usize>) -> impl Strategy<Value = Vec<AmbiNuc>> {
        prop_oneof![
            any_ambi_dna(size.clone()),
            any_dna(size).prop_map(|dna| dna.into_iter().map(Into::into).collect()),
        ]
    }

    prop_compose! {
        fn dna_and_index(size: Range<usize>)
                        (dna in any_dna(size))
                        (index in 0..dna.len(), dna in Just(dna))
                        -> (Vec<AmbiNuc>, usize) {
            (dna.into_iter().map(Into::into).collect(), index)
        }
    }

    proptest! {
        #![proptest_config(ProptestConfig {
            #[cfg(miri)]
            failure_persistence: None,
            .. Default::default()
        })]

        #[test]
        fn nucs_as_ambi_matches_conversion(dna in any_dna(0..50)) {
            let expected: Vec<_> = dna.iter().map(|&n| AmbiNuc::from(n)).collect();
            let ambi_dna = nucs_as_ambi(&dna);
            assert_eq!(ambi_dna, expected);
        }

        #[test]
        fn ambi_to_nucs_matches_conversion(ambi_dna in often_convertible_dna(0..50)) {
            let expected: Option<Vec<_>> =
                ambi_dna.iter().map(|&n| Nuc::try_from(n).ok()).collect();
            let dna = ambi_to_nucs(&ambi_dna);
            assert_eq!(dna, expected.as_deref());
        }

        #[test]
        fn ambi_to_nucs_mut_matches_conversion(mut ambi_dna in often_convertible_dna(0..50)) {
            let mut expected: Option<Vec<_>> =
                ambi_dna.iter().map(|&n| Nuc::try_from(n).ok()).collect();
            let dna = ambi_to_nucs_mut(&mut ambi_dna);
            assert_eq!(dna, expected.as_deref_mut());
        }

        #[test]
        fn ambi_to_nucs_mut_provides_mut_ref(
            (mut ambi_dna, index) in dna_and_index(1..50),
            nuc in any::<Nuc>(),
        ) {
            ambi_to_nucs_mut(&mut ambi_dna).unwrap()[index] = nuc;
            assert_eq!(ambi_dna[index], AmbiNuc::from(nuc));
        }
    }
}
