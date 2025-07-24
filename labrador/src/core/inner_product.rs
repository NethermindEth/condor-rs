use std::{
    borrow::Borrow,
    ops::{Add, Mul},
};

/// Computes the **linear combination** of two equally-sized slices, returning the
/// sum of pair-wise products.
///
/// Conceptually, given two slices  
/// `elements = [A1, A2, … , An]` and `challenges = [c1, c2, … , cn]`,  
/// the function evaluates  
///
/// ```text
/// Σ (Ai · ci)  for i = 1 .. n
/// ```
///
/// where `·` is multiplication for the underlying types and `+` is the
/// corresponding addition. The result’s type (`B`) is the same as the borrowed
/// type of each element in `elements`.
///
/// # Arguments
///
/// * `elements`   – Slice of items whose borrowed values will be multiplied by
///   the corresponding entries in `challenges`.
/// * `challenges` – Slice of multipliers—**must** be the same length as
///   `elements`.
///
/// # Returns
///
/// The accumulated sum **Σ `(elements[i] * challenges[i])`** as a value of type
/// `B`.
///
///
/// # Examples
///
/// Basic usage with primitive integers:
///
/// ```rust
/// use std::borrow::Borrow;
/// use rand::rng;
/// use labrador::ring::rq::Rq;
/// use labrador::core::inner_product::compute_linear_combination;
///
/// // The function is generic; we can call it directly.
/// let elems = [Rq::random(&mut rng()), Rq::random(&mut rng()), Rq::random(&mut rng()), Rq::random(&mut rng()), Rq::random(&mut rng())];
/// let challenges = [Rq::random(&mut rng()), Rq::random(&mut rng()), Rq::random(&mut rng()), Rq::random(&mut rng()), Rq::random(&mut rng())];
/// let sum = compute_linear_combination(&elems, &challenges);
/// ```
///
pub fn compute_linear_combination<E, B, C>(elements: &[E], challenges: &[C]) -> B
where
    E: Borrow<B>,
    for<'a> &'a B: Mul<&'a C, Output = B>,
    for<'a> &'a B: Add<&'a B, Output = B>,
{
    debug_assert_eq!(
        elements.len(),
        challenges.len(),
        "vectors must be the same length"
    );
    debug_assert!(!elements.is_empty(), "`elements` must not be empty");

    let mut zipped_iter = elements.iter().zip(challenges.iter());
    // Must do the following as the init value in fold requires size of B
    let (e0, c0) = zipped_iter.next().unwrap();
    let init = e0.borrow() * c0;

    zipped_iter.fold(init, |acc, (elem, c)| &acc + &(elem.borrow() * c))
}

#[cfg(test)]
mod tests {
    use rand::rng;

    use crate::ring::{rq::Rq, rq_matrix::RqMatrix, rq_vector::RqVector, zq::Zq};

    use super::*;

    #[test]
    #[should_panic]
    fn test_inputs_with_different_lengths_panic() {
        let vector_a = RqVector::random(&mut rng(), 100);
        let vector_b = RqVector::random(&mut rng(), 90);
        let _ = compute_linear_combination(vector_a.elements(), vector_b.elements());
    }

    #[test]
    fn test_zq_type_inputs() {
        let poly_a = Rq::random(&mut rng());
        let poly_b = Rq::random(&mut rng());

        let result = compute_linear_combination(poly_a.coeffs(), poly_b.coeffs());

        let mut expected = Zq::ZERO;
        for i in 0..poly_a.coeffs().len() {
            expected += poly_a.coeffs()[i] * poly_b.coeffs()[i];
        }
        assert_eq!(result, expected);
    }

    #[test]
    fn test_rq_type_inputs() {
        let poly_vec_a = RqVector::random(&mut rng(), 100);
        let poly_vec_b = RqVector::random(&mut rng(), 100);

        let result = compute_linear_combination(poly_vec_a.elements(), poly_vec_b.elements());

        let mut expected = Rq::zero();
        for i in 0..poly_vec_a.elements().len() {
            expected = &expected + &(&poly_vec_a.elements()[i] * &poly_vec_b.elements()[i]);
        }
        assert_eq!(result, expected);
    }

    #[test]
    fn test_zq_rq_type_inputs() {
        let poly_a = Rq::random(&mut rng());
        let poly_vec_b = RqVector::random(&mut rng(), 64);

        let result = compute_linear_combination(poly_vec_b.elements(), poly_a.coeffs());

        let mut expected = Rq::zero();
        for i in 0..poly_a.coeffs().len() {
            expected = &expected + &(&poly_vec_b.elements()[i] * &poly_a.coeffs()[i]);
        }
        assert_eq!(result, expected);
    }

    #[test]
    fn test_zq_rqvector_type_inputs() {
        let poly_a = Rq::random(&mut rng());
        let poly_vec_b = RqMatrix::random(&mut rng(), 64, 100);

        let result = compute_linear_combination(poly_vec_b.elements(), poly_a.coeffs());

        let mut expected = RqVector::zero(100);
        for i in 0..poly_a.coeffs().len() {
            expected = &expected + &(&poly_vec_b.elements()[i] * poly_a.coeffs()[i]);
        }
        assert_eq!(result, expected);
    }

    #[test]
    fn test_rq_rqvector_type_inputs() {
        let poly_vec_a = RqVector::random(&mut rng(), 80);
        let poly_vec_b = RqMatrix::random(&mut rng(), 80, 100);

        let result = compute_linear_combination(poly_vec_b.elements(), poly_vec_a.elements());

        let mut expected = RqVector::zero(100);
        for i in 0..poly_vec_a.elements().len() {
            expected = &expected + &(&poly_vec_b.elements()[i] * &poly_vec_a.elements()[i]);
        }
        assert_eq!(result, expected);
    }
}
