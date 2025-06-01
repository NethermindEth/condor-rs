use std::{
    borrow::Borrow,
    ops::{Add, Mul},
};

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
    fn test_inputs_with_different_legths_panic() {
        let vector_a = RqVector::random(&mut rng(), 100);
        let vector_b = RqVector::random(&mut rng(), 90);
        let _ = compute_linear_combination(vector_a.get_elements(), vector_b.get_elements());
    }

    #[test]
    fn test_zq_type_inputs() {
        let poly_a = Rq::random(&mut rng());
        let poly_b = Rq::random(&mut rng());

        let result =
            compute_linear_combination(poly_a.get_coefficients(), poly_b.get_coefficients());

        let mut expected = Zq::ZERO;
        for i in 0..poly_a.get_coefficients().len() {
            expected += poly_a.get_coefficients()[i] * poly_b.get_coefficients()[i];
        }
        assert_eq!(result, expected);
    }

    #[test]
    fn test_rq_type_inputs() {
        let poly_vec_a = RqVector::random(&mut rng(), 100);
        let poly_vec_b = RqVector::random(&mut rng(), 100);

        let result =
            compute_linear_combination(poly_vec_a.get_elements(), poly_vec_b.get_elements());

        let mut expected = Rq::zero();
        for i in 0..poly_vec_a.get_elements().len() {
            expected = &expected + &(&poly_vec_a.get_elements()[i] * &poly_vec_b.get_elements()[i]);
        }
        assert_eq!(result, expected);
    }

    #[test]
    fn test_zq_rq_type_inputs() {
        let poly_a = Rq::random(&mut rng());
        let poly_vec_b = RqVector::random(&mut rng(), 64);

        let result =
            compute_linear_combination(poly_vec_b.get_elements(), poly_a.get_coefficients());

        let mut expected = Rq::zero();
        for i in 0..poly_a.get_coefficients().len() {
            expected = &expected + &(&poly_vec_b.get_elements()[i] * &poly_a.get_coefficients()[i]);
        }
        assert_eq!(result, expected);
    }

    #[test]
    fn test_zq_rqvector_type_inputs() {
        let poly_a = Rq::random(&mut rng());
        let poly_vec_b = RqMatrix::random(&mut rng(), 64, 100);

        let result =
            compute_linear_combination(poly_vec_b.get_elements(), poly_a.get_coefficients());

        let mut expected = RqVector::zero(100);
        for i in 0..poly_a.get_coefficients().len() {
            expected = &expected + &(&poly_vec_b.get_elements()[i] * poly_a.get_coefficients()[i]);
        }
        assert_eq!(result, expected);
    }

    #[test]
    fn test_rq_rqvector_type_inputs() {
        let poly_vec_a = RqVector::random(&mut rng(), 80);
        let poly_vec_b = RqMatrix::random(&mut rng(), 80, 100);

        let result =
            compute_linear_combination(poly_vec_b.get_elements(), poly_vec_a.get_elements());

        let mut expected = RqVector::zero(100);
        for i in 0..poly_vec_a.get_elements().len() {
            expected = &expected + &(&poly_vec_b.get_elements()[i] * &poly_vec_a.get_elements()[i]);
        }
        assert_eq!(result, expected);
    }
}
