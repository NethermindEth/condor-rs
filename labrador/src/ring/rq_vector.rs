//! Dynamic vector of `Rq` elements.

use crate::ring::zq::Zq;
use crate::ring::{rq::Rq, Norms};
use core::ops::Mul;
use rand::{CryptoRng, Rng};
use std::ops::Add;

/// Vector of polynomials in Rq
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RqVector {
    elements: Vec<Rq>,
}

impl RqVector {
    pub fn new(elements: Vec<Rq>) -> Self {
        Self { elements }
    }

    /// Create a zero vector
    pub fn zero(length: usize) -> Self {
        Self {
            elements: vec![Rq::zero(); length],
        }
    }

    pub fn set(&mut self, index: usize, val: Rq) {
        self.elements[index] = val;
    }

    pub fn len(&self) -> usize {
        self.elements.len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn elements(&self) -> &Vec<Rq> {
        &self.elements
    }

    pub fn from_zq_vector(elements: Vec<Zq>) -> Self {
        assert!(
            elements.len() % Rq::DEGREE == 0,
            "slice length not multiple of DEGREE"
        );
        let mut result = Vec::new();
        elements.chunks_exact(Rq::DEGREE).for_each(|chunk| {
            result.push(Rq::new(chunk.try_into().unwrap()));
        });
        Self { elements: result }
    }

    /// Random vector of given length with coefficients in `(0, Zq::MAX)`.
    pub fn random<R: Rng + CryptoRng>(rng: &mut R, length: usize) -> Self {
        Self {
            elements: (0..length).map(|_| Rq::random(rng)).collect(),
        }
    }

    /// Random vector of given length with coefficients in `(-bound, bound)`.
    pub fn random_with_bound<R: Rng + CryptoRng>(rng: &mut R, length: usize, bound: u32) -> Self {
        Self {
            elements: (0..length)
                .map(|_| Rq::random_with_bound(rng, bound))
                .collect(),
        }
    }

    /// Flatten the coefficients of every polynomial **in order**.
    pub fn concatenate_coefficients(&self) -> Vec<Zq> {
        let mut concatenated_coeffs = Vec::with_capacity(self.len() * Rq::DEGREE);
        // Iterate over each Rq, extracting the coefficients and concatenating them
        for poly in &self.elements {
            concatenated_coeffs.extend_from_slice(poly.coeffs());
        }
        concatenated_coeffs
    }

    /// Decompose the element to #num_parts number of parts,
    /// where each part's infinity norm is less than or equal to bound/2
    pub fn decompose(&self, b: Zq, parts: usize) -> Vec<RqVector> {
        let mut result = vec![RqVector::zero(self.len()); parts];
        self.elements()
            .iter()
            .enumerate()
            .for_each(|(index, poly)| {
                poly.decompose(b, u64::try_from(parts).expect("parts does not fit in u64"))
                    .into_iter()
                    .zip(result.iter_mut())
                    .for_each(|(decomposed_poly, decomposed_vec)| {
                        decomposed_vec.set(index, decomposed_poly)
                    });
            });
        result
    }
}

impl Add<&RqVector> for &RqVector {
    type Output = RqVector;
    // add two poly vectors
    fn add(self, other: &RqVector) -> RqVector {
        self.elements()
            .iter()
            .zip(other.elements())
            .map(|(a, b)| a + b)
            .collect()
    }
}

impl FromIterator<Rq> for RqVector {
    fn from_iter<T: IntoIterator<Item = Rq>>(iter: T) -> Self {
        let mut elements = Vec::new();
        for item in iter {
            elements.push(item);
        }
        RqVector::new(elements)
    }
}

/// Create a new vector from a `Vec` of elements
impl From<Vec<Rq>> for RqVector {
    fn from(elements: Vec<Rq>) -> Self {
        Self { elements }
    }
}

impl Mul<&Rq> for &RqVector {
    type Output = RqVector;
    // A poly vector multiple by a PolyRing
    fn mul(self, other: &Rq) -> RqVector {
        self.elements().iter().map(|s| s * other).collect()
    }
}

impl Mul<&Zq> for &RqVector {
    type Output = RqVector;
    // A poly vector multiple by a PolyRing
    fn mul(self, other: &Zq) -> RqVector {
        self.elements().iter().map(|s| s * other).collect()
    }
}

impl Mul<Zq> for &RqVector {
    type Output = RqVector;
    // A poly vector multiple by a PolyRing
    fn mul(self, other: Zq) -> RqVector {
        self.elements().iter().map(|s| s * &other).collect()
    }
}

impl Norms for RqVector {
    type NormType = u128;

    // Compute the squared l2 norm of a vector of polynomials
    fn l2_norm_squared(&self) -> Self::NormType {
        self.elements
            .iter()
            .map(|poly| poly.l2_norm_squared())
            .sum()
    }

    fn linf_norm(&self) -> Self::NormType {
        self.elements
            .iter()
            .map(|poly| poly.linf_norm())
            .max()
            .unwrap()
    }
}

#[cfg(test)]
mod tests {
    use rand::rng;

    use super::*;
    use crate::{core::inner_product, ring::rq::tests::generate_rq_from_zq_vector};

    #[test]
    fn test_rqvector_from_iterator() {
        let expected = vec![
            Rq::random(&mut rng()),
            Rq::random(&mut rng()),
            Rq::random(&mut rng()),
        ];
        let vector_of_polynomials = expected.clone().into_iter();
        let result: RqVector = vector_of_polynomials.collect();

        assert_eq!(result.elements(), &expected);
    }

    #[test]
    fn test_rq_vector_multiplication_with_zq() {
        let poly1: Rq = generate_rq_from_zq_vector(vec![Zq::ONE, Zq::new(22)]);
        let poly2: Rq = generate_rq_from_zq_vector(vec![Zq::new(17), Zq::new(12)]);
        let rq_vector = RqVector::from(vec![poly1, poly2]);
        let result = &rq_vector * Zq::new(3);

        assert_eq!(
            result.elements()[0],
            generate_rq_from_zq_vector(vec![Zq::new(3), Zq::new(66)])
        );
        assert_eq!(
            result.elements()[1],
            generate_rq_from_zq_vector(vec![Zq::new(51), Zq::new(36)])
        );

        #[allow(clippy::op_ref)]
        let result2 = &rq_vector * &Zq::new(3);
        assert_eq!(
            result2.elements()[0],
            generate_rq_from_zq_vector(vec![Zq::new(3), Zq::new(66)])
        );
        assert_eq!(
            result2.elements()[1],
            generate_rq_from_zq_vector(vec![Zq::new(51), Zq::new(36)])
        );
    }

    #[test]
    fn test_rqvector_mul() {
        let poly1: Rq = generate_rq_from_zq_vector(vec![Zq::ONE, Zq::new(2)]);
        let poly2: Rq = generate_rq_from_zq_vector(vec![Zq::ONE, Zq::new(4)]);
        let vec_1: RqVector = RqVector::from(vec![poly1]);
        let vec_2: RqVector = RqVector::from(vec![poly2]);
        let result = inner_product::compute_linear_combination(vec_1.elements(), vec_2.elements());
        let poly_exp: Rq = generate_rq_from_zq_vector(vec![Zq::new(1), Zq::new(6), Zq::new(8)]);
        assert_eq!(result, poly_exp);

        let poly3: Rq = generate_rq_from_zq_vector(vec![Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE]);
        let poly4: Rq = generate_rq_from_zq_vector(vec![Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE]);
        let vec_3: RqVector = RqVector::from(vec![poly3]);
        let vec_4: RqVector = RqVector::from(vec![poly4]);
        let result_1 =
            inner_product::compute_linear_combination(vec_3.elements(), vec_4.elements());
        let poly_exp_1: Rq = generate_rq_from_zq_vector(vec![
            Zq::new(1),
            Zq::new(2),
            Zq::new(3),
            Zq::new(4),
            Zq::new(3),
            Zq::new(2),
            Zq::new(1),
        ]);
        assert_eq!(result_1, poly_exp_1);

        let poly5: Rq = generate_rq_from_zq_vector(vec![Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE]);
        let poly6: Rq = generate_rq_from_zq_vector(vec![Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE]);
        let poly7: Rq = generate_rq_from_zq_vector(vec![Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE]);
        let poly8: Rq = generate_rq_from_zq_vector(vec![Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE]);
        let vec_5: RqVector = RqVector::from(vec![poly5, poly6]);
        let vec_6: RqVector = RqVector::from(vec![poly7, poly8]);
        let result_2 =
            inner_product::compute_linear_combination(vec_5.elements(), vec_6.elements());
        let poly_exp_2: Rq = generate_rq_from_zq_vector(vec![
            Zq::new(2),
            Zq::new(4),
            Zq::new(6),
            Zq::new(8),
            Zq::new(6),
            Zq::new(4),
            Zq::new(2),
        ]);
        assert_eq!(result_2, poly_exp_2);
    }
}

#[cfg(test)]
mod norm_tests {
    use super::*;
    use crate::ring::{rq::tests::generate_rq_from_zq_vector, Norms};

    // Test the square of the norm
    #[test]
    fn test_l2_norm() {
        let poly1 = generate_rq_from_zq_vector(vec![
            Zq::ONE,
            Zq::ZERO,
            Zq::new(5),
            Zq::NEG_ONE - Zq::new(1),
        ]);
        let poly2 = generate_rq_from_zq_vector(vec![Zq::ZERO, Zq::ZERO, Zq::new(5), Zq::ONE]);
        let poly_vec1: RqVector = vec![poly1.clone(), poly2.clone()].into();
        assert_eq!(
            poly_vec1.l2_norm_squared(),
            poly1.l2_norm_squared() + poly2.l2_norm_squared()
        );

        let zero_vec: RqVector = RqVector::zero(4);
        assert_eq!(zero_vec.l2_norm_squared(), 0);
    }

    #[test]
    fn test_linf_norm() {
        let poly1 = generate_rq_from_zq_vector(vec![
            Zq::ONE,
            Zq::new(500),
            Zq::new(5),
            Zq::NEG_ONE - Zq::new(1000),
        ]);
        let poly2 =
            generate_rq_from_zq_vector(vec![Zq::new(5), Zq::new(5), Zq::new(5000), Zq::new(5)]);
        let poly_vec1: RqVector = vec![poly1.clone(), poly2.clone()].into();
        assert_eq!(poly_vec1.linf_norm(), 5000);
        let poly_vec2: RqVector = vec![poly2.clone(), poly1.clone()].into();
        assert_eq!(poly_vec2.linf_norm(), 5000);

        let zero_vec: RqVector = RqVector::zero(4);
        assert_eq!(zero_vec.linf_norm(), 0);
    }
}

#[cfg(test)]
mod decomposition_tests {

    use rand::rng;

    use super::*;

    fn default_base_and_parts() -> (Zq, usize) {
        let base = Zq::new(2u32.pow(16));
        let parts = 3; // 32 / 16 + (additional parts)
        (base, parts)
    }

    #[test]
    fn test_decomposition_number_of_parts() {
        let rq_vector = RqVector::random(&mut rng(), 10);
        let (base, parts) = default_base_and_parts();
        let decomposed = rq_vector.decompose(base, parts);
        assert_eq!(decomposed.len(), parts);
    }

    #[test]
    fn test_decompositions_linf_norm() {
        let rq_vector = RqVector::random(&mut rng(), 10);
        let (base, parts) = default_base_and_parts();
        let decomposed_vector = rq_vector.decompose(base, parts);

        for decomposition in decomposed_vector {
            assert!(
                decomposition.linf_norm() <= base.div_floor_by(2).to_u128(),
                "decomposition.linf_norm(): {} and base.div_floor_by(2): {}",
                decomposition.linf_norm(),
                base.div_floor_by(2)
            );
        }
    }

    #[test]
    fn test_decompositions_computation_is_correct() {
        let rq_vector = RqVector::random(&mut rng(), 10);
        let (base, parts) = default_base_and_parts();
        let decomposed_vector = rq_vector.decompose(base, parts);

        let mut combined_vector = RqVector::zero(10);

        let mut exponential_base = Zq::new(1);
        for element in decomposed_vector {
            combined_vector = &combined_vector + &(&element * exponential_base);
            exponential_base *= base;
        }

        assert_eq!(combined_vector, rq_vector);
    }
}
