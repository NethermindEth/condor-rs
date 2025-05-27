use crate::ring::rq::Rq;
use crate::ring::zq::Zq;
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

    pub fn new_from_zq_vector(elements: Vec<Zq>) -> Self {
        let mut result = Vec::new();
        debug_assert!(elements.len() % Rq::DEGREE == 0);
        elements.chunks_exact(Rq::DEGREE).for_each(|chunk| {
            result.push(Rq::new(chunk.try_into().unwrap()));
        });

        RqVector { elements: result }
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

    pub fn get_length(&self) -> usize {
        self.elements.len()
    }

    pub fn get_elements(&self) -> &Vec<Rq> {
        &self.elements
    }

    /// Create a random vector
    pub fn random<R: Rng + CryptoRng>(rng: &mut R, length: usize) -> Self {
        Self {
            elements: (0..length).map(|_| Rq::random(rng)).collect(),
        }
    }

    // Create a random vector
    pub fn random_ternary<R: Rng + CryptoRng>(rng: &mut R, length: usize) -> Self {
        Self {
            elements: (0..length).map(|_| Rq::random_ternary(rng)).collect(),
        }
    }

    pub fn random_with_bound<R: Rng + CryptoRng>(rng: &mut R, length: usize, bound: u32) -> Self {
        Self {
            elements: (0..length)
                .map(|_| Rq::random_with_bound(rng, bound))
                .collect(),
        }
    }

    /// Function to concatenate coefficients from multiple Rq into a Vec<Zq>
    pub fn concatenate_coefficients(&self) -> Vec<Zq> {
        let total_coeffs = self.elements.len() * Rq::DEGREE;
        let mut concatenated_coeffs: Vec<Zq> = Vec::with_capacity(total_coeffs);
        // Iterate over each Rq, extracting the coefficients and concatenating them
        for rq in &self.elements {
            let coeffs = rq.get_coefficients();
            concatenated_coeffs.extend_from_slice(coeffs);
        }
        concatenated_coeffs
    }

    // Compute the squared norm of a vector of polynomials
    pub fn l2_norm_squared(&self) -> Zq {
        self.elements
            .iter()
            .map(|poly| poly.l2_norm_squared()) // Collect coefficients from all polynomials
            .sum()
    }

    pub fn decompose(&self, b: Zq, parts: usize) -> Vec<RqVector> {
        self.get_elements()
            .iter()
            .map(|i| RqVector::new(Rq::decompose(i, b, parts)))
            .collect()
    }
}

impl Add<&RqVector> for &RqVector {
    type Output = RqVector;
    // add two poly vectors
    fn add(self, other: &RqVector) -> RqVector {
        self.get_elements()
            .iter()
            .zip(other.get_elements())
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
        self.get_elements().iter().map(|s| s * other).collect()
    }
}

impl Mul<&Zq> for &RqVector {
    type Output = RqVector;
    // A poly vector multiple by a PolyRing
    fn mul(self, other: &Zq) -> RqVector {
        self.get_elements().iter().map(|s| s * other).collect()
    }
}

impl Mul<Zq> for &RqVector {
    type Output = RqVector;
    // A poly vector multiple by a PolyRing
    fn mul(self, other: Zq) -> RqVector {
        self.get_elements().iter().map(|s| s * &other).collect()
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

        assert_eq!(result.get_elements(), &expected);
    }

    #[test]
    fn test_rq_vector_multiplication_with_zq() {
        let poly1: Rq = generate_rq_from_zq_vector(vec![Zq::ONE, Zq::new(22)]);
        let poly2: Rq = generate_rq_from_zq_vector(vec![Zq::new(17), Zq::new(12)]);
        let rq_vector = RqVector::from(vec![poly1, poly2]);
        let result = &rq_vector * Zq::new(3);

        assert_eq!(
            result.get_elements()[0],
            generate_rq_from_zq_vector(vec![Zq::new(3), Zq::new(66)])
        );
        assert_eq!(
            result.get_elements()[1],
            generate_rq_from_zq_vector(vec![Zq::new(51), Zq::new(36)])
        );

        #[allow(clippy::op_ref)]
        let result2 = &rq_vector * &Zq::new(3);
        assert_eq!(
            result2.get_elements()[0],
            generate_rq_from_zq_vector(vec![Zq::new(3), Zq::new(66)])
        );
        assert_eq!(
            result2.get_elements()[1],
            generate_rq_from_zq_vector(vec![Zq::new(51), Zq::new(36)])
        );
    }

    #[test]
    fn test_rqvector_mul() {
        let poly1: Rq = generate_rq_from_zq_vector(vec![Zq::ONE, Zq::new(2)]);
        let poly2: Rq = generate_rq_from_zq_vector(vec![Zq::ONE, Zq::new(4)]);
        let vec_1: RqVector = RqVector::from(vec![poly1]);
        let vec_2: RqVector = RqVector::from(vec![poly2]);
        let result =
            inner_product::compute_linear_combination(vec_1.get_elements(), vec_2.get_elements());
        let poly_exp: Rq = generate_rq_from_zq_vector(vec![Zq::new(1), Zq::new(6), Zq::new(8)]);
        assert_eq!(result, poly_exp);

        let poly3: Rq = generate_rq_from_zq_vector(vec![Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE]);
        let poly4: Rq = generate_rq_from_zq_vector(vec![Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE]);
        let vec_3: RqVector = RqVector::from(vec![poly3]);
        let vec_4: RqVector = RqVector::from(vec![poly4]);
        let result_1 =
            inner_product::compute_linear_combination(vec_3.get_elements(), vec_4.get_elements());
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
            inner_product::compute_linear_combination(vec_5.get_elements(), vec_6.get_elements());
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

    // Test the square of the norm
    #[test]
    fn test_norm() {
        let poly1 =
            generate_rq_from_zq_vector(vec![Zq::ONE, Zq::ZERO, Zq::new(5), Zq::MAX - Zq::new(1)]);
        let poly2 = generate_rq_from_zq_vector(vec![Zq::ZERO, Zq::ZERO, Zq::new(5), Zq::ONE]);
        let poly_vec1: RqVector = vec![poly1.clone(), poly2.clone()].into();
        assert_eq!(
            poly_vec1.l2_norm_squared(),
            poly1.l2_norm_squared() + poly2.l2_norm_squared()
        );

        let zero_vec: RqVector = RqVector::zero(4);
        assert_eq!(zero_vec.l2_norm_squared(), Zq::ZERO);
    }
}
