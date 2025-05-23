use crate::ring::rq::Rq;
use crate::ring::zq::Zq;
use core::ops::{Index, IndexMut, Mul};
use core::slice::Iter;
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

    pub fn set(&mut self, index: usize, value: Rq) {
        self.elements[index] = value;
    }

    pub fn into_inner(self) -> Vec<Rq> {
        self.elements
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

    /// Create a random vector
    pub fn random_ternary<R: Rng + CryptoRng>(rng: &mut R, length: usize) -> Self {
        Self {
            elements: (0..length).map(|_| Rq::random_ternary(rng)).collect(),
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

    pub fn inner_product_poly_vector(&self, other: &RqVector) -> Rq {
        self.iter()
            .zip(other.iter())
            .map(|(a, b)| a * b)
            .fold(Rq::zero(), |acc, val| &acc + &val)
    }

    /// Get the underlying vector as slice
    pub fn as_slice(&self) -> &[Rq] {
        &self.elements
    }

    pub fn iter(&self) -> Iter<'_, Rq> {
        self.elements.iter()
    }

    // Compute the squared norm of a vector of polynomials
    pub fn compute_norm_squared(&self) -> Zq {
        self.elements
            .iter()
            .flat_map(|poly| poly.get_coefficients()) // Collect coefficients from all polynomials
            .map(|coeff| *coeff * *coeff)
            .sum()
    }

    pub fn decompose(&self, b: Zq, parts: usize) -> Vec<RqVector> {
        self.iter().map(|i| Rq::decompose(i, b, parts)).collect()
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

impl Add<&RqVector> for &RqVector {
    type Output = RqVector;
    // add two poly vectors
    fn add(self, other: &RqVector) -> RqVector {
        self.iter().zip(other.iter()).map(|(a, b)| a + b).collect()
    }
}

// Enable array-like indexing
impl Index<usize> for RqVector {
    type Output = Rq;

    fn index(&self, index: usize) -> &Self::Output {
        &self.elements[index]
    }
}

impl IndexMut<usize> for RqVector {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.elements[index]
    }
}

/// Create a new vector from a `Vec` of elements
impl From<Vec<Rq>> for RqVector {
    fn from(elements: Vec<Rq>) -> Self {
        Self { elements }
    }
}

impl Mul<&Vec<RqVector>> for &RqVector {
    type Output = RqVector;
    fn mul(self, other: &Vec<RqVector>) -> RqVector {
        other
            .iter()
            .map(|o| o.inner_product_poly_vector(self))
            .collect()
    }
}

impl Mul<&Rq> for &RqVector {
    type Output = RqVector;
    // A poly vector multiple by a PolyRing
    fn mul(self, other: &Rq) -> RqVector {
        self.iter().map(|s| s * other).collect()
    }
}

impl Mul<&Zq> for &RqVector {
    type Output = RqVector;
    // A poly vector multiple by a PolyRing
    fn mul(self, other: &Zq) -> RqVector {
        self.iter().map(|s| s * other).collect()
    }
}

impl Mul<Zq> for &RqVector {
    type Output = RqVector;
    // A poly vector multiple by a PolyRing
    fn mul(self, other: Zq) -> RqVector {
        self.iter().map(|s| s * other).collect()
    }
}

// Dot product between two vectors
impl Mul for &RqVector {
    type Output = Rq;
    fn mul(self, rhs: Self) -> Self::Output {
        self.elements
            .iter()
            .zip(rhs.elements.iter())
            .map(|(a, b)| a * b)
            .fold(Rq::zero(), |acc, x| &acc + &x)
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_rqvector_mul() {
        let poly1: Rq = vec![Zq::ONE, Zq::new(2)].into();
        let poly2: Rq = vec![Zq::ONE, Zq::new(4)].into();
        let vec_1: RqVector = RqVector::from(vec![poly1]);
        let vec_2: RqVector = RqVector::from(vec![poly2]);
        let result = vec_1.mul(&vec_2);
        let poly_exp: Rq = vec![Zq::new(1), Zq::new(6), Zq::new(8)].into();
        assert_eq!(result, poly_exp);

        let poly3: Rq = vec![Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE].into();
        let poly4: Rq = vec![Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE].into();
        let vec_3: RqVector = RqVector::from(vec![poly3]);
        let vec_4: RqVector = RqVector::from(vec![poly4]);
        let result_1 = vec_3.mul(&vec_4);
        let poly_exp_1: Rq = vec![
            Zq::new(1),
            Zq::new(2),
            Zq::new(3),
            Zq::new(4),
            Zq::new(3),
            Zq::new(2),
            Zq::new(1),
        ]
        .into();
        assert_eq!(result_1, poly_exp_1);

        let poly5: Rq = vec![Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE].into();
        let poly6: Rq = vec![Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE].into();
        let poly7: Rq = vec![Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE].into();
        let poly8: Rq = vec![Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE].into();
        let vec_5: RqVector = RqVector::from(vec![poly5, poly6]);
        let vec_6: RqVector = RqVector::from(vec![poly7, poly8]);
        let result_2 = vec_5.mul(&vec_6);
        let poly_exp_2: Rq = vec![
            Zq::new(2),
            Zq::new(4),
            Zq::new(6),
            Zq::new(8),
            Zq::new(6),
            Zq::new(4),
            Zq::new(2),
        ]
        .into();
        assert_eq!(result_2, poly_exp_2);
    }

    // Test the square of the norm
    #[test]
    fn test_norm() {
        let poly: RqVector = vec![
            vec![Zq::ONE, Zq::ZERO, Zq::new(5), Zq::MAX].into(),
            vec![Zq::ZERO, Zq::ZERO, Zq::new(5), Zq::ONE].into(),
        ]
        .into();
        let result = Zq::new(53);
        assert!(poly.compute_norm_squared() == result);

        let poly2: RqVector = vec![vec![Zq::new(5), Zq::ONE, Zq::MAX, Zq::ZERO].into()].into();
        let result2 = Zq::new(27);
        assert!(poly2.compute_norm_squared() == result2);

        let poly_zero: RqVector = RqVector::zero(4);
        let result_zero = Zq::ZERO;
        assert!(poly_zero.compute_norm_squared() == result_zero);
    }
}
