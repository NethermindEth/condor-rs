use crate::rq::Rq;
use crate::zq::Zq;
use core::ops::{Index, IndexMut, Mul};
use core::slice::Iter;
use rand::{CryptoRng, Rng};

/// Vector of polynomials in Rq
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RqVector<const N: usize, const D: usize> {
    elements: Vec<Rq<D>>,
}

impl<const N: usize, const D: usize> RqVector<N, D> {
    /// Create a zero vector
    pub fn zero() -> Self {
        Self {
            elements: vec![Rq::zero(); N],
        }
    }

    /// Create a random vector
    pub fn random<R: Rng + CryptoRng>(rng: &mut R) -> Self {
        Self {
            elements: (0..N).map(|_| Rq::random(rng)).collect(),
        }
    }

    /// Create a random vector
    pub fn random_small<R: Rng + CryptoRng>(rng: &mut R) -> Self {
        Self {
            elements: (0..N).map(|_| Rq::random_ternary(rng)).collect(),
        }
    }

    /// Function to concatenate coefficients from multiple Rq into a Vec<Zq>
    pub fn concatenate_coefficients(&self) -> Vec<Zq> {
        let total_coeffs = self.elements.len() * D;
        let mut concatenated_coeffs: Vec<Zq> = Vec::with_capacity(total_coeffs);
        // Iterate over each Rq, extracting the coefficients and concatenating them
        for rq in &self.elements {
            let coeffs = rq.get_coefficients();
            concatenated_coeffs.extend_from_slice(coeffs);
        }

        concatenated_coeffs
    }

    /// Get the underlying vector as slice
    pub fn as_slice(&self) -> &[Rq<D>] {
        &self.elements
    }

    pub fn iter(&self) -> Iter<'_, Rq<D>> {
        self.elements.iter()
    }

    /// Convert into array
    pub fn into_array(self) -> [Rq<D>; N] {
        self.elements
            .try_into()
            .unwrap_or_else(|_| panic!("Vector length mismatch"))
    }

    // Compute the squared norm of a vector of polynomials
    pub fn compute_norm_squared(&self) -> Zq {
        self.elements
            .iter()
            .flat_map(|poly| poly.get_coefficients()) // Collect coefficients from all polynomials
            .map(|coeff| *coeff * *coeff)
            .sum()
    }
}

// Enable array-like indexing
impl<const N: usize, const D: usize> Index<usize> for RqVector<N, D> {
    type Output = Rq<D>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.elements[index]
    }
}

impl<const N: usize, const D: usize> IndexMut<usize> for RqVector<N, D> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.elements[index]
    }
}

/// Create a new vector from a `Vec` of elements
impl<const N: usize, const D: usize> From<Vec<Rq<D>>> for RqVector<N, D> {
    fn from(elements: Vec<Rq<D>>) -> Self {
        assert_eq!(elements.len(), N, "Vector length must match N");
        Self { elements }
    }
}

// Dot product between two vectors
impl<const N: usize, const D: usize> Mul for &RqVector<N, D> {
    type Output = Rq<D>;

    fn mul(self, rhs: Self) -> Self::Output {
        self.elements
            .iter()
            .zip(rhs.elements.iter())
            .map(|(a, b)| a.clone() * b.clone())
            .fold(Rq::zero(), |acc, x| acc + x)
    }
}

#[cfg(test)]
mod tests {
    use crate::{rq::Rq, zq::Zq};

    use super::*;

    #[test]
    fn test_rqvector_mul() {
        let poly1: Rq<2> = vec![Zq::ONE, Zq::new(2)].into();
        let poly2: Rq<2> = vec![Zq::ONE, Zq::new(4)].into();
        let vec_1: RqVector<1, 2> = RqVector::from(vec![poly1]);
        let vec_2: RqVector<1, 2> = RqVector::from(vec![poly2]);
        let result = vec_1.mul(&vec_2);
        let poly_exp: Rq<2> = vec![Zq::new(u32::MAX - 6), Zq::new(6)].into();
        assert_eq!(result, poly_exp);

        let poly3: Rq<4> = vec![Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE].into();
        let poly4: Rq<4> = vec![Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE].into();
        let vec_3: RqVector<1, 4> = RqVector::from(vec![poly3]);
        let vec_4: RqVector<1, 4> = RqVector::from(vec![poly4]);
        let result_1 = vec_3.mul(&vec_4);
        let poly_exp_1: Rq<4> =
            vec![Zq::new(u32::MAX - 1), Zq::ZERO, Zq::new(2), Zq::new(4)].into();
        assert_eq!(result_1, poly_exp_1);

        let poly5: Rq<4> = vec![Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE].into();
        let poly6: Rq<4> = vec![Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE].into();
        let poly7: Rq<4> = vec![Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE].into();
        let poly8: Rq<4> = vec![Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE].into();
        let vec_5: RqVector<2, 4> = RqVector::from(vec![poly5, poly6]);
        let vec_6: RqVector<2, 4> = RqVector::from(vec![poly7, poly8]);
        let result_2 = vec_5.mul(&vec_6);
        let poly_exp_2: Rq<4> =
            vec![Zq::new(u32::MAX - 3), Zq::ZERO, Zq::new(4), Zq::new(8)].into();
        assert_eq!(result_2, poly_exp_2);
    }

    #[test]
    fn test_rqvector_mul() {
        let poly1: Rq<2> = vec![Zq::ONE, Zq::new(2)].into();
        let poly2: Rq<2> = vec![Zq::ONE, Zq::new(4)].into();
        let vec_1: RqVector<1, 2> = RqVector::from(vec![poly1]);
        let vec_2: RqVector<1, 2> = RqVector::from(vec![poly2]);
        let result = vec_1.mul(&vec_2);
        let poly_exp: Rq<2> = vec![Zq::new(u32::MAX - 6), Zq::new(6)].into();
        assert_eq!(result, poly_exp);

        let poly3: Rq<4> = vec![Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE].into();
        let poly4: Rq<4> = vec![Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE].into();
        let vec_3: RqVector<1, 4> = RqVector::from(vec![poly3]);
        let vec_4: RqVector<1, 4> = RqVector::from(vec![poly4]);
        let result_1 = vec_3.mul(&vec_4);
        let poly_exp_1: Rq<4> =
            vec![Zq::new(u32::MAX - 1), Zq::ZERO, Zq::new(2), Zq::new(4)].into();
        assert_eq!(result_1, poly_exp_1);

        let poly5: Rq<4> = vec![Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE].into();
        let poly6: Rq<4> = vec![Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE].into();
        let poly7: Rq<4> = vec![Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE].into();
        let poly8: Rq<4> = vec![Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE].into();
        let vec_5: RqVector<2, 4> = RqVector::from(vec![poly5, poly6]);
        let vec_6: RqVector<2, 4> = RqVector::from(vec![poly7, poly8]);
        let result_2 = vec_5.mul(&vec_6);
        let poly_exp_2: Rq<4> =
            vec![Zq::new(u32::MAX - 3), Zq::ZERO, Zq::new(4), Zq::new(8)].into();
        assert_eq!(result_2, poly_exp_2);
    }
}
