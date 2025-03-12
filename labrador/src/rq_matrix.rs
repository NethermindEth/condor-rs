use crate::rq_vector::RqVector;
use rand::{CryptoRng, Rng};
use std::ops::Mul;

/// Matrix of polynomials in Rq
#[derive(Debug, Clone)]
pub struct RqMatrix<const M: usize, const N: usize, const D: usize> {
    elements: Vec<RqVector<N, D>>,
}

impl<const M: usize, const N: usize, const D: usize> RqMatrix<M, N, D> {
    /// Constructor for the Matrix of polynomials in Rq
    pub const fn new(elements: Vec<RqVector<N, D>>) -> Self {
        RqMatrix { elements }
    }

    /// Create a random matrix of polynomials
    pub fn random<R: Rng + CryptoRng>(rng: &mut R) -> Self {
        Self {
            elements: (0..M).map(|_| RqVector::random(rng)).collect(),
        }
    }

    /// Create a random matrix of polynomials with small coefficients
    pub fn random_small<R: Rng + CryptoRng>(rng: &mut R) -> Self {
        Self {
            elements: (0..M).map(|_| RqVector::random_small(rng)).collect(),
        }
    }
}

// Implement matrix-vector multiplication for reference to matrix
impl<const M: usize, const N: usize, const D: usize> Mul<&RqVector<N, D>> for &RqMatrix<M, N, D> {
    type Output = RqVector<M, D>;

    fn mul(self, rhs: &RqVector<N, D>) -> Self::Output {
        let mut result = RqVector::zero();

        for (i, row) in self.elements.iter().enumerate() {
            result[i] = row * rhs;
        }

        result
    }
}

// Implement matrix-vector multiplication for owned matrix by delegating to reference implementation
impl<const M: usize, const N: usize, const D: usize> Mul<&RqVector<N, D>> for RqMatrix<M, N, D> {
    type Output = RqVector<M, D>;

    fn mul(self, rhs: &RqVector<N, D>) -> Self::Output {
        &self * rhs
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{rq::Rq, zq::Zq};

    #[test]
    #[cfg(not(feature = "skip-slow-tests"))]
    fn rqmatrix_fits_stack() {
        let mut rng = rand::rng();
        let _: RqMatrix<256, { 1 << 10 }, 64> = RqMatrix::random(&mut rng);
    }

    #[test]
    fn test_rqmartrix_mul() {
        let poly1: Rq<2> = vec![Zq::new(8), Zq::new(6)].into();
        let poly2: Rq<2> = vec![Zq::new(u32::MAX - 4), Zq::new(u32::MAX - 4)].into();
        let poly3: Rq<2> = vec![Zq::ONE, Zq::ZERO].into();
        let poly4: Rq<2> = vec![Zq::ZERO, Zq::new(4)].into();
        let matrix_1: RqMatrix<1, 2, 2> = RqMatrix::new(vec![RqVector::from(vec![poly1, poly2])]);
        let vec_1: RqVector<2, 2> = RqVector::from(vec![poly3, poly4]);

        let result_1 = matrix_1.mul(&vec_1);
        let expected_poly_1 = vec![Zq::new(28), Zq::new(u32::MAX - 13)].into();
        let expected_1 = RqVector::from(vec![expected_poly_1]);
        assert_eq!(result_1, expected_1);

        let poly5: Rq<2> = vec![Zq::new(u32::MAX - 6), Zq::new(7)].into();
        let poly6: Rq<2> = vec![Zq::new(u32::MAX - 2), Zq::ZERO].into();
        let poly7: Rq<2> = vec![Zq::new(8), Zq::new(u32::MAX - 1)].into();
        let poly8: Rq<2> = vec![Zq::new(u32::MAX - 3), Zq::new(4)].into();
        let poly9: Rq<2> = vec![Zq::MAX, Zq::new(u32::MAX - 1)].into();
        let poly10: Rq<2> = vec![Zq::new(u32::MAX - 2), Zq::new(u32::MAX - 2)].into();
        let matrix_2: RqMatrix<2, 2, 2> = RqMatrix::new(vec![
            RqVector::from(vec![poly5, poly6]),
            RqVector::from(vec![poly7, poly8]),
        ]);
        let vec_2: RqVector<2, 2> = RqVector::from(vec![poly9, poly10]);

        let result_2 = matrix_2.mul(&vec_2);
        let expected_poly_2_1 = vec![Zq::new(30), Zq::new(16)].into();
        let expected_poly_2_2 = vec![Zq::new(12), Zq::new(u32::MAX - 13)].into();
        let expected_2 = RqVector::from(vec![expected_poly_2_1, expected_poly_2_2]);
        assert_eq!(result_2, expected_2);
    }
}
