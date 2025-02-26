use crate::rq::Rq;
use crate::rq_vector::RqVector;
use std::ops::Mul;

/// Matrix of polynomials in Rq
#[derive(Debug, Clone)]
pub struct RqMatrix<const M: usize, const N: usize, const D: usize> {
    elements: [[Rq<D>; N]; M],
}

impl<const M: usize, const N: usize, const D: usize> RqMatrix<M, N, D> {
    /// Create a random matrix of polynomials with small coefficients
    pub fn random() -> Self {
        let matrix = std::array::from_fn(|_| std::array::from_fn(|_| Rq::random_small()));
        Self { elements: matrix }
    }

    /// Matrix-vector multiplication that returns an array
    pub fn mul_vec(&self, vec: &[Rq<D>; N]) -> [Rq<D>; M] {
        (self * &RqVector::from(vec.clone())).into_array()
    }
}

// Implement matrix-vector multiplication for reference to matrix
impl<const M: usize, const N: usize, const D: usize> Mul<&RqVector<N, D>> for &RqMatrix<M, N, D> {
    type Output = RqVector<M, D>;

    fn mul(self, rhs: &RqVector<N, D>) -> Self::Output {
        let mut result = RqVector::zero();

        for (i, row) in self.elements.iter().enumerate() {
            result[i] = &RqVector::from(row.clone()) * rhs;
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
