use crate::ring::rq_vector::RqVector;
use rand::{CryptoRng, Rng};
use std::ops::Mul;

use super::{rq::Rq, zq::Zq};

/// Matrix of polynomials in Rq
#[derive(Debug, Clone)]
pub struct RqMatrix {
    elements: Vec<RqVector>,
}

impl RqMatrix {
    /// Constructor for the Matrix of polynomials in Rq
    pub const fn new(elements: Vec<RqVector>) -> Self {
        RqMatrix { elements }
    }

    pub fn zero(row_len: usize, col_len: usize) -> Self {
        RqMatrix::new(vec![RqVector::zero(col_len); row_len])
    }

    pub fn get_row_len(&self) -> usize {
        self.elements.len()
    }

    pub fn get_col_len(&self) -> usize {
        self.elements[0].get_length()
    }

    pub fn get_cell(&self, row: usize, col: usize) -> &Rq {
        &self.elements[row].get_elements()[col]
    }

    pub fn set_sell(&mut self, row: usize, col: usize, value: Rq) {
        self.elements[row].set(col, value);
    }

    /// Create a random matrix of polynomials
    pub fn random<R: Rng + CryptoRng>(rng: &mut R, row_len: usize, col_len: usize) -> Self {
        Self {
            elements: (0..row_len)
                .map(|_| RqVector::random(rng, col_len))
                .collect(),
        }
    }

    pub fn get_cell_symmetric(&self, row: usize, col: usize) -> &Rq {
        if row >= col {
            &self.elements[row].get_elements()[col]
        } else {
            &self.elements[col].get_elements()[row]
        }
    }

    /// Create a random matrix of polynomials with ternary coefficients
    pub fn random_ternary<R: Rng + CryptoRng>(rng: &mut R, row_len: usize, col_len: usize) -> Self {
        Self {
            elements: (0..row_len)
                .map(|_| RqVector::random_ternary(rng, col_len))
                .collect(),
        }
    }

    pub fn get_elements(&self) -> &Vec<RqVector> {
        &self.elements
    }

    pub fn decompose_each_cell(&self, base: Zq, num_parts: usize) -> RqVector {
        let mut decomposed_vec = Vec::new();
        for ring_vector in self.get_elements() {
            for ring in ring_vector.get_elements() {
                decomposed_vec.extend(ring.decompose(base, num_parts).into_inner());
            }
        }
        RqVector::new(decomposed_vec)
    }
}

impl FromIterator<RqVector> for RqMatrix {
    fn from_iter<T: IntoIterator<Item = RqVector>>(iter: T) -> Self {
        let mut elements = Vec::new();
        for item in iter {
            elements.push(item);
        }
        RqMatrix::new(elements)
    }
}

// Implement matrix-vector multiplication for reference to matrix
impl Mul<&RqVector> for &RqMatrix {
    type Output = RqVector;

    fn mul(self, rhs: &RqVector) -> Self::Output {
        let mut result = RqVector::zero(self.elements.len());

        for (i, row) in self.elements.iter().enumerate() {
            result[i] = row * rhs;
        }

        result
    }
}

// Implement matrix-vector multiplication for owned matrix by delegating to reference implementation
impl Mul<&RqVector> for RqMatrix {
    type Output = RqVector;

    fn mul(self, rhs: &RqVector) -> Self::Output {
        &self * rhs
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ring::rq::Rq;
    use crate::ring::zq::Zq;

    #[test]
    #[cfg(not(feature = "skip-slow-tests"))]
    fn rqmatrix_fits_stack() {
        let mut rng = rand::rng();
        let _: RqMatrix = RqMatrix::random(&mut rng, 256, 1 << 10);
    }

    #[test]
    fn test_zero_matrix() {
        let matrix = RqMatrix::zero(10, 20);
        assert_eq!(matrix.get_row_len(), 10);
        assert_eq!(matrix.get_col_len(), 20);
        for row in matrix.get_elements() {
            for cell in row.get_elements() {
                assert!(cell.is_zero());
            }
        }
    }

    #[test]
    fn test_rqmartrix_mul() {
        let poly1: Rq = vec![Zq::new(8), Zq::new(6)].into();
        let poly2: Rq = vec![Zq::new(u32::MAX - 4), Zq::new(u32::MAX - 4)].into();
        let poly3: Rq = vec![Zq::ONE, Zq::ZERO].into();
        let poly4: Rq = vec![Zq::ZERO, Zq::new(4)].into();
        let matrix_1: RqMatrix = RqMatrix::new(vec![RqVector::from(vec![poly1, poly2])]);
        let vec_1: RqVector = RqVector::from(vec![poly3, poly4]);

        let result_1 = matrix_1.mul(&vec_1);
        let expected_poly_1 =
            vec![Zq::new(8), Zq::new(u32::MAX - 13), Zq::new(u32::MAX - 19)].into();
        let expected_1 = RqVector::from(vec![expected_poly_1]);
        assert_eq!(result_1, expected_1);

        let poly5: Rq = vec![Zq::new(u32::MAX - 6), Zq::new(7)].into();
        let poly6: Rq = vec![Zq::new(u32::MAX - 2), Zq::ZERO].into();
        let poly7: Rq = vec![Zq::new(8), Zq::new(u32::MAX - 1)].into();
        let poly8: Rq = vec![Zq::new(u32::MAX - 3), Zq::new(4)].into();
        let poly9: Rq = vec![Zq::MAX, Zq::new(u32::MAX - 1)].into();
        let poly10: Rq = vec![Zq::new(u32::MAX - 2), Zq::new(u32::MAX - 2)].into();
        let matrix_2: RqMatrix = RqMatrix::new(vec![
            RqVector::from(vec![poly5, poly6]),
            RqVector::from(vec![poly7, poly8]),
        ]);
        let vec_2: RqVector = RqVector::from(vec![poly9, poly10]);

        let result_2 = matrix_2.mul(&vec_2);
        let expected_poly_2_1 = vec![Zq::new(16), Zq::new(16), Zq::new(u32::MAX - 13)].into();
        let expected_poly_2_2 =
            vec![Zq::new(4), Zq::new(u32::MAX - 13), Zq::new(u32::MAX - 7)].into();
        let expected_2 = RqVector::from(vec![expected_poly_2_1, expected_poly_2_2]);
        assert_eq!(result_2, expected_2);
    }
}
