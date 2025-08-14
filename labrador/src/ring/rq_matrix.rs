//! Dense matrix of `Rq` polynomials with an optional *lower‑triangular (symmetric)*
//! storage optimisation.  A symmetric matrix stores **only** the lower triangle
//! (including the diagonal).

use crate::{core::inner_product, ring::rq_vector::RqVector};
use rand::{CryptoRng, Rng};
use std::ops::Mul;

use super::{rq::Rq, zq::Zq};

/// Matrix of polynomials in Rq
#[derive(Debug, Clone)]
pub struct RqMatrix {
    elements: Vec<RqVector>,
    /// If `true`, `rows[i]` has length `i+1` and represents the lower triangle.
    is_symmetric: bool,
}

impl RqMatrix {
    /// Constructor for the Matrix of polynomials in Rq
    pub fn new(elements: Vec<RqVector>, is_symmetric: bool) -> Self {
        if is_symmetric {
            debug_assert!(
                elements.iter().enumerate().all(|(i, r)| r.len() == i + 1),
                "symmetric storage must be lower-triangular"
            );
        } else {
            let width = elements.first().map(|r| r.len()).unwrap_or(0);
            debug_assert!(
                elements.iter().all(|r| r.len() == width),
                "row lengths differ"
            );
        }
        Self {
            elements,
            is_symmetric,
        }
    }

    pub fn zero(row_len: usize, col_len: usize) -> Self {
        RqMatrix::new(vec![RqVector::zero(col_len); row_len], false)
    }

    pub fn symmetric_zero(size: usize) -> Self {
        Self {
            elements: (0..size).map(|row| RqVector::zero(row + 1)).collect(),
            is_symmetric: true,
        }
    }

    pub fn row_len(&self) -> usize {
        self.elements.len()
    }

    pub fn col_len(&self) -> usize {
        let last_row = self.row_len() - 1;
        self.elements[last_row].len()
    }

    pub fn get_cell(&self, row: usize, col: usize) -> &Rq {
        if !self.is_symmetric || row >= col {
            &self.elements[row].elements()[col]
        } else {
            &self.elements[col].elements()[row]
        }
    }

    pub fn set_cell(&mut self, row: usize, col: usize, value: Rq) {
        self.elements[row].set(col, value);
    }

    /// Create a random matrix of polynomials
    pub fn random<R: Rng + CryptoRng>(rng: &mut R, row_len: usize, col_len: usize) -> Self {
        Self {
            elements: (0..row_len)
                .map(|_| RqVector::random(rng, col_len))
                .collect(),
            is_symmetric: false,
        }
    }

    /// Create a random symmetric matrix of polynomials
    pub fn symmetric_random<R: Rng + CryptoRng>(rng: &mut R, row_len: usize) -> Self {
        Self {
            elements: (0..row_len)
                .map(|row| RqVector::random(rng, row + 1))
                .collect(),
            is_symmetric: true,
        }
    }

    pub fn elements(&self) -> &[RqVector] {
        &self.elements
    }

    /// Decompose every cell into `parts` low‑norm digits in base `base` and
    /// concat all results into **one** long vector.
    pub fn decompose_each_cell(&self, base: Zq, num_parts: usize) -> RqVector {
        let mut decomposed_vec = Vec::new();
        for ring_vector in self.elements() {
            for ring in ring_vector.elements() {
                decomposed_vec.extend(ring.decompose(
                    base,
                    u64::try_from(num_parts).expect("num_parts does not fit in u64"),
                ));
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
        RqMatrix::new(elements, false)
    }
}

// Implement matrix-vector multiplication for reference to matrix
impl Mul<&RqVector> for &RqMatrix {
    type Output = RqVector;

    fn mul(self, rhs: &RqVector) -> Self::Output {
        debug_assert_eq!(self.col_len(), rhs.len(), "dimension mismatch");
        let mut result = RqVector::zero(self.elements.len());

        for (i, row) in self.elements.iter().enumerate() {
            result.set(
                i,
                inner_product::compute_linear_combination(row.elements(), rhs.elements()),
            );
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use rand::rng;

    use super::*;
    use crate::ring::rq::tests::generate_rq_from_zq_vector;
    use crate::ring::rq::Rq;
    use crate::ring::zq::Zq;

    #[test]
    #[cfg(not(feature = "skip-slow-tests"))]
    fn rqmatrix_fits_stack() {
        let mut rng = rand::rng();
        let _: RqMatrix = RqMatrix::random(&mut rng, 256, 1 << 10);
    }

    #[test]
    fn test_set_sell() {
        let mut matrix = RqMatrix::zero(10, 18);
        matrix.set_cell(4, 9, Rq::new([Zq::new(10); Rq::DEGREE]));
        matrix.set_cell(8, 1, Rq::new([Zq::new(3); Rq::DEGREE]));

        for (i, vector) in matrix.elements().iter().enumerate() {
            for (j, poly) in vector.elements().iter().enumerate() {
                if (i == 4) && (j == 9) {
                    assert_eq!(poly, &Rq::new([Zq::new(10); Rq::DEGREE]))
                } else if (i == 8) && (j == 1) {
                    assert_eq!(poly, &Rq::new([Zq::new(3); Rq::DEGREE]))
                } else {
                    assert_eq!(poly, &Rq::zero())
                }
            }
        }
    }

    #[test]
    fn test_symmetric_matrix() {
        let symmetric_matrix = RqMatrix::symmetric_random(&mut rng(), 12);
        assert_eq!(symmetric_matrix.row_len(), 12);
        assert_eq!(symmetric_matrix.col_len(), 12);
        for i in 0..symmetric_matrix.row_len() {
            assert_eq!(symmetric_matrix.elements()[i].len(), i + 1);
        }
        for i in 0..symmetric_matrix.row_len() {
            for j in 0..symmetric_matrix.col_len() {
                assert_eq!(
                    symmetric_matrix.get_cell(i, j),
                    symmetric_matrix.get_cell(j, i)
                )
            }
        }
    }

    #[test]
    fn test_rq_matrix_from_iterator() {
        let expected = vec![
            RqVector::random(&mut rng(), 5),
            RqVector::random(&mut rng(), 5),
            RqVector::random(&mut rng(), 5),
            RqVector::random(&mut rng(), 5),
        ];
        let polynomial_matrix = expected.clone().into_iter();
        let result: RqMatrix = polynomial_matrix.collect();

        assert_eq!(result.elements(), &expected);
    }

    #[test]
    fn test_zero_matrix() {
        /// Check if Polynomial == 0
        pub fn is_polynomial_zero(poly: &Rq) -> bool {
            poly.coeffs().iter().all(|&coeff| coeff == Zq::ZERO)
        }

        let matrix = RqMatrix::zero(10, 20);
        assert_eq!(matrix.row_len(), 10);
        assert_eq!(matrix.col_len(), 20);
        for row in matrix.elements() {
            for cell in row.elements() {
                assert!(is_polynomial_zero(cell));
            }
        }
    }

    #[test]
    fn test_rqmartrix_mul() {
        let poly1: Rq = generate_rq_from_zq_vector(vec![Zq::new(8), Zq::new(6)]);
        let poly2: Rq = generate_rq_from_zq_vector(vec![-Zq::new(5), -Zq::new(5)]);
        let poly3: Rq = generate_rq_from_zq_vector(vec![Zq::ONE, Zq::ZERO]);
        let poly4: Rq = generate_rq_from_zq_vector(vec![Zq::ZERO, Zq::new(4)]);
        let matrix_1: RqMatrix = RqMatrix::new(vec![RqVector::from(vec![poly1, poly2])], false);
        let vec_1: RqVector = RqVector::from(vec![poly3, poly4]);

        let result_1 = matrix_1.mul(&vec_1);
        let expected_poly_1 =
            generate_rq_from_zq_vector(vec![Zq::new(8), -Zq::new(14), -Zq::new(20)]);
        let expected_1 = RqVector::from(vec![expected_poly_1]);
        assert_eq!(result_1, expected_1);

        let poly5: Rq = generate_rq_from_zq_vector(vec![-Zq::new(7), Zq::new(7)]);
        let poly6: Rq = generate_rq_from_zq_vector(vec![-Zq::new(3), Zq::ZERO]);
        let poly7: Rq = generate_rq_from_zq_vector(vec![Zq::new(8), -Zq::new(2)]);
        let poly8: Rq = generate_rq_from_zq_vector(vec![-Zq::new(4), Zq::new(4)]);
        let poly9: Rq = generate_rq_from_zq_vector(vec![Zq::NEG_ONE, -Zq::new(2)]);
        let poly10: Rq = generate_rq_from_zq_vector(vec![-Zq::new(3), -Zq::new(3)]);
        let matrix_2: RqMatrix = RqMatrix::new(
            vec![
                RqVector::from(vec![poly5, poly6]),
                RqVector::from(vec![poly7, poly8]),
            ],
            false,
        );
        let vec_2: RqVector = RqVector::from(vec![poly9, poly10]);

        let result_2 = matrix_2.mul(&vec_2);
        let expected_poly_2_1 =
            generate_rq_from_zq_vector(vec![Zq::new(16), Zq::new(16), -Zq::new(14)]);
        let expected_poly_2_2 =
            generate_rq_from_zq_vector(vec![Zq::new(4), -Zq::new(14), -Zq::new(8)]);
        let expected_2 = RqVector::from(vec![expected_poly_2_1, expected_poly_2_2]);
        assert_eq!(result_2, expected_2);
    }
}
