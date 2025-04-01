use std::ops::{Add, Mul};

use crate::{rq::Rq, rq_vector::RqVector, zq::Zq};
use rand::distr::{Distribution, Uniform};
use rand::{CryptoRng, Rng};
use rustfft::{num_complex::Complex, FftPlanner};

/// A PolyRing is a vector of Zq elements with a flexible degree that is less than or equal to DEGREE_BOUND.
#[derive(Debug, Clone, PartialEq, PartialOrd, Eq, Ord, Default)]
pub struct PolyRing {
    coeffs: Vec<Zq>,
}
impl PolyRing {
    // degree bound of the polynomial ring, set 64 from the paper
    pub const DEGREE_BOUND: usize = 64;

    pub fn new(coeffs: Vec<Zq>) -> Self {
        assert!(
            coeffs.len() <= Self::DEGREE_BOUND,
            "Polynomial degree should be less than {}",
            Self::DEGREE_BOUND
        );
        Self { coeffs }
    }
    pub fn zero(degree: usize) -> Self {
        Self::new(vec![Zq::ZERO; degree])
    }

    pub fn zero_poly() -> Self {
        Self::new(vec![Zq::ZERO; 1])
    }

    pub fn len(&self) -> usize {
        self.coeffs.len()
    }

    pub fn is_empty(&self) -> bool {
        self.coeffs.is_empty()
    }

    pub fn get_coeffs(&self) -> &Vec<Zq> {
        &self.coeffs
    }

    pub fn iter(&self) -> impl Iterator<Item = &Zq> {
        self.coeffs.iter()
    }

    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut Zq> {
        self.coeffs.iter_mut()
    }

    /// inner product of two polynomials
    pub fn inner_product(&self, other: &Self) -> Zq {
        self.coeffs
            .iter()
            .zip(other.coeffs.iter())
            .map(|(a, b)| *a * *b)
            .sum()
    }

    /// Generate random Zq vector with a provided cryptographically secure RNG
    pub fn random<R: Rng + CryptoRng>(rng: &mut R, n: usize) -> Self {
        let uniform = Uniform::new_inclusive(Zq::ZERO, Zq::MAX).unwrap();
        let mut coeffs = Vec::with_capacity(n);
        coeffs.extend((0..n).map(|_| uniform.sample(rng)));
        Self { coeffs }
    }

    /// Compute the conjugate automorphism \sigma_{-1} of vector based on B) Constraints..., Page 21.
    pub fn conjugate_automorphism(&self) -> PolyRing {
        let q_minus_1 = Zq::MAX;
        let mut new_coeffs = vec![Zq::ZERO; PolyRing::DEGREE_BOUND];
        for (i, new_coeff) in new_coeffs
            .iter_mut()
            .enumerate()
            .take(PolyRing::DEGREE_BOUND)
        {
            if i < self.get_coeffs().len() {
                if i == 0 {
                    *new_coeff = self.get_coeffs()[i];
                } else {
                    *new_coeff = self.get_coeffs()[i] * q_minus_1;
                }
            } else {
                *new_coeff = Zq::ZERO;
            }
        }
        let reversed_coefficients = new_coeffs
            .iter()
            .take(1)
            .cloned()
            .chain(new_coeffs.iter().skip(1).rev().cloned())
            .collect::<Vec<Zq>>();

        PolyRing::new(reversed_coefficients)
    }

    /// Compute the operator norm of a polynomial given its coefficients.
    /// The operator norm is defined as the maximum magnitude of the DFT (eigenvalues)
    /// of the coefficient vector.
    #[allow(clippy::as_conversions)]
    pub fn operator_norm(&self) -> f64 {
        let coeffs = self.get_coeffs();
        let n = coeffs.len();
        let mut planner = FftPlanner::new();
        let fft = planner.plan_fft_forward(n);

        // Convert coefficients into complex numbers (with zero imaginary parts)
        let mut buffer: Vec<Complex<f64>> = coeffs
            .iter()
            .map(|&x| {
                let half = Zq::MAX.scale_by(Zq::TWO);
                let converted_value = if x > half {
                    x.to_u128() as f64 - Zq::MAX.to_u128() as f64 - 1.0
                } else {
                    x.to_u128() as f64
                };
                Complex {
                    re: converted_value,
                    im: 0.0,
                }
            })
            .collect();

        // Compute the FFT (this gives the eigenvalues of the circulant matrix)
        fft.process(&mut buffer);

        // Return the maximum absolute value (norm) among the eigenvalues
        buffer
            .iter()
            .map(|c| c.norm())
            .fold(0.0, |max, x| max.max(x))
    }
}

impl<const D: usize> From<PolyRing> for Rq<D> {
    fn from(zqs: PolyRing) -> Self {
        zqs.get_coeffs().clone().into()
    }
}

impl FromIterator<Zq> for PolyRing {
    fn from_iter<T: IntoIterator<Item = Zq>>(iter: T) -> Self {
        let coeffs: Vec<Zq> = iter.into_iter().collect();
        PolyRing::new(coeffs)
    }
}

impl Add<&PolyRing> for &PolyRing {
    type Output = PolyRing;
    /// Add two polynomials with flexible degree
    fn add(self, other: &PolyRing) -> PolyRing {
        let max_degree = self.get_coeffs().len().max(other.get_coeffs().len());
        let mut coeffs = vec![Zq::ZERO; max_degree];
        for (i, coeff) in coeffs.iter_mut().enumerate().take(max_degree) {
            if i < self.get_coeffs().len() {
                *coeff += self.get_coeffs()[i];
            }
            if i < other.get_coeffs().len() {
                *coeff += other.get_coeffs()[i];
            }
        }
        PolyRing::new(coeffs)
    }
}

impl Mul<&PolyRing> for &PolyRing {
    type Output = PolyRing;
    /// Polynomial multiplication of two polynomials
    fn mul(self, other: &PolyRing) -> PolyRing {
        // Initialize a vector to hold the intermediate multiplication result
        let mut result_coefficients =
            vec![Zq::new(0); self.get_coeffs().len() + other.get_coeffs().len() - 1];
        for (i, &coeff1) in self.get_coeffs().iter().enumerate() {
            for (j, &coeff2) in other.get_coeffs().iter().enumerate() {
                result_coefficients[i + j] += coeff1 * coeff2;
            }
        }

        // Reduce modulo X^d + 1
        if result_coefficients.len() > PolyRing::DEGREE_BOUND {
            let q_minus_1 = Zq::MAX;
            let (left, right) = result_coefficients.split_at_mut(PolyRing::DEGREE_BOUND);
            for (i, &overflow) in right.iter().enumerate() {
                left[i] += overflow * q_minus_1;
            }
            result_coefficients.truncate(PolyRing::DEGREE_BOUND);
        }
        PolyRing::new(result_coefficients)
    }
}

impl Mul<&Zq> for &PolyRing {
    type Output = PolyRing;
    /// Scalar multiplication of a polynomial
    fn mul(self, other: &Zq) -> PolyRing {
        PolyRing::new(self.coeffs.iter().map(|c| c * *other).collect())
    }
}

#[derive(Debug, Clone, PartialEq, PartialOrd, Eq, Ord, Default)]
pub struct PolyVector {
    elements: Vec<PolyRing>,
}

impl PolyVector {
    pub fn new(elements: Vec<PolyRing>) -> Self {
        Self { elements }
    }

    pub fn zero() -> Self {
        Self {
            elements: vec![PolyRing::zero(0)],
        }
    }

    pub fn get_elements(&self) -> &Vec<PolyRing> {
        &self.elements
    }

    pub fn len(&self) -> usize {
        self.elements.len()
    }

    pub fn is_empty(&self) -> bool {
        self.elements.is_empty()
    }

    pub fn iter(&self) -> impl Iterator<Item = &PolyRing> {
        self.elements.iter()
    }

    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut PolyRing> {
        self.elements.iter_mut()
    }

    // Generate a random polynomial vector with n polynomials of degree m
    pub fn random(n: usize, m: usize) -> PolyVector {
        let mut vector = PolyVector::new(vec![]);
        vector.elements = (0..n)
            .map(|_| PolyRing::random(&mut rand::rng(), m))
            .collect();
        vector
    }

    pub fn inner_product_poly_vector(&self, other: &PolyVector) -> PolyRing {
        self.iter().zip(other.iter()).map(|(a, b)| a * b).fold(
            PolyRing::zero(self.get_elements()[0].get_coeffs().len()),
            |acc, val| &acc + &val,
        )
    }
}

impl<const N: usize, const D: usize> From<PolyVector> for RqVector<N, D> {
    fn from(polys: PolyVector) -> Self {
        let mut rq_vector = RqVector::zero();
        for (i, poly) in polys.elements.iter().enumerate() {
            rq_vector[i] = poly.get_coeffs().clone().into();
        }
        rq_vector
    }
}

impl FromIterator<PolyRing> for PolyVector {
    fn from_iter<T: IntoIterator<Item = PolyRing>>(iter: T) -> Self {
        let mut elements = Vec::new();
        for item in iter {
            elements.push(item);
        }
        PolyVector::new(elements)
    }
}

impl Add<&PolyVector> for &PolyVector {
    type Output = PolyVector;
    // add two poly vectors
    fn add(self, other: &PolyVector) -> PolyVector {
        self.iter().zip(other.iter()).map(|(a, b)| a + b).collect()
    }
}

impl Mul<&Zq> for &PolyVector {
    type Output = PolyVector;
    // A poly vector multiply by a Zq
    fn mul(self, other: &Zq) -> PolyVector {
        self.iter().map(|a| a * other).collect()
    }
}

/// A ZqVector is a vector of Zq elements with a flexible size.
/// Mainly used for store random Zq elements
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ZqVector {
    coeffs: Vec<Zq>,
}
impl ZqVector {
    pub fn new(coeffs: Vec<Zq>) -> Self {
        Self { coeffs }
    }

    pub fn get_coeffs(&self) -> &Vec<Zq> {
        &self.coeffs
    }

    pub fn len(&self) -> usize {
        self.coeffs.len()
    }

    pub fn is_empty(&self) -> bool {
        self.coeffs.is_empty()
    }

    pub fn random<R: Rng + CryptoRng>(rng: &mut R, n: usize) -> Self {
        let uniform = Uniform::new_inclusive(Zq::ZERO, Zq::MAX).unwrap();
        let mut coeffs = Vec::with_capacity(n);
        coeffs.extend((0..n).map(|_| uniform.sample(rng)));
        Self { coeffs }
    }

    pub fn iter(&self) -> impl Iterator<Item = &Zq> {
        self.coeffs.iter()
    }

    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut Zq> {
        self.coeffs.iter_mut()
    }
}

impl<const D: usize> From<ZqVector> for Rq<D> {
    fn from(zqs: ZqVector) -> Self {
        zqs.get_coeffs().clone().into()
    }
}

impl FromIterator<Zq> for ZqVector {
    fn from_iter<T: IntoIterator<Item = Zq>>(iter: T) -> Self {
        let coeffs: Vec<Zq> = iter.into_iter().collect();
        ZqVector::new(coeffs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_conjugate_automorphism() {
        let poly1: PolyRing = PolyRing::new(vec![Zq::ONE, Zq::TWO, Zq::new(3)]);
        let poly2: PolyRing = PolyRing::new(vec![Zq::new(4), Zq::new(5), Zq::new(6)]);
        let inner_12 = poly1.inner_product(&poly2);
        let conjugated_1 = poly1.conjugate_automorphism();
        let inner_conjugated_12 = &conjugated_1 * &poly2;

        assert_eq!(inner_conjugated_12.len(), PolyRing::DEGREE_BOUND);
        assert_eq!(inner_conjugated_12.get_coeffs()[0], Zq::from(32));
        assert_eq!(inner_conjugated_12.get_coeffs()[1], Zq::from(17));
        assert_eq!(inner_conjugated_12.get_coeffs()[2], Zq::new(6));

        // ct<\sigma_{-1}(poly1), poly2> ?= <poly1, poly2>
        let ct_inner_conjugated_12 = inner_conjugated_12.get_coeffs()[0];
        assert_eq!(ct_inner_conjugated_12, inner_12);
    }

    #[test]
    fn test_polyring_to_rq() {
        let polyring = PolyRing::new(vec![Zq::ONE, Zq::TWO, Zq::new(3)]);
        let rq_vector: Rq<3> = polyring.into();
        let expect_rq = Rq::new([Zq::ONE, Zq::TWO, Zq::new(3)]);
        assert_eq!(rq_vector, expect_rq);
    }

    #[test]
    fn test_zqvector_to_rq() {
        let zq_vector = ZqVector::new(vec![Zq::ONE, Zq::TWO, Zq::new(3)]);
        let rq_vector: Rq<3> = zq_vector.into();
        let expect_rq = Rq::new([Zq::ONE, Zq::TWO, Zq::new(3)]);
        assert_eq!(rq_vector, expect_rq);
    }

    #[test]
    fn test_polyvector_to_rqvector() {
        let poly_vector = PolyVector::new(vec![PolyRing::new(vec![Zq::ONE, Zq::TWO, Zq::new(3)])]);
        let rqs: RqVector<1, 3> = poly_vector.into();
        let expect_rq = Rq::new([Zq::ONE, Zq::TWO, Zq::new(3)]);
        let expect_rqvector: RqVector<1, 3> = RqVector::from(vec![expect_rq]);
        assert_eq!(rqs, expect_rqvector);
    }

    #[test]
    fn test_scalar_mul_vector_mul_vector() {
        let vector = PolyVector::new(vec![PolyRing::new(vec![Zq::ONE, Zq::TWO, Zq::new(3)])]);
        let zq = Zq::new(2);
        let result = &vector * &zq;
        let expect = PolyVector::new(vec![PolyRing::new(vec![Zq::TWO, Zq::new(4), Zq::new(6)])]);
        assert_eq!(result, expect)
    }

    #[test]
    fn test_add_poly_vector() {
        let vector1 = PolyVector::new(vec![PolyRing::new(vec![Zq::ONE, Zq::TWO, Zq::new(3)])]);
        let vector2 = PolyVector::new(vec![PolyRing::new(vec![Zq::new(4), Zq::new(5)])]);
        let result = &vector1 + &vector2;
        let expect = PolyVector::new(vec![PolyRing::new(vec![
            Zq::new(5),
            Zq::new(7),
            Zq::new(3),
        ])]);
        assert_eq!(result, expect)
    }
}
