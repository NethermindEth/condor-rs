use std::ops::{Add, Mul, Sub};

use crate::ring::rq::Rq;
use crate::ring::rq_vector::RqVector;
use crate::ring::zq::Zq;
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

    /// Generate random small polynomial with secure RNG implementation
    pub fn random_ternary<R: Rng + CryptoRng>(rng: &mut R, n: usize) -> Self {
        let mut coeffs = vec![Zq::ZERO; n];

        for coeff in coeffs.iter_mut() {
            // Explicitly sample from {-1, 0, 1} with equal probability
            let val = match rng.random_range(0..3) {
                0 => Zq::TWO,  // 2
                1 => Zq::ZERO, // 0
                2 => Zq::ONE,  // 1
                _ => unreachable!(),
            };
            *coeff = val;
        }

        Self::new(coeffs)
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
    ///
    /// Note that: The operator norm only affects the coefficients of the random PolyRings generated from the challenge space.
    /// Prover and Verifier will not do the operator norm check, because random PolyRings are determined after generation.
    /// Both party will have access to the same PolyRings through transcript,
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

    /// Decomposes a polynomial into base-B representation:
    /// p = p⁽⁰⁾ + p⁽¹⁾·B + p⁽²⁾·B² + ... + p⁽ᵗ⁻¹⁾·B^(t-1)
    /// Where each p⁽ⁱ⁾ has small coefficients, using centered representatives
    pub fn decompose(&self, base: Zq, num_parts: usize) -> PolyVector {
        let mut parts = Vec::with_capacity(num_parts);
        let mut current = self.clone();

        for i in 0..num_parts {
            if i == num_parts - 1 {
                parts.push(current.clone());
            } else {
                // Extract low part (mod base, centered around 0)
                let mut low_coeffs = vec![Zq::ZERO; self.len()];

                for (j, coeff) in current.get_coeffs().iter().enumerate() {
                    low_coeffs[j] = coeff.centered_mod(base);
                }

                let low_part = Self::new(low_coeffs);
                parts.push(low_part.clone());

                // Update current
                current = &current - &low_part;

                // Scale by base
                let mut scaled_coeffs = vec![Zq::ZERO; self.len()];
                for (j, coeff) in current.get_coeffs().iter().enumerate() {
                    scaled_coeffs[j] = coeff.scale_by(base);
                }
                current = Self::new(scaled_coeffs);
            }
        }

        PolyVector::new(parts)
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

impl Sub<&PolyRing> for &PolyRing {
    type Output = PolyRing;
    /// Sub two polynomials with flexible degree
    fn sub(self, other: &PolyRing) -> PolyRing {
        let max_degree = self.get_coeffs().len().max(other.get_coeffs().len());
        let mut coeffs = vec![Zq::ZERO; max_degree];
        for (i, coeff) in coeffs.iter_mut().enumerate().take(max_degree) {
            if i < self.get_coeffs().len() {
                *coeff += self.get_coeffs()[i];
            }
            if i < other.get_coeffs().len() {
                *coeff -= other.get_coeffs()[i];
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
    pub fn random(n: usize, m: usize) -> Self {
        let mut vector = PolyVector::new(vec![]);
        vector.elements = (0..n)
            .map(|_| PolyRing::random(&mut rand::rng(), m))
            .collect();
        vector
    }

    /// Generate random small polynomial with secure RNG implementation
    pub fn random_ternary(n: usize, m: usize) -> Self {
        let mut vector = PolyVector::new(vec![]);
        vector.elements = (0..n)
            .map(|_| PolyRing::random_ternary(&mut rand::rng(), m))
            .collect();
        vector
    }

    pub fn inner_product_poly_vector(&self, other: &PolyVector) -> PolyRing {
        self.iter().zip(other.iter()).map(|(a, b)| a * b).fold(
            PolyRing::zero(self.get_elements()[0].get_coeffs().len()),
            |acc, val| &acc + &val,
        )
    }

    // Compute the squared norm of a vector of polynomials
    pub fn compute_norm_squared(&self) -> Zq {
        self.elements
            .iter()
            .flat_map(|poly| poly.get_coeffs()) // Collect coefficients from all polynomials
            .map(|coeff| *coeff * *coeff)
            .sum()
    }

    /// Function to concatenate coefficients from multiple Rq into a Vec<Zq>
    pub fn concatenate_coefficients(&self, s: usize) -> ZqVector {
        let total_coeffs = self.get_elements().len() * s;
        let mut concatenated_coeffs: Vec<Zq> = Vec::with_capacity(total_coeffs);
        // Iterate over each Rq, extracting the coefficients and concatenating them
        for rq in self.get_elements() {
            let coeffs = rq.get_coeffs();
            concatenated_coeffs.extend_from_slice(coeffs);
        }

        ZqVector {
            coeffs: concatenated_coeffs,
        }
    }

    pub fn decompose(&self, b: Zq, parts: usize) -> Vec<PolyVector> {
        self.iter()
            .map(|i| PolyRing::decompose(i, b, parts))
            .collect()
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

impl Mul<&PolyRing> for &PolyVector {
    type Output = PolyVector;
    // A poly vector multiple by a PolyRing
    fn mul(self, other: &PolyRing) -> PolyVector {
        self.iter().map(|s| s * other).collect()
    }
}

impl Mul<&Vec<PolyVector>> for &PolyVector {
    type Output = PolyVector;
    // a poly vector mnultiple by a vec<PolyVector>
    fn mul(self, other: &Vec<PolyVector>) -> PolyVector {
        other
            .iter()
            .map(|o| o.inner_product_poly_vector(self))
            .collect()
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

    pub fn zero() -> Self {
        Self::new(vec![])
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

    /// Generate random small polynomial with secure RNG implementation
    pub fn random_ternary<R: Rng + CryptoRng>(rng: &mut R, n: usize) -> Self {
        let mut coeffs = vec![Zq::ZERO; n];

        for coeff in coeffs.iter_mut() {
            // Explicitly sample from {-1, 0, 1} with equal probability
            let val = match rng.random_range(0..3) {
                0 => Zq::TWO,  // 2
                1 => Zq::ZERO, // 0
                2 => Zq::ONE,  // 1
                _ => unreachable!(),
            };
            *coeff = val;
        }

        Self::new(coeffs)
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
