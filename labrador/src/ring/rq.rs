// This file is part of the polynomial ring operations module.
//
//
// Currently implemented functions include:
// - Polynomial addition:          +
// - Polynomial multiplication:    *
// - inner_product/ Dot product:   inner_product()
// - Polynomial subtraction:       -
// - Polynomial negation:          neg()
// - Scalar multiplication:        scalar_mul()
// - Polynomial evaluation:        eval()
// - Zero check:                   is_zero()
// - Polynomial equality check:    is_equal()
// - Get the Coefficients:         get_coefficients()
// - Random small norm vector:     random_small_vector()
// - Squared norm of coefficients: compute_norm_squared()
//
// Further operations and optimizations will be added in future versions.

// We use the Zq ring
use crate::ring::zq::Zq;
use core::ops::{Add, Mul, Neg, Sub};
use rand::distr::{Distribution, Uniform};
use rand::{CryptoRng, Rng};
use rustfft::num_complex::Complex;
use rustfft::FftPlanner;
use std::iter::Sum;

use super::rq_vector::RqVector;

/// This module provides implementations for various operations
/// in the polynomial ring R = Z_q\[X\] / (X^d + 1).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Rq {
    coeffs: [Zq; Self::DEGREE],
}

impl Rq {
    pub const DEGREE: usize = 64;
    /// Constructor for the polynomial ring
    pub const fn new(coeffs: [Zq; Self::DEGREE]) -> Self {
        Rq { coeffs }
    }

    /// Generate zero polynomial
    pub const fn zero() -> Self {
        Self {
            coeffs: [Zq::ZERO; Self::DEGREE],
        }
    }

    /// Get the coefficients as a vector
    pub fn get_coefficients(&self) -> &[Zq; Self::DEGREE] {
        &self.coeffs
    }

    /// Dot product between coefficients
    pub fn inner_product(&self, other: &Self) -> Zq {
        self.coeffs
            .iter()
            .zip(other.coeffs.iter())
            .map(|(&a, &b)| a * b)
            .fold(Zq::ZERO, |acc, x| acc + x)
    }

    /// Evaluate the polynomial at a specific point
    pub fn eval(&self, x: Zq) -> Zq {
        let mut result = Zq::ZERO;
        for coeff in self.coeffs.iter().rev() {
            result = result * x + *coeff;
        }

        result
    }

    pub fn add_coefficients(&self) -> Zq {
        self.get_coefficients()
            .iter()
            .fold(Zq::ZERO, |acc, coeff| acc + *coeff)
    }

    /// Check if Polynomial == 0
    pub fn is_zero(&self) -> bool {
        self.coeffs.iter().all(|&coeff| coeff == Zq::ZERO)
    }

    /// Check if two polynomials are equal
    pub fn is_equal(&self, other: &Self) -> bool {
        self.coeffs == other.coeffs
    }

    /// Generate random polynomial with a provided cryptographically secure RNG
    pub fn random<R: Rng + CryptoRng>(rng: &mut R) -> Self {
        let uniform = Uniform::new_inclusive(Zq::ZERO, Zq::MAX).unwrap();
        let mut coeffs = [Zq::ZERO; Self::DEGREE];
        coeffs.iter_mut().for_each(|c| *c = uniform.sample(rng));
        Self { coeffs }
    }

    /// Generate random small polynomial with secure RNG implementation
    pub fn random_ternary<R: Rng + CryptoRng>(rng: &mut R) -> Self {
        let mut coeffs = [Zq::ZERO; Self::DEGREE];

        for coeff in coeffs.iter_mut() {
            // Explicitly sample from {-1, 0, 1} with equal probability
            let val = match rng.random_range(0..3) {
                0 => Zq::MAX,  // -1 mod q
                1 => Zq::ZERO, // 0
                2 => Zq::ONE,  // 1
                _ => unreachable!(),
            };
            *coeff = val;
        }

        Rq::new(coeffs)
    }

    /// Decomposes a polynomial into base-B representation:
    /// p = p⁽⁰⁾ + p⁽¹⁾·B + p⁽²⁾·B² + ... + p⁽ᵗ⁻¹⁾·B^(t-1)
    /// Where each p⁽ⁱ⁾ has small coefficients, using centered representatives
    pub fn decompose(&self, base: Zq, num_parts: usize) -> RqVector {
        debug_assert!(num_parts > 0, "num_parts must be positive");
        let mut parts = Vec::with_capacity(num_parts);
        let mut initial_coeffs = *self.get_coefficients();

        for _ in 0..num_parts - 1 {
            // Extract low part (mod base, centered around 0)
            let mut low_coeffs = [Zq::ZERO; Self::DEGREE];

            for (low_c, coeff) in low_coeffs.iter_mut().zip(initial_coeffs.iter_mut()) {
                let centered = coeff.centered_mod(base);
                *low_c = centered;
                *coeff = (*coeff - centered).scale_by(base);
            }
            parts.push(Self::new(low_coeffs));
        }
        parts.push(Self::new(initial_coeffs));
        RqVector::new(parts)
    }

    /// Encode message into polynomial with small coefficients.
    ///
    /// # Arguments
    /// * `message` - A slice of booleans representing a binary message
    ///
    /// # Returns
    /// * `Some(Rq)` - A polynomial where each coefficient is 0 or 1 based on the message bits
    /// * `None` - If the message length exceeds the polynomial degree D
    ///
    /// # Format
    /// * Each boolean is encoded as a coefficient: false -> 0, true -> 1
    /// * Message bits are mapped to coefficients in order (index 0 -> constant term)
    /// * Remaining coefficients (if message is shorter than D) are set to 0
    pub fn encode_message(message: &[bool]) -> Option<Self> {
        if message.len() > Self::DEGREE {
            return None;
        }

        let mut coeffs = [Zq::ZERO; Self::DEGREE];
        for (i, &bit) in message.iter().enumerate() {
            coeffs[i] = Zq::new(u32::from(bit));
        }
        Some(Rq::new(coeffs))
    }

    /// Iterator over coefficients
    pub fn iter(&self) -> std::slice::Iter<'_, Zq> {
        self.coeffs.iter()
    }

    /// Check if polynomial coefficients are within bounds
    pub fn check_bounds(&self, bound: Zq) -> bool {
        self.iter().all(|coeff| coeff <= &bound || coeff >= &-bound)
    }

    /// Compute the conjugate automorphism \sigma_{-1} of vector based on B) Constraints..., Page 21.
    pub fn conjugate_automorphism(&self) -> Self {
        let q_minus_1 = Zq::MAX;
        let mut new_coeffs = [Zq::ZERO; Self::DEGREE];
        for (i, new_coeff) in new_coeffs.iter_mut().enumerate().take(Self::DEGREE) {
            if i < self.get_coefficients().len() {
                if i == 0 {
                    *new_coeff = self.get_coefficients()[i];
                } else {
                    *new_coeff = self.get_coefficients()[i] * q_minus_1;
                }
            } else {
                *new_coeff = Zq::ZERO;
            }
        }
        debug_assert_eq!(new_coeffs.len(), Self::DEGREE);
        let mut reversed_coefficients = new_coeffs; // copy / move the array
        reversed_coefficients[1..].reverse(); // reverse everything except index 0
        Self::new(reversed_coefficients)
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
        let coeffs = self.get_coefficients();
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

impl From<Vec<Zq>> for Rq {
    fn from(vec: Vec<Zq>) -> Self {
        let mut temp = [Zq::ZERO; Self::DEGREE];
        // Process excess terms with sign adjustment
        for i in (0..vec.len()).rev() {
            let m = i / Self::DEGREE;
            let r = i % Self::DEGREE;
            let sign = if m % 2 == 0 { 1 } else { -1 };
            if sign == 1 {
                temp[r] += vec[i];
            } else {
                temp[r] -= vec[i];
            }
        }
        Rq::new(temp)
    }
}

impl Mul<Zq> for &Rq {
    type Output = Rq;
    /// Scalar multiplication of a polynomial
    fn mul(self, other: Zq) -> Rq {
        let mut copied_coeffs = self.coeffs;
        for elem in copied_coeffs.iter_mut() {
            *elem *= other;
        }
        Rq::new(copied_coeffs)
    }
}

impl Mul<&Zq> for &Rq {
    type Output = Rq;
    /// Scalar multiplication of a polynomial
    fn mul(self, other: &Zq) -> Rq {
        let mut copied_coeffs = self.coeffs;
        for elem in copied_coeffs.iter_mut() {
            *elem *= *other;
        }
        Rq::new(copied_coeffs)
    }
}

impl Add<&Rq> for &Rq {
    type Output = Rq;
    /// Add two polynomials
    fn add(self, other: &Rq) -> Rq {
        let mut coeffs = [Zq::ZERO; Rq::DEGREE];
        for (r, (a, b)) in coeffs
            .iter_mut()
            .zip(self.coeffs.iter().zip(other.coeffs.iter()))
        {
            *r = *a + *b;
        }
        Rq::new(coeffs)
    }
}

impl Sub<&Rq> for &Rq {
    type Output = Rq;
    /// Add two polynomials
    fn sub(self, other: &Rq) -> Rq {
        let mut coeffs = [Zq::ZERO; Rq::DEGREE];
        for (r, (a, b)) in coeffs
            .iter_mut()
            .zip(self.coeffs.iter().zip(other.coeffs.iter()))
        {
            *r = *a - *b;
        }
        Rq::new(coeffs)
    }
}

impl Sum for Zq {
    // Accumulate using the addition operator
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = Zq>,
    {
        iter.fold(Zq::ZERO, |acc, x| acc + x)
    }
}

impl Mul<&Rq> for &Rq {
    type Output = Rq;
    /// Polynomial multiplication modulo x^D + 1
    fn mul(self, other: &Rq) -> Rq {
        let mut result = [Zq::ZERO; Rq::DEGREE];
        let mut out_of_field = [Zq::ZERO; Rq::DEGREE];
        for (i, &self_coeff) in self.coeffs.iter().enumerate() {
            for (j, &other_coeff) in other.coeffs.iter().enumerate() {
                if i + j < Rq::DEGREE {
                    result[i + j] += self_coeff * other_coeff;
                } else {
                    out_of_field[(i + j) % Rq::DEGREE] += self_coeff * other_coeff;
                }
            }
        }
        // Process excess terms with sign adjustment
        for i in (0..Rq::DEGREE).rev() {
            let m = i / Rq::DEGREE;
            let r = i % Rq::DEGREE;
            let sign = if (m + 1) % 2 == 0 { 1 } else { -1 };
            if sign == 1 {
                result[r] += out_of_field[i];
            } else {
                result[r] -= out_of_field[i];
            }
        }
        Rq::new(result)
    }
}

// Implementing the Neg trait
impl Neg for Rq {
    type Output = Self;

    /// Polynomial negation
    fn neg(self) -> Self {
        let mut result = [Zq::ZERO; Rq::DEGREE];
        for (i, &coeff) in self.coeffs.iter().enumerate() {
            result[i] = Zq::ZERO - coeff;
        }
        Rq::new(result)
    }
}

impl FromIterator<Zq> for Rq {
    fn from_iter<T: IntoIterator<Item = Zq>>(iter: T) -> Self {
        let coeffs_vec: Vec<Zq> = iter.into_iter().collect();
        assert_eq!(
            coeffs_vec.len(),
            Self::DEGREE,
            "Iterator must contain exactly {} elements to create Rq<{}>",
            Self::DEGREE,
            Self::DEGREE
        );

        // Convert the vector to an array
        let mut coeffs = [Zq::default(); Self::DEGREE];
        for (i, coeff) in coeffs_vec.into_iter().enumerate().take(Self::DEGREE) {
            coeffs[i] = coeff;
        }

        Rq::new(coeffs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    mod helper {
        use super::*;
        pub fn padded(prefix: &[Zq]) -> [Zq; Rq::DEGREE] {
            assert!(
                prefix.len() <= Rq::DEGREE,
                "too many coefficients for degree {}",
                Rq::DEGREE
            );

            let mut out = [Zq::ZERO; Rq::DEGREE];
            out[..prefix.len()].copy_from_slice(prefix);
            out
        }

        pub fn rq_from(prefix: &[Zq]) -> Rq {
            Rq {
                coeffs: padded(prefix),
            }
        }
    }

    // Test new() and polynomial creation
    #[test]
    fn test_new_and_create_poly() {
        let coeffs = [Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)];
        let poly = helper::rq_from(&coeffs);
        assert_eq!(poly.coeffs, helper::padded(&coeffs));

        // Direct conversion
        let coeffs2 = [Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)];
        let poly_from_vec_direct: Rq = coeffs.to_vec().into();
        assert_eq!(poly_from_vec_direct.coeffs, helper::padded(&coeffs2));
    }

    // Test addition of polynomials
    #[test]
    fn test_add() {
        // Within bounds
        let poly1: Rq = vec![Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)].into();
        let poly2: Rq = vec![Zq::new(4), Zq::new(3), Zq::new(2), Zq::ONE].into();
        let result = &poly1 + &poly2;
        assert_eq!(
            result.coeffs,
            helper::padded(&[Zq::new(5), Zq::new(5), Zq::new(5), Zq::new(5)])
        );

        // Outside of bounds
        let poly3: Rq = vec![Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)].into();
        let poly4: Rq = vec![Zq::MAX, Zq::new(3), Zq::MAX, Zq::ONE].into();
        let result2 = &poly3 + &poly4;
        assert_eq!(
            result2.coeffs,
            helper::padded(&[Zq::ZERO, Zq::new(5), Zq::new(2), Zq::new(5)])
        );
        // Addition with zero polynomial
        let poly5: Rq = vec![Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)].into();
        let poly6: Rq = vec![Zq::ZERO].into();
        let result3 = &poly5 + &poly6;
        assert_eq!(
            result3.coeffs,
            helper::padded(&[Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)])
        );
        // Addition with high coefficients
        let poly7: Rq = vec![Zq::ONE, Zq::new(2), Zq::new(3), Zq::MAX].into();
        let poly8: Rq = vec![Zq::MAX, Zq::MAX, Zq::MAX, Zq::MAX].into();
        let result3 = &poly7 + &poly8;
        assert_eq!(
            result3.coeffs,
            helper::padded(&[
                Zq::ZERO,
                Zq::ONE,
                Zq::new(2),
                Zq::new(u32::MAX.wrapping_add(u32::MAX))
            ])
        );
    }

    // Test multiplication of polynomials
    #[test]
    fn test_mul() {
        // Multiplication with wrapping
        let poly1: Rq = vec![Zq::ONE, Zq::ONE, Zq::new(2)].into();
        let poly2: Rq = vec![Zq::ONE, Zq::ONE].into();
        let result = &poly1 * &poly2;
        assert_eq!(
            result.coeffs,
            helper::padded(&[Zq::new(1), Zq::new(2), Zq::new(3), Zq::new(2)])
        );

        // Multiplication with zero polynomial
        let poly3: Rq = vec![Zq::ONE, Zq::ONE, Zq::new(2)].into();
        let poly4: Rq = vec![Zq::ZERO].into();
        let result2 = &poly3 * &poly4;
        assert_eq!(
            result2.coeffs,
            helper::padded(&[Zq::ZERO, Zq::ZERO, Zq::ZERO])
        );

        // Needs to be revised later
        // // Multiplication with wrapping higher order
        // let poly5: Rq = vec![Zq::ONE, Zq::ONE, Zq::new(2)].into();
        // let poly6: Rq = vec![Zq::ONE, Zq::ONE, Zq::new(7), Zq::new(5)].into();
        // let result3 = poly5 * poly6;
        // assert_eq!(
        //     result3.coeffs,
        //     helper::padded(&[Zq::new(u32::MAX - 12), Zq::new(u32::MAX - 16), Zq::ZERO])
        // );
    }

    // Test subtraction of polynomials
    #[test]
    fn test_sub() {
        // within bounds
        let poly1: Rq = vec![Zq::new(5), Zq::new(10), Zq::new(15), Zq::new(20)].into();
        let poly2: Rq = vec![Zq::new(2), Zq::new(4), Zq::new(6), Zq::new(8)].into();
        let result = &poly1 - &poly2;
        assert_eq!(
            result.coeffs,
            helper::padded(&[Zq::new(3), Zq::new(6), Zq::new(9), Zq::new(12)])
        );

        // Outside of bounds
        let poly3: Rq = vec![Zq::ONE, Zq::ONE, Zq::new(3), Zq::new(2)].into();
        let poly4: Rq = vec![Zq::new(2), Zq::new(4), Zq::new(6), Zq::new(8)].into();
        let result2 = &poly3 - &poly4;
        assert_eq!(
            result2.coeffs,
            helper::padded(&[
                Zq::MAX,
                Zq::new(u32::MAX - 2),
                Zq::new(u32::MAX - 2),
                Zq::new(u32::MAX - 5)
            ])
        );
        // Subtraction with zero polynomial
        let poly5: Rq = vec![Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)].into();
        let poly6: Rq = vec![Zq::ZERO].into();
        let result3 = &poly6 - &poly5;
        let result4 = &poly5 - &poly6;
        assert_eq!(
            result3.coeffs,
            helper::padded(&[
                Zq::MAX,
                Zq::new(u32::MAX - 1),
                Zq::new(u32::MAX - 2),
                Zq::new(u32::MAX - 3)
            ])
        );
        assert_eq!(
            result4.coeffs,
            helper::padded(&[Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)])
        );
    }

    // Test negation of polynomial
    #[test]
    fn test_neg() {
        let poly: Rq = vec![Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)].into();
        let result = -poly;
        assert_eq!(
            result.coeffs,
            helper::padded(&[
                Zq::MAX,
                Zq::new(u32::MAX - 1),
                Zq::new(u32::MAX - 2),
                Zq::new(u32::MAX - 3)
            ])
        );
    }

    // Test scalar multiplication
    #[test]
    fn test_scalar_mul() {
        let poly: Rq = vec![Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)].into();
        let result = &poly * Zq::new(2);
        assert_eq!(
            result.coeffs,
            helper::padded(&[Zq::new(2), Zq::new(4), Zq::new(6), Zq::new(8)])
        );
    }

    // Test polynomial evaluation
    #[test]
    fn test_eval() {
        let poly: Rq = vec![Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)].into();
        let result = poly.eval(Zq::new(2));
        assert_eq!(result, Zq::new(49));
    }

    // Test equality check
    #[test]
    fn test_is_equal() {
        let poly1: Rq = vec![Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)].into();
        let poly2: Rq = vec![Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)].into();
        let poly3: Rq = vec![Zq::new(4), Zq::new(3), Zq::new(2), Zq::ONE].into();
        assert!(poly1.is_equal(&poly2));
        assert!(!poly1.is_equal(&poly3));
    }

    // Test zero polynomial check
    #[test]
    fn test_is_zero_poly() {
        let zero_poly: Rq = vec![Zq::ZERO; 4].into();
        let non_zero_poly: Rq = vec![Zq::ONE, Zq::ZERO, Zq::ZERO, Zq::ZERO].into();
        assert!(zero_poly.is_zero());
        assert!(!non_zero_poly.is_zero());
    }

    #[test]
    fn test_encode_message() {
        // Test successful encoding
        let message = vec![true, false, true, false];
        let encoded = Rq::encode_message(&message).unwrap();
        assert_eq!(
            encoded.coeffs,
            helper::padded(&[Zq::ONE, Zq::ZERO, Zq::ONE, Zq::ZERO])
        );

        // Test message shorter than degree
        let short_message = vec![true, false];
        let encoded_short = Rq::encode_message(&short_message).unwrap();
        assert_eq!(
            encoded_short.coeffs,
            helper::padded(&[Zq::ONE, Zq::ZERO, Zq::ZERO, Zq::ZERO])
        );

        // Test message too long
        let long_message = vec![true; Rq::DEGREE + 2];
        assert!(Rq::encode_message(&long_message).is_none());

        // Test empty message
        let empty_message: Vec<bool> = vec![];
        let encoded_empty = Rq::encode_message(&empty_message).unwrap();
        assert!(encoded_empty.is_zero());
    }

    // Test coefficient extraction
    #[test]
    fn test_get_coefficient() {
        let poly: Rq = vec![Zq::ONE, Zq::ZERO, Zq::new(5), Zq::MAX].into();
        let vec = helper::padded(&[Zq::ONE, Zq::ZERO, Zq::new(5), Zq::MAX]).to_vec();
        assert!(poly.get_coefficients().to_vec() == vec);

        let poly_zero: Rq = vec![Zq::ZERO, Zq::ZERO, Zq::ZERO, Zq::ZERO].into();
        let vec_zero = helper::padded(&[Zq::ZERO, Zq::ZERO, Zq::ZERO, Zq::ZERO]).to_vec();
        assert!(poly_zero.get_coefficients().to_vec() == vec_zero);
    }

    #[test]
    fn test_base2_decomposition() {
        // Test case 1: Base 2 decomposition
        let poly: Rq = vec![Zq::new(5), Zq::new(3), Zq::new(7), Zq::new(1)].into();
        let parts = poly.decompose(Zq::TWO, 2);

        // Part 0: remainders mod 2 (no centering needed for base 2)
        assert_eq!(
            parts[0].coeffs,
            helper::padded(&[
                Zq::ONE, // 5 mod 2 = 1
                Zq::ONE, // 3 mod 2 = 1
                Zq::ONE, // 7 mod 2 = 1
                Zq::ONE, // 1 mod 2 = 1
            ])
        );

        // Part 1: quotients after division by 2
        assert_eq!(
            parts[1].coeffs,
            helper::padded(&[
                Zq::new(2), // 5 div 2 = 2
                Zq::ONE,    // 3 div 2 = 1
                Zq::new(3), // 7 div 2 = 3
                Zq::ZERO,   // 1 div 2 = 0
            ])
        );

        // Verify Base 2 reconstruction coefficient by coefficient
        for i in 0..4 {
            let expected = poly.coeffs[i];
            let actual = parts[0].coeffs[i] + parts[1].coeffs[i] * Zq::TWO;
            assert_eq!(actual, expected, "Base 2: Coefficient {} mismatch", i);
        }
    }

    #[test]
    fn test_base3_decomposition() {
        // Test case: Base 3 decomposition with centering
        let specific_poly: Rq = vec![Zq::new(8), Zq::new(11), Zq::new(4), Zq::new(15)].into();
        let parts = specific_poly.decompose(Zq::new(3), 2);

        // Part 0: centered remainders mod 3
        assert_eq!(
            parts[0].coeffs,
            helper::padded(&[
                Zq::MAX,  // 8 mod 3 = 2 -> -1 (centered)
                Zq::MAX,  // 11 mod 3 = 2 -> -1 (centered)
                Zq::ONE,  // 4 mod 3 = 1 -> 1 (centered)
                Zq::ZERO, // 15 mod 3 = 0 -> 0 (centered)
            ])
        );

        // Part 1: quotients
        assert_eq!(
            parts[1].coeffs,
            helper::padded(&[
                Zq::new(3), // (8 + 1) div 3 = 3
                Zq::new(4), // (11 + 1) div 3 = 4
                Zq::ONE,    // 4 div 3 = 1
                Zq::new(5), // 15 div 3 = 5
            ])
        );

        // Verify Base 3 reconstruction coefficient by coefficient
        for i in 0..4 {
            let expected = specific_poly.coeffs[i];
            let p0 = parts[0].coeffs[i];
            let p1 = parts[1].coeffs[i];
            let actual = p0 + p1 * Zq::new(3);
            assert_eq!(actual, expected, "Base 3: Coefficient {} mismatch", i);
        }
    }

    #[test]
    fn test_decomposition_edge_cases() {
        // Test zero polynomial
        let zero_poly: Rq = vec![Zq::ZERO; 4].into();
        let parts = zero_poly.decompose(Zq::TWO, 2);
        assert!(
            parts.iter().all(|p| p.is_zero()),
            "Zero polynomial decomposition failed"
        );

        // Test single part decomposition
        let simple_poly: Rq = vec![Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)].into();
        let parts = simple_poly.decompose(Zq::TWO, 1);
        assert_eq!(
            parts.get_elements().len(),
            1,
            "Single part decomposition length incorrect"
        );
        assert_eq!(
            parts[0], simple_poly,
            "Single part decomposition value incorrect"
        );
    }

    #[test]
    fn test_large_base_decomposition() {
        // Test decomposition with larger bases (8 and 16)
        let poly: Rq = vec![Zq::new(120), Zq::new(33), Zq::new(255), Zq::new(19)].into();

        // Base 8 decomposition
        let parts_base8 = poly.decompose(Zq::new(8), 2);

        // Part 0: centered remainders mod 8
        assert_eq!(
            parts_base8[0].coeffs,
            helper::padded(&[
                Zq::ZERO,   // 120 mod 8 = 0 -> 0 (centered)
                Zq::ONE,    // 33 mod 8 = 1 -> 1 (centered)
                Zq::MAX,    // 255 mod 8 = 7 -> -1 (centered)
                Zq::new(3), // 19 mod 8 = 3 -> 3 (centered)
            ])
        );

        // Part 1: quotients
        assert_eq!(
            parts_base8[1].coeffs,
            helper::padded(&[
                Zq::new(15), // 120 div 8 = 15
                Zq::new(4),  // 33 div 8 = 4
                Zq::new(32), // (255 + 1) div 8 = 32
                Zq::new(2),  // 19 div 8 = 2
            ])
        );

        // Verify reconstruction coefficient by coefficient
        for i in 0..4 {
            let expected = poly.coeffs[i];
            let p0 = parts_base8[0].coeffs[i];
            let p1 = parts_base8[1].coeffs[i];
            let actual = p0 + p1 * Zq::new(8);
            assert_eq!(actual, expected, "Base 8: Coefficient {} mismatch", i);
        }

        // Base 16 decomposition
        let parts_base16 = poly.decompose(Zq::new(16), 2);

        // Verify reconstruction for base 16
        for i in 0..4 {
            let expected = poly.coeffs[i];
            let p0 = parts_base16[0].coeffs[i];
            let p1 = parts_base16[1].coeffs[i];
            let actual = p0 + p1 * Zq::new(16);
            assert_eq!(actual, expected, "Base 16: Coefficient {} mismatch", i);
        }
    }

    #[test]
    fn test_multi_part_decomposition() {
        // Test with more than 2 parts
        let poly: Rq = vec![Zq::new(123), Zq::new(456), Zq::new(789), Zq::new(101112)].into();

        // Decompose into 3 parts with base 4
        let parts = poly.decompose(Zq::new(4), 3);
        assert_eq!(parts.get_elements().len(), 3, "Should have 3 parts");

        // Test reconstruction with all 3 parts
        let reconstructed =
            &(&parts[0] + &(&(&parts[1] * Zq::new(4)) + &(&parts[2] * Zq::new(16)))); // 4²

        // Verify reconstruction coefficient by coefficient
        for i in 0..4 {
            assert_eq!(
                reconstructed.coeffs[i], poly.coeffs[i],
                "3-part base 4: Coefficient {} mismatch",
                i
            );
        }
    }

    #[test]
    fn test_centering_properties() {
        // Test that centering works correctly for various values
        // Using base 5 which has half_base = 2
        let values = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let poly: Rq = values
            .iter()
            .map(|&v| Zq::new(v))
            .collect::<Vec<Zq>>()
            .into();

        let parts = poly.decompose(Zq::new(5), 2);

        // Expected centered values for each coefficient:
        // 0 mod 5 = 0 -> 0
        // 1 mod 5 = 1 -> 1
        // 2 mod 5 = 2 -> 2 (at threshold)
        // 3 mod 5 = 3 -> -2 (centered)
        // 4 mod 5 = 4 -> -1 (centered)
        // 5 mod 5 = 0 -> 0
        // ... and so on
        let expected_centered = [
            Zq::ZERO,    // 0 centered
            Zq::ONE,     // 1 centered
            Zq::new(2),  // 2 centered (at threshold)
            -Zq::new(2), // 3 centered to -2
            -Zq::ONE,    // 4 centered to -1
            Zq::ZERO,    // 5 centered
            Zq::ONE,     // 6 centered
            Zq::new(2),  // 7 centered
            -Zq::new(2), // 8 centered to -2
            -Zq::ONE,    // 9 centered to -1
            Zq::ZERO,    // 10 centered
        ];

        for (i, &expected) in expected_centered.iter().enumerate() {
            assert_eq!(
                parts[0].coeffs[i], expected,
                "Base 5 centering: Coefficient {} incorrectly centered",
                i
            );
        }
    }

    #[test]
    fn test_extreme_values() {
        // Test with values near the extremes of the Zq range
        let poly: Rq = vec![Zq::ZERO, Zq::MAX, Zq::MAX - Zq::ONE].into();

        // Decompose with base 3
        let parts = poly.decompose(Zq::new(3), 2);

        // Verify reconstruction
        let reconstructed = &parts[0] + &(&parts[1] * Zq::new(3));

        for i in 0..3 {
            assert_eq!(
                reconstructed.coeffs[i], poly.coeffs[i],
                "Extreme values: Coefficient {} mismatch",
                i
            );
        }

        // Corrected test for high value divisibility
        // u32::MAX = 4294967295, which equals 1431655765 * 3 + 0
        // So u32::MAX mod 3 = 0, which remains 0 (no centering needed)
        assert_eq!(parts[0].coeffs[1], Zq::ZERO); // Remainder after division by 3
        assert_eq!(parts[1].coeffs[1], Zq::new(1431655765)); // Quotient

        // Check u32::MAX - 1 as well
        // 4294967294 mod 3 = 1, which remains 1 (no centering needed since 1 <= half_base)
        assert_eq!(parts[0].coeffs[2], Zq::MAX); // u32::MAX - 1 is the third coefficient
        assert_eq!(parts[1].coeffs[2], Zq::new(1431655765)); // Should be same quotient
    }

    #[test]
    fn test_decomposition_properties() {
        // Test the algebraic property that all coefficients in first part should be small
        let poly: Rq = vec![
            Zq::new(100),
            Zq::new(200),
            Zq::new(300),
            Zq::new(400),
            Zq::new(500),
            Zq::new(600),
            Zq::new(700),
            Zq::new(800),
        ]
        .into();

        for base in [2, 3, 4, 5, 8, 10, 16].iter() {
            let parts = poly.decompose(Zq::new(*base), 2);
            let half_base = Zq::new(*base).scale_by(Zq::TWO);

            // Check that all coefficients in first part are properly "small"
            for coeff in parts[0].coeffs.iter() {
                // In centered representation, all coefficients should be <= half_base
                let abs_coeff = if *coeff > Zq::new(u32::MAX / 2) {
                    Zq::ZERO - *coeff // Handle negative values (represented as large positive ones)
                } else {
                    *coeff
                };

                assert!(
                    abs_coeff <= half_base,
                    "Base {}: First part coefficient {} exceeds half-base {}",
                    base,
                    coeff,
                    half_base
                );
            }

            // Verify reconstruction
            let reconstructed = &parts[0] + &(&parts[1] * Zq::new(*base));
            assert_eq!(reconstructed, poly, "Base {}: Reconstruction failed", base);
        }
    }

    #[test]
    fn test_conjugate_automorphism() {
        let poly1 = helper::rq_from(&[Zq::ONE, Zq::TWO, Zq::new(3)]);
        let poly2 = helper::rq_from(&[Zq::new(4), Zq::new(5), Zq::new(6)]);
        let inner_12 = poly1.inner_product(&poly2);
        let conjugated_1 = poly1.conjugate_automorphism();
        let inner_conjugated_12 = &conjugated_1 * &poly2;

        assert_eq!(inner_conjugated_12.coeffs.len(), Rq::DEGREE);
        assert_eq!(inner_conjugated_12.get_coefficients()[0], Zq::new(32));
        assert_eq!(inner_conjugated_12.get_coefficients()[1], Zq::new(17));
        assert_eq!(inner_conjugated_12.get_coefficients()[2], Zq::new(6));

        // ct<\sigma_{-1}(poly1), poly2> ?= <poly1, poly2>
        let ct_inner_conjugated_12 = inner_conjugated_12.get_coefficients()[0];
        assert_eq!(ct_inner_conjugated_12, inner_12);
    }
}
