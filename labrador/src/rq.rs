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
use crate::zq::Zq;
use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use rand::distr::{Distribution, Uniform};
use rand::{CryptoRng, Rng};
use std::iter::Sum;

/// This module provides implementations for various operations
/// in the polynomial ring R = Z_q\[X\] / (X^d + 1).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Rq<const D: usize> {
    coeffs: [Zq; D],
}

impl<const D: usize> Rq<D> {
    /// Constructor for the polynomial ring
    pub const fn new(coeffs: [Zq; D]) -> Self {
        Rq { coeffs }
    }
    /// Get the coefficients as a vector
    pub fn get_coefficients(&self) -> &[Zq; D] {
        &self.coeffs
    }

    /// Polynomial addition
    fn addition(&self, other: &Self) -> Self {
        let mut result = [Zq::ZERO; D];
        for (r, (a, b)) in result
            .iter_mut()
            .zip(self.coeffs.iter().zip(other.coeffs.iter()))
        {
            *r = *a + *b;
        }
        Rq::new(result)
    }

    /// Polynomial subtraction
    fn subtraction(&self, other: &Self) -> Self {
        let mut result = [Zq::ZERO; D];
        for (r, (a, b)) in result
            .iter_mut()
            .zip(self.coeffs.iter().zip(other.coeffs.iter()))
        {
            *r = *a - *b;
        }
        Rq::new(result)
    }

    /// Polynomial multiplication modulo x^D + 1
    fn multiplication(&self, other: &Self) -> Self {
        let mut result = [Zq::ZERO; D];
        let mut out_of_field = [Zq::ZERO; D];
        for (i, &self_coeff) in self.coeffs.iter().enumerate() {
            for (j, &other_coeff) in other.coeffs.iter().enumerate() {
                if i + j < D {
                    result[i + j] += self_coeff * other_coeff;
                } else {
                    out_of_field[(i + j) % D] += self_coeff * other_coeff;
                }
            }
        }
        // Process excess terms with sign adjustment
        for i in (0..D).rev() {
            let m = i / D;
            let r = i % D;
            let sign = if (m + 1) % 2 == 0 { 1 } else { -1 };
            if sign == 1 {
                result[r] += out_of_field[i];
            } else {
                result[r] -= out_of_field[i];
            }
        }
        Rq::new(result)
    }

    /// Dot product between coefficients
    pub fn inner_product(&self, other: &Self) -> Zq {
        self.coeffs
            .iter()
            .zip(other.coeffs.iter())
            .map(|(&a, &b)| a * b)
            .fold(Zq::ZERO, |acc, x| acc + x)
    }

    /// Scalar multiplication
    pub fn scalar_mul(&self, s: Zq) -> Self {
        let mut result = [Zq::ZERO; D];
        for (i, &coeff) in self.coeffs.iter().enumerate() {
            result[i] = s * (coeff);
        }
        Rq::new(result)
    }

    /// Evaluate the polynomial at a specific point
    pub fn eval(&self, x: Zq) -> Zq {
        let mut result = Zq::ZERO;
        for coeff in self.coeffs.iter().rev() {
            result = result * x + *coeff;
        }

        result
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
        let mut coeffs = [Zq::ZERO; D];
        coeffs.iter_mut().for_each(|c| *c = uniform.sample(rng));
        Self { coeffs }
    }

    /// Generate random small polynomial with secure RNG implementation
    pub fn random_ternary<R: Rng + CryptoRng>(rng: &mut R) -> Self {
        let mut coeffs = [Zq::ZERO; D];

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
    pub fn decompose(&self, base: Zq, num_parts: usize) -> Vec<Self> {
        let mut parts = Vec::with_capacity(num_parts);
        let mut current = self.clone();

        for i in 0..num_parts {
            if i == num_parts - 1 {
                parts.push(current.clone());
            } else {
                // Extract low part (mod base, centered around 0)
                let mut low_coeffs = [Zq::ZERO; D];

                for (j, coeff) in current.get_coefficients().iter().enumerate() {
                    low_coeffs[j] = coeff.centered_mod(base);
                }

                let low_part = Self::new(low_coeffs);
                parts.push(low_part.clone());

                // Update current
                current -= low_part;

                // Scale by base
                let mut scaled_coeffs = [Zq::ZERO; D];
                for (j, coeff) in current.get_coefficients().iter().enumerate() {
                    scaled_coeffs[j] = coeff.scale_by(base);
                }
                current = Self::new(scaled_coeffs);
            }
        }

        parts
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
        if message.len() > D {
            return None;
        }

        let mut coeffs = [Zq::ZERO; D];
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

    pub const fn zero() -> Self {
        Self::new([Zq::ZERO; D])
    }
}

macro_rules! impl_arithmetic {
    ($trait:ident, $assign_trait:ident, $method:ident, $assign_method:ident, $op_method:ident) => {
        impl<const D: usize> $trait for Rq<{ D }> {
            type Output = Self;

            fn $method(self, rhs: Self) -> Self::Output {
                self.$op_method(&rhs)
            }
        }

        impl<const D: usize> $assign_trait for Rq<{ D }> {
            fn $assign_method(&mut self, rhs: Self) {
                let result = self.$op_method(&rhs);
                self.coeffs = result.coeffs;
            }
        }
    };
}

impl_arithmetic!(Add, AddAssign, add, add_assign, addition);
impl_arithmetic!(Sub, SubAssign, sub, sub_assign, subtraction);
impl_arithmetic!(Mul, MulAssign, mul, mul_assign, multiplication);

impl<const D: usize> From<Vec<Zq>> for Rq<D> {
    fn from(vec: Vec<Zq>) -> Self {
        let mut temp = [Zq::ZERO; D];
        // Process excess terms with sign adjustment
        for i in (0..vec.len()).rev() {
            let m = i / D;
            let r = i % D;
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

impl Sum for Zq {
    // Accumulate using the addition operator
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = Zq>,
    {
        iter.fold(Zq::ZERO, |acc, x| acc + x)
    }
}

// Implementing the Neg trait
impl<const D: usize> Neg for Rq<D> {
    type Output = Self;

    /// Polynomial negation
    fn neg(self) -> Self {
        let mut result = [Zq::ZERO; D];
        for (i, &coeff) in self.coeffs.iter().enumerate() {
            result[i] = Zq::ZERO - coeff;
        }
        Rq::new(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Test new() and polynomial creation
    #[test]
    fn test_new_and_create_poly() {
        let poly = Rq::new([Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)]);
        assert_eq!(poly.coeffs, [Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)]);

        // Direct conversion
        let poly_from_vec_direct: Rq<4> = vec![Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)].into();
        assert_eq!(
            poly_from_vec_direct.coeffs,
            [Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)]
        );
        // Wrapping around
        let poly_from_vec_wrapping: Rq<4> =
            vec![Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4), Zq::ONE].into();
        assert_eq!(
            poly_from_vec_wrapping.coeffs,
            [Zq::ZERO, Zq::new(2), Zq::new(3), Zq::new(4)]
        );
        // Filling up with zeros
        let poly_from_vec_zeros: Rq<4> = vec![Zq::ONE, Zq::new(2)].into();
        assert_eq!(
            poly_from_vec_zeros.coeffs,
            [Zq::ONE, Zq::new(2), Zq::ZERO, Zq::ZERO]
        );
        // High-Degree Term Reduction
        let poly_high_degree_reduction: Rq<2> = vec![
            Zq::ONE,
            Zq::new(2),
            Zq::ZERO,
            Zq::ZERO,
            Zq::ZERO,
            Zq::ZERO,
            Zq::ONE,
        ]
        .into();
        assert_eq!(poly_high_degree_reduction.coeffs, [Zq::ZERO, Zq::new(2)]);
    }

    // Test addition of polynomials
    #[test]
    fn test_add() {
        // Within bounds
        let poly1: Rq<4> = vec![Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)].into();
        let poly2: Rq<4> = vec![Zq::new(4), Zq::new(3), Zq::new(2), Zq::ONE].into();
        let result = poly1 + poly2;
        assert_eq!(
            result.coeffs,
            [Zq::new(5), Zq::new(5), Zq::new(5), Zq::new(5)]
        );

        // Outside of bounds
        let poly3: Rq<4> = vec![Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)].into();
        let poly4: Rq<4> = vec![Zq::MAX, Zq::new(3), Zq::MAX, Zq::ONE].into();
        let result2 = poly3 + poly4;
        assert_eq!(
            result2.coeffs,
            [Zq::ZERO, Zq::new(5), Zq::new(2), Zq::new(5)]
        );
        // Addition with zero polynomial
        let poly5: Rq<4> = vec![Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)].into();
        let poly6: Rq<4> = vec![Zq::ZERO].into();
        let result3 = poly5 + poly6;
        assert_eq!(
            result3.coeffs,
            [Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)]
        );
        // Addition with high coefficients
        let poly7: Rq<4> = vec![Zq::ONE, Zq::new(2), Zq::new(3), Zq::MAX].into();
        let poly8: Rq<4> = vec![Zq::MAX, Zq::MAX, Zq::MAX, Zq::MAX].into();
        let result3 = poly7 + poly8;
        assert_eq!(
            result3.coeffs,
            [
                Zq::ZERO,
                Zq::ONE,
                Zq::new(2),
                Zq::new(u32::MAX.wrapping_add(u32::MAX))
            ]
        );
    }
    // Test multiplication of polynomials
    #[test]

    fn test_mul() {
        // Multiplication with wrapping
        let poly1: Rq<3> = vec![Zq::ONE, Zq::ONE, Zq::new(2)].into();
        let poly2: Rq<3> = vec![Zq::ONE, Zq::ONE].into();
        let result = poly1 * poly2;
        assert_eq!(result.coeffs, [Zq::MAX, Zq::new(2), Zq::new(3)]);

        // Multiplication with zero polynomial
        let poly3: Rq<3> = vec![Zq::ONE, Zq::ONE, Zq::new(2)].into();
        let poly4: Rq<3> = vec![Zq::ZERO].into();
        let result2 = poly3 * poly4;
        assert_eq!(result2.coeffs, [Zq::ZERO, Zq::ZERO, Zq::ZERO]);

        // Multiplication with wrapping higher order
        let poly5: Rq<3> = vec![Zq::ONE, Zq::ONE, Zq::new(2)].into();
        let poly6: Rq<3> = vec![Zq::ONE, Zq::ONE, Zq::new(7), Zq::new(5)].into();
        let result3 = poly5 * poly6;
        assert_eq!(
            result3.coeffs,
            [Zq::new(u32::MAX - 12), Zq::new(u32::MAX - 16), Zq::ZERO]
        );
    }

    // Test subtraction of polynomials
    #[test]
    fn test_sub() {
        // within bounds
        let poly1: Rq<4> = vec![Zq::new(5), Zq::new(10), Zq::new(15), Zq::new(20)].into();
        let poly2: Rq<4> = vec![Zq::new(2), Zq::new(4), Zq::new(6), Zq::new(8)].into();
        let result = poly1 - poly2;
        assert_eq!(
            result.coeffs,
            [Zq::new(3), Zq::new(6), Zq::new(9), Zq::new(12)]
        );

        // Outside of bounds
        let poly3: Rq<4> = vec![Zq::ONE, Zq::ONE, Zq::new(3), Zq::new(2)].into();
        let poly4: Rq<4> = vec![Zq::new(2), Zq::new(4), Zq::new(6), Zq::new(8)].into();
        let result2 = poly3 - poly4;
        assert_eq!(
            result2.coeffs,
            [
                Zq::MAX,
                Zq::new(u32::MAX - 2),
                Zq::new(u32::MAX - 2),
                Zq::new(u32::MAX - 5)
            ]
        );
        // Subtraction with zero polynomial
        let poly5: Rq<4> = vec![Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)].into();
        let poly6: Rq<4> = vec![Zq::ZERO].into();
        let result3 = poly6.clone() - poly5.clone();
        let result4 = poly5.clone() - poly6.clone();
        assert_eq!(
            result3.coeffs,
            [
                Zq::MAX,
                Zq::new(u32::MAX - 1),
                Zq::new(u32::MAX - 2),
                Zq::new(u32::MAX - 3)
            ]
        );
        assert_eq!(
            result4.coeffs,
            [Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)]
        );
    }

    // Test negation of polynomial
    #[test]
    fn test_neg() {
        let poly: Rq<4> = vec![Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)].into();
        let result = -poly;
        assert_eq!(
            result.coeffs,
            [
                Zq::MAX,
                Zq::new(u32::MAX - 1),
                Zq::new(u32::MAX - 2),
                Zq::new(u32::MAX - 3)
            ]
        );
    }

    // Test scalar multiplication
    #[test]
    fn test_scalar_mul() {
        let poly: Rq<4> = vec![Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)].into();
        let result = poly.scalar_mul(Zq::new(2));
        assert_eq!(
            result.coeffs,
            [Zq::new(2), Zq::new(4), Zq::new(6), Zq::new(8)]
        );
    }

    // Test polynomial evaluation
    #[test]
    fn test_eval() {
        let poly: Rq<4> = vec![Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)].into();
        let result = poly.eval(Zq::new(2));
        assert_eq!(result, Zq::new(49));
    }

    // Test equality check
    #[test]
    fn test_is_equal() {
        let poly1: Rq<4> = vec![Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)].into();
        let poly2: Rq<4> = vec![Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)].into();
        let poly3: Rq<4> = vec![Zq::new(4), Zq::new(3), Zq::new(2), Zq::ONE].into();
        assert!(poly1.is_equal(&poly2));
        assert!(!poly1.is_equal(&poly3));
    }

    // Test zero polynomial check
    #[test]
    fn test_is_zero_poly() {
        let zero_poly: Rq<4> = vec![Zq::ZERO; 4].into();
        let non_zero_poly: Rq<4> = vec![Zq::ONE, Zq::ZERO, Zq::ZERO, Zq::ZERO].into();
        assert!(zero_poly.is_zero());
        assert!(!non_zero_poly.is_zero());
    }

    #[test]
    fn test_encode_message() {
        // Test successful encoding
        let message = vec![true, false, true, false];
        let encoded = Rq::<4>::encode_message(&message).unwrap();
        assert_eq!(encoded.coeffs, [Zq::ONE, Zq::ZERO, Zq::ONE, Zq::ZERO]);

        // Test message shorter than degree
        let short_message = vec![true, false];
        let encoded_short = Rq::<4>::encode_message(&short_message).unwrap();
        assert_eq!(
            encoded_short.coeffs,
            [Zq::ONE, Zq::ZERO, Zq::ZERO, Zq::ZERO]
        );

        // Test message too long
        let long_message = vec![true; 5];
        assert!(Rq::<4>::encode_message(&long_message).is_none());

        // Test empty message
        let empty_message: Vec<bool> = vec![];
        let encoded_empty = Rq::<4>::encode_message(&empty_message).unwrap();
        assert!(encoded_empty.is_zero());
    }

    // Test coefficient extraction
    #[test]
    fn test_get_coefficient() {
        let poly: Rq<4> = vec![Zq::ONE, Zq::ZERO, Zq::new(5), Zq::MAX].into();
        let vec = vec![Zq::ONE, Zq::ZERO, Zq::new(5), Zq::MAX];
        assert!(poly.get_coefficients().to_vec() == vec);

        let poly_zero: Rq<4> = vec![Zq::ZERO, Zq::ZERO, Zq::ZERO, Zq::ZERO].into();
        let vec_zero = vec![Zq::ZERO, Zq::ZERO, Zq::ZERO, Zq::ZERO];
        assert!(poly_zero.get_coefficients().to_vec() == vec_zero);
    }

    #[test]
    fn test_base2_decomposition() {
        // Test case 1: Base 2 decomposition
        let poly: Rq<4> = vec![Zq::new(5), Zq::new(3), Zq::new(7), Zq::new(1)].into();
        let parts = poly.decompose(Zq::TWO, 2);

        // Part 0: remainders mod 2 (no centering needed for base 2)
        assert_eq!(
            parts[0].coeffs,
            [
                Zq::ONE, // 5 mod 2 = 1
                Zq::ONE, // 3 mod 2 = 1
                Zq::ONE, // 7 mod 2 = 1
                Zq::ONE, // 1 mod 2 = 1
            ]
        );

        // Part 1: quotients after division by 2
        assert_eq!(
            parts[1].coeffs,
            [
                Zq::new(2), // 5 div 2 = 2
                Zq::ONE,    // 3 div 2 = 1
                Zq::new(3), // 7 div 2 = 3
                Zq::ZERO,   // 1 div 2 = 0
            ]
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
        let specific_poly: Rq<4> = vec![Zq::new(8), Zq::new(11), Zq::new(4), Zq::new(15)].into();
        let parts = specific_poly.decompose(Zq::new(3), 2);

        // Part 0: centered remainders mod 3
        assert_eq!(
            parts[0].coeffs,
            [
                Zq::MAX,  // 8 mod 3 = 2 -> -1 (centered)
                Zq::MAX,  // 11 mod 3 = 2 -> -1 (centered)
                Zq::ONE,  // 4 mod 3 = 1 -> 1 (centered)
                Zq::ZERO, // 15 mod 3 = 0 -> 0 (centered)
            ]
        );

        // Part 1: quotients
        assert_eq!(
            parts[1].coeffs,
            [
                Zq::new(3), // (8 + 1) div 3 = 3
                Zq::new(4), // (11 + 1) div 3 = 4
                Zq::ONE,    // 4 div 3 = 1
                Zq::new(5), // 15 div 3 = 5
            ]
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
        let zero_poly: Rq<4> = vec![Zq::ZERO; 4].into();
        let parts = zero_poly.decompose(Zq::TWO, 2);
        assert!(
            parts.iter().all(|p| p.is_zero()),
            "Zero polynomial decomposition failed"
        );

        // Test single part decomposition
        let simple_poly: Rq<4> = vec![Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)].into();
        let parts = simple_poly.decompose(Zq::TWO, 1);
        assert_eq!(parts.len(), 1, "Single part decomposition length incorrect");
        assert_eq!(
            parts[0], simple_poly,
            "Single part decomposition value incorrect"
        );
    }

    #[test]
    fn test_large_base_decomposition() {
        // Test decomposition with larger bases (8 and 16)
        let poly: Rq<4> = vec![Zq::new(120), Zq::new(33), Zq::new(255), Zq::new(19)].into();

        // Base 8 decomposition
        let parts_base8 = poly.decompose(Zq::new(8), 2);

        // Part 0: centered remainders mod 8
        assert_eq!(
            parts_base8[0].coeffs,
            [
                Zq::ZERO,   // 120 mod 8 = 0 -> 0 (centered)
                Zq::ONE,    // 33 mod 8 = 1 -> 1 (centered)
                Zq::MAX,    // 255 mod 8 = 7 -> -1 (centered)
                Zq::new(3), // 19 mod 8 = 3 -> 3 (centered)
            ]
        );

        // Part 1: quotients
        assert_eq!(
            parts_base8[1].coeffs,
            [
                Zq::new(15), // 120 div 8 = 15
                Zq::new(4),  // 33 div 8 = 4
                Zq::new(32), // (255 + 1) div 8 = 32
                Zq::new(2),  // 19 div 8 = 2
            ]
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
        let poly: Rq<4> = vec![Zq::new(123), Zq::new(456), Zq::new(789), Zq::new(101112)].into();

        // Decompose into 3 parts with base 4
        let parts = poly.decompose(Zq::new(4), 3);
        assert_eq!(parts.len(), 3, "Should have 3 parts");

        // Test reconstruction with all 3 parts
        let reconstructed = parts[0].clone()
            + parts[1].clone().scalar_mul(Zq::new(4))
            + parts[2].clone().scalar_mul(Zq::new(16)); // 4²

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
        let poly: Rq<11> = values
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
        let poly: Rq<3> = vec![Zq::ZERO, Zq::MAX, Zq::MAX - Zq::ONE].into();

        // Decompose with base 3
        let parts = poly.decompose(Zq::new(3), 2);

        // Verify reconstruction
        let reconstructed = parts[0].clone() + parts[1].clone().scalar_mul(Zq::new(3));

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
        let poly: Rq<8> = vec![
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
            let reconstructed = parts[0].clone() + parts[1].clone().scalar_mul(Zq::new(*base));
            assert_eq!(reconstructed, poly, "Base {}: Reconstruction failed", base);
        }
    }
}
