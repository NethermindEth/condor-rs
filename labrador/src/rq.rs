// This file is part of the polynomial ring operations module.
//
//
// Currently implemented functions include:
// - Polynomial addition:         +
// - Polynomial multiplication:   *
// - inner_product/ Dot product:  inner_product()
// - Polynomial subtraction:      -
// - Polynomial negation:         neg()
// - Scalar multiplication:       scalar_mul()
// - Polynomial evaluation:       eval()
// - Zero check:                  is_zero()
// - Polynomial equality check:   is_equal()
//
// Further operations and optimizations will be added in future versions.

// We use the Zq ring
use crate::zq::Zq;
use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use rand::distr::{Distribution, Uniform};
use rand::{CryptoRng, Rng};

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
        self.iter().all(|coeff| {
            let val = coeff.value();
            // Check if value is within [-bound, bound]
            val <= bound.value() || val >= Zq::Q.wrapping_sub(bound.value())
        })
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
}
