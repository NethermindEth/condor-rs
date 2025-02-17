// This file is part of the polynomial ring operations module.
//
// This module provides implementations for various operations
// in the polynomial ring R = Z_q[X] / (X^d + 1).
//
// Currently implemented functions include:
// - Polynomial addition:         +
// - Polynomial multiplication:   *
// - inner_product/ Dot product:  inner_product()
// - Polynomial subtraction:      -
// - Polynomial negation:         neg()
// - Scalar multiplication:       scalar_mul()
// - Division by monomials:       div_by_monomial()
// - Polynomial evaluation:       eval()
// - Zero check:                  is_zero()
// - Polynomial equality check:   is_equal()
//
// Further operations and optimizations will be added in future versions.

// We use the Zq ring
use labrador::zq::Zq;
use std::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Rq<const D: usize> {
    coeffs: [Zq; D],
}

impl<const D: usize> Rq<D> {
    // Constructor for the polynomial ring
    pub fn new(coeffs: [Zq; D]) -> Self {
        Rq { coeffs }
    }

    // Polynomial addition
    pub fn _add(&self, other: &Self) -> Self {
        let mut result = [Zq::zero(); D];
        for (r, (a, b)) in result
            .iter_mut()
            .zip(self.coeffs.iter().zip(other.coeffs.iter()))
        {
            *r = *a + *b;
        }
        Rq::new(result)
    }

    // Polynomial subtraction
    pub fn _sub(&self, other: &Self) -> Self {
        let mut result = [Zq::zero(); D];
        for (r, (a, b)) in result
            .iter_mut()
            .zip(self.coeffs.iter().zip(other.coeffs.iter()))
        {
            *r = *a - *b;
        }
        Rq::new(result)
    }

    // Polynomial multiplication modulo x^D + 1
    pub fn _mul(&self, other: &Self) -> Self {
        // Replace conditional subtraction with full reduction:
        let mut temp = vec![Zq::zero(); 2 * D - 1];
        for i in 0..D {
            for j in 0..D {
                temp[i + j] += self.coeffs[i] * other.coeffs[j];
            }
        }
        let result: Rq<D> = temp.into();

        result
    }

    // Dot product between coefficients
    pub fn inner_product(&self, other: &Self) -> Zq {
        self.coeffs
            .iter()
            .zip(other.coeffs.iter())
            .map(|(&a, &b)| a * (b))
            .fold(Zq::zero(), |acc, x| acc + x)
    }

    // Polynomial negation
    pub fn neg(&self) -> Self {
        let mut result = [Zq::zero(); D];
        for (i, &coeff) in self.coeffs.iter().enumerate() {
            result[i] = Zq::zero() - coeff;
        }
        Rq::new(result)
    }

    // Scalar multiplication
    pub fn scalar_mul(&self, s: Zq) -> Self {
        let mut result = [Zq::zero(); D];
        for (i, &coeff) in self.coeffs.iter().enumerate() {
            result[i] = s * (coeff);
        }
        Rq::new(result)
    }

    // Evaluate the polynomial at a specific point
    pub fn eval(&self, x: Zq) -> Zq {
        let mut result = Zq::zero();
        let mut power = Zq::one();
        for &coeff in &self.coeffs {
            result += coeff * power;
            power *= x;
        }

        result
    }

    // Check if Polynomial == 0
    pub fn is_zero(&self) -> bool {
        self.coeffs.iter().all(|&coeff| coeff == Zq::zero())
    }

    // Check if two polynomials are equal
    pub fn is_equal(&self, other: &Self) -> bool {
        self.coeffs == other.coeffs
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

// Implement the traits using your custom methods
impl_arithmetic!(Add, AddAssign, add, add_assign, _add);
impl_arithmetic!(Sub, SubAssign, sub, sub_assign, _sub);
impl_arithmetic!(Mul, MulAssign, mul, mul_assign, _mul);

impl<const D: usize> From<Vec<Zq>> for Rq<D> {
    fn from(vec: Vec<Zq>) -> Self {
        let mut temp = [Zq::zero(); D];
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

#[cfg(test)]
mod tests {
    use super::*;

    // Test new() and polynomial creation
    #[test]
    fn test_new_and_create_poly() {
        let poly = Rq::new([Zq::new(1), Zq::new(2), Zq::new(3), Zq::new(4)]);
        assert_eq!(
            poly.coeffs,
            [Zq::new(1), Zq::new(2), Zq::new(3), Zq::new(4)]
        );

        // Positive coefficients
        let poly_from_vec: Rq<4> = vec![Zq::new(1), Zq::new(2), Zq::new(3), Zq::new(4)].into();
        assert_eq!(
            poly_from_vec.coeffs,
            [Zq::new(1), Zq::new(2), Zq::new(3), Zq::new(4)]
        );
    }

    // Test addition of polynomials
    #[test]
    fn test_add() {
        // within bounds
        let poly1: Rq<4> = vec![Zq::new(1), Zq::new(2), Zq::new(3), Zq::new(4)].into();
        let poly2: Rq<4> = vec![Zq::new(4), Zq::new(3), Zq::new(2), Zq::new(1)].into();
        let result = poly1 + poly2;
        assert_eq!(
            result.coeffs,
            [Zq::new(5), Zq::new(5), Zq::new(5), Zq::new(5)]
        );

        // Outside of bounds
        let poly3: Rq<4> = vec![Zq::new(1), Zq::new(2), Zq::new(3), Zq::new(4)].into();
        let poly4 = Rq::<4>::new([Zq::new(u32::MAX), Zq::new(3), Zq::new(u32::MAX), Zq::new(1)]);
        let result2 = poly3 + poly4;
        assert_eq!(
            result2.coeffs,
            [Zq::new(0), Zq::new(5), Zq::new(2), Zq::new(5)]
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
        let poly3: Rq<4> = vec![Zq::new(1), Zq::new(1), Zq::new(3), Zq::new(2)].into();
        let poly4: Rq<4> = vec![Zq::new(2), Zq::new(4), Zq::new(6), Zq::new(8)].into();
        let result2 = poly3 - poly4;
        assert_eq!(
            result2.coeffs,
            [
                Zq::new(u32::MAX),
                Zq::new(u32::MAX - 2),
                Zq::new(u32::MAX - 2),
                Zq::new(u32::MAX - 5)
            ]
        );
    }

    // Test negation of polynomial
    #[test]
    fn test_neg() {
        let poly: Rq<4> = vec![Zq::new(1), Zq::new(2), Zq::new(3), Zq::new(4)].into();
        let result = poly.neg();
        assert_eq!(
            result.coeffs,
            [
                Zq::new(u32::MAX),
                Zq::new(u32::MAX - 1),
                Zq::new(u32::MAX - 2),
                Zq::new(u32::MAX - 3)
            ]
        );
    }

    // Test scalar multiplication
    #[test]
    fn test_scalar_mul() {
        let poly: Rq<4> = vec![Zq::new(1), Zq::new(2), Zq::new(3), Zq::new(4)].into();
        let result = poly.scalar_mul(Zq::new(2));
        assert_eq!(
            result.coeffs,
            [Zq::new(2), Zq::new(4), Zq::new(6), Zq::new(8)]
        );
    }

    // Test division by monomials
    #[test]
    fn test_div_by_monomial() {
        let poly: Rq<4> = vec![Zq::new(1), Zq::new(2), Zq::new(3), Zq::new(4)].into();
        let result = poly.div_by_monomial(2);
        assert_eq!(
            result.coeffs,
            [Zq::new(3), Zq::new(4), Zq::new(1), Zq::new(2)]
        );
    }

    // Test polynomial evaluation
    #[test]
    fn test_eval() {
        let poly: Rq<4> = vec![Zq::new(1), Zq::new(2), Zq::new(3), Zq::new(4)].into();
        let result = poly.eval(Zq::new(2));
        assert_eq!(result, Zq::new(49));
    }

    // Test equality check
    #[test]
    fn test_is_equal() {
        let poly1: Rq<4> = vec![Zq::new(1), Zq::new(2), Zq::new(3), Zq::new(4)].into();
        let poly2: Rq<4> = vec![Zq::new(1), Zq::new(2), Zq::new(3), Zq::new(4)].into();
        let poly3: Rq<4> = vec![Zq::new(4), Zq::new(3), Zq::new(2), Zq::new(1)].into();
        assert!(poly1.is_equal(&poly2));
        assert!(!poly1.is_equal(&poly3));
    }

    // Test zero polynomial check
    #[test]
    fn test_is_zero_poly() {
        let zero_poly: Rq<4> = vec![Zq::zero(); 4].into();
        let non_zero_poly: Rq<4> = vec![Zq::new(1), Zq::zero(), Zq::zero(), Zq::zero()].into();
        assert!(zero_poly.is_zero());
        assert!(!non_zero_poly.is_zero());
    }
}
