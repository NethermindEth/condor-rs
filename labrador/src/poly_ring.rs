// This file is part of the polynomial ring operations module.
//
// This module provides implementations for various operations
// in the polynomial ring R = Z_q[X] / (X^d + 1).
//
// Currently implemented functions include:
// - Polynomial addition:         add()
// - Polynomial multiplication:   mul()
// - inner_product/ Dot product:  inner_product()
// - Polynomial subtraction:     sub()
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Poly<const D: usize> {
    coeffs: [Zq; D],
}

impl<const D: usize> Poly<D> {
    // Constructor for the polynomial ring
    pub fn new(coeffs: [Zq; D]) -> Self {
        Poly { coeffs }
    }

    // Polynomial addition
    pub fn add(&self, other: &Self) -> Self {
        let mut result = [Zq::zero(); D];
        for (r, (a, b)) in result
            .iter_mut()
            .zip(self.coeffs.iter().zip(other.coeffs.iter()))
        {
            *r = *a + *b;
        }
        Poly::new(result)
    }

    // Polynomial subtraction
    pub fn sub(&self, other: &Self) -> Self {
        let mut result = [Zq::zero(); D];
        for (r, (a, b)) in result
            .iter_mut()
            .zip(self.coeffs.iter().zip(other.coeffs.iter()))
        {
            *r = *a - *b;
        }
        Poly::new(result)
    }

    // Create a Polynomial from a vector
    pub fn create_poly(coeffs: Vec<i32>) -> Poly<D> {
        let mut arr = [Zq::zero(); D];
        let u32_coeffs: Vec<u32> = coeffs.iter().map(|&coeff| coeff as u32).collect();
        // First D elements are assigned directly
        for (i, &u32_coeffs) in u32_coeffs.iter().take(D).enumerate() {
            arr[i] = Zq::new(u32_coeffs);
        }

        // Handle additional elements by subtracting them at (index % D)
        for (i, &u32_coeffs) in u32_coeffs.iter().skip(D).enumerate() {
            let mod_index = i % D;
            arr[mod_index] -= Zq::new(u32_coeffs);
        }

        Poly::new(arr)
    }

    // Polynomial multiplication modulo x^D + 1
    pub fn mul(&self, other: &Self) -> Self {
        let mut result = [Zq::zero(); D];

        for i in 0..D {
            for j in 0..D {
                let degree = (i + j) % D;
                if (i + j) > D {
                    result[degree] -= self.coeffs[i] * other.coeffs[j];
                } else {
                    // normal multiplication
                    result[degree] -= self.coeffs[i] * other.coeffs[j];
                }
            }
        }

        Poly::new(result)
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
        Poly::new(result)
    }

    // Scalar multiplication
    pub fn scalar_mul(&self, s: Zq) -> Self {
        let mut result = [Zq::zero(); D];
        for (i, &coeff) in self.coeffs.iter().enumerate() {
            result[i] = s * (coeff);
        }
        Poly::new(result)
    }

    // (Division by monomial) X^k | performing a cyclic right shift
    pub fn div_by_monomial(&self, k: usize) -> Self {
        let mut result = [Zq::zero(); D];
        for (i, &coeff) in self.coeffs.iter().enumerate() {
            // (k-shift)
            let new_index = (i + k) % D;
            result[new_index] = coeff;
        }
        Poly::new(result)
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

#[cfg(test)]
mod tests {
    use super::*;

    // Test new() and polynomial creation
    #[test]
    fn test_new_and_create_poly() {
        let poly = Poly::new([Zq::new(1), Zq::new(2), Zq::new(3), Zq::new(4)]);
        assert_eq!(
            poly.coeffs,
            [Zq::new(1), Zq::new(2), Zq::new(3), Zq::new(4)]
        );

        // Positive coefficients
        let poly_from_vec = Poly::<4>::create_poly(vec![1, 2, 3, 4]);
        assert_eq!(
            poly_from_vec.coeffs,
            [Zq::new(1), Zq::new(2), Zq::new(3), Zq::new(4)]
        );

        // Negative coefficients
        let poly_from_vec_negative = Poly::<4>::create_poly(vec![-1, -2, -3, -4]);
        let u32_max = u32::MAX;
        assert_eq!(
            poly_from_vec_negative.coeffs,
            [
                Zq::new(u32_max),
                Zq::new(u32_max - 1),
                Zq::new(u32_max - 2),
                Zq::new(u32_max - 3)
            ]
        );
    }

    // Test addition of polynomials
    #[test]
    fn test_add() {
        // within bounds
        let poly1 = Poly::<4>::create_poly(vec![1, 2, 3, 4]);
        let poly2 = Poly::<4>::create_poly(vec![4, 3, 2, 1]);
        let result = poly1.add(&poly2);
        assert_eq!(
            result.coeffs,
            [Zq::new(5), Zq::new(5), Zq::new(5), Zq::new(5)]
        );

        // Outside of bounds
        let poly3 = Poly::<4>::create_poly(vec![1, 2, 3, 4]);
        let poly4 = Poly::<4>::new([Zq::new(u32::MAX), Zq::new(3), Zq::new(u32::MAX), Zq::new(1)]);
        let result2 = poly3.add(&poly4);
        assert_eq!(
            result2.coeffs,
            [Zq::new(0), Zq::new(5), Zq::new(2), Zq::new(5)]
        );
    }

    // Test subtraction of polynomials
    #[test]
    fn test_sub() {
        // within bounds
        let poly1 = Poly::<4>::create_poly(vec![5, 10, 15, 20]);
        let poly2 = Poly::<4>::create_poly(vec![2, 4, 6, 8]);
        let result = poly1.sub(&poly2);
        assert_eq!(
            result.coeffs,
            [Zq::new(3), Zq::new(6), Zq::new(9), Zq::new(12)]
        );

        // Outside of bounds
        let poly3 = Poly::<4>::create_poly(vec![1, 1, 3, 2]);
        let poly4 = Poly::<4>::create_poly(vec![2, 4, 6, 8]);
        let result2 = poly3.sub(&poly4);
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
        let poly = Poly::<4>::create_poly(vec![1, 2, 3, 4]);
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
        let poly = Poly::<4>::create_poly(vec![1, 2, 3, 4]);
        let result = poly.scalar_mul(Zq::new(2));
        assert_eq!(
            result.coeffs,
            [Zq::new(2), Zq::new(4), Zq::new(6), Zq::new(8)]
        );
    }

    // Test division by monomials
    #[test]
    fn test_div_by_monomial() {
        let poly = Poly::<4>::create_poly(vec![1, 2, 3, 4]);
        let result = poly.div_by_monomial(2);
        assert_eq!(
            result.coeffs,
            [Zq::new(3), Zq::new(4), Zq::new(1), Zq::new(2)]
        );
    }

    // Test polynomial evaluation
    #[test]
    fn test_eval() {
        let poly = Poly::<4>::create_poly(vec![1, 2, 3, 4]);
        let result = poly.eval(Zq::new(2));
        assert_eq!(result, Zq::new(49));
    }

    // Test equality check
    #[test]
    fn test_is_equal() {
        let poly1 = Poly::<4>::create_poly(vec![1, 2, 3, 4]);
        let poly2 = Poly::<4>::create_poly(vec![1, 2, 3, 4]);
        let poly3 = Poly::<4>::create_poly(vec![4, 3, 2, 1]);
        assert!(poly1.is_equal(&poly2));
        assert!(!poly1.is_equal(&poly3));
    }

    // Test zero polynomial check
    #[test]
    fn test_is_zero_poly() {
        let zero_poly = Poly::<4>::create_poly(vec![0, 0, 0, 0]);
        let non_zero_poly = Poly::<4>::create_poly(vec![1, 0, 0, 0]);
        assert!(zero_poly.is_zero());
        assert!(!non_zero_poly.is_zero());
    }
}
