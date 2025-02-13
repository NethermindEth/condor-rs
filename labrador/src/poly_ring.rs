// This file is part of the polynomial ring operations module.
//
// This module provides implementations for various operations
// in the polynomial ring R = Z_q[X] / (X^d + 1).
//
// Currently implemented functions include:
// - Polynomial addition:         add()
// - Polynomial multiplication:   mul()
// - inner_product/ Dot product:  inner_product()
// - Polynomial substraction:     sub()
// - Polynomial negation:         neg()
// - Scalar multiplication:       scalar_mul()
// - Division by monomials:       div_by_monomial()
// - Polynomial evaluation:       eval()
// - Zero check:                  is_zero()
// - Polynomial equality check:   is_equal()
//
// Further operations and optimizations will be added in future versions.

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Poly<const D: usize> {
    coeffs: [u32; D],
}

impl<const D: usize> Poly<D> {
    // Constructor for the polynomial ring
    pub fn new(coeffs: [u32; D]) -> Self {
        Poly { coeffs }
    }

    // Polynomial addition
    pub fn add(&self, other: &Self) -> Self {
        let mut result = [0u32; D];
        for (r, (a, b)) in result
            .iter_mut()
            .zip(self.coeffs.iter().zip(other.coeffs.iter()))
        {
            *r = a.wrapping_add(*b);
        }
        Poly::new(result)
    }

    // Polynomial substraction
    pub fn sub(&self, other: &Self) -> Self {
        let mut result = [0u32; D];
        for (r, (a, b)) in result
            .iter_mut()
            .zip(self.coeffs.iter().zip(other.coeffs.iter()))
        {
            *r = a.wrapping_sub(*b);
        }
        Poly::new(result)
    }

    // Create a Polynomial from a vector
    pub fn create_poly(coeffs: Vec<i32>) -> Poly<D> {
        let mut arr = [0u32; D];
        let u32_coeffs: Vec<u32> = coeffs.iter().map(|&coeff| coeff as u32).collect();
        // First D elements are assigned directly
        for (i, &u32_coeffs) in u32_coeffs.iter().take(D).enumerate() {
            arr[i] = u32_coeffs;
        }

        // Handle additional elements by subtracting them at (index % D)
        for (i, &u32_coeffs) in u32_coeffs.iter().skip(D).enumerate() {
            let mod_index = i % D;
            arr[mod_index] = arr[mod_index].wrapping_sub(u32_coeffs);
        }

        Poly::new(arr)
    }

    // Polynomial multiplication modulo x^D + 1
    pub fn mul(&self, other: &Self) -> Self {
        let mut result = [0u32; D];

        for i in 0..D {
            for j in 0..D {
                let degree = (i + j) % D;
                if (i + j) > D {
                    result[degree] =
                        result[degree].wrapping_sub(self.coeffs[i].wrapping_mul(other.coeffs[j]));
                } else {
                    // normal multiplication
                    result[degree] =
                        result[degree].wrapping_add(self.coeffs[i].wrapping_mul(other.coeffs[j]));
                }
            }
        }

        Poly::new(result)
    }

    // Dot product between coefficients
    pub fn inner_product(&self, other: &Self) -> u32 {
        self.coeffs
            .iter()
            .zip(other.coeffs.iter())
            .map(|(&a, &b)| a.wrapping_mul(b))
            .fold(0u32, u32::wrapping_add)
    }

    // Polynomial negation
    pub fn neg(&self) -> Self {
        let mut result = [0u32; D];
        for (i, &coeff) in self.coeffs.iter().enumerate() {
            result[i] = coeff.wrapping_neg();
        }
        Poly::new(result)
    }

    // Scalar multiplication
    pub fn scalar_mul(&self, s: u32) -> Self {
        let mut result = [0u32; D];
        for (i, &coeff) in self.coeffs.iter().enumerate() {
            result[i] = s.wrapping_mul(coeff);
        }
        Poly::new(result)
    }

    // (Division by monomials) X^k | performing a cyclic right shift
    pub fn div_by_monomial(&self, k: usize) -> Self {
        let mut result = [0u32; D];
        for (i, &coeff) in self.coeffs.iter().enumerate() {
            // Calculate the new index with wrap-around (k-shift)
            let new_index = (i + k) % D;
            result[new_index] = coeff;
        }
        Poly::new(result)
    }

    // Evaluate the polynomial at a specific point
    pub fn eval(&self, x: u32) -> u32 {
        let mut result = 0u32;
        let mut power = 1u32;
        for &coeff in &self.coeffs {
            result = result.wrapping_add(coeff.wrapping_mul(power));
            power = power.wrapping_mul(x);
        }

        result
    }

    // Check if Polynomial == 0
    pub fn is_zero(&self) -> bool {
        self.coeffs.iter().all(|&coeff| coeff == 0)
    }

    // Check if two polynomials are equal
    pub fn is_equal(&self, other: &Self) -> bool {
        self.coeffs == other.coeffs
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Test constructor and polynomial creation
    #[test]
    fn test_new_and_create_poly() {
        let poly = Poly::new([1, 2, 3, 4]);
        assert_eq!(poly.coeffs, [1, 2, 3, 4]);

        // vector positive coefficients
        let poly_from_vec = Poly::<4>::create_poly(vec![1, 2, 3, 4]);
        assert_eq!(poly_from_vec.coeffs, [1, 2, 3, 4]);

        // negative coefficients
        let poly_from_vec_negative = Poly::<4>::create_poly(vec![-1, -2, -3, -4]);
        let u32_max = u32::MAX;
        assert_eq!(
            poly_from_vec_negative.coeffs,
            [u32_max, u32_max - 1, u32_max - 2, u32_max - 3]
        );
    }

    // Test addition of polynomials
    #[test]
    fn test_add() {
        // within bounds
        let poly1 = Poly::new([1, 2, 3, 4]);
        let poly2 = Poly::new([4, 3, 2, 1]);
        let result = poly1.add(&poly2);
        assert_eq!(result.coeffs, [5, 5, 5, 5]);

        // Outside of bounds
        let poly3 = Poly::new([1, 2, 3, 4]);
        let poly4 = Poly::new([u32::MAX, 3, u32::MAX, 1]);
        let result2 = poly3.add(&poly4);
        assert_eq!(result2.coeffs, [0, 5, 2, 5]);
    }

    // Test subtraction of polynomials
    #[test]
    fn test_sub() {
        let poly1 = Poly::new([5, 10, 15, 20]);
        let poly2 = Poly::new([2, 4, 6, 8]);
        let result = poly1.sub(&poly2);
        assert_eq!(result.coeffs, [3, 6, 9, 12]);

        let poly3 = Poly::new([1, 1, 3, 2]);
        let poly4 = Poly::new([2, 4, 6, 8]);
        let result2 = poly3.sub(&poly4);
        assert_eq!(
            result2.coeffs,
            [u32::MAX, u32::MAX - 2, u32::MAX - 2, u32::MAX - 5]
        );
    }

    // Test negation of polynomial
    #[test]
    fn test_neg() {
        let poly = Poly::new([1, 2, 3, 4]);
        let result = poly.neg();
        assert_eq!(
            result.coeffs,
            [u32::MAX, u32::MAX - 1, u32::MAX - 2, u32::MAX - 3]
        );
    }

    // Test scalar multiplication
    #[test]
    fn test_scalar_mul() {
        let poly = Poly::new([1, 2, 3, 4]);
        let result = poly.scalar_mul(2);
        assert_eq!(result.coeffs, [2, 4, 6, 8]);
    }

    // Test division by monomials
    #[test]
    fn test_div_by_monomial() {
        let poly = Poly::new([1, 2, 3, 4]);
        let result = poly.div_by_monomial(2);
        assert_eq!(result.coeffs, [3, 4, 1, 2]);
    }

    // Test polynomial evaluation
    #[test]
    fn test_eval() {
        let poly = Poly::new([1, 2, 3, 4]);
        let result = poly.eval(2);
        assert_eq!(result, 49);
    }

    // Test equality check
    #[test]
    fn test_is_equal() {
        let poly1 = Poly::new([1, 2, 3, 4]);
        let poly2 = Poly::new([1, 2, 3, 4]);
        let poly3 = Poly::new([4, 3, 2, 1]);
        assert!(poly1.is_equal(&poly2));
        assert!(!poly1.is_equal(&poly3));
    }

    // Test zero polynomial check
    #[test]
    fn test_is_zero_poly() {
        let zero_poly = Poly::new([0, 0, 0, 0]);
        let non_zero_poly = Poly::new([1, 0, 0, 0]);
        assert_eq!(zero_poly.is_zero(), true);
        assert_eq!(non_zero_poly.is_zero(), false);
    }
}
