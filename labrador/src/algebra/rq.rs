use crate::algebra::zq::Zq;
use crate::algebra::{Ring, RingOps};
use std::cmp::Ordering;
use std::fmt::{Debug, Display, Formatter};
use std::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};

/// A quotient polynomial ring `R_q = Z_q[X] / (X^d + 1)` - i.e. polynomials with coefficients in `Z_q` modulo
/// `X^d + 1` (so `X^d + 1` is equivalent to zero).
/// In this ring, the degree of the polynomial is at most `d - 1`.
#[derive(Clone, PartialEq, Eq)]
struct Rq<const D: usize> {
    /// Coefficients of the polynomial, from the constant term to the highest degree term.
    pub coeffs: [Zq; D],
}

impl<const D: usize> Rq<D> {
    /// Creates a new polynomial with the given coefficients.
    #[allow(dead_code)]
    pub fn new(coeffs: &[Zq]) -> Self {
        let mut raw_coeffs: [Zq; D] = [Zq::ZERO; D];
        for (i, coeff) in coeffs.iter().enumerate() {
            raw_coeffs[i] = *coeff;
        }
        Self { coeffs: raw_coeffs }
    }
}

impl<const D: usize> Display for Rq<D> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "Rq({:?})", self.coeffs)
    }
}

impl<const D: usize> Add for &Rq<D> {
    type Output = Rq<D>;

    fn add(self, rhs: Self) -> Self::Output {
        // iterator_try_collect is unstable, so use plain old for loop
        let mut coeffs = [Zq::ZERO; D];
        for (i, coeff) in coeffs.iter_mut().enumerate().take(D) {
            *coeff = self.coeffs[i] + rhs.coeffs[i];
        }
        Rq { coeffs }
    }
}

impl<const D: usize> Sub for &Rq<D> {
    type Output = Rq<D>;

    fn sub(self, rhs: Self) -> Self::Output {
        // iterator_try_collect is unstable, so use plain old for loop
        let mut coeffs = [Zq::ZERO; D];
        for (i, coeff) in coeffs.iter_mut().enumerate().take(D) {
            *coeff = self.coeffs[i] - rhs.coeffs[i];
        }
        Rq { coeffs }
    }
}

impl<const D: usize> Mul for &Rq<D> {
    type Output = Rq<D>;

    // TODO: Use NTT
    fn mul(self, rhs: Self) -> Self::Output {
        // Initialize a vector to hold the intermediate multiplication result
        let mut coeffs = vec![Zq::ZERO; 2 * D - 1];
        for (i, &coeff1) in self.coeffs.iter().enumerate() {
            for (j, &coeff2) in rhs.coeffs.iter().enumerate() {
                coeffs[i + j] += coeff1 * coeff2;
            }
        }

        // Reduce modulo X^d + 1
        if coeffs.len() > D {
            let modulus_minus_one = Zq::new(Zq::Q - 1);
            let (front, back) = coeffs.split_at_mut(D);
            for (i, &overflow) in back.iter().enumerate() {
                front[i] += overflow * modulus_minus_one;
            }
            coeffs.truncate(D);
        }
        Rq {
            coeffs: coeffs.try_into().unwrap(),
        }
    }
}

impl<const D: usize> AddAssign for Rq<D> {
    fn add_assign(&mut self, rhs: Self) {
        for i in 0..D {
            self.coeffs[i] += rhs.coeffs[i];
        }
    }
}

impl<const D: usize> SubAssign for Rq<D> {
    fn sub_assign(&mut self, rhs: Self) {
        for i in 0..D {
            self.coeffs[i] -= rhs.coeffs[i];
        }
    }
}

impl<const D: usize> MulAssign for Rq<D> {
    fn mul_assign(&mut self, rhs: Self) {
        for i in 0..D {
            self.coeffs[i] *= rhs.coeffs[i];
        }
    }
}

impl<const D: usize> Ord for Rq<D> {
    fn cmp(&self, other: &Self) -> Ordering {
        for i in (0..D).rev() {
            match self.coeffs[i].cmp(&other.coeffs[i]) {
                Ordering::Equal => continue,
                non_equal => return non_equal,
            }
        }
        Ordering::Equal
    }
}

impl<const D: usize> PartialOrd for Rq<D> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<const D: usize> Debug for Rq<D> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        // No need for special formatting, just use Display
        <Self as Display>::fmt(self, f)
    }
}

impl<const D: usize> RingOps<Rq<D>> for &Rq<D> {}

impl<const D: usize> Ring for Rq<D> {
    const ZERO: Self = Rq {
        coeffs: [Zq::ZERO; D],
    };
    const ONE: Self = {
        let mut coeffs = [Zq::ZERO; D];
        coeffs[0] = Zq::ONE;
        Rq { coeffs }
    };
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rq_mul() {
        let poly1 = Rq::new(&[Zq::new(1), Zq::new(2)]);
        let poly2 = Rq::new(&[Zq::new(3), Zq::new(4)]);
        let result = &poly1 * &poly2;
        assert_eq!(result.coeffs, [Zq::new(3), Zq::new(10), Zq::new(8)]); // 1*3, 1*4 + 2*3, 2*4
    }

    #[test]
    fn test_rq_mul_overflow() {
        // Create two polynomials that will cause overflow when multiplied
        // For example, (X^63 + 1) * (X^63 + x) = X^126 + X^64 + X^63 + X
        // Modulo X^64 + 1, X^64 = -1, so X^126 = X^(2*64 -2) = X^-2 = X^62, X*X^63 = -1
        // Thus, X^126 + X^64 + X^63 + X mod X^64+1 = (-1)*X^62 + (-1) + X + X^63 = - 1 + X - X^62 + X^63

        // Initialize poly1 as X^63 + 1
        let mut poly1_coeffs = [Zq::ZERO; 64];
        poly1_coeffs[0] = Zq::ONE; // Constant term
        poly1_coeffs[63] = Zq::ONE; // X^63 term
        let poly1 = Rq::new(&poly1_coeffs);

        // Initialize poly1 as X^63 + X
        let mut poly2_coeffs = [Zq::ZERO; 64];
        poly2_coeffs[1] = Zq::ONE; // X term
        poly2_coeffs[63] = Zq::ONE; // X^63 term
        let poly2 = Rq::new(&poly2_coeffs);

        // Multiply poly1 by poly2
        let product = &poly1 * &poly2;

        // Expected coefficients after reduction modulo X^64 + 1:
        // coeffs[0] = 1
        // coeffs[62] = Zq::Q - 1  (since -1 mod q)
        // coeffs[63] = 2
        // All other coefficients should be 0
        let mut expected_coeffs = [Zq::ZERO; 64];
        expected_coeffs[0] = Zq::new(Zq::Q - 1); // Constant term
        expected_coeffs[1] = Zq::ONE; // X term
        expected_coeffs[62] = Zq::new(Zq::Q - 1); // X^62 term
        expected_coeffs[63] = Zq::ONE; // X^63 term

        // Assert that the coefficients match the expected values
        assert_eq!(
            product.coeffs, expected_coeffs,
            "Overflow handling in multiplication is incorrect"
        );
    }

    #[test]
    fn test_rq_add_mul() {
        let a = Rq::<5>::new(&[Zq::new(1), Zq::new(2), Zq::new(3)]);
        let b = Rq::<5>::new(&[Zq::new(4), Zq::new(5), Zq::new(6)]);
        let c = &a + &b;
        assert_eq!(
            c.coeffs,
            [Zq::new(5), Zq::new(7), Zq::new(9), Zq::ZERO, Zq::ZERO]
        );
        let d = &a * &b;
        assert_eq!(
            d.coeffs,
            [
                Zq::new(4),
                Zq::new(13),
                Zq::new(28),
                Zq::new(27),
                Zq::new(18)
            ]
        );
    }
}
