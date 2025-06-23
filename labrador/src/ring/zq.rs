use crate::ring::Norms;
use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use rand::distr::uniform::{Error, SampleBorrow, SampleUniform, UniformInt, UniformSampler};
use rand::prelude::*;
use std::fmt;
use std::iter::Sum;

/// Element of the group **Z/(2^32 − 1)**.
/// Uses native u32 operations with automatic modulo reduction through wrapping arithmetic.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Eq, Ord, Default)]
pub struct Zq {
    /// Values in the range `0..u32::MAX−1`.
    value: u32,
}

impl Zq {
    /// Modulus `q = 2^32 − 1`
    #[allow(clippy::as_conversions)]
    pub const Q: u32 = u32::MAX;

    // ------- constants -------
    pub const ZERO: Self = Self::new(0);
    pub const ONE: Self = Self::new(1);
    pub const TWO: Self = Self::new(2);
    // -1 or Maximum possible value. Equals `q - 1` or ` 2^32 − 2`
    pub const NEG_ONE: Self = Self::new(u32::MAX - 1);

    /// Creates a new Zq element from a raw u32 value.
    /// No explicit modulo needed as u32 automatically wraps
    pub const fn new(value: u32) -> Self {
        debug_assert!(value < Self::Q, "value not reduced modulo q");
        Self { value }
    }

    pub fn to_u128(&self) -> u128 {
        u128::from(self.value)
    }

    pub fn get_value(&self) -> u32 {
        self.value
    }

    pub const fn is_zero(&self) -> bool {
        self.value == 0
    }

    /// Returns `1` iff the element is in `(q-1/2, q)`
    #[allow(clippy::as_conversions)]
    pub fn is_larger_than_half(&self) -> bool {
        self.value > (Self::Q - 1) / 2
    }

    /// Centered representative in `(-q/2, q/2]`.
    #[allow(clippy::as_conversions)]
    pub(crate) fn centered_mod(&self) -> i128 {
        let bound = Self::Q as i128;
        let value = self.value as i128;

        if value > (bound - 1) / 2 {
            value - bound
        } else {
            value
        }
    }

    /// Floor division by another `Zq` value (*not* a field inverse!, just dividing the values).
    pub(crate) fn div_floor_by(&self, rhs: u32) -> Self {
        assert_ne!(rhs, 0, "division by zero");
        Self::new(self.value / rhs)
    }

    /// Decompose the element to #num_parts number of parts,
    /// where each part's infinity norm is less than or equal to bound/2
    pub(crate) fn decompose(&self, bound: Self, num_parts: usize) -> Vec<Zq> {
        assert!(bound >= Self::TWO, "base must be ≥ 2");
        assert_ne!(num_parts, 0, "num_parts cannot be zero");

        let mut parts = vec![Self::ZERO; num_parts];
        let half_bound = bound.div_floor_by(2);
        let mut abs_self = match self.is_larger_than_half() {
            true => -(*self),
            false => *self,
        };

        for part in &mut parts {
            let mut remainder = Self::new(abs_self.value % bound.value);
            if remainder > half_bound {
                remainder -= bound;
            }
            *part = match self.is_larger_than_half() {
                true => -remainder,
                false => remainder,
            };
            abs_self = Self::new((abs_self - remainder).value / bound.value);
            if abs_self == Self::ZERO {
                break;
            }
        }
        parts
    }

    #[allow(clippy::as_conversions)]
    fn add_op(self, rhs: Zq) -> Zq {
        let sum = (self.value as u64 + rhs.value as u64) % Zq::Q as u64;
        Zq::new(sum as u32)
    }

    #[allow(clippy::as_conversions)]
    fn sub_op(self, rhs: Zq) -> Zq {
        let sub = (self.value as u64 + Zq::Q as u64 - rhs.value as u64) % Zq::Q as u64;
        Zq::new(sub as u32)
    }

    #[allow(clippy::as_conversions)]
    fn mul_op(self, b: Zq) -> Zq {
        let prod = (self.value as u64 * b.value as u64) % Zq::Q as u64;
        Zq::new(prod as u32)
    }
}

// Macro to generate arithmetic trait implementations
macro_rules! impl_arithmetic {
    ($trait:ident, $assign_trait:ident, $method:ident, $assign_method:ident, $op:ident) => {
        impl $trait for Zq {
            type Output = Self;

            fn $method(self, rhs: Self) -> Self::Output {
                self.$op(rhs)
            }
        }

        impl $assign_trait for Zq {
            fn $assign_method(&mut self, rhs: Self) {
                *self = self.$op(rhs);
            }
        }

        impl $trait<Zq> for &Zq {
            type Output = Zq;

            fn $method(self, rhs: Zq) -> Self::Output {
                self.$op(rhs)
            }
        }

        impl $trait<&Zq> for &Zq {
            type Output = Zq;

            fn $method(self, rhs: &Zq) -> Self::Output {
                self.$op(*rhs)
            }
        }
    };
}

impl_arithmetic!(Add, AddAssign, add, add_assign, add_op);
impl_arithmetic!(Sub, SubAssign, sub, sub_assign, sub_op);
impl_arithmetic!(Mul, MulAssign, mul, mul_assign, mul_op);

// Implement the Neg trait for Zq.
impl Neg for Zq {
    type Output = Zq;

    /// Returns the additive inverse of the field element.
    ///
    /// Wrap around (q - a) mod q.
    fn neg(self) -> Zq {
        // If the value is zero, its inverse is itself.
        if self.value == 0 {
            self
        } else {
            #[allow(clippy::as_conversions)]
            Zq::new(Zq::Q - self.get_value())
        }
    }
}

impl fmt::Display for Zq {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Shows value with modulus for clarity
        write!(f, "{} (mod {})", self.value, Zq::Q)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct UniformZq(UniformInt<u32>);

impl UniformSampler for UniformZq {
    type X = Zq;

    fn new<B1, B2>(low: B1, high: B2) -> Result<Self, Error>
    where
        B1: SampleBorrow<Self::X> + Sized,
        B2: SampleBorrow<Self::X> + Sized,
    {
        UniformInt::<u32>::new(low.borrow().value, high.borrow().value).map(UniformZq)
    }
    fn new_inclusive<B1, B2>(low: B1, high: B2) -> Result<Self, Error>
    where
        B1: SampleBorrow<Self::X> + Sized,
        B2: SampleBorrow<Self::X> + Sized,
    {
        UniformInt::<u32>::new_inclusive(low.borrow().value, high.borrow().value).map(UniformZq)
    }
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Self::X {
        Self::X::new(self.0.sample(rng))
    }
}

impl SampleUniform for Zq {
    type Sampler = UniformZq;
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

/// Adds `rhs` into `lhs` component‑wise.
pub fn add_assign_two_zq_vectors(lhs: &mut [Zq], rhs: Vec<Zq>) {
    debug_assert_eq!(lhs.len(), rhs.len(), "vector length mismatch");
    lhs.iter_mut().zip(rhs).for_each(|(l, r)| *l += r);
}

// Implement l2 and infinity norms for a slice of Zq elements
impl Norms for [Zq] {
    type NormType = u128;

    #[allow(clippy::as_conversions)]
    fn l2_norm_squared(&self) -> Self::NormType {
        self.iter().fold(0u128, |acc, coeff| {
            let c = coeff.centered_mod();
            acc + (c * c) as u128
        })
    }

    #[allow(clippy::as_conversions)]
    fn linf_norm(&self) -> Self::NormType {
        self.iter()
            .map(|coeff| coeff.centered_mod().unsigned_abs())
            .max()
            .unwrap_or(0)
    }
}

#[cfg(test)]
mod norm_tests {
    use super::*;

    #[test]
    fn test_l2_norm() {
        let zq_vector = [
            Zq::new(1),
            Zq::new(2),
            Zq::new(3),
            Zq::new(4),
            Zq::new(5),
            Zq::new(6),
            Zq::new(7),
        ];
        let res = zq_vector.l2_norm_squared();

        assert_eq!(res, 140);
    }

    #[test]
    fn test_l2_norm_with_negative_values() {
        let zq_vector = [
            Zq::new(1),
            Zq::new(2),
            Zq::new(3),
            -Zq::new(4),
            -Zq::new(5),
            -Zq::new(6),
            -Zq::new(7),
        ];
        let res = zq_vector.l2_norm_squared();

        assert_eq!(res, 140);
    }

    #[test]
    fn test_linf_norm() {
        let zq_vector = [
            Zq::new(1),
            Zq::new(200),
            Zq::new(300),
            Zq::new(40),
            -Zq::new(5),
            -Zq::new(6),
            -Zq::new(700000),
        ];
        let res = zq_vector.linf_norm();
        assert_eq!(res, 700000);

        let zq_vector = [
            Zq::new(1000000),
            Zq::new(200),
            Zq::new(300),
            Zq::new(40),
            -Zq::new(5),
            -Zq::new(6),
            -Zq::new(999999),
        ];
        let res = zq_vector.linf_norm();
        assert_eq!(res, 1000000);

        let zq_vector = [
            Zq::new(1),
            Zq::new(2),
            Zq::new(3),
            -Zq::new(4),
            Zq::new(0),
            -Zq::new(3),
            -Zq::new(2),
            -Zq::new(1),
        ];
        let res = zq_vector.linf_norm();
        assert_eq!(res, 4);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_to_u128() {
        let a = Zq::new(10);
        let b = a.to_u128();
        assert_eq!(b, 10u128);
    }

    #[test]
    fn test_is_zero() {
        let a = Zq::new(0);
        let b = Zq::new(10);
        assert!(a.is_zero());
        assert!(!b.is_zero());
    }

    #[test]
    fn test_get_value() {
        let a = Zq::new(1000);
        assert_eq!(a.get_value(), 1000u32);
    }

    #[test]
    fn test_basic_arithmetic() {
        let a = Zq::new(5);
        let b = Zq::new(3);

        // Addition
        assert_eq!((a + b).value, 8, "5 + 3 should be 8");
        // Subtraction
        assert_eq!((a - b).value, 2, "5 - 3 should be 2");
        // Multiplication
        assert_eq!((a * b).value, 15, "5 * 3 should be 15");
    }

    #[test]
    fn test_wrapping_arithmetic() {
        let a = Zq::NEG_ONE;
        let b = Zq::ONE;

        assert_eq!((a + b).value, 0, "u32::MAX + 1 should wrap to 0");
        assert_eq!((b - a).value, 2, "1 - u32::MAX should wrap to 2 (mod 2^32)");
    }

    #[test]
    fn test_subtraction_edge_cases() {
        let max = Zq::NEG_ONE;
        let one = Zq::ONE;
        let two = Zq::TWO;

        assert_eq!((one - max).value, 2);
        assert_eq!((two - max).value, 3);
        assert_eq!((max - max).value, 0);
    }

    #[test]
    fn test_multiplication_wrapping() {
        let a = Zq::new(1 << 31);
        let two = Zq::TWO;

        // Multiplication wraps when exceeding u32 range
        assert_eq!((a * two).value, 1, "2^31 * 2 should wrap to 1");
    }

    #[test]
    fn test_assignment_operators() {
        let mut a = Zq::new(5);
        let b = Zq::new(3);

        a += b;
        assert_eq!(a.value, 8, "5 += 3 should be 8");

        a -= b;
        assert_eq!(a.value, 5, "8 -= 3 should be 5");

        a *= b;
        assert_eq!(a.value, 15, "5 *= 3 should be 15");
    }

    #[test]
    fn test_conversion_from_u32() {
        let a: Zq = Zq::new(5);
        assert_eq!(a.value, 5, "Conversion from u32 should preserve value");
    }

    #[test]
    fn test_negative_arithmetic() {
        let small = Zq::new(3);
        let large = Zq::new(5);

        // Test underflow handling (3 - 5 in u32 terms)
        let result = small - large;
        assert_eq!(result.value, u32::MAX - 2, "3 - 5 should wrap to 2^32 - 2");

        // Test compound negative operations
        let mut x = Zq::new(10);
        x -= Zq::new(15);
        assert_eq!(x.value, u32::MAX - 5, "10 -= 15 should wrap to 2^32 - 5");

        // Test negative equivalent value in multiplication
        let a = Zq::NEG_ONE; // Represents -1 in mod 2^32 arithmetic
        let b = Zq::TWO;
        assert_eq!(
            (a * b).value,
            u32::MAX - 2,
            "(-1) * 2 should be -2 ≡ 2^32 - 2"
        );
    }

    #[test]
    fn test_display_implementation() {
        let a = Zq::new(5);
        let max = Zq::NEG_ONE;
        assert_eq!(format!("{}", a), format!("5 (mod {})", Zq::Q));
        assert_eq!(format!("{}", max), format!("4294967294 (mod {})", Zq::Q));
    }

    #[test]
    fn test_maximum_element() {
        dbg!(Zq::NEG_ONE);
        dbg!(Zq::ZERO);
        dbg!(Zq::ONE);
        dbg!(Zq::ZERO - Zq::ONE);
        assert_eq!(Zq::NEG_ONE, Zq::ZERO - Zq::ONE);
    }

    #[test]
    fn test_ord() {
        let a = Zq::new(100);
        let b = Zq::new(200);
        let c = Zq::new(100);
        let d = Zq::new(400);

        let res_1 = a.cmp(&b);
        let res_2 = a.cmp(&c);
        let res_3 = d.cmp(&b);
        assert!(res_1.is_lt());
        assert!(res_2.is_eq());
        assert!(res_3.is_gt());
        assert_eq!(a, c);
        assert!(a < b);
        assert!(d > b);
    }

    #[test]
    fn test_neg() {
        let a = Zq::new(100);
        let b = Zq::ZERO;
        let neg_a: Zq = -a;
        let neg_b: Zq = -b;

        assert_eq!(neg_a + a, Zq::ZERO);
        assert_eq!(neg_b, Zq::ZERO);
    }

    #[test]
    fn test_centered_mod() {
        let a = -Zq::new(1);
        assert_eq!(-1, a.centered_mod());

        let a = Zq::new(4294967103);
        assert_eq!(a, -Zq::new(192));
        assert_eq!(-192, a.centered_mod());
    }
}

#[cfg(test)]
mod test_decomposition {
    use crate::ring::{zq::Zq, Norms};

    #[test]
    fn test_zq_decomposition() {
        let (base, parts) = (Zq::new(12), 10);
        let pos_zq = Zq::new(29);
        let neg_zq = -Zq::new(29);

        let pos_decomposed = pos_zq.decompose(base, parts);
        let neg_decomposed = neg_zq.decompose(base, parts);

        assert_eq!(
            pos_decomposed,
            vec![
                Zq::new(5),
                Zq::new(2),
                Zq::ZERO,
                Zq::ZERO,
                Zq::ZERO,
                Zq::ZERO,
                Zq::ZERO,
                Zq::ZERO,
                Zq::ZERO,
                Zq::ZERO
            ]
        );
        assert_eq!(
            neg_decomposed,
            vec![
                -Zq::new(5),
                -Zq::new(2),
                Zq::ZERO,
                Zq::ZERO,
                Zq::ZERO,
                Zq::ZERO,
                Zq::ZERO,
                Zq::ZERO,
                Zq::ZERO,
                Zq::ZERO
            ]
        );
    }

    #[test]
    fn test_zq_recompositoin() {
        let (base, parts) = (Zq::new(1802), 10);
        let pos_zq = -Zq::new(16200);

        let pos_decomposed = pos_zq.decompose(base, parts);
        let mut exponensial_base = Zq::new(1);
        let mut result = Zq::new(0);
        for decomposed_part in pos_decomposed {
            result += decomposed_part * exponensial_base;
            exponensial_base *= base;
        }
        assert_eq!(result, pos_zq)
    }

    #[test]
    fn test_zq_recompositoin_positive() {
        let (base, parts) = (Zq::new(1802), 10);
        let pos_zq = Zq::new(23071);

        let pos_decomposed = pos_zq.decompose(base, parts);
        let mut exponensial_base = Zq::new(1);
        let mut result = Zq::new(0);
        for decomposed_part in pos_decomposed {
            result += decomposed_part * exponensial_base;
            exponensial_base *= base;
        }
        assert_eq!(result, pos_zq)
    }

    #[test]
    fn test_linf_norm() {
        let (base, parts) = (Zq::new(1802), 10);
        let pos_zq = Zq::new(16200);

        let pos_decomposed = pos_zq.decompose(base, parts);
        dbg!(&pos_decomposed);
        assert!(pos_decomposed.linf_norm() <= 901);
    }
}
