use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use rand::distr::uniform::{Error, SampleBorrow, SampleUniform, UniformInt, UniformSampler};
use rand::prelude::*;
use std::fmt;
use std::iter::Sum;
/// Represents an element in the ring Z/qZ where q = 2^32.
/// Uses native u32 operations with automatic modulo reduction through wrapping arithmetic.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Eq, Ord, Default)]
pub struct Zq {
    /// Stored value is always in [0, q-1] due to u32's wrapping behavior
    value: u32,
}

impl Zq {
    /// Modulus q = 2^32 - 1
    #[allow(clippy::as_conversions)]
    pub const Q: u64 = u32::MAX as u64;
    /// Zero element (additive identity)
    pub const ZERO: Self = Self::new(0);
    /// Multiplicative identity
    pub const ONE: Self = Self::new(1);
    /// Two
    pub const TWO: Self = Self::new(2);
    /// Maximum element
    pub const NEG_ONE: Self = Self::new(u32::MAX - 1);

    /// Creates a new Zq element from a raw u32 value.
    /// No explicit modulo needed as u32 automatically wraps
    pub const fn new(value: u32) -> Self {
        Self { value }
    }

    pub fn to_u128(&self) -> u128 {
        u128::from(self.value)
    }

    pub const fn is_zero(&self) -> bool {
        self.value == 0
    }

    pub fn get_value(&self) -> u32 {
        self.value
    }

    /// Returns the centered representative modulo the given bound
    /// Result is guaranteed to be in (-bound/2, bound/2]
    ///
    /// # Panics
    ///
    /// Panics if `bound` is zero.
    pub(crate) fn centered_mod(&self, bound: Self) -> Self {
        assert!(
            bound != Zq::ZERO,
            "cannot get centered representative modulo for zero bound"
        );
        let bounded_coeff = Self::new(self.value % bound.value);
        let half_bound = bound.scale_by(Self::TWO);

        if bounded_coeff > half_bound {
            bounded_coeff - bound
        } else {
            bounded_coeff
        }
    }

    /// Scales by other Zq.
    ///
    /// Effectively it is a floor division of internal values.
    /// But for the ring of integers there is no defined division
    /// operation.
    ///
    /// # Panics
    ///
    /// Panics if `bound` is zero.
    pub(crate) fn scale_by(&self, rhs: Self) -> Self {
        assert!(rhs != Zq::ZERO, "cannot scale by zero");
        Self::new(self.value / rhs.value)
    }

    #[allow(clippy::as_conversions)]
    fn add_op(self, rhs: Zq) -> Zq {
        let sum = (self.value as u64 + rhs.value as u64) % Zq::Q;
        Zq::new(sum as u32)
    }

    #[allow(clippy::as_conversions)]
    fn sub_op(self, rhs: Zq) -> Zq {
        let sub = (self.value as u64 + Zq::Q - rhs.value as u64) % Zq::Q;
        Zq::new(sub as u32)
    }

    #[allow(clippy::as_conversions)]
    fn mul_op(self, b: Zq) -> Zq {
        let prod = (self.value as u64 * b.value as u64) % Zq::Q;
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
            Zq::new(Zq::Q as u32 - self.get_value())
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

pub fn add_assign_two_zq_vectors(first: &mut [Zq], second: Vec<Zq>) {
    first
        .iter_mut()
        .zip(second)
        .for_each(|(first_coeff, second_coeff)| *first_coeff += second_coeff);
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
}
