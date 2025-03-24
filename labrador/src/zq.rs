use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use rand::distr::uniform::{Error, SampleBorrow, SampleUniform, UniformInt, UniformSampler};
use rand::prelude::*;
use std::fmt;
/// Represents an element in the ring Z/qZ where q = 2^32.
/// Uses native u32 operations with automatic modulo reduction through wrapping arithmetic.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Eq, Ord, Default)]
pub struct Zq {
    /// Stored value is always in [0, q-1] due to u32's wrapping behavior
    value: u32,
}

impl Zq {
    /// Modulus q = 2^32 (stored as 0 in u32 due to wrapping behavior)
    pub const Q: u32 = u32::MAX.wrapping_add(1);
    /// Zero element (additive identity)
    pub const ZERO: Self = Self::new(0);
    /// Multiplicative identity
    pub const ONE: Self = Self::new(1);
    /// Maximum element
    pub const MAX: Self = Self::new(u32::MAX);

    /// Creates a new Zq element from a raw u32 value.
    /// No explicit modulo needed as u32 automatically wraps
    pub const fn new(value: u32) -> Self {
        Self { value }
    }

    pub fn to_u128(&self) -> u128 {
        u128::from(self.value)
    }

    pub fn to_i128(&self) -> i128 {
        i128::from(self.value)
    }

    pub fn to_f64(&self) -> f64 {
        let half = Self::MAX.value / 2;
        if self.value > half {
            f64::from(self.value) - f64::from(Self::MAX.value) - 1.0
        } else {
            f64::from(self.value)
        }
    }

    /// Converts an element from Zq (in the range [0, q-1])
    /// to its balanced (signed) representation in i128.
    pub fn to_signed_zq(&self) -> i128 {
        // Compute the threshold for conversion.
        let q: u128 = Self::MAX.to_u128() + 1;
        let half_q = q / 2;
        if self.to_u128() > half_q {
            // If x is above the threshold, subtract q to get a negative value.
            self.to_i128() - i128::try_from(q).unwrap()
        } else {
            self.to_i128()
        }
    }

    /// Converts a balanced (signed) representative `x` in Zq to its standard (unsigned)
    /// representative in the range [0, q-1].
    pub fn to_unsigned_zq(value: i128) -> Zq {
        // If x is negative, add q to bring it into the range [0, q-1]
        let q: i128 = i128::try_from(Self::MAX.to_u128() + 1).unwrap();
        if value < 0 {
            Zq::new(u32::try_from(value + q).unwrap())
        } else {
            Zq::new(u32::try_from(value).unwrap())
        }
    }

    pub const fn is_zero(&self) -> bool {
        self.value == 0
    }
}

// Macro to generate arithmetic trait implementations
macro_rules! impl_arithmetic {
    ($trait:ident, $assign_trait:ident, $method:ident, $assign_method:ident, $op:ident) => {
        impl $trait for Zq {
            type Output = Self;

            fn $method(self, rhs: Self) -> Self::Output {
                Self::new(self.value.$op(rhs.value))
            }
        }

        impl $assign_trait for Zq {
            fn $assign_method(&mut self, rhs: Self) {
                self.value = self.value.$op(rhs.value);
            }
        }
    };
}

impl_arithmetic!(Add, AddAssign, add, add_assign, wrapping_add);
impl_arithmetic!(Sub, SubAssign, sub, sub_assign, wrapping_sub);
impl_arithmetic!(Mul, MulAssign, mul, mul_assign, wrapping_mul);

impl From<u32> for Zq {
    fn from(value: u32) -> Self {
        Self::new(value)
    }
}

impl fmt::Display for Zq {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Shows value with modulus for clarity
        write!(f, "{} (mod 2^32)", self.value)
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
        self.0.sample(rng).into()
    }
}

impl SampleUniform for Zq {
    type Sampler = UniformZq;
}

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
            Zq::MAX + Zq::ONE - self
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
        let a = Zq::MAX;
        let b = Zq::ONE;

        assert_eq!((a + b).value, 0, "u32::MAX + 1 should wrap to 0");
        assert_eq!((b - a).value, 2, "1 - u32::MAX should wrap to 2 (mod 2^32)");
    }

    #[test]
    fn test_subtraction_edge_cases() {
        let max = Zq::MAX;
        let one = Zq::ONE;
        let two = Zq::new(2);

        assert_eq!((one - max).value, 2);
        assert_eq!((two - max).value, 3);
        assert_eq!((max - max).value, 0);
    }

    #[test]
    fn test_multiplication_wrapping() {
        let a = Zq::new(1 << 31);
        let two = Zq::new(2);

        // Multiplication wraps when exceeding u32 range
        assert_eq!((a * two).value, 0, "2^31 * 2 should wrap to 0");
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
        let a: Zq = 5_u32.into();
        assert_eq!(a.value, 5, "Conversion from u32 should preserve value");
    }

    #[test]
    fn test_negative_arithmetic() {
        let small = Zq::new(3);
        let large = Zq::new(5);

        // Test underflow handling (3 - 5 in u32 terms)
        let result = small - large;
        assert_eq!(result.value, u32::MAX - 1, "3 - 5 should wrap to 2^32 - 2");

        // Test compound negative operations
        let mut x = Zq::new(10);
        x -= Zq::new(15);
        assert_eq!(x.value, u32::MAX - 4, "10 -= 15 should wrap to 2^32 - 5");

        // Test negative equivalent value in multiplication
        let a = Zq::MAX; // Represents -1 in mod 2^32 arithmetic
        let b = Zq::new(2);
        assert_eq!(
            (a * b).value,
            u32::MAX - 1,
            "(-1) * 2 should be -2 ≡ 2^32 - 2"
        );
    }

    #[test]
    fn test_display_implementation() {
        let a = Zq::new(5);
        let max = Zq::MAX;

        assert_eq!(format!("{}", a), "5 (mod 2^32)");
        assert_eq!(format!("{}", max), "4294967295 (mod 2^32)");
    }

    #[test]
    fn test_maximum_element() {
        assert_eq!(Zq::MAX, Zq::ZERO - Zq::ONE);
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
    fn test_to_signed_zq() {
        let a = Zq::MAX - Zq::new(99);
        let b = Zq::ONE;

        let c: i128 = a.to_signed_zq();
        let d: i128 = b.to_signed_zq();

        assert_eq!(c, -100);
        assert_eq!(d, 1);
    }
}
