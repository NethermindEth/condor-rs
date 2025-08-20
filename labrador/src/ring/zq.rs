use crate::ring::Norms;
use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use rand::distr::uniform::{Error, SampleBorrow, SampleUniform, UniformInt, UniformSampler};
use rand::prelude::*;
use std::fmt;
use std::iter::Sum;
use std::marker::PhantomData;

/// Trait defining modular arithmetic parameters and operations
pub trait Mod: 'static + Copy + Clone + Send + Sync + PartialEq + Eq + PartialOrd + Ord {
    const MODULUS: u64;
    const INV: u64; // multiplicative inverse for Montgomery reduction, etc.
    const BITS: u32; // number of bits needed for this modulus
    const ROOT_OF_UNITY: Option<u64> = None; // for NTT operations
    const IS_PRIME: bool = true;
    
    /// Optimized modular reduction - can be overridden for specific moduli
    #[inline]
    fn reduce(value: u128) -> u64 {
        (value % Self::MODULUS as u128) as u64
    }
    
    /// Optimized modular addition
    #[inline]
    fn add_mod(a: u64, b: u64) -> u64 {
        let sum = a as u128 + b as u128;
        Self::reduce(sum)
    }
    
    /// Optimized modular subtraction
    #[inline]
    fn sub_mod(a: u64, b: u64) -> u64 {
        let diff = (a as u128 + Self::MODULUS as u128 - b as u128);
        Self::reduce(diff)
    }
    
    /// Optimized modular multiplication
    #[inline]
    fn mul_mod(a: u64, b: u64) -> u64 {
        let prod = a as u128 * b as u128;
        Self::reduce(prod)
    }
    
    /// Modular negation
    #[inline]
    fn neg_mod(a: u64) -> u64 {
        if a == 0 { 0 } else { Self::MODULUS - a }
    }
}

/// Element of the group **Z/M::MODULUS** using trait-based modulus definition.
/// Uses native u64 operations with automatic modulo reduction.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Eq, Ord, Default, Hash)]
pub struct Zq<M: Mod> {
    /// Values in the range `0..M::MODULUS`.
    value: u64,
    _phantom: PhantomData<M>,
}

impl<M: Mod> Zq<M> {
    /// Creates a new Zq element from a raw u64 value.
    /// Applies modulo M::MODULUS automatically
    #[inline]
    pub const fn new(value: u64) -> Self {
        Self { 
            value: if M::MODULUS == 0 { value } else { value % M::MODULUS },
            _phantom: PhantomData,
        }
    }
    
    /// Creates a Zq element from a pre-reduced value (unsafe - no modulo check)
    #[inline]
    pub const fn from_raw(value: u64) -> Self {
        Self {
            value,
            _phantom: PhantomData,
        }
    }

    // Constants - now using the trait's MODULUS
    pub const ZERO: Self = Self::from_raw(0);
    pub const ONE: Self = Self::from_raw(1);
    pub const TWO: Self = Self::from_raw(2);
    // -1 or Maximum possible value. Equals `M::MODULUS - 1`
    pub const NEG_ONE: Self = Self::from_raw(M::MODULUS - 1);

    #[inline]
    pub fn to_u128(&self) -> u128 {
        u128::from(self.value)
    }

    #[inline]
    pub fn get_value(&self) -> u64 {
        self.value
    }

    #[inline]
    pub const fn is_zero(&self) -> bool {
        self.value == 0
    }

    /// Returns `true` iff the element is in `(M::MODULUS-1/2, M::MODULUS)`
    #[inline]
    pub fn is_larger_than_half(&self) -> bool {
        self.value > (M::MODULUS - 1) / 2
    }

    /// Centered representative in `(-M::MODULUS/2, M::MODULUS/2]`.
    pub(crate) fn centered_mod(&self) -> i128 {
        let bound = M::MODULUS as i128;
        let value = self.value as i128;

        if value > (bound - 1) / 2 {
            value - bound
        } else {
            value
        }
    }

    /// Floor division by another value (*not* a field inverse!, just dividing the values).
    pub(crate) fn div_floor_by(&self, rhs: u64) -> Self {
        assert_ne!(rhs, 0, "division by zero");
        Self::new(self.value / rhs)
    }

    /// Decompose the element to #num_parts number of parts,
    /// where each part's infinity norm is less than or equal to bound/2
    pub(crate) fn decompose(&self, bound: Self, num_parts: u64) -> Vec<Self> {
        assert!(bound >= Self::TWO, "base must be ≥ 2");
        assert_ne!(num_parts, 0, "num_parts cannot be zero");

        let mut parts =
            vec![Self::ZERO; usize::try_from(num_parts).expect("num_parts does not fit in usize")];
        let half_bound = bound.div_floor_by(2);
        let mut abs_self = if self.is_larger_than_half() {
            -(*self)
        } else {
            *self
        };

        for part in &mut parts {
            let mut remainder = Self::new(abs_self.value % bound.value);
            if remainder > half_bound {
                remainder -= bound;
            }
            *part = if self.is_larger_than_half() {
                -remainder
            } else {
                remainder
            };
            abs_self = Self::new((abs_self - remainder).value / bound.value);
            if abs_self == Self::ZERO {
                break;
            }
        }
        parts
    }

    // Remove the internal arithmetic methods since we'll use the trait methods directly
}

// Arithmetic implementations using the trait methods
impl<M: Mod> Add for Zq<M> {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        Self::from_raw(M::add_mod(self.value, rhs.value))
    }
}

impl<M: Mod> Sub for Zq<M> {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        Self::from_raw(M::sub_mod(self.value, rhs.value))
    }
}

impl<M: Mod> Mul for Zq<M> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        Self::from_raw(M::mul_mod(self.value, rhs.value))
    }
}

impl<M: Mod> Neg for Zq<M> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        Self::from_raw(M::neg_mod(self.value))
    }
}

// Assignment operators
impl<M: Mod> AddAssign for Zq<M> {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl<M: Mod> SubAssign for Zq<M> {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl<M: Mod> MulAssign for Zq<M> {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

// Reference implementations
impl<M: Mod> Add<Zq<M>> for &Zq<M> {
    type Output = Zq<M>;
    
    #[inline]
    fn add(self, rhs: Zq<M>) -> Self::Output {
        *self + rhs
    }
}

impl<M: Mod> Add<&Zq<M>> for &Zq<M> {
    type Output = Zq<M>;
    
    #[inline]
    fn add(self, rhs: &Zq<M>) -> Self::Output {
        *self + *rhs
    }
}

impl<M: Mod> Sub<Zq<M>> for &Zq<M> {
    type Output = Zq<M>;
    
    #[inline]
    fn sub(self, rhs: Zq<M>) -> Self::Output {
        *self - rhs
    }
}

impl<M: Mod> Sub<&Zq<M>> for &Zq<M> {
    type Output = Zq<M>;
    
    #[inline]
    fn sub(self, rhs: &Zq<M>) -> Self::Output {
        *self - *rhs
    }
}

impl<M: Mod> Mul<Zq<M>> for &Zq<M> {
    type Output = Zq<M>;
    
    #[inline]
    fn mul(self, rhs: Zq<M>) -> Self::Output {
        *self * rhs
    }
}

impl<M: Mod> Mul<&Zq<M>> for &Zq<M> {
    type Output = Zq<M>;
    
    #[inline]
    fn mul(self, rhs: &Zq<M>) -> Self::Output {
        *self * *rhs
    }
}

impl<M: Mod> fmt::Display for Zq<M> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} (mod {})", self.value, M::MODULUS)
    }
}

// Random sampling
#[derive(Clone, Copy, Debug)]
pub struct UniformZq<M: Mod>(UniformInt<u64>, PhantomData<M>);

impl<M: Mod> UniformSampler for UniformZq<M> {
    type X = Zq<M>;

    fn new<B1, B2>(low: B1, high: B2) -> Result<Self, Error>
    where
        B1: SampleBorrow<Self::X> + Sized,
        B2: SampleBorrow<Self::X> + Sized,
    {
        UniformInt::<u64>::new(low.borrow().value, high.borrow().value)
            .map(|sampler| UniformZq(sampler, PhantomData))
    }
    
    fn new_inclusive<B1, B2>(low: B1, high: B2) -> Result<Self, Error>
    where
        B1: SampleBorrow<Self::X> + Sized,
        B2: SampleBorrow<Self::X> + Sized,
    {
        UniformInt::<u64>::new_inclusive(low.borrow().value, high.borrow().value)
            .map(|sampler| UniformZq(sampler, PhantomData))
    }
    
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Self::X {
        Self::X::new(self.0.sample(rng))
    }
}

impl<M: Mod> SampleUniform for Zq<M> {
    type Sampler = UniformZq<M>;
}

impl<M: Mod> Sum for Zq<M> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::ZERO, |a, b| a + b)
    }
}

/// Adds `rhs` into `lhs` component‑wise.
pub fn add_assign_two_zq_vectors<M: Mod>(lhs: &mut [Zq<M>], rhs: Vec<Zq<M>>) {
    debug_assert_eq!(lhs.len(), rhs.len(), "vector length mismatch");
    lhs.iter_mut().zip(rhs).for_each(|(l, r)| *l += r);
}

// Implement l2 and infinity norms for a slice of Zq elements
impl<M: Mod> Norms for [Zq<M>] {
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

// Define specific moduli for different use cases

/// The original u32::MAX modulus from your implementation
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct U32MaxMod;
impl Mod for U32MaxMod {
    const MODULUS: u64 = u64::MAX;
    const INV: u64 = 0; // Not used for this modulus
    const BITS: u32 = 32;
    const IS_PRIME: bool = false; // 2^32 - 1 is not prime
}

/// Standard LaBRADOR modulus (CRYSTALS-Dilithium compatible)
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct LabradorMod;
impl Mod for LabradorMod {
    const MODULUS: u64 = 8380417; // 2^23 - 2^13 + 1
    const INV: u64 = 58728449; // Precomputed inverse
    const BITS: u32 = 23;
    const ROOT_OF_UNITY: Option<u64> = Some(1753);
}

/// Falcon-512 single signature modulus (41 bits, but truncated to fit u32)
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Falcon512Mod;
impl Mod for Falcon512Mod {
    const MODULUS: u64 = (1u64 << 31) - 1;
    const INV: u64 = 0;
    const BITS: u32 = 31;
    const IS_PRIME: bool = true; // 2^31 - 1 is a Mersenne prime
}

// Type aliases for convenience and backward compatibility
pub type ZqU32Max = Zq<U32MaxMod>;
pub type ZqLabrador = Zq<LabradorMod>;
pub type ZqFalcon512 = Zq<Falcon512Mod>;

// For your existing code, you can use ZqU32Max as a drop-in replacement
// Or create a simple type alias:
// pub type Zq = ZqU32Max;

#[cfg(test)]
mod tests {
    use super::*;

    // Use the original modulus for backward compatibility
    type TestZq = ZqU32Max;

    #[test]
    fn test_to_u128() {
        let a = TestZq::new(10);
        let b = a.to_u128();
        assert_eq!(b, 10u128);
    }

    #[test]
    fn test_is_zero() {
        let a = TestZq::new(0);
        let b = TestZq::new(10);
        assert!(a.is_zero());
        assert!(!b.is_zero());
    }

    #[test]
    fn test_get_value() {
        let a = TestZq::new(1000);
        assert_eq!(a.get_value(), 1000u64);
    }

    #[test]
    fn test_basic_arithmetic() {
        let a = TestZq::new(5);
        let b = TestZq::new(3);

        // Addition
        assert_eq!((a + b).value, 8, "5 + 3 should be 8");
        // Subtraction
        assert_eq!((a - b).value, 2, "5 - 3 should be 2");
        // Multiplication
        assert_eq!((a * b).value, 15, "5 * 3 should be 15");
    }

    #[test]
    fn test_wrapping_arithmetic() {
        let a = TestZq::NEG_ONE;
        let b = TestZq::ONE;

        assert_eq!((a + b).value, 0, "u32::MAX-1 + 1 should wrap to 0");
        assert_eq!((b - a).value, 2, "1 - (u32::MAX-1) should wrap to 2");
    }

    #[test]
    fn test_different_moduli() {
        let labrador_a = ZqLabrador::new(5);
        let labrador_b = ZqLabrador::new(3);
        let result = labrador_a + labrador_b;
        assert_eq!(result.get_value(), 8);

        // Test that we can't accidentally mix different moduli
        // This would cause a compile-time error:
        // let mixed = labrador_a + TestZq::new(3); // Compile error!
    }

    #[test]
    fn test_assignment_operators() {
        let mut a = TestZq::new(5);
        let b = TestZq::new(3);

        a += b;
        assert_eq!(a.value, 8, "5 += 3 should be 8");

        a -= b;
        assert_eq!(a.value, 5, "8 -= 3 should be 5");

        a *= b;
        assert_eq!(a.value, 15, "5 *= 3 should be 15");
    }

    #[test]
    fn test_neg() {
        let a = TestZq::new(100);
        let b = TestZq::ZERO;
        let neg_a: TestZq = -a;
        let neg_b: TestZq = -b;

        assert_eq!(neg_a + a, TestZq::ZERO);
        assert_eq!(neg_b, TestZq::ZERO);
    }

    #[test]
    fn test_display_implementation() {
        let a = TestZq::new(5);
        let max = TestZq::NEG_ONE;
        assert_eq!(format!("{a}"), "5 (mod 4294967295)");
        assert_eq!(format!("{max}"), "4294967294 (mod 4294967295)");
    }
}

#[cfg(test)]
mod norm_tests {
    use super::*;
    
    type TestZq = ZqU32Max;

    #[test]
    fn test_l2_norm() {
        let zq_vector = [
            TestZq::new(1),
            TestZq::new(2),
            TestZq::new(3),
            TestZq::new(4),
            TestZq::new(5),
            TestZq::new(6),
            TestZq::new(7),
        ];
        let res = zq_vector.l2_norm_squared();
        assert_eq!(res, 140);
    }

    #[test]
    fn test_l2_norm_with_negative_values() {
        let zq_vector = [
            TestZq::new(1),
            TestZq::new(2),
            TestZq::new(3),
            -TestZq::new(4),
            -TestZq::new(5),
            -TestZq::new(6),
            -TestZq::new(7),
        ];
        let res = zq_vector.l2_norm_squared();
        assert_eq!(res, 140);
    }

    #[test]
    fn test_linf_norm() {
        let zq_vector = [
            TestZq::new(1),
            TestZq::new(200),
            TestZq::new(300),
            TestZq::new(40),
            -TestZq::new(5),
            -TestZq::new(6),
            -TestZq::new(700000),
        ];
        let res = zq_vector.linf_norm();
        assert_eq!(res, 700000);
    }
}

#[cfg(test)]
mod decomposition_tests {
    use crate::ring::{zq::Zq, Norms};
    use super::*;

    type TestZq = ZqU32Max;

    #[test]
    fn test_zq_decomposition() {
        let (base, parts) = (TestZq::new(12), 10);
        let pos_zq = TestZq::new(29);
        let neg_zq = -TestZq::new(29);

        let pos_decomposed = pos_zq.decompose(base, parts);
        let neg_decomposed = neg_zq.decompose(base, parts);

        assert_eq!(
            pos_decomposed,
            vec![
                TestZq::new(5),
                TestZq::new(2),
                TestZq::ZERO,
                TestZq::ZERO,
                TestZq::ZERO,
                TestZq::ZERO,
                TestZq::ZERO,
                TestZq::ZERO,
                TestZq::ZERO,
                TestZq::ZERO
            ]
        );
        assert_eq!(
            neg_decomposed,
            vec![
                -TestZq::new(5),
                -TestZq::new(2),
                TestZq::ZERO,
                TestZq::ZERO,
                TestZq::ZERO,
                TestZq::ZERO,
                TestZq::ZERO,
                TestZq::ZERO,
                TestZq::ZERO,
                TestZq::ZERO
            ]
        );
    }

    #[test]
    fn test_zq_recomposition() {
        let (base, parts) = (TestZq::new(1802), 10);
        let pos_zq = -TestZq::new(16200);

        let pos_decomposed = pos_zq.decompose(base, parts);
        let mut exponential_base = TestZq::new(1);
        let mut result = TestZq::new(0);
        for decomposed_part in pos_decomposed {
            result += decomposed_part * exponential_base;
            exponential_base *= base;
        }
        assert_eq!(result, pos_zq)
    }
}