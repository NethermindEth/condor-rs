use crate::algebra::{Ring, RingOps};
use std::fmt::{Debug, Display, Formatter};
use std::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};

type ValueType = u64;

/// A quotient ring `Z_q` of integers mod `q`.
/// Q is hardcoded to 2^32.
// TODO: Should we use a field from the crate ark_ff? Note that they are mostly focusing on prime fields
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Zq {
    pub value: ValueType,
}

impl Zq {
    pub const Q: ValueType = 2_u64.pow(32);

    pub fn new(value: ValueType) -> Self {
        Self {
            value: value % Self::Q,
        }
    }
}

impl Display for Zq {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.value)
    }
}

impl Add for Zq {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Zq::new(self.value + rhs.value)
    }
}

impl Add for &Zq {
    type Output = Zq;

    fn add(self, rhs: Self) -> Self::Output {
        Zq::new(self.value + rhs.value)
    }
}

impl Sub for Zq {
    type Output = Zq;

    fn sub(self, rhs: Self) -> Self::Output {
        Zq::new(self.value + Zq::Q - rhs.value)
    }
}

impl Sub for &Zq {
    type Output = Zq;

    fn sub(self, rhs: Self) -> Self::Output {
        Zq::new(self.value + Zq::Q - rhs.value)
    }
}

impl Mul for Zq {
    type Output = Zq;

    fn mul(self, rhs: Self) -> Self::Output {
        Zq::new(self.value * rhs.value)
    }
}

impl Mul for &Zq {
    type Output = Zq;

    fn mul(self, rhs: Self) -> Self::Output {
        Zq::new(self.value * rhs.value)
    }
}

impl AddAssign for Zq {
    fn add_assign(&mut self, rhs: Self) {
        self.value = (self.value + rhs.value) % Zq::Q;
    }
}

impl SubAssign for Zq {
    fn sub_assign(&mut self, rhs: Self) {
        self.value = (self.value + Self::Q - rhs.value) % Self::Q;
    }
}

impl MulAssign for Zq {
    fn mul_assign(&mut self, rhs: Self) {
        self.value = (self.value * rhs.value) % Self::Q;
    }
}

impl Debug for Zq {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        // No need for special formatting, just use Display
        <Self as Display>::fmt(self, f)
    }
}

impl RingOps<Zq> for Zq {}
impl RingOps<Zq> for &Zq {}

impl Ring for Zq {
    const ZERO: Self = Zq { value: 0 };
    const ONE: Self = Zq { value: 1 };
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zq_new() {
        let value = Zq::Q + 1;
        let zq = Zq::new(value);
        assert_eq!(zq.value, 1);
    }

    #[test]
    fn test_zq_add() {
        let a = Zq::new(Zq::Q - 1);
        assert_eq!((a + Zq::new(1)).value, 0);
        assert_eq!((a + Zq::new(2)).value, 1);
    }

    #[test]
    fn test_zq_sub() {
        let a = Zq::new(Zq::Q - 1);
        let b = Zq::new(Zq::Q - 2);
        let result = a - b;
        assert_eq!(result.value, 1);
    }

    #[test]
    fn test_zq_mul() {
        let a = Zq::new(6);
        let b = Zq::new(7);
        let result = a * b;
        assert_eq!(result.value, 42);
    }

    #[test]
    fn test_zq_mul_overflow() {
        let a = Zq::new(Zq::Q - 1);
        let b = Zq::new(2);
        let result = a * b;
        let expected = Zq::new(Zq::Q - 2); // -2
        assert_eq!(result, expected);
    }

    // #[test]
    // fn test_zq_rem() {
    //     let a = Zq::new(10);
    //     let b = Zq::new(3);
    //     let result = a % b;
    //     assert_eq!(result.value, 1);
    // }

    // #[test]
    // fn test_zq_ceil_div() {
    //     let mut a;
    //     let mut b;
    //     a = Zq::new(10);
    //     b = Zq::new(3);
    //     assert_eq!(a.div_ceil(b).value(), 4);
    //
    //     a = Zq::new(9);
    //     b = Zq::new(3);
    //     assert_eq!(a.div_ceil(b).value(), 3);
    // }
}
