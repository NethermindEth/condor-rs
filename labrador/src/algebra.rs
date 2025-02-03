pub mod rq;
pub mod zq;

use std::fmt::Debug;
use std::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};

// TODO: Do we also need operators for Rhs=Base?
pub trait RingOps<Output>:
    Sized + Add<Output = Output> + Sub<Output = Output> + Mul<Output = Output>
{
}

pub trait Ring:
    AddAssign + SubAssign + MulAssign + PartialOrd + Clone + PartialEq + Eq + Debug
where
    for<'a> &'a Self: RingOps<Self>,
{
    const ZERO: Self;
    const ONE: Self;
}
