use crate::rq::Rq;
use std::ops::{Index, IndexMut, Mul};

/// Vector of polynomials in Rq
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RqVector<const N: usize, const D: usize> {
    elements: [Rq<D>; N],
}

impl<const N: usize, const D: usize> RqVector<N, D> {
    /// Create a new vector from an array of elements
    pub fn new(elements: [Rq<D>; N]) -> Self {
        Self { elements }
    }

    /// Create a zero vector
    pub fn zero() -> Self {
        Self {
            elements: std::array::from_fn(|_| Rq::zero()),
        }
    }

    /// Create a random small vector
    pub fn random_small() -> Self {
        Self {
            elements: std::array::from_fn(|_| Rq::random_small()),
        }
    }

    /// Get the underlying array
    pub fn as_array(&self) -> &[Rq<D>; N] {
        &self.elements
    }

    /// Convert into the underlying array
    pub fn into_array(self) -> [Rq<D>; N] {
        self.elements
    }
}

// Enable array-like indexing
impl<const N: usize, const D: usize> Index<usize> for RqVector<N, D> {
    type Output = Rq<D>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.elements[index]
    }
}

impl<const N: usize, const D: usize> IndexMut<usize> for RqVector<N, D> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.elements[index]
    }
}

// From array conversion
impl<const N: usize, const D: usize> From<[Rq<D>; N]> for RqVector<N, D> {
    fn from(elements: [Rq<D>; N]) -> Self {
        Self::new(elements)
    }
}

// Into array conversion
impl<const N: usize, const D: usize> From<RqVector<N, D>> for [Rq<D>; N] {
    fn from(vector: RqVector<N, D>) -> Self {
        vector.into_array()
    }
}

// Dot product between two vectors
impl<const N: usize, const D: usize> Mul for &RqVector<N, D> {
    type Output = Rq<D>;

    fn mul(self, rhs: Self) -> Self::Output {
        self.elements
            .iter()
            .zip(rhs.elements.iter())
            .map(|(a, b)| a.clone() * b.clone())
            .fold(Rq::zero(), |acc, x| acc + x)
    }
}
