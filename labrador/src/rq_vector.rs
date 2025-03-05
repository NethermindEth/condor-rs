use crate::rq::Rq;
use core::ops::{Index, IndexMut, Mul};
use core::slice::Iter;
use rand::{CryptoRng, Rng};

/// Vector of polynomials in Rq
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RqVector<const N: usize, const D: usize> {
    elements: Vec<Rq<D>>,
}

impl<const N: usize, const D: usize> RqVector<N, D> {
    /// Create a zero vector
    pub fn zero() -> Self {
        Self {
            elements: vec![Rq::zero(); N],
        }
    }

    /// Create a random vector
    pub fn random<R: Rng + CryptoRng>(rng: &mut R) -> Self {
        Self {
            elements: (0..N).map(|_| Rq::random(rng)).collect(),
        }
    }

    /// Create a random vector
    pub fn random_small<R: Rng + CryptoRng>(rng: &mut R) -> Self {
        Self {
            elements: (0..N).map(|_| Rq::random_ternary(rng)).collect(),
        }
    }

    /// Get the underlying vector as slice
    pub fn as_slice(&self) -> &[Rq<D>] {
        &self.elements
    }

    pub fn iter(&self) -> Iter<'_, Rq<D>> {
        self.elements.iter()
    }

    /// Convert into array
    pub fn into_array(self) -> [Rq<D>; N] {
        self.elements
            .try_into()
            .unwrap_or_else(|_| panic!("Vector length mismatch"))
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

/// Create a new vector from a `Vec` of elements
impl<const N: usize, const D: usize> From<Vec<Rq<D>>> for RqVector<N, D> {
    fn from(elements: Vec<Rq<D>>) -> Self {
        assert_eq!(elements.len(), N, "Vector length must match N");
        Self { elements }
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
