use crate::{rq::Rq, rq_vector::RqVector, zq::Zq};
use rand::distr::{Distribution, Uniform};
use rand::{CryptoRng, Rng};

#[derive(Debug, Clone, PartialEq, PartialOrd, Eq, Ord, Default)]
pub struct PolyRing {
    coeffs: Vec<Zq>,
}
impl PolyRing {
    // degree bound of the polynomial ring, set 64 from the paper
    pub const DEGREE_BOUND: usize = 64;

    pub fn new(coeffs: Vec<Zq>) -> Self {
        Self { coeffs }
    }
    pub fn zero(length: usize) -> Self {
        Self {
            coeffs: vec![Zq::ZERO; length],
        }
    }

    pub fn zero_poly() -> Self {
        Self {
            coeffs: vec![Zq::ZERO; 1],
        }
    }

    pub fn len(&self) -> usize {
        self.coeffs.len()
    }

    pub fn is_empty(&self) -> bool {
        self.coeffs.is_empty()
    }

    pub fn get_coeffs(&self) -> &Vec<Zq> {
        &self.coeffs
    }

    pub fn iter(&self) -> impl Iterator<Item = &Zq> {
        self.coeffs.iter()
    }

    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut Zq> {
        self.coeffs.iter_mut()
    }

    /// inner product of two polynomials
    pub fn inner_product(&self, other: &Self) -> Zq {
        self.coeffs
            .iter()
            .zip(other.coeffs.iter())
            .map(|(a, b)| *a * *b)
            .sum()
    }

    /// Add two polynomials with flexible degree
    pub fn add(&self, other: &Self) -> Self {
        let max_degree = self.get_coeffs().len().max(other.get_coeffs().len());
        let mut coeffs = vec![Zq::ZERO; max_degree];
        for (i, coeff) in coeffs.iter_mut().enumerate().take(max_degree) {
            if i < self.get_coeffs().len() {
                *coeff += self.get_coeffs()[i];
            }
            if i < other.get_coeffs().len() {
                *coeff += other.get_coeffs()[i];
            }
        }
        Self { coeffs }
    }

    /// Scalar multiplication of a polynomial
    pub fn scalar_mul(&self, other: &Zq) -> Self {
        Self {
            coeffs: self.coeffs.iter().map(|c| *c * *other).collect(),
        }
    }

    /// Polynomial multiplication of two polynomials
    pub fn mul_poly(&self, other: &Self) -> PolyRing {
        // Initialize a vector to hold the intermediate multiplication result
        let mut result_coefficients =
            vec![Zq::new(0); self.get_coeffs().len() + other.get_coeffs().len() - 1];
        for (i, &coeff1) in self.get_coeffs().iter().enumerate() {
            for (j, &coeff2) in other.get_coeffs().iter().enumerate() {
                result_coefficients[i + j] += coeff1 * coeff2;
            }
        }

        // Reduce modulo X^d + 1
        if result_coefficients.len() > Self::DEGREE_BOUND {
            let q_minus_1 = Zq::MAX;
            let (left, right) = result_coefficients.split_at_mut(Self::DEGREE_BOUND);
            for (i, &overflow) in right.iter().enumerate() {
                left[i] += overflow * q_minus_1;
            }
            result_coefficients.truncate(Self::DEGREE_BOUND);
        }
        PolyRing::new(result_coefficients)
    }

    /// Generate random Zq vector with a provided cryptographically secure RNG
    pub fn random<R: Rng + CryptoRng>(rng: &mut R, n: usize) -> Self {
        let uniform = Uniform::new_inclusive(Zq::ZERO, Zq::MAX).unwrap();
        let mut coeffs = Vec::with_capacity(n);
        coeffs.extend((0..n).map(|_| uniform.sample(rng)));
        Self { coeffs }
    }

    pub fn to_rq<const D: usize>(&self) -> Rq<D> {
        Rq::from_vec(self.get_coeffs().clone())
    }

    /// Compute the conjugate automorphism \sigma_{-1} of vector based on B) Constraints..., Page 21.
    pub fn conjugate_automorphism(&self) -> PolyRing {
        let q_minus_1 = Zq::MAX;
        let mut new_coeffs = vec![Zq::ZERO; PolyRing::DEGREE_BOUND];
        for (i, new_coeff) in new_coeffs
            .iter_mut()
            .enumerate()
            .take(PolyRing::DEGREE_BOUND)
        {
            if i < self.get_coeffs().len() {
                if i == 0 {
                    *new_coeff = self.get_coeffs()[i];
                } else {
                    *new_coeff = self.get_coeffs()[i] * q_minus_1;
                }
            } else {
                *new_coeff = Zq::ZERO;
            }
        }
        let reversed_coefficients = new_coeffs
            .iter()
            .take(1)
            .cloned()
            .chain(new_coeffs.iter().skip(1).rev().cloned())
            .collect::<Vec<Zq>>();

        PolyRing::new(reversed_coefficients)
    }
}

impl FromIterator<Zq> for PolyRing {
    fn from_iter<T: IntoIterator<Item = Zq>>(iter: T) -> Self {
        let coeffs: Vec<Zq> = iter.into_iter().collect();
        PolyRing::new(coeffs)
    }
}

#[derive(Debug, Clone, PartialEq, PartialOrd, Eq, Ord, Default)]
pub struct PolyVector {
    elements: Vec<PolyRing>,
}

impl PolyVector {
    pub fn new(elements: Vec<PolyRing>) -> Self {
        Self { elements }
    }

    pub fn zero() -> Self {
        Self {
            elements: vec![PolyRing::zero(0)],
        }
    }

    pub fn get_elements(&self) -> &Vec<PolyRing> {
        &self.elements
    }

    pub fn len(&self) -> usize {
        self.elements.len()
    }

    pub fn is_empty(&self) -> bool {
        self.elements.is_empty()
    }

    pub fn iter(&self) -> impl Iterator<Item = &PolyRing> {
        self.elements.iter()
    }

    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut PolyRing> {
        self.elements.iter_mut()
    }

    // Generate a random polynomial vector with n polynomials of degree m
    pub fn random(n: usize, m: usize) -> PolyVector {
        let mut vector = PolyVector::new(vec![]);
        vector.elements = (0..n)
            .map(|_| PolyRing::random(&mut rand::rng(), m))
            .collect();
        vector
    }

    pub fn add(&self, other: &PolyVector) -> PolyVector {
        self.iter()
            .zip(other.iter())
            .map(|(a, b)| a.add(b))
            .collect()
    }

    pub fn scalar_mul_vector(&self, other: &Zq) -> PolyVector {
        self.iter().map(|a| a.scalar_mul(other)).collect()
    }

    pub fn inner_product_poly_vector(&self, other: &PolyVector) -> PolyRing {
        self.iter()
            .zip(other.iter())
            .map(|(a, b)| a.mul_poly(b))
            .fold(
                PolyRing::zero(self.get_elements()[0].get_coeffs().len()),
                |acc, val| acc.add(&val),
            )
    }

    pub fn to_rqvector<const N: usize, const D: usize>(&self) -> RqVector<N, D> {
        let mut rq_vector = RqVector::zero();
        for (i, poly) in self.elements.iter().enumerate() {
            rq_vector[i] = Rq::from_vec(poly.get_coeffs().clone());
        }
        rq_vector
    }
}

impl FromIterator<PolyRing> for PolyVector {
    fn from_iter<T: IntoIterator<Item = PolyRing>>(iter: T) -> Self {
        let mut elements = Vec::new();
        for item in iter {
            elements.push(item);
        }
        PolyVector::new(elements)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_conjugate_automorphism() {
        let poly1: PolyRing = PolyRing::new(vec![Zq::new(1), Zq::new(2), Zq::new(3)]);
        let poly2: PolyRing = PolyRing::new(vec![Zq::new(4), Zq::new(5), Zq::new(6)]);
        let inner_12 = poly1.inner_product(&poly2);
        let conjugated_1 = poly1.conjugate_automorphism();
        let inner_conjugated_12 = conjugated_1.mul_poly(&poly2);

        assert_eq!(inner_conjugated_12.len(), PolyRing::DEGREE_BOUND);
        assert_eq!(inner_conjugated_12.get_coeffs()[0], Zq::from(32));
        assert_eq!(inner_conjugated_12.get_coeffs()[1], Zq::from(17));
        assert_eq!(inner_conjugated_12.get_coeffs()[2], Zq::new(6));

        // ct<\sigma_{-1}(poly1), poly2> ?= <poly1, poly2>
        let ct_inner_conjugated_12 = inner_conjugated_12.get_coeffs()[0];
        assert_eq!(ct_inner_conjugated_12, inner_12);
    }

    #[test]
    fn test_to_rq() {
        let zq_vector = PolyRing::new(vec![Zq::new(1), Zq::new(2), Zq::new(3)]);
        let rq_vector = zq_vector.to_rq::<3>();
        assert_eq!(
            rq_vector,
            Rq::from_vec(vec![Zq::new(1), Zq::new(2), Zq::new(3)])
        );
    }

    #[test]
    fn test_to_rq_vector() {
        let vector = PolyVector::new(vec![PolyRing::new(vec![
            Zq::new(1),
            Zq::new(2),
            Zq::new(3),
        ])]);
        let rq_vector = vector.to_rqvector::<1, 3>();
        assert_eq!(
            rq_vector[0],
            Rq::from_vec(vec![Zq::new(1), Zq::new(2), Zq::new(3)])
        );
    }

    #[test]
    fn test_scalar_mul_vector_mul_vector() {
        let vector = PolyVector::new(vec![PolyRing::new(vec![
            Zq::new(1),
            Zq::new(2),
            Zq::new(3),
        ])]);
        let zq = Zq::new(2);
        let result = vector.scalar_mul_vector(&zq);
        println!("result: {:?}", result.get_elements());
    }

    #[test]
    fn test_add_poly_vector() {
        let vector1 = PolyVector::new(vec![PolyRing::new(vec![
            Zq::new(1),
            Zq::new(2),
            Zq::new(3),
        ])]);
        let vector2 = PolyVector::new(vec![PolyRing::new(vec![Zq::new(4), Zq::new(5)])]);
        let result = vector1.add(&vector2);
        println!("result: {:?}", result.get_elements());
    }
}
