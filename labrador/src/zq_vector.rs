use crate::{rq::Rq, rq_vector::RqVector, zq::Zq};
use rand::distr::{Distribution, Uniform};
use rand::{CryptoRng, Rng};

#[derive(Debug, Clone, PartialEq, PartialOrd, Eq, Ord, Default)]
pub struct ZqVector {
    coeffs: Vec<Zq>,
}

impl ZqVector {
    pub fn new(coeffs: Vec<Zq>) -> Self {
        Self { coeffs }
    }

    pub fn zero() -> Self {
        Self { coeffs: vec![] }
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

    /// Generate random Zq vector with a provided cryptographically secure RNG
    pub fn random<R: Rng + CryptoRng>(rng: &mut R, n: usize) -> Self {
        let uniform = Uniform::new_inclusive(Zq::ZERO, Zq::MAX).unwrap();
        let mut coeffs = Vec::with_capacity(n);
        coeffs.iter_mut().for_each(|c| *c = uniform.sample(rng));
        Self { coeffs }
    }
    pub fn to_rq<const D: usize>(&self) -> Rq<D> {
        Rq::from_vec(self.get_coeffs().clone())
    }
}

#[derive(Debug, Clone, PartialEq, PartialOrd, Eq, Ord, Default)]
pub struct ZqMatrix {
    elements: Vec<ZqVector>,
}

impl ZqMatrix {
    pub fn new(elements: Vec<ZqVector>) -> Self {
        Self { elements }
    }

    pub fn get_elements(&self) -> &Vec<ZqVector> {
        &self.elements
    }

    pub fn get_row(&self, i: usize) -> &ZqVector {
        &self.elements[i]
    }

    pub fn len(&self) -> usize {
        self.elements.len()
    }

    pub fn is_empty(&self) -> bool {
        self.elements.is_empty()
    }

    pub fn iter(&self) -> impl Iterator<Item = &ZqVector> {
        self.elements.iter()
    }

    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut ZqVector> {
        self.elements.iter_mut()
    }

    pub fn random_matrix(n: usize, m: usize) -> ZqMatrix {
        let mut matrix = ZqMatrix::new(vec![]);
        for _ in 0..n {
            matrix.elements.push(ZqVector::random(&mut rand::rng(), m));
        }
        matrix
    }

    // Compute the conjugate automorphism \sigma_{-1} of vector based on B) Constraints..., Page 21.
    pub fn conjugate_automorphism(&self) -> ZqMatrix {
        let mut new_matrix = ZqMatrix::new(vec![]);
        for poly in self.elements.iter() {
            let poly_len = poly.len();
            let mut new_coeffs = vec![Zq::ZERO; poly_len];
            new_coeffs[0] = poly.get_coeffs()[0];

            for i in 1..poly_len {
                let j = poly_len - i;
                new_coeffs[j] = -poly.get_coeffs()[i];
            }

            new_matrix.elements.push(ZqVector::new(new_coeffs));
        }

        new_matrix
    }

    pub fn to_rqvector<const N: usize, const D: usize>(&self) -> RqVector<N, D> {
        let mut rq_vector = RqVector::zero();
        for (i, poly) in self.elements.iter().enumerate() {
            rq_vector[i] = Rq::from_vec(poly.get_coeffs().clone());
        }
        rq_vector
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_conjugate_automorphism() {
        let matrix_1 = ZqMatrix::new(vec![ZqVector::new(vec![
            Zq::new(1),
            Zq::new(2),
            Zq::new(3),
        ])]);
        let conjugated_1 = matrix_1.conjugate_automorphism();
        let expected_1 = ZqMatrix::new(vec![ZqVector::new(vec![
            Zq::new(1),
            -Zq::new(3),
            -Zq::new(2),
        ])]);
        assert_eq!(conjugated_1.get_elements(), expected_1.get_elements());

        let matrix_2 = ZqMatrix::new(vec![
            ZqVector::new(vec![Zq::new(1), Zq::new(2), Zq::new(3), Zq::new(4)]),
            ZqVector::new(vec![Zq::new(5), Zq::new(6), Zq::new(7), Zq::new(8)]),
        ]);
        let conjugated_2 = matrix_2.conjugate_automorphism();
        let expected_2 = ZqMatrix::new(vec![
            ZqVector::new(vec![Zq::new(1), -Zq::new(4), -Zq::new(3), -Zq::new(2)]),
            ZqVector::new(vec![Zq::new(5), -Zq::new(8), -Zq::new(7), -Zq::new(6)]),
        ]);
        assert_eq!(conjugated_2.get_elements(), expected_2.get_elements());
    }

    #[test]
    fn test_to_rq() {
        let zq_vector = ZqVector::new(vec![Zq::new(1), Zq::new(2), Zq::new(3)]);
        let rq_vector = zq_vector.to_rq::<3>();
        assert_eq!(
            rq_vector,
            Rq::from_vec(vec![Zq::new(1), Zq::new(2), Zq::new(3)])
        );
    }

    #[test]
    fn test_to_rq_vector() {
        let matrix = ZqMatrix::new(vec![ZqVector::new(vec![
            Zq::new(1),
            Zq::new(2),
            Zq::new(3),
        ])]);
        let rq_vector = matrix.to_rqvector::<1, 3>();
        assert_eq!(
            rq_vector[0],
            Rq::from_vec(vec![Zq::new(1), Zq::new(2), Zq::new(3)])
        );
    }
}
