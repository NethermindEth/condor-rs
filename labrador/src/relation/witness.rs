use rand::rng;

use crate::ring::{rq::Rq, rq_vector::RqVector, Norms};

pub struct Witness {
    pub s: Vec<RqVector>,
}

impl Witness {
    pub fn new(rank: usize, multiplicity: usize, bound_sq: u128) -> Self {
        #[allow(clippy::as_conversions)]
        let std = ((bound_sq as f64) / ((rank * multiplicity * Rq::DEGREE) as f64)).sqrt();
        #[allow(clippy::as_conversions)]
        let std = std as u32;
        loop {
            let s: Vec<RqVector> = (0..multiplicity)
                .map(|_| RqVector::random_with_bound(&mut rng(), rank, std))
                .collect();
            if Self::validate_l2_norm(&s, bound_sq) {
                return Self { s };
            }
        }
    }

    fn validate_l2_norm(candidate: &[RqVector], bound_sq: u128) -> bool {
        for witness in candidate {
            if witness.l2_norm_squared() > bound_sq {
                return false;
            }
        }
        true
    }
}

#[cfg(test)]
mod tests {

    use crate::ring::{rq::Rq, rq_vector::RqVector, zq::Zq, Norms};

    use super::Witness;

    #[test]
    fn test_witness_vector_norm() {
        let bound_sq = 320000;
        let witness_vector = Witness::new(40, 100, bound_sq);
        assert_eq!(witness_vector.s.len(), 100);
        assert_eq!(witness_vector.s[0].len(), 40);
        for witness in witness_vector.s.iter() {
            let l2_norm = witness.l2_norm_squared();
            assert!(l2_norm < bound_sq * bound_sq)
        }
    }

    #[test]
    fn test_witness_with_larger_bound() {
        let poly1 = Rq::new([Zq::new(10000); Rq::DEGREE]);
        let poly2 = Rq::new([Zq::new(142310); Rq::DEGREE]);
        let poly3 = Rq::new([Zq::new(9310); Rq::DEGREE]);
        let poly_vector = vec![poly1, poly2, poly3];
        assert!(!Witness::validate_l2_norm(
            &[RqVector::new(poly_vector)],
            1000
        ));
    }
}
