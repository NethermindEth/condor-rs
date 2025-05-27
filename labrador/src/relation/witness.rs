use rand::rng;

use crate::ring::{rq::Rq, rq_vector::RqVector, zq::Zq};

pub struct Witness {
    pub s: Vec<RqVector>,
}

impl Witness {
    pub fn new(rank: usize, multiplicity: usize, bound: Zq) -> Self {
        #[allow(clippy::as_conversions)]
        let std = (bound.get_value() as f64) / f64::sqrt((rank * multiplicity * Rq::DEGREE) as f64);
        #[allow(clippy::as_conversions)]
        let std = std as u32;
        loop {
            let s: Vec<RqVector> = (0..multiplicity)
                .map(|_| RqVector::random_with_bound(&mut rng(), rank, std))
                .collect();
            if Self::validate_l2_norm(&s, bound) {
                return Self { s };
            }
        }
    }

    fn validate_l2_norm(candidate: &[RqVector], bound: Zq) -> bool {
        for witness in candidate {
            if witness.l2_norm_squared() > bound * bound {
                return false;
            }
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use crate::ring::{rq_vector::RqVector, zq::Zq};

    use super::Witness;

    #[test]
    fn test_witness_vector_norm() {
        let bound = Zq::new(320000);
        let witness_vector = Witness::new(40, 100, bound);
        assert_eq!(witness_vector.s.len(), 100);
        assert_eq!(witness_vector.s[0].get_length(), 40);
        for witness in witness_vector.s.iter() {
            let l2_norm = RqVector::l2_norm_squared(witness);
            assert!(l2_norm < bound * bound)
        }
    }
}
