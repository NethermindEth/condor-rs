use crate::ring::poly::PolyRing;
use crate::ring::zq::Zq;
use rand::prelude::*;
use rand::seq::SliceRandom;

pub struct ChallengeSet {
    challenges: PolyRing,
    norm: f64,
}

/// Challenge Space over R = Z_q\[X\] / (X^d + 1), paper page 6:
impl ChallengeSet {
    pub fn new(deg_bound_d: usize) -> Self {
        // threshold \in R, threshold is "T" in the paper(Challenge Space, page 6) which is at most 15.
        let op_norm = 15.0;
        // Sample challenges with a given norm.
        let challenges: PolyRing = sample_challenge_with_norm(deg_bound_d, op_norm);
        let norm = PolyRing::operator_norm(&challenges);
        ChallengeSet { challenges, norm }
    }

    /// Returns a reference to the challenge polynomial.
    pub fn get_challenges(&self) -> &PolyRing {
        &self.challenges
    }

    /// Returns the norm of the challenges.
    pub fn get_norm(&self) -> f64 {
        self.norm
    }
}

/// Generates a candidate challenge polynomial with 64 coefficients:
/// - 23 coefficients are 0
/// - 31 coefficients are chosen uniformly from {+1, -1}
/// - 10 coefficients are chosen uniformly from {+2, -2}
///
/// The coefficients are then shuffled randomly.
fn sample_challenge(deg_bound_d: usize) -> PolyRing {
    let mut rng = rand::rng();
    let mut challenge: Vec<Zq> = Vec::with_capacity(deg_bound_d);

    // Add 23 zeros.
    for _ in 0..23 {
        challenge.push(Zq::ZERO);
    }
    // Add 31 coefficients, each either +1 or -1.
    for _ in 0..31 {
        let value = if rng.random_bool(0.5) {
            Zq::ONE
        } else {
            -Zq::ONE
        };
        challenge.push(value);
    }
    // Add 10 coefficients, each either +2 or -2.
    for _ in 0..10 {
        let value = if rng.random_bool(0.5) {
            Zq::new(2)
        } else {
            -Zq::new(2)
        };
        challenge.push(value);
    }

    // Shuffle the vector to randomize the positions.
    challenge.shuffle(&mut rng);
    PolyRing::new(challenge)
}

/// Rejection sampling: repeatedly sample candidate challenges until one has an operator norm
/// less than the specified threshold. Returns the accepted challenge and the number of samples tried.
fn sample_challenge_with_norm(deg_bound_d: usize, threshold: f64) -> PolyRing {
    loop {
        let candidate = sample_challenge(deg_bound_d);
        let norm = PolyRing::operator_norm(&candidate);
        if norm < threshold {
            return candidate;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Test Challenge Set:
    /// const 71 and const 15 are from Challenge Space, paper page 6.
    /// l2 norm <= 71
    /// operator norm <= 15
    #[test]
    fn test_challenge_set() {
        let op_norm = 15.0;
        let deg_bound_d: usize = 8;
        let cs_1: ChallengeSet = ChallengeSet::new(deg_bound_d);
        let norm = cs_1.get_norm();
        let challenges_1 = cs_1.get_challenges();

        let cs_2: ChallengeSet = ChallengeSet::new(deg_bound_d);
        let challenges_2 = cs_2.get_challenges();

        // l2 norm 71 is from paper page 6, Challenge Space.
        assert_eq!(challenges_1.inner_product(challenges_1), Zq::new(71));
        assert!(norm <= op_norm);

        assert_ne!(challenges_1, challenges_2);
    }
}
