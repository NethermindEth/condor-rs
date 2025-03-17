use crate::zq::Zq;
use rand::prelude::*;
use rand::seq::SliceRandom;
use rustfft::{num_complex::Complex, FftPlanner};

/// Compute the operator norm of a polynomial given its coefficients.
/// The operator norm is defined as the maximum magnitude of the DFT (eigenvalues)
/// of the coefficient vector.
pub fn operator_norm(coeffs: &[f64]) -> f64 {
    let n = coeffs.len();
    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_forward(n);

    // Convert coefficients into complex numbers (with zero imaginary parts)
    let mut buffer: Vec<Complex<f64>> =
        coeffs.iter().map(|&x| Complex { re: x, im: 0.0 }).collect();

    // Compute the FFT (this gives the eigenvalues of the circulant matrix)
    fft.process(&mut buffer);

    // Return the maximum absolute value (norm) among the eigenvalues
    buffer
        .iter()
        .map(|c| c.norm())
        .fold(0.0, |max, x| max.max(x))
}

/// Challenge Space, page 6 in the paper:
/// Generates a candidate challenge polynomial with 64 coefficients:
/// - 23 coefficients are 0
/// - 31 coefficients are chosen uniformly from {+1, -1}
/// - 10 coefficients are chosen uniformly from {+2, -2}
///
/// The coefficients are then shuffled randomly.
pub fn sample_challenge() -> Vec<Zq> {
    let mut rng = rand::rng();
    let mut challenge: Vec<Zq> = Vec::new();

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
    challenge
}

/// Rejection sampling: repeatedly sample candidate challenges until one has an operator norm
/// less than the specified threshold. Returns the accepted challenge and the number of samples tried.
pub fn sample_challenge_with_norm(threshold: f64) -> Vec<Zq> {
    // let mut count = 0;
    loop {
        // count += 1;
        let candidate = sample_challenge();
        let candidate_f64 = zq_to_f64(candidate.clone());
        let norm = operator_norm(&candidate_f64);
        if norm < threshold {
            return candidate.clone();
        }
    }
}

pub fn zq_to_f64(zq: Vec<Zq>) -> Vec<f64> {
    let _zq = zq.iter().map(|z| z.to_f64()).collect();
    _zq
}

#[cfg(test)]
mod tests {
    use super::*;

    // Test Size correctness of projection matrix
    #[test]
    fn test_sample_challenge_with_norm() {
        // threshold \in R, threshold is "T" in the paper(Challenge Space, page 6) which is at most 15.
        let threshold = 15.0;
        let challenge = sample_challenge_with_norm(threshold);
        let candidate_f64 = zq_to_f64(challenge.clone());
        let norm = operator_norm(&candidate_f64);
        println!(
            "Accepted challenge polynomial coefficients: {:?}",
            challenge
        );
        assert!(norm <= threshold);
    }
}
