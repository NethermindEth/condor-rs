use crate::ring::rq_vector::RqVector;
use crate::ring::zq::Zq;

// LaBRADOR: Compact Proofs for R1CS from Module-SIS | Page 5 | Proving smallness section
const UPPER_BOUND_FACTOR: Zq = Zq::new(128);
const LOWER_BOUND_FACTOR: Zq = Zq::new(30);

pub struct Projection {
    random_linear_map_vector: Vec<Vec<Vec<Zq>>>,
    security_level: usize,
}

impl Projection {
    pub fn new(random_linear_map_vector: Vec<Vec<Vec<Zq>>>, security_level: usize) -> Self {
        Self {
            random_linear_map_vector,
            security_level,
        }
    }

    fn compute_projection(&self, index: usize, witness: &RqVector) -> Vec<Zq> {
        let mut projection = vec![Zq::ZERO; 2 * self.security_level];
        let coefficients = witness.concatenate_coefficients();
        for (i, item) in projection.iter_mut().enumerate() {
            *item = self.random_linear_map_vector[index][i]
                .iter()
                .zip(coefficients.iter())
                .map(|(m, s)| *m * *s)
                .sum::<Zq>();
        }
        projection
    }

    pub fn compute_batch_projection(&self, witness_vector: &[RqVector]) -> Vec<Zq> {
        let mut result = vec![Zq::ZERO; 2 * self.security_level];
        for (index_i, witness) in witness_vector.iter().enumerate() {
            for (index_j, element) in result.iter_mut().enumerate() {
                *element += self.compute_projection(index_i, witness)[index_j];
            }
        }
        result
    }

    pub fn get_projection_matrices(&self) -> &[Vec<Vec<Zq>>] {
        &self.random_linear_map_vector
    }

    fn norm_squared(projection: &[Zq]) -> Zq {
        projection.iter().map(|coeff| *coeff * *coeff).sum()
    }

    // Function to verify upper bound of projection
    pub fn verify_projection_upper_bound(projection: &[Zq], beta_squared: Zq) -> bool {
        Self::norm_squared(projection) < (UPPER_BOUND_FACTOR * beta_squared)
    }

    // Function to verify lower bound of projection
    pub fn verify_projection_lower_bound(projection: &[Zq], beta_squared: Zq) -> bool {
        Self::norm_squared(projection) > (LOWER_BOUND_FACTOR * beta_squared)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::transcript::sponges::shake::ShakeSponge;
    use crate::transcript::{LabradorTranscript, Sponge};
    use rand::rng;

    // Test that the probability of the inequality being true is close to 1/2
    #[test]
    #[cfg(not(feature = "skip-slow-tests"))]
    fn test_projection_is_smaller_than_upper_bound() {
        use crate::transcript::Sponge;

        let (security_parameter, rank, multiplicity) = (128, 5, 1);
        // 1000 was chosen to provide a reasonably large sample size

        let trials: f64 = 1000.0;
        let mut success_count: f64 = 0.0;
        for _ in 0..1000 {
            let mut transcript = LabradorTranscript::new(ShakeSponge::default());
            // This gives randomness to the transcript, to generate random projection matrices.
            transcript.set_u1(RqVector::random(&mut rng(), 1));
            let projections =
                transcript.generate_projections(security_parameter, rank, multiplicity);

            let witness = RqVector::random(&mut rand::rng(), rank);
            let result = projections.compute_projection(0, &witness);

            let beta = witness.compute_norm_squared();
            // Check if the norm of the projection is smaller than 128 * (squared norm of the projection of the random polynomial)
            let test: bool = Projection::verify_projection_upper_bound(&result, beta);
            if test {
                success_count += 1.0;
            }
        }

        let observed_probability = success_count / trials;

        // We allow some tolerance because of the statistical nature of the results.
        let tolerance = 0.05;
        assert!(
            (observed_probability - 0.5).abs() < tolerance,
            "Observed probability {} is not close to 0.5",
            observed_probability
        );
    }

    // On average the projected norm squared is the same as 128 * vector norm squared
    #[test]
    #[cfg(not(feature = "skip-slow-tests"))]
    fn test_projection_average_value() {
        use crate::transcript::Sponge;

        let (security_parameter, rank, multiplicity) = (128, 3, 1);
        let trials: u128 = 10000;

        let witness = RqVector::random_ternary(&mut rand::rng(), rank);
        let witness_norm = (Zq::new(security_parameter.try_into().unwrap())
            * RqVector::compute_norm_squared(&witness))
        .to_u128();

        let mut norm_sum = Zq::ZERO;
        // Run the test multiple times to simulate the probability
        for _ in 0..trials {
            let mut transcript = LabradorTranscript::new(ShakeSponge::default());
            // This gives randomness to the transcript, to generate random projection matrices.
            transcript.set_u1(RqVector::random(&mut rng(), 1));
            let projections =
                transcript.generate_projections(security_parameter, rank, multiplicity);
            let result = projections.compute_projection(0, &witness);
            norm_sum += Projection::norm_squared(&result);
        }

        // Calculate the observed probability
        let average = norm_sum.to_u128() / trials;
        let difference = if witness_norm <= average {
            average - witness_norm
        } else {
            witness_norm - average
        };

        // we choose a small tolerance value for possible statistical error
        let tolerance: u128 = 50;
        assert!(
            difference < tolerance,
            "Average norm value {} is not equal to {}.",
            average,
            Zq::new(security_parameter.try_into().unwrap()) * witness.compute_norm_squared(),
        );
    }

    // Test lower bound verification
    #[test]
    fn test_lower_bound() {
        let (security_parameter, rank, multiplicity) = (128, 5, 1);

        let mut transcript = LabradorTranscript::new(ShakeSponge::default());
        transcript.set_u1(RqVector::random(&mut rng(), 1));
        let projections = transcript.generate_projections(security_parameter, rank, multiplicity);
        let witness = RqVector::random_ternary(&mut rng(), rank);

        let beta = witness.compute_norm_squared();
        // Check if the norm of the projection is bigger than 30 * (squared norm of the projection of the random polynomial)
        assert!(Projection::verify_projection_lower_bound(
            &projections.compute_projection(0, &witness),
            beta
        ));
    }
}
