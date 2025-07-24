use crate::ring::rq_matrix::RqMatrix;
use crate::ring::rq_vector::RqVector;
use crate::ring::zq::Zq;
use crate::ring::{self, Norms};

// LaBRADOR: Compact Proofs for R1CS from Module-SIS | Page 5 | Proving smallness section
const UPPER_BOUND_FACTOR: u128 = 128;
const LOWER_BOUND_FACTOR: u128 = 30;

pub struct Projection {
    random_linear_map_vector: Vec<RqMatrix>,
    security_level: usize,
}

impl Projection {
    pub fn new(random_linear_map_vector: Vec<RqMatrix>, security_level: usize) -> Self {
        Self {
            random_linear_map_vector,
            security_level,
        }
    }

    fn compute_projection(&self, index: usize, witness: &RqVector) -> Vec<Zq> {
        let mut projection = vec![Zq::ZERO; 2 * self.security_level];
        let coefficients = witness.concatenate_coefficients();
        for (i, pi_ij) in self.random_linear_map_vector[index]
            .elements()
            .iter()
            .enumerate()
        {
            projection[i] = pi_ij
                .concatenate_coefficients()
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
            ring::zq::add_assign_two_zq_vectors(
                &mut result,
                self.compute_projection(index_i, witness),
            );
        }
        result
    }

    pub fn get_projection_matrices(&self) -> &[RqMatrix] {
        &self.random_linear_map_vector
    }

    pub fn get_conjugated_projection_matrices(&self) -> Vec<RqMatrix> {
        self.random_linear_map_vector
            .iter()
            .map(|pi_i| {
                pi_i.elements()
                    .iter()
                    .map(|pi_ij| {
                        pi_ij
                            .elements()
                            .iter()
                            .map(|polynomial| polynomial.conjugate_automorphism())
                            .collect()
                    })
                    .collect()
            })
            .collect()
    }

    // Function to verify upper bound of projection
    pub fn verify_projection_upper_bound(projection: &[Zq], beta_squared: u128) -> bool {
        projection.l2_norm_squared() <= (UPPER_BOUND_FACTOR * beta_squared)
    }

    // Function to verify lower bound of projection
    pub fn verify_projection_lower_bound(projection: &[Zq], beta_squared: u128) -> bool {
        projection.l2_norm_squared() > (LOWER_BOUND_FACTOR * beta_squared)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::relation::witness::Witness;
    use crate::ring::Norms;
    use crate::transcript::sponges::shake::ShakeSponge;
    use crate::transcript::LabradorTranscript;
    use rand::rng;

    // Test that the probability of the inequality being true is close to 1/2
    #[test]
    #[cfg(not(feature = "skip-slow-tests"))]
    fn test_projection_is_smaller_than_upper_bound() {
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

            let witness_vector = Witness::new(rank, multiplicity, 6400 * 6400).s;
            let result = projections.compute_projection(0, &witness_vector[0]);
            let beta = witness_vector[0].l2_norm_squared();
            // dbg!(&result);
            dbg!(result.l2_norm_squared());
            dbg!(UPPER_BOUND_FACTOR * beta);
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
        use crate::{relation::witness::Witness, ring::Norms};

        let (security_parameter, rank, multiplicity) = (128, 3, 1);
        let trials: u128 = 10000;

        // let witness = RqVector::random_ternary(&mut rand::rng(), rank);
        let witness = Witness::new(rank, multiplicity, Zq::NEG_ONE).s;
        let witness_norm = 128 * witness[0].l2_norm_squared();

        let mut norm_sum = 0u128;
        // Run the test multiple times to simulate the probability
        for _ in 0..trials {
            let mut transcript = LabradorTranscript::new(ShakeSponge::default());
            // This gives randomness to the transcript, to generate random projection matrices.
            transcript.set_u1(RqVector::random(&mut rng(), 1));
            let projections =
                transcript.generate_projections(security_parameter, rank, multiplicity);
            let result = projections.compute_projection(0, &witness[0]);
            norm_sum += result.l2_norm_squared();
        }

        // Calculate the observed probability
        let average = norm_sum as f64 / trials as f64;
        let ratio = if witness_norm as f64 <= average {
            average / witness_norm as f64
        } else {
            witness_norm as f64 / average
        };

        // we choose a small tolerance value for possible statistical error
        let tolerance_percent: f64 = 1.05;
        assert!(
            ratio < tolerance_percent,
            "Average norm value {} is not equal to {}.",
            average,
            witness_norm,
        );
    }

    // Test lower bound verification
    #[test]
    fn test_lower_bound() {
        let (security_parameter, rank, multiplicity) = (128, 5, 1);

        let mut transcript = LabradorTranscript::new(ShakeSponge::default());
        transcript.set_u1(RqVector::random(&mut rng(), 1));
        let projections = transcript.generate_projections(security_parameter, rank, multiplicity);
        let witness = Witness::new(rank, multiplicity, 6400).s;

        let beta = witness[0].l2_norm_squared();
        // Check if the norm of the projection is bigger than 30 * (squared norm of the projection of the random polynomial)
        assert!(Projection::verify_projection_lower_bound(
            &projections.compute_projection(0, &witness[0]),
            beta
        ));
    }
}
