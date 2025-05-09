use crate::ring::rq::Rq;
use crate::ring::rq_vector::RqVector;
use crate::ring::zq::Zq;
use rand::prelude::*;
use rand::rng;

/// The security level \lambda is 128
pub const SECURITY_LEVEL: usize = 128;
/// The size of the projection matrix is 2\lambda
pub const PROJECTION_MATRIX_SIZE: usize = 2 * SECURITY_LEVEL;
/// Projection matrix with values in {1,0,-1} mod q
pub struct ProjectionMatrix {
    matrix: Vec<Vec<Zq>>,
}

impl ProjectionMatrix {
    /// Defines a matrix of size 256xnxD
    /// n is the size of the vector of polynomials
    pub fn new(n: usize) -> Self {
        let mut matrix = vec![vec![Zq::ZERO; n * Rq::DEGREE]; PROJECTION_MATRIX_SIZE];
        let mut rng = rng();
        for row in matrix.iter_mut() {
            for elem in row.iter_mut() {
                let rand_val: f64 = rng.random();
                *elem = if rand_val < 0.25 {
                    Zq::MAX // -1 in Zq
                } else if rand_val < 0.75 {
                    Zq::ZERO
                } else {
                    Zq::ONE
                };
            }
        }

        ProjectionMatrix { matrix }
    }

    /// Returns the matrix
    pub fn get_matrix(&self) -> &Vec<Vec<Zq>> {
        &self.matrix
    }
}

/// Calculate projection vector
#[derive(Debug, Clone)]
pub struct ProjectionVector {
    projection: Vec<Zq>,
}

impl ProjectionVector {
    /// Calculates Projection
    pub fn new(matrix: &ProjectionMatrix, s_i: &RqVector) -> Self {
        let mut projection = vec![Zq::ZERO; PROJECTION_MATRIX_SIZE];
        let coefficients = s_i.concatenate_coefficients();
        for (i, item) in projection.iter_mut().enumerate() {
            *item = matrix.get_matrix()[i]
                .iter()
                .zip(coefficients.iter())
                .map(|(m, s)| *m * *s)
                .sum::<Zq>();
        }
        ProjectionVector { projection }
    }

    /// Obtain projection
    pub fn get_projection(&self) -> &Vec<Zq> {
        &self.projection
    }

    /// Euclidean norm
    pub fn norm_squared(&self) -> Zq {
        self.projection.iter().map(|coeff| *coeff * *coeff).sum()
    }
}

#[derive(Debug, Clone)]
pub struct Projections {
    pub projection: Vec<Zq>,
}

pub fn inner_product(first: &[Zq], second: &[Zq]) -> Zq {
    first
        .iter()
        .zip(second.iter())
        .map(|(&a, &b)| a * b)
        .fold(Zq::ZERO, |acc, x| acc + x)
}

impl Projections {
    /// Calculates Projection
    pub fn new(matrices: &[Vec<Vec<Zq>>], witness: &[RqVector]) -> Self {
        let mut temp = Vec::with_capacity(PROJECTION_MATRIX_SIZE);
        let coefficients: Vec<Vec<Zq>> = witness
            .iter()
            .map(|s_i| s_i.concatenate_coefficients())
            .collect();

        for j in 0..PROJECTION_MATRIX_SIZE {
            let mut sum = Zq::ZERO;
            for i in 0..matrices.len() {
                let pi = &matrices[i][j];
                let s_i = &coefficients[i];
                let inner = inner_product(pi, s_i);
                sum += inner;
            }
            temp.push(sum);
        }
        let projection = temp;

        Projections { projection }
    }

    /// Obtain projection
    pub fn get_projection(&self) -> &Vec<Zq> {
        &self.projection
    }

    /// Euclidean norm
    pub fn norm_squared(&self) -> Zq {
        self.projection.iter().map(|coeff| *coeff * *coeff).sum()
    }
}

// Bound verification

// LaBRADOR: Compact Proofs for R1CS from Module-SIS | Page 5 | Proving smallness section
const UPPER_BOUND_FACTOR: Zq = Zq::new(128);
const LOWER_BOUND_FACTOR: Zq = Zq::new(30);

// Function to verify upper bound of projection
pub fn verify_upper_bound(projection: ProjectionVector, beta_squared: Zq) -> bool {
    projection.norm_squared() < (UPPER_BOUND_FACTOR * beta_squared)
}
// Function to verify lower bound of projection
pub fn verify_lower_bound(projection: ProjectionVector, beta_squared: Zq) -> bool {
    projection.norm_squared() > (LOWER_BOUND_FACTOR * beta_squared)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Test Size correctness of projection matrix
    #[test]
    fn test_size_projection_matrix() {
        let n = 10;
        let matrix = ProjectionMatrix::new(n);

        assert_eq!(
            matrix.matrix.len(),
            PROJECTION_MATRIX_SIZE,
            "Matrix should have {} rows",
            PROJECTION_MATRIX_SIZE
        );
        assert_eq!(
            matrix.matrix[0].len(),
            n * Rq::DEGREE,
            "Matrix should have {} columns",
            n * Rq::DEGREE
        );

        let n2 = 1;
        let matrix = ProjectionMatrix::new(n2);

        assert_eq!(
            matrix.matrix.len(),
            PROJECTION_MATRIX_SIZE,
            "Matrix should have {} rows",
            PROJECTION_MATRIX_SIZE
        );
        assert_eq!(
            matrix.matrix[0].len(),
            n2 * Rq::DEGREE,
            "Matrix should have {} columns",
            n2 * Rq::DEGREE
        );
    }

    // Test the distribution of values in the random matrix
    #[test]
    #[allow(clippy::as_conversions)]
    fn test_random_distribution_matrix() {
        // 1000 was chosen to provide a reasonably large sample size
        let n = 1000;
        let matrix = ProjectionMatrix::new(n);
        let mut counts = [0.0, 0.0, 0.0]; // -1, 0, 1
        for row in matrix.matrix.iter() {
            for &elem in row.iter() {
                if elem == Zq::MAX {
                    counts[0] += 1.0;
                } else if elem == Zq::ZERO {
                    counts[1] += 1.0;
                } else if elem == Zq::ONE {
                    counts[2] += 1.0;
                }
            }
        }
        // Number of elements in the matrix as f64 (256x4x1000)
        let total: f64 = (256 * Rq::DEGREE * n) as f64;
        println!("this is the total amount of elements{}", total);
        let expected = [0.25, 0.5, 0.25];
        for i in 0..3 {
            let actual = counts[i] / total;
            println!("This is the actual value {}", actual);
            assert!(
                //Since its a statistical test some small error tolerance is allowed
                (actual - expected[i]).abs() < 0.005,
                "Values are not within expected proportions"
            );
        }
    }

    // Test that the probability of the inequality being true is close to 1/2
    #[test]
    #[cfg(not(feature = "skip-slow-tests"))]
    fn test_probability_is_close_to_half() {
        // 1000 was chosen to provide a reasonably large sample size
        let trials: f64 = 1000.0;
        let mut success_count: f64 = 0.0;
        let n = 5;
        for _ in 0..1000 {
            let mut rng = rng();
            // Generate the random polynomials
            let polynomials = RqVector::random_ternary(&mut rng, 5);
            // Generate projection matrix
            let matrix = ProjectionMatrix::new(n);
            // Generate Projection
            let projection = ProjectionVector::new(&matrix, &polynomials);
            let beta = polynomials.compute_norm_squared();
            // Check if the norm of the projection is smaller than 128 * (squared norm of the projection of the random polynomial)
            let test: bool = verify_upper_bound(projection, beta);
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
    fn average_value() {
        // 10000 was chosen to provide a reasonably large sample size
        let trials: u128 = 10000;
        let n = 3;
        let mut rng = rng();
        let polynomials = RqVector::random_ternary(&mut rng, 3);
        let mut matrix = ProjectionMatrix::new(n);
        let mut projection = ProjectionVector::new(&matrix, &polynomials);
        let mut norm_sum = projection.norm_squared();
        let norm_value = (Zq::new(SECURITY_LEVEL.try_into().unwrap())
            * RqVector::compute_norm_squared(&polynomials))
        .to_u128();
        // Run the test multiple times to simulate the probability
        for _ in 0..trials {
            matrix = ProjectionMatrix::new(n);
            projection = ProjectionVector::new(&matrix, &polynomials);
            norm_sum += projection.norm_squared();
        }

        // Calculate the observed probability
        let average = norm_sum.to_u128() / trials;
        let difference = if norm_value <= average {
            average - norm_value
        } else {
            norm_value - average
        };

        // we choose a small tolerance value for possible statistical error
        let tolerance: u128 = 20;
        assert!(
            difference < tolerance,
            "Average norm value {} is not equal to {}.",
            average,
            Zq::new(SECURITY_LEVEL.try_into().unwrap()) * polynomials.compute_norm_squared(),
        );
    }

    // Test lower bound verification
    #[test]
    fn test_lower_bound() {
        let n = 5;
        let mut rng = rng();
        // Generate random vector of polynomials of small norm
        let polynomials = RqVector::random_ternary(&mut rng, 5);
        // Generate projection matrix
        let matrix = ProjectionMatrix::new(n);
        // Generate Projection
        let projection = ProjectionVector::new(&matrix, &polynomials);
        let beta = polynomials.compute_norm_squared();
        // Check if the norm of the projection is bigger than 30 * (squared norm of the projection of the random polynomial)
        assert!(verify_lower_bound(projection, beta));
    }

    // Test vector concatenation
    #[test]
    fn test_vector_concatenation() {
        let mut vec1 = vec![Zq::new(10); Rq::DEGREE];
        let mut vec2 = vec![Zq::new(12); Rq::DEGREE];
        let polynomials = RqVector::from(vec![vec1.clone().into(), vec2.clone().into()]);

        vec1.append(&mut vec2);
        assert!(polynomials.concatenate_coefficients() == vec1);
    }
}
