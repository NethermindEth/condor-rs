use crate::rq::Rq;
use crate::zq::Zq;
use rand::prelude::*;
use rand::rng;

/// Projection matrix with values in {1,0,-1} mod q
pub struct ProjectionMatrix<const D: usize> {
    matrix: Vec<Vec<Zq>>,
}

impl<const D: usize> ProjectionMatrix<D> {
    /// Defines a matrix of size 256xnxD
    /// n is the size of the vector of polynomials
    pub fn new(n: usize) -> Self {
        let mut matrix = vec![vec![Zq::ZERO; n * D]; 256];
        let mut rng = rng();
        // Fill the matrix with random values from {-1, 0, 1}
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
#[derive(Debug)]
pub struct ProjectionVector<const D: usize> {
    projection: [Zq; 256], // 256-dimensional projection vector
}

impl<const D: usize> ProjectionVector<D> {
    /// Calculates Projection  
    pub fn new(matrix: &ProjectionMatrix<D>, s_i: &[Rq<D>]) -> Self {
        let mut projection = [Zq::ZERO; 256];
        let coefficients = Self::concatenate_coefficients(s_i);
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
    pub fn get_projection(&self) -> &[Zq; 256] {
        &self.projection
    }

    /// Function to concatenate coefficients from multiple Rq into a Vec<Zq>
    fn concatenate_coefficients(rqvect: &[Rq<D>]) -> Vec<Zq> {
        let mut concatenated_coeffs: Vec<Zq> = Vec::new();

        // Iterate over each Rq, extracting the coefficients and concatenating them
        for rq in rqvect {
            let coeffs = rq.get_coefficients();
            concatenated_coeffs.extend_from_slice(&coeffs);
        }

        concatenated_coeffs
    }

    /// Euclidean norm
    pub fn norm_squared(&self) -> Zq {
        self.projection.iter().map(|coeff| *coeff * *coeff).sum()
    }
}

// Function to verify upper bound of projection
pub fn verify_upper_bound<const D: usize>(
    projection: ProjectionVector<D>,
    beta_squared: Zq,
) -> bool {
    projection.norm_squared().value() < (Zq::new(128) * beta_squared).value()
}
// Function to verify lower bound of projection
pub fn verify_lower_bound<const D: usize>(
    projection: ProjectionVector<D>,
    beta_squared: Zq,
) -> bool {
    projection.norm_squared().value() > (Zq::new(30) * beta_squared).value()
}

#[cfg(test)]
mod tests {
    use super::*;

    // Test Size correctness of projection matrix
    #[test]
    fn test_size_projection_matrix() {
        let n = 10;
        let matrix = ProjectionMatrix::<4>::new(n);

        assert_eq!(matrix.matrix.len(), 256, "Matrix should have 256 rows");
        assert_eq!(
            matrix.matrix[0].len(),
            n * 4,
            "Matrix should have n * D columns"
        );

        let n2 = 1;
        let matrix = ProjectionMatrix::<4>::new(n2);

        assert_eq!(matrix.matrix.len(), 256, "Matrix should have 256 rows");
        assert_eq!(
            matrix.matrix[0].len(),
            n2 * 4,
            "Matrix should have n * D columns"
        );
    }

    // Test the distribution of values in the random matrix
    #[test]
    fn test_random_distribution_matrix() {
        let n = 1000;
        let matrix = ProjectionMatrix::<4>::new(n);
        let mut counts = [0.0, 0.0, 0.0]; // -1, 0, 1
        for row in matrix.matrix.iter() {
            for &elem in row.iter() {
                if elem == Zq::MAX {
                    counts[0] += 1.0;
                } else if elem == Zq::zero() {
                    counts[1] += 1.0;
                } else if elem == Zq::one() {
                    counts[2] += 1.0;
                }
            }
        }
        let total: f64 = 1024000.0;
        let expected = [0.25, 0.5, 0.25];
        for i in 0..3 {
            let actual = counts[i] / total;
            println!("This is the actual value {}", actual);
            assert!(
                (actual - expected[i]).abs() < 0.05,
                "Values are not within expected proportions"
            ); // we allow some tolerance
        }
    }

    // Test vector concatenation
    #[test]
    fn test_vector_concatenation() {
        let polynomials = vec![
            vec![Zq::new(1), Zq::zero(), Zq::zero(), Zq::MAX].into(),
            vec![Zq::new(6), Zq::zero(), Zq::new(5), Zq::new(3)].into(),
        ];
        let vector = vec![
            Zq::new(1),
            Zq::zero(),
            Zq::zero(),
            Zq::MAX,
            Zq::new(6),
            Zq::zero(),
            Zq::new(5),
            Zq::new(3),
        ];
        assert!(ProjectionVector::<4>::concatenate_coefficients(&polynomials) == vector);
    }

    // Test that the probability of the inequality being true is close to 1/2
    #[test]
    fn test_probability_is_close_to_half() {
        let trials: f64 = 1000.0;
        let mut success_count: f64 = 0.0;
        let n = 5;
        for _ in 0..1000 {
            // Generate the random polynomials
            let polynomials = Rq::<64>::random_small_vector(n);
            // Generate projection matrix
            let matrix = ProjectionMatrix::new(n);
            // Generate Projection
            let projection = ProjectionVector::new(&matrix, &polynomials);
            let beta = Rq::compute_norm_squared(&polynomials);
            // Check if the norm of the projection is smaller than 128 * (squared norm of the projection of the random polynomial)
            let test: bool = verify_upper_bound(projection, beta);
            if test {
                success_count += 1.0;
            }
        }

        let observed_probability = success_count / trials;

        // we allow some tolerance
        let tolerance = 0.05;
        assert!(
            (observed_probability - 0.5).abs() < tolerance,
            "Observed probability {} is not close to 0.5",
            observed_probability
        );
    }

    // On average the projected norm squared is the same as 128 * vector norm squared
    #[test]
    fn average_value() {
        let trials: u32 = 5000;
        let n = 5;
        let polynomials = Rq::<64>::random_small_vector(n);
        let mut matrix = ProjectionMatrix::new(n);
        let mut projection = ProjectionVector::new(&matrix, &polynomials);
        let mut norm_sum = projection.norm_squared();
        let norm_value = (Zq::new(128) * Rq::compute_norm_squared(&polynomials)).value();
        // Run the test multiple times to simulate the probability
        for _ in 0..trials {
            matrix = ProjectionMatrix::new(n);
            projection = ProjectionVector::new(&matrix, &polynomials);
            norm_sum += projection.norm_squared();
        }

        // Calculate the observed probability
        let average = norm_sum.value() / trials;

        // we allow some tolerance
        let tolerance: u32 = 100;
        let abs_difference = if average < norm_value {
            norm_value - average
        } else {
            average - norm_value
        };

        assert!(
            abs_difference < tolerance,
            "Average norm value {} is not close to {}. Difference: {}",
            average,
            (Zq::new(128) * Rq::compute_norm_squared(&polynomials)).value(),
            abs_difference
        );
    }

    // Test lower bound verification
    #[test]
    fn test_lower_bound() {
        let n = 5;
        // Generate random vector of polynomials of small norm
        let polynomials = Rq::<64>::random_small_vector(n);
        // Generate projection matrix
        let matrix = ProjectionMatrix::new(n);
        // Generate Projection
        let projection = ProjectionVector::new(&matrix, &polynomials);
        let beta = Rq::compute_norm_squared(&polynomials);
        // Check if the norm of the projection is bigger than 30 * (squared norm of the projection of the random polynomial)
        assert!(verify_lower_bound(projection, beta));
    }
}
