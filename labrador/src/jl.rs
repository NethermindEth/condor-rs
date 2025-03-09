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
        for row in matrix.iter_mut() {
            //let coefficients = Rq::<D>::random_ternary(&mut rng).get_coefficients().to_vec();
            //let vector = Rq::<D>::random_small_vector(n);
            //let mut r: Vec<Zq> = Vec::new();
            //for polynomial in vector{
            //    for element in polynomial.get_coefficients(){
            //        r.push(*element);
            //    }
            //}
            //*row = r;
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
        let total_coeffs = rqvect.len() * D;
        let mut concatenated_coeffs: Vec<Zq> = Vec::with_capacity(total_coeffs);
        // Iterate over each Rq, extracting the coefficients and concatenating them
        for rq in rqvect {
            let coeffs = rq.get_coefficients();
            concatenated_coeffs.extend_from_slice(coeffs);
        }

        concatenated_coeffs
    }

    /// Euclidean norm
    pub fn norm_squared(&self) -> Zq {
        self.projection.iter().map(|coeff| *coeff * *coeff).sum()
    }
}
// Bound verification
const UPPER_BOUND_FACTOR: Zq = Zq::new(128);
const LOWER_BOUND_FACTOR: Zq = Zq::new(30);
// Function to verify upper bound of projection
pub fn verify_upper_bound<const D: usize>(
    projection: ProjectionVector<D>,
    beta_squared: Zq,
) -> bool {
    projection.norm_squared().value() < (UPPER_BOUND_FACTOR * beta_squared).value()
}
// Function to verify lower bound of projection
pub fn verify_lower_bound<const D: usize>(
    projection: ProjectionVector<D>,
    beta_squared: Zq,
) -> bool {
    projection.norm_squared().value() > (LOWER_BOUND_FACTOR * beta_squared).value()
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
        // 1000 was chosen to provide a reasonably large sample size
        let n = 1000;
        let matrix = ProjectionMatrix::<4>::new(n);
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
        let total: f64 = 1024000.0;
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
    fn test_probability_is_close_to_half() {
        // 10.000 was chosen to provide a reasonably large sample size
        let trials: f64 = 10000.0;
        let mut success_count: f64 = 0.0;
        let n = 5;
        for _ in 0..10000 {
            // Generate the random polynomials
            let polynomials = Rq::<5>::random_small_vector(n);
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
    fn average_value() {
        // 100.000 was chosen to provide a reasonably large sample size
        let trials: u32 = 100000;
        let n = 3;
        let polynomials = Rq::<3>::random_small_vector(n);
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
        let difference = if norm_value <= average {
            average - norm_value
        } else {
            norm_value - average
        };

        // we choose a small tolerance value for possible statistical error
        let tolerance: u32 = 2;
        assert!(
            difference < tolerance,
            "Average norm value {} is not equal to {}.",
            average,
            (Zq::new(128) * Rq::compute_norm_squared(&polynomials)).value(),
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

    // Test vector concatenation
    #[test]
    fn test_vector_concatenation() {
        let polynomials = vec![
            vec![Zq::ONE, Zq::ZERO, Zq::ZERO, Zq::MAX].into(),
            vec![Zq::new(6), Zq::ZERO, Zq::new(5), Zq::new(3)].into(),
        ];
        let vector = vec![
            Zq::ONE,
            Zq::ZERO,
            Zq::ZERO,
            Zq::MAX,
            Zq::new(6),
            Zq::ZERO,
            Zq::new(5),
            Zq::new(3),
        ];
        assert!(ProjectionVector::<4>::concatenate_coefficients(&polynomials) == vector);
    }
}
