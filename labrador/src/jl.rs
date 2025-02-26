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
        let mut matrix = vec![vec![Zq::zero(); n * D]; 256];
        let mut rng = rng();
        // Fill the matrix with random values from {-1, 0, 1}
        for row in matrix.iter_mut() {
            for elem in row.iter_mut() {
                let rand_val: f64 = rng.random();
                *elem = if rand_val < 0.25 {
                    Zq::new(u32::MAX) // -1 in Zq
                } else if rand_val < 0.75 {
                    Zq::zero()
                } else {
                    Zq::one()
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
    /// Function to concatenate coefficients from multiple Rq into a Vec<Zq>
    fn concatenate_coefficients(rqs: Vec<Rq<D>>) -> Vec<Zq> {
        let mut concatenated_coeffs: Vec<Zq> = Vec::new();

        // Iterate over each Rq, extracting the coefficients and concatenating them
        for rq in rqs {
            let coeffs = rq.get_coefficients();
            concatenated_coeffs.extend_from_slice(&coeffs);
        }

        concatenated_coeffs
    }
    /// Euclidean norm
    pub fn norm(&self) -> Zq {
        self.projection.iter().map(|coeff| *coeff * *coeff).sum()
    }

    /// Calculates Projection  
    pub fn new(matrix: &ProjectionMatrix<D>, s_i: &[Rq<D>]) -> Self {
        let mut projection = [Zq::zero(); 256];
        let coefficients = Self::concatenate_coefficients(s_i.to_vec());
        for (i, item) in projection.iter_mut().enumerate().take(256) {
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
}

/// Returns a vector of `n` random polynomials, each of degree `d`, for testing purposes.
pub fn generate_random_polynomials<R: Rng, const D: usize>(
    n: usize,
    rng: &mut R,
    beta: Zq, // vector norm bound
) -> Vec<Rq<D>> {
    let mut polynomials = Vec::with_capacity(n);

    for _ in 0..n {
        loop {
            // Generate random coefficients
            let coeffs: [Zq; D] = std::array::from_fn(|_| {
                // Small values so we get a vector of 'small norm'
                let small_value: u32 = rng.random_range(0..100);
                Zq::new(small_value)
            });

            // Create the polynomial
            let polynomial = Rq::new(coeffs);

            // Compute the norm of all polynomials including the new one
            let mut all_polynomials = polynomials.clone();
            all_polynomials.push(polynomial.clone());
            let norm = Rq::compute_norm_squared(&all_polynomials);

            // If the norm is smaller than beta, add the polynomial to the list
            if norm.value() < beta.value() {
                polynomials.push(polynomial);
                break;
            }
        }
    }

    polynomials
}

/// Returns a boolean if the norm of the projection is smaller than 128 for a random polynomial
fn _smaller_than_128b<const D: usize>(n: usize) -> bool {
    // Generate the random polynomials
    let polynomials = Rq::<D>::random_small_vector(n);
    // Generate projection matrix
    let matrix = ProjectionMatrix::new(n);
    // Generate Projection
    let projection = ProjectionVector::new(&matrix, &polynomials);
    // Check if the norm of the projection is smaller than 128 * (projection of the random polynomial)
    let result =
        projection.norm().value() < (Zq::new(128) * Rq::compute_norm_squared(&polynomials)).value();
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    // Test that the probability of the inequality being true is close to 1/2
    fn test_probability_is_close_to_half() {
        let trials: f64 = 1000.0;
        let mut success_count: f64 = 0.0;
        for _ in 0..1000 {
            if _smaller_than_128b::<64>(5) {
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

    #[test]

    // on average the projected norm squared is the same as 128 * vector norm squared
    fn average_value() {
        let trials: u32 = 5000;
        let n = 5;
        let polynomials = Rq::<64>::random_small_vector(n);
        let mut matrix = ProjectionMatrix::new(n);
        let mut projection = ProjectionVector::new(&matrix, &polynomials);
        let mut norm_sum = projection.norm();
        let norm_value = (Zq::new(128) * Rq::compute_norm_squared(&polynomials)).value();
        // Run the test multiple times to simulate the probability
        for _ in 0..trials {
            matrix = ProjectionMatrix::new(n);
            projection = ProjectionVector::new(&matrix, &polynomials);
            norm_sum += projection.norm();
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
}

// define the dot product function
// implement verification of norm bounds for a guiven polynomial vector
