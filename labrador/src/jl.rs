use crate::rq::Rq;
use crate::zq::Zq;
use rand::prelude::*;
use rand::rng;

/// Projection matrix with values in {1,0,-1} mod q
pub struct ProjectionMatrix {
    matrix: Vec<Vec<Zq>>,
}

impl ProjectionMatrix {
    /// Defines a matrix of size 256xn
    pub fn new(n: usize) -> Self {
        let mut matrix = vec![vec![Zq::zero(); n]; 256];
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
    projection: Vec<Zq>, // 256-dimensional projection vector
}

impl<const D: usize> ProjectionVector<D> {
    // Function to concatenate coefficients from multiple Rq into a Vec<Zq>
    fn concatenate_coefficients(rqs: Vec<Rq<D>>) -> Vec<Zq> {
        let mut concatenated_coeffs: Vec<Zq> = Vec::new();

        // Iterate over each Rq, extracting the coefficients and concatenating them
        for rq in rqs {
            let coeffs = rq.get_coefficients();
            concatenated_coeffs.extend_from_slice(&coeffs);
        }

        concatenated_coeffs
    }
    // Euclidean norm
    pub fn norm(&self) -> f64 {
        self.projection
            .iter()
            .map(|coeff| coeff.to_f64().powi(2))
            .sum::<f64>()
            .sqrt()
    }

    /// Calculates Projection  
    pub fn new(matrix: &ProjectionMatrix, s_i: &[Rq<D>]) -> Self {
        let mut projection = vec![Zq::zero(); 256];
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
    pub fn get_projection(&self) -> &Vec<Zq> {
        &self.projection
    }
}

/// Returns a vector of `n` random polynomials, each of degree `d`, for testing purposes.
pub fn generate_random_polynomials<R: Rng, const D: usize>(
    n: usize,
    rng: &mut R,
    beta: f64, // vector norm bound
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
            let norm = compute_norm(&all_polynomials);

            // If the norm is smaller than beta, add the polynomial to the list
            if norm < beta {
                polynomials.push(polynomial);
                break; // Exit the loop as the polynomial is valid
            }
            // If the norm is greater than or equal to beta, try again with a new random polynomial
        }
    }

    polynomials
}

// Helper function to compute the norm of the polynomial
pub fn compute_norm<const D: usize>(polynomials: &[Rq<D>]) -> f64 {
    polynomials
        .iter()
        .flat_map(|poly| poly.get_coefficients()) // Collect coefficients from all polynomials
        .map(|coeff| coeff.to_f64().powi(2))
        .sum::<f64>()
        .sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::ThreadRng;
    use rand::Rng;

    // Norm size error
    #[test]
    fn norm_distance() {
        for _ in 0..5 {
            let mut rng = rand::rngs::ThreadRng::default();

            // Generate random values for n and beta at runtime
            const D: usize = 4;
            let n: usize = rng.random_range(3..5); // Random vector size between 3 and 10
            let beta: f64 = rng.random_range(200.0..500.0); // 'Small' Random value for beta

            // Generate random polynomials using d and n as runtime values
            let polynomials = generate_random_polynomials::<ThreadRng, D>(n, &mut rng, beta);
            let matrix = ProjectionMatrix::new(D * n);
            let projection = ProjectionVector::new(&matrix, &polynomials);

            assert!(
                projection.norm() > 128.0_f64.sqrt() * compute_norm(&polynomials),
                "This error message implies the Modular Johnson-Lindenstrauss Lemma is working"
            );
        }
    }
    // Average of norms error
    #[test]
    fn average_norm_distance() {
        let mut rng = rand::rngs::ThreadRng::default();
        const D: usize = 4;
        let n: usize = 3;
        let beta: f64 = 500_f64;
        let polynomials = generate_random_polynomials::<ThreadRng, D>(n, &mut rng, beta);
        let mut norm_sum = 0.0;

        for _ in 0..50 {
            // Generate random projections
            let matrix = ProjectionMatrix::new(D * n);
            let projection = ProjectionVector::new(&matrix, &polynomials);
            // Sum up the norms
            norm_sum += projection.norm();
        }

        // Compute the average of the norms
        let average_norm = norm_sum / 50.0;
        // Assert that the average norm is greater than the expected value
        assert!(
            average_norm > 128.0_f64.sqrt() * compute_norm(&polynomials),
            "This error message implies the Modular Johnson-Lindenstrauss Lemma is working"
        );
    }
}
