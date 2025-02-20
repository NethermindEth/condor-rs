use crate::rq::Rq;
use crate::zq::Zq;
use rand::rng;
//use labrador::rq::Rq;
use rand::prelude::*;

/// Projection matrix with values in {1,0,-1}
pub struct ProjectionMatrix {
    matrix: Vec<Vec<Zq>>, // 256 x (nxd) matrix, entries -1 or 1 in Zq
}

impl ProjectionMatrix {
    /// Defines a matrix of size 256xnd
    pub fn new(N: usize) -> Self {
        let mut matrix = vec![vec![Zq::zero(); N]; 256];
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
            concatenated_coeffs.extend_from_slice(&coeffs); // Extend the Vec with the coefficients
        }

        concatenated_coeffs
    }

    /// Calculates Projection  
    pub fn new(matrix: &ProjectionMatrix, s_i: &Vec<Rq<D>>) -> Self {
        let mut projection = vec![Zq::zero(); 256];
        let coefficients = Self::concatenate_coefficients(s_i.clone());
        for i in 0..256 {
            projection[i] = matrix.get_matrix()[i]
                .iter()
                .zip(coefficients.iter())
                .map(|(m, s)| *m * *s)
                .sum();
        }
        ProjectionVector { projection }
    }
    /// Obtain projection
    pub fn get_projection(&self) -> &Vec<Zq> {
        &self.projection
    }
}

/// Function to generate a random polynomial
/// Returns a vector of `n` random polynomials, each of degree `d`
pub fn generate_random_polynomials<R: Rng, const D: usize>(
    n: usize,
    rng: &mut R,
    beta: f64,
) -> Vec<Rq<D>> {
    let mut polynomials = Vec::with_capacity(n);

    for _ in 0..n {
        loop {
            // Generate random coefficients
            let coeffs: [Zq; D] = std::array::from_fn(|_| {
                let small_value: u32 = rng.random_range(0..100 as u32); // random value between 0 and 100 (as u322)
                match small_value.try_into() {
                    Ok(value) => Zq::new(value),
                    Err(_) => {
                        eprintln!("Failed to convert value: {}", small_value);
                        Zq::new(0) // Default value in case of failure
                    }
                }
            });

            // Create the polynomial
            let polynomial = Rq::new(coeffs);

            // Check the norm of the polynomial (assuming norm is some value you can compute)
            let norm = compute_norm(&polynomial); // This function should compute the norm of the polynomial

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
fn compute_norm<const D: usize>(polynomial: &Rq<D>) -> f64 {
    polynomial
        .get_coefficients()
        .iter()
        .map(|&coeff| coeff.to_f64())
        .sum()
}
