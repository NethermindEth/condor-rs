use rand::rng;

use crate::core::env_params::EnvironmentParameters;
use crate::ring::rq::Rq;
use crate::ring::rq_matrix::RqMatrix;

#[derive(Clone)]
pub struct PublicPrams {
    // A \in R_q^(k*n)
    pub matrix_a: RqMatrix,
    // B_{ik} \in R_q^(k_1*k), for i \in [1,r] and k \in [0, t_1-1]
    pub matrix_b: RqMatrix,
    // C_{ijk} \in R_q^(k_2*1), for i \in [1,r], j \in [i, r], and k \in [0, t_2-1]
    pub matrix_c: RqMatrix,
    // D_{ijk} \in R_q^(k_2*1), for i \in [1,r], j \in [i, r], and k \in [0, t_1-1]
    pub matrix_d: RqMatrix,
}

impl PublicPrams {
    pub fn new(ep: &EnvironmentParameters) -> Self {
        let matrix_a = Self::challenge_rq_matrix(ep.kappa, ep.n);
        let matrix_b = Self::challenge_rq_matrix(ep.kappa_1, ep.r * ep.t_1 * ep.kappa);
        // Todo: Need to confirm the row length is k_1 or k_2. For now, it is set to be k_1.
        let matrix_c = Self::challenge_rq_matrix(ep.kappa_1, ep.t_2 * ((ep.r.pow(2)) + ep.r) / 2);
        let matrix_d = Self::challenge_rq_matrix(ep.kappa_2, ep.t_1 * ((ep.r.pow(2)) + ep.r) / 2);

        Self {
            matrix_a,
            matrix_b,
            matrix_c,
            matrix_d,
        }
    }

    fn challenge_rq_matrix(row: usize, col: usize) -> RqMatrix {
        (0..row)
            .map(|_| (0..col).map(|_| Rq::random(&mut rng())).collect())
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Instant;

    #[test]
    fn test_crs() {
        // set up example environment parameters, use default set for testing.
        let ep_1 = EnvironmentParameters::default();

        let total_start = Instant::now();
        // generate the common reference string matrices A, B, C, D
        let _pp = PublicPrams::new(&ep_1);

        println!(
            "Total time for PublicPrams::new: {:?}",
            total_start.elapsed()
        );
    }
}
