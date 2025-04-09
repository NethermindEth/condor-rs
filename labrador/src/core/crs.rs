use crate::core::{challenge_set::ChallengeSet, env_params::EnvironmentParameters};
use crate::ring::poly::PolyVector;

pub struct PublicPrams {
    // A \in R_q^(k*n)
    pub matrix_a: Vec<PolyVector>,
    // B_{ik} \in R_q^(k_1*k), for i \in [1,r] and k \in [0, t_1-1]
    pub matrix_b: Vec<Vec<Vec<PolyVector>>>,
    // C_{ijk} \in R_q^(k_2*1), for i \in [1,r], j \in [i, r], and k \in [0, t_2-1]
    pub matrix_c: Vec<Vec<Vec<Vec<PolyVector>>>>,
    // D_{ijk} \in R_q^(k_2*1), for i \in [1,r], j \in [i, r], and k \in [0, t_1-1]
    pub matrix_d: Vec<Vec<Vec<Vec<PolyVector>>>>,
}

impl PublicPrams {
    pub fn new(ep: &EnvironmentParameters) -> Self {
        let d = ep.deg_bound_d;

        let matrix_a = Self::challenge_matrix(ep.k, ep.n, d);

        let matrix_b: Vec<Vec<Vec<PolyVector>>> = (0..ep.r)
            .map(|_| {
                (0..ep.t_1)
                    .map(|_| Self::challenge_matrix(ep.k_1, ep.k, d))
                    .collect()
            })
            .collect();

        let matrix_c: Vec<Vec<Vec<Vec<PolyVector>>>> = (0..ep.r)
            .map(|_| {
                (0..ep.r)
                    .map(|_| {
                        (0..ep.t_2)
                            .map(|_| Self::challenge_matrix(ep.k_2, 1, d))
                            .collect()
                    })
                    .collect()
            })
            .collect();

        let matrix_d: Vec<Vec<Vec<Vec<PolyVector>>>> = (0..ep.r)
            .map(|_| {
                (0..ep.r)
                    .map(|_| {
                        (0..ep.t_1)
                            .map(|_| Self::challenge_matrix(ep.k_2, 1, d))
                            .collect()
                    })
                    .collect()
            })
            .collect();

        Self {
            matrix_a,
            matrix_b,
            matrix_c,
            matrix_d,
        }
    }

    fn challenge_matrix(row: usize, col: usize, d: usize) -> Vec<PolyVector> {
        (0..row)
            .map(|_| {
                (0..col)
                    .map(|_| ChallengeSet::new(d).get_challenges().clone())
                    .collect()
            })
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
