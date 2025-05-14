use rand::rng;

use crate::commitments::ajtai_commitment::AjtaiScheme;
use crate::core::env_params::EnvironmentParameters;
use crate::ring::rq::Rq;
use crate::ring::rq_matrix::RqMatrix;

pub struct AjtaiInstances {
    // A \in R_q^(k*n)
    pub commitment_scheme_a: AjtaiScheme,
    // B_{ik} \in R_q^(k_1*k), for i \in [1,r] and k \in [0, t_1-1]
    pub commitment_scheme_b: AjtaiScheme,
    // C_{ijk} \in R_q^(k_2*1), for i \in [1,r], j \in [i, r], and k \in [0, t_2-1]
    pub commitment_scheme_c: AjtaiScheme,
    // D_{ijk} \in R_q^(k_2*1), for i \in [1,r], j \in [i, r], and k \in [0, t_1-1]
    pub commitment_scheme_d: AjtaiScheme,
}

impl AjtaiInstances {
    pub fn new(ep: &EnvironmentParameters) -> Self {
        Self {
            commitment_scheme_a: AjtaiScheme::new(
                ep.beta,
                ep.beta,
                Self::challenge_rq_matrix(ep.kappa, ep.n),
            )
            .unwrap(),
            commitment_scheme_b: AjtaiScheme::new(
                ep.gamma_1,
                ep.gamma_1,
                Self::challenge_rq_matrix(ep.kappa_1, ep.r * ep.t_1 * ep.kappa),
            )
            .unwrap(),
            // Todo: gamma_1 should be changed to a valid witness bound
            commitment_scheme_c: AjtaiScheme::new(
                ep.gamma_1,
                ep.gamma_1,
                Self::challenge_rq_matrix(ep.kappa_1, ep.t_2 * ((ep.r.pow(2)) + ep.r) / 2),
            )
            .unwrap(),
            commitment_scheme_d: AjtaiScheme::new(
                ep.gamma_2,
                ep.gamma_2,
                Self::challenge_rq_matrix(ep.kappa_2, ep.t_1 * ((ep.r.pow(2)) + ep.r) / 2),
            )
            .unwrap(),
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
        let _pp = AjtaiInstances::new(&ep_1);

        println!(
            "Total time for PublicPrams::new: {:?}",
            total_start.elapsed()
        );
    }
}
