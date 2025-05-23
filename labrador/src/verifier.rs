#![allow(clippy::result_large_err)]

use thiserror::Error;

use crate::commitments::common_instances::AjtaiInstances;
use crate::commitments::outer_commitments::{self, DecompositionParameters, OuterCommitment};
use crate::core::aggregate::{self, FunctionsAggregation, ZeroConstantFunctionsAggregation};
use crate::core::{env_params::EnvironmentParameters, statement::Statement};
use crate::ring::rq::Rq;
use crate::ring::rq_matrix::RqMatrix;
use crate::ring::rq_vector::RqVector;
use crate::ring::zq::{Zq, ZqVector};
use crate::transcript::{LabradorTranscript, Sponge};

#[derive(Debug, Error)]
pub enum VerifierError {
    #[error("matrix not symmetric at ({i},{j}): expected {expected:?}, found {found:?}")]
    NotSymmetric {
        i: usize,
        j: usize,
        expected: Rq,
        found: Rq,
    },
    #[error("B0 mismatch at index {index}: expected {expected}, computed {computed}")]
    B0Mismatch {
        index: usize,
        expected: Zq,
        computed: Zq,
    },
    #[error("‖z‖² = {norm} exceeds allowed bound {allowed}")]
    NormSumExceeded { norm: Zq, allowed: Zq },
    #[error("A·z check failed: expected {expected:?}, computed {computed:?}")]
    AzError {
        computed: RqVector,
        expected: RqVector,
    },
    #[error("⟨z,z⟩ mismatch: expected {expected:?}, computed {computed:?}")]
    ZInnerError { computed: Rq, expected: Rq },
    #[error("φ(z) mismatch: expected {expected:?}, computed {computed:?}")]
    PhiError { computed: Rq, expected: Rq },
    #[error("relation check failed")]
    RelationCheckFailed,
    #[error("outer commitment mismatch: expected {expected:?}, computed {computed:?}")]
    OuterCommitError {
        computed: RqVector,
        expected: RqVector,
    },
    #[error(transparent)]
    DecompositionError(#[from] outer_commitments::DecompositionError),
}
pub struct LabradorVerifier<'a> {
    pub pp: &'a AjtaiInstances,
    pub st: &'a Statement,
}

impl<'a> LabradorVerifier<'a> {
    pub fn new(pp: &'a AjtaiInstances, st: &'a Statement) -> Self {
        Self { pp, st }
    }

    /// All check conditions are from page 18
    pub fn verify<S: Sponge>(
        &mut self,
        proof: &LabradorTranscript<S>,
        ep: &EnvironmentParameters,
    ) -> Result<bool, VerifierError> {
        let mut transcript = LabradorTranscript::new(S::default());
        let mut constant_aggregation = ZeroConstantFunctionsAggregation::new(ep);
        let mut funcs_aggregator = FunctionsAggregation::new(ep);

        transcript.absorb_u1(&proof.u1);
        let projections =
            transcript.generate_projections(ep.security_parameter, ep.rank, ep.multiplicity);
        transcript.absorb_vector_p(&proof.vector_p);
        let size_of_psi = usize::div_ceil(ep.security_parameter, ep.log_q);
        let size_of_omega = size_of_psi;
        let psi = transcript.generate_vector_psi(size_of_psi, ep.constraint_l);
        let omega = transcript.generate_vector_omega(size_of_omega, ep.security_parameter);
        transcript.absorb_vector_b_ct_aggr(&proof.b_ct_aggr);
        let vector_alpha = transcript.generate_rq_vector(ep.constraint_k);
        let size_of_beta = size_of_psi;
        let vector_beta = transcript.generate_rq_vector(size_of_beta);
        transcript.absorb_u2(&proof.u2);
        let challenges = transcript.generate_challenges(ep.operator_norm, ep.multiplicity);

        // check b_0^{''(k)} ?= <omega^(k),p> + \sum(psi_l^(k) * b_0^{'(l)})
        Self::check_b_0_aggr(self, proof, ep, &psi, &omega)?;

        // 3. line 14: check norm_sum(z, t, g, h) <= (beta')^2

        // decompose z into z = z^(0) + z^(1) * b, only two parts.
        let z_ij = RqVector::decompose(&proof.z, ep.b, 2);
        let t_ij: Vec<Vec<RqVector>> = proof
            .t
            .get_elements()
            .iter()
            .map(|i| RqVector::decompose(i, ep.b, ep.t_1))
            .collect();
        let g_ij: Vec<Vec<RqVector>> = proof
            .g
            .get_elements()
            .iter()
            .map(|i| RqVector::decompose(i, ep.b, ep.t_2))
            .collect();
        let h_ij: Vec<Vec<RqVector>> = proof
            .h
            .get_elements()
            .iter()
            .map(|i| RqVector::decompose(i, ep.b, ep.t_1))
            .collect();
        let norm_z_ij = z_ij
            .iter()
            .fold(Zq::ZERO, |acc, p| acc + p.compute_norm_squared());
        let norm_t_ij = Self::norm_squared(&t_ij);
        let norm_g_ij = Self::norm_squared(&g_ij);
        let norm_h_ij = Self::norm_squared(&h_ij);
        let norm_sum = norm_z_ij + norm_t_ij + norm_g_ij + norm_h_ij;

        if norm_sum > ep.beta * ep.beta {
            return Err(VerifierError::NormSumExceeded {
                norm: norm_sum,
                allowed: ep.beta * ep.beta,
            });
        }

        // 4. line 15: check Az ?= c_1 * t_1 + ... + c_r * t_r
        let az = self.pp.commitment_scheme_a.matrix() * &proof.z;
        let ct_sum = aggregate::compute_linear_combination(
            proof.t.get_elements(),
            challenges.get_elements(),
        );
        if az != ct_sum {
            return Err(VerifierError::AzError {
                computed: az,
                expected: ct_sum,
            });
        }

        // 5. lne 16: check <z, z> ?= \sum(g_ij * c_i * c_j)

        let z_inner = proof.z.inner_product_poly_vector(&proof.z);
        let sum_gij_cij = Self::calculate_gh_ci_cj(&proof.g, &challenges, ep.multiplicity);

        if z_inner != sum_gij_cij {
            return Err(VerifierError::ZInnerError {
                computed: z_inner,
                expected: sum_gij_cij,
            });
        }

        // 6. line 17: check \sum(<\phi_i, z>c_i) ?= \sum(h_ij * c_i * c_j)
        constant_aggregation.calculate_agg_phi_double_prime(
            &self.st.phi_ct,
            &projections.get_conjugated_projection_matrices(),
            &psi,
            &omega,
        );
        funcs_aggregator.calculate_aggr_phi(
            &self.st.phi_constraint,
            constant_aggregation.get_phi_double_prime(),
            &vector_alpha,
            &vector_beta,
        );
        let sum_phi_z_c =
            Self::calculate_phi_z_c(funcs_aggregator.get_appr_phi(), &challenges, &proof.z);
        let sum_hij_cij = Self::calculate_gh_ci_cj(&proof.h, &challenges, ep.multiplicity);

        // Left side multiple by 2 because of when we calculate h_ij, we didn't apply the division (divided by 2)
        if &sum_phi_z_c * Zq::TWO != sum_hij_cij {
            return Err(VerifierError::PhiError {
                computed: &sum_phi_z_c * Zq::TWO,
                expected: sum_hij_cij,
            });
        }

        // 7. line 18: check \sum(a_ij * g_ij) + \sum(h_ii) - b ?= 0

        constant_aggregation.calculate_agg_a_double_prime(&psi, &self.st.a_ct);
        funcs_aggregator.calculate_agg_a(
            &self.st.a_constraint,
            constant_aggregation.get_alpha_double_prime(),
            &vector_alpha,
            &vector_beta,
        );

        funcs_aggregator.calculate_aggr_b(
            &self.st.b_constraint,
            &proof.b_ct_aggr,
            &vector_alpha,
            &vector_beta,
        );

        if !Self::check_relation(
            funcs_aggregator.get_agg_a(),
            funcs_aggregator.get_aggr_b(),
            &proof.g,
            &proof.h,
        ) {
            return Err(VerifierError::RelationCheckFailed);
        }

        // 8. line 19: u_1 ?= \sum(\sum(B_ik * t_i^(k))) + \sum(\sum(C_ijk * g_ij^(k)))

        let u_1 = &proof.u1;
        let mut outer_commitments = OuterCommitment::new(self.pp);
        let commitment_u1 = outer_commitments.compute_u1(
            &proof.t,
            DecompositionParameters::new(ep.b, ep.t_1)?,
            &proof.g,
            DecompositionParameters::new(ep.b, ep.t_2)?,
        );

        if proof.u1 != commitment_u1 {
            return Err(VerifierError::OuterCommitError {
                computed: u_1.clone(),
                expected: commitment_u1,
            });
        }

        // 9. line 20: u_2 ?= \sum(\sum(D_ijk * h_ij^(k)))
        let commitment_u2 =
            outer_commitments.compute_u2(&proof.h, DecompositionParameters::new(ep.b, ep.t_1)?);

        if proof.u2 != commitment_u2 {
            return Err(VerifierError::OuterCommitError {
                computed: commitment_u2,
                expected: proof.u2.clone(),
            });
        }

        Ok(true)
    }

    /// calculate the right hand side of line 16 or line 17, \sum(g_ij * c_i * c_j) or \sum(h_ij * c_i * c_j)
    fn calculate_gh_ci_cj(x_ij: &RqMatrix, random_c: &RqVector, r: usize) -> Rq {
        (0..r)
            .map(|i| {
                (0..r)
                    .map(|j| {
                        &(x_ij.get_cell(i, j) * &random_c.get_elements()[i])
                            * &random_c.get_elements()[j]
                    })
                    .fold(Rq::zero(), |acc, x| &acc + &x)
            })
            .fold(Rq::zero(), |acc, x| &acc + &x)
    }

    /// calculate the left hand side of line 17, \sum(<\phi_z, z> * c_i)
    fn calculate_phi_z_c(phi: &[RqVector], c: &RqVector, z: &RqVector) -> Rq {
        phi.iter()
            .zip(c.iter())
            .map(|(phi_i, c_i)| &(phi_i.inner_product_poly_vector(z)) * c_i)
            .fold(Rq::zero(), |acc, x| &acc + &x)
    }

    fn norm_squared(polys: &[Vec<RqVector>]) -> Zq {
        polys.iter().fold(Zq::ZERO, |acc, poly| {
            acc + poly
                .iter()
                .fold(Zq::ZERO, |acc, p| acc + p.compute_norm_squared())
        })
    }

    /// line 18, page 18: check if \sum(a_{ij} * g_{ij}) + \sum(h_{ii}) - b ?= 0
    /// in the verifier process, page 18 from the paper.
    ///
    /// param: a_primes: a_{ij}^{''(k)}
    /// param: b_primes: b^{''(k)}
    /// param: g: g_{ij}
    /// param: h: h_{ii}
    ///
    /// return: true if the relation holds, false otherwise
    pub fn check_relation(a_primes: &RqMatrix, b_primes: &Rq, g: &RqMatrix, h: &RqMatrix) -> bool {
        let r = a_primes.get_elements().len();

        let mut sum_a_primes_g = Rq::zero();
        // walk only over the stored half: i ≤ j
        for i in 0..r {
            for j in 0..r {
                sum_a_primes_g = &sum_a_primes_g + &(a_primes.get_cell(i, j) * g.get_cell(i, j));
            }
        }

        let sum_h_ii = (0..r).fold(Rq::zero(), |acc, i| &acc + h.get_cell(i, i));

        let b_primes2 = b_primes * Zq::TWO;
        let sum_a_primes_g2 = &sum_a_primes_g * Zq::TWO;

        &sum_a_primes_g2 + &sum_h_ii == b_primes2
    }

    fn check_b_0_aggr<S: Sponge>(
        &self,
        proof: &LabradorTranscript<S>,
        ep: &EnvironmentParameters,
        psi: &[Vec<Zq>],
        omega: &[Vec<Zq>],
    ) -> Result<bool, VerifierError> {
        for k in 0..ep.kappa {
            let b_0_poly = proof.b_ct_aggr.get_elements()[k].get_coefficients()[0];
            let mut b_0: Zq = (0..ep.constraint_l)
                .map(|l| psi[k][l] * self.st.b_0_ct[l])
                .sum();

            let inner_omega_p = omega[k].inner_product(&proof.vector_p);
            b_0 += inner_omega_p;
            if b_0 != b_0_poly {
                return Err(VerifierError::B0Mismatch {
                    index: k,
                    expected: b_0_poly,
                    computed: b_0,
                });
            }
        }

        Ok(true)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prover::{LabradorProver, Witness};
    use crate::transcript::sponges::shake::ShakeSponge;

    #[test]
    fn test_verify() {
        // set up example environment, use default set for testing.
        let ep_1 = EnvironmentParameters::default();
        // generate a random witness based on ep above
        let witness_1 = Witness::new(&ep_1);
        // generate public statements based on witness_1
        let st: Statement = Statement::new(&witness_1, &ep_1);
        // generate the common reference string matrices
        let pp = AjtaiInstances::new(&ep_1);

        // create a new prover
        let mut prover = LabradorProver::new(&pp, &witness_1, &st);
        let proof: LabradorTranscript<ShakeSponge> = prover.prove(&ep_1).unwrap();

        // create a new verifier
        let mut verifier = LabradorVerifier::new(&pp, &st);
        let result = verifier.verify(&proof, &ep_1);
        assert!(result.unwrap());
    }
}
