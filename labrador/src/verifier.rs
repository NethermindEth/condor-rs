#![allow(clippy::result_large_err)]

use thiserror::Error;

use crate::commitments::common_instances::AjtaiInstances;
use crate::commitments::outer_commitments::{self, DecompositionParameters};
use crate::core::aggregate::{FunctionsAggregation, ZeroConstantFunctionsAggregation};
use crate::core::inner_product;
use crate::core::jl::Projection;
use crate::relation::env_params;
use crate::relation::{env_params::EnvironmentParameters, statement::Statement};
use crate::ring::rq::Rq;
use crate::ring::rq_matrix::RqMatrix;
use crate::ring::rq_vector::RqVector;
use crate::ring::zq::Zq;
use crate::ring::Norms;
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
    #[error("Computed Norm: {norm} -- Maximum Bound {allowed}, Step: {step}")]
    NormSumExceeded {
        norm: u128,
        allowed: u128,
        step: String,
    },
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
    params: &'a EnvironmentParameters,
    crs: &'a AjtaiInstances,
    st: &'a Statement,
    // Aggregation instances
    constant_aggregator: ZeroConstantFunctionsAggregation<'a>,
    funcs_aggregator: FunctionsAggregation<'a>,
}

impl<'a> LabradorVerifier<'a> {
    pub fn new(
        params: &'a EnvironmentParameters,
        crs: &'a AjtaiInstances,
        st: &'a Statement,
    ) -> Self {
        Self {
            params,
            crs,
            st,
            constant_aggregator: ZeroConstantFunctionsAggregation::new(params),
            funcs_aggregator: FunctionsAggregation::new(params),
        }
    }

    /// All check conditions are from page 18
    pub fn verify<S: Sponge>(
        &mut self,
        proof: &LabradorTranscript<S>,
    ) -> Result<bool, VerifierError> {
        let mut transcript = LabradorTranscript::new(S::default());

        let projections = self.check_vector_p_norm_bound(proof, &mut transcript)?;

        // check b_0^{''(k)} ?= <omega^(k),p> + \sum(psi_l^(k) * b_0^{'(l)})
        let (psi, omega) = self.check_b_double_prime_constant(proof, &mut transcript)?;

        // 3. line 14: check norm_sum(z, t, g, h) <= (beta')^2
        self.check_final_norm_sum(proof)?;

        let vector_alpha = transcript.generate_rq_vector(self.params.constraint_k);
        let vector_beta = transcript.generate_rq_vector(usize::div_ceil(
            env_params::SECURITY_PARAMETER,
            self.params.log_q,
        ));
        transcript.absorb_u2(&proof.u2);

        // 4. line 15: check Az ?= c_1 * t_1 + ... + c_r * t_r
        let challenges = self.check_az_amortization_correctness(proof, &mut transcript)?;

        // 5. lne 16: check <z, z> ?= \sum(g_ij * c_i * c_j)
        self.check_g_correctness(proof, &challenges)?;

        // 6. line 17: check \sum(<\phi_i, z>c_i) ?= \sum(h_ij * c_i * c_j)
        self.constant_aggregator.calculate_agg_phi_double_prime(
            &self.st.phi_ct,
            &projections.get_conjugated_projection_matrices(),
            &psi,
            &omega,
        );
        self.funcs_aggregator.calculate_aggr_phi(
            &self.st.phi_constraint,
            self.constant_aggregator.get_phi_double_prime(),
            &vector_alpha,
            &vector_beta,
        );
        self.check_h_correctness(proof, &challenges)?;

        // 7. line 18: check \sum(a_ij * g_ij) + \sum(h_ii) - b ?= 0
        self.constant_aggregator
            .calculate_agg_a_double_prime(&psi, &self.st.a_ct);
        self.funcs_aggregator.calculate_agg_a(
            &self.st.a_constraint,
            self.constant_aggregator.get_alpha_double_prime(),
            &vector_alpha,
            &vector_beta,
        );

        self.funcs_aggregator.calculate_aggr_b(
            &self.st.b_constraint,
            &proof.b_ct_aggr,
            &vector_alpha,
            &vector_beta,
        );
        self.check_aggregated_relation(proof)?;

        // 8. line 19: u_1 ?= \sum(\sum(B_ik * t_i^(k))) + \sum(\sum(C_ijk * g_ij^(k)))
        self.check_u1(proof)?;

        // 9. line 20: u_2 ?= \sum(\sum(D_ijk * h_ij^(k)))
        self.check_u2(proof)?;

        Ok(true)
    }

    fn check_vector_p_norm_bound<S: Sponge>(
        &self,
        proof: &LabradorTranscript<S>,
        transcript: &mut LabradorTranscript<S>,
    ) -> Result<Projection, VerifierError> {
        transcript.absorb_u1(&proof.u1);

        let projections = transcript.generate_projections(
            env_params::SECURITY_PARAMETER,
            self.params.rank,
            self.params.multiplicity,
        );
        transcript.absorb_vector_p(&proof.vector_p);
        if proof.vector_p.l2_norm_squared().to_u128()
            > 128 * self.params.beta.to_u128() * self.params.beta.to_u128()
        {
            return Err(VerifierError::NormSumExceeded {
                norm: proof.vector_p.l2_norm_squared().to_u128(),
                allowed: 128 * self.params.beta.to_u128() * self.params.beta.to_u128(),
                step: String::from("vector p norm check"),
            });
        }
        Ok(projections)
    }

    fn check_b_double_prime_constant<S: Sponge>(
        &self,
        proof: &LabradorTranscript<S>,
        transcript: &mut LabradorTranscript<S>,
    ) -> Result<(Vec<Vec<Zq>>, Vec<Vec<Zq>>), VerifierError> {
        let size_of_psi = usize::div_ceil(env_params::SECURITY_PARAMETER, self.params.log_q);
        let size_of_omega = size_of_psi;
        let psi = transcript.generate_vector_psi(size_of_psi, self.params.constraint_l);
        let omega = transcript.generate_vector_omega(size_of_omega, env_params::SECURITY_PARAMETER);
        transcript.absorb_vector_b_ct_aggr(&proof.b_ct_aggr);

        for k in 0..self.params.kappa {
            let b_0_poly = proof.b_ct_aggr.get_elements()[k].get_coefficients()[0];
            let mut b_0: Zq = (0..self.params.constraint_l)
                .map(|l| psi[k][l] * self.st.b_0_ct[l])
                .sum();

            let inner_omega_p =
                inner_product::compute_linear_combination(&omega[k], &proof.vector_p);
            b_0 += inner_omega_p;
            if b_0 != b_0_poly {
                return Err(VerifierError::B0Mismatch {
                    index: k,
                    expected: b_0_poly,
                    computed: b_0,
                });
            }
        }

        Ok((psi, omega))
    }

    fn check_final_norm_sum<S: Sponge>(
        &self,
        proof: &LabradorTranscript<S>,
    ) -> Result<(), VerifierError> {
        // decompose z into z = z^(0) + z^(1) * b, only two parts.
        let z_ij = RqVector::decompose(&proof.z, self.params.b, 2);
        let t_ij: Vec<Vec<RqVector>> = proof
            .t
            .get_elements()
            .iter()
            .map(|i| RqVector::decompose(i, self.params.b, self.params.t_1))
            .collect();
        let g_ij: Vec<Vec<RqVector>> = proof
            .g
            .get_elements()
            .iter()
            .map(|i| RqVector::decompose(i, self.params.b, self.params.t_2))
            .collect();
        let h_ij: Vec<Vec<RqVector>> = proof
            .h
            .get_elements()
            .iter()
            .map(|i| RqVector::decompose(i, self.params.b, self.params.t_1))
            .collect();
        let norm_z_ij = z_ij
            .iter()
            .fold(Zq::ZERO, |acc, p| acc + p.l2_norm_squared());
        let norm_t_ij = Self::norm_squared(&t_ij);
        let norm_g_ij = Self::norm_squared(&g_ij);
        let norm_h_ij = Self::norm_squared(&h_ij);
        let norm_sum = norm_z_ij + norm_t_ij + norm_g_ij + norm_h_ij;

        if norm_sum > self.params.beta * self.params.beta {
            return Err(VerifierError::NormSumExceeded {
                norm: norm_sum.to_u128(),
                allowed: (self.params.beta * self.params.beta).to_u128(),
                step: String::from("Step 14 in verification"),
            });
        }
        Ok(())
    }

    fn check_az_amortization_correctness<S: Sponge>(
        &self,
        proof: &LabradorTranscript<S>,
        transcript: &mut LabradorTranscript<S>,
    ) -> Result<RqVector, VerifierError> {
        let challenges =
            transcript.generate_challenges(env_params::OPERATOR_NORM, self.params.multiplicity);
        let az = self.crs.commitment_scheme_a.matrix() * &proof.z;
        let ct_sum = inner_product::compute_linear_combination(
            proof.t.get_elements(),
            challenges.get_elements(),
        );
        if az != ct_sum {
            return Err(VerifierError::AzError {
                computed: az,
                expected: ct_sum,
            });
        }
        Ok(challenges)
    }

    fn check_g_correctness<S: Sponge>(
        &self,
        proof: &LabradorTranscript<S>,
        challenges: &RqVector,
    ) -> Result<(), VerifierError> {
        let z_inner = inner_product::compute_linear_combination(
            proof.z.get_elements(),
            proof.z.get_elements(),
        );
        let sum_gij_cij = Self::calculate_gh_ci_cj(&proof.g, challenges, self.params.multiplicity);

        if z_inner != sum_gij_cij {
            return Err(VerifierError::ZInnerError {
                computed: z_inner,
                expected: sum_gij_cij,
            });
        }
        Ok(())
    }

    fn check_h_correctness<S: Sponge>(
        &self,
        proof: &LabradorTranscript<S>,
        challenges: &RqVector,
    ) -> Result<(), VerifierError> {
        let sum_phi_z_c =
            Self::calculate_phi_z_c(self.funcs_aggregator.get_appr_phi(), challenges, &proof.z);
        let sum_hij_cij = Self::calculate_gh_ci_cj(&proof.h, challenges, self.params.multiplicity);

        // Left side multiple by 2 because of when we calculate h_ij, we didn't apply the division (divided by 2)
        if &sum_phi_z_c * &Zq::TWO != sum_hij_cij {
            return Err(VerifierError::PhiError {
                computed: &sum_phi_z_c * &Zq::TWO,
                expected: sum_hij_cij,
            });
        }
        Ok(())
    }

    fn check_aggregated_relation<S: Sponge>(
        &self,
        proof: &LabradorTranscript<S>,
    ) -> Result<(), VerifierError> {
        if !Self::check_relation(
            self.funcs_aggregator.get_agg_a(),
            self.funcs_aggregator.get_aggr_b(),
            &proof.g,
            &proof.h,
        ) {
            return Err(VerifierError::RelationCheckFailed);
        }
        Ok(())
    }

    fn check_u1<S: Sponge>(&self, proof: &LabradorTranscript<S>) -> Result<(), VerifierError> {
        let u_1 = &proof.u1;
        let commitment_u1 = outer_commitments::compute_u1(
            self.crs,
            &proof.t,
            DecompositionParameters::new(self.params.b, self.params.t_1)?,
            &proof.g,
            DecompositionParameters::new(self.params.b, self.params.t_2)?,
        );

        if proof.u1 != commitment_u1 {
            return Err(VerifierError::OuterCommitError {
                computed: u_1.clone(),
                expected: commitment_u1,
            });
        }
        Ok(())
    }

    fn check_u2<S: Sponge>(&self, proof: &LabradorTranscript<S>) -> Result<(), VerifierError> {
        let commitment_u2 = outer_commitments::compute_u2(
            self.crs,
            &proof.h,
            DecompositionParameters::new(self.params.b, self.params.t_1)?,
        );

        if proof.u2 != commitment_u2 {
            return Err(VerifierError::OuterCommitError {
                computed: commitment_u2,
                expected: proof.u2.clone(),
            });
        }
        Ok(())
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
            .zip(c.get_elements())
            .map(|(phi_i, c_i)| {
                &(inner_product::compute_linear_combination(phi_i.get_elements(), z.get_elements()))
                    * c_i
            })
            .fold(Rq::zero(), |acc, x| &acc + &x)
    }

    fn norm_squared(polys: &[Vec<RqVector>]) -> Zq {
        polys.iter().fold(Zq::ZERO, |acc, poly| {
            acc + poly
                .iter()
                .fold(Zq::ZERO, |acc, p| acc + p.l2_norm_squared())
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

        let b_primes2 = b_primes * &Zq::TWO;
        let sum_a_primes_g2 = &sum_a_primes_g * &Zq::TWO;

        &sum_a_primes_g2 + &sum_h_ii == b_primes2
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prover::LabradorProver;
    use crate::relation::witness::Witness;
    use crate::transcript::sponges::shake::ShakeSponge;

    #[test]
    fn test_verify() {
        let ep_1 = EnvironmentParameters::default();
        // generate a random witness based on ep above
        let witness_1 = Witness::new(ep_1.rank, ep_1.multiplicity, ep_1.beta);
        // generate public statements based on witness_1
        let st: Statement = Statement::new(&witness_1, &ep_1);
        // generate the common reference string matrices
        let crs = AjtaiInstances::new(&ep_1);

        // create a new prover
        let mut prover = LabradorProver::new(&ep_1, &crs, &witness_1, &st);
        let proof: LabradorTranscript<ShakeSponge> = prover.prove().unwrap();

        // create a new verifier
        let mut verifier = LabradorVerifier::new(&ep_1, &crs, &st);
        let result = verifier.verify(&proof);
        assert!(result.unwrap());
    }
}
