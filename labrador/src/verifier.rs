//! Labrador Verifier
//!
//! This module implements the **verifier** side of the protocol in
//! Fig. 2 (and the verification algorithm in Fig. 3) of the paper.  The prover
//! proof `proof` is accepted **iff** every arithmetic check listed in Fig. 2
//! and lines 8–20 of Fig. 3 succeeds.
//!
//! The verifier receives:
//! * a public statement `st`;
//! * the CRS containing the Ajtai commitment matrices \(\mathbf{A},\mathbf{B},\mathbf{C},\mathbf{D}\);
//! * a proof `proof` comprising
//!   \[(\mathbf{u}_1,\mathbf{u}_2,\,\vec p,\,b''^{(k)}\,\vec z,\tilde\mathbf t,\tilde\mathbf g,\tilde\mathbf h)\]  
//!   exactly as dispatched by the prover.
//!
//! ## Outline of the nine verifier checks
//! The implementation follows Fig. 3 literally:
//!
//! On failure the function returns [`VerifierError`] precisely identifying the
//! violated condition.

#![allow(clippy::result_large_err)]

use crate::commitments::{
    common_instances::AjtaiInstances,
    outer_commitments::{self, DecompositionParameters},
};
use crate::core::aggregate::{FunctionsAggregation, ZeroConstantFunctionsAggregation};
use crate::core::{inner_product, jl::Projection};
use crate::relation::env_params;
use crate::relation::{env_params::EnvironmentParameters, statement::Statement};
use crate::ring::{rq::Rq, rq_matrix::RqMatrix, rq_vector::RqVector, Norms};
use crate::ring::zq::ZqLabrador;
type Zq = ZqLabrador;
use crate::transcript::{LabradorTranscript, Sponge};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum VerifierError {
    #[error("matrix not symmetric at ({i},{j}): expected {expected:?}, found {found:?}")]
    NotSymmetric {
        i: u64,
        j: u64,
        expected: Rq,
        found: Rq,
    },
    #[error("B0 mismatch at index {index}: expected {expected}, computed {computed}")]
    B0Mismatch {
        index: u64,
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

/// Implements the algorithm executed by the verifier \(\mathcal{V}\) in the paper.
pub struct LabradorVerifier<'a> {
    /// System‑wide environment parameters like rank and multiplicity.
    params: &'a EnvironmentParameters,
    /// Common reference string (CRS) holding the commitment matrices
    /// \(\mathbf{A},\mathbf{B},\mathbf{C},\mathbf{D}\).
    crs: &'a AjtaiInstances,
    /// Public relation statement.
    st: &'a Statement,
    /* --- aggregation instances ---------------------------------------------------------- */
    constant_aggregator: ZeroConstantFunctionsAggregation<'a>,
    funcs_aggregator: FunctionsAggregation<'a>,
}

impl<'a> LabradorVerifier<'a> {
    /// Constructs a new verifier instance.
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

    /// Validate JL projection bound  (Fig. 2):
    /// \[
    ///   \|\vec p\|_2^2 = \sum_{j=1}^{2\lambda} p_j^2 \;\stackrel{?}{\le}\; 128\,\beta^2.
    /// \]
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
        if proof.vector_p.l2_norm_squared() > 128 * self.params.beta_sq {
            return Err(VerifierError::NormSumExceeded {
                norm: proof.vector_p.l2_norm_squared(),
                allowed: 128 * self.params.beta_sq,
                step: String::from("vector p norm check"),
            });
        }
        Ok(projections)
    }

    /// Validate b'' constant term (Fig. 2):
    /// For every `k`, verifies
    /// \[
    ///   b_{0}^{''(k)} \stackrel{?}{=} \sum_{l=1}^{L} \psi_{l}^{(k)}\,b_{0}^{'(l)}
    ///                    + \langle\,\vec\omega^{(k)},\vec p\rangle .
    /// \]
    #[allow(clippy::type_complexity)]
    fn check_b_double_prime_constant<S: Sponge>(
        &self,
        proof: &LabradorTranscript<S>,
        transcript: &mut LabradorTranscript<S>,
    ) -> Result<(Vec<Vec<Zq>>, Vec<Vec<Zq>>), VerifierError> {
        let size_of_psi = {
            let sec_param = u64::try_from(env_params::SECURITY_PARAMETER)
                .expect("SECURITY_PARAMETER does not fit in u64");
            let log_q = u64::try_from(self.params.log_q).expect("log_q does not fit in u64");
            usize::try_from(sec_param.div_ceil(log_q))
                .expect("div_ceil result does not fit in usize")
        };
        let size_of_omega = size_of_psi;
        let psi = transcript.generate_vector_psi(size_of_psi, self.params.constraint_l);
        let omega = transcript.generate_vector_omega(size_of_omega, env_params::SECURITY_PARAMETER);
        transcript.absorb_vector_b_ct_aggr(&proof.b_ct_aggr);

        let kappa_u64 = u64::try_from(self.params.kappa).expect("kappa does not fit in u64");
        for k in 0..kappa_u64 {
            let k_usize = usize::try_from(k).expect("k does not fit in usize");
            let b_0_poly = proof.b_ct_aggr.elements()[k_usize].coeffs()[0];
            let constraint_l_u64 =
                u64::try_from(self.params.constraint_l).expect("constraint_l does not fit in u64");
            let mut b_0: Zq = (0..constraint_l_u64)
                .map(|l| {
                    let l_usize = usize::try_from(l).expect("l does not fit in usize");
                    psi[k_usize][l_usize] * self.st.b_0_ct[l_usize]
                })
                .sum();

            let inner_omega_p =
                inner_product::compute_linear_combination(&omega[k_usize], &proof.vector_p);
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

    /// Validate consolidated norm bound  (Fig. 3 line 14):
    /// \[
    ///  \sum_{\ell=0}^{1}\|\mathbf z^{(\ell)}\|_2^2
    ///   + \sum_{i=1}^r\sum_{k=0}^{t_1-1}\|\mathbf t_i^{(k)}\|_2^2
    ///   + \sum_{1\le i\le j\le r}\sum_{k=0}^{t_2-1}\|\mathbf g_{ij}^{(k)}\|_2^2
    ///   + \sum_{1\le i\le j\le r}\sum_{k=0}^{t_1-1}\|\mathbf h_{ij}^{(k)}\|_2^2
    ///     \stackrel{?}{\le} (\beta')^2 .
    ///     \]
    fn check_final_norm_sum<S: Sponge>(
        &self,
        proof: &LabradorTranscript<S>,
    ) -> Result<(), VerifierError> {
        // decompose z into z = z^(0) + z^(1) * b, only two parts.
        let z_ij = proof.z.decompose(self.params.b, 2);
        let t_ij: Vec<Vec<RqVector>> = proof
            .t
            .elements()
            .iter()
            .map(|i| RqVector::decompose(i, self.params.b, self.params.t_1))
            .collect();
        let g_ij: Vec<Vec<RqVector>> = proof
            .g
            .elements()
            .iter()
            .map(|i| RqVector::decompose(i, self.params.b, self.params.t_2))
            .collect();
        let h_ij: Vec<Vec<RqVector>> = proof
            .h
            .elements()
            .iter()
            .map(|i| RqVector::decompose(i, self.params.b, self.params.t_1))
            .collect();
        let norm_z_ij = z_ij.iter().fold(0, |acc, p| acc + p.l2_norm_squared());
        let norm_t_ij = Self::norm_squared(&t_ij);
        let norm_g_ij = Self::norm_squared(&g_ij);
        let norm_h_ij = Self::norm_squared(&h_ij);
        let norm_sum = norm_z_ij + norm_t_ij + norm_g_ij + norm_h_ij;

        if norm_sum > self.params.beta_prime_sq {
            return Err(VerifierError::NormSumExceeded {
                norm: norm_sum,
                allowed: self.params.beta_prime_sq,
                step: String::from("Step 14 in verification"),
            });
        }
        Ok(())
    }

    /// Validate amortised Ajtai relation (Fig. 3 line 15):
    /// \[
    ///   \mathbf A\,\mathbf z \stackrel{?}{=} \sum_{i=1}^r c_i\,\tilde\mathbf t_i
    /// \]
    /// where the challenges `c_i` are re‑sampled locally from the sponge.
    fn check_az_amortization_correctness<S: Sponge>(
        &self,
        proof: &LabradorTranscript<S>,
        transcript: &mut LabradorTranscript<S>,
    ) -> Result<RqVector, VerifierError> {
        let challenges =
            transcript.generate_challenges(env_params::OPERATOR_NORM, self.params.multiplicity);
        let az = self.crs.commitment_scheme_a.matrix() * &proof.z;
        let ct_sum =
            inner_product::compute_linear_combination(proof.t.elements(), challenges.elements());
        if az != ct_sum {
            return Err(VerifierError::AzError {
                computed: az,
                expected: ct_sum,
            });
        }
        Ok(challenges)
    }

    /// Validate the following (Fig. 3 line 16):
    /// \[
    ///   \langle\mathbf z,\mathbf z\rangle \stackrel{?}{=}
    ///   \sum_{i,j} g_{ij}\,c_i\,c_j .
    /// \]
    fn check_g_correctness<S: Sponge>(
        &self,
        proof: &LabradorTranscript<S>,
        challenges: &RqVector,
    ) -> Result<(), VerifierError> {
        let z_inner =
            inner_product::compute_linear_combination(proof.z.elements(), proof.z.elements());
        let sum_gij_cij = Self::calculate_gh_ci_cj(&proof.g, challenges, self.params.multiplicity);

        if z_inner != sum_gij_cij {
            return Err(VerifierError::ZInnerError {
                computed: z_inner,
                expected: sum_gij_cij,
            });
        }
        Ok(())
    }

    /// Validate the following (Fig. 3 line 17):
    /// \[
    ///   2\,\sum_{i=1}^r \langle\varphi_i,\mathbf z\rangle c_i
    ///   \stackrel{?}{=} \sum_{i,j} h_{ij}\,c_i\,c_j .
    /// \]
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

    /// Validate the following (Fig. 3 line 18):
    /// \[
    ///   2\sum_{i,j} a_{ij} g_{ij} + \sum_i h_{ii}\;\stackrel{?}{=}\; 2b .
    /// \]
    fn check_aggregated_relation<S: Sponge>(
        &self,
        proof: &LabradorTranscript<S>,
    ) -> Result<(), VerifierError> {
        let r = self.funcs_aggregator.get_agg_a().elements().len();

        let mut sum_a_primes_g = Rq::zero();
        // walk only over the stored half: i ≤ j
        for i in 0..r {
            for j in 0..r {
                sum_a_primes_g = &sum_a_primes_g
                    + &(self.funcs_aggregator.get_agg_a().get_cell(i, j) * proof.g.get_cell(i, j));
            }
        }

        let sum_h_ii = (0..r).fold(Rq::zero(), |acc, i| &acc + proof.h.get_cell(i, i));

        let b_primes2 = self.funcs_aggregator.get_aggr_b() * &Zq::TWO;
        let sum_a_primes_g2 = &sum_a_primes_g * &Zq::TWO;

        if &sum_a_primes_g2 + &sum_h_ii != b_primes2 {
            return Err(VerifierError::RelationCheckFailed);
        }
        Ok(())
    }

    /// Validate uter commitment u1 (Fig. 3 line 19):
    /// \[
    ///   \mathbf u_1^{\text{re}} = \sum_{i=1}^r\sum_{k=0}^{t_1-1}\mathbf B_{ik}\,\mathbf t_i^{(k)}
    ///          + \sum_{1\le i\le j\le r}\sum_{k=0}^{t_2-1}\mathbf C_{ijk}\,\mathbf g_{ij}^{(k)}
    /// \]
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

    /// Validate outer commitment u2 (Fig. 3 line 20):
    /// \[
    ///   \mathbf u_2^{\text{re}} = \sum_{1\le i\le j\le r}\sum_{k=0}^{t_1-1}
    ///                              \mathbf D_{ijk}\,\mathbf h_{ij}^{(k)}
    /// \]
    /// and checks equality with the transcript.
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
        let c_elements = random_c.elements();
        let mut result = Rq::zero();
        for i in 0..r {
            let c_i = &c_elements[i];
            for (j, c_j) in c_elements.iter().enumerate().take(r) {
                // x_ij[i,j] * c_i * c_j
                let term = &(&(x_ij.get_cell(i, j) * c_i) * c_j);
                result = &result + term;
            }
        }
        result
    }

    /// calculate the left hand side of line 17, \sum(<\phi_z, z> * c_i)
    fn calculate_phi_z_c(phi: &[RqVector], c: &RqVector, z: &RqVector) -> Rq {
        phi.iter()
            .zip(c.elements())
            .map(|(phi_i, c_i)| {
                let inner_prod =
                    inner_product::compute_linear_combination(phi_i.elements(), z.elements());
                &inner_prod * c_i
            })
            .fold(Rq::zero(), |acc, term| &acc + &term)
    }

    fn norm_squared(polys: &[Vec<RqVector>]) -> u128 {
        polys.iter().fold(0, |acc, poly| {
            acc + poly.iter().fold(0, |acc, p| acc + p.l2_norm_squared())
        })
    }

    /// Verifies a prover transcript.
    /// Executes **all** verifier steps (Fig. 2 page 17, Fig. 3 page 18).
    ///
    /// Returns `Ok(true)` if *all* checks succeed, otherwise an explanatory
    /// [`VerifierError`].  The sponge seeded inside this routine **must** be
    /// identical to the prover’s so that the generated randomizers match.
    pub fn verify<S: Sponge>(
        &mut self,
        proof: &LabradorTranscript<S>,
    ) -> Result<bool, VerifierError> {
        // fresh sponge — deterministic (Todo: Add domain separator)
        let mut transcript = LabradorTranscript::new(S::default());

        // --- Fig. 2: JL Projections Norm Check ------------------------------------------------------------------
        let projections = self.check_vector_p_norm_bound(proof, &mut transcript)?;

        // --- Fig. 2: Constant of Constant Functions Aggregation Check -------------------------------------------
        let (psi, omega) = self.check_b_double_prime_constant(proof, &mut transcript)?;

        // --- Fig. 3, Line 14: Check Prover's Last Message Norm --------------------------------------------------
        self.check_final_norm_sum(proof)?;

        // Generate Randomness for Next Steps
        let vector_alpha = transcript.generate_rq_vector(self.params.constraint_k);
        let vector_beta = transcript.generate_rq_vector(usize::div_ceil(
            env_params::SECURITY_PARAMETER,
            self.params.log_q,
        ));
        transcript.absorb_u2(&proof.u2);

        // --- Fig. 3, Line 15: Check Correctness of Amortization --------------------------------------------------
        let challenges = self.check_az_amortization_correctness(proof, &mut transcript)?;

        // --- Fig. 3, Line 16: Check Relation of z, g, and challenges ---------------------------------------------
        self.check_g_correctness(proof, &challenges)?;

        // Generate Randomness for Next Steps
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
        // --- Fig. 3, Line 17: Check Relation of z, \varphi, h, and challenges -------------------------------------
        self.check_h_correctness(proof, &challenges)?;

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
        // --- Fig. 3, Line 18: Check Aggregated Elements Equal to 0 -----------------------------------------------
        self.check_aggregated_relation(proof)?;

        // --- Fig. 3, Line 19-20: Check Correntness of u1 and u2---------------------------------------------------
        self.check_u1(proof)?;
        self.check_u2(proof)?;

        Ok(true)
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
        let witness_1 = Witness::new(ep_1.rank, ep_1.multiplicity, ep_1.beta_sq);
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
