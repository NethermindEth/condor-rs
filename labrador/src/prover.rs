//! Labrador Prover
//!
//! This module implements the **prover** side of the protocol sketched in
//! Fig. 2 of the paper (page 17). All algebra is carried out over the ring
//! \(R_q = \mathbb{Z}_q[x]/(x^d+1)\), where `q` is 32-bit modulo and `d` is polynomial degree.  
//!
//! **Important — security parameters**:  The protocol is parameterised by
//! *rank* `n`, *multiplicity* `r`, *operator norm* `β`, decomposition depths `t_1,t_2`, and
//! security parameter `λ` = [`env_params::SECURITY_PARAMETER`].

use crate::commitments::{
    ajtai_commitment, common_instances::AjtaiInstances, outer_commitments,
    outer_commitments::DecompositionParameters, CommitError,
};
use crate::core::{
    aggregate::FunctionsAggregation, aggregate::ZeroConstantFunctionsAggregation,
    garbage_polynomials, inner_product, jl::Projection,
};
use crate::relation::{
    env_params::{self, EnvironmentParameters},
    statement::Statement,
    witness::Witness,
};
use crate::ring::{rq_matrix::RqMatrix, rq_vector::RqVector};
use crate::ring::zq::ZqLabrador;
type Zq = ZqLabrador;
use crate::transcript::{LabradorTranscript, Sponge};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum ProverError {
    /// Indicates that the L2 norm (squared) of the witness exceeded the allowed threshold.
    #[error("invalid witness size: norm_squared {norm_squared}, allowed {allowed}")]
    WitnessL2NormViolated { norm_squared: Zq, allowed: Zq },
    #[error("Invalid Projection of index {index}. Expected {expected}, got {computed}")]
    ProjectionError {
        index: usize,
        expected: Zq,
        computed: Zq,
    },
    #[error("commitment failure")]
    CommitError(#[from] ajtai_commitment::CommitError),
    #[error("decomposition failure")]
    DecompositionError(#[from] outer_commitments::DecompositionError),
}

/// Implements the steps executed by the prover \(\mathcal{P}\) in the paper.
pub struct LabradorProver<'a> {
    /// System‑wide environment parameters like rank and multiplicity.
    params: &'a EnvironmentParameters,
    /// Common reference string (CRS) holding the commitment matrices
    /// \(\mathbf{A},\mathbf{B},\mathbf{C},\mathbf{D}\).
    crs: &'a AjtaiInstances,
    /// Secret witness vectors \(\{\mathbf{s}_1,\dots,\mathbf{s}_r\}\).
    witness: &'a Witness,
    /// Public relation statement.
    st: &'a Statement,
    /* --- aggregation instances ---------------------------------------------------------- */
    constant_aggregator: ZeroConstantFunctionsAggregation<'a>,
    funcs_aggregator: FunctionsAggregation<'a>,
}

impl<'a> LabradorProver<'a> {
    /// Constructs a new prover instance.
    pub fn new(
        params: &'a EnvironmentParameters,
        crs: &'a AjtaiInstances,
        witness: &'a Witness,
        st: &'a Statement,
    ) -> Self {
        Self {
            params,
            crs,
            witness,
            st,
            constant_aggregator: ZeroConstantFunctionsAggregation::new(params),
            funcs_aggregator: FunctionsAggregation::new(params),
        }
    }

    // -----------------------------------------------------------------------------
    // Step 1.a — Ajtai commitments t_i
    // -----------------------------------------------------------------------------
    /// Computes the **Ajtai commitments**
    /// \[\mathbf{t}_i = \mathbf{A}\,\mathbf{s}_i\\].
    ///
    /// # Errors
    /// Returns [`CommitError`] if the commitment scheme fails.
    fn compute_vector_ti(&self) -> Result<RqMatrix, CommitError> {
        let commitments = self
            .witness
            .s
            .iter()
            .cloned()
            .map(|s_i| self.crs.commitment_scheme_a.commit(&s_i))
            .collect::<Result<Vec<_>, CommitError>>()?;

        Ok(RqMatrix::new(commitments, false))
    }

    // -----------------------------------------------------------------------------
    // Step 1.b — First outer commitment u_1
    // -----------------------------------------------------------------------------

    /// Computes the *first outer commitment* \(\mathbf{u}_1\).
    ///
    /// Internally this function
    /// 1. calls [`compute_vector_ti`] to obtain all \(\mathbf{t}_i\),
    /// 2. evaluates the garbage polynomials \(\mathbf{g}_{ij}=\langle\mathbf{s}_i,\mathbf{s}_j\rangle\),
    /// 4. finally aggregates and commit using [`compute_u1`].
    ///
    /// The resulting value is
    /// \[
    ///     \mathbf{u}_1
    ///      = \sum_{i=1}^r \sum_{k=0}^{t_1-1} \mathbf{B}_{ik}\,\mathbf{t}_i^{(k)}
    ///      + \sum_{i\le j} \sum_{k=0}^{t_2-1}
    ///        \mathbf{C}_{ijk}\,\mathbf{g}_{ij}^{(k)}.
    /// \]
    ///
    /// The commitment is recorded inside the transcript so the verifier can later
    /// retrieve it.
    fn compute_u1<S: Sponge>(
        &mut self,
        transcript: &mut LabradorTranscript<S>,
    ) -> Result<(RqMatrix, RqMatrix), ProverError> {
        let t_i = self.compute_vector_ti()?;
        let garbage_polynomial_g = garbage_polynomials::compute_g(&self.witness.s);
        let commitment_u1 = outer_commitments::compute_u1(
            self.crs,
            &t_i,
            DecompositionParameters::new(self.params.b, self.params.t_1)?,
            &garbage_polynomial_g,
            DecompositionParameters::new(self.params.b, self.params.t_2)?,
        );
        transcript.set_u1(commitment_u1);
        Ok((t_i, garbage_polynomial_g))
    }

    // -----------------------------------------------------------------------------
    // Step 2 — JL projections
    // -----------------------------------------------------------------------------

    /// Computes the Johnson–Lindenstrauss projections `p`.
    ///
    /// \[p_j = \sum_{i=1}^r \langle {\pi}_i^{(j)},\mathbf{s}_i\rangle\]
    ///
    /// |p| =  `2 * SECURITY_PARAMETER`
    ///
    /// p is recorded inside the transcript so the verifier can later
    /// retrieve it. Projection instance is returned for use in next steps of the prover.
    fn compute_p<S: Sponge>(&self, transcript: &mut LabradorTranscript<S>) -> Projection {
        let projections = transcript.generate_projections(
            env_params::SECURITY_PARAMETER,
            self.params.rank,
            self.params.multiplicity,
        );
        let vector_p = projections.compute_batch_projection(&self.witness.s);
        transcript.set_vector_p(vector_p);
        projections
    }

    // -----------------------------------------------------------------------------
    // Step 3 — Aggregation of zero‑constant functions
    // -----------------------------------------------------------------------------

    /// Computes the aggregated triples \((a''^{(k)},\varphi''_i,b''^{(k)})\).
    ///
    /// We sample vectors *ψ* and *ω* from the sponge and call the helper struct
    /// [`ZeroConstantFunctionsAggregation`] for the computation.
    fn compute_b_double_prime<S: Sponge>(
        &mut self,
        transcript: &mut LabradorTranscript<S>,
        projections: &Projection,
    ) {
        // --- sample randomisers ------------------------------------------------------
        let vector_psi =
            transcript.generate_vector_psi(self.params.const_agg_length, self.params.constraint_l);
        let vector_omega = transcript
            .generate_vector_omega(self.params.const_agg_length, env_params::SECURITY_PARAMETER);
        // --- perform aggregations -------------------------------------------
        self.constant_aggregator
            .calculate_agg_a_double_prime(&vector_psi, &self.st.a_ct);
        self.constant_aggregator.calculate_agg_phi_double_prime(
            &self.st.phi_ct,
            &projections.get_conjugated_projection_matrices(),
            &vector_psi,
            &vector_omega,
        );
        let b_ct_aggr = self
            .constant_aggregator
            .calculate_agg_b_double_prime(&self.witness.s);
        transcript.set_vector_b_ct_aggr(b_ct_aggr);
    }

    // -----------------------------------------------------------------------------
    // Step 4 — Second outer commitment u_2
    // -----------------------------------------------------------------------------

    /// Computes the *second outer commitment* \(\mathbf{u}_2\).
    ///
    /// First we create an *aggregated* vector of functions (via
    /// [`FunctionsAggregation`]). We then evaluate the garbage polynomials
    /// \(h_{ij} = \langle\tilde{\varphi}_i,\tilde{\mathbf{s}}_j\rangle
    ///                    + \langle\varphi_i^{\prime(j)},\tilde{\mathbf{s}}_i\rangle)
    /// and finally commit using the CRS matrices \(\mathbf{D}_{ijk}\):
    ///
    /// \[
    ///     \mathbf{u}_2 = \sum_{i\le j}\sum_{k=0}^{t_1-1}
    ///                     \mathbf{D}_{ijk}\,h_{ij}^{(k)}.
    /// \]
    /// * Note: As dividing by two is not supported in our scheme (Zq is not necessarily a filed), h_ij in the implementation does not have 1/2 factor in the paper. In subsequent verifications, this is considered.
    fn compute_u2<S: Sponge>(
        &mut self,
        transcript: &mut LabradorTranscript<S>,
    ) -> Result<RqMatrix, ProverError> {
        // --- sample randomisers ------------------------------------------------------
        let alpha_vector = transcript.generate_rq_vector(self.params.constraint_k);
        let beta_vector = transcript.generate_rq_vector(self.params.const_agg_length);
        // --- perform aggregations -------------------------------------------
        self.funcs_aggregator.calculate_aggr_phi(
            &self.st.phi_constraint,
            self.constant_aggregator.get_phi_double_prime(),
            &alpha_vector,
            &beta_vector,
        );

        let garbage_polynomial_h =
            garbage_polynomials::compute_h(&self.witness.s, self.funcs_aggregator.get_appr_phi());
        let commitment_u2 = outer_commitments::compute_u2(
            self.crs,
            &garbage_polynomial_h,
            DecompositionParameters::new(self.params.b, self.params.t_1)?,
        );
        transcript.set_u2(commitment_u2);
        Ok(garbage_polynomial_h)
    }

    // -----------------------------------------------------------------------------
    // Step 5 — Linear combination z
    // -----------------------------------------------------------------------------

    /// Computes the random linear combination
    /// \[\mathbf{z}=\sum_{i=1}^r c_i\,\mathbf{s}_i\] with coefficients
    /// `c_i ← C` (uniformly random challenges).
    fn compute_z<S: Sponge>(&mut self, transcript: &mut LabradorTranscript<S>) -> RqVector {
        let challenges =
            transcript.generate_challenges(env_params::OPERATOR_NORM, self.params.multiplicity);
        let z = inner_product::compute_linear_combination(&self.witness.s, challenges.elements());
        z
    }

    /// Executes **all** prover steps (Fig. 2, page 17).
    ///
    /// On success the returned [`LabradorTranscript`] contains every message that
    /// needs to be sent to the verifier.
    pub fn prove<S: Sponge>(&mut self) -> Result<LabradorTranscript<S>, ProverError> {
        // Initialise transcript and sponge -----------------------------------------------
        let mut transcript = LabradorTranscript::new(S::default());

        // --- Step 1: Outer Commitment u1 ------------------------------------------------
        let (t_i, garbage_polynomial_g) = self.compute_u1(&mut transcript)?;

        // --- Step 2: JL Projections -----------------------------------------------------
        let projections = self.compute_p(&mut transcript);

        // --- Step 3: Aggregations +  ----------------------------------------------------
        self.compute_b_double_prime(&mut transcript, &projections);

        // --- Step 4: Outer Commitment u2 ------------------------------------------------
        let garbage_polynomial_h = self.compute_u2(&mut transcript)?;

        // --- Step 5: Amortization -------------------------------------------------------
        let z = self.compute_z(&mut transcript);

        // --- Finalise Transcript: Add t_i, g_ij, h_ij to transcript----------------------
        transcript.set_recursive_part(z, t_i, garbage_polynomial_g, garbage_polynomial_h);

        Ok(transcript)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::transcript::sponges::shake::ShakeSponge;

    #[test]
    fn test_prove() {
        // set up example environment parameters, use default set for testing.
        let ep_1 = EnvironmentParameters::default();
        // generate a random witness based on environment parameters above
        let witness_1 = Witness::new(ep_1.rank, ep_1.multiplicity, ep_1.beta_sq);
        // generate public statement based on witness_1
        let st: Statement = Statement::new(&witness_1, &ep_1);
        // generate the common reference string matrices A, B, C, D
        let crs: AjtaiInstances = AjtaiInstances::new(&ep_1);

        // create a new prover
        let mut prover = LabradorProver::new(&ep_1, &crs, &witness_1, &st);
        let _: LabradorTranscript<ShakeSponge> = prover.prove().unwrap();
    }
}
