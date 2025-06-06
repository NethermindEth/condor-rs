use crate::commitments::ajtai_commitment;
use crate::commitments::common_instances::AjtaiInstances;
use crate::commitments::outer_commitments;
use crate::commitments::outer_commitments::DecompositionParameters;
use crate::commitments::CommitError;
use crate::core::aggregate::FunctionsAggregation;
use crate::core::aggregate::ZeroConstantFunctionsAggregation;
use crate::core::garbage_polynomials;
use crate::core::inner_product;
use crate::core::jl::Projection;
use crate::relation::env_params;
use crate::relation::witness::Witness;
use crate::relation::{env_params::EnvironmentParameters, statement::Statement};
use crate::ring::rq_matrix::RqMatrix;
use crate::ring::rq_vector::RqVector;
use crate::ring::zq::Zq;
use crate::transcript::LabradorTranscript;
use crate::transcript::Sponge;
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

pub struct LabradorProver<'a> {
    params: &'a EnvironmentParameters,
    crs: &'a AjtaiInstances,
    witness: &'a Witness,
    st: &'a Statement,
    // Aggregation instances
    constant_aggregator: ZeroConstantFunctionsAggregation<'a>,
    funcs_aggregator: FunctionsAggregation<'a>,
}

impl<'a> LabradorProver<'a> {
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

    fn compute_vector_ti(&self) -> Result<RqMatrix, CommitError> {
        // Ajtai Commitments t_i = A * s_i
        let commitments = self
            .witness
            .s
            .iter()
            .cloned()
            .map(|s_i| self.crs.commitment_scheme_a.commit(&s_i))
            .collect::<Result<Vec<_>, CommitError>>()?;

        Ok(RqMatrix::new(commitments, false))
    }

    fn compute_u1<S: Sponge>(
        &mut self,
        transcript: &mut LabradorTranscript<S>,
    ) -> Result<(RqMatrix, RqMatrix), ProverError> {
        let t_i = self.compute_vector_ti()?;
        // g_ij = <s_i, s_j>
        let garbage_polynomial_g = garbage_polynomials::compute_g(&self.witness.s);
        // calculate outer commitment u_1 = \sum(B_ik * t_i^(k)) + \sum(C_ijk * g_ij^(k))
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

    fn compute_b_double_prime<S: Sponge>(
        &mut self,
        transcript: &mut LabradorTranscript<S>,
        projections: &Projection,
    ) {
        let vector_psi =
            transcript.generate_vector_psi(self.params.const_agg_length, self.params.constraint_l);
        let vector_omega = transcript
            .generate_vector_omega(self.params.const_agg_length, env_params::SECURITY_PARAMETER);
        // first aggregation
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

    fn compute_u2<S: Sponge>(
        &mut self,
        transcript: &mut LabradorTranscript<S>,
    ) -> Result<RqMatrix, ProverError> {
        let alpha_vector = transcript.generate_rq_vector(self.params.constraint_k);
        let beta_vector = transcript.generate_rq_vector(self.params.const_agg_length);
        self.funcs_aggregator.calculate_aggr_phi(
            &self.st.phi_constraint,
            self.constant_aggregator.get_phi_double_prime(),
            &alpha_vector,
            &beta_vector,
        );

        // Step 4: Calculate h_ij, u_2, and z starts: ---------------------------------------
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

    // calculate z = c_1*s_1 + ... + c_r*s_r
    fn compute_z<S: Sponge>(&mut self, transcript: &mut LabradorTranscript<S>) -> RqVector {
        let challenges =
            transcript.generate_challenges(env_params::OPERATOR_NORM, self.params.multiplicity);
        let z =
            inner_product::compute_linear_combination(&self.witness.s, challenges.get_elements());
        z
    }

    /// all prove steps are from page 17
    pub fn prove<S: Sponge>(&mut self) -> Result<LabradorTranscript<S>, ProverError> {
        // Generate random challenges used between prover and verifier
        let mut transcript = LabradorTranscript::new(S::default());

        // Step 1: Outer commitments u_1 starts: --------------------------------------------
        let (t_i, garbage_polynomial_g) = self.compute_u1(&mut transcript)?;
        // Step 1: Outer commitments u_1 ends: ----------------------------------------------

        // Step 2: JL projection starts: ----------------------------------------------------
        let projections = self.compute_p(&mut transcript);
        // Step 2: JL projection ends: ------------------------------------------------------

        // Step 3: Aggregation starts: --------------------------------------------------------------
        self.compute_b_double_prime(&mut transcript, &projections);

        // second aggregation

        // Aggregation ends: ----------------------------------------------------------------
        let garbage_polynomial_h = self.compute_u2(&mut transcript)?;

        let z = self.compute_z(&mut transcript);

        transcript.set_recursive_part(z, t_i, garbage_polynomial_g, garbage_polynomial_h);

        // Step 4: Calculate h_ij, u_2, and z ends: -----------------------------------------

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
        let witness_1 = Witness::new(ep_1.rank, ep_1.multiplicity, ep_1.beta);
        // generate public statement based on witness_1
        let st: Statement = Statement::new(&witness_1, &ep_1);
        // generate the common reference string matrices A, B, C, D
        let crs: AjtaiInstances = AjtaiInstances::new(&ep_1);

        // create a new prover
        let mut prover = LabradorProver::new(&ep_1, &crs, &witness_1, &st);
        let _: LabradorTranscript<ShakeSponge> = prover.prove().unwrap();
    }
}
