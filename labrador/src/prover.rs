use crate::commitments::ajtai_commitment;
use crate::commitments::common_instances::AjtaiInstances;
use crate::commitments::outer_commitments;
use crate::commitments::outer_commitments::DecompositionParameters;
use crate::commitments::outer_commitments::OuterCommitment;
use crate::commitments::CommitError;
use crate::core::aggregate::FunctionsAggregation;
use crate::core::aggregate::ZeroConstantFunctionsAggregation;
use crate::core::garbage_polynomials::GarbagePolynomials;
use crate::core::inner_product;
use crate::core::{env_params::EnvironmentParameters, statement::Statement};
use crate::relation::witness::Witness;
use crate::ring::rq_matrix::RqMatrix;
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
    pub pp: &'a AjtaiInstances,
    pub witness: &'a Witness,
    pub st: &'a Statement,
}

impl<'a> LabradorProver<'a> {
    pub fn new(pp: &'a AjtaiInstances, witness: &'a Witness, st: &'a Statement) -> Self {
        Self { pp, witness, st }
    }

    fn compute_vector_ti(&self) -> Result<RqMatrix, CommitError> {
        // collect all commitments, propagating any CommitError with `?`
        let commitments = self
            .witness
            .s
            .iter()
            .cloned()
            .map(|s_i| self.pp.commitment_scheme_a.commit(&s_i))
            .collect::<Result<Vec<_>, CommitError>>()?;

        Ok(RqMatrix::new(commitments, false))
    }

    /// all prove steps are from page 17
    pub fn prove<S: Sponge>(
        &mut self,
        ep: &EnvironmentParameters,
    ) -> Result<LabradorTranscript<S>, ProverError> {
        // Generate random challenges used between prover and verifier
        let mut transcript = LabradorTranscript::new(S::default());
        // Compute garbage polynomials g and h
        let mut garbage_polynomials = GarbagePolynomials::new(&self.witness.s);
        // Compute outer commitments u1 and u2
        let mut outer_commitments = OuterCommitment::new(self.pp);
        // Aggregation instances
        let mut constant_aggregator = ZeroConstantFunctionsAggregation::new(ep);
        let mut funcs_aggregator = FunctionsAggregation::new(ep);

        // Step 1: Outer commitments u_1 starts: --------------------------------------------

        // Ajtai Commitments t_i = A * s_i
        let t_i = self.compute_vector_ti()?;
        // g_ij = <s_i, s_j>
        garbage_polynomials.compute_g();
        // calculate outer commitment u_1 = \sum(B_ik * t_i^(k)) + \sum(C_ijk * g_ij^(k))
        let commitment_u1 = outer_commitments.compute_u1(
            &t_i,
            DecompositionParameters::new(ep.b, ep.t_1)?,
            &garbage_polynomials.g,
            DecompositionParameters::new(ep.b, ep.t_2)?,
        );
        transcript.set_u1(commitment_u1);
        // Step 1: Outer commitments u_1 ends: ----------------------------------------------

        // Step 2: JL projection starts: ----------------------------------------------------
        let projections =
            transcript.generate_projections(ep.security_parameter, ep.rank, ep.multiplicity);
        let vector_p = projections.compute_batch_projection(&self.witness.s);
        transcript.set_vector_p(vector_p);
        // Step 2: JL projection ends: ------------------------------------------------------

        // Step 3: Aggregation starts: --------------------------------------------------------------
        let vector_psi = transcript.generate_vector_psi(ep.const_agg_length, ep.constraint_l);
        let vector_omega =
            transcript.generate_vector_omega(ep.const_agg_length, ep.security_parameter);
        // first aggregation
        constant_aggregator.calculate_agg_a_double_prime(&vector_psi, &self.st.a_ct);
        constant_aggregator.calculate_agg_phi_double_prime(
            &self.st.phi_ct,
            &projections.get_conjugated_projection_matrices(),
            &vector_psi,
            &vector_omega,
        );
        let b_ct_aggr = constant_aggregator.calculate_agg_b_double_prime(&self.witness.s);
        transcript.set_vector_b_ct_aggr(b_ct_aggr);

        // second aggregation
        let alpha_vector = transcript.generate_rq_vector(ep.constraint_k);
        let beta_vector = transcript.generate_rq_vector(ep.const_agg_length);
        funcs_aggregator.calculate_aggr_phi(
            &self.st.phi_constraint,
            constant_aggregator.get_phi_double_prime(),
            &alpha_vector,
            &beta_vector,
        );
        // Aggregation ends: ----------------------------------------------------------------

        // Step 4: Calculate h_ij, u_2, and z starts: ---------------------------------------
        garbage_polynomials.compute_h(funcs_aggregator.get_appr_phi());
        let commitment_u2 = outer_commitments.compute_u2(
            &garbage_polynomials.h,
            DecompositionParameters::new(ep.b, ep.t_1)?,
        );
        transcript.set_u2(commitment_u2);

        // calculate z = c_1*s_1 + ... + c_r*s_r
        let challenges = transcript.generate_challenges(ep.operator_norm, ep.multiplicity);
        let z =
            inner_product::compute_linear_combination(&self.witness.s, challenges.get_elements());

        transcript.set_recursive_part(z, t_i, garbage_polynomials.g, garbage_polynomials.h);

        // Step 4: Calculate h_ij, u_2, and z ends: -----------------------------------------

        Ok(transcript)
    }

    // The following is not part of the prover.
    // Todo: Add jl projection constraints to the statement.
    // /// check p_j? = ct(sum(<σ−1(pi_i^(j)), s_i>))
    // fn check_projection(&self, p: &[Zq], pi: &[Vec<Vec<Zq>>]) -> Result<bool, ProverError> {
    //     let s_coeffs: Vec<Vec<Zq>> = self
    //         .witness
    //         .s
    //         .iter()
    //         .map(|s_i| {
    //             s_i.iter()
    //                 .flat_map(|s_i_p| *s_i_p.get_coefficients())
    //                 .collect()
    //         })
    //         .collect();

    //     for (j, &p_j) in p.iter().enumerate() {
    //         let mut poly = vec![Zq::ZERO; p.len()];
    //         for (i, s_i) in s_coeffs.iter().enumerate() {
    //             let pi_ele = &pi[i][j];
    //             let pi_ele_ca = pi_ele.conjugate_automorphism();
    //             poly = poly.add(&(pi_ele_ca.multiply(s_i)));
    //         }

    //         if poly[0] != p_j {
    //             return Err(ProverError::ProjectionError {
    //                 index: j,
    //                 expected: p_j,
    //                 computed: poly[0],
    //             });
    //         }
    //     }

    //     Ok(true)
    // }
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
        let pp = AjtaiInstances::new(&ep_1);

        // create a new prover
        let mut prover = LabradorProver::new(&pp, &witness_1, &st);
        let _: LabradorTranscript<ShakeSponge> = prover.prove(&ep_1).unwrap();
    }
}
