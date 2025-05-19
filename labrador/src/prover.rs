use crate::commitments::ajtai_commitment;
use crate::commitments::common_instances::AjtaiInstances;
use crate::commitments::outer_commitments;
use crate::commitments::outer_commitments::DecompositionParameters;
use crate::commitments::outer_commitments::OuterCommitment;
use crate::core::garbage_polynomials::GarbagePolynomials;
use crate::ring::rq_matrix::RqMatrix;
use crate::ring::zq::Zq;
use crate::ring::zq::ZqVector;
use crate::transcript::{LabradorTranscript, Sponge};
use crate::{
    core::{aggregate, env_params::EnvironmentParameters, statement::Statement},
    ring::rq_vector::RqVector,
};
use rand::rng;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum ProverError {
    ProjectionError {
        index: usize,
        expected: Zq,
        computed: Zq,
    },
    #[error("commitment failure")]
    CommitmentError(#[from] ajtai_commitment::CommitError),
    #[error("decomposition failure")]
    DecompositionError(#[from] outer_commitments::DecompositionError),
}

pub struct Witness {
    pub s: Vec<RqVector>,
}

impl Witness {
    pub fn new(ep: &EnvironmentParameters) -> Self {
        loop {
            let s: Vec<RqVector> = (0..ep.multiplicity)
                .map(|_| RqVector::random_ternary(&mut rng(), ep.rank))
                .collect();
            if Self::validate_l2_norm(&s, ep) {
                return Self { s };
            }
        }
    }

    fn validate_l2_norm(candidate: &[RqVector], ep: &EnvironmentParameters) -> bool {
        let beta2 = ep.beta * ep.beta;
        for polys in candidate {
            let witness_l2norm_squared = RqVector::compute_norm_squared(polys);
            if witness_l2norm_squared > beta2 {
                return false;
            }
        }
        true
    }
}

pub struct LabradorProver<'a, S: Sponge> {
    pub pp: &'a AjtaiInstances,
    pub witness: &'a Witness,
    pub st: &'a Statement,
}

impl<'a> LabradorProver<'a> {
    pub fn new(pp: &'a AjtaiInstances, witness: &'a Witness, st: &'a Statement) -> Self {
        Self { pp, witness, st }
    }

    fn compute_vector_ti(&self) -> RqMatrix {
        RqMatrix::new(
            self.witness
                .s
                .iter()
                .map(|s_i| {
                    self.pp
                        .commitment_scheme_a
                        .commit(s_i)
                        .expect("Commitment error in committing to s_i")
                })
                .collect(),
        )
    }

    /// all prove steps are from page 17
    pub fn prove(
        &mut self,
        ep: &EnvironmentParameters,
    ) -> Result<LabradorTranscript<ShakeSponge>, ProverError> {
        // Generate random challenges used between prover and verifier
        let mut transcript = LabradorTranscript::new(ShakeSponge::default());
        // Compute garbage polynomials g and h
        let mut garbage_polynomials = GarbagePolynomials::new(&self.witness.s);
        // Compute outer commitments u1 and u2
        let mut outer_commitments = OuterCommitment::new(self.pp);

        // Step 1: Outer commitments u_1 starts: --------------------------------------------

        // Ajtai Commitments t_i = A * s_i
        let t_i = self.compute_vector_ti();
        // g_ij = <s_i, s_j>
        garbage_polynomials.compute_g();
        // calculate outer commitment u_1 = \sum(B_ik * t_i^(k)) + \sum(C_ijk * g_ij^(k))
        let commitment_u1 = outer_commitments.compute_u1(
            &t_i,
            DecompositionParameters::new(ep.b, ep.t_1)
                .expect("Decomposition error in decomposing t"),
            &garbage_polynomials.g,
            DecompositionParameters::new(ep.b, ep.t_2)
                .expect("Decomposition error in decomposing g"),
        );
        transcript.set_u1(commitment_u1);
        // Step 1: Outer commitments u_1 ends: ----------------------------------------------

        // Step 2: JL projection starts: ----------------------------------------------------
        let vector_of_projection_matrices = transcript.generate_vector_of_projection_matrices(
            ep.security_parameter,
            ep.rank,
            ep.multiplicity,
        );
        let jl_projection_instance =
            jl::Projection::new(vector_of_projection_matrices, ep.security_parameter);
        let vector_p = jl_projection_instance.compute_batch_projection(&self.witness.s);
        transcript.set_vector_p(vector_p);
        // Notice that this check is resource-intensive due to the multiplication of two ZqVector<256> instances,
        // followed by the removal of high-degree terms. It might not be a necessary check.
        // Omid's Note: This can be removed later. However, we need to ensure a correct projection matrix with correct upper-bound.
        Self::check_projection(
            self,
            &transcript.vector_p,
            jl_projection_instance.get_random_linear_map_vector(),
        )
        .expect("Projection check failed");
        // Step 2: JL projection ends: ------------------------------------------------------

        // Step 3: Aggregation starts: --------------------------------------------------------------
        let size_of_psi = usize::div_ceil(ep.security_parameter, ep.log_q);
        let size_of_omega = size_of_psi;
        let vector_psi = transcript.generate_vector_psi(size_of_psi, ep.constraint_l);
        let vector_omega = transcript.generate_vector_omega(size_of_omega, ep.security_parameter);
        // first aggregation
        let a_ct_aggr = aggregate::calculate_aggr_ct_a(&vector_psi, &self.st.a_ct, ep);
        let phi_ct_aggr = aggregate::calculate_aggr_ct_phi(
            &self.st.phi_ct,
            jl_projection_instance.get_random_linear_map_vector(),
            &vector_psi,
            &vector_omega,
            ep,
        );
        let b_ct_aggr =
            aggregate::calculate_aggr_ct_b(&a_ct_aggr, &phi_ct_aggr, &self.witness.s, ep);
        transcript.set_vector_b_ct_aggr(b_ct_aggr);

        // second aggregation
        let size_of_beta = size_of_psi;
        let alpha_vector = transcript.generate_rq_vector(ep.constraint_k);
        let beta_vector = transcript.generate_rq_vector(size_of_beta);
        let phi_i = aggregate::calculate_aggr_phi(
            &self.st.phi_constraint,
            &phi_ct_aggr,
            &alpha_vector,
            &beta_vector,
            ep,
        );
        // Aggregation ends: ----------------------------------------------------------------

        // Step 4: Calculate h_ij, u_2, and z starts: ---------------------------------------
        garbage_polynomials.compute_h(&phi_i);
        let commitment_u2 = outer_commitments.compute_u2(
            &garbage_polynomials.h,
            DecompositionParameters::new(ep.b, ep.t_1)
                .expect("Decomposition error in decomposing h"),
        );
        transcript.set_u2(commitment_u2);

        // calculate z = c_1*s_1 + ... + c_r*s_r
        let challenges = transcript.generate_challenges(ep.operator_norm, ep.multiplicity);
        let z = aggregate::calculate_z(&self.witness.s, &challenges);
        transcript.set_recursive_part(z, t_i, garbage_polynomials.g, garbage_polynomials.h);

        // Step 4: Calculate h_ij, u_2, and z ends: -----------------------------------------

        Ok(transcript)
    }

    // The following is not part of the prover.
    // Todo: Add jl projection constraints to the statement.
    // /// check p_j? = ct(sum(<σ−1(pi_i^(j)), s_i>))
    fn check_projection(&self, p: &[Zq], pi: &[Vec<Vec<Zq>>]) -> Result<bool, ProverError> {
        let s_coeffs: Vec<Vec<Zq>> = self
            .witness
            .s
            .iter()
            .map(|s_i| {
                s_i.iter()
                    .flat_map(|s_i_p| *s_i_p.get_coefficients())
                    .collect()
            })
            .collect();

        for (j, &p_j) in p.iter().enumerate() {
            let mut poly = vec![Zq::ZERO; p.len()];
            for (i, s_i) in s_coeffs.iter().enumerate() {
                let pi_ele = &pi[i][j];
                let pi_ele_ca = pi_ele.conjugate_automorphism();
                poly = poly.add(&(pi_ele_ca.multiply(s_i)));
            }

            if poly[0] != p_j {
                return Err(ProverError::ProjectionError {
                    index: j,
                    expected: p_j,
                    computed: poly[0],
                });
            }
        }

        Ok(true)
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
        let witness_1 = Witness::new(&ep_1);
        // generate public statement based on witness_1
        let st: Statement = Statement::new(&witness_1, &ep_1);
        // generate the common reference string matrices A, B, C, D
        let pp = AjtaiInstances::new(&ep_1);

        // create a new prover
        let mut prover = LabradorProver::new(&pp, &witness_1, &st);
        let _ = prover.prove(&ep_1).unwrap();
    }
}
