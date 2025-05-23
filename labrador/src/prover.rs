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
    CommitmentError(#[from] ajtai_commitment::CommitError),
    #[error("decomposition failure")]
    DecompositionError(#[from] outer_commitments::DecompositionError),
}

// Proof contains the parameters will be sent to verifier
// All parameters are from tr, line 2 on page 18
pub struct Proof {
    pub u_1: RqVector,
    pub p: Vec<Zq>,
    pub b_ct_aggr: RqVector,
    pub u_2: RqVector,
    pub z: RqVector,
    pub t_i: Vec<RqVector>,
    pub g_ij: RqMatrix,
    pub h_ij: RqMatrix,
}
pub struct Witness {
    pub s: Vec<RqVector>,
}

impl Witness {
    pub fn new(ep: &EnvironmentParameters) -> Self {
        let s = (0..ep.multiplicity)
            .map(|_| RqVector::random_ternary(&mut rng(), ep.rank))
            .collect();
        Self { s }
    }
}

pub struct LabradorProver<'a, S: Sponge> {
    pub pp: &'a AjtaiInstances,
    pub witness: &'a Witness,
    pub st: &'a Statement,
    pub transcript: LabradorTranscript<S>,
}

impl<'a, S: Sponge> LabradorProver<'a, S> {
    pub fn new(
        pp: &'a AjtaiInstances,
        witness: &'a Witness,
        st: &'a Statement,
        transcript: LabradorTranscript<S>,
    ) -> Self {
        Self {
            pp,
            witness,
            st,
            transcript,
        }
    }

    /// all prove steps are from page 17
    pub fn prove(&mut self, ep: &EnvironmentParameters) -> Result<Proof, ProverError> {
        // check the L2 norm of the witness
        // not sure whether this should be handled during the proving or managed by the witness generator.
        Self::check_witness_l2norm(self, ep)?;
        // Step 1: Outer commitments u_1 starts: --------------------------------------------

        // Ajtai Commitments t_i = A * s_i
        let t_i: Vec<RqVector> = self
            .witness
            .s
            .iter()
            .map(|s_i| self.pp.commitment_scheme_a.commit(s_i))
            .collect::<Result<_, _>>()?;

        // This replaces the following code
        let mut garbage_polynomials = GarbagePolynomials::new(self.witness.s.clone());
        garbage_polynomials.compute_g();
        // calculate outer commitment u_1 = \sum(B_ik * t_i^(k)) + \sum(C_ijk * g_ij^(k))
        let mut outer_commitments = OuterCommitment::new(self.pp);
        outer_commitments.compute_u1(
            RqMatrix::new(t_i.clone()),
            DecompositionParameters::new(ep.b, ep.t_1)?,
            garbage_polynomials.g.clone(),
            DecompositionParameters::new(ep.b, ep.t_2)?,
        );
        self.transcript.absorb_u1(outer_commitments.u_1.clone());
        // Step 1: Outer commitments u_1 ends: ----------------------------------------------

        // Step 2: JL projection starts: ----------------------------------------------------

        // JL projection p_j + check p_j = ct(sum(<\sigma_{-1}(pi_i^(j)), s_i>))
        let projections = self.transcript.generate_projections();
        let vector_p = projections.compute_batch_projection(&self.witness.s);
        self.transcript.absorb_vector_p(vector_p);
        // Projections::new(pi, &self.witness.s);

        // Notice that this check is resource-intensive due to the multiplication of two ZqVector<256> instances,
        // followed by the removal of high-degree terms. It might not be a necessary check.
        // Omid's Note: This can be removed later. However, we need to ensure a correct projection matrix with correct upper-bound.
        Self::check_projection(
            self,
            &self.transcript.vector_p,
            projections.get_projection_matrices(),
        )?;
        // Step 2: JL projection ends: ------------------------------------------------------

        // Step 3: Aggregation starts: --------------------------------------------------------------

        let size_of_psi = usize::div_ceil(ep.security_parameter, ep.log_q);
        let size_of_omega = size_of_psi;
        let vector_psi = self
            .transcript
            .generate_vector_psi(size_of_psi, ep.constraint_l);
        let vector_omega = self.transcript.generate_vector_omega(size_of_omega);
        // first aggregation
        let aggr_1 = aggregate::AggregationOne::new(
            self.witness,
            self.st,
            ep,
            projections.get_projection_matrices(),
            &vector_psi,
            &vector_omega,
        );
        self.transcript
            .absorb_vector_b_ct_aggr(aggr_1.b_ct_aggr.clone());

        // second aggregation
        let size_of_beta = size_of_psi;
        let alpha_vector = self.transcript.generate_rq_vector(ep.constraint_k);
        let beta_vector = self.transcript.generate_rq_vector(size_of_beta);
        let aggr_2 =
            aggregate::AggregationTwo::new(&aggr_1, self.st, ep, &alpha_vector, &beta_vector);
        // Aggregation ends: ----------------------------------------------------------------

        // Step 4: Calculate h_ij, u_2, and z starts: ---------------------------------------

        let phi_i = aggr_2.phi_i;
        garbage_polynomials.compute_h(&phi_i);
        outer_commitments.compute_u2(
            garbage_polynomials.h.clone(),
            DecompositionParameters::new(ep.b, ep.t_1)?,
        );
        self.transcript.absorb_u2(outer_commitments.u_2);

        // calculate z = c_1*s_1 + ... + c_r*s_r
        let challenges = self.transcript.generate_challenges(ep.operator_norm);
        let z = aggregate::calculate_z(&self.witness.s, &challenges);

        // Step 4: Calculate h_ij, u_2, and z ends: -----------------------------------------

        Ok(Proof {
            u_1: self.transcript.u1.clone(),
            p: self.transcript.vector_p.clone(),
            b_ct_aggr: self.transcript.b_ct_aggr.clone(),
            u_2: self.transcript.u2.clone(),
            z,
            t_i,
            g_ij: garbage_polynomials.g,
            h_ij: garbage_polynomials.h,
        })
    }

    /// check p_j? = ct(sum(<σ−1(pi_i^(j)), s_i>))
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

    /// check the L2 norm of the witness, || s_i || <= beta
    fn check_witness_l2norm(&self, ep: &EnvironmentParameters) -> Result<(), ProverError> {
        let beta2 = ep.beta * ep.beta;
        for polys in &self.witness.s {
            let witness_l2norm_squared = RqVector::compute_norm_squared(polys);
            if witness_l2norm_squared > beta2 {
                return Err(ProverError::WitnessL2NormViolated {
                    norm_squared: witness_l2norm_squared,
                    allowed: beta2,
                });
            }
        }
        Ok(())
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
        // generate random challenges used between prover and verifier.
        let transcript = LabradorTranscript::new(
            ShakeSponge::default(),
            ep_1.security_parameter,
            ep_1.rank,
            ep_1.multiplicity,
        );

        // create a new prover
        let mut prover = LabradorProver::new(&pp, &witness_1, &st, transcript);
        let _proof = prover.prove(&ep_1).unwrap();
    }
}
