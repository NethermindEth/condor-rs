use crate::commitments::ajtai_commitment::AjtaiCommitment;
use crate::commitments::outer_commitments::DecompositionParameters;
use crate::commitments::outer_commitments::OuterCommitment;
use crate::core::garbage_polynomials::GarbagePolynomials;
use crate::core::jl;
use crate::ring::rq_matrix::RqMatrix;
use crate::ring::zq::Zq;
use crate::ring::zq::ZqVector;
use crate::transcript::lib::LabradorTranscript;
use crate::transcript::shake_sponge::ShakeSponge;
use crate::{
    core::{
        aggregate, challenge_set::ChallengeSet, crs::PublicPrams,
        env_params::EnvironmentParameters, statement::Statement,
    },
    ring::rq_vector::RqVector,
};
use rand::rng;

#[derive(Debug)]
pub enum ProverError {
    /// Indicates that the L2 norm (squared) of the witness exceeded the allowed threshold.
    WitnessL2NormViolated { norm_squared: Zq, allowed: Zq },
    ProjectionError {
        index: usize,
        expected: Zq,
        computed: Zq,
    },
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

// pub struct Challenges just for testing, should be replaced by the Transcript
pub struct Challenges {
    pub psi: Vec<Vec<Zq>>,
    pub omega: Vec<Vec<Zq>>,
    pub random_alpha: RqVector,
    pub random_beta: RqVector,
    pub random_c: RqVector,
}

impl Challenges {
    pub fn new(ep: &EnvironmentParameters) -> Self {
        // generate random psi with size: k * constraint_l, each element is Zq
        let psi: Vec<Vec<Zq>> = (0..ep.kappa)
            .map(|_| Vec::<Zq>::random(&mut rng(), ep.constraint_l))
            .collect();

        // generate randm omega is with size: k * lambda2, each element is Zq
        let omega: Vec<Vec<Zq>> = (0..ep.kappa)
            .map(|_| Vec::<Zq>::random(&mut rng(), 2 * ep.lambda))
            .collect();

        // generate random alpha and beta from challenge set
        let cs_alpha: ChallengeSet = ChallengeSet::new();
        let random_alpha: RqVector = (0..ep.constraint_k)
            .map(|_| *cs_alpha.get_challenges())
            .collect();

        let cs_beta: ChallengeSet = ChallengeSet::new();
        let random_beta: RqVector = (0..ep.constraint_k)
            .map(|_| *cs_beta.get_challenges())
            .collect();

        let cs_c: ChallengeSet = ChallengeSet::new();
        let random_c: RqVector = (0..ep.r).map(|_| *cs_c.get_challenges()).collect();

        Self {
            psi,
            omega,
            random_alpha,
            random_beta,
            random_c,
        }
    }
}
pub struct Witness {
    pub s: Vec<RqVector>,
}

impl Witness {
    pub fn new(ep: &EnvironmentParameters) -> Self {
        let s = (0..ep.r)
            .map(|_| RqVector::random_ternary(&mut rng(), ep.n))
            .collect();
        Self { s }
    }
}

pub struct LabradorProver<'a> {
    pub pp: &'a PublicPrams,
    pub witness: &'a Witness,
    pub st: &'a Statement,
    pub tr: &'a Challenges,
    pub transcript: LabradorTranscript<ShakeSponge>,
}

impl<'a> LabradorProver<'a> {
    pub fn new(
        pp: &'a PublicPrams,
        witness: &'a Witness,
        st: &'a Statement,
        tr: &'a Challenges,
        transcript: LabradorTranscript<ShakeSponge>,
    ) -> Self {
        Self {
            pp,
            witness,
            st,
            tr,
            transcript,
        }
    }

    /// all prove steps are from page 17
    pub fn prove(&mut self, ep: &EnvironmentParameters) -> Result<Proof, ProverError> {
        // check the L2 norm of the witness
        // not sure whether this should be handled during the proving or managed by the witness generator.
        Self::check_witness_l2norm(self, ep).unwrap();
        // Step 1: Outer commitments u_1 starts: --------------------------------------------

        // Ajtai Commitments t_i = A * s_i
        let t_i: Vec<RqVector> = self
            .witness
            .s
            .iter()
            .map(|s_i| {
                AjtaiCommitment::new(ep.beta, ep.beta, self.pp.matrix_a.clone())
                    .unwrap()
                    .commit(s_i)
                    .unwrap()
            })
            .collect();

        // This replaces the following code
        let mut garbage_polynomials = GarbagePolynomials::new(self.witness.s.clone());
        garbage_polynomials.compute_g();
        // calculate outer commitment u_1 = \sum(B_ik * t_i^(k)) + \sum(C_ijk * g_ij^(k))
        let mut outer_commitments = OuterCommitment::new(self.pp.clone(), ep.clone());
        outer_commitments.compute_u1(
            RqMatrix::new(t_i.clone()),
            DecompositionParameters::new(ep.b, ep.t_1).unwrap(),
            garbage_polynomials.g.clone(),
            DecompositionParameters::new(ep.b, ep.t_2).unwrap(),
        );
        self.transcript.absorb_u1(outer_commitments.u_1.clone());
        // Step 1: Outer commitments u_1 ends: ----------------------------------------------

        // Step 2: JL projection starts: ----------------------------------------------------

        // JL projection p_j + check p_j = ct(sum(<\sigma_{-1}(pi_i^(j)), s_i>))
        let vector_of_projection_matrices =
            self.transcript.generate_vector_of_projection_matrices();
        let vector_p = jl::Projection::new(vector_of_projection_matrices.clone(), ep.lambda)
            .compute_batch_projection(&self.witness.s);
        self.transcript.absorb_vector_p(vector_p);
        // Projections::new(pi, &self.witness.s);

        // Notice that this check is resource-intensive due to the multiplication of two ZqVector<256> instances,
        // followed by the removal of high-degree terms. It might not be a necessary check.
        // Omid's Note: This can be removed later. However, we need to ensure a correct projection matrix with correct upper-bound.
        Self::check_projection(
            self,
            &self.transcript.vector_p,
            vector_of_projection_matrices.clone(),
        )
        .unwrap();
        // Step 2: JL projection ends: ------------------------------------------------------

        // Step 3: Aggregation starts: --------------------------------------------------------------

        // first aggregation
        let aggr_1 = aggregate::AggregationOne::new(
            self.witness,
            self.st,
            ep,
            self.tr,
            &vector_of_projection_matrices,
        );
        // second aggregation
        let aggr_2 = aggregate::AggregationTwo::new(&aggr_1, self.st, ep, self.tr);

        // Aggregation ends: ----------------------------------------------------------------

        // Step 4: Calculate h_ij, u_2, and z starts: ---------------------------------------

        let phi_i = aggr_2.phi_i;
        garbage_polynomials.compute_h(&phi_i);
        outer_commitments.compute_u2(
            garbage_polynomials.h.clone(),
            DecompositionParameters::new(ep.b, ep.t_1).unwrap(),
        );

        // calculate z = c_1*s_1 + ... + c_r*s_r
        let z = aggregate::calculate_z(&self.witness.s, &self.tr.random_c);

        // Step 4: Calculate h_ij, u_2, and z ends: -----------------------------------------

        Ok(Proof {
            u_1: outer_commitments.u_1,
            p: self.transcript.vector_p.clone(),
            b_ct_aggr: aggr_1.b_ct_aggr,
            u_2: outer_commitments.u_2,
            z,
            t_i,
            g_ij: garbage_polynomials.g,
            h_ij: garbage_polynomials.h,
        })
    }

    /// check p_j? = ct(sum(<σ−1(pi_i^(j)), s_i>))
    fn check_projection(&self, p: &[Zq], pi: Vec<Vec<Vec<Zq>>>) -> Result<bool, ProverError> {
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
    fn check_witness_l2norm(&self, ep: &EnvironmentParameters) -> Result<bool, ProverError> {
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
        Ok(true)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_prove() {
        // set up example environment parameters, use default set for testing.
        let ep_1 = EnvironmentParameters::default();
        // generate a random witness based on environment parameters above
        let witness_1 = Witness::new(&ep_1);
        // generate public statement based on witness_1
        let st: Statement = Statement::new(&witness_1, &ep_1);
        // generate the common reference string matrices A, B, C, D
        let pp = PublicPrams::new(&ep_1);
        // generate random challenges used between prover and verifier.
        let tr = Challenges::new(&ep_1);
        let transcript =
            LabradorTranscript::new(ShakeSponge::default(), ep_1.lambda, ep_1.n, ep_1.r);

        // create a new prover
        let mut prover = LabradorProver::new(&pp, &witness_1, &st, &tr, transcript);
        let _proof = prover.prove(&ep_1).unwrap();
    }
}
