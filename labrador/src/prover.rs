use crate::commitments::ajtai_commitment::AjtaiCommitment;
use crate::commitments::outer_commitments::DecompositionParameters;
use crate::commitments::outer_commitments::OuterCommitment;
use crate::core::garbage_polynomials::GarbagePolynomials;
use crate::ring::rq_matrix::RqMatrix;
use crate::ring::zq::Zq;
use crate::ring::zq::ZqVector;
use crate::{
    core::{
        aggregate,
        challenge_set::ChallengeSet,
        crs::PublicPrams,
        env_params::EnvironmentParameters,
        jl::{ProjectionMatrix, Projections},
        statement::Statement,
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
    pub p: Projections,
    pub b_ct_aggr: RqVector,
    pub u_2: RqVector,
    pub z: RqVector,
    pub t_i: Vec<RqVector>,
    pub g_ij: RqMatrix,
    pub h_ij: RqMatrix,
}

// pub struct Challenges just for testing, should be replaced by the Transcript
pub struct Challenges {
    pub pi: Vec<Vec<Vec<Zq>>>,
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

        // \pi is from JL projection, pi contains r matrices and each matrix: security_level2 * (n*d), (security_level2 is 256 in the paper).
        let pi: Vec<Vec<Vec<Zq>>> = Self::get_pi(ep.r, ep.n);

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
            pi,
            psi,
            omega,
            random_alpha,
            random_beta,
            random_c,
        }
    }

    pub fn get_pi(r: usize, n: usize) -> Vec<Vec<Vec<Zq>>> {
        (0..r)
            .map(|_| ProjectionMatrix::new(n).get_matrix().clone())
            .collect()
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
}

impl<'a> LabradorProver<'a> {
    pub fn new(
        pp: &'a PublicPrams,
        witness: &'a Witness,
        st: &'a Statement,
        tr: &'a Challenges,
    ) -> Self {
        Self {
            pp,
            witness,
            st,
            tr,
        }
    }

    /// all prove steps are from page 17
    pub fn prove(&self, ep: &EnvironmentParameters) -> Result<Proof, ProverError> {
        // check the L2 norm of the witness
        // not sure whether this should be handled during the proving or managed by the witness generator.
        Self::check_witness_l2norm(self, ep).unwrap();
        // Step 1: Outer commitments u_1 starts: --------------------------------------------

        // Ajtai Commitments t_i = A * s_i
        // let matrix_a = &self.pp.matrix_a;
        // let t_i: Vec<RqVector> = self.witness.s.iter().map(|s_i| s_i * matrix_a).collect();
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
        let mut outer_commitments = OuterCommitment::new(self.pp.clone(), ep.clone());
        outer_commitments.compute_u1(
            RqMatrix::new(t_i.clone()),
            DecompositionParameters::new(ep.b, ep.t_1),
            garbage_polynomials.g.clone(),
            DecompositionParameters::new(ep.b, ep.t_2),
        );
        // Probably do not need the following code
        // decompose t_i into t_i^(0) + ... + t_i^(t_1-1) * b_1^(t_1-1)
        // let t_ij: Vec<Vec<RqVector>> = t_i
        //     .iter()
        //     .map(|i| RqVector::decompose(i, ep.b, ep.t_1))
        //     .collect();
        // // calculate garbage polynomial g = <s_i, s_j>
        // let g_gp: Vec<RqVector> = aggregate::calculate_gij(&self.witness.s, ep.r);
        // // decompose g_gp into g_ij = g_ij^(0) + ... + g_ij^(t_2-1) * b_2^(t_2=1)
        // let g_ij: Vec<Vec<RqVector>> = g_gp
        //     .iter()
        //     .map(|i| RqVector::decompose(i, ep.b, ep.t_2))
        //     .collect();
        // let matrix_b = &self.pp.matrix_b;
        // let matrix_c = &self.pp.matrix_c;
        // calculate outer commitment u_1 = \sum(B_ik * t_i^(k)) + \sum(C_ijk * g_ij^(k))

        // let u_1 = aggregate::calculate_u_1(matrix_b, matrix_c, &t_ij, &g_ij, ep);

        // Step 1: Outer commitments u_1 ends: ----------------------------------------------

        // Step 2: JL projection starts: ----------------------------------------------------

        // JL projection p_j + check p_j = ct(sum(<\sigma_{-1}(pi_i^(j)), s_i>))
        let matrices = &self.tr.pi;
        let p = Projections::new(matrices, &self.witness.s);

        // Notice that this check is resource-intensive due to the multiplication of two ZqVector<256> instances,
        // followed by the removal of high-degree terms. It might not be a necessary check.
        Self::check_projection(self, p.get_projection()).unwrap();

        // Step 2: JL projection ends: ------------------------------------------------------

        // Step 3: Aggregation starts: --------------------------------------------------------------

        // first aggregation
        let aggr_1 = aggregate::AggregationOne::new(self.witness, self.st, ep, self.tr);
        // second aggregation
        let aggr_2 = aggregate::AggregationTwo::new(&aggr_1, self.st, ep, self.tr);

        // Aggregation ends: ----------------------------------------------------------------

        // Step 4: Calculate h_ij, u_2, and z starts: ---------------------------------------

        let phi_i = aggr_2.phi_i;
        garbage_polynomials.compute_h(&phi_i);
        // let h_gp = aggregate::calculate_hij(&phi_i, &self.witness.s, ep);
        // // decompose h_gp into h_ij = h_ij^(0) + ... + h_ij^(t_1-1) * b_1^(t_1-1)
        // let h_ij: Vec<Vec<RqVector>> = h_gp
        //     .iter()
        //     .map(|i| RqVector::decompose(i, ep.b, ep.t_1))
        //     .collect();
        // // Outer commitments: u_2
        // let matrix_d = &self.pp.matrix_d;
        // // calculate outer commitment u_2 = \sum(D_ijk * h_ij^(k))
        outer_commitments.compute_u2(
            garbage_polynomials.h.clone(),
            DecompositionParameters::new(ep.b, ep.t_1),
        );
        // let u_2 = aggregate::calculate_u_2(matrix_d, &h_ij, ep);
        // calculate z = c_1*s_1 + ... + c_r*s_r
        let z = aggregate::calculate_z(&self.witness.s, &self.tr.random_c);

        // Step 4: Calculate h_ij, u_2, and z ends: -----------------------------------------

        Ok(Proof {
            u_1: outer_commitments.u_1,
            p,
            b_ct_aggr: aggr_1.b_ct_aggr,
            u_2: outer_commitments.u_2,
            z,
            t_i,
            g_ij: garbage_polynomials.g,
            h_ij: garbage_polynomials.h,
        })
    }

    /// check p_j? = ct(sum(<σ−1(pi_i^(j)), s_i>))
    fn check_projection(&self, p: &[Zq]) -> Result<bool, ProverError> {
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
                let pi_ele = &self.tr.pi[i][j];
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

        // create a new prover
        let prover = LabradorProver::new(&pp, &witness_1, &st, &tr);
        let _proof = prover.prove(&ep_1).unwrap();
    }
}
