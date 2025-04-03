use crate::jl::ProjectionMatrix;
use crate::poly::{PolyVector, ZqVector};
use crate::utils::{
    aggregate, challenge_set::ChallengeSet, crs::PublicPrams, env_params::EnvironmentParameters,
    statement::Statement,
};
use rand::rng;

/// explicitly set the deg_bound_d to D == deg_bound_d, which is 4
const D: usize = 4;

// Proof contains the parameters will be sent to verifier
// All parameters are from tr, line 2 on page 18
pub struct Proof {
    pub u_1: PolyVector,
    pub p: ZqVector,
    pub b_ct_aggr: PolyVector,
    pub u_2: PolyVector,
    pub z: PolyVector,
    pub t_i: Vec<PolyVector>,
    pub g_ij: Vec<PolyVector>,
    pub h_ij: Vec<PolyVector>,
}

// pub struct Challenges just for testing, should be replaced by the Transcript
pub struct Challenges {
    pub pi: Vec<Vec<ZqVector>>,
    pub psi: Vec<ZqVector>,
    pub omega: Vec<ZqVector>,
    pub random_alpha: PolyVector,
    pub random_beta: PolyVector,
    pub random_c: PolyVector,
}

impl Challenges {
    pub fn new(ep: &EnvironmentParameters) -> Self {
        // generate random psi with size: k * constraint_l, each element is Zq
        let psi: Vec<ZqVector> = (0..ep.k)
            .map(|_| ZqVector::random(&mut rng(), ep.constraint_l))
            .collect();

        // generate randm omega is with size: k * lambda2, each element is Zq
        let omega: Vec<ZqVector> = (0..ep.k)
            .map(|_| ZqVector::random(&mut rng(), ep.lambda2))
            .collect();

        // \pi is from JL projection, pi contains r matrices and each matrix: security_level2 * (n*d), (security_level2 is 256 in the paper).
        let pi: Vec<Vec<ZqVector>> = (0..ep.r)
            .map(|_| ProjectionMatrix::<D>::new(ep.n).get_matrix().clone())
            .collect();

        // generate random alpha and beta from challenge set
        let cs_alpha: ChallengeSet = ChallengeSet::new(ep.deg_bound_d);
        let random_alpha: PolyVector = (0..ep.constraint_k)
            .map(|_| cs_alpha.get_challenges().clone())
            .collect();

        let cs_beta: ChallengeSet = ChallengeSet::new(ep.deg_bound_d);
        let random_beta: PolyVector = (0..ep.constraint_k)
            .map(|_| cs_beta.get_challenges().clone())
            .collect();

        let cs_c: ChallengeSet = ChallengeSet::new(ep.deg_bound_d);
        let random_c: PolyVector = (0..ep.r).map(|_| cs_c.get_challenges().clone()).collect();

        Self {
            pi,
            psi,
            omega,
            random_alpha,
            random_beta,
            random_c,
        }
    }
}
pub struct Witness {
    pub s: Vec<PolyVector>,
}

impl Witness {
    pub fn new(ep: &EnvironmentParameters) -> Self {
        let s = (0..ep.r)
            .map(|_| PolyVector::random_ternary(ep.n, ep.deg_bound_d))
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

    pub fn prove(&self, ep: &EnvironmentParameters) -> Proof {
        // check the L2 norm of the witness
        // not sure whether this should be handled during the proving or managed by the witness generator.
        let beta2 = ep.beta * ep.beta;
        self.witness.s.iter().for_each(|polys| {
            let witness_l2norm_squared = PolyVector::compute_norm_squared(polys);
            assert!(
                witness_l2norm_squared <= beta2,
                "witness l2-norm is not satisfied"
            )
        });

        // Step 1: Outer commitments u_1 starts: --------------------------------------------

        // Ajtai Commitments t_i = A * s_i
        let _matrix_a = &self.pp.matrix_a;
        let t_i = vec![PolyVector::zero()];
        let g_ij = vec![PolyVector::zero()];

        let _matrix_b = &self.pp.matrix_b;
        let _matrix_c = &self.pp.matrix_c;
        let u_1 = PolyVector::zero();

        // Step 1: Outer commitments u_1 ends: ----------------------------------------------

        // Step 2: JL projection starts: ----------------------------------------------------

        // JL projection p_j + check p_j = ct(sum(<\sigma_{-1}(pi_i^(j)), s_i>))
        let p = ZqVector::zero();

        // Step 2: JL projection ends: ------------------------------------------------------

        // Step 3: Aggregation starts: --------------------------------------------------------------

        // first aggregation
        let aggr_1 = aggregate::AggregationOne::new(self.witness, self.st, ep, self.tr);
        // second aggregation
        let aggr_2 = aggregate::AggregationTwo::new(&aggr_1, self.st, ep, self.tr);

        // Aggregation ends: ----------------------------------------------------------------

        // Step 3: Calculate h_ij, u_2, and z starts: ---------------------------------------

        let _phi_i = aggr_2.phi_i;
        let h_ij = vec![PolyVector::zero()];

        // Outer commitments: u_2
        let _matrix_d = &self.pp.matrix_d;
        let u_2 = PolyVector::zero();

        // calculate z = c_1*s_1 + ... + c_r*s_r
        let z = aggregate::calculate_z(&self.witness.s, &self.tr.random_c);

        // Step 3: Calculate h_ij, u_2, and z ends: -----------------------------------------

        Proof {
            u_1,
            p,
            b_ct_aggr: aggr_1.b_ct_aggr,
            u_2,
            z,
            t_i,
            g_ij,
            h_ij,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_prove() {
        // set up example environment, use set1 for testing.
        let ep_1 = EnvironmentParameters::set_1();
        // generate a random witness based on ep above
        let witness_1 = Witness::new(&ep_1);
        // generate public statements based on witness_1
        let st: Statement = Statement::new(&witness_1, &ep_1);
        // generate the common reference string matriices
        let pp = PublicPrams::new(&ep_1);
        // generate random challenges
        let tr = Challenges::new(&ep_1);

        // create a new prover
        let prover = LabradorProver::new(&pp, &witness_1, &st, &tr);
        let _proof = LabradorProver::prove(&prover, &ep_1);
    }
}
