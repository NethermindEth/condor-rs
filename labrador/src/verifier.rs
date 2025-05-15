#![allow(clippy::result_large_err)]

use crate::commitments::common_instances::AjtaiInstances;
use crate::commitments::outer_commitments::{DecompositionParameters, OuterCommitment};
use crate::core::{aggregate, env_params::EnvironmentParameters, statement::Statement};
use crate::prover::{Challenges, Proof};
use crate::ring::rq::Rq;
use crate::ring::rq_matrix::RqMatrix;
use crate::ring::rq_vector::RqVector;
use crate::ring::zq::{Zq, ZqVector};

#[derive(Debug)]
pub enum VerifierError {
    NotSymmetric {
        i: usize,
        j: usize,
        expected: Rq,
        found: Rq,
    },
    B0Mismatch {
        index: usize,
        expected: Zq,
        computed: Zq,
    },
    NormSumExceeded {
        norm: Zq,
        allowed: Zq,
    },
    AzError {
        computed: RqVector,
        expected: RqVector,
    },
    ZInnerError {
        computed: Rq,
        expected: Rq,
    },
    PhiError {
        computed: Rq,
        expected: Rq,
    },
    RelationCheckFailed,
    OuterCommitError {
        computed: RqVector,
        expected: RqVector,
    },
}
pub struct LabradorVerifier<'a> {
    pub pp: &'a AjtaiInstances,
    pub st: &'a Statement,
    pub tr: &'a Challenges,
}

impl<'a> LabradorVerifier<'a> {
    pub fn new(pp: &'a AjtaiInstances, st: &'a Statement, tr: &'a Challenges) -> Self {
        Self { pp, st, tr }
    }

    /// All check conditions are from page 18
    pub fn verify(&self, proof: &Proof, ep: &EnvironmentParameters) -> Result<bool, VerifierError> {
        // check b_0^{''(k)} ?= <omega^(k),p> + \sum(psi_l^(k) * b_0^{'(l)})
        Self::check_b_0_aggr(self, proof, ep).unwrap();

        // 3. line 14: check norm_sum(z, t, g, h) <= (beta')^2

        // decompose z into z = z^(0) + z^(1) * b, only two parts.
        let z_ij = RqVector::decompose(&proof.z, ep.b, 2);
        let t_ij: Vec<Vec<RqVector>> = proof
            .t_i
            .iter()
            .map(|i| RqVector::decompose(i, ep.b, ep.t_1))
            .collect();
        let g_ij: Vec<Vec<RqVector>> = proof
            .g_ij
            .get_elements()
            .iter()
            .map(|i| RqVector::decompose(i, ep.b, ep.t_2))
            .collect();
        let h_ij: Vec<Vec<RqVector>> = proof
            .h_ij
            .get_elements()
            .iter()
            .map(|i| RqVector::decompose(i, ep.b, ep.t_1))
            .collect();
        let norm_z_ij = z_ij
            .iter()
            .fold(Zq::ZERO, |acc, p| acc + p.compute_norm_squared());
        let norm_t_ij = Self::norm_squared(&t_ij);
        let norm_g_ij = Self::norm_squared(&g_ij);
        let norm_h_ij = Self::norm_squared(&h_ij);
        let norm_sum = norm_z_ij + norm_t_ij + norm_g_ij + norm_h_ij;

        if norm_sum > ep.beta * ep.beta {
            return Err(VerifierError::NormSumExceeded {
                norm: norm_sum,
                allowed: ep.beta * ep.beta,
            });
        }

        // 4. line 15: check Az ?= c_1 * t_1 + ... + c_r * t_r

        let az = self.pp.commitment_scheme_a.matrix() * &proof.z;
        let ct_sum = aggregate::calculate_z(&proof.t_i, &self.tr.random_c);

        if az != ct_sum {
            return Err(VerifierError::AzError {
                computed: az,
                expected: ct_sum,
            });
        }

        // 5. lne 16: check <z, z> ?= \sum(g_ij * c_i * c_j)

        let z_inner = proof.z.inner_product_poly_vector(&proof.z);
        let sum_gij_cij = Self::calculate_gh_ci_cj(&proof.g_ij, &self.tr.random_c, ep.r);

        if z_inner != sum_gij_cij {
            return Err(VerifierError::ZInnerError {
                computed: z_inner,
                expected: sum_gij_cij,
            });
        }

        // 6. line 17: check \sum(<\phi_i, z>c_i) ?= \sum(h_ij * c_i * c_j)

        let phi_ct_aggr = aggregate::AggregationOne::get_phi_ct_aggr(
            &self.st.phi_ct,
            &self.tr.pi,
            &self.tr.psi,
            &self.tr.omega,
            ep,
        );
        let phi_i = aggregate::AggregationTwo::get_phi_i(
            &self.st.phi_constraint,
            &phi_ct_aggr,
            &self.tr.random_alpha,
            &self.tr.random_beta,
            ep,
        );
        let sum_phi_z_c = Self::calculate_phi_z_c(&phi_i, &self.tr.random_c, &proof.z);
        let sum_hij_cij = Self::calculate_gh_ci_cj(&proof.h_ij, &self.tr.random_c, ep.r);

        // Left side multiple by 2 because of when we calculate h_ij, we didn't apply the division (divided by 2)
        if &sum_phi_z_c * &Zq::TWO != sum_hij_cij {
            return Err(VerifierError::PhiError {
                computed: &sum_phi_z_c * &Zq::TWO,
                expected: sum_hij_cij,
            });
        }

        // 7. line 18: check \sum(a_ij * g_ij) + \sum(h_ii) - b ?= 0

        let a_ct_aggr = aggregate::AggregationOne::get_a_ct_aggr(&self.tr.psi, &self.st.a_ct, ep);
        let a_primes = aggregate::AggregationTwo::get_a_i(
            &self.st.a_constraint,
            &a_ct_aggr,
            &self.tr.random_alpha,
            &self.tr.random_beta,
            ep,
        );
        let b_primes = aggregate::AggregationTwo::get_b_i(
            &self.st.b_constraint,
            &proof.b_ct_aggr,
            &self.tr.random_alpha,
            &self.tr.random_beta,
            ep,
        );

        if !Self::check_relation(
            &RqMatrix::new(a_primes),
            &b_primes,
            &proof.g_ij,
            &proof.h_ij,
        ) {
            return Err(VerifierError::RelationCheckFailed);
        }

        // 8. line 19: u_1 ?= \sum(\sum(B_ik * t_i^(k))) + \sum(\sum(C_ijk * g_ij^(k)))

        let u_1 = &proof.u_1;
        let mut outer_commitments = OuterCommitment::new(self.pp);
        outer_commitments.compute_u1(
            RqMatrix::new(proof.t_i.clone()),
            DecompositionParameters::new(ep.b, ep.t_1).unwrap(),
            proof.g_ij.clone(),
            DecompositionParameters::new(ep.b, ep.t_2).unwrap(),
        );

        if proof.u_1 != outer_commitments.u_1 {
            return Err(VerifierError::OuterCommitError {
                computed: u_1.clone(),
                expected: outer_commitments.u_1,
            });
        }

        // 9. line 20: u_2 ?= \sum(\sum(D_ijk * h_ij^(k)))
        outer_commitments.compute_u2(
            proof.h_ij.clone(),
            DecompositionParameters::new(ep.b, ep.t_1).unwrap(),
        );

        if proof.u_2 != outer_commitments.u_2 {
            return Err(VerifierError::OuterCommitError {
                computed: outer_commitments.u_2.clone(),
                expected: proof.u_2.clone(),
            });
        }

        Ok(true)
    }

    /// calculate the right hand side of line 16 or line 17, \sum(g_ij * c_i * c_j) or \sum(h_ij * c_i * c_j)
    fn calculate_gh_ci_cj(x_ij: &RqMatrix, random_c: &RqVector, r: usize) -> Rq {
        (0..r)
            .map(|i| {
                (0..r)
                    .map(|j| {
                        (x_ij.get_cell_symmetric(i, j) * random_c.get_elements()[i].clone())
                            * random_c.get_elements()[j].clone()
                    })
                    .fold(Rq::zero(), |acc, x| acc + x)
            })
            .fold(Rq::zero(), |acc, x| acc + x)
    }

    /// calculate the left hand side of line 17, \sum(<\phi_z, z> * c_i)
    fn calculate_phi_z_c(phi: &[RqVector], c: &RqVector, z: &RqVector) -> Rq {
        phi.iter()
            .zip(c.iter())
            .map(|(phi_i, c_i)| &(phi_i.inner_product_poly_vector(z)) * c_i)
            .fold(Rq::zero(), |acc, x| acc + x)
    }

    fn norm_squared(polys: &[Vec<RqVector>]) -> Zq {
        polys.iter().fold(Zq::ZERO, |acc, poly| {
            acc + poly
                .iter()
                .fold(Zq::ZERO, |acc, p| acc + p.compute_norm_squared())
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
                sum_a_primes_g += a_primes.get_elements()[i].get_elements()[j].clone()
                    * g.get_cell_symmetric(i, j);
            }
        }

        let sum_h_ii = (0..r).fold(Rq::zero(), |acc, i| acc + h.get_cell_symmetric(i, i));

        let b_primes2 = b_primes * &Zq::TWO;
        let sum_a_primes_g2 = &sum_a_primes_g * &Zq::TWO;

        sum_a_primes_g2 + sum_h_ii == b_primes2
    }

    fn check_b_0_aggr(
        &self,
        proof: &Proof,
        ep: &EnvironmentParameters,
    ) -> Result<bool, VerifierError> {
        for k in 0..ep.kappa {
            let b_0_poly = proof.b_ct_aggr.get_elements()[k].get_coefficients()[0];
            let mut b_0: Zq = (0..ep.constraint_l)
                .map(|l| self.tr.psi[k][l] * self.st.b_0_ct[l])
                .sum();
            let inner_omega_p = self.tr.omega[k].inner_product(proof.p.get_projection());
            b_0 += inner_omega_p;
            if b_0 != b_0_poly {
                return Err(VerifierError::B0Mismatch {
                    index: k,
                    expected: b_0_poly,
                    computed: b_0,
                });
            }
        }

        Ok(true)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prover::{LabradorProver, Witness};

    #[test]
    fn test_verify() {
        // set up example environment, use default set for testing.
        let ep_1 = EnvironmentParameters::default();
        // generate a random witness based on ep above
        let witness_1 = Witness::new(&ep_1);
        // generate public statements based on witness_1
        let st: Statement = Statement::new(&witness_1, &ep_1);
        // generate the common reference string matrices
        let pp = AjtaiInstances::new(&ep_1);
        // generate random challenges
        let tr = Challenges::new(&ep_1);

        // create a new prover
        let prover = LabradorProver::new(&pp, &witness_1, &st, &tr);
        let proof = prover.prove(&ep_1).unwrap();

        // create a new verifier
        let verifier = LabradorVerifier::new(&pp, &st, &tr);
        let result = verifier.verify(&proof, &ep_1);
        assert!(result.unwrap());
    }
}
