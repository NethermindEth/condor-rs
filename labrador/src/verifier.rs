use crate::core::{
    aggregate, crs::PublicPrams, env_params::EnvironmentParameters, statement::Statement,
};
use crate::prover::{Challenges, Proof};
use crate::ring::poly::{PolyRing, PolyVector};
use crate::ring::zq::Zq;

pub struct LabradorVerifier<'a> {
    pub pp: &'a PublicPrams,
    pub st: &'a Statement,
    pub tr: &'a Challenges,
}

impl<'a> LabradorVerifier<'a> {
    pub fn new(pp: &'a PublicPrams, st: &'a Statement, tr: &'a Challenges) -> Self {
        Self { pp, st, tr }
    }

    /// All check conditions are from page 18
    pub fn verify(&self, proof: &Proof, ep: &EnvironmentParameters) -> bool {
        // 1. line 08: check g_ij ?= g_ji
        // 2. line 09: check h_ij ?= h_ji

        for i in 0..ep.r {
            for j in (i + 1)..ep.r {
                assert_eq!(
                    proof.g_ij[i].get_elements()[j],
                    proof.g_ij[j].get_elements()[i]
                );
                assert_eq!(
                    proof.h_ij[i].get_elements()[j],
                    proof.h_ij[j].get_elements()[i]
                );
            }
        }

        // 3. line 14: check norm_sum(z, t, g, h) <= (beta')^2

        // decompose z into z = z^(0) + z^(1) * b, only two parts.

        let z_ij = PolyVector::decompose(&proof.z, ep.b, 2);
        let t_ij: Vec<Vec<PolyVector>> = proof
            .t_i
            .iter()
            .map(|i| PolyVector::decompose(i, ep.b, ep.t_1))
            .collect();
        let g_ij: Vec<Vec<PolyVector>> = proof
            .g_ij
            .iter()
            .map(|i| PolyVector::decompose(i, ep.b, ep.t_2))
            .collect();
        let h_ij: Vec<Vec<PolyVector>> = proof
            .h_ij
            .iter()
            .map(|i| PolyVector::decompose(i, ep.b, ep.t_1))
            .collect();
        let norm_z_ij = z_ij
            .iter()
            .fold(Zq::ZERO, |acc, p| acc + p.compute_norm_squared());
        let norm_t_ij = Self::norm_squared(&t_ij);
        let norm_g_ij = Self::norm_squared(&g_ij);
        let norm_h_ij = Self::norm_squared(&h_ij);
        println!("{:?}", norm_z_ij);
        println!("{:?}", norm_t_ij);
        println!("{:?}", norm_g_ij);
        println!("{:?}", norm_h_ij);
        let norm_sum = norm_z_ij + norm_t_ij + norm_g_ij + norm_h_ij;
        println!("{:?}", norm_sum);

        assert!(norm_sum <= ep.beta * ep.beta);

        // 4. line 15: check Az ?= c_1 * t_1 + ... + c_r * t_r

        let az = &proof.z * &self.pp.matrix_a;
        let ct_sum: PolyVector = aggregate::calculate_z(&proof.t_i, &self.tr.random_c);

        assert_eq!(az, ct_sum);

        // 5. lne 16: check <z, z> ?= \sum(g_ij * c_i * c_j)

        let z_inner = &proof.z.inner_product_poly_vector(&proof.z);
        let sum_gij_cij = Self::calculate_gh_ci_cj(&proof.g_ij, &self.tr.random_c, ep.r);

        assert_eq!(z_inner, &sum_gij_cij);

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
        assert_eq!(&sum_phi_z_c * &Zq::TWO, sum_hij_cij);

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

        assert!(Self::check_relation(
            &a_primes,
            &b_primes,
            &proof.g_ij,
            &proof.h_ij
        ));

        // 8. line 19: u_1 ?= \sum(\sum(B_ik * t_i^(k))) + \sum(\sum(C_ijk * g_ij^(k)))

        let u_1 = &proof.u_1;
        let outer_commit_u_1 =
            aggregate::calculate_u_1(&self.pp.matrix_b, &self.pp.matrix_c, &t_ij, &g_ij, ep);

        assert_eq!(u_1, &outer_commit_u_1);

        // 9. line 20: u_2 ?= \sum(\sum(D_ijk * h_ij^(k)))

        let u_2 = &proof.u_2;
        let outer_commit_u_2 = aggregate::calculate_u_2(&self.pp.matrix_d, &h_ij, ep);

        assert_eq!(u_2, &outer_commit_u_2);

        true
    }

    #[rustfmt::skip]
    fn calculate_gh_ci_cj(x_ij: &[PolyVector], random_c: &PolyVector, r: usize) -> PolyRing {
        (0..r).map(|i| {
            (0..r).map(|j| {
                &(&x_ij[i].get_elements()[j] * &random_c.get_elements()[i])
                    * &random_c.get_elements()[j]
                }).fold(PolyRing::zero_poly(), |acc, x| &acc + &x)
            }).fold(PolyRing::zero_poly(), |acc, x| &acc + &x)
    }

    fn calculate_phi_z_c(phi: &[PolyVector], c: &PolyVector, z: &PolyVector) -> PolyRing {
        phi.iter()
            .zip(c.iter())
            .map(|(phi_i, c_i)| &(phi_i.inner_product_poly_vector(z)) * c_i)
            .fold(PolyRing::zero_poly(), |acc, x| &acc + &x)
    }

    fn norm_squared(polys: &[Vec<PolyVector>]) -> Zq {
        polys.iter().fold(Zq::ZERO, |acc, poly| {
            acc + poly
                .iter()
                .fold(Zq::ZERO, |acc, p| acc + p.compute_norm_squared())
        })
    }

    /// line 18: check if \sum(a_{ij} * g_{ij}) + \sum(h_{ii}) - b ?= 0
    /// in the verifier process, page 18 from the paper.
    ///
    /// param: a_primes: a_{ij}^{''(k)}
    /// param: b_primes: b^{''(k)}
    /// param: g: g_{ij}
    /// param: h: h_{ii}
    ///
    /// return: true if the relation holds, false otherwise
    pub fn check_relation(
        a_primes: &[PolyVector],
        b_primes: &PolyRing,
        g: &[PolyVector],
        h: &[PolyVector],
    ) -> bool {
        let r = a_primes.len();
        let d = a_primes[0].get_elements()[0].get_coeffs().len();

        let sum_a_primes_g: PolyRing = a_primes
            .iter()
            .zip(g.iter())
            .map(|(a_i, g_i)| {
                a_i.iter()
                    .zip(g_i.iter())
                    .map(|(a_ij, g_ij)| a_ij * g_ij)
                    .fold(PolyRing::new(vec![Zq::ZERO; d]), |acc, val| &acc + &val)
            })
            .fold(PolyRing::new(vec![Zq::ZERO; d]), |acc, val| &acc + &val);

        let sum_h_ii: PolyRing = (0..r).fold(PolyRing::new(vec![Zq::ZERO; d]), |acc, i| {
            &acc + &h[i].get_elements()[i]
        });

        let b_primes2 = b_primes * &Zq::TWO;
        let sum_a_primes_g2 = &sum_a_primes_g * &Zq::TWO;

        &sum_a_primes_g2 + &sum_h_ii == b_primes2
    }
}

#[cfg(test)]
mod tests {
    // use super::*;
    // use crate::prover::{LabradorProver, Witness};

    // #[test]
    // fn test_verify() {
    //     // set up example environment, use set1 for testing.
    //     let ep_1 = EnvironmentParameters::default();
    //     // generate a random witness based on ep above
    //     let witness_1 = Witness::new(&ep_1);
    //     // generate public statements based on witness_1
    //     let st: Statement = Statement::new(&witness_1, &ep_1);
    //     // generate the common reference string matriices
    //     let pp = PublicPrams::new(&ep_1);
    //     // generate random challenges
    //     let tr = Challenges::new(&ep_1);

    //     // create a new prover
    //     let prover = LabradorProver::new(&pp, &witness_1, &st, &tr);
    //     let proof = prover.prove(&ep_1);

    //     // create a new verifier
    //     let verifier = LabradorVerifier::new(&pp, &st, &tr);
    //     let result = verifier.verify(&proof, &ep_1);
    //     assert!(result)
    // }
}
