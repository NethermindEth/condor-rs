use std::borrow::Borrow;
use std::ops::{Add, Mul};

use crate::core::env_params::EnvironmentParameters;
use crate::ring::rq::Rq;
use crate::ring::rq_matrix::RqMatrix;
use crate::ring::rq_vector::RqVector;
use crate::ring::zq::Zq;

pub fn compute_linear_combination<E, B, C>(elements: &[E], challenges: &[C]) -> B
where
    E: Borrow<B>,
    for<'a> &'a B: Mul<&'a C, Output = B>,
    for<'a> &'a B: Add<&'a B, Output = B>,
{
    debug_assert_eq!(
        elements.len(),
        challenges.len(),
        "vectors must be the same length"
    );
    debug_assert!(!elements.is_empty(), "`elements` must not be empty");

    let mut zipped_iter = elements.iter().zip(challenges.iter());
    // Must do the following as the init value in fold requires size of B
    let (e0, c0) = zipped_iter.next().unwrap();
    let init = e0.borrow() * c0;

    zipped_iter.fold(init, |acc, (elem, c)| &acc + &(elem.borrow() * c))
}

/// This struct serves as aggregation of functions with constant value 0.
pub struct ZeroConstantFunctionsAggregation<'a> {
    ep: &'a EnvironmentParameters,
    a_double_prime: Vec<RqMatrix>,
    phi_double_prime: Vec<Vec<RqVector>>,
}

impl<'a> ZeroConstantFunctionsAggregation<'a> {
    pub fn new(parameters: &'a EnvironmentParameters, k_range: usize) -> Self {
        Self {
            ep: parameters,
            a_double_prime: vec![RqMatrix::zero_symmetric(parameters.multiplicity); k_range],
            phi_double_prime: vec![
                vec![RqVector::zero(parameters.rank); parameters.multiplicity];
                k_range
            ],
        }
    }

    /// Calculate a_double_primes from a_prime, a_{i,j}^{''k} = \sum_{l=1}^{L}\psi_l^{k}a_{ij}^{'(l)}
    ///
    /// @param: vector_psi: \psi_l^k
    /// @param: a_prime: a_{ij}^{'(l)}, each a_{ij} is a ring element (PolyRing)
    ///
    /// @return: a_{ij}^{''(k)}, return a vector length k of matrix a_{ij}^{''}
    pub fn calculate_agg_a_double_prime(&mut self, vector_psi: &[Vec<Zq>], a_prime: &[RqMatrix]) {
        let mut a_prime_l_vector: Vec<&Rq> = Vec::new();
        for i in 0..self.ep.multiplicity {
            for j in 0..i + 1 {
                a_prime_l_vector.clear(); // Re-use a_prime_l to prevent repetitive heap allocations
                a_prime_l_vector = a_prime
                    .iter()
                    .map(|matrix| matrix.get_cell_symmetric(i, j))
                    .collect();

                for (k, matrix) in self.a_double_prime.iter_mut().enumerate() {
                    matrix.set_sell(
                        i,
                        j,
                        compute_linear_combination(&a_prime_l_vector, &vector_psi[k]),
                    );
                }
            }
        }
    }

    /// calculate \phi_{i}^{''(k)} = \sum_{l=1}^{L}\psi_l^{k}\phi_{i}^{'(l)} + \sum(\omega_j^{k} * \sigma_{-1} * pi_i^{j})
    /// in the prover process, page 17 from the paper.
    ///
    /// @param: phi_ct: \phi_{i}^{'(l)}
    /// @param: pi: pi_i^{j}
    /// @param: random_psi: \psi_l^{k}
    /// @param: random_omega: \omega_j^{k}
    ///
    /// return: \phi_{i}^{''(k)}
    pub fn calculate_agg_phi_double_prime(
        &mut self,
        phi_prime: &[Vec<RqVector>],
        pi: &[Vec<Vec<Zq>>],
        vector_psi: &[Vec<Zq>],
        vector_omega: &[Vec<Zq>],
    ) {
        let mut phi_prime_l_vector: Vec<&RqVector> = Vec::new();
        for i in 0..self.ep.multiplicity {
            phi_prime_l_vector.clear();
            phi_prime_l_vector = phi_prime.iter().map(|elems| &elems[i]).collect();
            for (k, phi_k) in self.phi_double_prime.iter_mut().enumerate() {
                phi_k[i] = compute_linear_combination(&phi_prime_l_vector, &vector_psi[k]);
            }
        }

        let mut conjugated_pi_j: Vec<RqVector> = Vec::new();
        for (i, pi_i) in pi.iter().enumerate() {
            conjugated_pi_j.clear();
            // conjugated_pi_j =
            conjugated_pi_j = pi_i
                .iter()
                .map(|zq_vector| {
                    zq_vector
                        .chunks_exact(Rq::DEGREE)
                        .map(|coeffs| Rq::new(coeffs.try_into().unwrap()).conjugate_automorphism())
                        .collect::<RqVector>()
                })
                .collect();
            for (k, phi_k) in self.phi_double_prime.iter_mut().enumerate() {
                phi_k[i] =
                    &phi_k[i] + &compute_linear_combination(&conjugated_pi_j, &vector_omega[k]);
            }
        }
    }

    /// calculate b^{''(k)} = \sum_{i,j=1}^{r} a_{ij}^{''(k)} * <s_i, s_j> + \sum_{i=1}^{r} <\phi_{i}^{''(k)} * s_i>
    ///
    /// @param: a_ct_aggr: a_{ij}^{''(k)}
    /// @param: phi_ct_aggr: \phi_{i}^{''(k)}
    /// @param: witness: s_i
    ///
    /// @return: b^{''(k)}
    pub fn calculate_agg_b_double_prime(&mut self, witness: &[RqVector]) -> RqVector {
        (0..self.ep.kappa)
            .map(|k| {
                (0..self.ep.multiplicity)
                    .map(|i| {
                        &(0..self.ep.multiplicity).map(|j| {
                    // calculate a_{ij}^{''(k)} * <s_i, s_j>
                    self.a_double_prime[k].get_cell_symmetric(i, j)
                        * &witness[i].inner_product_poly_vector(&witness[j])
                })
                .fold(
                    // sum over all i,j
                    Rq::zero(),
                    |acc, val| &acc + &val,
                )
                // add \phi_{i}^{''(k)} * s[i]
                + &self.phi_double_prime[k][i].inner_product_poly_vector(&witness[i])
                    }) // sum over all i,j
                    .fold(Rq::zero(), |acc, val| &acc + &val)
            })
            .collect()
    }

    pub fn get_alpha_double_prime(&self) -> &[RqMatrix] {
        &self.a_double_prime
    }

    pub fn get_phi_double_prime(&self) -> &[Vec<RqVector>] {
        &self.phi_double_prime
    }
}

pub struct FunctionsAggregation<'a> {
    ep: &'a EnvironmentParameters,
    aggregated_a: RqMatrix,
    aggregated_phi: Vec<RqVector>,
    aggregated_b: Rq,
}

impl<'a> FunctionsAggregation<'a> {
    pub fn new(parameters: &'a EnvironmentParameters) -> Self {
        Self {
            ep: parameters,
            aggregated_a: RqMatrix::zero_symmetric(parameters.multiplicity),
            aggregated_phi: vec![RqVector::zero(parameters.rank); parameters.multiplicity],
            aggregated_b: Rq::zero(),
        }
    }

    /// calculate a_i = \sum(alpha_k * a_{ij}) + \sum(beta_k * a_{ij}^{''(k)})
    /// equation 5, in the verifier process, page 18 from the paper.
    ///
    /// @param: a_constraint: a_{ij}
    /// @param: a_ct_aggr: a_{ij}^{''(k)}
    /// @param: random_alpha: alpha_k
    /// @param: random_beta: beta_k
    /// @param: ep: struct SizeParams
    ///
    /// @return: a_i
    pub fn calculate_agg_a(
        &mut self,
        a_constraint: &[RqMatrix],
        a_double_prime: &[RqMatrix],
        vector_alpha: &RqVector,
        vector_beta: &RqVector,
    ) {
        let mut a_constraint_k: Vec<&Rq> = Vec::new();
        let mut a_double_prime_k: Vec<&Rq> = Vec::new();
        for i in 0..self.ep.multiplicity {
            for j in 0..i + 1 {
                a_constraint_k.clear();
                a_constraint_k = a_constraint
                    .iter()
                    .map(|matrix| matrix.get_cell_symmetric(i, j))
                    .collect();
                a_double_prime_k.clear();
                a_double_prime_k = a_double_prime
                    .iter()
                    .map(|matrix| matrix.get_cell_symmetric(i, j))
                    .collect();
                self.aggregated_a.set_sell(
                    i,
                    j,
                    &compute_linear_combination::<&Rq, Rq, Rq>(
                        &a_constraint_k,
                        vector_alpha.get_elements(),
                    ) + &compute_linear_combination(&a_double_prime_k, vector_beta.get_elements()),
                );
            }
        }
    }

    /// calculate phi_i = \sum(alpha_k * \phi_{i}^{k}) + \sum(beta_k * \phi_{i}^{''(k)})
    /// equation 6, in the verifier process, page 18 from the paper.
    ///
    /// param: phi_constraint: \phi_{i}^{k}
    /// param: phi_ct_aggr: \phi_{i}^{''(k)}
    /// param: random_alpha: alpha_k
    /// param: random_beta: beta_k
    /// param: ep: struct SizeParams
    ///
    /// return: phi_i
    pub fn calculate_aggr_phi(
        &mut self,
        phi_constraint: &[Vec<RqVector>],
        phi_double_prime: &[Vec<RqVector>],
        vector_alpha: &RqVector,
        vector_beta: &RqVector,
    ) {
        let mut phi_constraint_k: Vec<&RqVector> = Vec::new();
        let mut phi_double_prime_k: Vec<&RqVector> = Vec::new();

        for i in 0..self.ep.multiplicity {
            phi_constraint_k.clear();
            phi_constraint_k = phi_constraint.iter().map(|element| &element[i]).collect();
            phi_double_prime_k.clear();
            phi_double_prime_k = phi_double_prime.iter().map(|element| &element[i]).collect();
            self.aggregated_phi[i] =
                &compute_linear_combination::<&RqVector, RqVector, Rq>(
                    &phi_constraint_k,
                    vector_alpha.get_elements(),
                ) + &compute_linear_combination(&phi_double_prime_k, vector_beta.get_elements());
        }
    }

    /// calculate b_i = \sum(alpha_k * b^{k}) + \sum(beta_k * b^{''(k})
    /// equation 7, in the verifier process, page 18 from the paper.
    ///
    /// @param: b_constraint: b^{k}
    /// @param: b_ct_aggr: b^{''(k)}
    /// @param: random_alpha: alpha_k
    /// @param: random_beta: beta_k
    /// @param: ep: struct SizeParams
    ///
    /// @return: b_i
    pub fn calculate_aggr_b(
        &mut self,
        b_constraint: &RqVector,
        b_double_prime: &RqVector,
        vector_alpha: &RqVector,
        vector_beta: &RqVector,
    ) {
        self.aggregated_b =
            &compute_linear_combination(b_constraint.get_elements(), vector_alpha.get_elements())
                + &compute_linear_combination(
                    b_double_prime.get_elements(),
                    vector_beta.get_elements(),
                )
    }

    pub fn get_agg_a(&self) -> &RqMatrix {
        &self.aggregated_a
    }

    pub fn get_appr_phi(&self) -> &[RqVector] {
        &self.aggregated_phi
    }

    pub fn get_aggr_b(&self) -> &Rq {
        &self.aggregated_b
    }
}

// #[cfg(test)]
// mod tests {
//     use super::*;
//     use crate::prover::Challenges;
//     use crate::verifier::LabradorVerifier;
//     #[test]
//     fn test_check_relation_full() {
//         // set up example environment, use set1 for testing.
//         let ep = EnvironmentParameters::default();
//         // generate a random witness based on ep above
//         let witness_1 = Witness::new(&ep);
//         // generate public statements based on witness_1
//         let st = Statement::new(&witness_1, &ep);
//         // generate random challenges
//         let tr = Challenges::new(&ep);
//         // first aggregation
//         let aggr_1 = AggregationOne::new(&witness_1, &st, &ep, &tr);
//         // second aggregation
//         let aggr_2 = AggregationTwo::new(&aggr_1, &st, &ep, &tr);

//         // calculate garbage polynomial g_{ij} = <s_i, s_j>
//         let g = (0..ep.multiplicity)
//             .map(|i| {
//                 (0..ep.multiplicity)
//                     .map(|j| witness_1.s[i].inner_product_poly_vector(&witness_1.s[j]))
//                     .collect()
//             })
//             .collect();

//         // calculate h_{ii}
//         let h = (0..ep.multiplicity)
//             .map(|i| {
//                 (0..ep.multiplicity)
//                     .map(|j| {
//                         let inner_phii_sj =
//                             aggr_2.phi_i[i].inner_product_poly_vector(&witness_1.s[j]);
//                         let inner_phij_si =
//                             aggr_2.phi_i[j].inner_product_poly_vector(&witness_1.s[i]);
//                         &inner_phii_sj + &inner_phij_si
//                     })
//                     .collect()
//             })
//             .collect();

//         // check aggregation relation
//         let relation = LabradorVerifier::check_relation(&aggr_2.a_i, &aggr_2.b_i, &g, &h);

//         assert!(relation);
//     }
// }
