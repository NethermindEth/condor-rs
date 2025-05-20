use crate::core::env_params::EnvironmentParameters;
use crate::ring::rq::Rq;
use crate::ring::rq_matrix::RqMatrix;
use crate::ring::rq_vector::RqVector;
use crate::ring::zq::Zq;

/// This struct serves as aggregation of functions with constant value 0.
pub struct ZeroConstantFunctionsAggregation<'a> {
    ep: &'a EnvironmentParameters,
}

impl<'a> ZeroConstantFunctionsAggregation<'a> {
    pub fn new(parameters: &'a EnvironmentParameters) -> Self {
        Self { ep: parameters }
    }

    /// Calculate a_double_primes from a_prime, a_{i,j}^{''k} = \sum_{l=1}^{L}\psi_l^{k}a_{ij}^{'(l)}
    ///
    /// @param: vector_psi: \psi_l^k
    /// @param: a_prime: a_{ij}^{'(l)}, each a_{ij} is a ring element (PolyRing)
    ///
    /// @return: a_{ij}^{''(k)}, return a vector length k of matrix a_{ij}^{''}
    pub fn calculate_agg_a_double_prime(
        &mut self,
        vector_psi: &[Vec<Zq>],
        a_prime: &[RqMatrix],
    ) -> Vec<RqMatrix> {
        let a_double_prime: Vec<Vec<RqVector>> = (0..vector_psi.len())
            .map(|k| {
                (0..self.ep.multiplicity)
                    .map(|i| {
                        (0..self.ep.multiplicity)
                            .map(|j| {
                                // calculate a_{ij}^{'(l)} * \psi_l^k
                                (0..self.ep.constraint_l)
                                    .map(|l: usize| {
                                        &a_prime[l].get_elements()[i].get_elements()[j]
                                            * &vector_psi[k][l]
                                    })
                                    .fold(
                                        // sum over all l
                                        Rq::zero(),
                                        |acc, val| &acc + &val,
                                    )
                            })
                            .collect::<RqVector>()
                    })
                    .collect::<Vec<RqVector>>()
            })
            .collect();

        a_double_prime.into_iter().map(RqMatrix::new).collect()
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
    ) -> Vec<Vec<RqVector>> {
        let phi_double_prime: Vec<Vec<RqVector>> = (0..self.ep.kappa)
            .map(|k| {
                (0..self.ep.multiplicity)
                    .map(|i| {
                        // \sum_{l=1}^{L}\psi_l^{k}\phi_{i}^{'(l)}
                        let left_side = (0..self.ep.constraint_l)
                            .map(|l| {
                                phi_prime[l][i]
                                    .iter()
                                    .map(|phi| phi * &vector_psi[k][l])
                                    .collect::<RqVector>()
                            })
                            .fold(RqVector::new(vec![Rq::zero(); self.ep.rank]), |acc, val| {
                                acc.iter().zip(val.iter()).map(|(a, b)| a + b).collect()
                            });

                        // Calculate the right side: \sum(\omega_j^{k} * \sigma_{-1} * pi_i^{j})
                        // Because the length of pi is n*d, so we need to split it into n parts, each part has d elements to do the conjugate automorphism.
                        let right_side = (0..(2 * self.ep.security_parameter))
                            .map(|j| {
                                let omega_j = vector_omega[k][j];
                                (0..self.ep.rank)
                                    .map(|chunk_index| {
                                        let start = chunk_index * Rq::DEGREE;
                                        let end = start + Rq::DEGREE;

                                        let pi_poly =
                                            Rq::new(pi[i][j][start..end].try_into().unwrap());
                                        let pi_poly_conjugate = pi_poly.conjugate_automorphism();
                                        &pi_poly_conjugate * &omega_j
                                    })
                                    .collect::<RqVector>()
                            })
                            .fold(RqVector::new(vec![Rq::zero(); self.ep.rank]), |acc, val| {
                                acc.iter().zip(val.iter()).map(|(a, b)| a + b).collect()
                            });

                        &left_side + &right_side
                    })
                    .collect::<Vec<RqVector>>()
            })
            .collect::<Vec<Vec<RqVector>>>();

        phi_double_prime
    }

    /// calculate b^{''(k)} = \sum_{i,j=1}^{r} a_{ij}^{''(k)} * <s_i, s_j> + \sum_{i=1}^{r} <\phi_{i}^{''(k)} * s_i>
    ///
    /// @param: a_ct_aggr: a_{ij}^{''(k)}
    /// @param: phi_ct_aggr: \phi_{i}^{''(k)}
    /// @param: witness: s_i
    ///
    /// @return: b^{''(k)}
    pub fn calculate_agg_b_double_prime(
        &mut self,
        a_double_prime: &[RqMatrix],
        phi_ct_aggr: &[Vec<RqVector>],
        witness: &[RqVector],
    ) -> RqVector {
        (0..self.ep.kappa)
            .map(|k| {
                (0..self.ep.multiplicity)
                    .map(|i| {
                        &(0..self.ep.multiplicity).map(|j| {
                    // calculate a_{ij}^{''(k)} * <s_i, s_j>
                    &a_double_prime[k].get_elements()[i].get_elements()[j]
                        * &witness[i].inner_product_poly_vector(&witness[j])
                })
                .fold(
                    // sum over all i,j
                    Rq::zero(),
                    |acc, val| &acc + &val,
                )
                // add \phi_{i}^{''(k)} * s[i]
                + &phi_ct_aggr[k][i].inner_product_poly_vector(&witness[i])
                    }) // sum over all i,j
                    .fold(Rq::zero(), |acc, val| &acc + &val)
            })
            .collect()
    }
}

pub struct FunctionsAggregation<'a> {
    ep: &'a EnvironmentParameters,
}

impl<'a> FunctionsAggregation<'a> {
    pub fn new(parameters: &'a EnvironmentParameters) -> Self {
        Self { ep: parameters }
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
    pub fn calculate_aggr_a(
        &self,
        a_constraint: &[RqMatrix],
        a_double_prime: &[RqMatrix],
        vector_alpha: &RqVector,
        vector_beta: &RqVector,
    ) -> Vec<RqVector> {
        let a_i: Vec<RqVector> = (0..self.ep.multiplicity)
            .map(|i| {
                (0..self.ep.multiplicity)
                    .map(|j| {
                        // calculate \sum(alpha_k * a_{ij}), k is constraint_k
                        let left_side = (0..self.ep.constraint_k)
                            .map(|k| {
                                &a_constraint[k].get_elements()[i].get_elements()[j]
                                    * &vector_alpha.get_elements()[k]
                            })
                            .fold(Rq::zero(), |acc, val| &acc + &val);

                        // calculate \sum(beta_k * a_{ij}^{''(k)}), k is size k
                        let right_side = (0..self.ep.kappa)
                            .map(|k| {
                                &a_double_prime[k].get_elements()[i].get_elements()[j]
                                    * &vector_beta.get_elements()[k]
                            })
                            .fold(Rq::zero(), |acc, val| &acc + &val);

                        &left_side + &right_side
                    })
                    .collect::<RqVector>()
            })
            .collect::<Vec<RqVector>>();

        a_i
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
        &self,
        phi_constraint: &[Vec<RqVector>],
        phi_ct_aggr: &[Vec<RqVector>],
        random_alpha: &RqVector,
        random_beta: &RqVector,
    ) -> Vec<RqVector> {
        let phi_i: Vec<RqVector> = (0..self.ep.multiplicity)
            .map(|i| {
                // calculate \sum(alpha_k * \phi_{i}^{k})
                let left_side: RqVector = (0..self.ep.constraint_k)
                    .map(|k| {
                        phi_constraint[k][i]
                            .iter()
                            .map(|phi| phi * &random_alpha.get_elements()[k])
                            .collect::<RqVector>()
                    })
                    .fold(RqVector::new(vec![Rq::zero(); self.ep.rank]), |acc, val| {
                        acc.iter().zip(val.iter()).map(|(a, b)| a + b).collect()
                    });

                // calculate \sum(beta_k * \phi_{i}^{''(k)})
                let right_side: RqVector = (0..self.ep.kappa)
                    .map(|k| {
                        phi_ct_aggr[k][i]
                            .iter()
                            .map(|phi| phi * &random_beta.get_elements()[k])
                            .collect::<RqVector>()
                    })
                    .fold(RqVector::new(vec![Rq::zero(); self.ep.rank]), |acc, val| {
                        acc.iter().zip(val.iter()).map(|(a, b)| a + b).collect()
                    });

                &left_side + &right_side
            })
            .collect::<Vec<RqVector>>();

        phi_i
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
        &self,
        b_constraint: &RqVector,
        b_ct_aggr: &RqVector,
        random_alpha: &RqVector,
        random_beta: &RqVector,
    ) -> Rq {
        let left_side = (0..self.ep.constraint_k)
            .map(|k| &b_constraint.get_elements()[k] * &random_alpha.get_elements()[k])
            .fold(Rq::zero(), |acc, val| &acc + &val);

        let right_side = (0..self.ep.kappa)
            .map(|k| &b_ct_aggr.get_elements()[k] * &random_beta.get_elements()[k])
            .fold(Rq::zero(), |acc, val| &acc + &val);

        &left_side + &right_side
    }
}

/// calculate z = c_1*s_1 + ... + c_r*s_r
/// or calculate Az = c_1*t_1 + ... + c_r*t_r
///
/// @param: x: witness s_i or Ajtai commitments t_i
/// @param: random_c: c_i from challenge set
///
/// return z
pub fn calculate_z(x: &[RqVector], random_c: &RqVector) -> RqVector {
    x.iter()
        .zip(random_c.iter())
        .map(|(s_row, c_element)| s_row * c_element)
        .fold(
            RqVector::new(vec![Rq::zero(); x[0].get_elements().len()]),
            |acc, x| &acc + &x,
        )
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
