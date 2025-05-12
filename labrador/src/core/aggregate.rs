use crate::core::{env_params::EnvironmentParameters, statement::Statement};
use crate::prover::Witness;
use crate::ring::rq::Rq;
use crate::ring::rq_vector::RqVector;
use crate::ring::zq::Zq;

/// First step of aggregation
pub struct AggregationOne {
    // a_{ij}^{''(k)}
    pub a_ct_aggr: Vec<Vec<RqVector>>,
    // \phi_{i}^{''(k)}
    pub phi_ct_aggr: Vec<Vec<RqVector>>,
    // b^{''(k)}
    pub b_ct_aggr: RqVector,
}

impl AggregationOne {
    pub fn new(
        witness: &Witness,
        st: &Statement,
        ep: &EnvironmentParameters,
        pi: &[Vec<Vec<Zq>>],
        psi: &[Vec<Zq>],
        omega: &[Vec<Zq>],
    ) -> Self {
        Self::aggregate(witness, st, ep, pi, psi, omega)
    }

    fn aggregate(
        witness: &Witness,
        st: &Statement,
        ep: &EnvironmentParameters,
        pi: &[Vec<Vec<Zq>>],
        psi: &[Vec<Zq>],
        omega: &[Vec<Zq>],
    ) -> Self {
        // calculate a_{ij}^{''(k)}
        let a_ct_aggr: Vec<Vec<RqVector>> = Self::get_a_ct_aggr(psi, &st.a_ct, ep);

        // calculate \phi_{i}^{''(k)}
        let phi_ct_aggr: Vec<Vec<RqVector>> = Self::get_phi_ct_aggr(&st.phi_ct, pi, psi, omega, ep);

        // calculate b^{''(k)}
        let b_ct_aggr: RqVector = Self::get_b_ct_aggr(&a_ct_aggr, &phi_ct_aggr, witness, ep);

        Self {
            a_ct_aggr,
            phi_ct_aggr,
            b_ct_aggr,
        }
    }

    pub fn get_a_ct_aggr(
        psi: &[Vec<Zq>],
        a_ct: &[Vec<RqVector>],
        ep: &EnvironmentParameters,
    ) -> Vec<Vec<RqVector>> {
        calculate_aggr_ct_a(psi, a_ct, ep)
    }

    pub fn get_phi_ct_aggr(
        phi_ct: &[Vec<RqVector>],
        pi: &[Vec<Vec<Zq>>],
        psi: &[Vec<Zq>],
        omega: &[Vec<Zq>],
        ep: &EnvironmentParameters,
    ) -> Vec<Vec<RqVector>> {
        calculate_aggr_ct_phi(phi_ct, pi, psi, omega, ep)
    }

    pub fn get_b_ct_aggr(
        a_ct_aggr: &[Vec<RqVector>],
        phi_ct_aggr: &[Vec<RqVector>],
        witness: &Witness,
        ep: &EnvironmentParameters,
    ) -> RqVector {
        calculate_aggr_ct_b(a_ct_aggr, phi_ct_aggr, &witness.s, ep)
    }
}

/// Second step of aggregation
pub struct AggregationTwo {
    // a_{ij}
    pub a_i: Vec<RqVector>,
    // \phi_{i}
    pub phi_i: Vec<RqVector>,
    // b
    pub b_i: Rq,
}

impl AggregationTwo {
    pub fn new(
        aggr_one: &AggregationOne,
        st: &Statement,
        ep: &EnvironmentParameters,
        alpha_vector: &RqVector,
        beta_vector: &RqVector,
    ) -> Self {
        Self::aggregate(aggr_one, st, ep, alpha_vector, beta_vector)
    }

    fn aggregate(
        aggr_one: &AggregationOne,
        st: &Statement,
        ep: &EnvironmentParameters,
        alpha_vector: &RqVector,
        beta_vector: &RqVector,
    ) -> Self {
        // calculate a_i
        let a_i = Self::get_a_i(
            &st.a_constraint,
            &aggr_one.a_ct_aggr,
            alpha_vector,
            beta_vector,
            ep,
        );

        // calculate phi_i
        let phi_i = Self::get_phi_i(
            &st.phi_constraint,
            &aggr_one.phi_ct_aggr,
            alpha_vector,
            beta_vector,
            ep,
        );

        // calculate b_i
        let b_i = Self::get_b_i(
            &st.b_constraint,
            &aggr_one.b_ct_aggr,
            alpha_vector,
            beta_vector,
            ep,
        );

        Self { a_i, phi_i, b_i }
    }

    pub fn get_a_i(
        a_constraint: &[Vec<RqVector>],
        a_ct_aggr: &[Vec<RqVector>],
        random_alpha: &RqVector,
        random_beta: &RqVector,
        ep: &EnvironmentParameters,
    ) -> Vec<RqVector> {
        calculate_aggr_a(a_constraint, a_ct_aggr, random_alpha, random_beta, ep)
    }

    pub fn get_phi_i(
        phi_constraint: &[Vec<RqVector>],
        phi_ct_aggr: &[Vec<RqVector>],
        random_alpha: &RqVector,
        random_beta: &RqVector,
        ep: &EnvironmentParameters,
    ) -> Vec<RqVector> {
        calculate_aggr_phi(phi_constraint, phi_ct_aggr, random_alpha, random_beta, ep)
    }

    pub fn get_b_i(
        b_constraint: &RqVector,
        b_ct_aggr: &RqVector,
        random_alpha: &RqVector,
        random_beta: &RqVector,
        ep: &EnvironmentParameters,
    ) -> Rq {
        calculate_aggr_b(b_constraint, b_ct_aggr, random_alpha, random_beta, ep)
    }
}

/// Calculate aprimes from aprime_l, a_{i,j}^{''k} = \sum_{l=1}^{L}\psi_l^{k}a_{ij}^{'(l)}
///
/// @param: random_psi: \psi_l^k
/// @param: a_ct: a_{ij}^{'(l)}, each a_{ij} is a ring element (PolyRing)
/// @param: ep: struct SizeParams
///
/// @return: a_{ij}^{'(k)}, return a vector length k of matrix a_{ij}
#[rustfmt::skip]
fn calculate_aggr_ct_a(
    random_psi: &[Vec<Zq>],
    a_ct: &[Vec<RqVector>],
    ep: &EnvironmentParameters,
) -> Vec<Vec<RqVector>> {
    let aprimes: Vec<Vec<RqVector>> = (0..ep.kappa).map(|k| {
        let psi_k = &random_psi[k];
        (0..ep.r).map(|i| {
            (0..ep.r).map(|j| {
                // calculate a_{ij}^{'(l)} * \psi_l^k
                (0..ep.constraint_l).map(|l| {
                    &a_ct[l][i].get_elements()[j]
                        * &psi_k[l]
                })
                .fold(
                    // sum over all l
                    Rq::zero(),
                    |acc, val| acc + val,
                )
            }).collect::<RqVector>()
        }).collect::<Vec<RqVector>>()
    }).collect();

    aprimes
}

/// calculate \phi_{i}^{''(k)} = \sum_{l=1}^{L}\psi_l^{k}\phi_{i}^{'(l)} + \sum(\omega_j^{k} * \sigma_{-1} * pi_i^{j})
/// in the prover process, page 17 from the paper.
///
/// @param: phi_ct: \phi_{i}^{'(l)}
/// @param: pi: pi_i^{j}
/// @param: random_psi: \psi_l^{k}
/// @param: random_omega: \omega_j^{k}
/// @param: ep: struct SizeParams
/// @param: security_level2: 256 in the paper
///
/// return: \phi_{i}^{''(k)}
fn calculate_aggr_ct_phi(
    phi_ct: &[Vec<RqVector>],
    pi: &[Vec<Vec<Zq>>],
    random_psi: &[Vec<Zq>],
    random_omega: &[Vec<Zq>],
    ep: &EnvironmentParameters,
) -> Vec<Vec<RqVector>> {
    let phi_ct_aggr: Vec<Vec<RqVector>> = (0..ep.kappa)
        .map(|k| {
            (0..ep.r)
                .map(|i| {
                    // \sum_{l=1}^{L}\psi_l^{k}\phi_{i}^{'(l)}
                    let left_side = (0..ep.constraint_l)
                        .map(|l| {
                            phi_ct[l][i]
                                .iter()
                                .map(|phi| phi * &random_psi[k][l])
                                .collect::<RqVector>()
                        })
                        .fold(RqVector::new(vec![Rq::zero(); ep.n]), |acc, val| {
                            acc.iter().zip(val.iter()).map(|(a, b)| a + b).collect()
                        });

                    // Calculate the right side: \sum(\omega_j^{k} * \sigma_{-1} * pi_i^{j})
                    // Because the length of pi is n*d, so we need to split it into n parts, each part has d elements to do the conjugate automorphism.
                    let right_side = (0..(2 * ep.lambda))
                        .map(|j| {
                            let omega_j = random_omega[k][j];
                            (0..ep.n)
                                .map(|chunk_index| {
                                    let start = chunk_index * Rq::DEGREE;
                                    let end = start + Rq::DEGREE;

                                    let pi_poly = Rq::new(pi[i][j][start..end].try_into().unwrap());
                                    let pi_poly_conjugate = pi_poly.conjugate_automorphism();
                                    &pi_poly_conjugate * &omega_j
                                })
                                .collect::<RqVector>()
                        })
                        .fold(RqVector::new(vec![Rq::zero(); ep.n]), |acc, val| {
                            acc.iter().zip(val.iter()).map(|(a, b)| a + b).collect()
                        });

                    &left_side + &right_side
                })
                .collect::<Vec<RqVector>>()
        })
        .collect::<Vec<Vec<RqVector>>>();

    phi_ct_aggr
}

/// calculate b^{''(k)} = \sum_{i,j=1}^{r} a_{ij}^{''(k)} * <s_i, s_j> + \sum_{i=1}^{r} <\phi_{i}^{''(k)} * s_i>
///
/// @param: a_ct_aggr: a_{ij}^{''(k)}
/// @param: phi_ct_aggr: \phi_{i}^{''(k)}
/// @param: witness: s_i
/// @param: ep: struct SizeParams
///
/// @return: b^{''(k)}
#[rustfmt::skip]
fn calculate_aggr_ct_b(
    a_ct_aggr: &[Vec<RqVector>],
    phi_ct_aggr: &[Vec<RqVector>],
    witness: &[RqVector],
    ep: &EnvironmentParameters,
) -> RqVector {
    (0..ep.kappa).map(|k| {
        (0..ep.r).map(|i| {
            (0..ep.r).map(|j| {
                // calculate a_{ij}^{''(k)} * <s_i, s_j>
                a_ct_aggr[k][i].get_elements()[j]
                    * witness[i].inner_product_poly_vector(&witness[j])
            })
            .fold(
                // sum over all i,j
                Rq::zero(),
                |acc, val| acc + val,
            )
            // add \phi_{i}^{''(k)} * s[i]
            + phi_ct_aggr[k][i].inner_product_poly_vector(&witness[i])
        }) // sum over all i,j
        .fold(Rq::zero(), |acc, val| {
            acc + val
        })
    }).collect()
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
fn calculate_aggr_a(
    a_constraint: &[Vec<RqVector>],
    a_ct_aggr: &[Vec<RqVector>],
    random_alpha: &RqVector,
    random_beta: &RqVector,
    ep: &EnvironmentParameters,
) -> Vec<RqVector> {
    let a_i: Vec<RqVector> = (0..ep.r)
        .map(|i| {
            (0..ep.r)
                .map(|j| {
                    // calculate \sum(alpha_k * a_{ij}), k is constraint_k
                    let left_side = (0..ep.constraint_k)
                        .map(|k| {
                            a_constraint[k][i].get_elements()[j] * random_alpha.get_elements()[k]
                        })
                        .fold(Rq::zero(), |acc, val| acc + val);

                    // calculate \sum(beta_k * a_{ij}^{''(k)}), k is size k
                    let right_side = (0..ep.kappa)
                        .map(|k| a_ct_aggr[k][i].get_elements()[j] * random_beta.get_elements()[k])
                        .fold(Rq::zero(), |acc, val| acc + val);

                    left_side + right_side
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
fn calculate_aggr_phi(
    phi_constraint: &[Vec<RqVector>],
    phi_ct_aggr: &[Vec<RqVector>],
    random_alpha: &RqVector,
    random_beta: &RqVector,
    ep: &EnvironmentParameters,
) -> Vec<RqVector> {
    let phi_i: Vec<RqVector> = (0..ep.r)
        .map(|i| {
            // calculate \sum(alpha_k * \phi_{i}^{k})
            let left_side: RqVector = (0..ep.constraint_k)
                .map(|k| {
                    phi_constraint[k][i]
                        .iter()
                        .map(|phi| phi * &random_alpha.get_elements()[k])
                        .collect::<RqVector>()
                })
                .fold(RqVector::new(vec![Rq::zero(); ep.n]), |acc, val| {
                    acc.iter().zip(val.iter()).map(|(a, b)| a + b).collect()
                });

            // calculate \sum(beta_k * \phi_{i}^{''(k)})
            let right_side: RqVector = (0..ep.kappa)
                .map(|k| {
                    phi_ct_aggr[k][i]
                        .iter()
                        .map(|phi| phi * &random_beta.get_elements()[k])
                        .collect::<RqVector>()
                })
                .fold(RqVector::new(vec![Rq::zero(); ep.n]), |acc, val| {
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
fn calculate_aggr_b(
    b_constraint: &RqVector,
    b_ct_aggr: &RqVector,
    random_alpha: &RqVector,
    random_beta: &RqVector,
    ep: &EnvironmentParameters,
) -> Rq {
    let left_side = (0..ep.constraint_k)
        .map(|k| b_constraint.get_elements()[k] * random_alpha.get_elements()[k])
        .fold(Rq::zero(), |acc, val| acc + val);

    let right_side = (0..ep.kappa)
        .map(|k| b_ct_aggr.get_elements()[k] * random_beta.get_elements()[k])
        .fold(Rq::zero(), |acc, val| acc + val);

    left_side + right_side
}

/// calculate h_{ij} = 1/2 * (<\phi_i, s_j> + <\phi_j, s_i>), then use base b to decompose the polynomial
///
/// @param: phi_i: phi_i
/// @param: s: witness s_i
/// @param: ep: struct SizeParams
///
/// return h_{ij}
pub fn calculate_hij(
    phi_i: &[RqVector],
    s: &[RqVector],
    ep: &EnvironmentParameters,
) -> Vec<RqVector> {
    (0..ep.r)
        .map(|i| {
            (0..ep.r)
                .map(|j| {
                    let left_side = &phi_i[i].inner_product_poly_vector(&s[j]);
                    let right_side = &phi_i[j].inner_product_poly_vector(&s[i]);
                    left_side + right_side
                })
                .collect()
        })
        .collect()
}

/// calculate garbage polynomial g = <s_i, s_j>
///
/// @param: phi_i: phi_i
/// @param: s: witness s_i
/// @param: r: size_r
///
/// return g_{ij}
pub fn calculate_gij(s: &[RqVector], r: usize) -> Vec<RqVector> {
    (0..r)
        .map(|i| {
            (0..r)
                .map(|j| s[i].inner_product_poly_vector(&s[j]))
                .collect()
        })
        .collect()
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

/// line 19, page 18: calculate u_1 = \sum(B_ik * t_i^(k)) + sum(C_ijk * g_ij^(k))
pub fn calculate_u_1(
    b: &[Vec<Vec<RqVector>>],
    c: &[Vec<Vec<Vec<RqVector>>>],
    t_i: &[Vec<RqVector>],
    g_ij: &[Vec<RqVector>],
    ep: &EnvironmentParameters,
) -> RqVector {
    let mut u_1 = vec![Rq::zero(); ep.kappa_1];
    // calculate left side
    for i in 0..ep.r {
        for k in 0..ep.t_1 {
            let b_ik_t_ik = &t_i[i][k] * &b[i][k];
            u_1 = u_1
                .iter()
                .zip(b_ik_t_ik.iter())
                .map(|(a, b)| a + b)
                .collect();
        }
    }
    // calculate right side
    for i in 0..ep.r {
        for j in i..ep.r {
            for k in 0..ep.t_2 {
                let c_ijk_g_ij = &g_ij[i][j] * &c[i][j][k];
                u_1 = u_1
                    .iter()
                    .zip(c_ijk_g_ij.iter())
                    .map(|(a, b)| a + b)
                    .collect();
            }
        }
    }

    RqVector::new(u_1)
}

/// line 20, page 18: calculate u_2 = \sum(D_ijk * h_ij^(k))
pub fn calculate_u_2(
    d: &[Vec<Vec<Vec<RqVector>>>],
    h_ij: &[Vec<RqVector>],
    ep: &EnvironmentParameters,
) -> RqVector {
    // Pre-collect the iterator over (i, j, k) into a vector.
    let flat_vec: Vec<(usize, usize, usize)> = (0..ep.r)
        .flat_map(|i| (i..ep.r).flat_map(move |j| (0..ep.t_2).map(move |k| (i, j, k))))
        .collect();

    flat_vec
        .into_iter()
        .fold(RqVector::new(vec![Rq::zero(); ep.n]), |acc, (i, j, k)| {
            let d_ijk_h_ij = &h_ij[i][j] * &d[i][j][k];
            acc.iter()
                .zip(d_ijk_h_ij.iter())
                .map(|(a, b)| a + b)
                .collect::<RqVector>()
        })
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
//         let g = (0..ep.r)
//             .map(|i| {
//                 (0..ep.r)
//                     .map(|j| witness_1.s[i].inner_product_poly_vector(&witness_1.s[j]))
//                     .collect()
//             })
//             .collect();

//         // calculate h_{ii}
//         let h = (0..ep.r)
//             .map(|i| {
//                 (0..ep.r)
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
