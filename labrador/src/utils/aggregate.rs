use crate::{
    poly::{PolyRing, PolyVector, ZqVector},
    prover::{Challenges, Witness},
    utils::{env_params::EnvironmentParameters, statement::Statement},
    zq::Zq,
};

/// First step of aggregation
pub struct AggregationOne {
    // a_{ij}^{''(k)}
    pub a_ct_aggr: Vec<Vec<PolyVector>>,
    // \phi_{i}^{''(k)}
    pub phi_ct_aggr: Vec<Vec<PolyVector>>,
    // b^{''(k)}
    pub b_ct_aggr: PolyVector,
}

impl AggregationOne {
    pub fn new(
        witness: &Witness,
        st: &Statement,
        ep: &EnvironmentParameters,
        tr: &Challenges,
    ) -> Self {
        Self::aggregate(witness, st, ep, tr)
    }

    fn aggregate(
        witness: &Witness,
        st: &Statement,
        ep: &EnvironmentParameters,
        tr: &Challenges,
    ) -> Self {
        // calculate a_{ij}^{''(k)}
        let a_ct_aggr: Vec<Vec<PolyVector>> = calculate_aggr_ct_a(&tr.psi, &st.a_ct, ep);

        // calculate \phi_{i}^{''(k)}
        let phi_ct_aggr: Vec<Vec<PolyVector>> =
            calculate_aggr_ct_phi(&st.phi_ct, &tr.pi, &tr.psi, &tr.omega, ep);

        // calculate b^{''(k)}
        let b_ct_aggr: PolyVector = calculate_aggr_ct_b(&a_ct_aggr, &phi_ct_aggr, &witness.s, ep);

        Self {
            a_ct_aggr,
            phi_ct_aggr,
            b_ct_aggr,
        }
    }
}

/// Second step of aggregation
pub struct AggregationTwo {
    // a_{ij}
    pub a_i: Vec<PolyVector>,
    // \phi_{i}
    pub phi_i: Vec<PolyVector>,
    // b
    pub b_i: PolyRing,
}

impl AggregationTwo {
    pub fn new(
        aggr_one: &AggregationOne,
        st: &Statement,
        ep: &EnvironmentParameters,
        tr: &Challenges,
    ) -> Self {
        Self::aggregate(aggr_one, st, ep, tr)
    }

    fn aggregate(
        aggr_one: &AggregationOne,
        st: &Statement,
        ep: &EnvironmentParameters,
        tr: &Challenges,
    ) -> Self {
        // calculate a_i
        let a_i = calculate_aggr_a(
            &st.a_constraint,
            &aggr_one.a_ct_aggr,
            &tr.random_alpha,
            &tr.random_beta,
            ep,
        );

        // calculate phi_i
        let phi_i = calculate_aggr_phi(
            &st.phi_constraint,
            &aggr_one.phi_ct_aggr,
            &tr.random_alpha,
            &tr.random_beta,
            ep,
        );

        // calculate b_i
        let b_i = calculate_aggr_b(
            &st.b_constraint,
            &aggr_one.b_ct_aggr,
            &tr.random_alpha,
            &tr.random_beta,
            ep,
        );

        Self { a_i, phi_i, b_i }
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
    random_psi: &[ZqVector],
    a_ct: &[Vec<PolyVector>],
    ep: &EnvironmentParameters,
) -> Vec<Vec<PolyVector>> {
    let aprimes: Vec<Vec<PolyVector>> = (0..ep.k).map(|k| {
        let psi_k = &random_psi[k];
        (0..ep.r).map(|i| {
            (0..ep.r).map(|j| {
                // calculate a_{ij}^{'(l)} * \psi_l^k
                (0..ep.constraint_l).map(|l| {
                    &a_ct[l][i].get_elements()[j]
                        * &psi_k.get_coeffs()[l]
                })
                .fold(
                    // sum over all l
                    PolyRing::new(vec![Zq::ZERO; ep.deg_bound_d]),
                    |acc, val| &acc + &val,
                )
            }).collect::<PolyVector>()
        }).collect::<Vec<PolyVector>>()
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
#[rustfmt::skip]
fn calculate_aggr_ct_phi(
    phi_ct: &[Vec<PolyVector>],
    pi: &[Vec<ZqVector>],
    random_psi: &[ZqVector],
    random_omega: &[ZqVector],
    ep: &EnvironmentParameters,
) -> Vec<Vec<PolyVector>> {
    let phi_ct_aggr: Vec<Vec<PolyVector>> = (0..ep.k).map(|k| {
        (0..ep.r).map(|i| {
            // \sum_{l=1}^{L}\psi_l^{k}\phi_{i}^{'(l)}
            let left_side: PolyVector = (0..ep.constraint_l).map(|l| {
                phi_ct[l][i]
                    .iter()
                    .map(|phi| {
                        phi * &random_psi[k].get_coeffs()[l]
                    }).collect::<PolyVector>()
            })
            .fold(
                PolyVector::new(vec![PolyRing::new(vec![Zq::ZERO; ep.deg_bound_d]); ep.n]),
                |acc, val| acc.iter().zip(val.iter()).map(|(a, b)| a + b).collect(),
            );

            // Calculate the right side: \sum(\omega_j^{k} * \sigma_{-1} * pi_i^{j})
            // Because the length of pi is n*d, so we need to split it into n parts, each part has d elements to do the conjugate automorphism.
            let right_side: PolyVector = (0..ep.lambda2).map(|j| {
                let omega_j = random_omega[k].get_coeffs()[j];

                let poly_vec: PolyVector = (0..ep.n).map(|chunk_index| {
                    let start = chunk_index * ep.deg_bound_d;
                    let end = start + ep.deg_bound_d;

                    let pi_poly = PolyRing::new(
                        pi[i][j].get_coeffs()[start..end].to_vec(),
                    );
                    let pi_poly_conjugate = pi_poly.conjugate_automorphism();
                    &pi_poly_conjugate * &omega_j
                }).collect::<PolyVector>();

                poly_vec
            })
            .fold(
                PolyVector::new(vec![PolyRing::new(vec![Zq::ZERO; ep.deg_bound_d]); ep.n]),
                |acc, val| acc.iter().zip(val.iter()).map(|(a, b)| a + b).collect(),
            );

            &left_side + &right_side
        }).collect::<Vec<PolyVector>>()
    }).collect::<Vec<Vec<PolyVector>>>();

    phi_ct_aggr
}

/// calculate b^{''(k)} = \sum_{i,j=1}^{r} a_{ij}^{''(k)} * <s_i, s_j> + \sum_{i=1}^{r} \phi_{i}^{''(k)} * w[i]
///
/// @param: a_ct_aggr: a_{ij}^{''(k)}
/// @param: phi_ct_aggr: \phi_{i}^{''(k)}
/// @param: witness: s_i
/// @param: ep: struct SizeParams
///
/// @return: b^{''(k)}
#[rustfmt::skip]
fn calculate_aggr_ct_b(
    a_ct_aggr: &[Vec<PolyVector>],
    phi_ct_aggr: &[Vec<PolyVector>],
    witness: &[PolyVector],
    ep: &EnvironmentParameters,
) -> PolyVector {
    let b_primes: PolyVector = (0..ep.k).map(|k| {
        (0..ep.r).map(|i| {
            &(0..ep.r).map(|j| {
                // calculate a_{ij}^{''(k)} * <s_i, s_j>
                &a_ct_aggr[k][i].get_elements()[j]
                    * &witness[i].inner_product_poly_vector(&witness[j])
            })
            .fold(
                // sum over all i,j
                PolyRing::new(vec![Zq::ZERO; ep.deg_bound_d]),
                |acc, val| &acc + &val,
            )
            // add \phi_{i}^{''(k)} * w[i]
            + &phi_ct_aggr[k][i].inner_product_poly_vector(&witness[i])
        }) // sum over all i,j
        .fold(PolyRing::new(vec![Zq::ZERO; ep.deg_bound_d]), |acc, val| {
            &acc + &val
        })
    }).collect();

    b_primes
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
#[rustfmt::skip]
fn calculate_aggr_a(
    a_constraint: &[Vec<PolyVector>],
    a_ct_aggr: &[Vec<PolyVector>],
    random_alpha: &PolyVector,
    random_beta: &PolyVector,
    ep: &EnvironmentParameters,
) -> Vec<PolyVector> {
    let a_i: Vec<PolyVector> = (0..ep.r).map(|i| {
        (0..ep.r).map(|j| {
            // calculate \sum(alpha_k * a_{ij}), k is constraint_k
            let left_side = (0..ep.constraint_k).map(|k| {
                &a_constraint[k][i].get_elements()[j]
                    * &random_alpha.get_elements()[k]
            })
            .fold(PolyRing::zero_poly(), |acc, val| &acc + &val);

            // calculate \sum(beta_k * a_{ij}^{''(k)}), k is size k
            let right_side = (0..ep.k).map(|k| {
                &a_ct_aggr[k][i].get_elements()[j]
                    * &random_beta.get_elements()[k]
            })
            .fold(PolyRing::zero_poly(), |acc, val| &acc + &val);

            &left_side + &right_side
        }).collect::<PolyVector>()
    }).collect::<Vec<PolyVector>>();

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
#[rustfmt::skip]
fn calculate_aggr_phi(
    phi_constraint: &[Vec<PolyVector>],
    phi_ct_aggr: &[Vec<PolyVector>],
    random_alpha: &PolyVector,
    random_beta: &PolyVector,
    ep: &EnvironmentParameters,
) -> Vec<PolyVector> {
    let phi_i: Vec<PolyVector> = (0..ep.r).map(|i| {
        // calculate \sum(alpha_k * \phi_{i}^{k})
        let left_side: PolyVector = (0..ep.constraint_k).map(|k| {
            phi_constraint[k][i]
                .iter()
                .map(|phi| phi * &random_alpha.get_elements()[k])
                .collect::<PolyVector>()
        })
        .fold(
            PolyVector::new(vec![PolyRing::new(vec![Zq::ZERO; ep.deg_bound_d]); ep.n]),
            |acc, val| acc.iter().zip(val.iter()).map(|(a, b)| a + b ).collect(),
        );

        // calculate \sum(beta_k * \phi_{i}^{''(k)})
        let right_side: PolyVector = (0..ep.k).map(|k| {
            phi_ct_aggr[k][i]
                .iter()
                .map(|phi| phi * &random_beta.get_elements()[k])
                .collect::<PolyVector>()
        })
        .fold(
            PolyVector::new(vec![PolyRing::new(vec![Zq::ZERO; ep.deg_bound_d]); ep.n]),
            |acc, val| acc.iter().zip(val.iter()).map(|(a, b)| a + b ).collect(),
        );

        &left_side + &right_side
    }).collect::<Vec<PolyVector>>();

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
    b_constraint: &PolyVector,
    b_ct_aggr: &PolyVector,
    random_alpha: &PolyVector,
    random_beta: &PolyVector,
    ep: &EnvironmentParameters,
) -> PolyRing {
    let left_side = (0..ep.constraint_k)
        .map(|k| &b_constraint.get_elements()[k] * &random_alpha.get_elements()[k])
        .fold(PolyRing::zero_poly(), |acc, val| &acc + &val);

    let right_side = (0..ep.k)
        .map(|k| &b_ct_aggr.get_elements()[k] * &random_beta.get_elements()[k])
        .fold(PolyRing::zero_poly(), |acc, val| &acc + &val);

    &left_side + &right_side
}

/// equation 18: check if \sum(a_{ij} * g_{ij}) + \sum(h_{ii}) - b ?= 0
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

/// calculate h_{ij} = 1/2 * (<\phi_i, s_j> + <\phi_j, s_i>), then use base b to decompose the polynomial
///
/// @param: phi_i: phi_i
/// @param: s: witness s_i
/// @param: ep: struct SizeParams
///
/// return h_{ij}
pub fn calculate_hij(
    phi_i: &[PolyVector],
    s: &[PolyVector],
    ep: &EnvironmentParameters,
) -> Vec<PolyVector> {
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
    // Todo: decompose h: h{ij} = h_{ij}^0 + ... + h_{ij}^(t_1-1)*b_1^(t_1-1)
}

/// calculate z = c_1*s_1 + ... + c_r*s_r
///
/// @param: s: witness s_i
/// @param: random_c: c_i from challenge set
///
/// return z
pub fn calculate_z(s: &[PolyVector], random_c: &PolyVector) -> PolyVector {
    s.iter()
        .zip(random_c.iter())
        .map(|(s_row, c_element)| s_row * c_element)
        .fold(
            PolyVector::new(vec![PolyRing::zero_poly(); s[0].len()]),
            |acc, x| &acc + &x,
        )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prover::Challenges;

    #[test]
    fn test_check_relation_full() {
        // set up example environment, use set1 for testing.
        let ep = EnvironmentParameters::set_1();
        // generate a random witness based on ep above
        let witness_1 = Witness::new(&ep);
        // generate public statements based on witness_1
        let st: Statement = Statement::new(&witness_1, &ep);
        // generate random challenges
        let tr = Challenges::new(&ep);
        // first aggregation
        let aggr_1 = AggregationOne::new(&witness_1, &st, &ep, &tr);
        // second aggregation
        let aggr_2 = AggregationTwo::new(&aggr_1, &st, &ep, &tr);

        // calculate garbage polynomial g_{ij} = <s_i, s_j>
        let g: Vec<PolyVector> = (0..ep.r)
            .map(|i| {
                (0..ep.r)
                    .map(|j| witness_1.s[i].inner_product_poly_vector(&witness_1.s[j]))
                    .collect()
            })
            .collect();

        // calculate h_{ii}
        let h: Vec<PolyVector> = (0..ep.r)
            .map(|i| {
                (0..ep.r)
                    .map(|j| {
                        let inner_phii_sj =
                            aggr_2.phi_i[i].inner_product_poly_vector(&witness_1.s[j]);
                        let inner_phij_si =
                            aggr_2.phi_i[j].inner_product_poly_vector(&witness_1.s[i]);
                        &inner_phii_sj + &inner_phij_si
                    })
                    .collect()
            })
            .collect();

        // check aggregation relation
        let relation = check_relation(&aggr_2.a_i, &aggr_2.b_i, &g, &h);

        assert!(relation);

        // let z = calculate_z(&witness_1, &random_c);
    }
}
