use crate::ring::poly::{PolyRing, PolyVector, ZqVector};
use crate::ring::zq::Zq;

pub struct SizeParams {
    size_r: usize,
    size_n: usize,
    deg_bound_d: usize,
    size_k: usize,
    constraint_l: usize,
    constraint_k: usize,
}

/// Calculate aprimes from aprime_l, a_{i,j}^{''k} = \sum_{l=1}^{L}\psi_l^{k}a_{ij}^{'(l)}
///
/// @param: random_psi: \psi_l^k
/// @param: a_ct: a_{ij}^{'(l)}, each a_{ij} is a ring element (PolyRing)
/// @param: sp: struct SizeParams
///
/// @return: a_{ij}^{'(k)}, return a vector length k of matrix a_{ij}
#[rustfmt::skip]
pub fn calculate_aggr_ct_a(
    random_psi: &[ZqVector],
    a_ct: &[Vec<PolyVector>],
    sp: &SizeParams,
) -> Vec<Vec<PolyVector>> {
    let aprimes: Vec<Vec<PolyVector>> = (0..sp.size_k).map(|k| {
        let psi_k = &random_psi[k];
        (0..sp.size_r).map(|i| {
            (0..sp.size_r).map(|j| {
                // calculate a_{ij}^{'(l)} * \psi_l^k
                (0..sp.constraint_l).map(|l| {
                    &a_ct[l][i].get_elements()[j]
                        * &psi_k.get_coeffs()[l]
                })
                .fold(
                    // sum over all l
                    PolyRing::new(vec![Zq::ZERO; sp.deg_bound_d]),
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
/// @param: sp: struct SizeParams
/// @param: security_level2: 256 in the paper
///
/// return: \phi_{i}^{''(k)}
#[rustfmt::skip]
pub fn calculate_aggr_ct_phi(
    phi_ct: &[Vec<PolyVector>],
    pi: &[Vec<ZqVector>],
    random_psi: &[ZqVector],
    random_omega: &[ZqVector],
    sp: &SizeParams,
    security_level2: usize,
) -> Vec<Vec<PolyVector>> {
    let phi_ct_aggr: Vec<Vec<PolyVector>> = (0..sp.size_k).map(|k| {
        (0..sp.size_r).map(|i| {
            // \sum_{l=1}^{L}\psi_l^{k}\phi_{i}^{'(l)}
            let left_side: PolyVector = (0..sp.constraint_l).map(|l| {
                phi_ct[l][i]
                    .iter()
                    .map(|phi| {
                        phi * &random_psi[k].get_coeffs()[l]
                    }).collect::<PolyVector>()
            })
            .fold(
                PolyVector::new(vec![PolyRing::new(vec![Zq::ZERO; sp.deg_bound_d]); sp.size_n]),
                |acc, val| acc.iter().zip(val.iter()).map(|(a, b)| a + b).collect(),
            );

            // Calculate the right side: \sum(\omega_j^{k} * \sigma_{-1} * pi_i^{j})
            // Because the length of pi is n*d, so we need to split it into n parts, each part has d elements to do the conjugate automorphism.
            let right_side: PolyVector = (0..security_level2).map(|j| {
                let omega_j = random_omega[k].get_coeffs()[j];

                let poly_vec: PolyVector = (0..sp.size_n).map(|chunk_index| {
                    let start = chunk_index * sp.deg_bound_d;
                    let end = start + sp.deg_bound_d;

                    let pi_poly = PolyRing::new(
                        pi[i][j].get_coeffs()[start..end].to_vec(),
                    );
                    let pi_poly_conjugate = pi_poly.conjugate_automorphism();
                    &pi_poly_conjugate * &omega_j
                }).collect::<PolyVector>();

                poly_vec
            })
            .fold(
                PolyVector::new(vec![PolyRing::new(vec![Zq::ZERO; sp.deg_bound_d]); sp.size_n]),
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
/// @param: sp: struct SizeParams
///
/// @return: b^{''(k)}
#[rustfmt::skip]
pub fn calculate_aggr_ct_b(
    a_ct_aggr: &[Vec<PolyVector>],
    phi_ct_aggr: &[Vec<PolyVector>],
    witness: &[PolyVector],
    sp: &SizeParams,
) -> PolyVector {
    let b_primes: PolyVector = (0..sp.size_k).map(|k| {
        (0..sp.size_r).map(|i| {
            &(0..sp.size_r).map(|j| {
                // calculate a_{ij}^{''(k)} * <s_i, s_j>
                &a_ct_aggr[k][i].get_elements()[j]
                    * &witness[i].inner_product_poly_vector(&witness[j])
            })
            .fold(
                // sum over all i,j
                PolyRing::new(vec![Zq::ZERO; sp.deg_bound_d]),
                |acc, val| &acc + &val,
            )
            // add \phi_{i}^{''(k)} * w[i]
            + &phi_ct_aggr[k][i].inner_product_poly_vector(&witness[i])
        }) // sum over all i,j
        .fold(PolyRing::new(vec![Zq::ZERO; sp.deg_bound_d]), |acc, val| {
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
/// @param: sp: struct SizeParams
///
/// @return: a_i
#[rustfmt::skip]
pub fn calculate_aggr_a(
    a_constraint: &[Vec<PolyVector>],
    a_ct_aggr: &[Vec<PolyVector>],
    random_alpha: &PolyVector,
    random_beta: &PolyVector,
    sp: &SizeParams,
) -> Vec<PolyVector> {
    let a_i: Vec<PolyVector> = (0..sp.size_r).map(|i| {
        (0..sp.size_r).map(|j| {
            // calculate \sum(alpha_k * a_{ij}), k is constraint_k
            let left_side = (0..sp.constraint_k).map(|k| {
                &a_constraint[k][i].get_elements()[j]
                    * &random_alpha.get_elements()[k]
            })
            .fold(PolyRing::zero_poly(), |acc, val| &acc + &val);

            // calculate \sum(beta_k * a_{ij}^{''(k)}), k is size k
            let right_side = (0..sp.size_k).map(|k| {
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
/// param: sp: struct SizeParams
///
/// return: phi_i
#[rustfmt::skip]
pub fn calculate_aggr_phi(
    phi_constraint: &[Vec<PolyVector>],
    phi_ct_aggr: &[Vec<PolyVector>],
    random_alpha: &PolyVector,
    random_beta: &PolyVector,
    sp: &SizeParams,
) -> Vec<PolyVector> {
    let phi_i: Vec<PolyVector> = (0..sp.size_r).map(|i| {
        // calculate \sum(alpha_k * \phi_{i}^{k})
        let left_side: PolyVector = (0..sp.constraint_k).map(|k| {
            phi_constraint[k][i]
                .iter()
                .map(|phi| phi * &random_alpha.get_elements()[k])
                .collect::<PolyVector>()
        })
        .fold(
            PolyVector::new(vec![PolyRing::new(vec![Zq::ZERO; sp.deg_bound_d]); sp.size_n]),
            |acc, val| acc.iter().zip(val.iter()).map(|(a, b)| a + b ).collect(),
        );

        // calculate \sum(beta_k * \phi_{i}^{''(k)})
        let right_side: PolyVector = (0..sp.size_k).map(|k| {
            phi_ct_aggr[k][i]
                .iter()
                .map(|phi| phi * &random_beta.get_elements()[k])
                .collect::<PolyVector>()
        })
        .fold(
            PolyVector::new(vec![PolyRing::new(vec![Zq::ZERO; sp.deg_bound_d]); sp.size_n]),
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
/// @param: sp: struct SizeParams
///
/// @return: b_i
pub fn calculate_aggr_b(
    b_constraint: PolyVector,
    b_ct_aggr: PolyVector,
    random_alpha: PolyVector,
    random_beta: PolyVector,
    sp: &SizeParams,
) -> PolyRing {
    let left_side = (0..sp.constraint_k)
        .map(|k| &b_constraint.get_elements()[k] * &random_alpha.get_elements()[k])
        .fold(PolyRing::zero_poly(), |acc, val| &acc + &val);

    let right_side = (0..sp.size_k)
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

/// calculate b^{k} = \sum(a_{ij}^{k}<s_i, s_j>) + \sum(<phi_{i}^{k}, s_i>), k \in [K]
/// in Prover initialization process, page 17 from the paper.
///
/// @param: s: s_i
/// @param: a_constraint: a_{ij}^{k}
/// @param: phi_constraint: \phi_{i}^{k}
///
/// @return: b^{k}
#[rustfmt::skip]
pub fn calculate_b_constraint(
    s: &[PolyVector],
    a_constraint: &[PolyVector],
    phi_constraint: &[PolyVector],
) -> PolyRing {
    let size_s = s.len();
    // calculate \sum(a_{ij}^{k}<s_i, s_j>)
    let left_side = (0..size_s).map(|i| {
        (0..size_s).map(|j| {
            &a_constraint[i].get_elements()[j]
                * &s[i].inner_product_poly_vector(&s[j])
        })
        .fold(PolyRing::zero_poly(), |acc, val| &acc + &val )
    })
    .fold(PolyRing::zero_poly(), |acc, val| &acc + &val );

    // calculate \sum(<phi_{i}^{k}, s_i>)
    let right_side = (0..size_s).fold(PolyRing::zero_poly(), |acc, i| {
        &acc + &phi_constraint[i].inner_product_poly_vector(&s[i])
    });

    &left_side + &right_side
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rng;
    use crate::core::challenge_set::ChallengeSet;
    use crate::core::jl::ProjectionMatrix;

    #[test]
    /// dummy test for checking the relation, random generated a_aggr and phi_aggr instead of using the corresponding calculation functions.
    fn test_check_relation_dummy() {
        // number of witnesses
        let size_r: usize = 3;
        let size_n: usize = 5;
        let deg_bound_d: usize = 2;
        // generate example witness_1 with size: r * s
        let witness_1: Vec<PolyVector> = (0..size_r)
            .map(|_| PolyVector::random(size_n, deg_bound_d))
            .collect();

        // generate random a_aggr with size: r * r
        let a_aggr: Vec<PolyVector> = (0..size_r)
            .map(|_| PolyVector::random(size_r, deg_bound_d))
            .collect();

        // generate random phi_aggr with size: r * n
        let phi_aggr: Vec<PolyVector> = (0..size_r)
            .map(|_| PolyVector::random(size_n, deg_bound_d))
            .collect();

        // calculate b_1, which is one of the b_aggr
        let b_1 = calculate_b_constraint(&witness_1, &a_aggr, &phi_aggr);

        // calculate garbage polynomial g_{ij} = <s_i, s_j>
        let g: Vec<PolyVector> = (0..size_r)
            .map(|i| {
                (0..size_r)
                    .map(|j| witness_1[i].inner_product_poly_vector(&witness_1[j]))
                    .collect()
            })
            .collect();

        // calculate h_{ii}
        let h: Vec<PolyVector> = (0..size_r)
            .map(|i| {
                (0..size_r)
                    .map(|j| {
                        let inner_phii_sj = phi_aggr[i].inner_product_poly_vector(&witness_1[j]);
                        let inner_phij_si = phi_aggr[j].inner_product_poly_vector(&witness_1[i]);
                        &inner_phii_sj + &inner_phij_si
                    })
                    .collect()
            })
            .collect();

        //check aggregation relation
        let relation = check_relation(&a_aggr, &b_1, &g, &h);
        assert!(relation);
    }

    #[test]
    fn test_check_relation_full() {
        // ---setup parameters starts---
        // number of witnesses
        let r: usize = 3;
        let n: usize = 5;
        let deg_bound_d: usize = 4;
        // set constant D for JL projection matrix, D = deg_bound_d
        const D: usize = 4;

        // security level is \lambda in the paper
        let security_level: usize = 128;
        let security_level2: usize = 2 * security_level;

        // constraint_l is L
        let constraint_l: usize = 5;
        // constraint_k is K
        let constraint_k: usize = 5;

        // aggregated to only 128/log_q functions, q is 2^32 in the paper
        let log_q: usize = 32;
        let k = security_level / log_q;
        // ---setup parameters ends---

        // size_params contains size_r, size_n, deg_bound_d, size_k, constraint_l to reduce the size of the input parameters.
        let size_params = SizeParams {
            size_r: r,
            size_n: n,
            deg_bound_d,
            size_k: k,
            constraint_l,
            constraint_k,
        };

        // generate random witness_1 with size: r * s
        let witness_1: Vec<PolyVector> =
            (0..r).map(|_| PolyVector::random(n, deg_bound_d)).collect();

        // generate random a_constraint with size: constraint_k * r * n
        let a_constraint: Vec<Vec<PolyVector>> = (0..constraint_k)
            .map(|_| (0..r).map(|_| PolyVector::random(n, deg_bound_d)).collect())
            .collect();

        // generate random phi_constraint with size: constraint_k * r * n
        let phi_constraint: Vec<Vec<PolyVector>> = (0..constraint_k)
            .map(|_| (0..r).map(|_| PolyVector::random(n, deg_bound_d)).collect())
            .collect();

        // generate random b_constraint with size: constraint_k
        let b_constraint: PolyVector = (0..constraint_k)
            .map(|k| calculate_b_constraint(&witness_1, &a_constraint[k], &phi_constraint[k]))
            .collect();

        // generate example a_ct with size: constraint_l * r * n
        let a_ct: Vec<Vec<PolyVector>> = (0..constraint_l)
            .map(|_| (0..r).map(|_| PolyVector::random(n, deg_bound_d)).collect())
            .collect();

        // generate random phi_ct with size: constraint_k * r * n
        // it is a k length vector of matrix with size: r * n
        let phi_ct: Vec<Vec<PolyVector>> = (0..constraint_k)
            .map(|_| (0..r).map(|_| PolyVector::random(n, deg_bound_d)).collect())
            .collect();

        // generate random psi with size: k * constraint_l, each element is Zq
        let random_psi: Vec<ZqVector> = (0..k)
            .map(|_| ZqVector::random(&mut rng(), constraint_l))
            .collect();

        // generate randm omega is with size: k * lambda2, each element is Zq
        let random_omega: Vec<ZqVector> = (0..k)
            .map(|_| ZqVector::random(&mut rng(), security_level2))
            .collect();

        // calculate a_{ij}^{''(k)}
        let a_ct_aggr: Vec<Vec<PolyVector>> = calculate_aggr_ct_a(&random_psi, &a_ct, &size_params);

        // calculate \phi_{i}^{''(k)}
        // \pi is from JL projection, pi contains r matrices and each matrix: security_level2 * (n*d), (security_level2 is 256 in the paper).
        // explicitly set the deg_bound_d to D == deg_bound_d
        let pi: Vec<Vec<ZqVector>> = (0..r)
            .map(|_| ProjectionMatrix::<D>::new(n).get_matrix().clone())
            .collect();

        // calculate \phi_{i}^{''(k)}
        // only this function needs "security_level2", so I didn't add it to size_params.
        let phi_ct_aggr: Vec<Vec<PolyVector>> = calculate_aggr_ct_phi(
            &phi_ct,
            &pi,
            &random_psi,
            &random_omega,
            &size_params,
            security_level2,
        );

        // calculate b^{''(k)}
        let b_ct_aggr: PolyVector =
            calculate_aggr_ct_b(&a_ct_aggr, &phi_ct_aggr, &witness_1, &size_params);

        // generate random alpha and beta from challenge set
        let cs_alpha: ChallengeSet = ChallengeSet::new(deg_bound_d);
        let random_alpha: PolyVector = (0..constraint_k)
            .map(|_| cs_alpha.get_challenges().clone())
            .collect();

        let cs_beta: ChallengeSet = ChallengeSet::new(deg_bound_d);
        let random_beta: PolyVector = (0..constraint_k)
            .map(|_| cs_beta.get_challenges().clone())
            .collect();

        // calculate a_i
        let a_i: Vec<PolyVector> = calculate_aggr_a(
            &a_constraint,
            &a_ct_aggr,
            &random_alpha,
            &random_beta,
            &size_params,
        );

        // calculate phi_i
        let phi_i: Vec<PolyVector> = calculate_aggr_phi(
            &phi_constraint,
            &phi_ct_aggr,
            &random_alpha,
            &random_beta,
            &size_params,
        );

        // calculate b_i
        let b_i: PolyRing = calculate_aggr_b(
            b_constraint,
            b_ct_aggr,
            random_alpha,
            random_beta,
            &size_params,
        );

        // calculate garbage polynomial g_{ij} = <s_i, s_j>
        let g: Vec<PolyVector> = (0..r)
            .map(|i| {
                (0..r)
                    .map(|j| witness_1[i].inner_product_poly_vector(&witness_1[j]))
                    .collect()
            })
            .collect();

        // calculate h_{ii}
        let h: Vec<PolyVector> = (0..r)
            .map(|i| {
                (0..r)
                    .map(|j| {
                        let inner_phii_sj = phi_i[i].inner_product_poly_vector(&witness_1[j]);
                        let inner_phij_si = phi_i[j].inner_product_poly_vector(&witness_1[i]);
                        &inner_phii_sj + &inner_phij_si
                    })
                    .collect()
            })
            .collect();

        // check aggregation relation
        let relation = check_relation(&a_i, &b_i, &g, &h);
        assert!(relation);
    }
}
