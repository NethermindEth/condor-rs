use crate::core::env_params::EnvironmentParameters;
use crate::prover::Witness;
use crate::ring::poly::{PolyRing, PolyVector, ZqVector};

/// Statement is the input of the prover, which contains the constraints.
/// All parameters are from line 1, st, in the verifier process, page 18 from the paper.
pub struct Statement {
    // $a_{ij}^k$
    pub a_constraint: Vec<Vec<PolyVector>>,
    // $\varphi_i^k$
    pub phi_constraint: Vec<Vec<PolyVector>>,
    // $b^{(k)}$
    pub b_constraint: PolyVector,
    // $a_{ij}^{'(l)}$
    pub a_ct: Vec<Vec<PolyVector>>,
    // $\varphi_i^{'(l)}$
    pub phi_ct: Vec<Vec<PolyVector>>,
    // $b_0^{'(l)}$
    pub b_0_ct: ZqVector,
}

impl Statement {
    pub fn new(witness: &Witness, ep: &EnvironmentParameters) -> Self {
        // generate random a_constraint with size: constraint_k * r * n
        let a_constraint: Vec<Vec<PolyVector>> = (0..ep.constraint_k)
            .map(|_| {
                (0..ep.r)
                    .map(|_| PolyVector::random(ep.n, ep.deg_bound_d))
                    .collect()
            })
            .collect();

        // generate random phi_constraint with size: constraint_k * r * n
        let phi_constraint: Vec<Vec<PolyVector>> = (0..ep.constraint_k)
            .map(|_| {
                (0..ep.r)
                    .map(|_| PolyVector::random(ep.n, ep.deg_bound_d))
                    .collect()
            })
            .collect();

        // calculate b_constraint b^k with size: constraint_k
        let b_constraint: PolyVector = (0..ep.constraint_k)
            .map(|k| calculate_b_constraint(&witness.s, &a_constraint[k], &phi_constraint[k]))
            .collect();

        // generate example a_ct with size: constraint_l * r * n
        let a_ct: Vec<Vec<PolyVector>> = (0..ep.constraint_l)
            .map(|_| {
                (0..ep.r)
                    .map(|_| PolyVector::random(ep.n, ep.deg_bound_d))
                    .collect()
            })
            .collect();

        // generate random phi_ct with size: constraint_k * r * n
        // it is a k length vector of matrix with size: r * n
        let phi_ct: Vec<Vec<PolyVector>> = (0..ep.constraint_k)
            .map(|_| {
                (0..ep.r)
                    .map(|_| PolyVector::random(ep.n, ep.deg_bound_d))
                    .collect()
            })
            .collect();

        // calculate b^l with size: constraint_l
        let b_constraint_l: PolyVector = (0..ep.constraint_l)
            .map(|l| calculate_b_constraint(&witness.s, &a_constraint[l], &phi_constraint[l]))
            .collect();

        // calculate b_0^l
        let b_0_ct: ZqVector = (0..ep.constraint_l)
            .map(|l| b_constraint_l.get_elements()[l].get_coeffs()[0])
            .collect();

        Self {
            a_constraint,
            phi_constraint,
            b_constraint,
            a_ct,
            phi_ct,
            b_0_ct,
        }
    }
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
