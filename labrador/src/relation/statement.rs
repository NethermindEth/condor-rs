use super::env_params::EnvironmentParameters;
use crate::relation::witness::Witness;
use crate::ring::rq::Rq;
use crate::ring::rq_matrix::RqMatrix;
use crate::ring::rq_vector::RqVector;
use crate::ring::zq::ZqLabrador;
type Zq = ZqLabrador;

use crate::core::inner_product;

/// Statement is the input of the prover, which contains the constraints.
/// All parameters are from line 1, st, in the verifier process, page 18 from the paper.
pub struct Statement {
    // $a_{ij}^k$
    pub a_constraint: Vec<RqMatrix>,
    // $\varphi_i^k$
    pub phi_constraint: Vec<Vec<RqVector>>,
    // $b^{(k)}$
    pub b_constraint: RqVector,
    // $a_{ij}^{'(l)}$
    pub a_ct: Vec<RqMatrix>,
    // $\varphi_i^{'(l)}$
    pub phi_ct: Vec<Vec<RqVector>>,
    // $b_0^{'(l)}$
    pub b_0_ct: Vec<Zq>,
}

impl Statement {
    pub fn new(witness: &Witness, ep: &EnvironmentParameters) -> Self {
        // generate random a_constraint with size: constraint_k * r * n
        let a_constraint: Vec<RqMatrix> = (0..ep.constraint_k)
            .map(|_| RqMatrix::symmetric_random(&mut rand::rng(), ep.multiplicity))
            .collect();

        // generate random phi_constraint with size: constraint_k * r * n
        let phi_constraint: Vec<Vec<RqVector>> = (0..ep.constraint_k)
            .map(|_| {
                (0..ep.multiplicity)
                    .map(|_| RqVector::random(&mut rand::rng(), ep.rank))
                    .collect()
            })
            .collect();

        // calculate b_constraint b^k with size: constraint_k
        let b_constraint: RqVector = (0..ep.constraint_k)
            .map(|k| calculate_b_constraint(&witness.s, &a_constraint[k], &phi_constraint[k]))
            .collect();

        // generate example a_ct with size: constraint_l * r * n
        let a_ct: Vec<RqMatrix> = (0..ep.constraint_l)
            .map(|_| RqMatrix::symmetric_random(&mut rand::rng(), ep.multiplicity))
            .collect();

        // generate random phi_ct with size: constraint_k * r * n
        // it is a k length vector of matrix with size: r * n
        let phi_ct: Vec<Vec<RqVector>> = (0..ep.constraint_l)
            .map(|_| {
                (0..ep.multiplicity)
                    .map(|_| RqVector::random(&mut rand::rng(), ep.rank))
                    .collect()
            })
            .collect();

        // calculate b^l with size: constraint_l
        let b_constraint_l: RqVector = (0..ep.constraint_l)
            .map(|l| calculate_b_constraint(&witness.s, &a_ct[l], &phi_ct[l]))
            .collect();

        // calculate b_0^l
        let b_0_ct: Vec<Zq> = (0..ep.constraint_l)
            .map(|l| b_constraint_l.elements()[l].coeffs()[0])
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
    s: &[RqVector],
    a_constraint: &RqMatrix,
    phi_constraint: &[RqVector],
) -> Rq {
    let size_s = s.len();
    // calculate \sum(a_{ij}^{k}<s_i, s_j>)
    let left_side = (0..size_s).map(|i| {
        (0..size_s).map(|j: usize| {
            a_constraint.get_cell(i, j)
                * &inner_product::compute_linear_combination(s[i].elements(), s[j].elements())
        })
        .fold(Rq::zero(), |acc, val| &acc + &val )
    })
    .fold(Rq::zero(), |acc, val| &acc + &val );

    // calculate \sum(<phi_{i}^{k}, s_i>)
    let right_side = (0..size_s).fold(Rq::zero(), |acc, i| {
        &acc + &inner_product::compute_linear_combination(phi_constraint[i].elements(), s[i].elements())
    });

    &left_side + &right_side
}
