use crate::relation::env_params::EnvironmentParameters;
use crate::ring::rq::Rq;
use crate::ring::rq_matrix::RqMatrix;
use crate::ring::rq_vector::RqVector;
use crate::ring::zq::Zq;

use super::inner_product;

/// This struct serves as aggregation of functions with constant value 0.
pub struct ZeroConstantFunctionsAggregation<'a> {
    ep: &'a EnvironmentParameters,
    a_double_prime: Vec<RqMatrix>,
    phi_double_prime: Vec<Vec<RqVector>>,
}

impl<'a> ZeroConstantFunctionsAggregation<'a> {
    pub fn new(parameters: &'a EnvironmentParameters) -> Self {
        Self {
            ep: parameters,
            a_double_prime: vec![
                RqMatrix::symmetric_zero(parameters.multiplicity);
                parameters.const_agg_length
            ],
            phi_double_prime: vec![
                vec![RqVector::zero(parameters.rank); parameters.multiplicity];
                parameters.const_agg_length
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
        for i in 0..self.ep.multiplicity {
            for j in 0..i + 1 {
                let a_prime_l_vector: Vec<&Rq> =
                    a_prime.iter().map(|matrix| matrix.get_cell(i, j)).collect();

                for (k, matrix) in self.a_double_prime.iter_mut().enumerate() {
                    matrix.set_cell(
                        i,
                        j,
                        inner_product::compute_linear_combination(
                            &a_prime_l_vector,
                            &vector_psi[k],
                        ),
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
        conjugated_pi: &[RqMatrix],
        vector_psi: &[Vec<Zq>],
        vector_omega: &[Vec<Zq>],
    ) {
        for i in 0..self.ep.multiplicity {
            let phi_prime_l_vector: Vec<&RqVector> =
                phi_prime.iter().map(|elems| &elems[i]).collect();
            for (k, phi_k) in self.phi_double_prime.iter_mut().enumerate() {
                phi_k[i] =
                    inner_product::compute_linear_combination(&phi_prime_l_vector, &vector_psi[k]);
            }
        }

        for (i, pi_i) in conjugated_pi.iter().enumerate() {
            for (k, phi_k) in self.phi_double_prime.iter_mut().enumerate() {
                phi_k[i] = &phi_k[i]
                    + &inner_product::compute_linear_combination(
                        pi_i.get_elements(),
                        &vector_omega[k],
                    );
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
                    self.a_double_prime[k].get_cell(i, j)
                        * &inner_product::compute_linear_combination(witness[i].get_elements(), witness[j].get_elements())
                })
                .fold(
                    // sum over all i,j
                    Rq::zero(),
                    |acc, val| &acc + &val,
                )
                // add \phi_{i}^{''(k)} * s[i]
                + &inner_product::compute_linear_combination(self.phi_double_prime[k][i].get_elements(), witness[i].get_elements())
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
            aggregated_a: RqMatrix::symmetric_zero(parameters.multiplicity),
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
        for i in 0..self.ep.multiplicity {
            for j in 0..i + 1 {
                let a_constraint_k: Vec<&Rq> = a_constraint
                    .iter()
                    .map(|matrix| matrix.get_cell(i, j))
                    .collect();
                let a_double_prime_k: Vec<&Rq> = a_double_prime
                    .iter()
                    .map(|matrix| matrix.get_cell(i, j))
                    .collect();
                self.aggregated_a.set_cell(
                    i,
                    j,
                    &inner_product::compute_linear_combination::<&Rq, Rq, Rq>(
                        &a_constraint_k,
                        vector_alpha.get_elements(),
                    ) + &inner_product::compute_linear_combination(
                        &a_double_prime_k,
                        vector_beta.get_elements(),
                    ),
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
        for i in 0..self.ep.multiplicity {
            let phi_constraint_k: Vec<&RqVector> =
                phi_constraint.iter().map(|element| &element[i]).collect();
            let phi_double_prime_k: Vec<&RqVector> =
                phi_double_prime.iter().map(|element| &element[i]).collect();
            self.aggregated_phi[i] =
                &inner_product::compute_linear_combination::<&RqVector, RqVector, Rq>(
                    &phi_constraint_k,
                    vector_alpha.get_elements(),
                ) + &inner_product::compute_linear_combination(
                    &phi_double_prime_k,
                    vector_beta.get_elements(),
                );
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
        self.aggregated_b = &inner_product::compute_linear_combination(
            b_constraint.get_elements(),
            vector_alpha.get_elements(),
        ) + &inner_product::compute_linear_combination(
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

#[cfg(test)]
#[allow(clippy::needless_range_loop)]
mod constant_agg_tests {
    use crate::relation::env_params;

    use super::*;

    mod variable_generator {
        use super::*;
        use rand::{
            distr::{Distribution, Uniform},
            rng,
        };
        fn sample_zq_vector(length: usize) -> Vec<Zq> {
            let uniform = Uniform::new_inclusive(Zq::ZERO, Zq::NEG_ONE).unwrap();
            let mut coeffs = vec![Zq::ZERO; length];
            coeffs
                .iter_mut()
                .for_each(|c| *c = uniform.sample(&mut rng()));
            coeffs
        }

        pub fn generate_vector_psi(vec_length: usize, inner_vec_size: usize) -> Vec<Vec<Zq>> {
            let mut vector_psi = Vec::new();
            for _ in 0..vec_length {
                vector_psi.push(sample_zq_vector(inner_vec_size));
            }
            vector_psi
        }

        pub fn generate_a_prime(vec_length: usize, matrix_size: usize) -> Vec<RqMatrix> {
            let mut a_prime = Vec::new();
            for _ in 0..vec_length {
                a_prime.push(RqMatrix::symmetric_random(&mut rng(), matrix_size));
            }
            a_prime
        }

        pub fn generate_phi_prime(
            vec_length: usize,
            inner_vec_length: usize,
            inner_inner_vec_length: usize,
        ) -> Vec<Vec<RqVector>> {
            let mut phi_prime = Vec::new();
            for _ in 0..vec_length {
                let mut inner_vec = Vec::new();
                for _ in 0..inner_vec_length {
                    inner_vec.push(RqVector::random(&mut rng(), inner_inner_vec_length));
                }
                phi_prime.push(inner_vec);
            }
            phi_prime
        }

        pub fn generate_conjugated_pi(
            vec_length: usize,
            matrix_row: usize,
            matrix_col: usize,
        ) -> Vec<RqMatrix> {
            let mut conjugated_pi = Vec::new();
            for _ in 0..vec_length {
                conjugated_pi.push(RqMatrix::random(&mut rng(), matrix_row, matrix_col));
            }
            conjugated_pi
        }

        pub fn generate_omega(vec_length: usize, inner_vec_length: usize) -> Vec<Vec<Zq>> {
            let mut vector_omega = Vec::new();
            for _ in 0..vec_length {
                vector_omega.push(variable_generator::sample_zq_vector(inner_vec_length));
            }
            vector_omega
        }

        pub fn generate_witness(vec_length: usize, inner_vec_length: usize) -> Vec<RqVector> {
            let mut witness = Vec::new();
            for _ in 0..vec_length {
                witness.push(RqVector::random(&mut rng(), inner_vec_length));
            }
            witness
        }
    }

    #[test]
    fn test_calculate_agg_a_double_prime() {
        let params = EnvironmentParameters::default();
        let mut aggregator = ZeroConstantFunctionsAggregation::new(&params);

        let vector_psi =
            variable_generator::generate_vector_psi(params.const_agg_length, params.constraint_l);
        let a_prime =
            variable_generator::generate_a_prime(params.constraint_l, params.multiplicity);

        aggregator.calculate_agg_a_double_prime(&vector_psi, &a_prime);

        for k in 0..params.const_agg_length {
            for i in 0..params.multiplicity {
                for j in 0..params.multiplicity {
                    let mut rhs = Rq::zero();
                    for l in 0..params.constraint_l {
                        rhs = &rhs + &(a_prime[l].get_cell(i, j) * &vector_psi[k][l])
                    }
                    assert_eq!(*aggregator.a_double_prime[k].get_cell(i, j), rhs);
                }
            }
        }
    }

    #[test]
    fn test_calculate_agg_phi_double_prime() {
        let params = EnvironmentParameters::default();
        let mut aggregator = ZeroConstantFunctionsAggregation::new(&params);

        let phi_prime = variable_generator::generate_phi_prime(
            params.constraint_l,
            params.multiplicity,
            params.rank,
        );
        let conjugated_pi = variable_generator::generate_conjugated_pi(
            params.multiplicity,
            2 * env_params::SECURITY_PARAMETER,
            params.rank,
        );
        let vector_psi =
            variable_generator::generate_vector_psi(params.const_agg_length, params.constraint_l);
        let vector_omega = variable_generator::generate_omega(
            params.const_agg_length,
            2 * env_params::SECURITY_PARAMETER,
        );

        aggregator.calculate_agg_phi_double_prime(
            &phi_prime,
            &conjugated_pi,
            &vector_psi,
            &vector_omega,
        );

        for k in 0..params.const_agg_length {
            for i in 0..params.multiplicity {
                let mut rhs = RqVector::zero(params.rank);
                for l in 0..params.constraint_l {
                    rhs = &rhs + &(&phi_prime[l][i] * vector_psi[k][l]);
                }
                let mut lhs = RqVector::zero(params.rank);
                for j in 0..2 * env_params::SECURITY_PARAMETER {
                    lhs = &lhs + &(&conjugated_pi[i].get_elements()[j] * vector_omega[k][j]);
                }
                assert_eq!(aggregator.phi_double_prime[k][i], &rhs + &lhs);
            }
        }
    }

    #[test]
    fn test_calculate_agg_b_double_prime() {
        let params = EnvironmentParameters::default();
        let mut aggregator = ZeroConstantFunctionsAggregation::new(&params);

        let phi_prime = variable_generator::generate_phi_prime(
            params.constraint_l,
            params.multiplicity,
            params.rank,
        );
        let conjugated_pi = variable_generator::generate_conjugated_pi(
            params.multiplicity,
            2 * env_params::SECURITY_PARAMETER,
            params.rank,
        );
        let vector_psi =
            variable_generator::generate_vector_psi(params.const_agg_length, params.constraint_l);
        let vector_omega = variable_generator::generate_omega(
            params.const_agg_length,
            2 * env_params::SECURITY_PARAMETER,
        );
        let a_prime =
            variable_generator::generate_a_prime(params.constraint_l, params.multiplicity);
        let witness_vector = variable_generator::generate_witness(params.multiplicity, params.rank);

        aggregator.calculate_agg_a_double_prime(&vector_psi, &a_prime);
        aggregator.calculate_agg_phi_double_prime(
            &phi_prime,
            &conjugated_pi,
            &vector_psi,
            &vector_omega,
        );
        let b_double_prime = aggregator.calculate_agg_b_double_prime(&witness_vector);

        for k in 0..params.const_agg_length {
            let mut rhs = Rq::zero();
            for i in 0..params.multiplicity {
                for j in 0..params.multiplicity {
                    rhs = &rhs
                        + &(aggregator.a_double_prime[k].get_cell(i, j)
                            * &inner_product::compute_linear_combination(
                                witness_vector[i].get_elements(),
                                witness_vector[j].get_elements(),
                            ));
                }
            }
            let mut lhs = Rq::zero();
            for i in 0..params.multiplicity {
                lhs = &lhs
                    + (&inner_product::compute_linear_combination(
                        aggregator.phi_double_prime[k][i].get_elements(),
                        witness_vector[i].get_elements(),
                    ));
            }
            assert_eq!(b_double_prime.get_elements()[k], &rhs + &lhs);
        }
    }
}

#[cfg(test)]
#[allow(clippy::needless_range_loop)]
mod func_agg_tests {
    use super::*;

    mod variable_generator {
        use super::*;
        use rand::rng;

        pub fn generate_rq_vector(vec_length: usize) -> RqVector {
            let mut vector_alpha = Vec::new();
            for _ in 0..vec_length {
                vector_alpha.push(Rq::random(&mut rng()));
            }
            RqVector::new(vector_alpha)
        }

        pub fn generate_matrix_vector(vec_length: usize, matrix_size: usize) -> Vec<RqMatrix> {
            let mut a_constraint = Vec::new();
            for _ in 0..vec_length {
                a_constraint.push(RqMatrix::symmetric_random(&mut rng(), matrix_size));
            }
            a_constraint
        }

        pub fn generate_rqvector_vector(
            vec_length: usize,
            inner_vec_length: usize,
            inner_inner_vec_length: usize,
        ) -> Vec<Vec<RqVector>> {
            let mut phi_constraint = Vec::new();
            for _ in 0..vec_length {
                let mut inner_vec = Vec::new();
                for _ in 0..inner_vec_length {
                    inner_vec.push(RqVector::random(&mut rng(), inner_inner_vec_length));
                }
                phi_constraint.push(inner_vec);
            }
            phi_constraint
        }
    }

    #[test]
    fn test_calculate_agg_a() {
        let params = EnvironmentParameters::default();
        let mut aggregator = FunctionsAggregation::new(&params);

        let vector_alpha = variable_generator::generate_rq_vector(params.constraint_k);
        let vector_beta = variable_generator::generate_rq_vector(params.const_agg_length);
        let a_constraint =
            variable_generator::generate_matrix_vector(params.constraint_k, params.multiplicity);
        let a_double_prime = variable_generator::generate_matrix_vector(
            params.const_agg_length,
            params.multiplicity,
        );

        aggregator.calculate_agg_a(&a_constraint, &a_double_prime, &vector_alpha, &vector_beta);

        for i in 0..params.multiplicity {
            for j in 0..params.multiplicity {
                let mut rhs = Rq::zero();
                for k in 0..params.constraint_k {
                    rhs = &rhs + &(a_constraint[k].get_cell(i, j) * &vector_alpha.get_elements()[k])
                }
                let mut lhs = Rq::zero();
                for k in 0..params.const_agg_length {
                    lhs =
                        &lhs + &(a_double_prime[k].get_cell(i, j) * &vector_beta.get_elements()[k])
                }
                assert_eq!(*aggregator.aggregated_a.get_cell(i, j), &rhs + &lhs);
            }
        }
    }

    #[test]
    fn test_calculate_agg_phi() {
        let params = EnvironmentParameters::default();
        let mut aggregator = FunctionsAggregation::new(&params);

        let vector_alpha = variable_generator::generate_rq_vector(params.constraint_k);
        let vector_beta = variable_generator::generate_rq_vector(params.const_agg_length);
        let phi_constraint = variable_generator::generate_rqvector_vector(
            params.constraint_k,
            params.multiplicity,
            params.rank,
        );
        let phi_double_prime = variable_generator::generate_rqvector_vector(
            params.const_agg_length,
            params.multiplicity,
            params.rank,
        );

        aggregator.calculate_aggr_phi(
            &phi_constraint,
            &phi_double_prime,
            &vector_alpha,
            &vector_beta,
        );

        for i in 0..params.multiplicity {
            let mut rhs = RqVector::zero(params.rank);
            for k in 0..params.constraint_k {
                rhs = &rhs + &(&phi_constraint[k][i] * &vector_alpha.get_elements()[k])
            }
            let mut lhs = RqVector::zero(params.rank);
            for k in 0..params.const_agg_length {
                lhs = &lhs + &(&phi_double_prime[k][i] * &vector_beta.get_elements()[k])
            }
            assert_eq!(aggregator.aggregated_phi[i], &rhs + &lhs);
        }
    }

    #[test]
    fn test_calculate_agg_b() {
        let params = EnvironmentParameters::default();
        let mut aggregator = FunctionsAggregation::new(&params);

        let vector_alpha = variable_generator::generate_rq_vector(params.constraint_k);
        let vector_beta = variable_generator::generate_rq_vector(params.const_agg_length);
        let b_constraint = variable_generator::generate_rq_vector(params.constraint_k);
        let b_double_prime = variable_generator::generate_rq_vector(params.const_agg_length);

        aggregator.calculate_aggr_b(&b_constraint, &b_double_prime, &vector_alpha, &vector_beta);

        let mut rhs = Rq::zero();
        for k in 0..params.constraint_k {
            rhs = &rhs + &(&b_constraint.get_elements()[k] * &vector_alpha.get_elements()[k])
        }
        let mut lhs = Rq::zero();
        for k in 0..params.const_agg_length {
            lhs = &lhs + &(&b_double_prime.get_elements()[k] * &vector_beta.get_elements()[k])
        }
        assert_eq!(aggregator.aggregated_b, &rhs + &lhs);
    }
}
