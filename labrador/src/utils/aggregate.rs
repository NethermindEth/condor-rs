use crate::{
    jl::{ProjectionVector, SECURITY_LEVEL},
    rq::Rq,
    rq_vector::RqVector,
    zq::Zq,
    zq_vector::ZqMatrix,
};
pub struct Aggregation<const N: usize, const D: usize> {
    polynomials: RqVector<N, D>,
}

impl<const N: usize, const D: usize> Aggregation<N, D> {
    // modulus of the ring, 2^32
    pub const Q: u128 = 42949672962;
    pub const SECURITY_LEVEL: u32 = 128;

    pub fn new() -> Self {
        Self {
            polynomials: RqVector::zero(),
        }
    }

    pub fn get_polynomials(&self) -> &RqVector<N, D> {
        &self.polynomials
    }

    pub fn aggregate(self, projection: ProjectionVector<N, D>) {
        let k = Self::SECURITY_LEVEL.div_ceil(Self::Q.ilog2());
        let l = self.polynomials.as_slice().len();
        let mut aggregated_polynomials_1: Vec<Rq<D>> = Vec::with_capacity(k.try_into().unwrap());
        let mut aggregated_polynomials_2: Zq = Zq::ZERO;

        for i in 0..k.try_into().unwrap() {
            let random_phi: ZqMatrix = ZqMatrix::random_matrix(k.try_into().unwrap(), l);
            for j_1 in 0..l {
                let mut aggregated_polynomial_1 = Rq::new([Zq::ZERO; D]);
                aggregated_polynomial_1 = Self::linear_combination(
                    aggregated_polynomial_1,
                    random_phi.get_row(j_1).to_rq(),
                );
                aggregated_polynomials_1.push(aggregated_polynomial_1);
            }

            let random_omega: ZqMatrix =
                ZqMatrix::random_matrix(k.try_into().unwrap(), 2 * SECURITY_LEVEL);
            let conjugate_pi: ZqMatrix = random_omega.conjugate_automorphism();
            for j_2 in 0..2 * SECURITY_LEVEL {
                let project_constraints = Self::generate_project_constraints(
                    conjugate_pi.get_row(j_2).to_rq(),
                    projection.clone(),
                    self.get_polynomials().as_slice()[i].clone(),
                );
                let aggregated_polynomial_2 =
                    project_constraints * random_omega.get_row(i).get_coeffs()[j_2];
                aggregated_polynomials_2 += aggregated_polynomial_2;
            }
        }
    }

    fn linear_combination(poly: Rq<D>, challenges: Rq<D>) -> Rq<D> {
        let mut result = [Zq::ZERO; D];
        for (i, val) in result.iter_mut().enumerate() {
            *val = poly.get_coefficients()[i] * challenges.get_coefficients()[i];
        }

        Rq::new(result)
    }

    /// \sum_{i=1}^r<\sigma_{-1}(\pi_i^j), \omega_i)> - p_j. Page 12, Aggregating
    fn generate_project_constraints(
        random_matrix_row: Rq<D>,
        projection: ProjectionVector<N, D>,
        s_i: Rq<D>,
    ) -> Zq {
        let projection_coefficients = projection.get_projection().get_coefficients();
        let mut result = Zq::ZERO;
        for p_j in projection_coefficients.iter() {
            let coeff = random_matrix_row.inner_product(&s_i);
            result = coeff - *p_j;
        }
        result
    }
}

impl<const N: usize, const D: usize> Default for Aggregation<N, D> {
    fn default() -> Self {
        Self::new()
    }
}
