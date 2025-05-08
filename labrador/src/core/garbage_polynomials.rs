use crate::ring::{rq_matrix::RqMatrix, rq_vector::RqVector};

pub struct GarbagePolynomials {
    witness_vector: Vec<RqVector>,
    pub g: RqMatrix,
    pub h: RqMatrix,
}

impl GarbagePolynomials {
    pub fn new(witness_vector: Vec<RqVector>) -> Self {
        Self {
            witness_vector,
            g: RqMatrix::new(Vec::new()),
            h: RqMatrix::new(Vec::new()),
        }
    }

    pub fn compute_g(&mut self) {
        let mut g_i = Vec::new();
        for i in 0..self.witness_vector.len() {
            let mut g_ij = Vec::new();
            for j in 0..=i {
                // Only calculate for j ≤ i (upper triangular)
                g_ij.push(&self.witness_vector[i] * &self.witness_vector[j]);
            }
            g_i.push(RqVector::new(g_ij));
        }
        self.g = RqMatrix::new(g_i);
    }

    pub fn compute_h(&mut self, phi: &[RqVector]) {
        let r = self.witness_vector.len();
        let mut h_i = Vec::with_capacity((r * (r + 1)) / 2);

        for i in 0..r {
            let mut h_ij = Vec::new();
            for j in 0..=i {
                // Only calculate for j ≤ i (upper triangular)
                let inner_phi_i_s_j = phi[i].inner_product_poly_vector(&self.witness_vector[j]);
                let inner_phi_j_s_i = phi[j].inner_product_poly_vector(&self.witness_vector[i]);
                h_ij.push(inner_phi_i_s_j + inner_phi_j_s_i);
            }
            h_i.push(RqVector::new(h_ij));
        }
        self.h = RqMatrix::new(h_i);
    }
}
