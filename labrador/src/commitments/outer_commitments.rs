use thiserror::Error;

use super::ajtai_commitment::AjtaiCommitment;
use crate::{
    core::crs::PublicPrams,
    ring::{rq_matrix::RqMatrix, rq_vector::RqVector, zq::Zq},
};

#[derive(Debug, Error)]
pub enum DecompositionError {
    #[error("invalid decomposition base: {0}")]
    InvalidBase(Zq),
    #[error("invalid number of parts: {0}")]
    InvalidPartCount(usize),
}

/// Parameters for polynomial decomposition in hierarchical commitments
/// The base parameter controls how coefficients are decomposed
/// The num_parts parameter determines how many parts each coefficient is split into
#[derive(Debug, Clone)]
pub struct DecompositionParameters {
    base: Zq,
    num_parts: usize,
}

impl DecompositionParameters {
    /// Creates new decomposition parameters with validation
    /// - base must be greater than 1 for meaningful decomposition
    /// - num_parts must be positive to ensure decomposition occurs
    pub fn new(base: Zq, num_parts: usize) -> Result<Self, DecompositionError> {
        if base <= Zq::ONE {
            return Err(DecompositionError::InvalidBase(base));
        }
        if num_parts == 0 {
            return Err(DecompositionError::InvalidPartCount(num_parts));
        }

        Ok(Self { base, num_parts })
    }

    /// Returns the decomposition base
    pub fn base(&self) -> Zq {
        self.base
    }

    /// Returns the number of decomposition parts
    pub fn num_parts(&self) -> usize {
        self.num_parts
    }
}

pub struct OuterCommitment {
    crs: PublicPrams,
    pub u_1: RqVector,
    pub u_2: RqVector,
}

impl OuterCommitment {
    pub fn new(crs: PublicPrams) -> Self {
        Self {
            crs,
            u_1: RqVector::new(Vec::new()),
            u_2: RqVector::new(Vec::new()),
        }
    }

    pub fn compute_u1(
        &mut self,
        t: RqMatrix,
        t_decomposition_params: DecompositionParameters,
        g: RqMatrix,
        g_decomposition_params: DecompositionParameters,
        ajtai_witness_bound: Zq,
    ) {
        let matrix_b = self.crs.matrix_b.clone();
        let commitment_scheme =
            AjtaiCommitment::new(ajtai_witness_bound, ajtai_witness_bound, matrix_b).unwrap();
        let decomposed_t = t.decompose_each_cell(
            t_decomposition_params.base,
            t_decomposition_params.num_parts,
        );
        let u1_left = commitment_scheme.commit(&decomposed_t).unwrap();

        let matrix_c = self.crs.matrix_c.clone();
        // Todo: gamma_1 should be changed to a valid witness bound
        let commitment_scheme =
            AjtaiCommitment::new(ajtai_witness_bound, ajtai_witness_bound, matrix_c).unwrap();
        let decomposed_g = g.decompose_each_cell(
            g_decomposition_params.base,
            g_decomposition_params.num_parts,
        );
        let u2_left = commitment_scheme.commit(&decomposed_g).unwrap();

        self.u_1 = &u1_left + &u2_left;
    }

    pub fn compute_u2(
        &mut self,
        h: RqMatrix,
        h_decomposition_params: DecompositionParameters,
        ajtai_witness_bound: Zq,
    ) {
        let matrix_d = self.crs.matrix_d.clone();
        // Todo: gamma_1 should be changed to a valid witness bound
        let commitment_scheme =
            AjtaiCommitment::new(ajtai_witness_bound, ajtai_witness_bound, matrix_d).unwrap();
        let decomposed_h = h.decompose_each_cell(
            h_decomposition_params.base,
            h_decomposition_params.num_parts,
        );
        self.u_2 = commitment_scheme.commit(&decomposed_h).unwrap();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decomposition_parameters() {
        assert!(DecompositionParameters::new(Zq::ZERO, 2).is_err());
        assert!(DecompositionParameters::new(Zq::TWO, 0).is_err());
        let params = DecompositionParameters::new(Zq::new(8), 3).unwrap();
        assert_eq!(params.base(), Zq::new(8));
        assert_eq!(params.num_parts(), 3);
    }
}
