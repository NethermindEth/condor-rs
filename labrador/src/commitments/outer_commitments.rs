use thiserror::Error;

use crate::{
    commitments::common_instances::AjtaiInstances,
    ring::{rq_matrix::RqMatrix, rq_vector::RqVector, zq::Zq},
};

use super::ajtai_commitment::AjtaiScheme;

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

fn decompose_and_commit(
    commitment_matrix: &AjtaiScheme,
    input: &RqMatrix,
    params: &DecompositionParameters,
) -> RqVector {
    let decomposed_input = input.decompose_each_cell(params.base, params.num_parts);
    commitment_matrix
        .commit(&decomposed_input)
        .expect("Commitment error in committing to decomposed input")
}

pub fn compute_u1(
    crs: &AjtaiInstances,
    t: &RqMatrix,
    t_decomposition_params: DecompositionParameters,
    g: &RqMatrix,
    g_decomposition_params: DecompositionParameters,
) -> RqVector {
    &decompose_and_commit(&crs.commitment_scheme_b, t, &t_decomposition_params)
        + &decompose_and_commit(&crs.commitment_scheme_c, g, &g_decomposition_params)
}

pub fn compute_u2(
    crs: &AjtaiInstances,
    h: &RqMatrix,
    h_decomposition_params: DecompositionParameters,
) -> RqVector {
    decompose_and_commit(&crs.commitment_scheme_d, h, &h_decomposition_params)
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
