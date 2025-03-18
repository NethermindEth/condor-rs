use crate::{
    ajtai_commitment::{
        AjtaiCommitment, AjtaiParameters, CommitError, Opening, ParameterError, VerificationError,
    },
    rq::Rq,
    rq_matrix::RqMatrix,
    rq_vector::RqVector,
    zq::Zq,
};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum HierarchicalError {
    #[error("parameter error: {0}")]
    ParameterError(#[from] ParameterError),
    #[error("commit error: {0}")]
    CommitError(#[from] CommitError),
    #[error("verification error: {0}")]
    VerificationError(#[from] VerificationError),
    #[error("invalid decomposition base: {0}")]
    InvalidBase(Zq),
    #[error("invalid number of parts: {0}")]
    InvalidPartCount(usize),
    #[error("empty witness batch")]
    EmptyWitnessBatch,
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
    pub fn new(base: Zq, num_parts: usize) -> Result<Self, HierarchicalError> {
        if base <= Zq::ONE {
            return Err(HierarchicalError::InvalidBase(base));
        }
        if num_parts == 0 {
            return Err(HierarchicalError::InvalidPartCount(num_parts));
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

/// Proof structure for hierarchical commitments
/// Contains all information needed to verify the commitment structure:
/// - Outer commitment and its opening
/// - All inner commitments and their openings
/// - The flattened decomposed parts for consistency checking
#[derive(Debug, Clone)]
pub struct HierarchicalProof<const M: usize, const N: usize, const D: usize> {
    /// Outer commitment that summarizes inner commitments
    pub outer_commitment: RqVector<M, D>,
    /// Opening information for the outer commitment
    pub outer_opening: Opening<N, D>,
    /// Inner commitments for individual witnesses
    pub inner_commitments: Vec<RqVector<M, D>>,
    /// Opening information for each inner commitment
    pub inner_openings: Vec<Opening<N, D>>,
    /// Flattened parts of all inner commitments
    pub flattened_parts: Vec<Rq<D>>,
}

/// Two-layer commitment structure with inner and outer commitments
/// - M: Size of commitment output
/// - N: Size of witness vector
/// - D: Degree of the polynomials
#[derive(Debug)]
pub struct HierarchicalCommitment<const M: usize, const N: usize, const D: usize> {
    // Inner commitment scheme for individual witnesses
    inner_scheme: AjtaiCommitment<M, N, D>,
    // Outer commitment scheme for summarizing inner commitments
    outer_scheme: AjtaiCommitment<M, N, D>,
    // Parameters for decomposing polynomials
    decomp_params: DecompositionParameters,
}

impl<const M: usize, const N: usize, const D: usize> HierarchicalCommitment<M, N, D> {
    /// Creates a new hierarchical commitment scheme with separate matrices for inner and outer commitments
    /// Both schemes use the same AjtaiParameters for consistent security
    pub fn new(
        params: AjtaiParameters,
        inner_matrix: RqMatrix<M, N, D>,
        outer_matrix: RqMatrix<M, N, D>,
        decomp_params: DecompositionParameters,
    ) -> Result<Self, HierarchicalError> {
        let inner_scheme = AjtaiCommitment::new(params.clone(), inner_matrix)?;
        let outer_scheme = AjtaiCommitment::new(params, outer_matrix)?;

        Ok(Self {
            inner_scheme,
            outer_scheme,
            decomp_params,
        })
    }

    fn decompose_polynomial(&self, poly: &Rq<D>) -> Vec<Rq<D>> {
        poly.decompose(self.decomp_params.base(), self.decomp_params.num_parts())
    }

    /// Helper method to decompose all inner commitments into a single flattened vector
    /// This collects all the decomposed parts from all polynomials in all commitments
    fn decompose_all_commitments(&self, commitments: &[RqVector<M, D>]) -> Vec<Rq<D>> {
        let capacity = commitments.len() * M * self.decomp_params.num_parts();
        let mut all_parts = Vec::with_capacity(capacity);
        for commitment in commitments {
            for poly in commitment.iter() {
                let decomposed_parts = self.decompose_polynomial(poly);
                all_parts.extend(decomposed_parts);
            }
        }
        all_parts
    }

    /// Creates a valid outer witness from decomposed parts
    /// This involves:
    /// 1. Selecting N polynomials from the decomposed parts
    /// 2. Ensuring all coefficients are within the witness bound
    /// 3. Using centered representatives to minimize coefficient size
    fn create_outer_witness(&self, parts: &[Rq<D>]) -> RqVector<N, D> {
        let witness_bound = self.outer_scheme.witness_bound();
        let mut combined_coeffs = Vec::with_capacity(N);

        // Handle the case where we have more or fewer parts than N
        if parts.len() >= N {
            // Take the first N parts
            for part in parts.iter().take(N) {
                combined_coeffs.push(part.clone());
            }
        } else {
            // Use all available parts and pad with zeros
            for part in parts {
                combined_coeffs.push(part.clone());
            }

            // Pad with zeros to reach size N
            while combined_coeffs.len() < N {
                combined_coeffs.push(Rq::zero());
            }
        }

        // Apply witness bounds to ensure all coefficients are valid
        let mut bounded_coeffs = Vec::with_capacity(N);
        for poly in &combined_coeffs {
            // Create a new polynomial with bounded coefficients
            let mut bounded_coeffs_array = [Zq::ZERO; D];
            for (i, coeff) in poly.get_coefficients().iter().enumerate() {
                bounded_coeffs_array[i] = coeff.centered_mod(witness_bound);
            }
            bounded_coeffs.push(Rq::new(bounded_coeffs_array));
        }

        RqVector::from(bounded_coeffs)
    }

    /// Commits to a batch of witnesses using the hierarchical structure
    /// This creates both inner commitments for each witness and an outer commitment
    /// that summarizes them all, providing a more compact representation
    pub fn commit_batch(
        &self,
        witnesses: &[RqVector<N, D>],
    ) -> Result<HierarchicalProof<M, N, D>, HierarchicalError> {
        if witnesses.is_empty() {
            return Err(HierarchicalError::EmptyWitnessBatch);
        }

        // Generate inner commitments and openings
        let mut inner_commitments = Vec::with_capacity(witnesses.len());
        let mut inner_openings = Vec::with_capacity(witnesses.len());

        for witness in witnesses {
            let (commitment, opening) = self.inner_scheme.commit(witness.clone())?;
            inner_commitments.push(commitment);
            inner_openings.push(opening);
        }

        // Decompose inner commitments and create a flattened vector for the outer commitment
        let flattened_parts = self.decompose_all_commitments(&inner_commitments);

        // Create a valid outer witness
        let outer_witness = self.create_outer_witness(&flattened_parts);

        // Create outer commitment
        let (outer_commitment, outer_opening) = self.outer_scheme.commit(outer_witness)?;

        Ok(HierarchicalProof {
            outer_commitment,
            outer_opening,
            inner_commitments,
            inner_openings,
            flattened_parts,
        })
    }

    /// Verifies a hierarchical proof
    /// This ensures:
    /// 1. All inner commitments match their openings
    /// 2. The decomposed parts match what's stored in the proof
    /// 3. The outer commitment correctly summarizes the inner commitments
    /// 4. The outer commitment matches its opening
    pub fn verify(&self, proof: &HierarchicalProof<M, N, D>) -> Result<(), HierarchicalError> {
        // Verify each inner commitment
        for (commitment, opening) in proof
            .inner_commitments
            .iter()
            .zip(proof.inner_openings.iter())
        {
            self.inner_scheme.verify(commitment, opening)?;
        }

        // Verify outer commitment consistency
        // First, decompose inner commitments
        let flattened_parts = self.decompose_all_commitments(&proof.inner_commitments);

        // Compare with stored flattened parts
        if flattened_parts != proof.flattened_parts {
            return Err(HierarchicalError::VerificationError(
                VerificationError::CommitmentMismatch,
            ));
        }

        // Create the same outer witness
        let outer_witness = self.create_outer_witness(&flattened_parts);

        // Create expected outer commitment
        let (expected_outer, _) = self.outer_scheme.commit(outer_witness)?;

        // Check if outer commitment matches expected
        if proof.outer_commitment != expected_outer {
            return Err(HierarchicalError::VerificationError(
                VerificationError::CommitmentMismatch,
            ));
        }

        // Verify outer opening
        self.outer_scheme
            .verify(&proof.outer_commitment, &proof.outer_opening)?;

        Ok(())
    }

    /// Returns a reference to the inner commitment scheme
    pub fn inner_scheme(&self) -> &AjtaiCommitment<M, N, D> {
        &self.inner_scheme
    }

    /// Returns a reference to the outer commitment scheme
    pub fn outer_scheme(&self) -> &AjtaiCommitment<M, N, D> {
        &self.outer_scheme
    }

    /// Returns a reference to the decomposition parameters
    pub fn decomp_params(&self) -> &DecompositionParameters {
        &self.decomp_params
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rng;

    const TEST_M: usize = 8;
    const TEST_N: usize = 8;
    const TEST_D: usize = 4;

    // Test helpers
    fn create_test_scheme() -> HierarchicalCommitment<TEST_M, TEST_N, TEST_D> {
        let beta = Zq::ONE;
        let witness_bound = Zq::ONE;
        let params = AjtaiParameters::new(beta, witness_bound).unwrap();

        let mut rng = rng();
        let inner_matrix = RqMatrix::<TEST_M, TEST_N, TEST_D>::random(&mut rng);
        let outer_matrix = RqMatrix::<TEST_M, TEST_N, TEST_D>::random(&mut rng);

        let decomp_params = DecompositionParameters::new(Zq::new(8), 2).unwrap();

        HierarchicalCommitment::new(params, inner_matrix, outer_matrix, decomp_params).unwrap()
    }

    fn create_random_witnesses(count: usize) -> Vec<RqVector<TEST_N, TEST_D>> {
        let mut rng = rng();
        (0..count)
            .map(|_| RqVector::random_ternary(&mut rng))
            .collect()
    }

    #[test]
    fn test_decomposition_parameters() {
        assert!(DecompositionParameters::new(Zq::ZERO, 2).is_err());
        assert!(DecompositionParameters::new(Zq::TWO, 0).is_err());
        let params = DecompositionParameters::new(Zq::new(8), 3).unwrap();
        assert_eq!(params.base(), Zq::new(8));
        assert_eq!(params.num_parts(), 3);
    }

    #[test]
    fn test_hierarchical_creation() {
        let _ = create_test_scheme();
    }

    #[test]
    fn test_commit_and_verify() {
        let scheme = create_test_scheme();
        let witnesses = create_random_witnesses(3);

        let proof = scheme.commit_batch(&witnesses).unwrap();
        assert!(scheme.verify(&proof).is_ok());
    }

    #[test]
    fn test_tampered_inner_commitment() {
        let scheme = create_test_scheme();
        let witnesses = create_random_witnesses(3);

        let mut proof = scheme.commit_batch(&witnesses).unwrap();

        // Tamper with an inner commitment
        let mut rng = rng();
        proof.inner_commitments[0] = RqVector::random(&mut rng);

        assert!(scheme.verify(&proof).is_err());
    }

    #[test]
    fn test_tampered_outer_commitment() {
        let scheme = create_test_scheme();
        let witnesses = create_random_witnesses(3);

        let mut proof = scheme.commit_batch(&witnesses).unwrap();

        // Tamper with the outer commitment
        let mut rng = rng();
        proof.outer_commitment = RqVector::random(&mut rng);

        assert!(scheme.verify(&proof).is_err());
    }
}
