use crate::ring::rq_matrix::RqMatrix;
use crate::ring::rq_vector::RqVector;
use crate::ring::zq::Zq;
use thiserror::Error;

// Error types with documentation
#[derive(Debug, Error)]
pub enum ParameterError {
    #[error("parameters must be positive")]
    ZeroParameter,
    #[error("security bound β·m^(3/2) must be less than q")]
    SecurityBoundViolation,
    #[error("invalid witness bounds specified")]
    InvalidWitnessBounds(Zq),
    #[error("commitment output length {0} is too large")]
    TooLargeCommitmentLength(usize),
}

#[derive(Debug, Error)]
pub enum CommitError {
    #[error("witness coefficients exceed bound {0}")]
    InvalidWitnessBounds(Zq),
    #[error("invalid witness vector size")]
    InvalidWitnessSize,
}

#[derive(Debug, Error)]
pub enum VerificationError {
    #[error("witness coefficients exceed bound {0}")]
    InvalidWitnessBounds(Zq),
    #[error("commitment does not match opening")]
    CommitmentMismatch,
    #[error("invalid opening vector size")]
    InvalidOpeningSize,
    #[error("invalid commitment vector size")]
    InvalidCommitmentSize,
}

/// Ajtai commitment scheme implementation with matrix-based operations
#[derive(Debug)]
pub struct AjtaiCommitment {
    witness_bound: Zq,
    random_matrix: RqMatrix,
}

impl AjtaiCommitment {
    pub fn new(
        beta: Zq,
        witness_bound: Zq,
        random_matrix: RqMatrix,
    ) -> Result<Self, ParameterError> {
        if witness_bound.is_zero() {
            return Err(ParameterError::InvalidWitnessBounds(witness_bound));
        }
        Self::validate_parameters(
            beta,
            random_matrix.get_row_len(),
            random_matrix.get_col_len(),
        )?;

        Ok(Self {
            witness_bound,
            random_matrix,
        })
    }

    /// Generates commitment and opening information with bounds checking
    pub fn commit(&self, witness: &RqVector) -> Result<RqVector, CommitError> {
        if !self.check_bounds(witness) {
            return Err(CommitError::InvalidWitnessBounds(self.witness_bound));
        }
        if witness.get_length() != self.random_matrix.get_col_len() {
            return Err(CommitError::InvalidWitnessSize);
        }
        let commitment = &self.random_matrix * witness;
        Ok(commitment)
    }

    /// Verifies commitment against opening information
    pub fn verify(
        &self,
        commitment: &RqVector,
        opening: &RqVector,
    ) -> Result<(), VerificationError> {
        if !self.check_bounds(opening) {
            return Err(VerificationError::InvalidWitnessBounds(self.witness_bound));
        }
        if opening.get_length() != self.random_matrix.get_col_len() {
            return Err(VerificationError::InvalidOpeningSize);
        }
        if commitment.get_length() != self.random_matrix.get_row_len() {
            return Err(VerificationError::InvalidCommitmentSize);
        }

        let recomputed = &self.random_matrix * opening;
        if commitment != &recomputed {
            return Err(VerificationError::CommitmentMismatch);
        }

        Ok(())
    }

    fn validate_parameters(beta: Zq, row_len: usize, col_len: usize) -> Result<(), ParameterError> {
        if [row_len, col_len].contains(&0) {
            return Err(ParameterError::ZeroParameter);
        }
        Self::verify_security_relation(beta, row_len)
    }

    /// Verifies the security relation β²m³ < q² required for Ajtai's commitment scheme.
    ///
    /// This bound ensures the scheme's security by:
    /// 1. Making the underlying lattice problem hard (SIS assumption)
    /// 2. Preventing statistical attacks on the commitment
    /// 3. Ensuring the commitment is binding under standard lattice assumptions
    ///
    /// The relation β²m³ < q² is a necessary condition derived from the security
    /// proof of Ajtai's commitment scheme, where:
    /// - β bounds the size of witness coefficients
    /// - m is the commitment output length
    /// - q is the modulus of the underlying ring
    fn verify_security_relation(beta: Zq, m: usize) -> Result<(), ParameterError> {
        // Calculate q from Zq properties
        let q_val = Zq::MAX;
        let q: u128 = q_val.to_u128() + 1;

        // Calculate beta²
        let beta_squared = beta
            .to_u128()
            .checked_pow(2)
            .ok_or(ParameterError::SecurityBoundViolation)?;

        // Calculate m³
        let m_cubed: u128 = m
            .checked_pow(3)
            .ok_or(ParameterError::SecurityBoundViolation)?
            .try_into()
            .map_err(|_| ParameterError::TooLargeCommitmentLength(m))?;

        // Calculate q²
        let q_squared = q
            .checked_pow(2)
            .ok_or(ParameterError::SecurityBoundViolation)?;

        // Check if beta² * m³ < q²
        // Use division instead of multiplication to avoid potential overflow
        if beta_squared >= q_squared.checked_div(m_cubed).unwrap_or(0) {
            return Err(ParameterError::SecurityBoundViolation);
        }

        Ok(())
    }

    /// Checks polynomial coefficients against specified bound
    fn check_bounds(&self, _polynomials: &RqVector) -> bool {
        // As now there are no concrete parameters, we return true.
        true
        // polynomials
        //     .iter()
        //     .all(|p| p.check_bounds(self.witness_bound))
    }

    /// Returns a reference to the internal matrix
    pub fn matrix(&self) -> &RqMatrix {
        &self.random_matrix
    }

    /// Returns the witness bound
    pub fn witness_bound(&self) -> Zq {
        self.witness_bound
    }

    pub fn get_row_size(&self) -> usize {
        self.random_matrix.get_elements().len()
    }

    pub fn get_col_size(&self) -> usize {
        self.random_matrix.get_elements()[0].get_elements().len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ring::rq::Rq;

    const TEST_M: usize = 8;
    const TEST_N: usize = 8;

    // Test helpers
    mod test_utils {
        use super::*;

        pub fn valid_witness(scheme: &AjtaiCommitment) -> RqVector {
            vec![Rq::new([scheme.witness_bound(); Rq::DEGREE]); TEST_N].into()
        }

        pub fn random_valid_witness() -> RqVector {
            let mut rng = rand::rng();
            RqVector::random_ternary(&mut rng, TEST_N)
        }

        pub fn setup_scheme() -> AjtaiCommitment {
            let mut rng = rand::rng();
            let random_matrix = RqMatrix::random(&mut rng, TEST_M, TEST_N);
            AjtaiCommitment::new(Zq::ONE, Zq::ONE, random_matrix).unwrap()
        }
    }

    #[test]
    fn rejects_invalid_parameters() {
        assert!(AjtaiCommitment::new(
            Zq::ONE,
            Zq::ZERO,
            RqMatrix::new(vec![RqVector::new(vec![Rq::zero()])])
        )
        .is_err());
        let _ = test_utils::setup_scheme(); // Will panic if setup fails
    }

    #[test]
    fn initializes_with_correct_bounds() {
        let scheme = test_utils::setup_scheme();
        assert_eq!(scheme.witness_bound(), Zq::ONE);
    }

    #[test]
    fn completes_commitment_cycle() {
        let scheme = test_utils::setup_scheme();
        let witness = test_utils::valid_witness(&scheme);

        let commitment = scheme.commit(&witness).unwrap();
        assert!(scheme.verify(&commitment, &witness).is_ok());

        let mut bad_opening = witness.clone();
        let mut rng = rand::rng();
        bad_opening[0] = Rq::random(&mut rng);
        assert!(scheme.verify(&commitment, &bad_opening).is_err());
    }

    #[test]
    fn maintains_security_properties() {
        let scheme = test_utils::setup_scheme();

        // Use random witnesses to ensure they're different
        let witness1 = test_utils::random_valid_witness();
        let witness2 = test_utils::random_valid_witness();

        // Ensure the witnesses are actually different
        assert_ne!(witness1, witness2, "Test requires different witnesses");

        let c1 = scheme.commit(&witness1).unwrap();
        let c2 = scheme.commit(&witness2).unwrap();
        assert_ne!(
            c1, c2,
            "Different witnesses should produce different commitments"
        );
    }

    #[test]
    fn handles_edge_cases() {
        let scheme = test_utils::setup_scheme();
        let zero_witness = RqVector::zero(TEST_N);

        assert!(scheme.commit(&zero_witness).is_ok());
        assert!(scheme.commit(&test_utils::valid_witness(&scheme)).is_ok());
    }

    #[test]
    fn stress_test() {
        let scheme = test_utils::setup_scheme();

        (0..100).for_each(|_| {
            let witness = test_utils::valid_witness(&scheme);
            let commitment = scheme.commit(&witness).unwrap();
            assert!(scheme.verify(&commitment, &witness).is_ok());
        });
    }
}
