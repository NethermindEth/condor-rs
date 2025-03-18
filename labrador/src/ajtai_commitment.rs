use crate::{rq_matrix::RqMatrix, rq_vector::RqVector, zq::Zq};
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
}

#[derive(Debug, Error)]
pub enum VerificationError {
    #[error("witness coefficients exceed bound {0}")]
    InvalidWitnessBounds(Zq),
    #[error("commitment does not match opening")]
    CommitmentMismatch,
}

/// Configuration parameters for Ajtai commitment scheme with validation invariants
#[derive(Debug, Clone)]
pub struct AjtaiParameters {
    beta: Zq,
    witness_bound: Zq,
}

impl AjtaiParameters {
    /// Creates new parameters with validation
    pub const fn new(beta: Zq, witness_bound: Zq) -> Result<Self, ParameterError> {
        if witness_bound.is_zero() {
            return Err(ParameterError::InvalidWitnessBounds(witness_bound));
        }

        Ok(Self {
            beta,
            witness_bound,
        })
    }

    /// Returns the beta value
    pub const fn beta(&self) -> Zq {
        self.beta
    }

    /// Returns the witness bound
    pub const fn witness_bound(&self) -> Zq {
        self.witness_bound
    }
}

/// Cryptographic opening containing witness
#[derive(Clone, Debug)]
pub struct Opening<const N: usize, const D: usize> {
    pub witness: RqVector<N, D>,
}

impl<const N: usize, const D: usize> Opening<N, D> {
    /// Creates a new opening from a witness
    pub const fn new(witness: RqVector<N, D>) -> Self {
        Self { witness }
    }
}

/// Ajtai commitment scheme implementation with matrix-based operations
#[derive(Debug)]
pub struct AjtaiCommitment<const M: usize, const N: usize, const D: usize> {
    matrix_a: RqMatrix<M, N, D>,
    witness_bound: Zq,
}

// Core implementation with security checks
impl<const M: usize, const N: usize, const D: usize> AjtaiCommitment<M, N, D> {
    /// Creates new commitment scheme with validated parameters
    pub fn new(
        params: AjtaiParameters,
        matrix_a: RqMatrix<M, N, D>,
    ) -> Result<Self, ParameterError> {
        Self::validate_parameters(&params)?;

        Ok(Self {
            matrix_a,
            witness_bound: params.witness_bound,
        })
    }

    /// Generates commitment and opening information with bounds checking
    pub fn commit(
        &self,
        witness: RqVector<N, D>,
    ) -> Result<(RqVector<M, D>, Opening<N, D>), CommitError> {
        if !Self::check_bounds(&witness, self.witness_bound) {
            return Err(CommitError::InvalidWitnessBounds(self.witness_bound));
        }

        let commitment = &self.matrix_a * &witness;
        let opening = Opening::new(witness);

        Ok((commitment, opening))
    }

    /// Verifies commitment against opening information
    pub fn verify(
        &self,
        commitment: &RqVector<M, D>,
        opening: &Opening<N, D>,
    ) -> Result<(), VerificationError> {
        if !Self::check_bounds(&opening.witness, self.witness_bound) {
            return Err(VerificationError::InvalidWitnessBounds(self.witness_bound));
        }

        let recomputed = &self.matrix_a * &opening.witness;
        if commitment != &recomputed {
            return Err(VerificationError::CommitmentMismatch);
        }

        Ok(())
    }

    /// Validates scheme parameters against cryptographic security requirements
    fn validate_parameters(params: &AjtaiParameters) -> Result<(), ParameterError> {
        if [M, N, D].contains(&0) {
            return Err(ParameterError::ZeroParameter);
        }

        Self::verify_security_relation(params.beta, M)
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
    fn check_bounds<const SIZE: usize>(polynomials: &RqVector<SIZE, D>, bound: Zq) -> bool {
        polynomials.iter().all(|p| p.check_bounds(bound))
    }

    /// Returns a reference to the internal matrix
    pub fn matrix(&self) -> &RqMatrix<M, N, D> {
        &self.matrix_a
    }

    /// Returns the witness bound
    pub fn witness_bound(&self) -> Zq {
        self.witness_bound
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rq::Rq;

    const TEST_M: usize = 8;
    const TEST_N: usize = 8;
    const TEST_D: usize = 4;
    type TestAjtai = AjtaiCommitment<TEST_M, TEST_N, TEST_D>;

    // Test helpers
    mod test_utils {
        use super::*;

        pub fn valid_witness(scheme: &TestAjtai) -> RqVector<TEST_N, TEST_D> {
            vec![Rq::new([scheme.witness_bound(); TEST_D]); TEST_N].into()
        }

        pub fn random_valid_witness() -> RqVector<TEST_N, TEST_D> {
            let mut rng = rand::rng();
            RqVector::random_ternary(&mut rng)
        }

        pub fn setup_scheme() -> TestAjtai {
            let mut rng = rand::rng();
            let matrix_a = RqMatrix::random(&mut rng);
            TestAjtai::new(AjtaiParameters::new(Zq::ONE, Zq::ONE).unwrap(), matrix_a).unwrap()
        }
    }

    #[test]
    fn rejects_invalid_parameters() {
        assert!(AjtaiParameters::new(Zq::ONE, Zq::ZERO).is_err());
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

        let (commitment, opening) = scheme.commit(witness).unwrap();
        assert!(scheme.verify(&commitment, &opening).is_ok());

        let mut bad_opening = opening.clone();
        let mut rng = rand::rng();
        bad_opening.witness[0] = Rq::random(&mut rng);
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

        let (c1, _) = scheme.commit(witness1).unwrap();
        let (c2, _) = scheme.commit(witness2).unwrap();
        assert_ne!(
            c1, c2,
            "Different witnesses should produce different commitments"
        );
    }

    #[test]
    fn handles_edge_cases() {
        let scheme = test_utils::setup_scheme();
        let zero_witness = RqVector::zero();

        assert!(scheme.commit(zero_witness).is_ok());
        assert!(scheme.commit(test_utils::valid_witness(&scheme)).is_ok());
    }

    #[test]
    fn stress_test() {
        let scheme = test_utils::setup_scheme();

        (0..100).for_each(|_| {
            let witness = test_utils::valid_witness(&scheme);
            let (commitment, opening) = scheme.commit(witness).unwrap();
            assert!(scheme.verify(&commitment, &opening).is_ok());
        });
    }
}
