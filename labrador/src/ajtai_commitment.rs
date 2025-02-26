use crate::{rq::Rq, rq_matrix::RqMatrix, zq::Zq};
use thiserror::Error;

// Error types with clear documentation
#[derive(Debug, Error)]
pub enum ParameterError {
    #[error("parameters must be positive")]
    ZeroParameter,
    #[error("security bound β·m^(3/2) must be less than q")]
    SecurityBoundViolation,
    #[error("invalid witness bounds specified")]
    InvalidWitnessBounds,
}

#[derive(Debug, Error)]
pub enum CommitError {
    #[error("witness coefficients exceed bound")]
    WitnessBoundViolation,
}

/// Configuration parameters for Ajtai commitment scheme with validation invariants
#[derive(Debug, Clone)]
pub struct AjtaiParameters {
    beta: Zq,
    witness_bound: Zq,
}

impl AjtaiParameters {
    /// Creates new parameters with validation
    pub fn new(beta: Zq, witness_bound: Zq) -> Result<Self, ParameterError> {
        if witness_bound.is_zero() {
            return Err(ParameterError::InvalidWitnessBounds);
        }

        Ok(Self {
            beta,
            witness_bound,
        })
    }

    /// Default parameters for ternary distribution {-1, 0, 1}
    pub fn ternary() -> Self {
        Self {
            beta: Zq::one(),
            witness_bound: Zq::one(),
        }
    }
}

/// Ajtai commitment scheme implementation with matrix-based operations
#[derive(Debug)]
pub struct AjtaiCommitment<const M: usize, const N: usize, const D: usize> {
    matrix_a: RqMatrix<M, N, D>,
    witness_bound: Zq,
}

/// Cryptographic opening containing witness
#[derive(Clone, Debug)]
pub struct Opening<const N: usize, const D: usize> {
    pub witness: [Rq<D>; N],
}

// Implement Default trait for more idiomatic Rust
impl<const M: usize, const N: usize, const D: usize> Default for AjtaiCommitment<M, N, D> {
    fn default() -> Self {
        Self::new(AjtaiParameters::ternary()).expect("Default parameters should always be valid")
    }
}

// Core implementation with security checks
impl<const M: usize, const N: usize, const D: usize> AjtaiCommitment<M, N, D> {
    /// Creates new commitment scheme with validated parameters
    pub fn new(params: AjtaiParameters) -> Result<Self, ParameterError> {
        Self::validate_parameters(&params)?;
        Ok(Self {
            matrix_a: RqMatrix::random(),
            witness_bound: params.witness_bound,
        })
    }

    /// Generates commitment and opening information with bounds checking
    pub fn commit(&self, witness: [Rq<D>; N]) -> Result<([Rq<D>; M], Opening<N, D>), CommitError> {
        if !Self::check_bounds(&witness, self.witness_bound) {
            return Err(CommitError::WitnessBoundViolation);
        }

        let commitment = self.matrix_a.mul_vec(&witness);

        Ok((commitment, Opening { witness }))
    }

    /// Verifies commitment against opening information
    pub fn verify(&self, commitment: &[Rq<D>; M], opening: &Opening<N, D>) -> bool {
        let bounds_valid = Self::check_bounds(&opening.witness, self.witness_bound);

        bounds_valid && self.verify_commitment_calculation(commitment, opening)
    }

    /// Validates scheme parameters against cryptographic security requirements
    fn validate_parameters(params: &AjtaiParameters) -> Result<(), ParameterError> {
        if [M, N, D].iter().any(|&v| v == 0) {
            return Err(ParameterError::ZeroParameter);
        }

        Self::verify_security_relation(
            params.beta.value(),
            u128::try_from(M).unwrap(),
            Self::modulus_u128(),
        )
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
    /// - β bounds the size of randomness/witness coefficients
    /// - m is the commitment output length
    /// - q is the modulus of the underlying ring
    fn verify_security_relation(beta: u32, m: u128, q: u128) -> Result<(), ParameterError> {
        let beta = u128::from(beta);
        let m_cubed = m
            .checked_pow(3)
            .ok_or(ParameterError::SecurityBoundViolation)?;
        let beta_squared = beta
            .checked_pow(2)
            .ok_or(ParameterError::SecurityBoundViolation)?;
        let q_squared = q
            .checked_pow(2)
            .ok_or(ParameterError::SecurityBoundViolation)?;

        if beta_squared
            .checked_mul(m_cubed)
            .map(|left| left >= q_squared)
            .unwrap_or(true)
        {
            Err(ParameterError::SecurityBoundViolation)
        } else {
            Ok(())
        }
    }

    /// Checks polynomial coefficients against specified bound
    fn check_bounds<const SIZE: usize>(polynomials: &[Rq<D>; SIZE], bound: Zq) -> bool {
        polynomials.iter().all(|p| p.check_bounds(bound))
    }

    /// Recomputes commitment from opening and verifies match
    fn verify_commitment_calculation(
        &self,
        commitment: &[Rq<D>; M],
        opening: &Opening<N, D>,
    ) -> bool {
        let recomputed = self.matrix_a.mul_vec(&opening.witness);
        commitment == &recomputed
    }

    fn modulus_u128() -> u128 {
        let q_val = (Zq::zero() - Zq::one()).value();
        u128::from(q_val) + 1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TEST_M: usize = 8;
    const TEST_N: usize = 8;
    const TEST_D: usize = 4;
    type TestAjtai = AjtaiCommitment<TEST_M, TEST_N, TEST_D>;

    // Test helpers
    mod test_utils {
        use super::*;

        pub fn valid_witness(scheme: &TestAjtai) -> [Rq<TEST_D>; TEST_N] {
            std::array::from_fn(|_| Rq::new(std::array::from_fn(|_| scheme.witness_bound)))
        }

        pub fn random_valid_witness() -> [Rq<TEST_D>; TEST_N] {
            std::array::from_fn(|_| Rq::random_small())
        }

        pub fn setup_scheme() -> TestAjtai {
            TestAjtai::new(AjtaiParameters::ternary()).unwrap()
        }
    }

    #[test]
    fn rejects_invalid_parameters() {
        assert!(AjtaiParameters::new(Zq::one(), Zq::zero()).is_err());
        let _ = test_utils::setup_scheme(); // Will panic if setup fails
    }

    #[test]
    fn initializes_with_correct_bounds() {
        let scheme = TestAjtai::default();
        assert_eq!(scheme.witness_bound.value(), 1);
    }

    #[test]
    fn completes_commitment_cycle() {
        let scheme = test_utils::setup_scheme();
        let witness = test_utils::valid_witness(&scheme);

        let (commitment, opening) = scheme.commit(witness).unwrap();
        assert!(scheme.verify(&commitment, &opening));

        let mut bad_opening = opening.clone();
        bad_opening.witness[0] = Rq::random_small();
        assert!(!scheme.verify(&commitment, &bad_opening));
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
        let zero_witness = std::array::from_fn(|_| Rq::zero());

        assert!(scheme.commit(zero_witness).is_ok());
        assert!(scheme.commit(test_utils::valid_witness(&scheme)).is_ok());
    }

    #[test]
    fn stress_test() {
        let scheme = TestAjtai::default();

        (0..100).for_each(|_| {
            let witness = test_utils::valid_witness(&scheme);
            let (commitment, opening) = scheme.commit(witness).unwrap();
            assert!(scheme.verify(&commitment, &opening));
        });
    }
}
