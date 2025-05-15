use crate::ring::zq::Zq;

// Example Environment parameters used for LaBRADOR, can be expanded as required by testing.
#[derive(Clone)]
pub struct EnvironmentParameters {
    /// Relation R Parameters
    pub n: usize, // rank (size of each witness s_i)
    pub r: usize, // multiplicity (number of witness elements)

    /// Decomposition Parameters
    pub b: Zq, // z decomposition base
    pub b_1: Zq,    // t_i decomposition base
    pub t_1: usize, // t_i number of parts
    pub b_2: Zq,    // g_ij decomposition base
    pub t_2: usize, // g_ij number of parts

    /// Norm Bounds
    pub beta: Zq, // Bound for witness s_i
    pub gamma: Zq,   // Bound for z
    pub gamma_1: Zq, // Bound for t
    pub gamma_2: Zq, // Bound for g and h

    /// Commitment Matrices Sizes
    pub kappa: usize, // Number of rows in A
    pub kappa_1: usize, // Number of rows in B
    pub kappa_2: usize, // Number of rows in C

    /// Security Parameter
    pub lambda: usize, // security parameter

    /// Function Families Sizes
    pub constraint_k: usize, // Number of constraints of the form f
    pub constraint_l: usize, // Number of constraints of the form f'

    /// Other Parameters
    pub log_q: usize, // Size of log(q) in bits, where q is the modulo
}

#[allow(clippy::too_many_arguments)]
impl EnvironmentParameters {
    pub fn new(
        n: usize,
        r: usize,
        b: Zq,
        b_1: Zq,
        t_1: usize,
        b_2: Zq,
        t_2: usize,
        beta: Zq,
        gamma: Zq,
        gamma_1: Zq,
        gamma_2: Zq,
        kappa: usize,
        kappa_1: usize,
        kappa_2: usize,
        lambda: usize,
        constraint_k: usize,
        constraint_l: usize,
        log_q: usize,
    ) -> Self {
        Self {
            n,
            r,
            b,
            b_1,
            t_1,
            b_2,
            t_2,
            beta,
            gamma,
            gamma_1,
            gamma_2,
            kappa,
            kappa_1,
            kappa_2,
            lambda,
            constraint_k,
            constraint_l,
            log_q,
        }
    }
}

impl Default for EnvironmentParameters {
    fn default() -> Self {
        Self {
            n: 5,
            r: 3,
            b: Zq::new(2),
            b_1: Zq::new(16),
            t_1: 4,
            b_2: Zq::new(16),
            t_2: 4,
            beta: Zq::new(65535),
            gamma: Zq::new(16),
            gamma_1: Zq::new(16),
            gamma_2: Zq::new(16),
            kappa: 4,
            kappa_1: 5,
            kappa_2: 5,
            lambda: 128,
            constraint_k: 5,
            constraint_l: 5,
            log_q: 32,
        }
    }
}

// Todo: Revise and complete the following functionality
// /// Calculate optimal decomposition parameters based on Section 5.4 of the paper
// #[allow(clippy::as_conversions)]
// pub fn calculate_optimal_parameters(
//     n: usize,  // Dimension of the witness vector
//     s: f64,    // Standard deviation of witness coefficients
//     beta: f64, // Ajtai commitment parameter
// ) -> GarbageParameters {
//     // Calculate estimated standard deviations based on Section 5.4
//     // For g_{ij}, std dev is approximately s_g = sqrt(n*d) * s^2
//     // For h_{ij}, std dev is approximately s_h = beta * sqrt(n*d) * s
//     let n_d_sqrt = (n * d) as f64;
//     let s_g = n_d_sqrt * s * s;
//     let s_h = beta * n_d_sqrt * s;

//     // Calculate optimal bases using formula from Section 5.4
//     // B_i ≈ sqrt(12 * s_i)
//     let b2 = ((12.0 * s_g).sqrt() as u32).max(2);
//     let b1 = ((12.0 * s_h).sqrt() as u32).max(2);

//     // Calculate optimal number of parts
//     // t_i ≈ log_B_i(q)
//     let t2 = ((32.0 / (b2 as f64).log2()).ceil() as usize).max(2);
//     let t1 = ((32.0 / (b1 as f64).log2()).ceil() as usize).max(2);

//     GarbageParameters {
//         g_base: Zq::new(b2),
//         g_parts: t2,
//         h_base: Zq::new(b1),
//         h_parts: t1,
//     }
// }

// #[test]
// #[allow(clippy::as_conversions)]
// fn test_optimal_parameters_accuracy() {
//     // Choose small, simple values for manual calculation
//     let n = 4; // Small dimension
//     let d = 4; // Small degree
//     let s = 2.0; // Simple standard deviation
//     let beta = 1.0; // Simple beta value

//     // Manually calculate the expected values according to the paper's formulas
//     let n_d_sqrt = (n * d) as f64; // = 4.0

//     // s_g = sqrt(n*d) * s^2 = 4.0 * 4.0 = 16.0
//     let s_g = n_d_sqrt * s * s;

//     // s_h = beta * sqrt(n*d) * s = 1.0 * 4.0 * 2.0 = 8.0
//     let s_h = beta * n_d_sqrt * s;

//     // b2 ≈ sqrt(12 * s_g) = sqrt(12 * 16.0) = sqrt(192) ≈ 13.856... => 13
//     let expected_b2 = (12.0 * s_g).sqrt() as u32;

//     // b1 ≈ sqrt(12 * s_h) = sqrt(12 * 8.0) = sqrt(96) ≈ 9.798... => 9
//     let expected_b1 = (12.0 * s_h).sqrt() as u32;

//     // For q = 2^32:
//     // t2 ≈ log_b2(q) = log_13(2^32) = 32/log_2(13) ≈ 32/3.7 ≈ 8.65 => 9
//     let expected_t2 = ((32.0 / (expected_b2 as f64).log2()).ceil() as usize).max(2);

//     // t1 ≈ log_b1(q) = log_9(2^32) = 32/log_2(9) ≈ 32/3.17 ≈ 10.09 => 11
//     let expected_t1 = ((32.0 / (expected_b1 as f64).log2()).ceil() as usize).max(2);

//     // Call the function under test
//     let params =
//         GarbagePolynomialCommitment::<TEST_M, TEST_N, TEST_D>::calculate_optimal_parameters(
//             n, d, s, beta,
//         );

//     // Check results
//     assert!(
//         params.g_base.to_u128() >= expected_b2 as u128 - 1
//             && params.g_base.to_u128() <= expected_b2 as u128 + 1,
//         "g_base {}, expected {}",
//         params.g_base.to_u128(),
//         expected_b2
//     );

//     assert!(
//         params.h_base.to_u128() >= expected_b1 as u128 - 1
//             && params.h_base.to_u128() <= expected_b1 as u128 + 1,
//         "h_base {}, expected {}",
//         params.h_base.to_u128(),
//         expected_b1
//     );

//     // For part counts, use approximate comparison due to potential floating point differences
//     assert!(
//         params.g_parts >= expected_t2 - 1 && params.g_parts <= expected_t2 + 1,
//         "g_parts {}, expected {}",
//         params.g_parts,
//         expected_t2
//     );

//     assert!(
//         params.h_parts >= expected_t1 - 1 && params.h_parts <= expected_t1 + 1,
//         "h_parts {}, expected {}",
//         params.h_parts,
//         expected_t1
//     );
// }
