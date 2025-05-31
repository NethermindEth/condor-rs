#![allow(clippy::as_conversions)]
use crate::ring::{rq::Rq, zq::Zq};

/// Security Parameter
pub const SECURITY_PARAMETER: usize = 128;
pub const OPERATOR_NORM: f64 = 15.0;

// Example Environment parameters used for LaBRADOR, can be expanded as required by testing.
#[derive(Clone)]
pub struct EnvironmentParameters {
    /// Relation R Parameters
    pub rank: usize, // size of each witness s_i
    pub multiplicity: usize, // number of witness elements

    /// Decomposition Parameters
    pub b: Zq, // z decomposition base
    pub b_1: usize, // t_i decomposition base
    pub t_1: usize, // t_i number of parts
    pub b_2: usize, // g_ij decomposition base
    pub t_2: usize, // g_ij number of parts

    /// Norm Bounds
    pub beta: Zq, // Bound for witness s_i
    pub gamma: Zq,   // Bound for z
    pub gamma_1: Zq, // Bound for t
    pub gamma_2: Zq, // Bound for g and h
    pub beta_prime: Zq,

    /// Commitment Matrices Sizes
    pub kappa: usize, // Number of rows in A
    pub kappa_1: usize, // Number of rows in B
    pub kappa_2: usize, // Number of rows in C

    /// Function Families Sizes
    pub constraint_k: usize, // Number of constraints of the form f
    pub constraint_l: usize,     // Number of constraints of the form f'
    pub const_agg_length: usize, // Number of functions in the first aggregation step

    /// Other Parameters
    pub log_q: usize, // Size of log(q) in bits, where q is the modulo
}

impl EnvironmentParameters {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        rank: usize,
        multiplicity: usize,
        beta: f64,
        kappa: usize,
        kappa_1: usize,
        kappa_2: usize,
        constraint_k: usize,
        constraint_l: usize,
        tau: f64,
        modulo_q: usize,
    ) -> Self {
        let std_s = Self::compute_std_s(rank, multiplicity, beta);

        // balanced radix for z‑split
        let base_b = Self::compute_base_b(std_s, multiplicity, tau);

        // Decomposition parameters for t and h
        let parts_t1 = Self::compute_t1(modulo_q, base_b);
        let base_b1 = Self::compute_b1(modulo_q, parts_t1);

        // Decomposition parameters for g
        let parts_t2 = Self::compute_t2(rank, std_s, base_b);
        let base_b2 = Self::compute_b2(rank, std_s, parts_t2);

        let gamma = Self::compute_gamma(beta, tau);
        let gamma_1 =
            Self::compute_gamma1(base_b1, parts_t1, multiplicity, kappa, base_b2, parts_t2);
        let gamma_2 = Self::compute_gamma2(base_b1, parts_t1, multiplicity);
        let beta_prime = Self::compute_beta_prime(base_b, gamma, gamma_1, gamma_2);
        let log_q = (modulo_q as f64).log2() as usize;
        EnvironmentParameters {
            rank,
            multiplicity,
            b: Zq::new(base_b as u32),
            b_1: base_b1,
            t_1: parts_t1,
            b_2: base_b2,
            t_2: parts_t2,
            beta: Zq::new(beta as u32),
            gamma: Zq::new(gamma as u32),
            gamma_1: Zq::new(gamma_1 as u32),
            gamma_2: Zq::new(gamma_2 as u32),
            beta_prime: Zq::new(beta_prime as u32),
            kappa,
            kappa_1,
            kappa_2,
            constraint_k,
            constraint_l,
            const_agg_length: usize::div_ceil(SECURITY_PARAMETER, log_q),
            log_q, // Size of log(q) in bits, where q is the modulo
        }
    }

    // s = β / √(r·n·d)
    fn compute_std_s(n: usize, r: usize, beta: f64) -> f64 {
        beta / ((r * n * Rq::DEGREE) as f64).sqrt()
    }

    /// b = round( √(√(12 r τ) · s) ), but at least 2
    fn compute_base_b(s: f64, multiplicity: usize, tau: f64) -> usize {
        let raw = ((12.0 * multiplicity as f64 * tau).sqrt() * s).sqrt();
        let mut b = raw.round() as usize;
        if b < 2 {
            b = 2;
        }
        b
    }

    /// t₁ = round( log_q / log_b )
    fn compute_t1(q: usize, b: usize) -> usize {
        let log_q = (q as f64).ln();
        let log_b = (b as f64).ln();
        (log_q / log_b).round() as usize
    }

    /// b₁ = ceil( q^(1/t₁) )
    fn compute_b1(q: usize, t1: usize) -> usize {
        let root = (q as f64).powf(1.0 / t1 as f64);
        root.ceil() as usize
    }

    /// t₂ = ceil( log(√(24 n d) · s^2) / log(b) )
    fn compute_t2(n: usize, std_s: f64, b: usize) -> usize {
        (((24.0 * n as f64 * Rq::DEGREE as f64).sqrt() * std_s * std_s).ln() / (b as f64).ln())
            .round() as usize
    }

    /// t₂ = ceil( log(√(24 n d) · s^2) / log(b) )
    fn compute_b2(n: usize, std_s: f64, t2: usize) -> usize {
        ((24.0 * n as f64 * Rq::DEGREE as f64).sqrt() * std_s * std_s)
            .powf(1.0 / t2 as f64)
            .round() as usize
    }

    /// γ = β·√τ
    fn compute_gamma(beta: f64, tau: f64) -> f64 {
        beta * (tau).sqrt()
    }

    /// γ₁ = √( (b₁² t₁ /12) · r κ d  +  (b₂² t₂ /12) · (r²+r)/2 · d )
    #[allow(clippy::too_many_arguments)]
    fn compute_gamma1(
        b1: usize,
        t1: usize,
        multiplicity: usize,
        kappa: usize,
        b2: usize,
        t2: usize,
    ) -> f64 {
        let term1 =
            (b1.pow(2) * t1) as f64 / 12.0 * multiplicity as f64 * kappa as f64 * Rq::DEGREE as f64;
        let term2 = (b2.pow(2) * t2) as f64 / 12.0
            * ((multiplicity * multiplicity + multiplicity) as f64)
            / 2.0
            * Rq::DEGREE as f64;
        (term1 + term2).sqrt()
    }

    /// γ₂ = √( (b₁² t₁ /12) · (r²+r)/2 · d )
    fn compute_gamma2(b1: usize, t1: usize, multiplicity: usize) -> f64 {
        let term = (b1.pow(2) * t1) as f64 / 12.0
            * ((multiplicity * multiplicity + multiplicity) as f64)
            / 2.0
            * Rq::DEGREE as f64;
        term.sqrt()
    }

    /// β' = √( 2 γ² / b²  +  γ₁²  +  γ₂² )
    fn compute_beta_prime(b: usize, gamma: f64, gamma1: f64, gamma2: f64) -> f64 {
        let part_z = 2.0 * (gamma * gamma) / (b * b) as f64;
        (part_z + gamma1 * gamma1 + gamma2 * gamma2).sqrt()
    }
}

impl Default for EnvironmentParameters {
    fn default() -> Self {
        Self::new(5, 3, 65535.0, 4, 5, 5, 5, 5, 65535.0, (1u64 << 32) as usize)
    }
}

#[cfg(test)]
mod tests {
    use crate::ring::zq::Zq;

    use super::EnvironmentParameters;
    use super::Rq;

    #[test]
    fn test_std_s_matches_formula() {
        let (n, r, beta) = (12, 5, 42.0);
        let result = EnvironmentParameters::compute_std_s(n, r, beta);
        let expected = beta / ((r * n * Rq::DEGREE) as f64).sqrt();
        assert_eq!(result, expected);
    }

    #[test]
    fn test_compute_base_b_not_below_two() {
        let b = EnvironmentParameters::compute_base_b(0.001, 1, 0.0001);
        assert!(b >= 2);
    }

    #[test]
    fn test_base_b_formula_and_floor() {
        let s = 0.17;
        let (r, tau) = (3, 12.0);
        let result = EnvironmentParameters::compute_base_b(s, r, tau);
        let expected = ((12.0 * r as f64 * tau).sqrt() * s).sqrt().round().max(2.0) as usize;
        assert_eq!(result, expected, "balanced radix b");
    }

    #[test]
    fn test_t1_rounding() {
        let (q, b) = (257usize, 4usize); // ln(q)/ln(b)=~4.006
        let result = EnvironmentParameters::compute_t1(q, b);
        assert_eq!(result, 4, "t1 = round(log_q / log_b)");
    }

    #[test]
    fn test_b1_is_ceil_root() {
        let (q, t1) = (1usize << 20, 5usize);
        let result = EnvironmentParameters::compute_b1(q, t1);
        let expected = (q as f64).powf(1.0 / t1 as f64).ceil() as usize;
        assert_eq!(result, expected);
    }

    #[test]
    fn test_t1_and_b1_are_consistent() {
        let q = 1usize << 23;
        let b = 137usize;
        let t1 = EnvironmentParameters::compute_t1(q, b);
        let b1 = EnvironmentParameters::compute_b1(q, t1);
        assert!(b1.pow(t1 as u32) >= q);
    }

    #[test]
    fn test_t2_b2_pair_match() {
        let (n, s, b) = (10usize, 0.12, 150usize);
        let t2: usize = EnvironmentParameters::compute_t2(n, s, b);
        let b2 = EnvironmentParameters::compute_b2(n, s, t2);
        // recompute tmp and double-check injectivity analogue
        let tmp_ln = ((24.0 * n as f64 * Rq::DEGREE as f64).sqrt() * s * s).ln();
        let expected_t2 = (tmp_ln / (b as f64).ln()).round() as usize;
        assert_eq!(t2, expected_t2, "t2 formula");
        let lhs = (b2 as f64).powf(t2 as f64);
        let rhs = (24.0 * n as f64 * Rq::DEGREE as f64).sqrt() * s * s;
        assert!(lhs >= rhs - 1.0);
    }

    #[test]
    fn test_gamma_functions() {
        let (beta, tau) = (32.0, 49.0);
        let gamma = EnvironmentParameters::compute_gamma(beta, tau);
        assert_eq!(gamma, beta * tau.sqrt());

        // simple numeric spot-check for γ₁, γ₂
        let (b1, t1, r, kappa, b2, t2) = (7, 3, 2, 4, 5, 2);
        let manual1 = ((b1 * b1 * t1) as f64 / 12.0 * r as f64 * kappa as f64 * Rq::DEGREE as f64
            + (b2 * b2 * t2) as f64 / 12.0 * ((r * r + r) as f64) / 2.0 * Rq::DEGREE as f64)
            .sqrt();
        let got1 = EnvironmentParameters::compute_gamma1(b1, t1, r, kappa, b2, t2);
        assert_eq!(got1, manual1);

        let manual2 = ((b1.pow(2) * t1) as f64 / 12.0 * ((r * r + r) as f64) / 2.0
            * Rq::DEGREE as f64)
            .sqrt();
        let got2 = EnvironmentParameters::compute_gamma2(b1, t1, r);
        assert_eq!(got2, manual2);
    }

    #[test]
    fn beta_prime_formula() {
        let (b, gamma, g1, g2) = (11usize, 3.3, 7.7, 2.2);
        let expected = ((2.0 * gamma * gamma) / (b * b) as f64 + g1 * g1 + g2 * g2).sqrt();
        let result = EnvironmentParameters::compute_beta_prime(b, gamma, g1, g2);
        assert_eq!(expected, result);
    }

    #[test]
    fn test_gamma_and_beta_prime_are_positive() {
        let beta = 100.0;
        let tau = 64.0;
        let gamma = EnvironmentParameters::compute_gamma(beta, tau);
        assert!(gamma > 0.0);
        let beta_p = EnvironmentParameters::compute_beta_prime(3, gamma, 2.0, 1.0);
        assert!(beta_p > 0.0);
    }

    #[test]
    fn default_instance_is_consistent() {
        let p = EnvironmentParameters::default();

        assert!(p.b >= Zq::TWO && p.b_1 >= 2 && p.b_2 >= 2);
        assert!(p.t_1 > 0 && p.t_2 > 0);

        assert!((p.b_1 as u64).pow(p.t_1 as u32) >= p.log_q as u64);
        assert!((p.b_2 as u64).pow(p.t_2 as u32) > 1);

        // c)  norm bounds are non-negative
        assert!(p.gamma > Zq::ZERO && p.gamma_1 > Zq::ZERO && p.gamma_2 > Zq::ZERO);
        assert!(p.beta_prime >= p.gamma_1.max(p.gamma_2));
    }

    #[test]
    fn random_parameter_sets_hold_invariants() {
        for seed in 1..=20 {
            let rank = 2 + seed as usize % 8; // 2..9
            let multiplicity = 1 + seed as usize % 5; // 1..5
            let beta = 10.0 + seed as f64;
            let tau = 32.0 + seed as f64;
            let q = (1usize << 20) + (seed as usize * 12345);
            let params =
                EnvironmentParameters::new(rank, multiplicity, beta, 4, 4, 4, 3, 3, tau, q);

            // main invariants
            assert!(params.b >= Zq::TWO);
            assert!(params.t_1 >= 1 && params.t_2 >= 1);
            assert!((params.b_1 as u64).pow(params.t_1 as u32) >= q as u64);
            assert!(params.gamma > Zq::ZERO);
            assert!(params.beta_prime >= params.gamma_1.max(params.gamma_2));
        }
    }
}
