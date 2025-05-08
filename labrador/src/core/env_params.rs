use crate::ring::zq::Zq;

// Example Environment parameters used for LaBRADOR, can be expanded as required by testing.
pub struct EnvironmentParameters {
    // r is the number of witness elements, multiplicity
    pub r: usize,
    // n is the witness sizes, rank
    pub n: usize,
    // witness norm bound beta
    pub beta: Zq,
    // decompose base b
    pub b: Zq,
    // the parts of decomposition.
    // t_1
    pub t_1: usize,
    // t_2
    pub t_2: usize,
    // kappa k, the example size, 128/log_q
    pub k: usize,
    pub k_1: usize,
    pub k_2: usize,
    // security level, \lambda = 128, lambda2 = 2 * lambda
    pub lambda2: usize,
    // the log of modulus q, q = 2^(32)
    pub log_q: usize,
    // L: F' functions family size, K: F functions family size
    pub constraint_l: usize,
    pub constraint_k: usize,
}

#[allow(clippy::too_many_arguments)]
impl EnvironmentParameters {
    pub fn new(
        r: usize,
        n: usize,
        beta: Zq,
        b: Zq,
        t_1: usize,
        t_2: usize,
        k: usize,
        k_1: usize,
        k_2: usize,
        lambda2: usize,
        log_q: usize,
        constraint_l: usize,
        constraint_k: usize,
    ) -> Self {
        Self {
            r,
            n,
            beta,
            b,
            t_1,
            t_2,
            k,
            k_1,
            k_2,
            lambda2,
            log_q,
            constraint_l,
            constraint_k,
        }
    }
}

impl Default for EnvironmentParameters {
    fn default() -> Self {
        Self {
            r: 3,
            n: 5,
            beta: Zq::new(65535),
            b: Zq::new(16),
            t_1: 4,
            t_2: 4,
            k: 4,
            k_1: 5,
            k_2: 5,
            lambda2: 256,
            log_q: 32,
            constraint_l: 5,
            constraint_k: 5,
        }
    }
}
