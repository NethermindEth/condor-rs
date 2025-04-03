use crate::zq::Zq;

// Example Environment parameters used for LaBRADOR, can be expanded as required by testing.
pub struct EnvironmentParameters {
    // r is the number of witness elements, multiplicity
    pub r: usize,
    // n is the witness sizes, rank
    pub n: usize,
    // witness norm bound beta
    pub beta: Zq,
    // decompose base b
    pub b: usize,
    // the parts of decomposition.
    pub digits: usize,
    // t_1
    pub t_1: usize,
    // t_2
    pub t_2: usize,
    // kappa k, the example size, 128/log_q
    pub k: usize,
    pub k_1: usize,
    pub k_2: usize,
    // security level, \lambda
    pub lambda: usize,
    pub lambda2: usize,
    // the log of modulus q, q = 2^(32)
    pub log_q: usize,
    // random Eqlynomial degree bound
    pub deg_bound_d: usize,

    // L: F' functions family size, K: F functions family size
    pub constraint_l: usize,
    pub constraint_k: usize,
}

impl EnvironmentParameters {
    pub fn set_1() -> Self {
        Self {
            r: 4,
            n: 4,
            beta: Zq::new(70),
            b: 4,
            digits: 3,
            t_1: 3,
            t_2: 3,
            k: 4,
            k_1: 5,
            k_2: 5,
            lambda: 128,
            lambda2: 256,
            log_q: 32,
            deg_bound_d: 4,
            constraint_l: 5,
            constraint_k: 5,
        }
    }
}
