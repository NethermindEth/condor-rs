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
