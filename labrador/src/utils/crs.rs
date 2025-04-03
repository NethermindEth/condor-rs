use crate::poly::PolyVector;
use crate::utils::{challenge_set::ChallengeSet, env_params::EnvironmentParameters};

pub struct PublicPrams {
    // A \in R_q^(k*n)
    pub matrix_a: Vec<PolyVector>,
    // B_{ik} \in R_q^(k_1*k), for i \in [1,r] and k \in [0, t_1-1]
    pub matrix_b: Vec<Vec<Vec<PolyVector>>>,
    // C_{ijk} \in R_q^(k_2*1), for i \in [1,r], j \in [i, r], and k \in [0, t_2-1]
    pub matrix_c: Vec<Vec<Vec<Vec<PolyVector>>>>,
    // D_{ijk} \in R_q^(k_2*1), for i \in [1,r], j \in [i, r], and k \in [0, t_1-1]
    pub matrix_d: Vec<Vec<Vec<Vec<PolyVector>>>>,
}

#[rustfmt::skip]
impl PublicPrams {
    pub fn new(ep: &EnvironmentParameters) -> Self {
        let cs_matrix_a = ChallengeSet::new(ep.deg_bound_d);
        let cs_matrix_b = ChallengeSet::new(ep.deg_bound_d);
        let cs_matrix_c = ChallengeSet::new(ep.deg_bound_d);
        let cs_matrix_d = ChallengeSet::new(ep.deg_bound_d);

        let matrix_a: Vec<PolyVector> = (0..ep.k).map(|_| {
            (0..ep.n).map(|_| 
                cs_matrix_a.get_challenges().clone()
            ).collect()
        }).collect();

        let matrix_b:Vec<Vec<Vec<PolyVector>>> = (0..ep.r).map(|_|{
            (0..ep.t_1).map(|_|{
                (0..ep.k_1).map(|_|{
                    (0..ep.k).map(|_|{
                        cs_matrix_b.get_challenges().clone()
                    }).collect()
                }).collect()
            }).collect()
        }).collect();

        let matrix_c:Vec<Vec<Vec<Vec<PolyVector>>>> = (0..ep.r).map(|_|{
            (0..ep.r).map(|_|{
                (0..ep.t_2).map(|_|{
                    (0..ep.k_2).map(|_|{
                        (0..1).map(|_|{
                            cs_matrix_c.get_challenges().clone()
                        }).collect()
                    }).collect()
                }).collect()
            }).collect()
        }).collect();

        let matrix_d:Vec<Vec<Vec<Vec<PolyVector>>>> = (0..ep.r).map(|_|{
            (0..ep.r).map(|_|{
                (0..ep.t_1).map(|_|{
                    (0..ep.k_2).map(|_|{
                        (0..1).map(|_|{
                            cs_matrix_d.get_challenges().clone()
                        }).collect()
                    }).collect()
                }).collect()
            }).collect()
        }).collect();

        Self {
            matrix_a,
            matrix_b,
            matrix_c,
            matrix_d,
        }
    }
}
