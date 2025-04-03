use crate::prover::Challenges;
use crate::utils::{crs::PublicPrams, env_params::EnvironmentParameters, statement::Statement};

pub struct LabradorVerifier<'a> {
    pub pp: &'a PublicPrams,
    pub st: &'a Statement,
    pub tr: &'a Challenges,
}

impl<'a> LabradorVerifier<'a> {
    pub fn new(pp: &'a PublicPrams, st: &'a Statement, tr: &'a Challenges) -> Self {
        Self { pp, st, tr }
    }

    pub fn verify(&self, ep: &EnvironmentParameters) -> bool {
        let _ep = ep;

        // Todo
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prover::{LabradorProver, Witness};

    #[test]
    fn test_verify() {
        // set up example environment, use set1 for testing.
        let ep_1 = EnvironmentParameters::set_1();
        // generate a random witness based on ep above
        let witness_1 = Witness::new(&ep_1);
        // generate public statements based on witness_1
        let st: Statement = Statement::new(&witness_1, &ep_1);
        // generate the common reference string matriices
        let pp = PublicPrams::new(&ep_1);
        // generate random challenges
        let tr = Challenges::new(&ep_1);

        // create a new prover
        let prover = LabradorProver::new(&pp, &witness_1, &st, &tr);
        let _proof = prover.prove(&ep_1);

        // create a new verifier
        let verifier = LabradorVerifier::new(&pp, &st, &tr);
        let result = verifier.verify(&ep_1);
        assert!(result)
    }
}
