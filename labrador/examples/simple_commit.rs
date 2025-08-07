extern crate labrador;

use labrador::commitments::common_instances::AjtaiInstances;
use labrador::relation::env_params::EnvironmentParameters;
use labrador::relation::witness::Witness;

fn main() {
    let params = EnvironmentParameters::default();
    let scheme = AjtaiInstances::new(&params);
    let witness = Witness::new(params.rank, params.multiplicity, params.beta_sq);

    let commitment = scheme.commitment_scheme_a.commit(&witness.s[0]).unwrap();
    let verification_result = scheme.commitment_scheme_a.verify(&commitment, &witness.s[0]);
    println!("Verification result: {verification_result:?}");
}
