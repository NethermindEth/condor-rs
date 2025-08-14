extern crate labrador;

use labrador::commitments::common_instances::AjtaiInstances;
use labrador::prover::LabradorProver;
use labrador::relation::env_params::EnvironmentParameters;
use labrador::relation::statement::Statement;
use labrador::relation::witness::Witness;
use labrador::transcript::sponges::shake::ShakeSponge;
use labrador::verifier::LabradorVerifier;

fn main() {
    let params = EnvironmentParameters::default();
    let witness = Witness::new(params.rank, params.multiplicity, params.beta_sq);
    let statement = Statement::new(&witness, &params);
    let crs = AjtaiInstances::new(&params);

    // Generate the proof using the prover
    let mut prover = LabradorProver::new(&params, &crs, &witness, &statement);
    let proof: labrador::transcript::LabradorTranscript<ShakeSponge> = prover.prove().unwrap();

    // Verify the proof using the verifier
    let mut verifier = LabradorVerifier::new(&params, &crs, &statement);
    let result = verifier.verify(&proof);

    println!("Verification result: {result:?}");
}
