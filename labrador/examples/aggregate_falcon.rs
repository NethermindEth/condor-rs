use labrador::{
    commitments::common_instances::AjtaiInstances,
    prover::LabradorProver,
    relation::{env_params::EnvironmentParameters, statement::Statement, witness::Witness},
    transcript::sponges::shake::ShakeSponge,
    verifier::LabradorVerifier,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Step 1: Set up LABRADOR parameters
    let params = EnvironmentParameters::default();

    // Step 2: Encode Falcon signatures as witnesses
    let mut witnesses = Vec::new();
    for _ in 0..4 {
        let witness = Witness::new(params.rank, params.multiplicity, params.beta_sq);
        witnesses.push(witness);
    }

    // Step 3: Formulate verification constraints as a Statement
    let statement = Statement::new(&witnesses[0], &params);

    // Step 4: Initialize the Ajtai commitment scheme instances
    let ajtai_instances = AjtaiInstances::new(&params);

    // Step 5: Create and run the LABRADOR prover
    let mut prover = LabradorProver::new(&params, &ajtai_instances, &witnesses[0], &statement);

    // Generate the proof
    let proof = prover.prove::<ShakeSponge>()?;

    // Step 6: Verify the aggregate signature using the LABRADOR verifier
    let mut verifier = LabradorVerifier::new(&params, &ajtai_instances, &statement);

    // Verify the proof
    let is_valid = verifier.verify(&proof)?;

    println!(
        "Aggregate signature verification result: {}",
        if is_valid { "Valid" } else { "Invalid" }
    );

    Ok(())
}
