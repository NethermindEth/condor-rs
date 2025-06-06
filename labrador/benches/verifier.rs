use criterion::{black_box, criterion_group, criterion_main, Criterion};
use labrador::{
    commitments::common_instances::AjtaiInstances,
    prover::LabradorProver,
    relation::{env_params::EnvironmentParameters, statement::Statement, witness::Witness},
    transcript::sponges::shake::ShakeSponge,
    verifier::LabradorVerifier,
};

pub fn proof_verification_benchmark(c: &mut Criterion) {
    let ep_1 = EnvironmentParameters::default();
    let witness_1 = Witness::new(ep_1.rank, ep_1.multiplicity, ep_1.beta);
    let st: Statement = Statement::new(&witness_1, &ep_1);
    let crs: AjtaiInstances = AjtaiInstances::new(&ep_1);
    let mut prover = LabradorProver::new(&ep_1, &crs, &witness_1, &st);
    let proof: labrador::transcript::LabradorTranscript<ShakeSponge> = prover.prove().unwrap();

    let mut verifier = LabradorVerifier::new(&ep_1, &crs, &st);

    c.bench_function("Proof Verification", |b| {
        b.iter(|| verifier.verify(black_box(&proof)))
    });
}

criterion_group!(
    name = benches;
    config = Criterion::default();
    targets = proof_verification_benchmark
);
criterion_main!(benches);
