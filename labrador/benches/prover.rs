use criterion::{criterion_group, criterion_main, Criterion};
use labrador::{
    commitments::common_instances::AjtaiInstances,
    prover::LabradorProver,
    relation::{env_params::EnvironmentParameters, statement::Statement, witness::Witness},
    transcript::sponges::shake::ShakeSponge,
};

pub fn proof_generation_benchmark(c: &mut Criterion) {
    let ep_1 = EnvironmentParameters::default();
    let witness_1 = Witness::new(ep_1.rank, ep_1.multiplicity, ep_1.beta_sq);
    let st: Statement = Statement::new(&witness_1, &ep_1);
    let crs: AjtaiInstances = AjtaiInstances::new(&ep_1);
    let mut prover = LabradorProver::new(&ep_1, &crs, &witness_1, &st);

    c.bench_function("Proof Generation", |b| {
        b.iter(|| prover.prove::<ShakeSponge>())
    });
}

criterion_group!(
    name = benches;
    config = Criterion::default();
    targets = proof_generation_benchmark
);
criterion_main!(benches);
