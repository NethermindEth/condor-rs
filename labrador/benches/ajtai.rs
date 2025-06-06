use criterion::{black_box, criterion_group, criterion_main, Criterion};
use labrador::{
    commitments::common_instances::AjtaiInstances,
    relation::{env_params::EnvironmentParameters, witness::Witness},
};

pub fn ajtai_commit_generation(c: &mut Criterion) {
    let params = EnvironmentParameters::default();
    let scheme = AjtaiInstances::new(&params);
    let witness = Witness::new(params.rank, params.multiplicity, params.beta);

    c.bench_function("Ajtai Commitment Generation", |b| {
        b.iter(|| scheme.commitment_scheme_a.commit(black_box(&witness.s[0])))
    });
}

pub fn ajtai_commit_verification(c: &mut Criterion) {
    let params = EnvironmentParameters::default();
    let scheme = AjtaiInstances::new(&params);
    let witness = Witness::new(params.rank, params.multiplicity, params.beta);
    let commitment = scheme.commitment_scheme_a.commit(&witness.s[0]).unwrap();

    c.bench_function("Ajtai Commitment Verification", |b| {
        b.iter(|| {
            scheme
                .commitment_scheme_a
                .verify(black_box(&commitment), black_box(&witness.s[0]))
        })
    });
}

criterion_group!(
    name = benches;
    config = Criterion::default();
    targets = ajtai_commit_generation, ajtai_commit_verification,
);
criterion_main!(benches);
