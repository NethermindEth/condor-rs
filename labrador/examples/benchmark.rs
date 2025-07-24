extern crate labrador;
extern crate criterion;

use labrador::commitments::common_instances::AjtaiInstances;
use labrador::relation::env_params::EnvironmentParameters;
use labrador::relation::witness::Witness;
use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn bench_ajtai_commit(c: &mut Criterion) {
    let params = EnvironmentParameters::default();
    let scheme = AjtaiInstances::new(&params);
    let witness = Witness::new(params.rank, params.multiplicity, params.beta);

    c.bench_function("Ajtai Commitment Generation", |b| {
        b.iter(|| scheme.commitment_scheme_a.commit(black_box(&witness.s[0])))
    });
}

criterion_group!(benches, bench_ajtai_commit);
criterion_main!(benches);