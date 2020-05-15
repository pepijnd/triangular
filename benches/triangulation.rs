use criterion::{black_box, criterion_group, criterion_main, Benchmark, Criterion};

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use triangular::Triangulation;

fn triangulate(points: &[(f64, f64)]) {
    let t =
        Triangulation::new(points)
            .with_bounds((0.0, 0.0), (50.0, 50.0));
    t.build::<usize>();
}

fn benchmarks(c: &mut Criterion) {
    c.bench(
        "triangulation",
        Benchmark::new("triangulate", move |b| {
            let mut rng = StdRng::from_seed([0; 32]);
            let points: Vec<(f64, f64)> = (0..2500)
                .into_iter()
                .map(|_| (rng.gen_range(0.0, 50.0), rng.gen_range(0.0, 50.0)))
                .collect();
            b.iter(|| {
                triangulate(black_box(&points));
            })
        })
    );
}

criterion_group!(benches, benchmarks);
criterion_main!(benches);
