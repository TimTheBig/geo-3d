use criterion::{criterion_group, criterion_main, Criterion};
use geo_3d::prelude::*;
use geo_3d::Polygon;

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("extremes f32", |bencher| {
        let norway = geo_test_fixtures::norway_main::<f32>();
        let polygon = Polygon::new(norway, vec![]);

        bencher.iter(|| {
            criterion::black_box(criterion::black_box(&polygon).extremes());
        });
    });

    c.bench_function("extremes f64", |bencher| {
        let norway = geo_test_fixtures::norway_main::<f64>();
        let polygon = Polygon::new(norway, vec![]);

        bencher.iter(|| {
            criterion::black_box(criterion::black_box(&polygon).extremes());
        });
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
