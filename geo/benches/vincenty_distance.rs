use criterion::{criterion_group, criterion_main};
use geo_3d::algorithm::VincentyDistance;
use geo_3d::Point;

fn criterion_benchmark(c: &mut criterion::Criterion) {
    c.bench_function("vincenty distance f32", |bencher| {
        let a = Point::<f32>::new(17.107558, 48.148636, 17.107558);
        let b = Point::<f32>::new(16.372477, 48.20881, 16.372477);

        bencher.iter(|| {
            let _ = criterion::black_box(
                criterion::black_box(&a).vincenty_distance(criterion::black_box(&b)),
            );
        });
    });

    c.bench_function("vincenty distance f64", |bencher| {
        let a = Point::new(17.107558, 48.148636, 17.107558);
        let b = Point::new(16.372477, 48.208810, 16.372477);

        bencher.iter(|| {
            let _ = criterion::black_box(
                criterion::black_box(&a).vincenty_distance(criterion::black_box(&b)),
            );
        });
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
