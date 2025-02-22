use criterion::{criterion_group, criterion_main, Criterion};
use geo_3d::prelude::*;
use geo_3d::{Coord, CoordNum};

use num_traits::Signed;
use rand::distr::uniform::SampleUniform;
use rand::Rng;

pub fn uniform_points_in_range<S: CoordNum + SampleUniform + Signed, R: Rng>(
    range: S,
    size: usize,
    rng: &mut R,
) -> Vec<Coord<S>> {
    (0..size)
        .map(|_| (rng.random_range(-range..=range), rng.random_range(-range..=range), rng.random_range(-range..=range)).into())
        .collect()
}

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("convex hull f64", |bencher| {
        let line_string = geo_test_fixtures::norway_main::<f64>();

        bencher.iter(|| {
            criterion::black_box(unsafe {
                criterion::black_box(&line_string)
                .convex_hull()
                .unwrap_unchecked()
            });
        });
    });

    c.bench_function("convex hull with collinear random f64", |bencher| {
        let mut points = uniform_points_in_range(10_000_f64, 1_000_000, &mut rand::rng());
        use geo_3d::convex_hull::graham_hull;
        bencher.iter(|| {
            criterion::black_box(graham_hull(
                criterion::black_box(&mut points),
                criterion::black_box(true),
            ));
        });
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
