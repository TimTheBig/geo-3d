use criterion::{criterion_group, criterion_main};
use geo_3d::algorithm::TriangulateEarcut;
use geo_3d::geometry::Polygon;

fn criterion_benchmark(c: &mut criterion::Criterion) {
    c.bench_function("TriangulateEarcut - small polys", |bencher| {
        let multi_poly = geo_test_fixtures::nl_zones::<f64>();
        bencher.iter(|| {
            for poly in &multi_poly {
                let triangulation = TriangulateEarcut::earcut_triangles(poly);
                assert!(triangulation.len() > 1);
                criterion::black_box(triangulation);
            }
        });
    });

    c.bench_function("TriangulateEarcut - large_poly", |bencher| {
        let poly = Polygon::new(geo_test_fixtures::norway_main::<f64>(), vec![]);
        bencher.iter(|| {
            let triangulation = TriangulateEarcut::earcut_triangles(&poly);
            assert!(triangulation.len() > 1);
            criterion::black_box(triangulation);
        });
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
