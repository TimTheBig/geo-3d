use criterion::{criterion_group, criterion_main, Criterion};
use geo_3d::intersects::Intersects;
use geo_3d::MultiPolygon;

fn multi_polygon_intersection(c: &mut Criterion) {
    let plot_polygons: MultiPolygon = geo_test_fixtures::nl_plots_wgs84();
    let zone_polygons: MultiPolygon = geo_test_fixtures::nl_zones();

    c.bench_function("MultiPolygon intersects", |bencher| {
        bencher.iter(|| {
            let mut intersects = 0;
            let mut non_intersects = 0;

            for a in &plot_polygons {
                for b in &zone_polygons {
                    if criterion::black_box(a.intersects(b)) {
                        intersects += 1;
                    } else {
                        non_intersects += 1;
                    }
                }
            }

            assert_eq!(intersects, 974);
            assert_eq!(non_intersects, 27782);
        });
    });
}

fn rect_intersection(c: &mut Criterion) {
    use geo_3d::algorithm::BoundingRect;
    use geo_3d::geometry::Rect;
    let plot_bbox: Vec<Rect> = geo_test_fixtures::nl_plots_wgs84()
        .iter()
        .map(|plot| plot.bounding_rect().unwrap())
        .collect();
    let zone_bbox: Vec<Rect> = geo_test_fixtures::nl_zones()
        .iter()
        .map(|plot| plot.bounding_rect().unwrap())
        .collect();

    c.bench_function("Rect intersects", |bencher| {
        bencher.iter(|| {
            let mut intersects = 0;
            let mut non_intersects = 0;

            for a in &plot_bbox {
                for b in &zone_bbox {
                    if criterion::black_box(a.intersects(b)) {
                        intersects += 1;
                    } else {
                        non_intersects += 1;
                    }
                }
            }

            assert_eq!(intersects, 3054);
            assert_eq!(non_intersects, 25702);
        });
    });
}

fn point_rect_intersection(c: &mut Criterion) {
    use geo_3d::algorithm::{BoundingRect, Centroid};
    use geo_3d::geometry::{Point, Rect};
    let plot_centroids: Vec<Point> = geo_test_fixtures::nl_plots_wgs84()
        .iter()
        .map(|plot| plot.centroid().unwrap())
        .collect();
    let zone_bbox: Vec<Rect> = geo_test_fixtures::nl_zones()
        .iter()
        .map(|plot| plot.bounding_rect().unwrap())
        .collect();

    c.bench_function("Point intersects rect", |bencher| {
        bencher.iter(|| {
            let mut intersects = 0;
            let mut non_intersects = 0;

            for a in &plot_centroids {
                for b in &zone_bbox {
                    if criterion::black_box(a.intersects(b)) {
                        intersects += 1;
                    } else {
                        non_intersects += 1;
                    }
                }
            }

            assert_eq!(intersects, 2246);
            assert_eq!(non_intersects, 26510);
        });
    });
}

fn point_triangle_intersection(c: &mut Criterion) {
    use geo_3d::{Centroid, Triangulate};
    use geo_types::{Point, Triangle};
    let plot_centroids: Vec<Point> = geo_test_fixtures::nl_plots_wgs84()
        .iter()
        .map(|plot| plot.centroid().unwrap())
        .collect();
    let zone_triangles: Vec<Triangle> = geo_test_fixtures::nl_zones()
        .iter()
        .flat_map(|plot| plot.triangles_iter())
        .collect();

    c.bench_function("Point intersects triangle", |bencher| {
        bencher.iter(|| {
            let mut intersects = 0;
            let mut non_intersects = 0;

            for a in &plot_centroids {
                for b in &zone_triangles {
                    if criterion::black_box(a.intersects(b)) {
                        intersects += 1;
                    } else {
                        non_intersects += 1;
                    }
                }
            }

            assert_eq!(intersects, 533);
            assert_eq!(non_intersects, 5450151);
        });
    });

    c.bench_function("Triangle intersects point", |bencher| {
        let triangle = Triangle::from([(0., 0., 0.), (10., 0., 10.), (5., 10., 5.)]);
        let point = Point::new(5., 5., 5.);

        bencher.iter(|| {
            assert!(criterion::black_box(&triangle).intersects(criterion::black_box(&point)));
        });
    });

    c.bench_function("Triangle intersects point on edge", |bencher| {
        let triangle = Triangle::from([(0., 0., 0.), (10., 0., 10.), (6., 10., 6.)]);
        let point = Point::new(3., 5., 3.);

        bencher.iter(|| {
            assert!(criterion::black_box(&triangle).intersects(criterion::black_box(&point)));
        });
    });
}

criterion_group! {
    name = bench_multi_polygons;
    config = Criterion::default().sample_size(10);
    targets = multi_polygon_intersection
}
criterion_group!(bench_rects, rect_intersection);
criterion_group! {
    name = bench_point_rect;
    config = Criterion::default().sample_size(50);
    targets = point_rect_intersection
}
criterion_group! {
    name = bench_point_triangle;
    config = Criterion::default().sample_size(50);
    targets = point_triangle_intersection
}

criterion_main!(
    bench_multi_polygons,
    bench_rects,
    bench_point_rect,
    bench_point_triangle
);
