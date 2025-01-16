[![geo](https://avatars1.githubusercontent.com/u/10320338?v=4&s=50)](https://github.com/georust)

<!-- todo update -->
[![Coverage Status](https://img.shields.io/coverallsCoverage/github/georust/geo.svg)](https://coveralls.io/github/georust/geo?branch=trying)
[![Documentation](https://img.shields.io/docsrs/geo/latest.svg)](https://docs.rs/geo)

# geo 3D

## 3D Geospatial Primitives, Algorithms, and Utilities

### TODO
- move all types to 3d
- make all methuds work with 3d
- fix doctests/tests
- valume calculations

The `geo` crate provides geospatial primitive types such as `Point`, `LineString`, and `Polygon`, and provides algorithms and operations such as:
- Area and centroid calculation
- Simplification and convex hull operations
- Euclidean and Haversine distance measurement
- Intersection checks
- Affine transforms such as rotation and translation
- All DE-9IM spatial predicates such as contains, crosses, and touches.

Please refer to [the documentation](https://docs.rs/geo-3d) for a complete list.

The primitive types also provide the basis for other functionality in the `Geo` ecosystem, including:

- [Coordinate transformation and projection](https://github.com/georust/proj)
- Serialization to and from [GeoJSON](https://github.com/georust/geojson) and [WKT](https://github.com/georust/wkt)

## Example

<!-- todo update -->
```rust
// primitives
use geo::{line_string, polygon};

// algorithms
use geo::ConvexHull;

// An L shape
let poly = polygon![
    (x: 0.0, y: 0.0),
    (x: 4.0, y: 0.0),
    (x: 4.0, y: 1.0),
    (x: 1.0, y: 1.0),
    (x: 1.0, y: 4.0),
    (x: 0.0, y: 4.0),
    (x: 0.0, y: 0.0),
];

// Calculate the polygon's convex hull
let hull = poly.convex_hull();

assert_eq!(
    hull.exterior(),
    &line_string![
        (x: 4.0, y: 0.0),
        (x: 4.0, y: 1.0),
        (x: 1.0, y: 4.0),
        (x: 0.0, y: 4.0),
        (x: 0.0, y: 0.0),
        (x: 4.0, y: 0.0),
    ]
);
```

## Contributing

Contributions are welcome! Have a look at the [issues](https://github.com/TimTheBig/geo-3d/issues), and open a pull request if you'd like to add an algorithm or some functionality.

## License

Licensed under either of

 * Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or https://www.apache.org/licenses/LICENSE-2.0)
 * MIT license ([LICENSE-MIT](LICENSE-MIT) or https://opensource.org/licenses/MIT)

at your option.
