[![geo](https://avatars1.githubusercontent.com/u/10320338?v=4&s=50)](https://github.com/TimTheBig/geo-3d)

<!-- todo update -->
[![Coverage Status](https://img.shields.io/coverallsCoverage/github/georust/geo.svg)](https://coveralls.io/github/georust/geo?branch=trying)
[![Documentation](https://img.shields.io/docsrs/geo-3d/latest.svg)](https://docs.rs/geo-3d)

# geo-3D

## 3D Geospatial Primitives, Algorithms, and Utilities

The `geo-3d` crate provides geospatial primitive types such as `Point`, `LineString`, and `Polygon`, and provides algorithms and operations such as:
- Volume, area and centroid calculation
- Simplification and convex hull operations
- Euclidean distance measurement
- Intersection checks
- Affine transforms such as rotation and translation

Please refer to [the documentation](https://docs.rs/geo-3d) for a complete list.

## Example

<!-- todo update -->
```rust
// primitives
use geo::{line_string, polygon};

// algorithms
use geo::ConvexHull;

// An L shape
let poly = polygon![
    (x: 0.0, y: 0.0, z: 0.0),
    (x: 4.0, y: 0.0, z: 4.0),
    (x: 4.0, y: 1.0, z: 2.5),
    (x: 1.0, y: 1.0, z: 1.0),
    (x: 1.0, y: 4.0, z: 7.0),
    (x: 0.0, y: 4.0, z: 0.0),
    (x: 0.0, y: 0.0, z: 0.0),
];

// Calculate the polygon's convex hull
let hull = poly.convex_hull().unwrap();

assert_eq!(
    hull.exterior(),
    &line_string![
        (x: 4.0, y: 0.0, z: 4.0),
        (x: 4.0, y: 1.0, z: 2.5),
        (x: 1.0, y: 4.0, z: 7.0),
        (x: 0.0, y: 4.0, z: 0.0),
        (x: 0.0, y: 0.0, z: 0.0),
        (x: 4.0, y: 0.0, z: 4.0),
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
