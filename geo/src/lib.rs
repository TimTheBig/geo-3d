#![doc(html_logo_url = "https://raw.githubusercontent.com/georust/meta/master/logo/logo.png")]

//! The `geo-3d` crate provides geospatial primitive types and algorithms.
//!
//! # Types
//!
//! - **[`Coord`]**: A three-dimensional coordinate. All geometry types are composed of [`Coord`]s, though [`Coord`] itself is not a [`Geometry`] type
//! - **[`Point`]**: A single point represented by one [`Coord`]
//! - **[`MultiPoint`]**: A collection of [`Point`]s
//! - **[`Line`]**: A line segment represented by two [`Coord`]s
//! - **[`LineString`]**: A series of contiguous line segments represented by two or more [`Coord`]s
//! - **[`MultiLineString`]**: A collection of [`LineString`]s
//! - **[`Polygon`]**: A bounded area represented by one [`LineString`] exterior ring, and zero or
//!   more [`LineString`] interior rings
//! - **[`MultiPolygon`]**: A collection of [`Polygon`]s
//! - **[`Rect`]**: An axis-aligned bounded rectangle represented by minimum and maximum [`Coord`]s
//! - **[`Triangle`]**: A bounded area represented by three [`Coord`] vertices
//! - **[`GeometryCollection`]**: A collection of [`Geometry`]s
//! - **[`Geometry`]**: An enumeration of all geometry types, excluding [`Coord`]
//!
//! The preceding types are reexported from the [`geo-types`] crate. Consider using that crate
//! if you only need access to these types and no other `geo-3d` functionality.
//!
//! # Algorithms
//!
//! ## Measures
//!
//! Algorithms for measures along a line, and how a line is measured.
//!
//! ### Metric Spaces
//!
//! - **[`Euclidean`]**: The [Euclidean plane] measures distance with the pythagorean formula.
//!
//! ### Operations on Metric Spaces
//!
//! - **[`Distance`]**: Calculate the minimum distance between two geometries.
//! - **[`Length`]**: Calculate the length of a `Line`, `LineString`, or `MultiLineString`.
//! - **[`Bearing`]**: Calculate the bearing between two points.
//!
//! - **[`Destination`]**: Calculate the destination point from an origin point, given a bearing and a distance.
//! - **[`InterpolatePoint`]**: Interpolate points along a line.
//! - **[`Densify`]**: Insert points into a geometry so there is never more than `max_segment_length` between points.
//!
//! ### Misc measures
//!
//! - **[`HausdorffDistance`]**: Calculate "the maximum of the distances from a point in any of the sets to the nearest point in the other set." (Rote, 1991)
//! - **[`FrechetDistance`]**: Calculate the similarity between [`LineString`]s using the Fréchet distance
//!
//! ## Area/Volume
//!
//! - **[`Area`]**: Calculate the area of a geometry
//! - **[`Volume`]**: Calculate the volume of a geometry
//!
//! ## Outlier Detection
//!
//! - **[`OutlierDetection`]**: Detect outliers in a group of points using [LOF](https://en.wikipedia.org/wiki/Local_outlier_factor)
//!
//! ## Simplification
//!
//! - **[`Simplify`]**: Simplify a geometry using the Ramer–Douglas–Peucker algorithm
//! - **[`SimplifyIdx`]**: Calculate a simplified geometry using the Ramer–Douglas–Peucker algorithm, returning coordinate indices
//! - **[`SimplifyVw`]**: Simplify a geometry using the Visvalingam-Whyatt algorithm
//! - **[`SimplifyVwPreserve`]**: Simplify a geometry using a topology-preserving variant of the Visvalingam-Whyatt algorithm
//! - **[`SimplifyVwIdx`]**: Calculate a simplified geometry using the Visvalingam-Whyatt algorithm, returning coordinate indices
//!
//! ## Query
//!
//! - **[`ClosestPoint`]**: Find the point on a geometry
//!   closest to a given point
//! - **[`IsConvex`]**:
//!   Calculate the convexity of a [`LineString`]
//! - **[`LineInterpolatePoint`]**:
//!   Generates a point that lies a given fraction along the line
//! - **[`LineLocatePoint`]**: Calculate the
//!   fraction of a line’s total length representing the location of the closest point on the
//!   line to the given point
//! - **[`InteriorPoint`]**:
//!     Calculates a representative point inside a `Geometry`
//!
//! ## Topology
//!
//! - **[`Contains`]**: Calculate if a geometry contains another
//!   geometry
//! - **[`CoordinatePosition`]**: Calculate
//!   the position of a coordinate relative to a geometry
//! - **[`HasDimensions`]**: Determine the dimensions of a geometry
//! - **[`Intersects`]**: Calculate if a geometry intersects
//!   another geometry
//! - **[`line_intersection`]**: Calculates the
//!   intersection, if any, between two lines
//! - **[`Relate`]**: Topologically relate two geometries based on
//!   [DE-9IM](https://en.wikipedia.org/wiki/DE-9IM) semantics
//! - **[`Within`]**: Calculate if a geometry lies completely within another geometry
//!
//! ## Triangulation
//!
//! - **[`TriangulateEarcut`](triangulate_earcut)**: Triangulate polygons using the earcut algorithm. Requires the `"earcutr"` feature, which is enabled by default
//!
//! ## Winding
//!
//! - **[`Orient`]**: Apply a specified winding [`Direction`](orient::Direction) to a [`Polygon`]’s interior and exterior rings
//! - **[`Winding`]**: Calculate and manipulate the [`WindingOrder`](winding_order::WindingOrder) of a [`LineString`]
//!
//! ## Iteration
//!
//! - **[`CoordsIter`]**: Iterate over the coordinates of a geometry
//! - **[`MapCoords`]**: Map a function over all the coordinates
//!   in a geometry, returning a new geometry
//! - **[`MapCoordsInPlace`]**: Map a function over all the
//!   coordinates in a geometry in-place
//! - **[`LinesIter`]**: Iterate over lines of a geometry
//!
//! ## Boundary
//!
//! - **[`BoundingRect`]**: Calculate the axis-aligned
//!   bounding rectangle of a geometry
//! - **[`MinimumRotatedRect`]**: Calculate the
//!   minimum bounding box of a geometry
//! - **[`ConcaveHull`]**: Calculate the concave hull of a
//!   geometry
//! - **[`ConvexHull`]**: Calculate the convex hull of a
//!   geometry
//! - **[`Extremes`]**: Calculate the extreme coordinates and
//!   indices of a geometry
//!
//! ## Affine transformations
//!
//! - **[`Rotate`]**: Rotate a geometry around an axis or its centroid
//! - **[`Scale`]**: Scale a geometry up or down by a factor
//! - **[`Skew`]**: Skew a geometry by shearing angles along the `x`, `y`, and `z` dimension
//! - **[`Translate`]**: Translate a geometry along its axis
//! - **[`AffineOps`]**: generalised composable affine operations
//!
//! ## Conversion
//!
//! - **[`Convert`]**: Convert (infalliby) the numeric type of a geometry’s coordinate value
//! - **[`TryConvert`]**: Convert (falliby) the numeric type of a geometry’s coordinate value
//! - **[`ToDegrees`]**: Radians to degrees coordinate transforms for a given geometry
//! - **[`ToRadians`]**: Degrees to radians coordinate transforms for a given geometry
//!
//! ## Miscellaneous
//!
//! - **[`Centroid`]**: Calculate the centroid of a geometry
//! - **[`ChaikinSmoothing`]**: Smoothen `LineString`, `Polygon`, `MultiLineString` and `MultiPolygon` using Chaikin's algorithm
//! - **[`proj`]**: Project geometries with the `proj` crate (requires the `use-proj` feature)
//! - **[`LineStringSegmentize`]**: Segment a LineString into `n` segments
//! - **[`Transform`]**: Transform a geometry using Proj
//! - **[`RemoveRepeatedPoints`]**: Remove repeated points from a geometry
//! - **[`Validation`]**: Checks if the geometry is well formed. Some algorithms may not work correctly with invalid geometries
//! - **[`TriangulateEarcut`]**: Uses the [`earcutr`] crate, which provides triangulation of polygons using the earcut algorithm
//!
//! # Spatial Indexing
//!
//! `geo-3d` geometries ([`Point`], [`Line`], [`LineString`], [`Polygon`], [`MultiPolygon`]) can be used with the [rstar](https://docs.rs/rstar/0.12.0/rstar/struct.RTree.html#usage)
//! R*-tree crate for fast distance and nearest-neighbour queries. Multi- geometries can be added to the tree by iterating over
//! their members and adding them. Note in particular the availability of the [`bulk_load`](https://docs.rs/rstar/0.12.0/rstar/struct.RTree.html#method.bulk_load)
//! method and [`GeomWithData`](https://docs.rs/rstar/0.12.0/rstar/primitives/struct.GeomWithData.html) struct.
//!
//! # Features
//!
//! The following optional [Cargo features] are available:
//!
//! - `proj-network`:
//!     - Enables [network grid] support for the [`proj` crate]
//!     - After enabling this feature, [further configuration][proj crate file download] is required to use the network grid.
//!     - ☐ Disabled by default
//! - `use-proj`:
//!     - Enables coordinate conversion and transformation of `Point` geometries using the [`proj` crate]
//!     - ☐ Disabled by default
//! - `use-serde`:
//!     - Allows geometry types to be serialized and deserialized with [Serde]
//!     - ☐ Disabled by default
//! - `multithreading`:
//!     - Enables multithreading support (via Rayon), and activates the `multithreading` flag
//!       in `geo-types`, enabling multi-threaded iteration over `Multi*` geometries
//!     - ☑ Enabled by default
//!
//! [Euclidean plane]: https://en.wikipedia.org/wiki/Euclidean_plane
//! [`geo-types`]: https://crates.io/crates/geo-types
//! [`proj` crate]: https://github.com/georust/proj
//! [Karney (2013)]:  https://arxiv.org/pdf/1109.4448.pdf
//! [shapefile crate]: https://crates.io/crates/shapefile
//! [georust website]: https://georust.org
//! [Cargo features]: https://doc.rust-lang.org/cargo/reference/features.html
//! [network grid]: https://proj.org/usage/network.html
//! [proj crate file download]: https://docs.rs/proj/*/proj/#grid-file-download
//! [Serde]: https://serde.rs/

#![deny(unused)]

#[cfg(feature = "use-serde")]
#[macro_use]
extern crate serde;

pub use crate::algorithm::*;
pub use crate::types::Closest;
use std::cmp::Ordering;

pub use crate::relate::PreparedGeometry;
pub use geo_types::{coord, line_string, point, polygon, wkt, CoordNum};

pub mod geometry;
pub use geometry::*;

/// This module includes all the functions of geometric calculations
pub mod algorithm;
mod geometry_cow;
mod types;
mod utils;
use crate::kernels::RobustKernel;
pub(crate) use geometry_cow::GeometryCow;

#[cfg(test)]
#[macro_use]
extern crate approx;

#[macro_use]
extern crate log;

/// Radius of Earth at the equator in meters (derived from the WGS-84 ellipsoid)
const EQUATORIAL_EARTH_RADIUS: f64 = 6_378_137.0;

/// Radius of Earth at the poles in meters (derived from the WGS-84 ellipsoid)
const POLAR_EARTH_RADIUS: f64 = 6_356_752.314_245;

/// Flattening of the WGS-84 ellipsoid - [wikipedia: Flattening](https://en.wikipedia.org/wiki/Flattening)
const EARTH_FLATTENING: f64 =
    (EQUATORIAL_EARTH_RADIUS - POLAR_EARTH_RADIUS) / EQUATORIAL_EARTH_RADIUS;

/// A prelude which re-exports the traits for manipulating objects in this
/// crate. Typically imported with `use geo::prelude::*`.
pub mod prelude {
    pub use crate::algorithm::*;
}

/// A common numeric trait used for geo algorithms\
/// A trait for methods which work for floating point, using a Kernel
///
/// Different numeric types have different tradeoffs. `geo-3d` strives to utilize generics to allow
/// users to choose their numeric types.
///
/// # Examples
///
/// ```
/// use geo_3d::{GeoNum, MultiPolygon, Polygon, Point};
///
/// // An admittedly silly method implementation, but the signature shows how to use the GeoNum trait
/// fn farthest_from<'a, T: GeoNum>(point: &Point<T>, polygons: &'a MultiPolygon<T>) -> Option<&'a Polygon<T>> {
///     polygons.iter().fold(None, |accum, next| {
///         match accum {
///             None => Some(next),
///             Some(farthest) => {
///                 use geo_3d::Distance;
///                 if next.distance(point) > farthest.distance(point) {
///                     Some(next)
///                 } else {
///                     Some(farthest)
///                 }
///             }
///         }
///     })
/// }
/// ```
pub trait GeoNum: CoordNum + num_traits::Float + num_traits::Signed + num_traits::Bounded + float_next_after::NextAfter {
    type Ker: Kernel<Self>;

    /// Return the ordering between self and other.
    ///
    /// For integers, this should behave just like [`Ord`].
    ///
    /// For floating point numbers, unlike the standard partial comparison between floating point numbers, this comparison
    /// always produces an ordering.
    ///
    /// See [f64::total_cmp](https://doc.rust-lang.org/src/core/num/f64.rs.html#1432) for details.
    fn total_cmp(&self, other: &Self) -> Ordering;
}

macro_rules! impl_geo_num_for_float {
    ($t: ident) => {
        impl GeoNum for $t {
            type Ker = RobustKernel;
            fn total_cmp(&self, other: &Self) -> Ordering {
                self.total_cmp(other)
            }
        }
    };
}

// This is the list of primitives that we support.
impl_geo_num_for_float!(f32);
impl_geo_num_for_float!(f64);
// todo for when f128 gets stabilized
// impl_geo_num_for_float!(f128);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn total_ord_f64() {
        assert_eq!(GeoNum::total_cmp(&3.0f64, &2.0f64), Ordering::Greater);
        assert_eq!(GeoNum::total_cmp(&2.0f64, &2.0f64), Ordering::Equal);
        assert_eq!(GeoNum::total_cmp(&1.0f64, &2.0f64), Ordering::Less);
        assert_eq!(GeoNum::total_cmp(&1.0f64, &f64::NAN), Ordering::Less);
        assert_eq!(GeoNum::total_cmp(&f64::NAN, &f64::NAN), Ordering::Equal);
        assert_eq!(GeoNum::total_cmp(&f64::INFINITY, &f64::NAN), Ordering::Less);
    }

    #[test]
    fn total_ord_f32() {
        assert_eq!(GeoNum::total_cmp(&3.0f32, &2.0f32), Ordering::Greater);
        assert_eq!(GeoNum::total_cmp(&2.0f32, &2.0f32), Ordering::Equal);
        assert_eq!(GeoNum::total_cmp(&1.0f32, &2.0f32), Ordering::Less);
        assert_eq!(GeoNum::total_cmp(&1.0f32, &f32::NAN), Ordering::Less);
        assert_eq!(GeoNum::total_cmp(&f32::NAN, &f32::NAN), Ordering::Equal);
        assert_eq!(GeoNum::total_cmp(&f32::INFINITY, &f32::NAN), Ordering::Less);
    }

    #[test]
    fn numeric_types() {
        let _n_f32 = Point::new(1.0f32, 2.0f32, 3.0f32);
        let _n_f64 = Point::new(1.0f64, 2.0f64, 3.0f64);
        // todo for when f128 gets stabilized
        // let _n_f128 = Point::new(1.0f128, 2.0f128, 3.0f128);
    }
}
