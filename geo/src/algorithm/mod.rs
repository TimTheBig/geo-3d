/// Kernels to compute various predicates
pub mod kernels;
pub use kernels::{Kernel, Orientation};

/// Calculate the area of the surface of a `Geometry`.
pub mod area;
pub use area::Area;

/// Calculate the bounding rectangle of a `Geometry`.
pub mod bounding_rect;
pub use bounding_rect::BoundingRect;

/// Calculate the minimum rotated rectangle of a `Geometry`.
pub mod minimum_rotated_rect;
pub use minimum_rotated_rect::MinimumRotatedRect;

/// Calculate the centroid of a `Geometry`.
pub mod centroid;
pub use centroid::Centroid;

/// Smoothen `LineString`, `Polygon`, `MultiLineString` and `MultiPolygon` using Chaikins algorithm.
pub mod chaikin_smoothing;
pub use chaikin_smoothing::ChaikinSmoothing;

/// Calculate the closest `Point` between a `Geometry` and an input `Point`.
pub mod closest_point;
pub use closest_point::ClosestPoint;

/// Calculate the concave hull of a `Geometry`.
pub mod concave_hull;
pub use concave_hull::ConcaveHull;

/// Determine whether `Geometry` `A` completely encloses `Geometry` `B`.
pub mod contains;
pub use contains::Contains;

/// Convert the type of a geometry’s coordinate value.
pub mod convert;
pub use convert::{Convert, TryConvert};

/// Convert coordinate angle units between radians and degrees.
pub mod convert_angle_unit;
pub use convert_angle_unit::{ToDegrees, ToRadians};

/// Calculate the convex hull of a `Geometry`.
pub mod convex_hull;
pub use convex_hull::ConvexHull;

/// Determine whether a `Coord` lies inside, outside, or on the boundary of a geometry.
pub mod coordinate_position;
pub use coordinate_position::CoordinatePosition;

/// Iterate over geometry coordinates.
pub mod coords_iter;
pub use coords_iter::CoordsIter;

/// Dimensionality of a geometry and its boundary, based on OGC-SFA.
pub mod dimensions;
pub use dimensions::HasDimensions;

/// Calculate the extreme coordinates and indices of a geometry.
pub mod extremes;
pub use extremes::Extremes;

/// Calculate the Frechet distance between two `LineStrings`.
pub mod frechet_distance;
pub use frechet_distance::FrechetDistance;

/// Calculate the Hausdorff distance between two geometries.
pub mod hausdorff_distance;
pub use hausdorff_distance::HausdorffDistance;

/// Calculate a representative `Point` inside a `Geometry`
pub mod interior_point;
pub use interior_point::InteriorPoint;

/// Determine whether `Geometry` `A` intersects `Geometry` `B`.
pub mod intersects;
pub use intersects::Intersects;

/// Determines whether a `LineString` is convex.
pub mod is_convex;
pub use is_convex::IsConvex;

/// Calculate concave hull using k-nearest algorithm
pub mod k_nearest_concave_hull;
pub use k_nearest_concave_hull::KNearestConcaveHull;

/// Interpolate a point along a `Line` or `LineString`.
pub mod line_interpolate_point;
pub use line_interpolate_point::LineInterpolatePoint;

/// Computes the intersection of two Lines.
pub mod line_intersection;
pub use line_intersection::LineIntersection;

/// Locate a point along a `Line` or `LineString`.
pub mod line_locate_point;
pub use line_locate_point::LineLocatePoint;

/// Iterate over the lines in a geometry.
pub mod lines_iter;
pub use lines_iter::LinesIter;

pub mod line_measures;
pub use line_measures::Euclidean;
pub use line_measures::{Bearing, Densify, Destination, Distance, InterpolatePoint, Length};

/// Split a LineString into n segments
pub mod linestring_segment;
pub use linestring_segment::LineStringSegmentize;

/// Apply a function to all `Coord`s of a `Geometry`.
pub mod map_coords;
pub use map_coords::{MapCoords, MapCoordsInPlace};

/// Orient a `Polygon`'s exterior and interior rings.
pub mod orient;
pub use orient::Orient;

/// Coordinate projections and transformations using the current stable version of [PROJ](http://proj.org).
#[cfg(feature = "use-proj")]
pub use proj::*;

/// Relate two geometries based on DE-9IM
pub mod relate;
pub use relate::Relate;

/// Remove (consecutive) repeated points
pub mod remove_repeated_points;
pub use remove_repeated_points::RemoveRepeatedPoints;

/// Rotate a `Geometry` by an angle given in degrees.
pub mod rotate;
pub use rotate::Rotate;

/// Scale a `Geometry` up or down by a factor
pub mod scale;
pub use scale::Scale;

/// Composable affine operations such as rotate, scale, skew, and translate
pub mod affine_ops;
pub use affine_ops::{AffineOps, AffineTransform};

/// Simplify `Geometries` using the Ramer-Douglas-Peucker algorithm.
pub mod simplify;
pub use simplify::{Simplify, SimplifyIdx};

/// Simplify `Geometries` using the Visvalingam-Whyatt algorithm. Includes a topology-preserving variant.
pub mod simplify_vw;
pub use simplify_vw::{SimplifyVw, SimplifyVwIdx, SimplifyVwPreserve};

/// Transform a geometry using PROJ.
#[cfg(feature = "use-proj")]
pub mod transform;
#[cfg(feature = "use-proj")]
pub use transform::Transform;

/// Translate a `Geometry` along the given offsets.
pub mod translate;
pub use translate::Translate;

/// Triangulate polygons using an [ear-cutting algorithm](https://www.geometrictools.com/Documentation/TriangulationByEarClipping.pdf).
pub mod triangulate_delaunay;
pub use triangulate_delaunay::TriangulateDelaunay;

/// Vector Operations for 3D coordinates
mod vector_ops;
pub use vector_ops::Vector3DOps;

/// Calculate the Vincenty distance between two `Point`s.
pub mod vincenty_distance;
pub use vincenty_distance::VincentyDistance;

/// Calculate the Vincenty length of a `LineString`.
pub mod vincenty_length;
pub use vincenty_length::VincentyLength;

/// Calculate and work with the winding order of `Linestring`s.
pub mod winding_order;
pub use winding_order::Winding;

/// Determine whether `Geometry` `A` is completely within by `Geometry` `B`.
pub mod within;
pub use within::Within;

/// Planar sweep algorithm and related utils
pub mod sweep;

/// Detect outliers in a group of points using [LOF](https://en.wikipedia.org/wiki/Local_outlier_factor)
pub mod outlier_detection;

pub use outlier_detection::OutlierDetection;

/// Monotonic polygon subdivision
pub mod monotone;
pub use monotone::{monotone_subdivision, MonoPoly, MonotonicPolygons};

pub mod validation;
pub use validation::Validation;

/// Volume of geometry
pub mod volume;
pub use volume::Volume;

/// Project geometry on to a 2d plane
pub mod plane_project;
pub use plane_project::ProjectToPlane;
