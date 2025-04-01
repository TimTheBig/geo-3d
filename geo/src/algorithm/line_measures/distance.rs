/// Calculate the minimum distance between two geometries.
pub trait Distance<F, Destination> {
    /// Note that not all implementations support all geometry combinations, but at least `Point` to `Point`
    /// is supported.
    /// See [specific implementations](#implementors) for details.
    ///
    /// # Units
    ///
    /// - `origin`, `destination`: geometry where the units of x/y/z depend on the trait implementation.
    /// - returns: depends on the trait implementation.
    fn distance(self, destination: Destination) -> F;
}

use crate::algorithm::Intersects;
use crate::coordinate_position::CoordPos;
use crate::{geometry::*, CoordinatePosition};
use crate::{CoordNum, GeoNum};
use num_traits::{Bounded, Float};
use rstar::primitives::CachedEnvelope;
use rstar::RTree;

// Distance is a symmetric operation, so we can implement it once for both
macro_rules! symmetric_distance_impl {
    ($t:ident, $a:ty, $b:ty) => {
        impl<F> $crate::Distance<F, $b> for $a
        where
            F: $t,
        {
            fn distance(self, b: $b) -> F {
                b.distance(self)
            }
        }
    };
}

// ┌───────────────────────────┐
// │ Implementations for Coord │
// └───────────────────────────┘

impl<F: CoordNum> Distance<F, Coord<F>> for Coord<F> {
    // todo check z
    fn distance(self, destination: Coord<F>) -> F {
        let delta = self - destination;
        delta.x.hypot(delta.y).hypot(delta.z)
    }
}
impl<F: CoordNum> Distance<F, &Line<F>> for Coord<F> {
    fn distance(self, line: &Line<F>) -> F {
        ::geo_types::private_utils::point_line_euclidean_distance(Point(self), *line)
    }
}
impl<F: CoordNum> Distance<F, Line<F>> for Point<F> {
    fn distance(self, line: Line<F>) -> F {
        ::geo_types::private_utils::point_line_euclidean_distance(self, line)
    }
}

// ┌───────────────────────────┐
// │ Implementations for Point │
// └───────────────────────────┘

/// Calculate the Euclidean distance (a.k.a. pythagorean distance) between two Points
impl<F: CoordNum> Distance<F, Point<F>> for Point<F> {
    /// Calculate the Euclidean distance (a.k.a. pythagorean distance) between two Points
    ///
    /// # Units
    /// - `origin`, `destination`: Point where the units of x/y/z represent non-angular units
    ///    — e.g. meters or miles
    /// - returns: distance in the same units as the `origin` and `destination` points
    ///
    /// # Example
    /// ```
    /// use geo_3d::{Euclidean, Distance};
    /// use geo_3d::Point;
    /// // web mercator
    /// let new_york_city = Point::new(-8238310.24, 4942194.78);
    /// // web mercator
    /// let london = Point::new(-14226.63, 6678077.70);
    /// let distance: f64 = new_york_city.distance(london);
    ///
    /// assert_eq!(
    ///     8_405_286., // meters in web mercator
    ///     distance.round()
    /// );
    /// ```
    ///
    /// [`Haversine`]: crate::line_measures::metric_spaces::Haversine
    /// [`Geodesic`]: crate::line_measures::metric_spaces::Geodesic
    /// [metric spaces]: crate::line_measures::metric_spaces
    fn distance(self, destination: Point<F>) -> F {
        self.0.distance(destination.0)
    }
}

impl<F: CoordNum> Distance<F, &Point<F>> for &Point<F> {
    fn distance(self, destination: &Point<F>) -> F {
        (*self).distance(*destination)
    }
}

impl<F: CoordNum> Distance<F, &Line<F>> for &Point<F> {
    fn distance(self, destination: &Line<F>) -> F {
        geo_types::private_utils::point_line_euclidean_distance(*self, *destination)
    }
}

impl<F: CoordNum> Distance<F, &LineString<F>> for &Point<F> {
    fn distance(self, destination: &LineString<F>) -> F {
        geo_types::private_utils::point_line_string_euclidean_distance(*self, destination)
    }
}

impl<F: GeoNum> Distance<F, &Polygon<F>> for &Point<F> {
    fn distance(self, polygon: &Polygon<F>) -> F {
        // No need to continue if the polygon intersects the point, or is zero-length
        if polygon.exterior().0.is_empty() || polygon.intersects(self) {
            return F::zero();
        }
        // fold the minimum interior ring distance if any, followed by the exterior
        // shell distance, returning the minimum of the two distances
        polygon
            .interiors()
            .iter()
            .map(|ring| Self::distance(self, ring))
            .fold(Bounded::max_value(), |accum: F, val| accum.min(val))
            .min(
                polygon
                    .exterior()
                    .lines()
                    .map(|line| {
                        ::geo_types::private_utils::line_segment_distance(
                            self.0, line.start, line.end,
                        )
                    })
                    .fold(Bounded::max_value(), |accum, val| accum.min(val)),
            )
    }
}

// ┌──────────────────────────┐
// │ Implementations for Line │
// └──────────────────────────┘

symmetric_distance_impl!(CoordNum, &Line<F>, Coord<F>);
symmetric_distance_impl!(CoordNum, Line<F>, Point<F>);
symmetric_distance_impl!(CoordNum, &Line<F>, &Point<F>);

impl<F: GeoNum> Distance<F, &Line<F>> for &Line<F> {
    fn distance(self, line_b: &Line<F>) -> F {
        if self.intersects(line_b) {
            return F::zero();
        }
        // minimum of all Point-Line distances
        (&self.start_point()).distance(line_b)
            .min((&self.end_point()).distance(line_b))
            .min((&line_b.start_point()).distance(self))
            .min((&line_b.end_point()).distance(self))
    }
}

impl<F: GeoNum> Distance<F, &LineString<F>> for &Line<F> {
    fn distance(self, line_string: &LineString<F>) -> F {
        line_string
            .lines()
            .fold(Bounded::max_value(), |acc, segment| {
                acc.min(Self::distance(self, &segment))
            })
    }
}

impl<F: GeoNum> Distance<F, &Polygon<F>> for &Line<F> {
    fn distance(self, polygon: &Polygon<F>) -> F {
        if self.intersects(polygon) {
            return F::zero();
        }

        // REVIEW: This impl changed slightly.
        std::iter::once(polygon.exterior())
            .chain(polygon.interiors().iter())
            .fold(Bounded::max_value(), |acc, line_string| {
                acc.min(Self::distance(self, line_string))
            })
    }
}

// ┌────────────────────────────────┐
// │ Implementations for LineString │
// └────────────────────────────────┘

symmetric_distance_impl!(CoordNum, &LineString<F>, &Point<F>);
symmetric_distance_impl!(GeoNum, &LineString<F>, &Line<F>);

impl<F: GeoNum> Distance<F, &LineString<F>> for &LineString<F> {
    fn distance(self, line_string_b: &LineString<F>) -> F {
        if self.intersects(line_string_b) {
            F::zero()
        } else {
            nearest_neighbour_distance(self, line_string_b)
        }
    }
}

impl<F: GeoNum> Distance<F, &Polygon<F>> for &LineString<F> {
    fn distance(self, polygon: &Polygon<F>) -> F {
        if self.intersects(polygon) {
            F::zero()
        } else if !polygon.interiors().is_empty()
            // FIXME: Explodes on empty line_string
            && ring_contains_coord(polygon.exterior(), self.0[0])
        {
            // check each ring distance, returning the minimum
            let mut mindist: F = Float::max_value();
            for ring in polygon.interiors() {
                mindist = mindist.min(nearest_neighbour_distance(self, ring))
            }
            mindist
        } else {
            nearest_neighbour_distance(self, polygon.exterior())
        }
    }
}

// ┌─────────────────────────────┐
// │ Implementations for Polygon │
// └─────────────────────────────┘

symmetric_distance_impl!(GeoNum, &Polygon<F>, &Point<F>);
symmetric_distance_impl!(GeoNum, &Polygon<F>, &Line<F>);
symmetric_distance_impl!(GeoNum, &Polygon<F>, &LineString<F>);
impl<F: GeoNum> Distance<F, &Polygon<F>> for &Polygon<F> {
    fn distance(self, polygon_b: &Polygon<F>) -> F {
        if self.intersects(polygon_b) {
            return F::zero();
        }
        // FIXME: explodes when polygon_b.exterior() is empty
        // Containment check
        if !self.interiors().is_empty()
            && ring_contains_coord(self.exterior(), polygon_b.exterior().0[0])
        {
            // check each ring distance, returning the minimum
            let mut mindist: F = Float::max_value();
            for ring in self.interiors() {
                mindist = mindist.min(nearest_neighbour_distance(polygon_b.exterior(), ring))
            }
            return mindist;
        } else if !polygon_b.interiors().is_empty()
            // FIXME: explodes when self.exterior() is empty
            && ring_contains_coord(polygon_b.exterior(), self.exterior().0[0])
        {
            let mut mindist: F = Float::max_value();
            for ring in polygon_b.interiors() {
                mindist = mindist.min(nearest_neighbour_distance(self.exterior(), ring))
            }
            return mindist;
        }
        nearest_neighbour_distance(self.exterior(), polygon_b.exterior())
    }
}

// ┌────────────────────────────────────────┐
// │ Implementations for Rect and Triangle  │
// └────────────────────────────────────────┘

/// Implements Euclidean distance for Triangles and Rects by converting them to polygons.
macro_rules! impl_euclidean_distance_for_polygonlike_geometry {
  ($polygonlike:ty,  [$($geometry_b:ty),*]) => {
      impl<F: GeoNum> Distance<F, $polygonlike> for $polygonlike
      {
          fn distance(self, destination: $polygonlike) -> F {
              self.to_polygon().distance(destination)
          }
      }
      $(
          impl<F: GeoNum> Distance<F, $geometry_b> for $polygonlike
          {
              fn distance(self, geometry_b: $geometry_b) -> F {
                    self.to_polygon().distance(geometry_b)
              }
          }
          symmetric_distance_impl!(GeoNum, $geometry_b, $polygonlike);
      )*
  };
}

impl_euclidean_distance_for_polygonlike_geometry!(&Triangle<F>, [&Point<F>, &Line<F>, &LineString<F>, &Polygon<F>, &Rect<F>]);
impl_euclidean_distance_for_polygonlike_geometry!(&Rect<F>,     [&Point<F>, &Line<F>, &LineString<F>, &Polygon<F>]);

// ┌───────────────────────────────────────────┐
// │ Implementations for multi geometry types  │
// └───────────────────────────────────────────┘

/// Euclidean distance implementation for multi geometry types.
macro_rules! impl_euclidean_distance_for_iter_geometry {
    ($iter_geometry:ty,  [$($to_geometry:ty),*]) => {
        impl<F: GeoNum> Distance<F, $iter_geometry> for $iter_geometry {
            fn distance(self, destination: $iter_geometry) -> F {
                self
                    .iter()
                    .fold(Bounded::max_value(), |accum: F, member| {
                        accum.min(member.distance(destination))
                    })
             }
        }
        $(
            impl<F: GeoNum> Distance<F, $to_geometry> for $iter_geometry {
                fn distance(self, to_geometry: $to_geometry) -> F {
                    self
                        .iter()
                        .fold(Bounded::max_value(), |accum: F, member| {
                            accum.min(member.distance(to_geometry))
                        })
                }
            }
            symmetric_distance_impl!(GeoNum, $to_geometry, $iter_geometry);
        )*
  };
}

impl_euclidean_distance_for_iter_geometry!(&MultiPoint<F>,         [&Point<F>, &Line<F>, &LineString<F>, &MultiLineString<F>, &Polygon<F>, &MultiPolygon<F>, &GeometryCollection<F>, &Rect<F>, &Triangle<F>]);
impl_euclidean_distance_for_iter_geometry!(&MultiLineString<F>,    [&Point<F>, &Line<F>, &LineString<F>,                      &Polygon<F>, &MultiPolygon<F>, &GeometryCollection<F>, &Rect<F>, &Triangle<F>]);
impl_euclidean_distance_for_iter_geometry!(&MultiPolygon<F>,       [&Point<F>, &Line<F>, &LineString<F>,                      &Polygon<F>,                   &GeometryCollection<F>, &Rect<F>, &Triangle<F>]);
impl_euclidean_distance_for_iter_geometry!(&GeometryCollection<F>, [&Point<F>, &Line<F>, &LineString<F>,                      &Polygon<F>,                                           &Rect<F>, &Triangle<F>]);

// ┌──────────────────────────────┐
// │ Implementation for Geometry  │
// └──────────────────────────────┘

/// Euclidean distance implementation for every specific Geometry type to Geometry<T>.
macro_rules! impl_euclidean_distance_for_geometry_and_variant {
  ([$($target:ty),*]) => {
      $(
          impl<F: GeoNum> Distance<F, &Geometry<F>> for $target {
              fn distance(self, destination: &Geometry<F>) -> F {
                  match destination {
                      Geometry::Point(point) => Self::distance(self, point),
                      Geometry::Line(line) => Self::distance(self, line),
                      Geometry::LineString(line_string) => Self::distance(self, line_string),
                      Geometry::Polygon(polygon) => Self::distance(self, polygon),
                      Geometry::MultiPoint(multi_point) => Self::distance(self, multi_point),
                      Geometry::MultiLineString(multi_line_string) => Self::distance(self, multi_line_string),
                      Geometry::MultiPolygon(multi_polygon) => Self::distance(self, multi_polygon),
                      Geometry::GeometryCollection(geometry_collection) => Self::distance(self, geometry_collection),
                      Geometry::Rect(rect) => Self::distance(self, rect),
                      Geometry::Triangle(triangle) => Self::distance(self, triangle),
                  }
              }
          }
          symmetric_distance_impl!(GeoNum, &Geometry<F>, $target);
      )*
  };
}

impl_euclidean_distance_for_geometry_and_variant!([&Point<F>, &MultiPoint<F>, &Line<F>, &LineString<F>, &MultiLineString<F>, &Polygon<F>, &MultiPolygon<F>, &Triangle<F>, &Rect<F>, &GeometryCollection<F>]);

impl<F: GeoNum> Distance<F, &Geometry<F>> for &Geometry<F> {
    fn distance(self, destination: &Geometry<F>) -> F {
        match self {
            Geometry::Point(point) => point.distance(destination),
            Geometry::Line(line) => line.distance(destination),
            Geometry::LineString(line_string) => line_string.distance(destination),
            Geometry::Polygon(polygon) => polygon.distance(destination),
            Geometry::MultiPoint(multi_point) => multi_point.distance(destination),
            Geometry::MultiLineString(multi_line_string) => {
                multi_line_string.distance(destination)
            }
            Geometry::MultiPolygon(multi_polygon) => multi_polygon.distance(destination),
            Geometry::GeometryCollection(geometry_collection) => {
                geometry_collection.distance(destination)
            }
            Geometry::Rect(rect) => rect.distance(destination),
            Geometry::Triangle(triangle) => triangle.distance(destination),
        }
    }
}

// ┌───────────────────────────┐
// │ Implementations utilities │
// └───────────────────────────┘

/// Uses an R* tree and nearest-neighbour lookups to calculate minimum distances
// This is somewhat slow and memory-inefficient, but certainly better than quadratic time
fn nearest_neighbour_distance<F: GeoNum>(geom1: &LineString<F>, geom2: &LineString<F>) -> F {
    let tree_a = RTree::bulk_load(geom1.lines().map(CachedEnvelope::new).collect());
    let tree_b = RTree::bulk_load(geom2.lines().map(CachedEnvelope::new).collect());
    // Return minimum distance between all geom a points and geom b lines, and all geom b points and geom a lines
    geom2
        .points()
        .fold(Bounded::max_value(), |acc: F, point| {
            let nearest = tree_a.nearest_neighbor(&point).unwrap();
            acc.min((nearest as &Line<F>).distance(&point))
        })
        .min(geom1.points().fold(Bounded::max_value(), |acc, point| {
            let nearest = tree_b.nearest_neighbor(&point).unwrap();
            acc.min((nearest as &Line<F>).distance(&point))
        }))
}

fn ring_contains_coord<T: GeoNum>(ring: &LineString<T>, c: Coord<T>) -> bool {
    match ring.coordinate_position(&c) {
        CoordPos::Inside => true,
        CoordPos::OnBoundary | CoordPos::Outside => false,
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::orient::{Direction, Orient};
    use crate::{Line, LineString, MultiLineString, MultiPoint, MultiPolygon, Point, Polygon};
    use geo_types::{coord, polygon, private_utils::line_segment_distance};

    #[test]
    fn line_segment_distance_test() {
        let o1 = Point::new(8.0, 0.0, 0.0);
        let o2 = Point::new(5.5, 0.0, 0.0);
        let o3 = Point::new(5.0, 0.0, 0.0);
        let o4 = Point::new(4.5, 1.5, 0.0);

        let p1 = Point::new(7.2, 2.0, 0.0);
        let p2 = Point::new(6.0, 1.0, 0.0);

        let dist = line_segment_distance(o1, p1, p2);
        let dist2 = line_segment_distance(o2, p1, p2);
        let dist3 = line_segment_distance(o3, p1, p2);
        let dist4 = line_segment_distance(o4, p1, p2);
        // Results agree with Shapely
        assert_relative_eq!(dist, 2.0485900789263356);
        assert_relative_eq!(dist2, 1.118033988749895);
        assert_relative_eq!(dist3, std::f64::consts::SQRT_2);
        assert_relative_eq!(dist4, 1.5811388300841898);
        // Point is on the line
        let zero_dist = line_segment_distance(p1, p1, p2);
        assert_relative_eq!(zero_dist, 0.0);
    }

    #[test]
    // Point to Polygon, outside point
    fn point_polygon_distance_outside_test() {
        // an octagon
        let points = vec![
            (5., 1., 0.0),
            (4., 2., 0.0),
            (4., 3., 4.0),
            (5., 4., 0.0),
            (6., 4., 0.0),
            (7., 3., 0.0),
            (7., 2., 0.0),
            (6., 1., 0.0),
            (5., 1., 0.0),
        ];
        let ls = LineString::from(points);
        let poly = Polygon::new(ls, vec![]);
        // A Random point outside the octagon
        let p = Point::new(2.5, 0.5, 0.0);
        let dist = (&p).distance(&poly);
        assert_relative_eq!(dist, 3.1213203435596424);
    }

    #[test]
    // Point to Polygon, inside point
    fn point_polygon_distance_inside_test() {
        // an octagon
        let points = vec![
            (5., 1., 0.0),
            (4., 2., 0.0),
            (4., 3., 0.0),
            (5., 4., 0.0),
            (6., 4., 0.0),
            (7., 3., 0.0),
            (7., 2., 0.0),
            (6., 1., 0.0),
            (5., 1., 0.0),
        ];
        let ls = LineString::from(points);
        let poly = Polygon::new(ls, vec![]);
        // A Random point inside the octagon
        let p = Point::new(5.5, 2.1, 0.0);
        let dist = (&p).distance(&poly);
        assert_relative_eq!(dist, 0.0);
    }

    #[test]
    // Point to Polygon, on boundary
    fn point_polygon_distance_boundary_test() {
        // an octagon
        let points = vec![
            (5., 1., 0.0),
            (4., 2., 0.0),
            (4., 3., 0.0),
            (5., 4., 0.0),
            (6., 4., 0.0),
            (7., 3., 0.0),
            (7., 2., 0.0),
            (6., 1., 0.0),
            (5., 1., 0.0),
        ];
        let ls = LineString::from(points);
        let poly = Polygon::new(ls, vec![]);
        // A point on the octagon
        let p = Point::new(5.0, 1.0, 0.0);
        let dist = (&p).distance(&poly);
        assert_relative_eq!(dist, 0.0);
    }

    #[test]
    // Point to Polygon, on boundary (second variant)
    fn point_polygon_boundary_test2() {
        let exterior = LineString::from(vec![
            (0., 0., 0.0),
            (0., 0.0004, 0.0),
            (0.0004, 0.0004, 0.0),
            (0.0004, 0., 0.0),
            (0., 0., 0.0),
        ]);

        let poly = Polygon::new(exterior, vec![]);
        let bugged_point = Point::new(0.0001, 0., 0.);
        assert_relative_eq!((&poly).distance(&bugged_point), 0.);
    }

    #[test]
    // Point to Polygon, empty Polygon
    fn point_polygon_empty_test() {
        // an empty Polygon
        let points = vec![];
        let ls = LineString::new(points);
        let poly = Polygon::new(ls, vec![]);
        // A point outside the empty polygon
        let p = Point::new(2.5, 0.5, 0.3);
        let dist = (&p).distance(&poly);
        assert_relative_eq!(dist, 0.0);
    }

    #[test]
    // Point to Polygon with an interior ring
    fn point_polygon_interior_cutout_test() {
        // an octagon
        let ext_points = vec![
            (4., 1., 0.0),
            (5., 2., 0.0),
            (5., 3., 0.0),
            (4., 4., 0.0),
            (3., 4., 0.0),
            (2., 3., 0.0),
            (2., 2., 0.0),
            (3., 1., 0.0),
            (4., 1., 0.0),
        ];
        // cut out a triangle inside octagon
        let int_points = vec![
            (3.5, 3.5, 0.0),
            (4.4, 1.5, 0.0),
            (2.6, 1.5, 0.0),
            (3.5, 3.5, 0.0),
        ];
        let ls_ext = LineString::from(ext_points);
        let ls_int = LineString::from(int_points);
        let poly = Polygon::new(ls_ext, vec![ls_int]);
        // A point inside the cutout triangle
        let p = Point::new(3.5, 2.5, 0.0);
        let dist = (&p).distance(&poly);

        // 0.41036467732879783 <-- Shapely
        assert_relative_eq!(dist, 0.41036467732879767);
    }

    #[test]
    fn line_distance_multipolygon_do_not_intersect_test() {
        // checks that the distance from the multipolygon
        // is equal to the distance from the closest polygon
        // taken in isolation, whatever that distance is
        let ls1 = LineString::from(vec![
            (0.0, 0.0, 0.0),
            (10.0, 0.0, 0.0),
            (10.0, 10.0, 0.0),
            (5.0, 15.0, 0.0),
            (0.0, 10.0, 0.0),
            (0.0, 0.0, 0.0),
        ]);
        let ls2 = LineString::from(vec![
            (0.0, 30.0, 0.4),
            (0.0, 25.0, 0.0),
            (10.0, 25.0, 0.0),
            (10.0, 30.0, 0.0),
            (0.0, 30.0, 0.0),
        ]);
        let ls3 = LineString::from(vec![
            (15.0, 30.0, 0.0),
            (15.0, 25.0, 0.0),
            (20.0, 25.0, 6.0),
            (20.0, 30.0, 0.0),
            (15.0, 30.0, 0.0),
        ]);
        let pol1 = Polygon::new(ls1, vec![]);
        let pol2 = Polygon::new(ls2, vec![]);
        let pol3 = Polygon::new(ls3, vec![]);
        let mp = MultiPolygon::new(vec![pol1.clone(), pol2, pol3]);
        let pnt1 = Point::new(0.0, 15.0, 0.0);
        let pnt2 = Point::new(10.0, 20.0, 0.0);
        let ln = Line::new(pnt1.0, pnt2.0);
        let dist_mp_ln = (&ln).distance(&mp);
        let dist_pol1_ln = (&ln).distance(&pol1);
        assert_relative_eq!(dist_mp_ln, dist_pol1_ln);
    }

    #[test]
    fn point_distance_multipolygon_test() {
        let ls1 = LineString::from(vec![
            (0.0, 0.0, 0.0),
            (1.0, 10.0, 0.0),
            (2.0, 0.0, 0.0),
            (0.0, 0.0, 0.0),
        ]);
        let ls2 = LineString::from(vec![
            (3.0, 0.0, 0.0),
            (4.0, 10.0, 0.0),
            (5.0, 0.0, 0.0),
            (3.0, 0.0, 0.0),
        ]);
        let p1 = Polygon::new(ls1, vec![]);
        let p2 = Polygon::new(ls2, vec![]);
        let mp = MultiPolygon::new(vec![p1, p2]);
        let p = Point::new(50.0, 50.0, 0.0);
        assert_relative_eq!((&p).distance(&mp), 60.959002616512684);
    }

    #[test]
    // Point to LineString
    fn point_linestring_distance_test() {
        // like an octagon, but missing the lowest horizontal segment
        let points = vec![
            (5., 1., 7.),
            (4., 2., 0.0),
            (4., 3., 0.0),
            (5., 4., 0.0),
            (6., 4., 0.0),
            (7., 3., 0.0),
            (7., 2., 0.0),
            (6., 1., 0.0),
        ];
        let ls = LineString::from(points);
        // A Random point "inside" the LineString
        let p = Point::new(5.5, 2.1, 0.0);
        let dist = (&p).distance(&ls);
        assert_relative_eq!(dist, 1.1313708498984762);
    }

    #[test]
    // Point to LineString, point lies on the LineString
    fn point_linestring_contains_test() {
        let points = vec![
            (5., 1., 0.0),
            (4., 2., 0.0),
            (4., 3., 2.),
            (5., 4., 3.),
            (6., 4., 0.0),
            (7., 3., 0.0),
            (7., 2., 0.0),
            (6., 1., 0.0),
        ];
        let ls = LineString::from(points);
        let p = Point::new(5.0, 4.0, 0.0);
        let dist = (&p).distance(&ls);
        assert_relative_eq!(dist, 0.0);
    }

    #[test]
    // Point to LineString, closed triangle
    fn point_linestring_triangle_test() {
        let points = vec![
            (3.5, 3.5, 3.5),
            (4.4, 2.0, 0.0),
            (2.6, 2.0, 0.0),
            (3.5, 3.5, 3.5),
        ];
        let ls = LineString::from(points);
        let p = Point::new(3.5, 2.5, 0.5);
        let dist = (&p).distance(&ls);
        assert_relative_eq!(dist, 0.5);
    }

    #[test]
    // Point to LineString, empty LineString
    fn point_linestring_empty_test() {
        let points = vec![];
        let ls = LineString::new(points);
        let p = Point::new(5.0, 4.0, 0.0);
        let dist = (&p).distance(&ls);
        assert_relative_eq!(dist, 0.0);
    }

    #[test]
    fn distance_multilinestring_test() {
        let v1 = LineString::from(vec![
            (0.0, 0.0, 0.0),
            (1.0, 10.0, 0.0),
        ]);
        let v2 = LineString::from(vec![
            (1.0, 10.0, 0.0),
            (2.0, 0.0, 0.0),
            (3.0, 1.0, 0.0),
        ]);
        let mls = MultiLineString::new(vec![v1, v2]);
        let p = Point::new(50.0, 50.0, 0.0);
        assert_relative_eq!((&p).distance(&mls), 63.25345840347388);
    }

    #[test]
    fn distance1_test() {
        assert_relative_eq!(
            (&Point::new(0., 0., 0.)).distance(&Point::new(1., 0., 0.)),
            1.
        );
    }

    #[test]
    fn distance2_test() {
        let dist = (&Point::new(-72.1235, 42.3521, 0.0))
            .distance(&Point::new(72.1260, 70.612, 0.0));
        assert_relative_eq!(dist, 146.99163308930207);
    }

    #[test]
    fn distance_multipoint_test() {
        let v = vec![
            Point::new(0.0, 10.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
            Point::new(10.0, 0.0, 0.0),
            Point::new(1.0, -1.0, 0.0),
            Point::new(0.0, -10.0, 0.0),
            Point::new(-1.0, -1.0, 0.0),
            Point::new(-10.0, 0.0, 0.0),
            Point::new(-1.0, 1.0, 0.0),
            Point::new(0.0, 10.0, 0.0),
        ];
        let mp = MultiPoint::new(v);
        let p = Point::new(50.0, 50.0, 50.);
        assert_relative_eq!((&p).distance(&mp), 64.03124237432849)
    }

    #[test]
    fn distance_line_test() {
        let line0 = Line::from([(0., 0., 0.), (5., 0., 5.)]);
        let p0 = Point::new(2., 3., 4.);
        let p1 = Point::new(3., 0., 3.);
        let p2 = Point::new(6., 0., 6.);
        assert_relative_eq!((&line0).distance(&p0), 3.);
        assert_relative_eq!((&p0).distance(&line0), 3.);

        assert_relative_eq!((&line0).distance(&p1), 0.);
        assert_relative_eq!((&p1).distance(&line0), 0.);

        assert_relative_eq!((&line0).distance(&p2), 1.);
        assert_relative_eq!((&p2).distance(&line0), 1.);
    }

    #[test]
    fn distance_line_line_test() {
        let line0 = Line::from([(0., 0., 0.), (5., 0., 5.)]);
        let line1 = Line::from([(2., 1., 0.), (7., 2., 0.)]);
        assert_relative_eq!((&line0).distance(&line1), 1.);
        assert_relative_eq!((&line1).distance(&line0), 1.);
    }

    #[test]
    // See https://github.com/georust/geo/issues/476
    fn distance_line_polygon_test() {
        let line = Line::new(
            coord! {
                x: -0.17084137691985102,
                y: 0.8748085493016657,
                z: 0.0,
            },
            coord! {
                x: -0.17084137691985102,
                y: 0.09858870312437906,
                z: 3.1213203435596424,
            },
        );
        let poly: Polygon<f64> = polygon![
            coord! {
                x: -0.10781391405721802,
                y: -0.15433610862574643,
                z: 0.0,
            },
            coord! {
                x: -0.7855276236615211,
                y: 0.23694208404779793,
                z: 0.0,
            },
            coord! {
                x: -0.7855276236615214,
                y: -0.5456143012992907,
                z: 0.0,
            },
            coord! {
                x: -0.10781391405721802,
                y: -0.15433610862574643,
                z: 0.0,
            },
        ];
        assert_eq!((&line).distance(&poly), 0.18752558079168907);
    }

    #[test]
    // test edge-vertex minimum distance
    fn test_minimum_polygon_distance() {
        let points_raw = [
            (126., 232., 0.0),
            (126., 212., 0.0),
            (112., 202., 0.0),
            (97., 204., 0.0),
            (87., 215., 0.0),
            (87., 232., 0.0),
            (100., 246., 0.0),
            (118., 247., 0.0),
        ];
        let points = points_raw
            .iter()
            .map(|e| Point::new(e.0, e.1, e.2))
            .collect::<Vec<_>>();
        let poly1 = Polygon::new(LineString::from(points), vec![]);

        let points_raw_2 = [
            (188., 231., 0.),
            (189., 207., 0.),
            (174., 196., 0.),
            (164., 196., 0.),
            (147., 220., 0.),
            (158., 242., 0.),
            (177., 242., 0.),
        ];
        let points2 = points_raw_2
            .iter()
            .map(|e| Point::new(e.0, e.1, e.2))
            .collect::<Vec<_>>();
        let poly2 = Polygon::new(LineString::from(points2), vec![]);
        let dist = nearest_neighbour_distance(poly1.exterior(), poly2.exterior());
        assert_relative_eq!(dist, 21.0);
    }

    #[test]
    // test vertex-vertex minimum distance
    fn test_minimum_polygon_distance_2() {
        let points_raw = [
            (118., 200., 2.),
            (153., 179., 10.),
            (106., 155., 20.),
            (88., 190., 30.),
            (118., 200., 37.),
        ];
        let points = points_raw
            .iter()
            .map(|e| Point::new(e.0, e.1, e.2))
            .collect::<Vec<_>>();
        let poly1 = Polygon::new(LineString::from(points), vec![]);

        let points_raw_2 = [
            (242., 186., 0.),
            (260., 146., 157.),
            (182., 175., 0.),
            (216., 193., 236.),
            (242., 186., 0.),
        ];
        let points2 = points_raw_2
            .iter()
            .map(|e| Point::new(e.0, e.1, e.2))
            .collect::<Vec<_>>();
        let poly2 = Polygon::new(LineString::from(points2), vec![]);
        let dist = nearest_neighbour_distance(poly1.exterior(), poly2.exterior());
        assert_relative_eq!(dist, 29.274562336608895);
    }

    #[test]
    // test edge-edge minimum distance
    fn test_minimum_polygon_distance_3() {
        let points_raw = [
            (182., 182., 0.0),
            (182., 168., 0.0),
            (138., 160., 0.0),
            (136., 193., 0.0),
            (182., 182., 0.0),
        ];
        let points = points_raw
            .iter()
            .map(|e| Point::new(e.0, e.1, e.2))
            .collect::<Vec<_>>();
        let poly1 = Polygon::new(LineString::from(points), vec![]);

        let points_raw_2 = [
            (232., 196., 0.0),
            (234., 150., 0.0),
            (194., 165., 0.0),
            (194., 191., 0.0),
            (232., 196., 0.0),
        ];
        let points2 = points_raw_2
            .iter()
            .map(|e| Point::new(e.0, e.1, e.2))
            .collect::<Vec<_>>();
        let poly2 = Polygon::new(LineString::from(points2), vec![]);
        let dist = nearest_neighbour_distance(poly1.exterior(), poly2.exterior());
        assert_relative_eq!(dist, 12.0);
    }

    #[test]
    fn test_large_polygon_distance() {
        let ls = geo_test_fixtures::norway_main::<f64>();
        let poly1 = Polygon::new(ls, vec![]);
        let vec2 = vec![
            (4.921875, 66.33750501996518, 31.213206424),
            (3.69140625, 65.21989393613207, 32.32096424),
            (6.15234375, 65.07213008560697, 0.0),
            (4.921875, 66.33750501996518, 31.213206424),
        ];
        let poly2 = Polygon::new(vec2.into(), vec![]);
        let distance = (&poly1).distance(&poly2);
        // GEOS says 2.2864896295566055
        assert_relative_eq!(distance, 2.2864896295566055);
    }

    #[test]
    // A polygon inside another polygon's ring; they're disjoint in the DE-9IM sense:
    // FF2FF1212
    fn test_poly_in_ring() {
        let shell = geo_test_fixtures::shell::<f64>();
        let ring = geo_test_fixtures::ring::<f64>();
        let poly_in_ring = geo_test_fixtures::poly_in_ring::<f64>();
        // inside is "inside" outside's ring, but they are disjoint
        let outside = Polygon::new(shell, vec![ring]);
        let inside = Polygon::new(poly_in_ring, vec![]);
        assert_relative_eq!((&outside).distance(&inside), 5.992772737231033);
    }

    #[test]
    // two ring LineStrings; one encloses the other but they neither touch nor intersect
    fn test_linestring_distance() {
        let ring = geo_test_fixtures::ring::<f64>();
        let poly_in_ring = geo_test_fixtures::poly_in_ring::<f64>();
        assert_relative_eq!((&ring).distance(&poly_in_ring), 5.992772737231033);
    }

    #[test]
    // Line-Polygon test: closest point on Polygon is NOT nearest to a Line end-point
    fn test_line_polygon_simple() {
        let line = Line::from([(0.0, 0.0, 0.0), (0.0, 3.0, 0.0)]);
        let v = vec![
            (5.0, 1.0, 0.0),
            (5.0, 2.0, 0.0),
            (0.25, 1.5, 0.0),
            (5.0, 1.0, 0.0),
        ];
        let poly = Polygon::new(v.into(), vec![]);
        assert_relative_eq!((&line).distance(&poly), 0.25);
    }

    #[test]
    // Line-Polygon test: Line intersects Polygon
    fn test_line_polygon_intersects() {
        let line = Line::from([(0.5, 0.0, 0.0), (0.0, 3.0, 0.0)]);
        let v = vec![
            (5.0, 1.0, 0.0),
            (5.0, 2.0, 6.0),
            (0.25, 1.5, 0.5),
            (5.0, 1.0, 0.0),
        ];
        let poly = Polygon::new(v.into(), vec![]);
        assert_relative_eq!((&line).distance(&poly), 0.0);
    }

    #[test]
    // Line-Polygon test: Line contained by interior ring
    fn test_line_polygon_inside_ring() {
        let line = Line::from([(4.4, 1.5, 0.0), (4.45, 1.5, 0.0)]);
        let v = vec![
            (5.0, 1.0, 0.0),
            (5.0, 2.0, 0.0),
            (0.25, 1.0, 0.0),
            (5.0, 1.0, 0.0),
        ];
        let v2 = vec![
            (4.5, 1.2, 0.0),
            (4.5, 1.8, 0.0),
            (3.5, 1.2, 0.0),
            (4.5, 1.2, 0.0),
        ];
        let poly = Polygon::new(v.into(), vec![v2.into()]);
        assert_relative_eq!((&line).distance(&poly), 0.04999999999999982);
    }

    #[test]
    // LineString-Line test
    fn test_linestring_line_distance() {
        let line = Line::from([(0.0, 0.0, 0.0), (0.0, 2.0, 0.0)]);
        let ls: LineString<_> = vec![(3.0, 0.0, 0.0), (1.0, 1.0, 0.0), (3.0, 2.0, 0.0)].into();
        assert_relative_eq!((&ls).distance(&line), 1.0);
    }

    #[test]
    // Triangle-Point test: point on vertex
    fn test_triangle_point_on_vertex_distance() {
        let triangle = Triangle::from([(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (2.0, 2.0, 0.0)]);
        let point = Point::new(0.0, 0.0, 0.0);
        assert_relative_eq!((&triangle).distance(&point), 0.0);
    }

    #[test]
    // Triangle-Point test: point on edge
    fn test_triangle_point_on_edge_distance() {
        let triangle = Triangle::from([(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (2.0, 2.0, 0.0)]);
        let point = Point::new(1.5, 0.0, 0.0);
        assert_relative_eq!((&triangle).distance(&point), 0.0);
    }

    #[test]
    // Triangle-Point test
    fn test_triangle_point_distance() {
        let triangle = Triangle::from([(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (2.0, 2.0, 0.0)]);
        let point = Point::new(2.0, 3.0, 0.0);
        assert_relative_eq!((&triangle).distance(&point), 1.0);
    }

    #[test]
    // Triangle-Point test: point within triangle
    fn test_triangle_point_inside_distance() {
        let triangle = Triangle::from([(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (2.0, 2.0, 0.0)]);
        let point = Point::new(1.0, 0.5, 0.0);
        assert_relative_eq!((&triangle).distance(&point), 0.0);
    }

    #[test]
    fn convex_and_nearest_neighbour_comparison() {
        let ls1: LineString<f64> = vec![
            Coord::from((57.39453770777941, 307.60533608924663, 0.0)),
            Coord::from((67.1800355576469, 309.6654408997451, 37.1213203435596424)),
            Coord::from((84.89693692793338, 225.5101593908847, 0.0)),
            Coord::from((75.1114390780659, 223.45005458038628, 0.0)),
            Coord::from((57.39453770777941, 307.60533608924663, 0.0)),
        ]
        .into();
        let first_polygon: Polygon<f64> = Polygon::new(ls1, vec![]);
        let ls2: LineString<f64> = vec![
            Coord::from((138.11769866645008, -45.75134112915392, 0.0)),
            Coord::from((130.50230476949187, -39.270154833870336, 0.0)),
            Coord::from((184.94426964987397, 24.699153900578573, 0.0)),
            Coord::from((192.55966354683218, 18.217967605294987, 0.0)),
            Coord::from((138.11769866645008, -45.75134112915392, 0.0)),
        ]
        .into();
        let second_polygon = Polygon::new(ls2, vec![]);

        assert_relative_eq!(
            (&first_polygon).distance(&second_polygon),
            224.35357967013238
        );
    }

    #[test]
    fn fast_path_regression() {
        // this test will fail if the fast path algorithm is reintroduced without being fixed
        let p1 = polygon!(
            (x: 0_f64, y: 0_f64, z: 0.0),
            (x: 300_f64, y: 0_f64, z: 0.0),
            (x: 300_f64, y: 100_f64, z: 100.),
            (x: 0_f64, y: 100_f64, z: 0.0),
        )
        .orient(Direction::Default);
        let p2 = polygon!(
            (x: 100_f64, y: 150_f64, z: 0.0),
            (x: 150_f64, y: 200_f64, z: 0.0),
            (x: 50_f64, y: 200_f64, z: 0.0),
        )
        .orient(Direction::Default);
        let p3 = polygon!(
            (x: 0_f64, y: 0_f64, z: 0.0),
            (x: 300_f64, y: 0_f64, z: 0.0),
            (x: 300_f64, y: 100_f64, z: 0.0),
            (x: 0_f64, y: 100_f64, z: 0.0),
        )
        .orient(Direction::Reversed);
        let p4 = polygon!(
            (x: 100_f64, y: 150_f64, z: 0.0),
            (x: 150_f64, y: 200_f64, z: 0.0),
            (x: 50_f64, y: 200_f64, z: 0.0),
        )
        .orient(Direction::Reversed);
        assert_eq!((&p1).distance(&p2), 50.0f64);
        assert_eq!((&p3).distance(&p4), 50.0f64);
        assert_eq!((&p1).distance(&p4), 50.0f64);
        assert_eq!((&p2).distance(&p3), 50.0f64);
    }

    #[test]
    fn all_types_geometry_collection_test() {
        let p = Point::new(0.0, 0.0, 0.0);
        let line = Line::from([(-1.0, -1.0, 0.0), (-2.0, -2.0, 0.0)]);
        let ls = LineString::from(vec![(0.0, 0.0, 0.0), (1.0, 10.0, 0.0), (2.0, 0.0, 0.0)]);
        let poly = Polygon::new(
            LineString::from(vec![(0.0, 0.0, 0.0), (1.0, 10.0, 0.0), (2.0, 0.0, 0.0), (0.0, 0.0, 0.0)]),
            vec![],
        );
        let tri = Triangle::from([(0.0, 0.0, 0.0), (1.0, 10.0, 0.0), (2.0, 0.0, 0.0)]);
        let rect = Rect::new((0.0, 0.0, 0.0), (-1.0, -1.0, 0.0));

        let ls1 = LineString::from(vec![(0.0, 0.0, 0.0), (1.0, 10.0, 0.0), (2.0, 0.0, 0.0), (0.0, 0.0, 0.0)]);
        let ls2 = LineString::from(vec![(3.0, 0.0, 0.0), (4.0, 10.0, 0.0), (5.0, 0.0, 0.0), (3.0, 0.0, 0.0)]);
        let p1 = Polygon::new(ls1, vec![]);
        let p2 = Polygon::new(ls2, vec![]);
        let mpoly = MultiPolygon::new(vec![p1, p2]);

        let v = vec![
            Point::new(0.0, 10.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
            Point::new(10.0, 0.0, 0.0),
            Point::new(1.0, -1.0, 0.0),
            Point::new(0.0, -10.0, 0.0),
            Point::new(-1.0, -1.0, 0.0),
            Point::new(-10.0, 0.0, 0.0),
            Point::new(-1.0, 1.0, 0.0),
            Point::new(0.0, 10.0, 0.0),
        ];
        let mpoint = MultiPoint::new(v);

        let v1 = LineString::from(vec![(0.0, 0.0, 0.0), (1.0, 10.0, 0.0)]);
        let v2 = LineString::from(vec![(1.0, 10.0, 0.0), (2.0, 0.0, 0.0), (3.0, 1.0, 0.0)]);
        let mls = MultiLineString::new(vec![v1, v2]);

        let gc = GeometryCollection(vec![
            Geometry::Point(p),
            Geometry::Line(line),
            Geometry::LineString(ls),
            Geometry::Polygon(poly),
            Geometry::MultiPoint(mpoint),
            Geometry::MultiLineString(mls),
            Geometry::MultiPolygon(mpoly),
            Geometry::Triangle(tri),
            Geometry::Rect(rect),
        ]);

        let test_p = Point::new(50., 50., 50.);
        assert_relative_eq!((&test_p).distance(&gc), 60.959002616512684);

        let test_multipoint = MultiPoint::new(vec![test_p]);
        assert_relative_eq!(
            (&test_multipoint).distance(&gc),
            60.959002616512684
        );

        let test_line = Line::from([(50., 50., 50.), (60., 60., 60.)]);
        assert_relative_eq!((&test_line).distance(&gc), 60.959002616512684);

        let test_ls = LineString::from(vec![(50., 50., 50.), (60., 60., 60.), (70., 70., 70.)]);
        assert_relative_eq!((&test_ls).distance(&gc), 60.959002616512684);

        let test_mls = MultiLineString::new(vec![test_ls]);
        assert_relative_eq!((&test_mls).distance(&gc), 60.959002616512684);

        let test_poly = Polygon::new(
            LineString::from(vec![
                (50., 50., 50.),
                (60., 50., 60.),
                (60., 60., 60.),
                (55., 55., 55.),
                (50., 50., 50.),
            ]),
            vec![],
        );
        assert_relative_eq!((&test_poly).distance(&gc), 60.959002616512684);

        let test_multipoly = MultiPolygon::new(vec![test_poly]);
        assert_relative_eq!(
            (&test_multipoly).distance(&gc),
            60.959002616512684
        );

        let test_tri = Triangle::from([(50., 50., 50.), (60., 50., 40.), (55., 55., 55.)]);
        assert_relative_eq!((&test_tri).distance(&gc), 60.959002616512684);

        let test_rect = Rect::new(coord! { x: 50., y: 50., z: 50. }, coord! { x: 60., y: 60., z: 60. });
        assert_relative_eq!((&test_rect).distance(&gc), 60.959002616512684);

        let test_gc = GeometryCollection(vec![Geometry::Rect(test_rect)]);
        assert_relative_eq!((&test_gc).distance(&gc), 60.959002616512684);
    }
}
