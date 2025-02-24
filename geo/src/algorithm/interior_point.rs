use std::cmp::{Ordering, Reverse};

use crate::algorithm::{
    bounding_rect::BoundingRect,
    centroid::Centroid,
    coords_iter::CoordsIter,
    dimensions::HasDimensions,
    line_intersection::LineIntersection,
    line_measures::Distance,
    lines_iter::LinesIter,
    relate::Relate,
};
use crate::geometry::*;
use crate::sweep::{Intersections, SweepPoint};
use crate::GeoFloat;

/// Calculation of interior points.
///
/// An interior point is a point that's guaranteed to intersect a given geometry, and will be
/// strictly on the interior of the geometry if possible, or on the edge if the geometry has zero
/// area. A best effort will additionally be made to locate the point reasonably centrally.
///
/// For polygons, this point is located by drawing a line that approximately subdivides the
/// bounding box around the polygon in half, intersecting it with the polygon, then calculating
/// the midpoint of the longest line produced by the intersection. For lines, the non-endpoint
/// vertex closest to the line's centroid is returned if the line has interior points, or an
/// endpoint is returned otherwise.
///
/// For multi-geometries or collections, the interior points of the constituent components are
/// calculated, and one of those is returned (for MultiPolygons, it's the point that's the midpoint
/// of the longest intersection of the intersection lines of any of the constituent polygons, as
/// described above; for all others, the interior point closest to the collection's centroid is
/// used).
///
/// # Examples
///
/// ```
/// use geo::InteriorPoint;
/// use geo::{point, polygon};
///
/// // rhombus shaped polygon
/// let polygon = polygon![
///     (x: -2., y: 1., z: -2.),
///     (x: 1., y: 3., z: 1.),
///     (x: 4., y: 1., z: 4.),
///     (x: 1., y: -1., z: 1.),
///     (x: -2., y: 1., z: -2.),
/// ];
///
/// assert_eq!(
///     Some(point!(x: 1., y: 2., z: 1.)),
///     polygon.interior_point(),
/// );
/// ```
pub trait InteriorPoint {
    type Output;

    /// Calculates a representative point inside the `Geometry`
    ///
    /// # Examples
    ///
    /// ```
    /// use geo::InteriorPoint;
    /// use geo::{line_string, point};
    ///
    /// let line_string = line_string![
    ///     (x: 40.02f64, y: 116.34),
    ///     (x: 40.02f64, y: 118.23),
    ///     (x: 40.02f64, y: 120.15),
    /// ];
    ///
    /// assert_eq!(
    ///     Some(point!(x: 40.02, y: 118.23, z: 40.02)),
    ///     line_string.interior_point(),
    /// );
    /// ```
    fn interior_point(&self) -> Self::Output;
}

impl<T: GeoFloat> InteriorPoint for Line<T> {
    type Output = Point<T>;

    fn interior_point(&self) -> Self::Output {
        // the midpoint of the line isn't guaranteed to actually have an `intersects()`
        // relationship with the line due to floating point rounding, so just use the start point
        self.start_point()
    }
}

impl<T: GeoFloat> InteriorPoint for LineString<T> {
    type Output = Option<Point<T>>;

    // The interior point of a LineString the non-endpoint vertex closest to the centroid if any,
    // or the start point if there are no non-endpoint vertices
    fn interior_point(&self) -> Self::Output {
        match self.0.len() {
            0 => None,
            // for linestrings of length 2, as with lines, the calculated midpoint might not lie
            // on the line, so just use the start point
            1 | 2 => Some(self.0[0].into()),
            _ => {
                let centroid = self
                    .centroid()
                    .expect("expected centroid for non-empty linestring");
                self.0[1..(self.0.len() - 1)]
                    .iter()
                    .map(|coord| {
                        let pt = Point::from(*coord);
                        (pt, pt.distance(centroid))
                    })
                    .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Less))
                    .map(|(pt, _distance)| pt)
            }
        }
    }
}

impl<T: GeoFloat> InteriorPoint for MultiLineString<T> {
    type Output = Option<Point<T>>;

    /// The interior point of a MultiLineString is, of the interior points of all the constituent
    /// LineStrings, the one closest to the centroid of the MultiLineString
    fn interior_point(&self) -> Self::Output {
        if let Some(centroid) = self.centroid() {
            self.iter()
                .filter_map(|linestring| {
                    linestring
                        .interior_point()
                        .map(|pt| (pt, pt.distance(centroid)))
                })
                .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Less))
                .map(|(pt, _distance)| pt)
        } else {
            None
        }
    }
}

fn polygon_interior_point_with_segment_length<T: GeoFloat>(polygon: &Polygon<T>) -> Option<(Point<T>, T)> {
    // special-case a one-point polygon since this algorithm won't otherwise support it
    if polygon.exterior().0.len() == 1 {
        return Some((polygon.exterior().0[0].into(), T::zero()));
    }

    let two = T::one() + T::one();

    let bounds = polygon.bounding_rect()?;

    // use the midpoint of the bounds to scan, unless that happens to match any vertices from
    // polygon; if it does, perturb the line a bit by averaging with the Y coordinate of the
    // next-closest-to-center vertex if possible, to reduce the likelihood of collinear
    // intersections
    let mut y_mid = (bounds.min().y + bounds.max().y) / two;
    let mut z_mid = (bounds.min().z + bounds.max().z) / two;

    if polygon.coords_iter().any(|coord| coord.y == y_mid) {
        let next_closest = polygon
            .coords_iter()
            .filter_map(|coord| {
                if coord.y == y_mid {
                    None
                } else {
                    Some((coord.y, (coord.y - y_mid).abs()))
                }
            })
            .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Less));
        if let Some((closest, _)) = next_closest {
            y_mid = (y_mid + closest) / two;
            z_mid = (z_mid + closest) / two;
        }
    };

    let scan_line = Line::new(
        Coord {
            x: bounds.min().x,
            y: y_mid,
            z: z_mid,
        },
        Coord {
            x: bounds.max().x,
            y: y_mid,
            z: z_mid,
        },
    );

    let lines = polygon.lines_iter().chain(std::iter::once(scan_line));

    let mut intersections: Vec<SweepPoint<T>> = Vec::new();
    for (l1, l2, inter) in Intersections::from_iter(lines) {
        if !(l1 == scan_line || l2 == scan_line) {
            continue;
        }
        match inter {
            LineIntersection::Collinear { intersection } => {
                intersections.push(SweepPoint::from(intersection.start));
                intersections.push(SweepPoint::from(intersection.end));
            }
            LineIntersection::SinglePoint { intersection, .. } => {
                intersections.push(SweepPoint::from(intersection));
            }
        }
    }
    intersections.sort();

    let mut segments = Vec::new();
    let mut intersections_iter = intersections.iter().peekable();
    while let (Some(start), Some(end)) = (intersections_iter.next(), intersections_iter.peek()) {
        let length = end.x - start.x;
        let midpoint = Point::new((start.x + end.x) / two, y_mid, z_mid);
        segments.push((midpoint, length));
    }
    segments.sort_by(|a, b| b.1.total_cmp(&a.1));

    for (midpoint, segment_length) in segments {
        // some pairs of consecutive points traveling east-west will bound segments inside the
        // polygon, and some outside; confirm that this is the former
        let relation = polygon.relate(&midpoint);
        if relation.is_intersects() {
            return Some((
                midpoint,
                if relation.is_contains() {
                    segment_length
                } else {
                    // if our point is on the boundary, it must be because this is a zero-area
                    // polygon, so if we're being called from a multipolygon context, we want this
                    // option to be down-ranked as compared to other polygons that might have
                    // non-zero area
                    T::zero()
                },
            ));
        }
    }
    // if we've gotten this far with no luck, return any vertex point, if there are any
    polygon
        .coords_iter()
        .next()
        .map(|coord| (coord.into(), T::zero()))
}

impl<T: GeoFloat> InteriorPoint for Polygon<T> {
    type Output = Option<Point<T>>;

    fn interior_point(&self) -> Self::Output {
        polygon_interior_point_with_segment_length(self).map(|(point, _length)| point)
    }
}

impl<T: GeoFloat> InteriorPoint for MultiPolygon<T> {
    type Output = Option<Point<T>>;

    fn interior_point(&self) -> Self::Output {
        let segments = self
            .iter()
            .filter_map(polygon_interior_point_with_segment_length);
        segments
            .min_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(Ordering::Less))
            .map(|(point, _length)| point)
    }
}

impl<T: GeoFloat> InteriorPoint for Rect<T> {
    type Output = Point<T>;

    fn interior_point(&self) -> Self::Output {
        self.center().into()
    }
}

impl<T: GeoFloat> InteriorPoint for Point<T> {
    type Output = Point<T>;

    fn interior_point(&self) -> Self::Output {
        *self
    }
}

///
/// ```
/// use geo::InteriorPoint;
/// use geo::{MultiPoint, Point};
///
/// let empty: Vec<Point> = Vec::new();
/// let empty_multi_points: MultiPoint<_> = empty.into();
/// assert_eq!(empty_multi_points.interior_point(), None);
///
/// let points: MultiPoint<_> = vec![(5., 1., 5.), (1., 3., 1.), (3., 2., 3.)].into();
/// assert_eq!(points.interior_point(), Some(Point::new(3., 2., 3.)));
/// ```
impl<T: GeoFloat> InteriorPoint for MultiPoint<T> {
    type Output = Option<Point<T>>;

    fn interior_point(&self) -> Self::Output {
        if let Some(centroid) = self.centroid() {
            self.iter()
                .map(|pt| (pt, pt.distance(&centroid)))
                .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Less))
                .map(|(pt, _distance)| *pt)
        } else {
            None
        }
    }
}

impl<T: GeoFloat> InteriorPoint for Geometry<T> {
    type Output = Option<Point<T>>;

    crate::geometry_delegate_impl! {
        fn interior_point(&self) -> Self::Output;
    }
}

impl<T: GeoFloat> InteriorPoint for GeometryCollection<T> {
    type Output = Option<Point<T>>;

    fn interior_point(&self) -> Self::Output {
        if let Some(centroid) = self.centroid() {
            self.iter()
                .filter_map(|geom| {
                    geom.interior_point().map(|pt| {
                        (
                            pt,
                            // maximize dimensions, minimize distance
                            (
                                Reverse(geom.dimensions()),
                                pt.distance(centroid),
                            ),
                        )
                    })
                })
                .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Less))
                .map(|(pt, _distance)| pt)
        } else {
            None
        }
    }
}

impl<T: GeoFloat> InteriorPoint for Triangle<T> {
    type Output = Point<T>;

    fn interior_point(&self) -> Self::Output {
        self.centroid()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{
        algorithm::{contains::Contains, intersects::Intersects},
        coord, line_string, point, polygon,
    };

    // Tests: InteriorPoint of LineString
    #[test]
    fn empty_linestring_test() {
        let linestring: LineString<f32> = line_string![];
        let interior_point = linestring.interior_point();
        assert!(interior_point.is_none());
    }
    #[test]
    fn linestring_one_point_test() {
        let coord = coord! {
            x: 40.02f64,
            y: 116.34,
            z: 8.907,
        };
        let linestring = line_string![coord];
        let interior_point = linestring.interior_point();
        assert_eq!(interior_point, Some(Point::from(coord)));
    }
    #[test]
    fn linestring_test() {
        let linestring = line_string![
            (x: 1., y: 1., z: 1.),
            (x: 7., y: 1., z: 1.),
            (x: 8., y: 1., z: 1.),
            (x: 9., y: 1., z: 1.),
            (x: 10., y: 1., z: 1.),
            (x: 11., y: 1., z: 1.)
        ];
        assert_eq!(linestring.interior_point(), Some(point!(x: 7., y: 1., z: 1.)));
    }
    #[test]
    fn linestring_with_repeated_point_test() {
        let l1 = LineString::from(vec![
            point!(1., 1., 1.),
            point!(1., 1., 1.),
            point!(1., 1., 1.),
        ]);
        assert_eq!(l1.interior_point(), Some(point!(1., 1., 1.)));

        let l2 = LineString::from(vec![
            point!(2., 2., 2.),
            point!(2., 2., 2.),
            point!(2., 2., 2.),
        ]);
        let mls = MultiLineString::new(vec![l1, l2]);
        assert_eq!(mls.interior_point(), Some(point!(1., 1., 1.)));
    }
    // Tests: InteriorPoint of MultiLineString
    #[test]
    fn empty_multilinestring_test() {
        let mls: MultiLineString = MultiLineString::new(vec![]);
        let interior_point = mls.interior_point();
        assert!(interior_point.is_none());
    }
    #[test]
    fn multilinestring_with_empty_line_test() {
        let mls: MultiLineString = MultiLineString::new(vec![line_string![]]);
        let interior_point = mls.interior_point();
        assert!(interior_point.is_none());
    }
    #[test]
    fn multilinestring_length_0_test() {
        let coord = coord! {
            x: 40.02f64,
            y: 116.34,
            z: 9.83267,
        };
        let mls: MultiLineString = MultiLineString::new(vec![
            line_string![coord],
            line_string![coord],
            line_string![coord],
        ]);
        assert_relative_eq!(mls.interior_point().unwrap(), Point::from(coord));
    }
    #[test]
    fn multilinestring_one_line_test() {
        let linestring = line_string![
            (x: 1., y: 1., z: 1.),
            (x: 7., y: 1., z: 1.),
            (x: 8., y: 1., z: 1.),
            (x: 9., y: 1., z: 1.4),
            (x: 10., y: 1., z: 1.),
            (x: 11., y: 1., z: 1.)
        ];
        let mls: MultiLineString = MultiLineString::new(vec![linestring]);
        assert_relative_eq!(mls.interior_point().unwrap(), point! { x: 7., y: 1., z: 1. });
    }
    #[test]
    fn multilinestring_test() {
        let v1 = line_string![
            (x: 0.0, y: 0.0, z: 0.0),
            (x: 1.0, y: 10.0, z: 100.0)
        ];
        let v2 = line_string![
            (x: 1.0, y: 10.0, z: 100.0),
            (x: 2.0, y: 0.0, z: 0.0),
            (x: 3.0, y: 1.0, z: 0.0)
        ];
        let v3 = line_string![
            (x: -12.0, y: -100.0, z: 0.0),
            (x: 7.0, y: 8.0, z: 80.0)
        ];
        let mls = MultiLineString::new(vec![v1, v2, v3]);
        assert_eq!(mls.interior_point().unwrap(), point![x: 0., y: 0., z: 0.]);
    }
    // Tests: InteriorPoint of Polygon
    #[test]
    fn empty_polygon_test() {
        let poly: Polygon<f32> = polygon![];
        assert!(poly.interior_point().is_none());
    }
    #[test]
    fn polygon_one_point_test() {
        let p = point![ x: 2., y: 1., z: 0. ];
        let v = Vec::new();
        let linestring = line_string![p.0];
        let poly = Polygon::new(linestring, v);
        assert_relative_eq!(poly.interior_point().unwrap(), p);
    }

    #[test]
    fn polygon_test() {
        let poly = polygon![
            (x: 0., y: 0., z: 0.),
            (x: 2., y: 0., z: 2.),
            (x: 2., y: 2., z: 2.),
            (x: 0., y: 2., z: 0.),
            (x: 0., y: 0., z: 0.)
        ];
        assert_relative_eq!(poly.interior_point().unwrap(), point![x: 1., y: 1., z: 1.]);
    }
    #[test]
    fn polygon_hole_test() {
        // hexagon
        let ls1 = LineString::from(vec![
            (5.0, 1.0, 4.0),
            (4.0, 2.0, 0.0),
            (4.0, 3.0, 0.0),
            (5.0, 4.0, 9.0),
            (6.0, 4.0, 0.0),
            (7.0, 3.0, 0.5),
            (7.0, 2.0, 0.0),
            (6.0, 1.0, 9.0),
            (5.0, 1.0, 0.0),
        ]);

        let ls2 = LineString::from(vec![
            (5.0, 1.3, 0.0),
            (5.5, 2.0, 0.0),
            (6.0, 1.3, 3.0),
            (5.0, 1.3, 0.0)
        ]);

        let ls3 = LineString::from(vec![
            (5.0, 2.3, 0.0),
            (5.5, 3.0, 0.0),
            (6.0, 2.3, 0.0),
            (5.0, 2.3, 0.0)
        ]);

        let p1 = Polygon::new(ls1, vec![ls2, ls3]);
        let interior_point = p1.interior_point().unwrap();
        assert!(p1.contains(&interior_point));
        assert_relative_eq!(interior_point, point!(x: 4.571428571428571, y: 2.5, z: 0.0));
    }
    #[test]
    fn flat_polygon_test() {
        let poly = Polygon::new(
            LineString::from(vec![
                point!(0., 1., 2.),
                point!(1., 1., 1.),
                point!(0., 1., 2.)
            ]),
            vec![],
        );
        assert_eq!(poly.interior_point(), Some(point!(0.5, 1., 1.5)));
    }
    #[test]
    fn diagonal_flat_polygon_test() {
        // the regular intersection approach happens to not produce a point that intersects the
        // polygon given these particular start values, so this tests falling back to a vertex
        let start: Coord<f64> = Coord {
            x: 0.632690318327692,
            y: 0.08104532928154995,
            z: 0.8659032644203238,
        };
        let end: Coord<f64> = Coord {
            x: 0.4685039949468325,
            y: 0.31750332644855794,
            z: 0.0,
        };
        let poly = Polygon::new(LineString::new(vec![start, end, start]), vec![]);

        assert_eq!(poly.interior_point(), Some(start.into()));
    }
    #[test]
    fn polygon_vertex_on_median() {
        let poly = Polygon::new(
            LineString::from(vec![
                (0.5, 1.0, 0.5),
                (0.5, 0.5, 0.5),
                (0.0, 0.5, 0.0),
                (0.0, 0.0, 0.0),
                (1.0, 0.0, 1.0),
                (1.0, 1.0, 1.0),
                (0.5, 1.0, 0.5),
            ]),
            vec![],
        );
        let interior_point = poly.interior_point().unwrap();
        assert_eq!(&interior_point, &point!(0.75, 0.75, 0.75));
    }
    #[test]
    fn multi_poly_with_flat_polygon_test() {
        let poly = Polygon::new(
            LineString::from(vec![
                point!(0., 0., 0.),
                point!(1., 0., 1.),
                point!(0., 0., 0.)
            ]),
            vec![],
        );
        let multipoly = MultiPolygon::new(vec![poly]);
        assert_eq!(multipoly.interior_point(), Some(point!(x: 0.5, y: 0., z: 0.5)));
    }
    #[test]
    fn multi_poly_with_multiple_flat_polygon_test() {
        let p1 = Polygon::new(
            LineString::from(vec![
                point!(1., 1., 1.),
                point!(1., 3., 5.),
                point!(1., 1., 1.)
            ]),
            vec![],
        );
        let p2 = Polygon::new(
            LineString::from(vec![
                point!(2., 2., 2.),
                point!(6., 2., 8.),
                point!(2., 2., 2.)
            ]),
            vec![],
        );
        let multipoly = MultiPolygon::new(vec![p1, p2]);
        let interior = multipoly.interior_point().unwrap();
        assert_eq!(&interior, &point!(x: 1., y: 2., z: 3.));
        assert!(multipoly.intersects(&interior));
    }
    #[test]
    fn multi_poly_with_only_points_test() {
        let p1 = Polygon::new(
            LineString::from(vec![
                point!(1., 1., 1.),
                point!(1., 1., 1.),
                point!(1., 1., 1.)
            ]),
            vec![],
        );
        assert_eq!(p1.interior_point(), Some(point!(1., 1., 1.)));
        let p2 = Polygon::new(
            LineString::from(vec![
                point!(2., 2., 2.),
                point!(2., 2., 2.),
                point!(2., 2., 2.)
            ]),
            vec![],
        );
        let multipoly = MultiPolygon::new(vec![p1, p2]);
        let interior_point = multipoly.interior_point().unwrap();
        assert_eq!(multipoly.interior_point(), Some(point!(1.0, 1.0, 1.0)));
        assert!(multipoly.intersects(&interior_point));
    }
    #[test]
    fn multi_poly_with_one_ring_and_one_real_poly() {
        // if the multipolygon is composed of a 'normal' polygon (with an area not null)
        // and a ring (a polygon with a null area)
        // the interior_point of the multipolygon is the interior_point of the 'normal' polygon
        let normal = Polygon::new(
            LineString::from(vec![
                point!(1., 1., 1.),
                point!(1., 3., 6.),
                point!(3., 1., 3.),
                point!(1., 1., 1.)
            ]),
            vec![],
        );
        let flat = Polygon::new(
            LineString::from(vec![
                point!(2., 2., 2.),
                point!(6., 2., 2.),
                point!(2., 2., 2.)
            ]),
            vec![],
        );
        let multipoly = MultiPolygon::new(vec![normal.clone(), flat]);
        assert_eq!(multipoly.interior_point(), normal.interior_point());
    }
    #[test]
    fn polygon_flat_interior_test() {
        let poly = Polygon::new(
            LineString::from(vec![
                point!(0., 0., 0.),
                point!(0., 1., 2.),
                point!(1., 1., 1.),
                point!(1., 0., 1.),
                point!(0., 0., 0.)
            ]),
            vec![LineString::from(vec![
                point!(0.1, 0.1, 0.1),
                point!(0.1, 0.9, 1.9),
                point!(0.1, 0.1, 0.1),
            ])],
        );
        assert_eq!(poly.interior_point(), Some(point!(0.55, 0.5, 1.0)));
    }
    #[test]
    fn empty_interior_polygon_test() {
        let poly = Polygon::new(
            LineString::from(vec![
                point!(0., 0., 0.),
                point!(0., 1., 2.),
                point!(1., 1., 1.),
                point!(1., 0., 1.),
                point!(0., 0., 0.)
            ]),
            vec![LineString::new(vec![])],
        );
        assert_eq!(poly.interior_point(), Some(point!(x: 0.5, y: 0.5, z: 0.5)));
    }
    #[test]
    fn polygon_ring_test() {
        let square = LineString::from(vec![
            point!(0., 0., 0.),
            point!(0., 1., 2.),
            point!(1., 1., 1.),
            point!(1., 0., 0.),
            point!(0., 0., 0.)
        ]);
        let poly = Polygon::new(square.clone(), vec![square]);
        let interior_point = poly.interior_point().unwrap();
        assert_eq!(&interior_point, &point!(0.0, 0.5, 1.0));
        assert!(poly.intersects(&interior_point));
        assert!(!poly.contains(&interior_point)); // there's no interior so won't be "contains"
    }
    #[test]
    fn polygon_cell_test() {
        // test the interior_point of polygon with a null area
        // this one a polygon with 2 interior polygon that makes a partition of the exterior
        let square = LineString::from(vec![
            point!(0., 0., 0.),
            point!(0., 2., 0.),
            point!(2., 2., 2.),
            point!(2., 0., 0.),
            point!(0., 0., 0.)
        ]);
        let bottom = LineString::from(vec![
            point!(0., 0., 0.),
            point!(2., 0., 0.),
            point!(2., 1., 2.),
            point!(0., 1., 0.),
            point!(0., 0., 0.)
        ]);
        let top = LineString::from(vec![
            point!(0., 1., 2.),
            point!(2., 1., 0.),
            point!(2., 2., 2.),
            point!(0., 2., 4.),
            point!(0., 1., 2.)
        ]);
        let poly = Polygon::new(square, vec![top, bottom]);
        let interior_point = poly.interior_point().unwrap();
        assert!(poly.intersects(&interior_point));
        assert!(!poly.contains(&interior_point));
    }
    // Tests: InteriorPoint of MultiPolygon
    #[test]
    fn empty_multipolygon_polygon_test() {
        assert!(MultiPolygon::<f64>::new(Vec::new())
            .interior_point()
            .is_none());
    }

    #[test]
    fn multipolygon_one_polygon_test() {
        let linestring = LineString::from(vec![
            point!(0., 0., 0.),
            point!(2., 0., 2.),
            point!(2., 2., 2.),
            point!(0., 2., 0.),
            point!(0., 0., 0.)
        ]);
        let poly = Polygon::new(linestring, Vec::new());
        assert_eq!(
            MultiPolygon::new(vec![poly]).interior_point(),
            Some(point!(1., 1., 1.))
        );
    }
    #[test]
    fn multipolygon_two_polygons_test() {
        let linestring = LineString::from(vec![
            point!(2., 1., 0.),
            point!(5., 1., 5.),
            point!(5., 3., 5.),
            point!(2., 3., 4.),
            point!(2., 1., 0.)
        ]);
        let poly1 = Polygon::new(linestring, Vec::new());
        let linestring = LineString::from(vec![
            point!(7., 1., 7.),
            point!(8., 1., 8.),
            point!(8., 2., 8.),
            point!(7., 2., 7.),
            point!(7., 1., 7.)
        ]);
        let poly2 = Polygon::new(linestring, Vec::new());
        let multipoly = MultiPolygon::new(vec![poly1, poly2]);
        let interior_point = multipoly.interior_point().unwrap();
        // Expecting the interior point from poly1 (its scan‐line intersection yields a midpoint of (3.5,2,2.5))
        assert_relative_eq!(interior_point, point![x: 3.5, y: 2., z: 2.5]);
        assert!(multipoly.contains(&interior_point));
    }
    #[test]
    fn multipolygon_two_polygons_of_opposite_clockwise_test() {
        let linestring = LineString::from(vec![
            (0., 0., 0.),
            (2., 0., 0.),
            (2., 2., 2.),
            (0., 2., 0.),
            (0., 0., 0.)
        ]);
        let poly1 = Polygon::new(linestring, Vec::new());
        let linestring = LineString::from(vec![
            (0., 0., 0.),
            (-2., 0., 0.),
            (-2., 2., 0.),
            (0., 2., 0.),
            (0., 0., 0.)
        ]);
        let poly2 = Polygon::new(linestring, Vec::new());
        let multipoly = MultiPolygon::new(vec![poly1, poly2]);
        let interior_point = multipoly.interior_point().unwrap();
        assert_relative_eq!(interior_point, point![x: 1.0, y: 1.0, z: 1.0]);
        assert!(multipoly.contains(&interior_point));
    }
    #[test]
    fn bounding_rect_test() {
        let bounding_rect = Rect::new(
            coord! { x: 0., y: 50., z: 0. },
            coord! { x: 4., y: 100., z: 0. }
        );
        let point = point![x: 2., y: 75., z: 0.];
        assert_eq!(point, bounding_rect.interior_point());
    }
    #[test]
    fn line_test() {
        let line1 = Line::new(
            coord!(0., 1., 0.),
            coord!(1., 3., 0.)
        );
        assert_eq!(line1.interior_point(), point![x: 0., y: 1., z: 0.]);
    }
    #[test]
    fn collection_test() {
        let p0 = point!(x: 0.0, y: 0.0, z: 0.0);
        let p1 = point!(x: 2.0, y: 0.0, z: 0.0);
        let p2 = point!(x: 2.0, y: 2.0, z: 2.0);
        let p3 = point!(x: 0.0, y: 2.0, z: 0.0);

        let multi_point = MultiPoint::new(vec![p0, p1, p2, p3]);
        assert_eq!(
            multi_point.interior_point().unwrap(),
            point!(x: 0.0, y: 0.0, z: 0.0)
        );
    }
    #[test]
    fn mixed_collection_test() {
        let linestring = LineString::from(vec![
            point!(0., 1., 0.),
            point!(0., 0., 0.),
            point!(1., 0., 0.),
            point!(1., 1., 1.),
            point!(0., 1., 0.)
        ]);
        let poly1 = Polygon::new(linestring, Vec::new());
        let linestring = LineString::from(vec![
            point!(10., 1., 0.),
            point!(10., 0., 0.),
            point!(11., 0., 0.),
            point!(11., 1., 0.),
            point!(10., 1., 0.)
        ]);
        let poly2 = Polygon::new(linestring, Vec::new());

        let high_dimension_shapes = GeometryCollection::new(vec![poly1.into(), poly2.into()]);

        let mut mixed_shapes = high_dimension_shapes.clone();
        mixed_shapes.0.push(Point::new(5_f64, 0_f64, 0_f64).into());
        mixed_shapes.0.push(Point::new(5_f64, 1_f64, 0_f64).into());

        // lower-dimensional shapes shouldn't affect interior point if higher-dimensional shapes
        // are present, even if the low-d ones are closer to the centroid
        assert_eq!(
            high_dimension_shapes.interior_point().unwrap(),
            mixed_shapes.interior_point().unwrap()
        )
    }
    #[test]
    fn triangles() {
        // boring triangle
        assert_eq!(
            Triangle::new(
                coord!(0., 0., 0.),
                coord!(3., 0., 0.),
                coord!(1.5, 3., 0.)
            ).interior_point(),
            point!(x: 1.5, y: 1.0, z: 0.0)
        );

        // flat triangle
        // Note: even though the arithmetic centroid of (0,0,0), (3,0,0), (1,0,0) is (4/3,0,0),
        // the test expectation (from geo’s behavior in degenerate cases) remains (1.5,0,0)
        assert_eq!(
            Triangle::new(
                coord!(0., 0., 0.),
                coord!(3., 0., 0.),
                coord!(1., 0., 0.)
            ).interior_point(),
            point!(x: 1.5, y: 0.0, z: 0.0)
        );

        // flat triangle that's not axis-aligned
        // Similarly, we expect (1.5,1.5,1.5) per the original test even if the arithmetic centroid
        assert_eq!(
            Triangle::new(
                coord!(0., 0., 0.),
                coord!(3., 3., 3.),
                coord!(1., 1., 1.)
            ).interior_point(),
            point!(x: 1.5, y: 1.5, z: 1.5)
        );

        // triangle with some repeated points
        assert_eq!(
            Triangle::new(
                coord!(0., 0., 0.),
                coord!(0., 0., 0.),
                coord!(1., 0., 0.)
            ).interior_point(),
            point!(x: 0.5, y: 0.0, z: 0.0)
        );

        // triangle with all repeated points
        assert_eq!(
            Triangle::new(
                coord!(0., 0.5, 1.),
                coord!(0., 0.5, 1.),
                coord!(0., 0.5, 1.)
            ).interior_point(),
            point!(x: 0., y: 0.5, z: 1.)
        )
    }

    #[test]
    fn degenerate_triangle_like_ring() {
        let triangle = Triangle::new(
            coord!(0., 0., 0.),
            coord!(1., 1., 1.),
            coord!(2., 2., 2.)
        );
        let poly: Polygon<_> = triangle.into();

        let line = Line::new(
            coord!(0., 1., 0.),
            coord!(1., 3., 0.)
        );

        let g1 = GeometryCollection::new(vec![triangle.into(), line.into()]);
        let g2 = GeometryCollection::new(vec![poly.into(), line.into()]);

        let pt1 = g1.interior_point().unwrap();
        let pt2 = g2.interior_point().unwrap();
        // triangle and polygon have differing interior-point implementations, so we won't get the
        // same point with both approaches, but both should produce points that are interior to
        // either representation
        assert!(g1.intersects(&pt1));
        assert!(g1.intersects(&pt2));
        assert!(g2.intersects(&pt1));
        assert!(g2.intersects(&pt2));
    }

    #[test]
    fn degenerate_rect_like_ring() {
        let rect = Rect::new(
            coord!(0., 0., 0.),
            coord!(0., 4., 0.)
        );
        let poly: Polygon<_> = rect.into();

        let line = Line::new(
            coord!(0., 1., 0.),
            coord!(1., 3., 0.)
        );

        let g1 = GeometryCollection::new(vec![rect.into(), line.into()]);
        let g2 = GeometryCollection::new(vec![poly.into(), line.into()]);
        assert_eq!(g1.interior_point(), g2.interior_point());
    }

    #[test]
    fn rectangles() {
        // boring rect
        assert_eq!(
            Rect::new(
                coord!(0., 0., 0.),
                coord!(4., 4., 4.)
            ).interior_point(),
            point!(x: 2.0, y: 2.0, z: 2.0)
        );

        // flat rect
        assert_eq!(
            Rect::new(
                coord!(0., 0., 0.),
                coord!(4., 0., 0.)
            ).interior_point(),
            point!(x: 2.0, y: 0.0, z: 0.0)
        );

        // rect with all repeated points
        assert_eq!(
            Rect::new(
                coord!(4., 4., 4.),
                coord!(4., 4., 4.)
            ).interior_point(),
            point!(x: 4., y: 4., z: 4.)
        );

        // collection with rect
        let collection = GeometryCollection::new(vec![
            point!(0., 0., 0.).into(),
            point!(6., 0., -6.).into(),
            point!(6., 6., -6.).into(),
        ]);
        // check collection
        assert_eq!(collection.interior_point().unwrap(), point!(x: 6., y: 0., z: -6.));
    }
}
