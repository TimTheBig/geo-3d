use super::{Distance, Euclidean};
use crate::{
    CoordNum, CoordsIter, InterpolatePoint, Line, LineString, MultiLineString, MultiPolygon, Point,
    Polygon, Rect, Triangle,
};
use num_traits::FromPrimitive;

/// Creates a copy of the geometry with additional points inserted as necessary to ensure there
/// is never more than `max_segment_length` between points.
///
/// ## Units
/// - `max_segment_length` units depend on the implementing [metric space]. It must be greater than 0.
///
/// # Examples
/// ```
/// # use approx::assert_relative_eq;
/// use geo::{wkt, Densify};
/// use geo::line_measures::Euclidean;
///
/// let line_string = wkt!(LINESTRING(0.0 0.0 0.0,0.0 6.0 9.0,1.0 7.0 9.0));
///
/// // For Euclidean calculations, the unit of distance is the same as the units
/// // of your coordinates.
/// let max_dist = 2.0;
/// let densified = line_string.densify(max_dist);
/// let expected_output = wkt!(LINESTRING(
///     0.0 0.0 0.0,
///     0.0 2.0 2.0,
///     0.0 4.0 6.0,
///     0.0 6.0 9.0,
///     1.0 7.0 9.0
/// ));
/// assert_relative_eq!(densified, expected_output);
///```
///
/// [metric space]: crate::line_measures::metric_spaces
pub trait Densify<F: CoordNum> {
    type Output;
    fn densify(&self, max_segment_length: F) -> Self::Output;
}

pub(crate) fn densify_between<F>(
    line_start: Point<F>,
    line_end: Point<F>,
    container: &mut Vec<Point<F>>,
    max_segment_length: F,
) where
    F: CoordNum + FromPrimitive,
{
    assert!(max_segment_length > F::zero());
    let num_segments = (line_start.distance(line_end) / max_segment_length)
        .ceil()
        .to_u64()
        .expect("unreasonable number of segments");

    // distance "unit" for this line segment
    let frac = F::one() / F::from(num_segments).unwrap();

    for segment_num in 1..num_segments {
        let ratio = frac * F::from(segment_num).unwrap();

        // PERF TODO: We recompute "total_distance" every step of this loop.
        // If we impl point_at_distance_between, we could compute it once and use it here.
        // At that point, I think this function could be a good candidate to be *the single* basis
        // for a unified generic of points_along_line for all metric spaces.
        let interpolated_point = Euclidean::point_at_ratio_between(line_start, line_end, ratio);
        container.push(interpolated_point);
    }
}

impl<F: CoordNum + FromPrimitive> Densify<F> for Line<F> {
    type Output = LineString<F>;

    fn densify(&self, max_segment_length: F) -> Self::Output {
        let mut points = vec![self.start_point()];
        densify_between::<F>(
            self.start_point(),
            self.end_point(),
            &mut points,
            max_segment_length,
        );
        points.push(self.end_point());
        LineString::from(points)
    }
}

impl<F: CoordNum + FromPrimitive> Densify<F> for LineString<F> {
    type Output = Self;

    fn densify(&self, max_segment_length: F) -> LineString<F> {
        if self.coords_count() == 0 {
            return LineString::new(vec![]);
        }

        let mut points = vec![];
        self.lines().for_each(|line| {
            points.push(line.start_point());
            densify_between::<F>(
                line.start_point(),
                line.end_point(),
                &mut points,
                max_segment_length,
            )
        });

        // we're done, push the last coordinate on to finish
        let final_coord = *self
            .0
            .last()
            .expect("we already asserted the line string is not empty");
        points.push(final_coord.into());

        LineString::from(points)
    }
}

impl<F: CoordNum + FromPrimitive> Densify<F> for MultiLineString<F> {
    type Output = Self;

    fn densify(&self, max_segment_length: F) -> Self::Output {
        MultiLineString::new(
            self.iter()
                .map(|line_string| line_string.densify(max_segment_length))
                .collect(),
        )
    }
}

impl<F: CoordNum + FromPrimitive> Densify<F> for Polygon<F> {
    type Output = Self;

    fn densify(&self, max_segment_length: F) -> Self::Output {
        Polygon::new(
            self.exterior().densify(max_segment_length),
            self.interiors()
                .iter()
                .map(|interior| interior.densify(max_segment_length))
                .collect(),
        )
    }
}

impl<F: CoordNum + FromPrimitive> Densify<F> for MultiPolygon<F> {
    type Output = Self;

    /// `max_segment_length` is the max distance between points
    fn densify(&self, max_segment_length: F) -> Self::Output {
        MultiPolygon::new(
            self.iter()
                .map(|polygon| polygon.densify(max_segment_length))
                .collect(),
        )
    }
}

impl<F: CoordNum + FromPrimitive> Densify<F> for Rect<F> {
    type Output = Polygon<F>;

    fn densify(&self, max_segment_length: F) -> Self::Output {
        self.to_polygon().densify(max_segment_length)
    }
}

impl<F: CoordNum + FromPrimitive> Densify<F> for Triangle<F> {
    type Output = Polygon<F>;

    fn densify(&self, max_segment_length: F) -> Self::Output {
        self.to_polygon().densify(max_segment_length)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{coord, polygon, wkt};

    #[test]
    fn densify_polygon() {
        let polygon = polygon![
            (x: -58.3816f64, y: -34.6037, z: 0.0), // Buenos Aires
            (x: -77.0428, y: -12.0464, z: 0.0),    // Lima
            (x: -47.9292, y: -15.7801, z: 0.0),    // Brasília
        ];

        let densified_polygon = polygon.densify(500_000.0); // 500 km max segment length
        assert!(densified_polygon.exterior().coords_count() > polygon.exterior().coords_count());
    }

    // ported from the old Deprecated trait, which only worked with Euclidean measures
    mod euclidean {
        use super::*;

        #[test]
        fn test_polygon_densify() {
            let polygon = wkt!(POLYGON(
                (-5.0 0.0 -5.0,0.0 5.0 0.0,5.0 0.0 5.0,-5.0 0.0 -5.0),
                (-3.0 0.0 -3.0,0.0 3.0 0.0,3.0 0.0 3.0,-3.0 0.0 -3.0)
            ));

            let expected = wkt!(POLYGON(
                (-5.0 0.0 -5.0,-3.75 1.25 -3.75,-2.5 2.5 -2.5,-1.25 3.75 -1.25,0.0 5.0 0.0,1.25 3.75 1.25,2.5 2.5 2.5,3.75 1.25 3.75,5.0 0.0 5.0,3.0 0.0 3.0,1.0 0.0 1.0,-1.0000000000000009 0.0 -1.0,-3.0 0.0 -3.0, -5.0 0.0 -5.0),
                (-3.0 0.0 -3.0,-2.0 1.0 -2.0,-1.0 2.0 -1.0,0.0 3.0 0.0,1.0 2.0 1.0,2.0 1.0 2.0,3.0 0.0 3.0,1.0 0.0 1.0,-1.0 0.0 -1.0,-3.0 0.0 -3.0)
            ));

            let max_dist = 2.0;
            let densified = polygon.densify(max_dist);
            assert_eq!(densified, expected);
        }

        #[test]
        fn test_empty_linestring_densify() {
            let linestring = LineString::<f64>::new(vec![]);
            let max_dist = 2.0;
            let densified = linestring.densify(max_dist);
            assert!(densified.0.is_empty());
        }

        #[test]
        fn test_linestring_densify() {
            let linestring = wkt!(LINESTRING(
               -1.0 0.0 -1.0,
                0.0 0.0 0.0,
                0.0 6.0 0.0,
                1.0 8.0 1.0
            ));
            let expected = wkt!(LINESTRING(
               -1.0 0.0 -1.0,
                0.0 0.0 0.0,
                0.0 2.0 0.0,
                0.0 4.0 0.0,
                0.0 6.0 0.0,
                0.5 7.0 0.5,
                1.0 8.0 1.0
            ));
            let max_dist = 2.0;
            let densified = linestring.densify(max_dist);
            assert_eq!(densified, expected);
        }

        #[test]
        fn test_line_densify() {
            let line: Line<f64> = Line::new(
                coord! { x: 0.0, y: 6.0, z: 0.0 },
                coord! { x: 1.0, y: 8.0, z: 0.0 },
            );
            let correct: LineString<f64> =
                vec![[0.0, 6.0, 0.0], [0.5, 7.0, 0.5], [1.0, 8.0, 1.0]].into();
            let max_dist = 2.0;
            let densified = line.densify(max_dist);
            assert_eq!(densified, correct);
        }
    }

    mod degenerate {
        use super::*;

        #[test]
        fn test_empty_linestring() {
            let input = wkt!(LINESTRING EMPTY);
            let dense = input.densify(1.0);
            assert_eq!(0, dense.coords_count());
            assert_eq!(input, dense);
        }

        #[test]
        fn test_one_point_linestring() {
            let input = wkt!(LINESTRING(1.0 1.0 1.0));
            let dense = input.densify(1.0);
            assert_eq!(1, dense.coords_count());
            assert_eq!(input, dense);
        }
    }
}
