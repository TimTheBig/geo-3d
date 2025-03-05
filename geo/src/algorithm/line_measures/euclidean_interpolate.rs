use super::super::InterpolatePoint;
use crate::line_measures::densify::densify_between;
use crate::{CoordNum, Point};
use num_traits::FromPrimitive;

/// Operations on the [Euclidean plane] measure distance with the pythagorean formula -
/// what you'd measure with a ruler.
///
/// [Euclidean plane]: https://en.wikipedia.org/wiki/Euclidean_plane
pub struct Euclidean;

/// Interpolate Point(s) along a line in [Euclidean space].
///
/// [Euclidean plane]: https://en.wikipedia.org/wiki/Euclidean_plane
impl<F: CoordNum + FromPrimitive> InterpolatePoint<F> for Euclidean {
    /// Returns the point at the given distance along the line between `start` and `end`.\
    ///
    /// # Units
    /// - `distance`: Measured in whatever units your `start` and `end` points use.
    ///
    ///   `distance` and your `start` and `end` points should have non-angular
    ///   units, like meters or miles.
    ///
    /// [metric spaces]: crate::line_measures::metric_spaces
    fn point_at_distance_between(
        start: Point<F>,
        end: Point<F>,
        distance_from_start: F,
    ) -> Point<F> {
        let diff = end - start;
        let total_distance = diff.x().hypot(diff.y()).hypot(diff.z());
        let offset = diff * distance_from_start / total_distance;
        start + offset
    }

    /// Returns the point at the given ratio along the line between `start` and `end`.
    ///
    /// # Units
    /// - `distance`: Measured in whatever units your `start` and `end` points use.
    ///
    ///   `distance` and your `start` and `end` points should have non-angular
    ///   units, like meters or miles.
    ///
    /// [metric spaces]: crate::line_measures::metric_spaces
    fn point_at_ratio_between(start: Point<F>, end: Point<F>, ratio_from_start: F) -> Point<F> {
        let diff = end - start;
        start + diff * ratio_from_start
    }

    /// Interpolates `Point`s along a line between `start` and `end`.
    ///
    /// As many points as necessary will be added such that the distance between points
    /// never exceeds `max_distance`. If the distance between start and end is less than
    /// `max_distance`, no additional points will be included in the output.
    ///
    /// `include_ends`: Should the start and end points be included in the output?
    ///
    /// # Units
    /// - `max_distance`: Measured in whatever units your `start` and `end` points use.
    ///
    ///   `max_distance` and your `start` and `end` points should have non-angular
    ///   units, like meters or miles.
    ///
    /// [metric spaces]: crate::line_measures::metric_spaces
    fn points_along_line(
        start: Point<F>,
        end: Point<F>,
        max_distance: F,
        include_ends: bool,
    ) -> impl Iterator<Item = Point<F>> {
        let mut container = vec![];
        if include_ends {
            container.push(start);
        }
        densify_between::<F>(start, end, &mut container, max_distance);
        if include_ends {
            container.push(end);
        }
        container.into_iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    type MetricSpace = Euclidean;

    mod distance {
        use crate::Distance;

        use super::*;

        #[test]
        fn new_york_to_london() {
            // web mercator
            let new_york_city = Point::new(-8238310.24, 4942194.78, 0.0);
            // web mercator
            let london = Point::new(-14226.63, 6678077.70, 0.0);
            let distance: f64 = new_york_city.distance(london);

            assert_relative_eq!(
                8_405_286., // meters in web mercator
                distance.round()
            );
        }

        #[test]
        fn test_point_at_distance_between() {
            let new_york_city = Point::new(-8_238_310.24, 4_942_194.78, 6_678_077.70);
            // web mercator
            let london = Point::new(-14_226.63, 6_678_077.70, -8_238_310.24);
            let start = MetricSpace::point_at_distance_between(new_york_city, london, 0.0);
            assert_relative_eq!(new_york_city, start);

            let midway =
                MetricSpace::point_at_distance_between(new_york_city, london, 8_405_286.0 / 2.0);
            assert_relative_eq!(
                Point::new(-4_126_268., 5_810_136., -8_238_310.24),
                midway,
                epsilon = 1.0
            );

            let end = MetricSpace::point_at_distance_between(new_york_city, london, 8_405_286.0);
            assert_relative_eq!(london, end, epsilon = 1.0);
        }
    }
}
