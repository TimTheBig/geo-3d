use crate::{CoordNum, Point};

/// Interpolate a `Point` along a line between two existing points
pub trait InterpolatePoint<F: CoordNum> {
    /// Returns a new Point along a line between two existing points.
    ///
    /// See [specific implementations](#implementors) for details.
    fn point_at_distance_between(
        start: Point<F>,
        end: Point<F>,
        distance_from_start: F,
    ) -> Point<F>;

    /// Returns a new Point along a line between two existing points.
    ///
    /// See [specific implementations](#implementors) for details.
    fn point_at_ratio_between(start: Point<F>, end: Point<F>, ratio_from_start: F) -> Point<F>;

    /// Interpolates `Point`s along a line between `start` and `end`.
    ///
    /// See [specific implementations](#implementors) for details.
    ///
    /// As many points as necessary will be added such that the distance between points
    /// never exceeds `max_distance`. If the distance between start and end is less than
    /// `max_distance`, no additional points will be included in the output.
    ///
    /// `include_ends`: Should the start and end points be included in the output?
    fn points_along_line(
        start: Point<F>,
        end: Point<F>,
        max_distance: F,
        include_ends: bool,
    ) -> impl Iterator<Item = Point<F>>;
}

#[cfg(test)]
mod tests {
    use crate::{Euclidean, InterpolatePoint, Point};

    #[test]
    fn point_at_ratio_between_line_ends() {
        let start = Point::new(0.0, 0.0, 0.0);
        let end = Point::new(1.0, 1.0, 1.0);

        let ratio = 0.0;
        assert_eq!(Euclidean::point_at_ratio_between(start, end, ratio), start);

        let ratio = 1.0;
        assert_eq!(Euclidean::point_at_ratio_between(start, end, ratio), end);
    }

    mod degenerate {
        use super::*;

        #[test]
        fn point_at_ratio_between_collapsed_line() {
            let start = Point::new(1.0, 1.0, 1.0);

            let ratio = 0.0;
            assert_eq!(
                Euclidean::point_at_ratio_between(start, start, ratio),
                start
            );

            let ratio = 0.5;
            assert_eq!(
                Euclidean::point_at_ratio_between(start, start, ratio),
                start
            );

            let ratio = 1.0;
            assert_eq!(
                Euclidean::point_at_ratio_between(start, start, ratio),
                start
            );
        }

        #[test]
        fn point_at_distance_between_collapsed_line() {
            // This method just documents existing behavior. I don't think our current behavior
            // is especially useful, but we might consider handling it uniformly one day.
            let start: Point = Point::new(1.0, 1.0, 1.0);

            let distance = 0.0;

            let euclidean_result = Euclidean::point_at_distance_between(start, start, distance);
            assert!(euclidean_result.x().is_nan());
            assert!(euclidean_result.y().is_nan());
            assert!(euclidean_result.z().is_nan());

            let distance = 100000.0;

            let euclidean_result = Euclidean::point_at_distance_between(start, start, distance);
            assert!(euclidean_result.x().is_nan());
            assert!(euclidean_result.y().is_nan());
            assert!(euclidean_result.z().is_nan());
        }

        #[test]
        fn points_along_collapsed_line() {
            let start = Point::new(1.0, 1.0, 1.0);

            let max_distance = 1.0;

            let include_ends = true;

            let points: Vec<_> =
                Euclidean::points_along_line(start, start, max_distance, include_ends).collect();
            assert_eq!(points, vec![start, start]);

            let include_ends = false;

            let points: Vec<_> =
                Euclidean::points_along_line(start, start, max_distance, include_ends).collect();
            assert_eq!(points, vec![]);
        }
    }
}
