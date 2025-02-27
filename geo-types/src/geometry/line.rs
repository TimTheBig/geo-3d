use crate::{Coord, CoordNum, Point};
#[cfg(any(feature = "approx", test))]
use approx::{AbsDiffEq, RelativeEq};

/// A line segment made up of exactly two
/// [`Coord`]s.
///
/// # Semantics
///
/// The _interior_ and _boundary_ are defined as with a
/// `LineString` with the two end points.
#[derive(Eq, PartialEq, Clone, Copy, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Line<T: CoordNum = f64> {
    pub start: Coord<T>,
    pub end: Coord<T>,
}

impl<T: CoordNum> Line<T> {
    /// Creates a new line segment.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_types::{coord, Line};
    ///
    /// let line = Line::new(coord! { x: 0., y: 0., z: 0. }, coord! { x: 1., y: 2., z: 3. });
    ///
    /// assert_eq!(line.start, coord! { x: 0., y: 0., z: 0. });
    /// assert_eq!(line.end, coord! { x: 1., y: 2., z: 3. });
    /// ```
    pub fn new<C>(start: C, end: C) -> Self
    where
        C: Into<Coord<T>>,
    {
        Self {
            start: start.into(),
            end: end.into(),
        }
    }

    /// Calculate the difference in coordinates (Δx, Δy, Δz).
    /// Equivalent to `self.end - self.start`
    pub fn delta(&self) -> Coord<T> {
        self.end - self.start
    }

    /// Calculate the difference in ‘x’ components (Δx).
    ///
    /// Equivalent to:
    ///
    /// ```rust
    /// # use geo_types::{Line, point};
    /// # let line = Line::new(
    /// #     point! { x: 4., y: -12., z: 0. },
    /// #     point! { x: 0., y: 9., z: 0. },
    /// # );
    /// # assert_eq!(
    /// #     line.dx(),
    /// line.end.x - line.start.x
    /// # );
    /// ```
    pub fn dx(&self) -> T {
        self.delta().x
    }

    /// Calculate the difference in ‘y’ components (Δy).
    ///
    /// Equivalent to:
    ///
    /// ```rust
    /// # use geo_types::{Line, point};
    /// # let line = Line::new(
    /// #     point! { x: 4., y: -12., z: 0. },
    /// #     point! { x: 0., y: 9., z: 0. },
    /// # );
    /// # assert_eq!(
    /// #     line.dy(),
    /// line.end.y - line.start.y
    /// # );
    /// ```
    pub fn dy(&self) -> T {
        self.delta().y
    }

    /// Calculate the difference in ‘z’ components (Δz).
    ///
    /// Equivalent to:
    ///
    /// ```rust
    /// # use geo_types::{Line, point};
    /// # let line = Line::new(
    /// #     point! { x: 4., y: -12., z: 0.3 },
    /// #     point! { x: 0., y: 9., z: 1. },
    /// # );
    /// # assert_eq!(
    /// #     line.dz(),
    /// line.end.z - line.start.z
    /// # );
    /// ```
    pub fn dz(&self) -> T {
        self.delta().z
    }

    /// Calculate the generalized 3D slope.
    /// 
    /// The slope in 3D is defined as the magnitude of the direction vector (Δx, Δy, Δz),
    /// normalized by the horizontal displacement (√(Δx² + Δy²)).
    /// This generalization ensures a scalar value representing the steepness of the line.
    /// 
    /// # Examples
    /// 
    /// ```rust
    /// use geo_types::{Line, point};
    /// 
    /// // A line with a steep vertical displacement
    /// let line1 = Line::new(
    ///     point! { x: 0., y: 0., z: 0. },
    ///     point! { x: 3., y: 4., z: 5. },
    /// );
    /// assert_eq!(line1.slope(), 5. / 5.); // dz / sqrt(dx² + dy²)
    /// 
    /// // A purely horizontal line
    /// let line2 = Line::new(
    ///     point! { x: 0., y: 0., z: 0. },
    ///     point! { x: 6., y: 8., z: 0. },
    /// );
    /// assert_eq!(line2.slope(), 0.); // dz = 0
    /// 
    /// // A purely vertical line
    /// let line3 = Line::new(
    ///     point! { x: 1., y: 1., z: 0. },
    ///     point! { x: 1., y: 1., z: 10. },
    /// );
    /// assert!(line3.slope().is_infinite()); // dz = 10, dx and dy = 0
    /// ```
    pub fn slope(&self) -> T {
        // The 3D slope is the ratio of the vertical magnitude to the horizontal magnitude.
        // If Δx² + Δy² == 0, this implies a vertical line, so the slope is infinite (or undefined).
        let horizontal_magnitude = T::hypot(self.dx(), self.dy());
        if horizontal_magnitude == T::zero() {
            T::infinity()
        } else {
            self.dz() / horizontal_magnitude
        }
    }

    /// Calculate the [determinant](https://en.wikipedia.org/wiki/Determinant) of the line in 3D.
    ///
    /// In 3D, the determinant is generalized as the magnitude of the cross product of the
    /// two points' position vectors (relative to the origin). This represents the signed
    /// "volume" or "area" spanned by the two points in 3D space.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use geo_types::{Line, point};
    /// # let line = Line::new(
    /// #     point! { x: 4., y: -12., z: 4. },
    /// #     point! { x: 0., y: 9., z: 0. },
    /// # );
    /// # assert_eq!(
    /// #     line.determinant(),
    /// #     ((4f64 * 9. - (-12.) * 0.).powi(2) + (4f64 * 0. - 4. * 0.).powi(2) + ((-12.) * 0. - 9. * 4f64).powi(2)).sqrt()
    /// # );
    /// ```
    pub fn determinant(&self) -> T {
        // Calculate the cross product of the two vectors
        let cross = self.start.cross(self.end);

        // Return the magnitude of the resulting 3D vector (the "volume" determinant)
        T::hypot(T::hypot(cross.x, cross.y), cross.z)
    }

    pub fn start_point(&self) -> Point<T> {
        Point::from(self.start)
    }

    pub fn end_point(&self) -> Point<T> {
        Point::from(self.end)
    }

    pub fn points(&self) -> (Point<T>, Point<T>) {
        (self.start_point(), self.end_point())
    }
}

impl<T: CoordNum> From<[(T, T, T); 2]> for Line<T> {
    fn from(coord: [(T, T, T); 2]) -> Self {
        Line::new(coord[0], coord[1])
    }
}
#[cfg(any(feature = "approx", test))]
impl<T> RelativeEq for Line<T>
where
    T: AbsDiffEq<Epsilon = T> + CoordNum + RelativeEq,
{
    #[inline]
    fn default_max_relative() -> Self::Epsilon {
        T::default_max_relative()
    }

    /// Equality assertion within a relative limit.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_types::{coord, Line};
    ///
    /// let a = Line::new(coord! { x: 0., y: 0., z: 0. }, coord! { x: 1., y: 1., z: 1. });
    /// let b = Line::new(coord! { x: 0., y: 0., z: 0. }, coord! { x: 1.001, y: 1., z: 1. });
    ///
    /// approx::assert_relative_eq!(a, b, max_relative=0.1);
    /// ```
    #[inline]
    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        self.start.relative_eq(&other.start, epsilon, max_relative)
            && self.end.relative_eq(&other.end, epsilon, max_relative)
    }
}

#[cfg(any(feature = "approx", test))]
impl<T: AbsDiffEq<Epsilon = T> + CoordNum> AbsDiffEq for Line<T> {
    type Epsilon = T;

    #[inline]
    fn default_epsilon() -> Self::Epsilon {
        T::default_epsilon()
    }

    /// Equality assertion with an absolute limit.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_types::{coord, Line};
    ///
    /// let a = Line::new(coord! { x: 0., y: 0., z: 0. }, coord! { x: 1., y: 1., z: 1. });
    /// let b = Line::new(coord! { x: 0., y: 0., z: 0. }, coord! { x: 1.001, y: 1., z: 1. });
    ///
    /// approx::assert_abs_diff_eq!(a, b, epsilon=0.1);
    /// ```
    #[inline]
    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.start.abs_diff_eq(&other.start, epsilon) && self.end.abs_diff_eq(&other.end, epsilon)
    }
}

#[cfg(feature = "rstar")]
impl <T> rstar::RTreeObject for Line<T>
    where T: num_traits::Float + rstar::RTreeNum
{
    type Envelope = rstar::AABB<Point<T>>;

    fn envelope(&self) -> Self::Envelope {
        let bounding_rect = crate::private_utils::line_bounding_rect(*self);
        rstar::AABB::from_corners(bounding_rect.min().into(), bounding_rect.max().into())
    }
}

#[cfg(feature = "rstar")]
impl <T> rstar::PointDistance for Line<T>
    where T: num_traits::Float + rstar::RTreeNum
{
    fn distance_2(&self, point: &Point<T>) -> T {
        let d = crate::private_utils::point_line_euclidean_distance(*point, *self);
        d.powi(2)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{coord, point};

    #[test]
    fn test_abs_diff_eq() {
        let delta = 1e-6;
        let line = Line::new(coord! { x: 0., y: 0., z: 0. }, coord! { x: 1., y: 1., z: 1. });
        let line_start_x = Line::new(
            point! {
                x: 0. + delta,
                y: 0.,
                z: 0.,
            },
            point! { x: 1., y: 1., z: 1. },
        );
        assert!(line.abs_diff_eq(&line_start_x, 1e-2));
        assert!(line.abs_diff_ne(&line_start_x, 1e-12));

        let line_start_y = Line::new(
            coord! {
                x: 0.,
                y: 0. + delta,
                z: 0.,
            },
            coord! { x: 1., y: 1., z: 1. },
        );
        assert!(line.abs_diff_eq(&line_start_y, 1e-2));
        assert!(line.abs_diff_ne(&line_start_y, 1e-12));

        let line_end_x = Line::new(
            coord! { x: 0., y: 0., z: 0. },
            coord! {
                x: 1. + delta,
                y: 1.,
                z: 1.,
            },
        );

        assert!(line.abs_diff_eq(&line_end_x, 1e-2));
        assert!(line.abs_diff_ne(&line_end_x, 1e-12));

        let line_end_y = Line::new(
            coord! { x: 0., y: 0., z: 0. },
            coord! {
                x: 1.,
                y: 1. + delta,
                z: 1. - delta,
            },
        );

        assert!(line.abs_diff_eq(&line_end_y, 1e-2));
        assert!(line.abs_diff_ne(&line_end_y, 1e-12));
    }

    #[test]
    fn test_relative_eq() {
        let delta = 1e-6;

        let line = Line::new(coord! { x: 0., y: 0., z: 0. }, coord! { x: 1., y: 1., z: 1. });
        let line_start_x = Line::new(
            point! {
                x: 0. + delta,
                y: 0.,
                z: 0.
            },
            point! { x: 1., y: 1., z: 1. },
        );
        let line_start_y = Line::new(
            coord! {
                x: 0.,
                y: 0. + delta,
                z: 0. - delta,
            },
            coord! { x: 1., y: 1., z: 1. },
        );

        assert!(line.relative_eq(&line_start_x, 1e-2, 1e-2));
        assert!(line.relative_ne(&line_start_x, 1e-12, 1e-12));

        assert!(line.relative_eq(&line_start_y, 1e-2, 1e-2));
        assert!(line.relative_ne(&line_start_y, 1e-12, 1e-12));
    }
}
