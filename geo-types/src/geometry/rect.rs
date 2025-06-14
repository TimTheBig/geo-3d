use crate::{coord, polygon, Coord, CoordNum, Line, Polygon};

#[cfg(any(feature = "approx", test))]
use approx::{AbsDiffEq, RelativeEq};

/// An _axis-aligned_ bounded 3D rectangle whose area is
/// defined by minimum and maximum `Coord`s.
///
/// The constructors and setters ensure the maximum
/// `Coord` is greater than or equal to the minimum.
/// Thus, a `Rect`s width, depth, height, and area is guaranteed to
/// be greater than or equal to zero.
///
/// **Note.** While `Rect` implements `MapCoords` and
/// `RotatePoint` algorithmic traits, the usage is expected
/// to maintain the axis alignment. In particular, only
/// rotation by integer multiples of 90 degrees, will
/// preserve the original shape. In other cases, the min,
/// and max points are rotated or transformed, and a new
/// rectangle is created (with coordinate swaps to ensure
/// min < max).
///
/// # Examples
///
/// ```
/// use geo_3d_types::{coord, Rect};
///
/// let rect = Rect::new(
///     coord! { x: 0.0, y: 4.0, z: 1.0 },
///     coord! { x: 3.0, y: 10.0, z: 2.0 },
/// );
///
/// assert_eq!(3.0, rect.width());
/// assert_eq!(6.0, rect.depth());
/// assert_eq!(1.0, rect.height());
/// assert_eq!(
///     coord! { x: 1.5, y: 7., z: 1.5 },
///     rect.center()
/// );
/// ```
#[derive(Eq, PartialEq, Clone, Copy, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Rect<T: CoordNum = f64> {
    min: Coord<T>,
    max: Coord<T>,
}

impl<T: CoordNum> Rect<T> {
    /// Creates a new rectangle from two corner coordinates.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_3d_types::{coord, Rect};
    ///
    /// let rect = Rect::new(
    ///     coord! { x: 10., y: 20., z: 30. },
    ///     coord! { x: 30., y: 20., z: 10. }
    /// );
    /// assert_eq!(rect.min(), coord! { x: 10., y: 20., z: 10. });
    /// assert_eq!(rect.max(), coord! { x: 30., y: 20., z: 30. });
    /// ```
    pub fn new<C>(c1: C, c2: C) -> Self
    where
        C: Into<Coord<T>>,
    {
        let c1 = c1.into();
        let c2 = c2.into();

        let (min_x, max_x) = if c1.x < c2.x {
            (c1.x, c2.x)
        } else {
            (c2.x, c1.x)
        };
        let (min_y, max_y) = if c1.y < c2.y {
            (c1.y, c2.y)
        } else {
            (c2.y, c1.y)
        };
        let (min_z, max_z) = if c1.z < c2.z {
            (c1.z, c2.z)
        } else {
            (c2.z, c1.z)
        };

        Self {
            min: coord! { x: min_x, y: min_y, z: min_z },
            max: coord! { x: max_x, y: max_y, z: max_z },
        }
    }

    /// Returns the minimum `Coord` of the `Rect`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use geo_3d_types::{coord, Rect};
    ///
    /// let rect = Rect::new(
    ///     coord! { x: 5., y: 5., z: 5. },
    ///     coord! { x: 15., y: 15., z: 15. },
    /// );
    ///
    /// assert_eq!(rect.min(), coord! { x: 5., y: 5., z: 5. });
    /// ```
    pub const fn min(self) -> Coord<T> {
        self.min
    }

    /// Set the `Rect`’s minimum coordinate.
    ///
    /// # Panics
    ///
    /// Panics if `min`’s x/y/z is greater than the maximum coordinate’s x/y/z.
    pub fn set_min<C>(&mut self, min: C)
    where
        C: Into<Coord<T>>,
    {
        self.min = min.into();
        self.assert_valid_bounds();
    }

    /// Returns the maximum `Coord` of the `Rect`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use geo_3d_types::{coord, Rect};
    ///
    /// let rect = Rect::new(
    ///     coord! { x: 5., y: 5., z: 5. },
    ///     coord! { x: 15., y: 15., z: 15. },
    /// );
    ///
    /// assert_eq!(rect.max(), coord! { x: 15., y: 15., z: 15. });
    /// ```
    pub const fn max(self) -> Coord<T> {
        self.max
    }

    /// Set the `Rect`’s maximum coordinate.
    ///
    /// # Panics
    ///
    /// Panics if `max`’s x/y/z is less than the minimum coordinate’s x/y/z.
    pub fn set_max<C>(&mut self, max: C)
    where
        C: Into<Coord<T>>,
    {
        self.max = max.into();
        self.assert_valid_bounds();
    }

    /// Returns the width of the `Rect`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_3d_types::{coord, Rect};
    ///
    /// let rect = Rect::new(
    ///     coord! { x: 5., y: 5., z: 5. },
    ///     coord! { x: 15., y: 10., z: 10. },
    /// );
    ///
    /// assert_eq!(rect.width(), 10.);
    /// ```
    #[doc(alias = "dx")]
    pub fn width(self) -> T {
        self.max().x - self.min().x
    }

    /// Returns the depth of the `Rect`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_3d_types::{coord, Rect};
    ///
    /// let rect = Rect::new(
    ///     coord! { x: 5., y: 5., z: 5. },
    ///     coord! { x: 10., y: 15., z: 10. },
    /// );
    ///
    /// assert_eq!(rect.depth(), 10.);
    /// ```
    #[doc(alias = "dy")]
    pub fn depth(self) -> T {
        self.max().y - self.min().y
    }

    /// Returns the height of the `Rect`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_3d_types::{coord, Rect};
    ///
    /// let rect = Rect::new(
    ///     coord! { x: 5., y: 5., z: 5. },
    ///     coord! { x: 10., y: 10., z: 15. },
    /// );
    ///
    /// assert_eq!(rect.height(), 10.);
    /// ```
    #[doc(alias = "dz")]
    pub fn height(self) -> T {
        self.max().z - self.min().z
    }

    /// Create a `Polygon` from the `Rect`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use geo_3d_types::{coord, Rect, polygon};
    ///
    /// let rect = Rect::new(
    ///     coord! { x: 0., y: 0., z: 0. },
    ///     coord! { x: 1., y: 2., z: 3. },
    /// );
    ///
    /// assert_eq!(
    ///     rect.to_polygon(),
    ///     polygon![
    ///         (x: 0., y: 0., z: 0.),
    ///         (x: 0., y: 2., z: 2.),
    ///         (x: 1., y: 2., z: 3.),
    ///         (x: 1., y: 0., z: 2.),
    ///         (x: 0., y: 0., z: 2.),
    ///     ],
    /// );
    /// ```
    pub fn to_polygon(self) -> Polygon<T> {
        polygon![
            // top
            (x: self.min.x, y: self.min.y, z: self.max.z),
            (x: self.min.x, y: self.max.y, z: self.max.z),
            (x: self.max.x, y: self.max.y, z: self.max.z),
            (x: self.max.x, y: self.min.y, z: self.max.z),
            // bottom
            (x: self.max.x, y: self.min.y, z: self.min.z),
            (x: self.min.x, y: self.min.y, z: self.min.z),
            (x: self.min.x, y: self.max.y, z: self.min.z),
            (x: self.max.x, y: self.max.y, z: self.min.z),
            (x: self.min.x, y: self.min.y, z: self.max.z),
        ]
    }

    /// Returns the spatial representation of the `Rect`
    pub const fn to_coords(self) -> [Coord<T>; 8] {
        [
            // top
            coord!(x: self.min.x, y: self.min.y, z: self.max.z),
            coord!(x: self.min.x, y: self.max.y, z: self.max.z),
            coord!(x: self.max.x, y: self.max.y, z: self.max.z),
            coord!(x: self.max.x, y: self.min.y, z: self.max.z),
            // bottom
            coord!(x: self.max.x, y: self.min.y, z: self.min.z),
            coord!(x: self.min.x, y: self.min.y, z: self.min.z),
            coord!(x: self.min.x, y: self.max.y, z: self.min.z),
            coord!(x: self.max.x, y: self.max.y, z: self.min.z),
        ]
    }

    /// Return the lines of a `Rect`
    pub const fn to_lines(&self) -> [Line<T>; 12] {
        [
            // top
            Line::new(
                coord! { x: self.min.x, y: self.min.y, z: self.max.z },
                coord! { x: self.min.x, y: self.max.y, z: self.max.z },
            ),
            Line::new(
                coord! { x: self.min.x, y: self.max.y, z: self.max.z },
                coord! { x: self.max.x, y: self.max.y, z: self.max.z },
            ),
            Line::new(
                coord! { x: self.max.x, y: self.max.y, z: self.max.z },
                coord! { x: self.max.x, y: self.min.y, z: self.max.z },
            ),
            Line::new(
                coord! { x: self.max.x, y: self.min.y, z: self.max.z },
                coord! { x: self.min.x, y: self.min.y, z: self.max.z },
            ),
            // z connecting lines
            Line::new(
                coord! { x: self.min.x, y: self.min.y, z: self.max.z },
                coord! { x: self.min.x, y: self.min.y, z: self.min.z },
            ),
            Line::new(
                coord! { x: self.min.x, y: self.max.y, z: self.max.z },
                coord! { x: self.min.x, y: self.max.y, z: self.min.z },
            ),
            Line::new(
                coord! { x: self.max.x, y: self.max.y, z: self.max.z },
                coord! { x: self.max.x, y: self.max.y, z: self.min.z },
            ),
            Line::new(
                coord! { x: self.max.x, y: self.min.y, z: self.max.z },
                coord! { x: self.max.x, y: self.min.y, z: self.min.z },
            ),
            // bottom
            Line::new(
                coord! { x: self.min.x, y: self.min.y, z: self.min.z },
                coord! { x: self.min.x, y: self.max.y, z: self.min.z },
            ),
            Line::new(
                coord! { x: self.min.x, y: self.max.y, z: self.min.z },
                coord! { x: self.max.x, y: self.max.y, z: self.min.z },
            ),
            Line::new(
                coord! { x: self.max.x, y: self.max.y, z: self.min.z },
                coord! { x: self.max.x, y: self.min.y, z: self.min.z },
            ),
            Line::new(
                coord! { x: self.max.x, y: self.min.y, z: self.min.z },
                coord! { x: self.min.x, y: self.min.y, z: self.min.z },
            ),
        ]
    }

    /// Split a rectangle into two rectangles along the X-axis with equal widths.
    ///
    /// # Examples
    ///
    /// ```
    /// let rect = geo_3d_types::Rect::new(
    ///     geo_3d_types::coord! { x: 0., y: 0., z: 0. },
    ///     geo_3d_types::coord! { x: 4., y: 4., z: 4. },
    /// );
    ///
    /// let [rect1, rect2] = rect.split_x();
    ///
    /// assert_eq!(
    ///     geo_3d_types::Rect::new(
    ///         geo_3d_types::coord! { x: 0., y: 0., z: 0. },
    ///         geo_3d_types::coord! { x: 2., y: 4., z: 4. },
    ///     ),
    ///     rect1,
    /// );
    /// assert_eq!(
    ///     geo_3d_types::Rect::new(
    ///         geo_3d_types::coord! { x: 2., y: 0., z: 0. },
    ///         geo_3d_types::coord! { x: 4., y: 4., z: 4. },
    ///     ),
    ///     rect2,
    /// );
    /// ```
    pub fn split_x(self) -> [Rect<T>; 2] {
        let two = T::one() + T::one();
        let mid_x = self.min().x + self.width() / two;
        [
            Rect::new(self.min(), coord! { x: mid_x, y: self.max().y, z: self.max().z }),
            Rect::new(coord! { x: mid_x, y: self.min().y, z: self.min().z }, self.max()),
        ]
    }

    /// Split a rectangle into two rectangles along the Y-axis with equal depths.
    ///
    /// # Examples
    ///
    /// ```
    /// let rect = geo_3d_types::Rect::new(
    ///     geo_3d_types::coord! { x: 0., y: 0., z: 0. },
    ///     geo_3d_types::coord! { x: 4., y: 4., z: 4. },
    /// );
    ///
    /// let [rect1, rect2] = rect.split_y();
    ///
    /// assert_eq!(
    ///     geo_3d_types::Rect::new(
    ///         geo_3d_types::coord! { x: 0., y: 0., z: 0. },
    ///         geo_3d_types::coord! { x: 4., y: 2., z: 4. },
    ///     ),
    ///     rect1,
    /// );
    /// assert_eq!(
    ///     geo_3d_types::Rect::new(
    ///         geo_3d_types::coord! { x: 0., y: 2., z: 0. },
    ///         geo_3d_types::coord! { x: 4., y: 4., z: 4. },
    ///     ),
    ///     rect2,
    /// );
    /// ```
    pub fn split_y(self) -> [Rect<T>; 2] {
        let two = T::one() + T::one();
        let mid_y = self.min().y + self.depth() / two;
        [
            Rect::new(self.min(), coord! { x: self.max().x, y: mid_y, z: self.max().z }),
            Rect::new(coord! { x: self.min().x, y: mid_y, z: self.min().z }, self.max()),
        ]
    }

    /// Split a rectangle into two rectangles along the Z-axis with equal heights.
    ///
    /// # Examples
    ///
    /// ```
    /// let rect = geo_3d_types::Rect::new(
    ///     geo_3d_types::coord! { x: 0., y: 0., z: 0. },
    ///     geo_3d_types::coord! { x: 4., y: 4., z: 4. },
    /// );
    ///
    /// let [rect1, rect2] = rect.split_z();
    ///
    /// assert_eq!(
    ///     geo_3d_types::Rect::new(
    ///         geo_3d_types::coord! { x: 0., y: 0., z: 0. },
    ///         geo_3d_types::coord! { x: 4., y: 4., z: 2. },
    ///     ),
    ///     rect1,
    /// );
    /// assert_eq!(
    ///     geo_3d_types::Rect::new(
    ///         geo_3d_types::coord! { x: 0., y: 0., z: 2. },
    ///         geo_3d_types::coord! { x: 4., y: 4., z: 4. },
    ///     ),
    ///     rect2,
    /// );
    /// ```
    pub fn split_z(self) -> [Rect<T>; 2] {
        let two = T::one() + T::one();
        let mid_z = self.min().z + self.height() / two;
        [
            Rect::new(self.min(), coord! { x: self.max().x, y: self.max().y, z: mid_z }),
            Rect::new(coord! { x: self.min().x, y: self.min().y, z: mid_z }, self.max()),
        ]
    }

    fn assert_valid_bounds(&self) {
        if !self.has_valid_bounds() {
            panic!("{}", RECT_INVALID_BOUNDS_ERROR);
        }
    }

    fn has_valid_bounds(&self) -> bool {
        self.min.x <= self.max.x && self.min.y <= self.max.y && self.min.z <= self.max.z
    }
}

impl<T: CoordNum> Rect<T> {
    /// Returns the center `Coord` of the `Rect`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use geo_3d_types::{coord, Rect};
    ///
    /// let rect = Rect::new(
    ///     coord! { x: 5., y: 5., z: 5. },
    ///     coord! { x: 15., y: 15., z: 15. },
    /// );
    ///
    /// assert_eq!(rect.center(), coord! { x: 10., y: 10., z: 10. });
    /// ```
    pub fn center(self) -> Coord<T> {
        let two = T::one() + T::one();
        coord! {
            x: (self.max.x + self.min.x) / two,
            y: (self.max.y + self.min.y) / two,
            z: (self.max.z + self.min.z) / two,
        }
    }
}

static RECT_INVALID_BOUNDS_ERROR: &str =
    "Failed to create Rect: 'min' coordinate's x/y/z value must be smaller or equal to the 'max' x/y/z value";

#[cfg(any(feature = "approx", test))]
impl<T> RelativeEq for Rect<T>
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
    /// use geo_3d_types::Rect;
    ///
    /// let a = Rect::new((0.0, 0.0, 0.0), (10.0, 10.0, 10.0));
    /// let b = Rect::new((0.0, 0.0, 0.0), (10.01, 10.0, 10.0));
    ///
    /// approx::assert_relative_eq!(a, b, max_relative=0.1);
    /// approx::assert_relative_ne!(a, b, max_relative=0.0001);
    /// ```
    #[inline]
    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        if !self.min.relative_eq(&other.min, epsilon, max_relative) {
            return false;
        }

        if !self.max.relative_eq(&other.max, epsilon, max_relative) {
            return false;
        }

        true
    }
}

#[cfg(any(feature = "approx", test))]
impl<T> AbsDiffEq for Rect<T>
where
    T: AbsDiffEq<Epsilon = T> + CoordNum,
    T::Epsilon: Copy,
{
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
    /// use geo_3d_types::{point, Rect};
    ///
    /// let a = Rect::new((0.0, 0.0, 0.0), (10.0, 10.0, 10.0));
    /// let b = Rect::new((0.0, 0.0, 0.0), (10.01, 10.0, 10.0));
    ///
    /// approx::abs_diff_eq!(a, b, epsilon=0.1);
    /// approx::abs_diff_ne!(a, b, epsilon=0.001);
    /// ```
    #[inline]
    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        if !self.min.abs_diff_eq(&other.min, epsilon) {
            return false;
        }

        if !self.max.abs_diff_eq(&other.max, epsilon) {
            return false;
        }

        true
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use approx::assert_relative_eq;
    use crate::coord;

    #[test]
    fn rect() {
        let rect = Rect::new((10., 10., 10.), (20., 20., 20.));
        assert_eq!(rect.min, coord! { x: 10., y: 10., z: 10. });
        assert_eq!(rect.max, coord! { x: 20., y: 20., z: 20. });

        let rect = Rect::new((20., 20., 20.), (10., 10., 10.));
        assert_eq!(rect.min, coord! { x: 10., y: 10., z: 10. });
        assert_eq!(rect.max, coord! { x: 20., y: 20., z: 20. });

        let rect = Rect::new((10., 20., 20.), (20., 10., 10.));
        assert_eq!(rect.min, coord! { x: 10., y: 10., z: 10. });
        assert_eq!(rect.max, coord! { x: 20., y: 20., z: 20. });
    }

    #[test]
    fn rect_width() {
        let rect = Rect::new((10., 10., 10.), (20., 20., 20.));
        assert_eq!(rect.width(), 10.);
    }

    #[test]
    fn rect_depth() {
        let rect = Rect::new((10., 10., 10.1), (20., 20., 20.2));
        assert_relative_eq!(rect.depth(), 10.);
    }

    #[test]
    fn rect_height() {
        let rect = Rect::new((10., 10., 10.1), (20., 20., 20.2));
        assert_relative_eq!(rect.height(), 10.1);
    }

    #[test]
    fn rect_center() {
        assert_relative_eq!(
            Rect::new((0., 10., 0.), (10., 90., 10.)).center(),
            Coord::from((5., 50., 5.))
        );
        assert_relative_eq!(
            Rect::new((-42., -42., -5.), (42., 42., 5.)).center(),
            Coord::from((0., 0., 0.))
        );
        assert_relative_eq!(
            Rect::new((0., 0., 0.), (0., 0., 0.)).center(),
            Coord::from((0., 0., 0.))
        );
    }
}
