use crate::{polygon, Coord, CoordNum, Line, Polygon};

#[cfg(any(feature = "approx", test))]
use approx::{AbsDiffEq, RelativeEq};

/// A bounded 3D area whose three co-planar vertices are defined by `Coord`s.\
/// The semantics and validity are that of
/// the equivalent [`Polygon`]; in addition, the three
/// vertices must not be collinear and they must be distinct.
#[derive(Copy, Clone, Hash, Eq, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Triangle<T: CoordNum = f64>(pub Coord<T>, pub Coord<T>, pub Coord<T>);

impl<T: CoordNum> Triangle<T> {
    /// Instantiate Self from the raw content value
    pub const fn new(v1: Coord<T>, v2: Coord<T>, v3: Coord<T>) -> Self {
        Self(v1, v2, v3)
    }

    pub const fn to_array(&self) -> [Coord<T>; 3] {
        [self.0, self.1, self.2]
    }

    pub const fn to_lines(&self) -> [Line<T>; 3] {
        [
            Line::new(self.0, self.1),
            Line::new(self.1, self.2),
            Line::new(self.2, self.0),
        ]
    }

    /// Create a `Polygon` from the `Triangle`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use geo_3d_types::{coord, Triangle, polygon};
    ///
    /// let triangle = Triangle::new(
    ///     coord! { x: 0., y: 0., z: 0. },
    ///     coord! { x: 10., y: 20., z: 30. },
    ///     coord! { x: 20., y: -10., z: 25. },
    /// );
    ///
    /// assert_eq!(
    ///     triangle.to_polygon(),
    ///     polygon![
    ///         (x: 0., y: 0., z: 0.),
    ///         (x: 10., y: 20., z: 30.),
    ///         (x: 20., y: -10., z: 25.),
    ///         (x: 0., y: 0., z: 0.),
    ///     ],
    /// );
    /// ```
    pub fn to_polygon(self) -> Polygon<T> {
        polygon![self.0, self.1, self.2, self.0]
    }
}

impl<IC: Into<Coord<T>> + Copy, T: CoordNum> From<[IC; 3]> for Triangle<T> {
    fn from(array: [IC; 3]) -> Self {
        Self(array[0].into(), array[1].into(), array[2].into())
    }
}

impl<T: CoordNum> From<&[Coord<T>; 3]> for Triangle<T> {
    fn from(coords: &[Coord<T>; 3]) -> Self {
        Triangle(coords[0], coords[1], coords[2])
    }
}

impl<T: CoordNum> From<(Coord<T>, Coord<T>, Coord<T>)> for Triangle<T> {
    fn from(coords: (Coord<T>, Coord<T>, Coord<T>)) -> Self {
        Triangle(coords.0, coords.1, coords.2)
    }
}

#[cfg(any(feature = "approx", test))]
impl<T> RelativeEq for Triangle<T>
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
    /// use geo_3d_types::{point, Triangle};
    ///
    /// let a = Triangle::new((0.0, 0.0, 0.0).into(), (10.0, 10.0, 10.0).into(), (0.0, 5.0, 0.0).into());
    /// let b = Triangle::new((0.0, 0.0, 0.0).into(), (10.01, 10.0, 10.0).into(), (0.0, 5.0, 0.0).into());
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
        if !self.0.relative_eq(&other.0, epsilon, max_relative) {
            return false;
        }
        if !self.1.relative_eq(&other.1, epsilon, max_relative) {
            return false;
        }
        if !self.2.relative_eq(&other.2, epsilon, max_relative) {
            return false;
        }

        true
    }
}

#[cfg(any(feature = "approx", test))]
impl<T> AbsDiffEq for Triangle<T>
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
    /// use geo_3d_types::{point, Triangle};
    ///
    /// let a = Triangle::new((0.0, 0.0, 0.0).into(), (10.0, 10.0, 10.0).into(), (0.0, 5.0, 0.0).into());
    /// let b = Triangle::new((0.0, 0.0, 0.0).into(), (10.01, 10.0, 10.0).into(), (0.0, 5.0, 0.0).into());
    ///
    /// approx::abs_diff_eq!(a, b, epsilon=0.1);
    /// approx::abs_diff_ne!(a, b, epsilon=0.001);
    /// ```
    #[inline]
    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        if !self.0.abs_diff_eq(&other.0, epsilon) {
            return false;
        }
        if !self.1.abs_diff_eq(&other.1, epsilon) {
            return false;
        }
        if !self.2.abs_diff_eq(&other.2, epsilon) {
            return false;
        }

        true
    }
}
