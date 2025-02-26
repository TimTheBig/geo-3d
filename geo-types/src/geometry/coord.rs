use crate::{coord, CoordNum, Point};

#[cfg(any(feature = "approx", test))]
use approx::{AbsDiffEq, RelativeEq, UlpsEq};

/// A lightweight struct used to store coordinates in 3-dimensional space.
///
/// Unlike `Point` (which in the future may contain additional information such
/// as an envelope, a precision model, and spatial reference system
/// information), a `Coord` only contains ordinate values and accessor
/// methods.
///
/// This type implements the [vector space] operations:
/// [`Add`], [`Sub`], [`Neg`], [`Zero`],
/// [`Mul<T>`][`Mul`], and [`Div<T>`][`Div`] traits.
///
/// # Semantics
///
/// This type does not represent any geospatial primitive,
/// but is used in their definitions. The only requirement
/// is that the coordinates it contains are valid numbers
/// (for eg. not `f64::NAN`).
///
/// [vector space]: //en.wikipedia.org/wiki/Vector_space
#[derive(Eq, PartialEq, Clone, Copy, Hash, Default)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Coord<T: CoordNum = f64> {
    pub x: T,
    pub y: T,
    #[cfg_attr(feature = "serde", serde(default = "T::zero"), serde(skip_serializing_if = "not_serializable"))]
    pub z: T,
}

/// This help serde to see if the z should be serialized, courntly this only skips zero
#[cfg(feature = "serde")]
fn not_serializable<T: CoordNum>(num: &T) -> bool {
    !num.is_zero()
}

impl<T: CoordNum> From<(T, T, T)> for Coord<T> {
    #[inline]
    fn from(coords: (T, T, T)) -> Self {
        coord! {
            x: coords.0,
            y: coords.1,
            z: coords.2,
        }
    }
}

impl<T: CoordNum> From<[T; 3]> for Coord<T> {
    #[inline]
    fn from(coords: [T; 3]) -> Self {
        coord! {
            x: coords[0],
            y: coords[1],
            z: coords[2],
        }
    }
}

impl<T: CoordNum> From<Point<T>> for Coord<T> {
    #[inline]
    fn from(point: Point<T>) -> Self {
        coord! {
            x: point.x(),
            y: point.y(),
            z: point.z(),
        }
    }
}

impl<T: CoordNum> From<Coord<T>> for (T, T, T) {
    #[inline]
    fn from(coord: Coord<T>) -> Self {
        (coord.x, coord.y, coord.z)
    }
}

impl<T: CoordNum> From<Coord<T>> for [T; 3] {
    #[inline]
    fn from(coord: Coord<T>) -> Self {
        [coord.x, coord.y, coord.z]
    }
}

impl<T: CoordNum> Coord<T> {
    // todo undeprecate
    #[deprecated = "coord is now 3d try `Coord::x_y_z`"]
    /// Returns a tuple that contains the x/horizontal & y/vertical component of the coordinate.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_types::coord;
    ///
    /// let c = coord! {
    ///     x: 40.02f64,
    ///     y: 116.34,
    ///     z: 397.01,
    /// };
    /// let (x, y) = c.x_y();
    ///
    /// assert_eq!(y, 116.34);
    /// assert_eq!(x, 40.02f64);
    /// ```
    #[inline]
    pub fn x_y(&self) -> (T, T) {
        (self.x, self.y)
    }

    /// Returns a tuple that contains the x/horizontal & y/vertical component of the coordinate.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_types::coord;
    ///
    /// let c = coord! {
    ///     x: 40.02f64,
    ///     y: 116.34,
    ///     z: 397.01,
    /// };
    /// let (x, y, z) = c.x_y_z();
    ///
    /// assert_eq!(y, 116.34);
    /// assert_eq!(x, 40.02f64);
    /// assert_eq!(z, 397.01);
    /// ```
    #[inline]
    pub fn x_y_z(&self) -> (T, T, T) {
        (self.x, self.y, self.z)
    }

    // todo document + example
    /// Get cross product of two `Coord`s
    /// 
    /// # Examples
    /// 
    /// ```
    /// let c1 = coord! {
    ///     x: 40.02f64,
    ///     y: 113.34,
    ///     z: 367.01,
    /// };
    /// 
    /// let c2 = coord! {
    ///     x: 55.0,
    ///     y: 116.0,
    ///     z: 497.21,
    /// };
    /// 
    /// let c3 = c1.cross(c2);
    /// 
    /// assert_eq!(c3.x, )
    /// assert_eq!(c3.y, )
    /// assert_eq!(c3.z, )
    /// ```
    pub fn cross(self, other: Coord<T>) -> Coord<T> {
        coord! {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }

    /// Computes the squared magnitude (length) of a 3D vector.
    /// This avoids the overhead of calculating the square root compared to `magnitude()`.
    pub fn magnitude_squared(self) -> T {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    /// Returns the dot product of the two coords:
    /// `dot = x1 * x2 + y1 * y2 + z1 * z2`
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_types::{coord, Coord};
    ///
    /// let coord = coord! { x: 1.5, y: 0.5, z: 8.5 };
    /// let dot = coord.dot(coord! { x: 2.0, y: 4.5, z: 6.3 });
    ///
    /// assert_eq!(dot, 58.8);
    /// ```
    pub fn dot(self, other: Self) -> T {
        self.x * other.x + self.y * other.y + self.z * other.z
    }
}

use core::ops::{Add, Div, Mul, Neg, Sub};

/// Negate a coordinate.
///
/// # Examples
///
/// ```
/// use geo_types::coord;
///
/// let p = coord! { x: 1.25, y: 2.5, z: 5.0 };
/// let q = -p;
///
/// assert_eq!(q.x, -p.x);
/// assert_eq!(q.y, -p.y);
/// assert_eq!(q.z, -p.z);
/// ```
impl<T> Neg for Coord<T>
where
    T: CoordNum + Neg<Output = T>,
{
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        coord! {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

/// Add two coordinates.
///
/// # Examples
///
/// ```
/// use geo_types::coord;
///
/// let p = coord! { x: 1.25, y: 2.5, z: 5.0 };
/// let q = coord! { x: 1.5, y: 2.5, z: 5.5 };
/// let sum = p + q;
///
/// assert_eq!(sum.x, 2.75);
/// assert_eq!(sum.y, 5.0);
/// assert_eq!(sum.z, 10.5);
/// ```
impl<T: CoordNum> Add for Coord<T> {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self {
        coord! {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

/// Subtract a coordinate from another.
///
/// # Examples
///
/// ```
/// use geo_types::coord;
///
/// let p = coord! { x: 1.5, y: 2.5, z: 5.0 };
/// let q = coord! { x: 1.25, y: 2.5, z: 5.5 };
/// let diff = p - q;
///
/// assert_eq!(diff.x, 0.25);
/// assert_eq!(diff.y, 0.);
/// assert_eq!(diff.z, 0.5);
/// ```
impl<T: CoordNum> Sub for Coord<T> {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        coord! {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

/// Multiply coordinate wise by a scalar.
///
/// # Examples
///
/// ```
/// use geo_types::coord;
///
/// let p = coord! { x: 1.25, y: 2.5, z: 5.0 };
/// let q = p * 4.;
///
/// assert_eq!(q.x, 5.0);
/// assert_eq!(q.y, 10.0);
/// assert_eq!(q.z, 20.0);
/// ```
impl<T: CoordNum> Mul<T> for Coord<T> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: T) -> Self {
        coord! {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

/// Divide coordinate wise by a scalar.
///
/// # Examples
///
/// ```
/// use geo_types::coord;
///
/// let p = coord! { x: 5., y: 10., z: 15. };
/// let q = p / 4.;
///
/// assert_eq!(q.x, 1.25);
/// assert_eq!(q.y, 2.5);
/// assert_eq!(q.z, 3.75);
/// ```
impl<T: CoordNum> Div<T> for Coord<T> {
    type Output = Self;

    #[inline]
    fn div(self, rhs: T) -> Self {
        coord! {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
        }
    }
}

use num_traits::Zero;
/// Create a coordinate at the origin.
///
/// # Examples
///
/// ```
/// use geo_types::Coord;
/// use num_traits::Zero;
///
/// let p: Coord = Zero::zero(); // or Coord::zero()
///
/// assert_eq!(p.x, 0.0);
/// assert_eq!(p.y, 0.0);
/// assert_eq!(p.z, 0.0);
/// ```
impl<T: CoordNum> Coord<T> {
    #[inline]
    pub fn zero() -> Self {
        coord! {
            x: T::zero(),
            y: T::zero(),
            z: T::zero(),
        }
    }
}

impl<T: CoordNum> Zero for Coord<T> {
    #[inline]
    fn zero() -> Self {
        Self::zero()
    }
    #[inline]
    fn is_zero(&self) -> bool {
        self.x.is_zero() && self.y.is_zero() && self.z.is_zero()
    }
}

#[cfg(any(feature = "approx", test))]
impl<T: CoordNum + AbsDiffEq> AbsDiffEq for Coord<T>
where
    T::Epsilon: Copy,
{
    type Epsilon = T::Epsilon;

    #[inline]
    fn default_epsilon() -> T::Epsilon {
        T::default_epsilon()
    }

    #[inline]
    fn abs_diff_eq(&self, other: &Self, epsilon: T::Epsilon) -> bool {
        T::abs_diff_eq(&self.x, &other.x, epsilon)
            && T::abs_diff_eq(&self.y, &other.y, epsilon)
            && T::abs_diff_eq(&self.z, &other.z, epsilon)
    }
}

#[cfg(any(feature = "approx", test))]
impl<T: CoordNum + RelativeEq> RelativeEq for Coord<T>
where
    T::Epsilon: Copy,
{
    #[inline]
    fn default_max_relative() -> T::Epsilon {
        T::default_max_relative()
    }

    #[inline]
    fn relative_eq(&self, other: &Self, epsilon: T::Epsilon, max_relative: T::Epsilon) -> bool {
        T::relative_eq(&self.x, &other.x, epsilon, max_relative)
            && T::relative_eq(&self.y, &other.y, epsilon, max_relative)
            && T::relative_eq(&self.z, &other.z, epsilon, max_relative)
    }
}

#[cfg(any(feature = "approx", test))]
impl<T: CoordNum + UlpsEq> UlpsEq for Coord<T>
where
    T::Epsilon: Copy,
{
    #[inline]
    fn default_max_ulps() -> u32 {
        T::default_max_ulps()
    }

    #[inline]
    fn ulps_eq(&self, other: &Self, epsilon: T::Epsilon, max_ulps: u32) -> bool {
        T::ulps_eq(&self.x, &other.x, epsilon, max_ulps)
            && T::ulps_eq(&self.y, &other.y, epsilon, max_ulps)
            && T::ulps_eq(&self.z, &other.z, epsilon, max_ulps)
    }
}

#[cfg(feature = "rstar")]
impl<T> ::rstar::Point for Coord<T>
where
    T: ::num_traits::Float + ::rstar::RTreeNum,
{
    type Scalar = T;

    const DIMENSIONS: usize = 3;

    #[inline]
    fn generate(mut generator: impl FnMut(usize) -> Self::Scalar) -> Self {
        coord! {
            x: generator(0),
            y: generator(1),
            z: generator(2),
        }
    }

    #[inline]
    fn nth(&self, index: usize) -> Self::Scalar {
        match index {
            0 => self.x,
            1 => self.y,
            2 => self.z,
            _ => unreachable!(),
        }
    }

    #[inline]
    fn nth_mut(&mut self, index: usize) -> &mut Self::Scalar {
        match index {
            0 => &mut self.x,
            1 => &mut self.y,
            2 => &mut self.z,
            _ => unreachable!(),
        }
    }
}
