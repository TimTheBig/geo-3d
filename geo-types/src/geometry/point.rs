use crate::{point, Coord, CoordNum};

#[cfg(any(feature = "approx", test))]
use approx::{AbsDiffEq, RelativeEq};

use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// A single point in 2D space.
///
/// Points can be created using the [`Point::new`] constructor,
/// the [`point!`] macro, or from a `Coord`, two-element
/// tuples, or arrays – see the `From` impl section for a
/// complete list.
///
/// # Semantics
///
/// The _interior_ of the point is itself (a singleton set),
/// and its _boundary_ is empty. A point is _valid_ if and
/// only if the `Coord` is valid.
///
/// # Examples
///
/// ```
/// use geo_types::{coord, Point};
/// let p1: Point = (0., 1.).into();
/// let c = coord! { x: 10., y: 20. };
/// let p2: Point = c.into();
/// ```
#[derive(Eq, PartialEq, Clone, Copy, Debug, Hash, Default)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Point<T: CoordNum = f64>(pub Coord<T>);

impl<T: CoordNum> From<Coord<T>> for Point<T> {
    fn from(x: Coord<T>) -> Self {
        Point(x)
    }
}

impl<T: CoordNum> From<(T, T, T)> for Point<T> {
    fn from(coords: (T, T, T)) -> Self {
        Point::new(coords.0, coords.1, coords.2)
    }
}

impl<T: CoordNum> From<[T; 3]> for Point<T> {
    fn from(coords: [T; 3]) -> Self {
        Point::new(coords[0], coords[1], coords[2])
    }
}

impl<T: CoordNum> From<Point<T>> for (T, T, T) {
    fn from(point: Point<T>) -> Self {
        point.0.into()
    }
}

impl<T: CoordNum> From<Point<T>> for [T; 3] {
    fn from(point: Point<T>) -> Self {
        point.0.into()
    }
}

impl<T: CoordNum> Point<T> {
    /// Creates a new point.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_types::Point;
    ///
    /// let p = Point::new(1.234, 2.345, 8.735);
    ///
    /// assert_eq!(p.x(), 1.234);
    /// assert_eq!(p.y(), 2.345);
    /// assert_eq!(p.z(), 8.735);
    /// ```
    pub fn new(x: T, y: T, z: T) -> Self {
        point! { x: x, y: y, z: z }
    }

    /// Returns the x/horizontal component of the point.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_types::Point;
    ///
    /// let p = Point::new(1.234, 2.345, 8.735);
    ///
    /// assert_eq!(p.x(), 1.234);
    /// ```
    pub fn x(self) -> T {
        self.0.x
    }

    /// Sets the x/horizontal component of the point.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_types::Point;
    ///
    /// let mut p = Point::new(1.234, 2.345, 8.735);
    /// p.set_x(9.876);
    ///
    /// assert_eq!(p.x(), 9.876);
    /// ```
    pub fn set_x(&mut self, x: T) -> &mut Self {
        self.0.x = x;
        self
    }

    /// Returns a mutable reference to the x/horizontal component of the point
    ///
    /// # Examples
    ///
    /// ```
    /// use approx::assert_relative_eq;
    /// use geo_types::Point;
    /// let mut p = Point::new(1.234, 2.345, 8.735);
    /// let mut p_x = p.x_mut();
    /// *p_x += 1.0;
    /// assert_relative_eq!(p.x(), 2.234);
    /// ```
    pub fn x_mut(&mut self) -> &mut T {
        &mut self.0.x
    }
    /// Returns the y/vertical component of the point.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_types::Point;
    ///
    /// let p = Point::new(1.234, 2.345, 8.735);
    ///
    /// assert_eq!(p.y(), 2.345);
    /// ```
    pub fn y(self) -> T {
        self.0.y
    }

    /// Sets the y/vertical component of the point.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_types::Point;
    ///
    /// let mut p = Point::new(1.234, 2.345, 8.735);
    /// p.set_y(9.876);
    ///
    /// assert_eq!(p.y(), 9.876);
    /// ```
    pub fn set_y(&mut self, y: T) -> &mut Self {
        self.0.y = y;
        self
    }

    /// Returns a mutable reference to the y/vertical component of the point
    ///
    /// # Examples
    ///
    /// ```
    /// use approx::assert_relative_eq;
    /// use geo_types::Point;
    /// let mut p = Point::new(1.234, 2.345, 8.735);
    /// let mut p_y = p.y_mut();
    /// *p_y += 1.0;
    /// assert_relative_eq!(p.y(), 3.345);
    /// ```
    pub fn y_mut(&mut self) -> &mut T {
        &mut self.0.y
    }

    /// Returns the z/hight component of the point.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_types::Point;
    ///
    /// let p = Point::new(1.234, 2.345, 4.382);
    ///
    /// assert_eq!(p.z(), 4.382);
    /// ```
    pub fn z(self) -> T {
        self.0.z
    }

    /// Sets the z/hight component of the point.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_types::Point;
    ///
    /// let mut p = Point::new(1.234, 2.345, 4.380);
    /// p.set_z(9.876);
    ///
    /// assert_eq!(p.z(), 9.876);
    /// ```
    pub fn set_z(&mut self, z: T) -> &mut Self {
        self.0.z = z;
        self
    }

    /// Returns a mutable reference to the z/hight component of the point
    ///
    /// # Examples
    ///
    /// ```
    /// use approx::assert_relative_eq;
    /// use geo_types::Point;
    /// let mut p = Point::new(1.234, 2.345, 4.380);
    /// let mut p_y = p.z_mut();
    /// *p_y += 1.0;
    /// assert_relative_eq!(p.y(), 5.380);
    /// ```
    pub fn z_mut(&mut self) -> &mut T {
        &mut self.0.z
    }

    /// Returns a tuple that contains the x/horizontal & y/vertical component of the point.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_types::Point;
    ///
    /// let mut p = Point::new(1.234, 2.345, 3.872);
    /// let (x, y) = p.x_y();
    ///
    /// assert_eq!(y, 2.345);
    /// assert_eq!(x, 1.234);
    /// ```
    pub fn x_y(self) -> (T, T) {
        (self.0.x, self.0.y)
    }

    /// Returns a tuple that contains the x/horizontal, y/vertical & z/hight component of the point.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_types::Point;
    ///
    /// let mut p = Point::new(1.234, 2.345, 7.509);
    /// let (x, y, z) = p.x_y_z();
    ///
    /// assert_eq!(y, 2.345);
    /// assert_eq!(x, 1.234);
    /// assert_eq!(z, 7.509);
    /// ```
    pub fn x_y_z(self) -> (T, T, T) {
        (self.0.x, self.0.y, self.0.z)
    }
}

impl<T: CoordNum> Point<T> {
    /// Returns the dot product of the two points:
    /// `dot = x1 * x2 + y1 * y2 + z1 * z2`
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_types::{point, Point};
    ///
    /// let point = point! { x: 1.5, y: 0.5, z: 8.5 };
    /// let dot = point.dot(point! { x: 2.0, y: 4.5, z: 6.3 });
    ///
    /// assert_eq!(dot, 58.8);
    /// ```
    pub fn dot(self, other: Self) -> T {
        self.x() * other.x() + self.y() * other.y() + self.z() * other.z()
    }

    /// Returns the cross product of 3 points. A positive value implies
    /// `self` → `point_b` → `point_c` is counter-clockwise, negative implies
    /// clockwise.
    ///
    /// # Note on Robustness
    ///
    /// This function is **not** robust against floating-point errors.
    /// The [`geo`](https://docs.rs/geo) crate
    /// offers robust predicates for standard numeric types using the
    /// [`Kernel`](https://docs.rs/geo/algorithm/kernels/trait.Kernel.html)
    /// trait, and these should be preferred if possible.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_types::point;
    ///
    /// let point_a = point! { x: 1., y: 2., z: 3. };
    /// let point_b = point! { x: 3., y: 5., z: 3. };
    /// let point_c = point! { x: 7., y: 12., z: 7. };
    ///
    /// let cross = point_a.cross_prod(point_b, point_c);
    ///
    /// assert_eq!(cross, 2.0)
    /// ```
    pub fn cross_prod_2d(self, point_b: Self, point_c: Self) -> T {
        (point_b.x() - self.x()) * (point_c.y() - self.y())
            - (point_b.y() - self.y()) * (point_c.x() - self.x())
    }

    /// Returns the cross product of 3 points in 3D space.
    /// The result is a 3D vector, which is perpendicular to the plane formed
    /// by the vectors `self` → `point_b` and `self` → `point_c`.
    ///
    /// # Note
    ///
    /// - This function is **not robust** against floating-point errors.
    /// - For 2D use cases, use `cross_prod_2d`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_types::point;
    ///
    /// let point_a = point! { x: 1., y: 2., z: 3. };
    /// let point_b = point! { x: 3., y: 5., z: 3. };
    /// let point_c = point! { x: 7., y: 12., z: 7. };
    ///
    /// let cross = point_a.cross_prod(point_b, point_c);
    ///
    /// assert_eq!(cross, point! { x: 0.0, y: 12.0, z: -9.0 });
    /// ```
    pub fn cross_prod(self, point_b: Self, point_c: Self) -> Self {
        let ux = point_b.x() - self.x();
        let uy = point_b.y() - self.y();
        let uz = point_b.z() - self.z();
        
        let vx = point_c.x() - self.x();
        let vy = point_c.y() - self.y();
        let vz = point_c.z() - self.z();

        // Compute the cross product:
        let x = uy * vz - uz * vy;
        let y = uz * vx - ux * vz;
        let z = ux * vy - uy * vx;

        point! { x: x, y: y, z: z }
    }
}

impl<T: CoordNum> Point<T> {
    /// Converts the (x, y, z) components of Point to degrees
    ///
    /// # Example
    /// ```
    /// use geo_types::Point;
    ///
    /// let p = Point::new(1.234, 2.345, 8.591);
    /// let (x, y, z): (f32, f32, f32) = p.to_degrees().x_y_z();
    /// assert_eq!(x.round(), 71.0);
    /// assert_eq!(y.round(), 134.0);
    /// assert_eq!(z.round(), 0.15);
    /// ```
    pub fn to_degrees(self) -> Self {
        let (x, y, z) = self.x_y_z();
        let x = x.to_degrees();
        let y = y.to_degrees();
        let z = z.to_degrees();
        Point::new(x, y, z)
    }

    /// Converts the (x, y, z) components of Point to radians
    ///
    /// # Example
    /// ```
    /// use geo_types::Point;
    ///
    /// let p = Point::new(180.0, 341.5, 0.15);
    /// let (x, y, z): (f32, f32, f32) = p.to_radians().x_y_z();
    /// assert_eq!(x.round(), 3.0);
    /// assert_eq!(y.round(), 6.0);
    /// assert_eq!(z.round(), 8.59);
    /// ```
    pub fn to_radians(self) -> Self {
        let (x, y, z) = self.x_y_z();
        let x = x.to_radians();
        let y = y.to_radians();
        let z = z.to_radians();
        Point::new(x, y, z)
    }
}

impl<T> Neg for Point<T>
where
    T: CoordNum + Neg<Output = T>,
{
    type Output = Self;

    /// Returns a point with the x and y components negated.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_types::Point;
    ///
    /// let p = -Point::new(-1.25, 2.5, -7.3);
    ///
    /// assert_eq!(p.x(), 1.25);
    /// assert_eq!(p.y(), -2.5);
    /// assert_eq!(p.z(), 7.3);
    /// ```
    fn neg(self) -> Self::Output {
        Point::from(-self.0)
    }
}

impl<T: CoordNum> Add for Point<T> {
    type Output = Self;

    /// Add a point to the given point.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_types::Point;
    ///
    /// let p = Point::new(1.25, 2.5, 5.0) + Point::new(1.5, 2.5, 5.0);
    ///
    /// assert_eq!(p.x(), 2.75);
    /// assert_eq!(p.y(), 5.0);
    /// assert_eq!(p.z(), 10.0);
    /// ```
    fn add(self, rhs: Self) -> Self::Output {
        Point::from(self.0 + rhs.0)
    }
}

impl<T: CoordNum> AddAssign for Point<T> {
    /// Add a point to the given point and assign it to the original point.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_types::Point;
    ///
    /// let mut p = Point::new(1.25, 2.5, 5.0);
    /// p += Point::new(1.5, 2.5, 5.0);
    ///
    /// assert_eq!(p.x(), 2.75);
    /// assert_eq!(p.y(), 5.0);
    /// assert_eq!(p.z(), 10.0);
    /// ```
    fn add_assign(&mut self, rhs: Self) {
        self.0 = self.0 + rhs.0;
    }
}

impl<T: CoordNum> Sub for Point<T> {
    type Output = Self;

    /// Subtract a point from the given point.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_types::Point;
    ///
    /// let p = Point::new(1.25, 3.0, 8.0) - Point::new(1.5, 2.5, 3.5);
    ///
    /// assert_eq!(p.x(), -0.25);
    /// assert_eq!(p.y(), 0.5);
    /// assert_eq!(p.z(), 4.5);
    /// ```
    fn sub(self, rhs: Self) -> Self::Output {
        Point::from(self.0 - rhs.0)
    }
}

impl<T: CoordNum> SubAssign for Point<T> {
    /// Subtract a point from the given point and assign it to the original point.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_types::Point;
    ///
    /// let mut p = Point::new(1.25, 2.5, 8.0);
    /// p -= Point::new(1.5, 2.5, 3.5);
    ///
    /// assert_eq!(p.x(), -0.25);
    /// assert_eq!(p.y(), 0.0);
    /// assert_eq!(p.z(), 4.5);
    /// ```
    fn sub_assign(&mut self, rhs: Self) {
        self.0 = self.0 - rhs.0;
    }
}

impl<T: CoordNum> Mul<T> for Point<T> {
    type Output = Self;

    /// Scaler multiplication of a point
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_types::Point;
    ///
    /// let p = Point::new(2.0, 3.0, 5.0) * 2.0;
    ///
    /// assert_eq!(p.x(), 4.0);
    /// assert_eq!(p.y(), 6.0);
    /// assert_eq!(p.z(), 10.0);
    /// ```
    fn mul(self, rhs: T) -> Self::Output {
        Point::from(self.0 * rhs)
    }
}

impl<T: CoordNum> MulAssign<T> for Point<T> {
    /// Scaler multiplication of a point in place
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_types::Point;
    ///
    /// let mut p = Point::new(2.0, 3.0, 5.0);
    /// p *= 2.0;
    ///
    /// assert_eq!(p.x(), 4.0);
    /// assert_eq!(p.y(), 6.0);
    /// assert_eq!(p.z(), 10.0);
    /// ```
    fn mul_assign(&mut self, rhs: T) {
        self.0 = self.0 * rhs
    }
}

impl<T: CoordNum> Div<T> for Point<T> {
    type Output = Self;

    /// Scaler division of a point
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_types::Point;
    ///
    /// let p = Point::new(2.0, 3.0, 5.0) / 2.0;
    ///
    /// assert_eq!(p.x(), 1.0);
    /// assert_eq!(p.y(), 1.5);
    /// assert_eq!(p.z(), 2.5);
    /// ```
    fn div(self, rhs: T) -> Self::Output {
        Point::from(self.0 / rhs)
    }
}

impl<T: CoordNum> DivAssign<T> for Point<T> {
    /// Scaler division of a point in place
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_types::Point;
    ///
    /// let mut p = Point::new(2.0, 3.0, 10.0);
    /// p /= 2.0;
    ///
    /// assert_eq!(p.x(), 1.0);
    /// assert_eq!(p.y(), 1.5);
    /// assert_eq!(p.z(), 2.5);
    /// ```
    fn div_assign(&mut self, rhs: T) {
        self.0 = self.0 / rhs
    }
}

#[cfg(any(feature = "approx", test))]
impl<T> RelativeEq for Point<T>
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
    /// use geo_types::Point;
    ///
    /// let a = Point::new(2.0, 3.0, 1.01);
    /// let b = Point::new(2.0, 3.01, 1.02);
    ///
    /// approx::assert_relative_eq!(a, b, max_relative=0.1)
    /// ```
    #[inline]
    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        self.0.relative_eq(&other.0, epsilon, max_relative)
    }
}

#[cfg(any(feature = "approx", test))]
impl<T> AbsDiffEq for Point<T>
where
    T: AbsDiffEq<Epsilon = T> + CoordNum,
    T::Epsilon: Copy,
{
    type Epsilon = T::Epsilon;

    #[inline]
    fn default_epsilon() -> Self::Epsilon {
        T::default_epsilon()
    }

    /// Equality assertion with an absolute limit.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_types::Point;
    ///
    /// let a = Point::new(2.0, 3.0, 4.0);
    /// let b = Point::new(2.0, 3.0000001, 4.1);
    ///
    /// approx::assert_relative_eq!(a, b, epsilon=0.1)
    /// ```
    #[inline]
    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.0.abs_diff_eq(&other.0, epsilon)
    }
}

#[cfg(feature = "rstar")]
impl<T> ::rstar::Point for Point<T>
where
    T: ::num_traits::Float + ::rstar::RTreeNum,
{
    type Scalar = T;

    const DIMENSIONS: usize = 3;

    fn generate(mut generator: impl FnMut(usize) -> Self::Scalar) -> Self {
        Point::new(generator(0), generator(1), generator(2))
    }

    fn nth(&self, index: usize) -> Self::Scalar {
        match index {
            0 => self.0.x,
            1 => self.0.y,
            2 => self.0.z,
            _ => unreachable!(),
        }
    }
    fn nth_mut(&mut self, index: usize) -> &mut Self::Scalar {
        match index {
            0 => &mut self.0.x,
            1 => &mut self.0.y,
            2 => &mut self.0.z,
            _ => unreachable!(),
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use approx::AbsDiffEq;

    #[test]
    fn test_abs_diff_eq() {
        let delta = 1e-6;
        let p = Point::new(1.0, 1.0, 1.0);

        let p_x = Point::new(1.0 - delta, 1.0, 1.0);
        assert!(p.abs_diff_eq(&p_x, 1e-2));
        assert!(p.abs_diff_ne(&p_x, 1e-12));

        let p_y = Point::new(1.0, 1.0 + delta, 1.0);
        assert!(p.abs_diff_eq(&p_y, 1e-2));
        assert!(p.abs_diff_ne(&p_y, 1e-12));

        let p_z = Point::new(1.0, 1.0, 1.0 + delta);
        assert!(p.abs_diff_eq(&p_z, 1e-2));
        assert!(p.abs_diff_ne(&p_z, 1e-12));

        let p_xy = Point::new(1.0 + delta, 1.0 - delta, 1.0);
        assert!(p.abs_diff_eq(&p_xy, 1e-2));
        assert!(p.abs_diff_ne(&p_xy, 1e-12));

        let p_inf = Point::new(f64::INFINITY, 1., 1.);
        assert!(p.abs_diff_ne(&p_inf, 1e-2));
    }

    #[test]
    fn test_relative_eq() {
        let delta = 1e-6;
        let p = Point::new(1.0, 1.0, 1.0);

        let p_x = Point::new(1.0 - delta, 1.0, 1.0);
        assert!(p.relative_eq(&p_x, 1e-2, 1e-2));
        assert!(p.relative_ne(&p_x, 1e-12, 1e-12));

        let p_y = Point::new(1.0, 1.0 + delta, 1.0);
        assert!(p.relative_eq(&p_y, 1e-2, 1e-2));
        assert!(p.relative_ne(&p_y, 1e-12, 1e-12));

        let p_z = Point::new(1.0, 1.0, 1.0 + delta);
        assert!(p.relative_eq(&p_z, 1e-2, 1e-2));
        assert!(p.relative_ne(&p_z, 1e-12, 1e-12));

        let p_xy = Point::new(1.0 + delta, 1.0 - delta, 1.0);
        assert!(p.relative_eq(&p_xy, 1e-2, 1e-2));
        assert!(p.relative_ne(&p_xy, 1e-12, 1e-12));

        let p_inf = Point::new(f64::INFINITY, 1., 1.);
        assert!(p.relative_ne(&p_inf, 1e-2, 1e-2));
    }
}
