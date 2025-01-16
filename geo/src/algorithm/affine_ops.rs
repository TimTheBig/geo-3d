use num_traits::ToPrimitive;

use crate::{Coord, CoordFloat, CoordNum, MapCoords, MapCoordsInPlace};
use std::{fmt, ops::Mul, ops::Neg};

/// Apply an [`AffineTransform`] like [`scale`](AffineTransform::scale),
/// [`skew`](AffineTransform::skew), or [`rotate`](AffineTransform::rotate) to a
/// [`Geometry`](crate::geometry::Geometry).
///
/// Multiple transformations can be composed in order to be efficiently applied in a single
/// operation. See [`AffineTransform`] for more on how to build up a transformation.
///
/// If you are not composing operations, traits that leverage this same machinery exist which might
/// be more readable. See: [`Scale`](crate::algorithm::Scale),
/// [`Translate`](crate::algorithm::Translate), [`Rotate`](crate::algorithm::Rotate),
/// and [`Skew`](crate::algorithm::Skew).
///
/// # Examples
/// ## Build up transforms by beginning with a constructor, then chaining mutation operations
/// ```
/// use geo::{AffineOps, AffineTransform};
/// use geo::{point, line_string, BoundingRect};
/// use approx::assert_relative_eq;
///
/// let line_string = line_string![(x: 0.0, y: 0.0, z: 0.0),(x: 1.0, y: 1.0, z: 1.0)];
///
/// let transform = AffineTransform::translate(1.0, 1.0).scaled(2.0, 2.0, point!(x: 0.0, y: 0.0, z: 0.0));
///
/// let transformed_line_string = line_string.affine_transform(&transform);
///
/// assert_relative_eq!(
///     transformed_line_string,
///     line_string![(x: 2.0, y: 2.0, z: 2.0),(x: 4.0, y: 4.0, z: 4.0)]
/// );
/// ```
pub trait AffineOps<T: CoordNum> {
    /// Apply `transform` immutably, outputting a new geometry.
    #[must_use]
    fn affine_transform(&self, transform: &AffineTransform<T>) -> Self;

    /// Apply `transform` to mutate `self`.
    fn affine_transform_mut(&mut self, transform: &AffineTransform<T>);
}

impl<T: CoordNum, M: MapCoordsInPlace<T> + MapCoords<T, T, Output = Self>> AffineOps<T> for M {
    fn affine_transform(&self, transform: &AffineTransform<T>) -> Self {
        self.map_coords(|c| transform.apply(c))
    }

    fn affine_transform_mut(&mut self, transform: &AffineTransform<T>) {
        self.map_coords_in_place(|c| transform.apply(c))
    }
}

/// A general affine transformation matrix, and associated operations.
///
/// Note that affine ops are **already implemented** on most `geo-types` primitives, using this module.
///
/// Affine transforms using the same numeric type (e.g. [`CoordFloat`]) can be **composed**,
/// and the result can be applied to geometries using e.g. [`MapCoords`]. This allows the
/// efficient application of transforms: an arbitrary number of operations can be chained.
/// These are then composed, producing a final transformation matrix which is applied to the geometry coordinates.
///
/// `AffineTransform` is a row-major matrix.
/// 2D affine transforms require six matrix parameters:
///
/// `[a, b, xoff, d, e, yoff]`
///
/// these map onto the `AffineTransform` rows as follows:
/// ```ignore
/// [[a, b, xoff],
/// [d, e, yoff],
/// [0, 0, 1]]
/// ```
/// The equations for transforming coordinates `(x, y) -> (x', y')` are given as follows:
///
/// `x' = ax + by + xoff`
///
/// `y' = dx + ey + yoff`
///
/// # Usage
///
/// Two types of operation are provided: construction and mutation. **Construction** functions create a *new* transform
/// and are denoted by the use of the **present tense**: `scale()`, `translate()`, `rotate()`, and `skew()`.
///
/// **Mutation** methods *add* a transform to the existing `AffineTransform`, and are denoted by the use of the past participle:
/// `scaled()`, `translated()`, `rotated()`, and `skewed()`.
///
/// # Examples
/// ## Build up transforms by beginning with a constructor, then chaining mutation operations
/// ```
/// use geo::{AffineOps, AffineTransform};
/// use geo::{point, line_string, BoundingRect};
/// use approx::assert_relative_eq;
///
/// let line_string = line_string![(x: 0.0, y: 0.0),(x: 1.0, y: 1.0)];
///
/// let transform = AffineTransform::translate(1.0, 1.0).scaled(2.0, 2.0, point!(x: 0.0, y: 0.0));
///
/// let transformed_line_string = line_string.affine_transform(&transform);
///
/// assert_relative_eq!(
///     transformed_line_string,
///     line_string![(x: 2.0, y: 2.0),(x: 4.0, y: 4.0)]
/// );
/// ```
///
/// ## Create affine transform manually, and access elements using getter methods
/// ```
/// use geo::AffineTransform;
///
/// let transform = AffineTransform::new(10.0, 0.0, 400_000.0, 0.0, -10.0, 500_000.0);
///
/// let a: f64 = transform.a();
/// let b: f64 = transform.b();
/// let xoff: f64 = transform.xoff();
/// let d: f64 = transform.d();
/// let e: f64 = transform.e();
/// let yoff: f64 = transform.yoff();
/// assert_eq!(transform, AffineTransform::new(a, b, xoff, d, e, yoff))
/// ```

#[derive(Copy, Clone, PartialEq, Eq)]
pub struct AffineTransform<T: CoordNum = f64>([[T; 4]; 4]);

impl<T: CoordNum> Default for AffineTransform<T> {
    fn default() -> Self {
        // identity matrix
        Self::identity()
    }
}

/// Makes an fn for a field at specified pos
macro_rules! affine_field {
    ($field_name: ident[$pos_0: literal][$pos_1: literal]) => {
        /// See [AffineTransform::new] for this value's role in the affine transformation.
        pub fn $field_name(&self) -> T {
            self.0[$pos_0][$pos_1]
        }
    };
}

impl<T: CoordNum> AffineTransform<T> {
    /// Create a new affine transformation by composing two `AffineTransform`s.
    ///
    /// This is a **cumulative** operation; the new transform is *added* to the existing transform.
    #[must_use]
    pub fn compose(&self, other: &Self) -> Self {
        Self([
            [
                (other.0[0][0] * self.0[0][0])
                    + (other.0[0][1] * self.0[1][0])
                    + (other.0[0][2] * self.0[2][0])
                    + (other.0[0][3] * self.0[3][0]),
                (other.0[0][0] * self.0[0][1])
                    + (other.0[0][1] * self.0[1][1])
                    + (other.0[0][2] * self.0[2][1])
                    + (other.0[0][3] * self.0[3][1]),
                (other.0[0][0] * self.0[0][2])
                    + (other.0[0][1] * self.0[1][2])
                    + (other.0[0][2] * self.0[2][2])
                    + (other.0[0][3] * self.0[3][2]),
                (other.0[0][0] * self.0[0][3])
                    + (other.0[0][1] * self.0[1][3])
                    + (other.0[0][2] * self.0[2][3])
                    + (other.0[0][3] * self.0[3][3]),
            ],
            [
                (other.0[1][0] * self.0[0][0])
                    + (other.0[1][1] * self.0[1][0])
                    + (other.0[1][2] * self.0[2][0])
                    + (other.0[1][3] * self.0[3][0]),
                (other.0[1][0] * self.0[0][1])
                    + (other.0[1][1] * self.0[1][1])
                    + (other.0[1][2] * self.0[2][1])
                    + (other.0[1][3] * self.0[3][1]),
                (other.0[1][0] * self.0[0][2])
                    + (other.0[1][1] * self.0[1][2])
                    + (other.0[1][2] * self.0[2][2])
                    + (other.0[1][3] * self.0[3][2]),
                (other.0[1][0] * self.0[0][3])
                    + (other.0[1][1] * self.0[1][3])
                    + (other.0[1][2] * self.0[2][3])
                    + (other.0[1][3] * self.0[3][3]),
            ],
            [
                (other.0[2][0] * self.0[0][0])
                    + (other.0[2][1] * self.0[1][0])
                    + (other.0[2][2] * self.0[2][0])
                    + (other.0[2][3] * self.0[3][0]),
                (other.0[2][0] * self.0[0][1])
                    + (other.0[2][1] * self.0[1][1])
                    + (other.0[2][2] * self.0[2][1])
                    + (other.0[2][3] * self.0[3][1]),
                (other.0[2][0] * self.0[0][2])
                    + (other.0[2][1] * self.0[1][2])
                    + (other.0[2][2] * self.0[2][2])
                    + (other.0[2][3] * self.0[3][2]),
                (other.0[2][0] * self.0[0][3])
                    + (other.0[2][1] * self.0[1][3])
                    + (other.0[2][2] * self.0[2][3])
                    + (other.0[2][3] * self.0[3][3]),
            ],
            [
                // this section isn't technically necessary since the last row is invariant: [0, 0, 0, 1]
                (other.0[3][0] * self.0[0][0])
                    + (other.0[3][1] * self.0[1][0])
                    + (other.0[3][2] * self.0[2][0])
                    + (other.0[3][3] * self.0[3][0]),
                (other.0[3][0] * self.0[0][1])
                    + (other.0[3][1] * self.0[1][1])
                    + (other.0[3][2] * self.0[2][1])
                    + (other.0[3][3] * self.0[3][1]),
                (other.0[3][0] * self.0[0][2])
                    + (other.0[3][1] * self.0[1][2])
                    + (other.0[3][2] * self.0[2][2])
                    + (other.0[3][3] * self.0[3][2]),
                (other.0[3][0] * self.0[0][3])
                    + (other.0[3][1] * self.0[1][3])
                    + (other.0[3][2] * self.0[2][3])
                    + (other.0[3][3] * self.0[3][3]),
            ],
        ])
    }

    /// Create a new affine transformation by composing an arbitrary number of `AffineTransform`s.
    ///
    /// This is a **cumulative** operation; the new transform is *added* to the existing transform.
    /// ```
    /// use geo::AffineTransform;
    /// let mut transform = AffineTransform::identity();
    ///
    /// // create two transforms that cancel each other
    /// let transform1 = AffineTransform::translate(1.0, 2.0);
    /// let transform2 = AffineTransform::translate(-1.0, -2.0);
    /// let transforms = vec![transform1, transform2];
    ///
    /// // apply them
    /// let outcome = transform.compose_many(&transforms);
    /// // we should be back to square one
    /// assert!(outcome.is_identity());
    /// ```
    #[must_use]
    pub fn compose_many(&self, transforms: &[Self]) -> Self {
        self.compose(&transforms.iter().fold(
            AffineTransform::default(),
            |acc: AffineTransform<T>, transform| acc.compose(transform),
        ))
    }

    /// Create the identity matrix
    ///
    /// The matrix is:
    /// ```ignore
    /// [[1, 0, 0, 0],
    ///  [0, 1, 0, 0],
    ///  [0, 0, 1, 0],
    ///  [0, 0, 0, 1]]
    /// ```
    pub fn identity() -> Self {
        Self([
            [T::one(), T::zero(), T::zero(), T::zero()],
            [T::zero(), T::one(), T::zero(), T::zero()],
            [T::zero(), T::zero(), T::one(), T::zero()],
            [T::zero(), T::zero(), T::zero(), T::one()],
        ])
    }

    /// Whether the transformation is equivalent to the [identity matrix](Self::identity),
    /// that is, whether it's application will be a a no-op.
    ///
    /// ```
    /// use geo::AffineTransform;
    /// let mut transform = AffineTransform::identity();
    /// assert!(transform.is_identity());
    ///
    /// // mutate the transform a bit
    /// transform = transform.translated(1.0, 2.0, 3.0);
    /// assert!(!transform.is_identity());
    ///
    /// // put it back
    /// transform = transform.translated(-1.0, -2.0, -3.0);
    /// assert!(transform.is_identity());
    /// ```
    pub fn is_identity(&self) -> bool {
        self == &Self::identity()
    }

    /// **Create** a new affine transform for scaling, scaled by factors along the `x` and `y` dimensions.
    /// The point of origin is *usually* given as the 2D bounding box centre of the geometry, but
    /// any coordinate may be specified.
    /// Negative scale factors will mirror or reflect coordinates.
    ///
    /// The matrix is:
    /// ```ignore
    /// [[xfact, 0, xoff],
    /// [0, yfact, yoff],
    /// [0, 0, 1]]
    ///
    /// xoff = origin.x - (origin.x * xfact)
    /// yoff = origin.y - (origin.y * yfact)
    /// ```
    /// Create a new affine transform for scaling in 3D
    pub fn scale(xfact: T, yfact: T, zfact: T, origin: impl Into<Coord<T>>) -> Self {
        let (x0, y0, z0) = origin.into().x_y_z();
        let xoff = x0 - x0 * xfact;
        let yoff = y0 - y0 * yfact;
        let zoff = z0 - z0 * zfact;

        Self::new(
            xfact, T::zero(), T::zero(), xoff,
            T::zero(), yfact, T::zero(), yoff,
            T::zero(), T::zero(), zfact, zoff
        )
    }

    /// **Add** an affine transform for scaling, scaled by factors along the `x` and `y` dimensions.
    /// The point of origin is *usually* given as the 2D bounding box centre of the geometry, but
    /// any coordinate may be specified.
    /// Negative scale factors will mirror or reflect coordinates.
    /// This is a **cumulative** operation; the new transform is *added* to the existing transform.
    #[must_use]
    pub fn scaled(mut self, xfact: T, yfact: T, zfact: T, origin: impl Into<Coord<T>>) -> Self {
        self.0 = self.compose(&Self::scale(xfact, yfact, zfact, origin)).0;
        self
    }

    /// **Create** an affine transform for translation, shifted by offsets along the `x` and `y` dimensions.
    ///
    /// The matrix is:
    /// ```ignore
    /// [[1, 0, 0, xoff],
    /// [0, 1, 0, yoff],
    /// [0, 0, 1, zoff],
    /// [0, 0, 0, 1]]
    /// ```
    pub fn translate(xoff: T, yoff: T, zoff: T) -> Self {
        Self::new(
            T::one(), T::zero(), T::zero(), xoff,
            T::zero(), T::one(), T::zero(), yoff,
            T::zero(), T::zero(), T::one(), zoff
        )
    }

    /// **Add** an affine transform for translation, shifted by offsets along the `x` and `y` dimensions
    ///
    /// This is a **cumulative** operation; the new transform is *added* to the existing transform.
    #[must_use]
    pub fn translated(mut self, xoff: T, yoff: T, zoff: T) -> Self {
        self.0 = self.compose(&Self::translate(xoff, yoff, zoff)).0;
        self
    }

    /// Apply the current transform to a coordinate
    pub fn apply(&self, coord: Coord<T>) -> Coord<T> {
        Coord {
            x: self.0[0][0] * coord.x + self.0[0][1] * coord.y + self.0[0][2] * coord.z + self.0[0][3],
            y: self.0[1][0] * coord.x + self.0[1][1] * coord.y + self.0[1][2] * coord.z + self.0[1][3],
            z: self.0[2][0] * coord.x + self.0[2][1] * coord.y + self.0[2][2] * coord.z + self.0[2][3],
        }
    }

    /// Create a new custom transform matrix
    ///
    /// The argument order matches that of the affine transform matrix:
    ///```ignore
    /// [[a, b, f, xoff],
    ///  [d, e, g, yoff],
    ///  [h, i, j, zoff],
    ///  [0, 0, 0, 1]] <-- not part of the input arguments
    /// ```
    pub fn new(a: T, b: T, f: T, xoff: T, d: T, e: T, g: T, yoff: T, h: T, i: T, j: T, zoff: T) -> Self {
        Self([
            [a, b, f, xoff],
            [d, e, g, yoff],
            [h, i, j, zoff],
            [T::zero(), T::zero(), T::zero(), T::one()]
        ])
    }

    // See [AffineTransform::new] for these value's roles in the affine transformation.
    affine_field!(a[0][0]);
    affine_field!(b[0][1]);
    affine_field!(f[0][2]);
    affine_field!(xoff[0][3]);

    affine_field!(d[1][0]);
    affine_field!(e[1][1]);
    affine_field!(g[1][2]);
    affine_field!(yoff[1][3]);

    affine_field!(h[2][0]);
    affine_field!(i[2][1]);
    affine_field!(j[2][2]);
    affine_field!(zoff[2][3]);
}

impl<T: CoordNum + Neg> AffineTransform<T> {
    /// Return the inverse of a given transform. Composing a transform with its inverse yields
    /// the [identity matrix](Self::identity)
    #[must_use]
    pub fn inverse(&self) -> Option<Self>
    where
        <T as Neg>::Output: Mul<T>,
        <<T as Neg>::Output as Mul<T>>::Output: ToPrimitive,
    {
        let a = self.0[0][0];
        let b = self.0[0][1];
        let xoff = self.0[0][2];
        let d = self.0[1][0];
        let e = self.0[1][1];
        let yoff = self.0[1][2];

        let determinant = a * e - b * d;

        if determinant == T::zero() {
            return None; // The matrix is not invertible
        }
        let inv_det = T::one() / determinant;

        // If conversion of either the b or d matrix value fails, bail out
        Some(Self::new(
            e * inv_det,
            T::from(-b * inv_det)?,
            (b * yoff - e * xoff) * inv_det,
            T::from(-d * inv_det)?,
            a * inv_det,
            (d * xoff - a * yoff) * inv_det,
        ))
    }
}

impl<T: CoordNum> fmt::Debug for AffineTransform<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("AffineTransform")
            .field("a", &self.0[0][0])
            .field("b", &self.0[0][1])
            .field("f", &self.0[0][2])
            .field("xoff", &self.0[0][3])
            .field("d", &self.0[1][0])
            .field("e", &self.0[1][1])
            .field("g", &self.0[1][2])
            .field("yoff", &self.0[1][3])
            .field("h", &self.0[2][0])
            .field("i", &self.0[2][1])
            .field("j", &self.0[2][2])
            .field("zoff", &self.0[2][3])
            .finish()
    }
}

impl<T: CoordNum> From<[T; 12]> for AffineTransform<T> {
    fn from(arr: [T; 12]) -> Self {
        Self::new(
            arr[0], arr[1], arr[2], arr[3],
            arr[4], arr[5], arr[6], arr[7],
            arr[8], arr[9], arr[10], arr[11]
        )
    }
}

impl<T: CoordNum> From<(T, T, T, T, T, T, T, T, T, T, T, T)> for AffineTransform<T> {
    fn from(tup: (T, T, T, T, T, T, T, T, T, T, T, T)) -> Self {
        Self::new(
            tup.0, tup.1, tup.2, tup.3,
            tup.4, tup.5, tup.6, tup.7,
            tup.8, tup.9, tup.10, tup.11
        )
    }
}


impl<U: CoordFloat> AffineTransform<U> {
    /// **Create** an affine transform for rotation, using an arbitrary point as its center.
    ///
    /// Note that this operation is only available for geometries with floating point coordinates.
    ///
    /// `angle` is given in **degrees**.
    ///
    /// The matrix (angle denoted as theta) is:
    /// ```ignore
    /// [[cos_theta, -sin_theta, 0, xoff],
    ///  [sin_theta,  cos_theta, 0, yoff],
    ///  [0,          0,         1, zoff],
    ///  [0,          0,         0, 1]]
    ///
    /// xoff = origin.x - (origin.x * cos(theta)) + (origin.y * sin(theta))
    /// yoff = origin.y - (origin.x * sin(theta)) - (origin.y * cos(theta))
    /// zoff = origin.z
    /// ```
    pub fn rotate(degrees: U, origin: impl Into<Coord<U>>) -> Self {
        let (sin_theta, cos_theta) = degrees.to_radians().sin_cos();
        let Coord { x: x0, y: y0, z: z0 } = origin.into();
        let xoff = x0 - (x0 * cos_theta) + (y0 * sin_theta);
        let yoff = y0 - (x0 * sin_theta) - (y0 * cos_theta);
        let zoff = z0; // No change in the z-axis for a 2D rotation in 3D space.

        Self::new(
            cos_theta, -sin_theta, U::zero(), xoff,
            sin_theta,  cos_theta, U::zero(), yoff,
            U::zero(),  U::zero(), U::one(), zoff,
        )
    }

    /// **Add** an affine transform for rotation, using an arbitrary point as its centre.
    ///
    /// Note that this operation is only available for geometries with floating point coordinates.
    ///
    /// `angle` is given in **degrees**.
    ///
    /// This is a **cumulative** operation; the new transform is *added* to the existing transform.
    #[must_use]
    pub fn rotated(mut self, angle: U, origin: impl Into<Coord<U>>) -> Self {
        self.0 = self.compose(&Self::rotate(angle, origin)).0;
        self
    }

    /// **Create** an affine transform for skewing.
    ///
    /// Note that this operation is only available for geometries with floating point coordinates.
    ///
    /// Geometries are sheared by angles along x (`xs`) and y (`ys`) dimensions.
    /// The point of origin is *usually* given as the 2D bounding box centre of the geometry, but
    /// any coordinate may be specified. Angles are given in **degrees**.
    /// The matrix is:
    /// ```ignore
    /// [[1, tan(x), 0, xoff],
    /// [tan(y), 1, 0, yoff],
    /// [0, 0, 1, zoff],
    /// [0, 0, 0, 1]]
    ///
    /// xoff = -origin.y * tan(xs)
    /// yoff = -origin.x * tan(ys)
    /// zoff = origin.z
    /// ```
    pub fn skew(xs: U, ys: U, origin: impl Into<Coord<U>>) -> Self {
        let Coord { x: x0, y: y0, z: z0 } = origin.into();
        let mut tanx = xs.to_radians().tan();
        let mut tany = ys.to_radians().tan();
        // These checks are stolen from Shapely's implementation -- may not be necessary
        if tanx.abs() < U::from::<f64>(2.5e-16).unwrap() {
            tanx = U::zero();
        }
        if tany.abs() < U::from::<f64>(2.5e-16).unwrap() {
            tany = U::zero();
        }
        let xoff = -y0 * tanx;
        let yoff = -x0 * tany;
        let zoff = z0; // Z-offset remains unchanged during a 2D skew operation.
        Self::new(U::one(), tanx, U::zero(), xoff, tany, U::one(), U::zero(), yoff, U::zero(), U::zero(), U::one(), zoff)
    }

    /// **Add** an affine transform for skewing.
    ///
    /// Note that this operation is only available for geometries with floating point coordinates.
    ///
    /// Geometries are sheared by angles along x (`xs`) and y (`ys`) dimensions.
    /// The point of origin is *usually* given as the 2D bounding box centre of the geometry, but
    /// any coordinate may be specified. Angles are given in **degrees**.
    ///
    /// This is a **cumulative** operation; the new transform is *added* to the existing transform.
    #[must_use]
    pub fn skewed(mut self, xs: U, ys: U, origin: impl Into<Coord<U>>) -> Self {
        self.0 = self.compose(&Self::skew(xs, ys, origin)).0;
        self
    }
}

#[cfg(test)]
mod tests {
    use approx::{AbsDiffEq, RelativeEq};

    impl<T> RelativeEq for AffineTransform<T>
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
        /// use geo_types::AffineTransform;
        /// use geo_types::point;
        ///
        /// let a = AffineTransform::new(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
        /// let b = AffineTransform::new(1.01, 2.02, 3.03, 4.04, 5.05, 6.06);
        ///
        /// approx::assert_relative_eq!(a, b, max_relative=0.1)
        /// approx::assert_relative_ne!(a, b, max_relative=0.055)
        /// ```
        #[inline]
        fn relative_eq(
            &self,
            other: &Self,
            epsilon: Self::Epsilon,
            max_relative: Self::Epsilon,
        ) -> bool {
            let mut mp_zipper = self.0.iter().flatten().zip(other.0.iter().flatten());
            mp_zipper.all(|(lhs, rhs)| lhs.relative_eq(rhs, epsilon, max_relative))
        }
    }

    impl<T> AbsDiffEq for AffineTransform<T>
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
        /// use geo_types::MultiPoint;
        /// use geo_types::point;
        ///
        /// let a = AffineTransform::new(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
        /// let b = AffineTransform::new(1.01, 2.02, 3.03, 4.04, 5.05, 6.06);
        ///
        /// approx::abs_diff_eq!(a, b, epsilon=0.1)
        /// approx::abs_diff_ne!(a, b, epsilon=0.055)
        /// ```
        #[inline]
        fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
            let mut mp_zipper = self.0.iter().flatten().zip(other.0.iter().flatten());
            mp_zipper.all(|(lhs, rhs)| lhs.abs_diff_eq(rhs, epsilon))
        }
    }

    use super::*;
    use crate::{wkt, Point};

    // given a matrix with the shape
    // [[a, b, xoff],
    // [d, e, yoff],
    // [0, 0, 1]]
    #[test]
    fn matrix_multiply() {
        let a = AffineTransform::new(1, 2, 5, 3, 4, 6);
        let b = AffineTransform::new(7, 8, 11, 9, 10, 12);
        let composed = a.compose(&b);
        assert_eq!(composed.0[0][0], 31);
        assert_eq!(composed.0[0][1], 46);
        assert_eq!(composed.0[0][2], 94);
        assert_eq!(composed.0[1][0], 39);
        assert_eq!(composed.0[1][1], 58);
        assert_eq!(composed.0[1][2], 117);
    }
    #[test]
    fn test_transform_composition() {
        let p0 = Point::new(0.0f64, 0.0);
        // scale once
        let mut scale_a = AffineTransform::default().scaled(2.0, 2.0, p0);
        // rotate
        scale_a = scale_a.rotated(45.0, p0);
        // rotate back
        scale_a = scale_a.rotated(-45.0, p0);
        // scale up again, doubling
        scale_a = scale_a.scaled(2.0, 2.0, p0);
        // scaled once
        let scale_b = AffineTransform::default().scaled(2.0, 2.0, p0);
        // scaled once, but equal to 2 + 2
        let scale_c = AffineTransform::default().scaled(4.0, 4.0, p0);
        assert_ne!(&scale_a.0, &scale_b.0);
        assert_relative_eq!(&scale_a, &scale_c);
    }

    #[test]
    fn affine_transformed() {
        let transform = AffineTransform::translate(1.0, 1.0).scaled(2.0, 2.0, (0.0, 0.0));
        let mut poly = wkt! { POLYGON((0.0 0.0,0.0 2.0,1.0 2.0)) };
        poly.affine_transform_mut(&transform);

        let expected = wkt! { POLYGON((2.0 2.0,2.0 6.0,4.0 6.0)) };
        assert_eq!(expected, poly);
    }
    #[test]
    fn affine_transformed_inverse() {
        let transform = AffineTransform::translate(1.0, 1.0).scaled(2.0, 2.0, (0.0, 0.0));
        let tinv = transform.inverse().unwrap();
        let identity = transform.compose(&tinv);
        // test really only needs this, but let's be sure
        assert!(identity.is_identity());

        let mut poly = wkt! { POLYGON((0.0 0.0,0.0 2.0,1.0 2.0)) };
        let expected = poly.clone();
        poly.affine_transform_mut(&identity);
        assert_eq!(expected, poly);
    }
    #[test]
    fn test_affine_transform_getters() {
        let transform = AffineTransform::new(10.0, 0.0, 400_000.0, 0.0, -10.0, 500_000.0);
        assert_eq!(transform.a(), 10.0);
        assert_eq!(transform.b(), 0.0);
        assert_eq!(transform.xoff(), 400_000.0);
        assert_eq!(transform.d(), 0.0);
        assert_eq!(transform.e(), -10.0);
        assert_eq!(transform.yoff(), 500_000.0);
    }
    #[test]
    fn test_compose() {
        let point = Point::new(1., 0.);

        let translate = AffineTransform::translate(1., 0.);
        let scale = AffineTransform::scale(4., 1., [0., 0.]);
        let composed = translate.compose(&scale);

        assert_eq!(point.affine_transform(&translate), Point::new(2., 0.));
        assert_eq!(point.affine_transform(&scale), Point::new(4., 0.));
        assert_eq!(
            point.affine_transform(&translate).affine_transform(&scale),
            Point::new(8., 0.)
        );

        assert_eq!(point.affine_transform(&composed), Point::new(8., 0.));
    }
}
