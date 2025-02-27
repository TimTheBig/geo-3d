//! This module defines the [Vector3DOps] trait and implements it for the
//! [Coord] struct.

use crate::{Coord, CoordNum};

/// Defines vector operations for 3D coordinate types which implement CoordNum
///
/// This trait is intended for internal use within the geo crate as a way to
/// bring together the various hand-crafted linear algebra operations used
/// throughout other algorithms and attached to various structs.
pub trait Vector3DOps<Rhs = Self>
where
    Self: Sized,
{
    type Scalar: CoordNum;

    /// The euclidean distance between this coordinate and the origin
    ///
    /// `sqrt(x² + y² + z²)`
    ///
    fn magnitude(self) -> Self::Scalar;

    /// The squared distance between this coordinate and the origin.
    /// (Avoids the square root calculation when it is not needed)
    ///
    /// `x² + y² + z²`
    ///
    fn magnitude_squared(self) -> Self::Scalar;

    /// Rotate this coordinate around the origin by 90 degrees clockwise.
    ///
    /// `a.left() => (-a.y, a.x)`
    ///
    /// Assumes a coordinate system where positive `y` is up and positive `x` is
    /// to the right. The described rotation direction is consistent with the
    /// documentation for [crate::algorithm::rotate::Rotate].
    /// This assumes a 2D rotation in the XY plane, leaving Z unchanged.
    fn left(self) -> Self;

    /// Rotate this coordinate around the origin by 90 degrees anti-clockwise.
    ///
    /// `a.right() => (a.y, -a.x)`
    ///
    /// Assumes a coordinate system where positive `y` is up and positive `x` is
    /// to the right. The described rotation direction is consistent with the
    /// documentation for [crate::algorithm::rotate::Rotate].
    /// This assumes a 2D rotation in the XY plane, leaving Z unchanged.
    fn right(self) -> Self;

    /// The calculates the `wedge product` between two vectors.\
    /// Note: This is 2D only
    ///
    /// `a ∧ b = a.x * b.y - a.y * b.x`
    ///
    /// Also known as:
    ///
    ///  - `exterior product`
    ///    - because the wedge product comes from 'Exterior Algebra'
    ///  - `perpendicular product`
    ///    -  because it is equivalent to `a.dot(b.right())`
    ///  - `2D cross product`
    ///    - because it is equivalent to the signed magnitude of the
    ///      conventional 3D cross product assuming `z` ordinates are zero
    ///  - `determinant`
    ///    - because it is equivalent to the `determinant` of the 2x2 matrix
    ///      formed by the column-vector inputs.
    ///
    /// ## Examples
    ///
    /// The following list highlights some examples in geo which might be
    /// brought together to use this function:
    ///
    /// 1. [geo_types::Point::cross_prod()] is already defined on
    ///    [geo_types::Point]... but that it seems to be some other
    ///    operation on 3 points??
    /// 2. [geo_types::Line] struct also has a [geo_types::Line::determinant()]
    ///    function which is the same as `line.start.wedge_product(line.end)`
    /// 3. The [crate::algorithm::Kernel::orient2d()] trait default
    ///    implementation uses cross product to compute orientation. It returns
    ///    an enum, not the numeric value which is needed for line segment
    ///    intersection.
    ///
    /// ## Properties
    ///
    /// - The absolute value of the cross product is the area of the
    ///   parallelogram formed by the operands
    /// - Anti-commutative: The sign of the output is reversed if the operands
    ///   are reversed
    /// - If the operands are collinear with the origin, the value is zero
    /// - The sign can be used to check if the operands are clockwise with
    ///   respect to the origin, or phrased differently:
    ///   "is a to the left of the line between the origin and b"?
    ///   - If this is what you are using it for, then please use
    ///     [crate::algorithm::Kernel::orient2d()] instead as this is more
    ///     explicit and has a `RobustKernel` option for extra precision.
    fn wedge_product(self, other: Rhs) -> Self::Scalar;

    /// Try to find a vector of unit length in the same direction as this
    /// vector.
    ///
    /// Returns `None` if the result is not finite. This can happen when
    ///
    /// - the vector is really small (or zero length) and the `.magnitude()`
    ///   calculation has rounded-down to `0.0`
    /// - the vector is really large and the `.magnitude()` has rounded-up
    ///   or 'overflowed' to `f64::INFINITY`
    /// - Either x or y are `f64::NAN` or `f64::INFINITY`
    fn try_normalize(self) -> Option<Self>;

    /// Returns true if the x, y, and z components are finite
    // Annotation to disable bad clippy lint; It is not good to use
    // `&self` as clippy suggests since Coord is Copy
    #[allow(clippy::wrong_self_convention)]
    fn is_finite(self) -> bool;
}

impl<T> Vector3DOps for Coord<T>
where
    T: CoordNum,
{
    type Scalar = T;

    fn wedge_product(self, other: Coord<T>) -> Self::Scalar {
        self.x * other.y - self.y * other.x
    }

    fn magnitude(self) -> Self::Scalar {
        self.dot(self).sqrt()
    }

    fn magnitude_squared(self) -> Self::Scalar {
        self.dot(self)
    }

    /// This assumes a 2D rotation in the XY plane, leaving Z unchanged.
    fn left(self) -> Self {
        Self {
            x: -self.y,
            y: self.x,
            z: self.z,
        }
    }

    /// This assumes a 2D rotation in the XY plane, leaving Z unchanged.
    fn right(self) -> Self {
        Self {
            x: self.y,
            y: -self.x,
            z: self.z,
        }
    }

    fn try_normalize(self) -> Option<Self> {
        let magnitude = self.magnitude();
        let result = self / magnitude;
        // Both the result AND the magnitude must be finite they are finite
        // Otherwise very large vectors overflow magnitude to Infinity,
        // and the after the division the result would be coord!{x:0.0, y:0.0, z:0.0}
        // Note we don't need to check if magnitude is zero, because after the division
        // that would have made result non-finite or NaN anyway.
        if result.is_finite() && magnitude.is_finite() {
            Some(result)
        } else {
            None
        }
    }

    fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }
}

#[cfg(test)]
mod test {
    use super::Vector3DOps;
    use crate::coord;

    #[test]
    fn test_wedge_product() {
        // perpendicular unit length
        let a = coord! { x: 1f64, y: 0f64, z: 0.0 };
        let b = coord! { x: 0f64, y: 1f64, z: 0.0 };

        // expect the area of the parallelogram
        assert_eq!(a.wedge_product(b), 1f64);
        // swapping the operands reverses the sign
        assert_eq!(b.wedge_product(a), -1f64);

        // Add skew; the results should be the same
        let a = coord! { x: 1f64, y: 0f64, z: 0.0 };
        let b = coord! { x: 1f64, y: 1f64, z: 0.0 };

        assert_eq!(a.wedge_product(b), 1f64);
        assert_eq!(b.wedge_product(a), -1f64);

        // Collinear case should yield zero
        let a = coord! { x: 2f64, y: 2f64, z: 3f64 };
        let b = coord! { x: 1f64, y: 1f64, z: 1f64 };
        assert_eq!(a.wedge_product(b), 0f64);
    }

    #[test]
    fn test_dot_product() {
        // Perpendicular unit vectors in the XY plane.
        let a = coord! { x: 1f64, y: 0f64, z: 0.0 };
        let b = coord! { x: 0f64, y: 1f64, z: 0.0 };
        assert_eq!(a.dot(b), 0f64);

        // Parallel vectors in the same direction.
        let a = coord! { x: 1f64, y: 0f64, z: 0.0 };
        let b = coord! { x: 2f64, y: 0f64, z: 0.0 };
        assert_eq!(a.dot(b), 2f64);
        assert_eq!(b.dot(a), 2f64);

        // Parallel but opposite vectors.
        let a = coord! { x: 3f64, y: 4f64, z: 0.0 };
        let b = coord! { x: -3f64, y: -4f64, z: 0.0 };
        assert_eq!(a.dot(b), -25f64);
        assert_eq!(b.dot(a), -25f64);
    }

    #[test]
    fn test_magnitude() {
        let a = coord! { x: 1.0, y: 0.0, z: 0.0 };
        assert_eq!(a.magnitude(), 1f64);

        let a = coord! { x: 0.0, y: 0.0, z: 1.0 };
        assert_eq!(a.magnitude(), 1f64);

        let a = coord! { x: 0.0, y: 0.0, z: 0.0 };
        assert_eq!(a.magnitude(), 0f64);

        let a = coord! { x: -3.0, y: 4.10008, z: 2.10708 };
        assert_relative_eq!(a.magnitude(), 5.5, epsilon = 0.00005);
    }

    #[test]
    fn test_magnitude_squared() {
        let a = coord! { x: 1f64, y: 0f64, z: 0.0 };
        assert_eq!(a.magnitude_squared(), 1f64);

        let a = coord! { x: 0f64, y: 0f64, z: 0.0 };
        assert_eq!(a.magnitude_squared(), 0f64);

        let a = coord! { x: -3f64, y: 4f64, z: 0.0 };
        assert_eq!(a.magnitude_squared(), 25f64);
    }

    #[test]
    fn test_left_right() {
        // Use a coordinate with a nonzero z so that we can confirm that left/right
        // rotations leave z unchanged.
        let a = coord! { x: 1f64, y: 0f64, z: 2.0 };
        let a_left = coord! { x: -0f64, y: 1f64, z: 2.0 };
        let a_right = coord! { x: 0f64, y: -1f64, z: 2.0 };

        assert_eq!(a.left(), a_left);
        assert_eq!(a.right(), a_right);
        assert_eq!(a.left(), -a.right());
    }

    #[test]
    fn test_left_right_match_rotate() {
        use crate::algorithm::rotate::Rotate;
        use crate::Point;
        // The aim of this test is to confirm that wording in documentation is
        // consistent.

        // when the user is in a coordinate system where the y axis is flipped
        // (eg screen coordinates in a HTML canvas), then rotation directions
        // will be different to those described in the documentation.

        // The documentation for the Rotate trait says: 'Positive angles are
        // counter-clockwise, and negative angles are clockwise rotations'

        let counter_clockwise_rotation_degrees = 90.0;
        let clockwise_rotation_degrees = -counter_clockwise_rotation_degrees;

        let a: Point = coord! { x: 1.0, y: 0.0, z: 0.0 }.into();
        let origin: Point = coord! { x: 0.0, y: 0.0, z: 0.0 }.into();

        // left() is equivalent to a 90° counter-clockwise rotation.
        assert_relative_eq!(
            Point::from(a.0.left()),
            a.rotate_around_point(counter_clockwise_rotation_degrees, origin),
        );
        // right() is equivalent to a 90° clockwise rotation.
        assert_relative_eq!(
            Point::from(a.0.right()),
            a.rotate_around_point(clockwise_rotation_degrees, origin),
        );
    }

    #[test]
    fn test_try_normalize() {
        // Already normalized.
        let a = coord! {
            x: 1.0,
            y: 0.0,
            z: 0.0
        };
        assert_relative_eq!(a.try_normalize().unwrap(), a);

        // Already normalized.
        let a = coord! {
            x: 1.0 / f64::sqrt(2.0),
            y: -1.0 / f64::sqrt(2.0),
            z: 1.0 / f64::sqrt(2.0),
        };
        assert_relative_eq!(a.try_normalize().unwrap(), a);

        // A nontrivial example.
        let a = coord! { x: -10.0, y: 8.0, z: 0.0 };
        assert_relative_eq!(
            a.try_normalize().unwrap(),
            coord! { x: -10.0, y: 8.0, z: 0.0 } / f64::sqrt(10.0 * 10.0 + 8.0 * 8.0)
        );
    }

    #[test]
    fn test_try_normalize_edge_cases_1() {
        use float_next_after::NextAfter;
        // Very small input: normalization should still succeed.
        let a = coord! {
            x: 0.0,
            y: 1e-301_f64,
            z: 0.0
        };
        assert_eq!(
            a.try_normalize(),
            Some(coord! {
                x: 0.0,
                y: 1.0,
                z: 0.0,
            })
        );

        // A large vector that does not overflow.
        let a = coord! {
            x: f64::sqrt(f64::MAX/2.0),
            y: f64::sqrt(f64::MAX/2.0),
            z: f64::sqrt(f64::MAX/2.0),
        };
        assert_relative_eq!(
            a.try_normalize().unwrap(),
            coord! {
                x: 1.0 / f64::sqrt(2.0),
                y: 1.0 / f64::sqrt(2.0),
                z: 1.0 / f64::sqrt(2.0),
            }
        );

        // A large vector that is barely under the overflow threshold.
        let a = coord! {
            x: f64::sqrt(f64::MAX / 2.0),
            y: f64::sqrt(f64::MAX / 2.0).next_after(f64::INFINITY),
            z: f64::sqrt(f64::MAX / 2.0)
        };
        assert_relative_eq!(
            a.try_normalize().unwrap(),
            coord! {
                x: 1.0 / f64::sqrt(2.0),
                y: 1.0 / f64::sqrt(2.0),
                z: 1.0 / f64::sqrt(2.0),
            }
        );
    }

    #[test]
    fn test_try_normalize_edge_cases_2() {
        // Zero vector: normalization should fail.
        let a = coord! { x: 0.0, y: 0.0, z: 0.0 };
        assert_eq!(a.try_normalize(), None);

        // Component is NaN: normalization should fail.
        let a = coord! { x: f64::NAN, y: 0.0, z: 0.0 };
        assert_eq!(a.try_normalize(), None);

        // Component is Infinite: normalization should fail.
        let a = coord! { x: f64::INFINITY, y: 0.0, z: 0.0 };
        assert_eq!(a.try_normalize(), None);
    }
}
