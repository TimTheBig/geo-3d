//! This module defines the [Vector3DOps] trait and implements it for the
//! [Coord] struct.

use crate::{Coord, Point, CoordNum};

/// Defines vector operations for 3D coordinate types which implement `CoordNum`
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

    /// The squared distance between this coordinate and the origin.\
    /// (Avoids the square root calculation when it is not needed)
    ///
    /// `x² + y² + z²`
    ///
    fn magnitude_squared(self) -> Self::Scalar;

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

impl<T: CoordNum> Vector3DOps for Coord<T> {
    type Scalar = T;

    fn magnitude(self) -> Self::Scalar {
        self.dot(self).sqrt()
    }

    fn magnitude_squared(self) -> Self::Scalar {
        self.dot(self)
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

impl<T: CoordNum> Vector3DOps for Point<T> {
    type Scalar = T;

    fn magnitude(self) -> Self::Scalar {
        self.0.magnitude()
    }

    fn magnitude_squared(self) -> Self::Scalar {
        self.0.magnitude_squared()
    }

    fn try_normalize(self) -> Option<Self> {
        self.0.try_normalize().map(|c| Point(c))
    }

    fn is_finite(self) -> bool {
        self.0.is_finite()
    }
}

#[cfg(test)]
mod test {
    use super::Vector3DOps;
    use crate::coord;

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
