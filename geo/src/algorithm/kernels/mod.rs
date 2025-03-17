use ::robust::Coord3D;
use num_traits::NumCast;
use std::cmp::Ordering;

use crate::{Coord, CoordNum};

/// A 2d orientation of points
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Clone, Copy)]
pub enum Orientation {
    CounterClockwise,
    Clockwise,
    Collinear,
}

impl Orientation {
    /// Helper to convert orientation-3d into an ordering.
    #[inline]
    pub(crate) const fn as_ordering(&self) -> Ordering {
        match self {
            Orientation::CounterClockwise => Ordering::Less,
            Orientation::Clockwise => Ordering::Greater,
            Orientation::Collinear => Ordering::Equal,
        }
    }
}

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Clone, Copy)]
pub enum Orientation3D {
    /// `CounterClockwise`,
    /// a point lies `below` a plane if the plane appears in clockwise order when viewed from below
    Below,
    /// `Clockwise`,
    /// a point lies `above` a plane if the plane appears in counterclockwise order when viewed from above
    Above,
    CoPlanar,
}

/// Kernel trait to provide predicates to operate on
/// different scalar types.
pub trait Kernel<T: CoordNum> {
    /// Gives the orientation of 3 2-dimensional points:
    /// ccw, cw or collinear (None)
    fn orient2d(p: Coord<T>, q: Coord<T>, r: Coord<T>) -> Orientation {
        let res = (q.x - p.x) * (r.y - q.y) - (q.y - p.y) * (r.x - q.x);
        if res > T::zero() {
            Orientation::CounterClockwise
        } else if res < T::zero() {
            Orientation::Clockwise
        } else {
            Orientation::Collinear
        }
    }

    /// Returns a positive value if the point `d` lies below the plane passing through `a`, `b`, and `c`
    /// ("below" is defined so that `a`, `b`, and `c` appear in counterclockwise order when viewed from above the plane).  
    /// Returns a negative value if `d` lies above the plane.  
    /// Returns `0` if they are **coplanar**.
    ///
    /// # Example
    /// ```
    /// # use geo_3d::{coord, Kernel};
    /// # use geo_3d::kernels::RobustKernel;
    /// // plane
    /// let pa = coord!(1., 0., 1.);
    /// let pb = coord!(-1., 0., -1.);
    /// let pc = coord!(-1., 0., 0.);
    ///
    /// // above plane - negative value expected
    /// let p1 = coord!(::core::f64::MIN_POSITIVE, ::core::f64::MIN_POSITIVE, ::core::f64::MIN_POSITIVE);
    /// // below plane - positive value expected
    /// let p2 = coord!(-::core::f64::MIN_POSITIVE, -::core::f64::MIN_POSITIVE, -::core::f64::MIN_POSITIVE);
    // // collinear to plane - zero expected
    /// let p3 = coord! { x: 0., y: 0., z: 0. };
    ///
    /// for &(p, sign) in &[(p1, -1.0), (p2, 1.0), (p3, 0.0)] {
    ///     let det = <RobustKernel as Kernel<f64>>::orient3d(pa, pb, pc, p);
    ///     assert!(det == sign || det.signum() == sign.signum());
    /// }
    /// ```
    fn orient3d(a: Coord<T>, b: Coord<T>, c: Coord<T>, d: Coord<T>) -> T {
        <T as NumCast>::from(::robust::orient3d(
            to_robust_coord3d(a),
            to_robust_coord3d(b),
            to_robust_coord3d(c),
            to_robust_coord3d(d),
        )).unwrap()
    }

    /// Get the 3D orientation of point `d` on a plane passing through `a`, `b`, and `c`
    fn orientation_3d(a: Coord<T>, b: Coord<T>, c: Coord<T>, d: Coord<T>) -> Orientation3D {
        match <RobustKernel as Kernel<T>>::orient3d(a, b, c, d) {
            o if o < T::zero() => Orientation3D::Above,
            o if o > T::zero() => Orientation3D::Below,
            o if o == T::zero() => Orientation3D::CoPlanar,
            _ => unreachable!("All numbers are covered"),
        }
    }

    fn square_euclidean_distance(p: Coord<T>, q: Coord<T>) -> T {
        (p.x - q.x) * (p.x - q.x) + (p.y - q.y) * (p.y - q.y) + (p.z - q.z) * (p.z - q.z)
    }

    /// Compute the sign of the dot product of `u` and `v` using
    /// robust predicates. The output is `CounterClockwise` if
    /// the sign is positive, `Clockwise` if negative, and
    /// `Collinear` if zero
    fn dot_product_sign_2d(u: Coord<T>, v: Coord<T>) -> Orientation {
        let zero = Coord::zero();
        let vdash = Coord {
            x: T::zero() - v.y,
            y: v.x,
            z: T::zero(),
        };
        Self::orient2d(zero, u, vdash)
    }
}

pub mod robust;
pub use self::robust::RobustKernel;

#[inline]
fn to_robust_coord3d<T: CoordNum>(coord: Coord<T>) -> Coord3D<f64> {
    Coord3D {
        x: <f64 as NumCast>::from(coord.x).unwrap(),
        y: <f64 as NumCast>::from(coord.y).unwrap(),
        z: <f64 as NumCast>::from(coord.z).unwrap(),
    }
}
