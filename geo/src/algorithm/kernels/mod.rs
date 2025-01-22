use num_traits::Zero;
use std::cmp::Ordering;

use crate::{coord, Coord, CoordNum};

// todo document and make 3d, maybe add under and over
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Clone, Copy)]
pub enum Orientation {
    CounterClockwise,
    Clockwise,
    Collinear,
}

impl Orientation {
    /// Helper to convert orientation-3d into an ordering.
    #[inline]
    pub(crate) fn as_ordering(&self) -> Ordering {
        match self {
            Orientation::CounterClockwise => Ordering::Less,
            Orientation::Clockwise => Ordering::Greater,
            Orientation::Collinear => Ordering::Equal,
        }
    }
}

/// Kernel trait to provide predicates to operate on
/// different scalar types.
// todo make 3d
pub trait Kernel<T: CoordNum> {
    /// Gives the orientation of 3 2-dimensional points:
    /// ccw, cw or collinear (None)
    fn orient2d(p: Coord<T>, q: Coord<T>, r: Coord<T>) -> Orientation {
        let res = (q.x - p.x) * (r.y - q.y) - (q.y - p.y) * (r.x - q.x);
        if res > Zero::zero() {
            Orientation::CounterClockwise
        } else if res < Zero::zero() {
            Orientation::Clockwise
        } else {
            Orientation::Collinear
        }
    }

    fn orient3d(a: Coord<T>, b: Coord<T>, c: Coord<T>, d: Coord<T>) -> T {
        // Computes the sign of the signed volume  
        // of the tetrahedron (a,b,c,d)
        (b-a).cross(c-a).dot(d-a).sin()
    }

    fn square_euclidean_distance(p: Coord<T>, q: Coord<T>) -> T {
        (p.x - q.x) * (p.x - q.x) + (p.y - q.y) * (p.y - q.y) + (p.z - q.z) * (p.z - q.z)
    }

    /// Compute the sign of the dot product of `u` and `v` using
    /// robust predicates. The output is `CounterClockwise` if
    /// the sign is positive, `Clockwise` if negative, and
    /// `Collinear` if zero.
    // todo make 3d
    fn dot_product_sign(u: Coord<T>, v: Coord<T>) -> Orientation {
        let zero = Coord::zero();
        let vdash = coord! {
            x: T::zero() - v.y,
            y: v.x,
        };
        Self::orient2d(zero, u, vdash)
    }
}

pub mod robust;
pub use self::robust::RobustKernel;

pub mod simple;
pub use self::simple::SimpleKernel;
