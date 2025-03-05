use crate::geometry::{Coord, LineString, Polygon};
use crate::kernels::*;
use crate::winding_order::triangle_winding_order;
use crate::GeoNum;

/// Returns the convex hull of a geometry. The hull is always oriented counter-clockwise.
///
/// This implementation uses the QuickHull algorithm,
/// based on [Barber, C. Bradford; Dobkin, David P.; Huhdanpaa, Hannu (1 December 1996)](https://dx.doi.org/10.1145%2F235815.235821)
/// Original paper here: <http://www.cs.princeton.edu/~dpd/Papers/BarberDobkinHuhdanpaa.pdf>
///
/// # Examples
///
/// ```
/// use geo::{line_string, polygon};
/// use geo::ConvexHull;
///
/// // an L shape
/// let poly = polygon![
///     (x: 0.0, y: 0.0),
///     (x: 4.0, y: 0.0),
///     (x: 4.0, y: 1.0),
///     (x: 1.0, y: 1.0),
///     (x: 1.0, y: 4.0),
///     (x: 0.0, y: 4.0),
///     (x: 0.0, y: 0.0),
/// ];
///
/// // The correct convex hull coordinates
/// let correct_hull = line_string![
///     (x: 4.0, y: 0.0),
///     (x: 4.0, y: 1.0),
///     (x: 1.0, y: 4.0),
///     (x: 0.0, y: 4.0),
///     (x: 0.0, y: 0.0),
///     (x: 4.0, y: 0.0),
/// ];
///
/// let res = poly.convex_hull();
/// assert_eq!(res.exterior(), &correct_hull);
/// assert_eq!(res.interiors(), &[]);
/// ```
pub trait ConvexHull<'a, T> {
    type Scalar: GeoNum;
    fn convex_hull(&'a self) -> Result<Polygon<Self::Scalar>, ErrorKind>;
}

use crate::algorithm::CoordsIter;
use crate::utils::lex_cmp;

impl<'a, T, G> ConvexHull<'a, T> for G
where
    T: CoordNum + GeoNum + From<f64>,
    G: CoordsIter<Scalar = T>,
{
    type Scalar = T;

    fn convex_hull(&'a self) -> Result<Polygon<T>, ErrorKind> {
        let mut exterior: Vec<_> = self.exterior_coords_iter().collect();
        Ok(Polygon::new(quick_hull(&mut exterior)?, vec![]))
    }
}

pub mod qhull;
use geo_types::{CoordNum, Triangle};
pub use qhull::{quick_hull, ConvexQHull};

pub mod graham;
pub use graham::graham_hull;
use quickhull::ErrorKind;

/// Helper function that outputs the convex hull in the
/// trivial case: input with at most 3 points.\
/// It ensures the output is ccw, and does not repeat points unless required.
fn trivial_hull<T: GeoNum>(points: &mut [Coord<T>], include_on_hull: bool) -> LineString<T> {
    assert!(points.len() < 4);

    let mut ls = points.to_vec();
    // Remove repeated coplanear points
    if !include_on_hull {
        ls.sort_unstable_by(lex_cmp);

        if ls.len() == 3 {
            let mut remove = None;
            for (i, d) in ls.iter().enumerate() {
                if T::Ker::orient3d(ls[0], ls[1], ls[2], *d) == T::zero() {
                    remove = Some(i)
                }
            }
            if let Some(i) = remove {
                ls.remove(i);
            }
        }
    }

    let mut ls: LineString<T> = ls.into();
    ls.close();

    // Maintain the CCW invariance
    use super::winding_order::WindingOrder;
    if ls.0.len() == 3 {
        if let Some(w_o) = triangle_winding_order(&Triangle::new(ls[0], ls[1], ls[2])) {
            if w_o == WindingOrder::Clockwise {
                ls.0.reverse();
            }
        }
    }

    ls
}

/// Utility function for convex hull ops
///
/// 1. _swap_ the element at `idx` with the element at `head` (0th position)
/// 2. remove the _new_ `head` element (modifying the slice)
/// 3. return a _mutable ref_ to the removed head element
fn swap_with_first_and_remove<'a, T>(slice: &mut &'a mut [T], idx: usize) -> &'a mut T {
    // temporarily replace `slice` with an empty value
    let tmp = std::mem::take(slice);
    tmp.swap(0, idx);
    let (h, t) = tmp.split_first_mut().unwrap();
    *slice = t;
    h
}

#[cfg(test)]
mod test;
