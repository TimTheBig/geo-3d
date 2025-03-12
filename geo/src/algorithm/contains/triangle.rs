use super::{impl_contains_from_relate, impl_contains_geometry_for, Contains};
use crate::{geometry::*, GeoNum, Vector3DOps};

// ┌──────────────────────────────┐
// │ Implementations for Triangle │
// └──────────────────────────────┘

// from maths-rs
/// Returns the [barycentric](https://en.wikipedia.org/wiki/Barycentric_coordinate_system) coordinate `(u, v, w)` of coord inside triangle
pub(crate) fn barycentric<T: GeoNum>(p: Coord<T>, tri: &Triangle<T>) -> (T, T, T) {
    let x13 = tri.0 - tri.2;
    let x23 = tri.1 - tri.2;
    let x03 = p - tri.2;
    let m13 = x13.magnitude_squared();
    let m23 = x23.magnitude_squared();

    let d = x13.dot(x23);
    let invdet = T::one() / T::max(m13 * m23 - d * d, T::epsilon());
    let a = x13.dot(x03);
    let b = x23.dot(x03);

    let u = invdet * (m23 * a - d * b);
    let v = invdet * (m13 * b - d * a);
    let w = T::one() - u - v;
    (u, v, w)
}

impl<T: GeoNum> Contains<Coord<T>> for Triangle<T> {
    fn contains(&self, coord: &Coord<T>) -> bool {
        // from [maths-rs](https://github.com/polymonster/maths-rs)
        // True if the coord is inside the triangle
        let (u, v, w) = barycentric(*coord, self);
        u > T::zero() && v > T::zero() && w > T::zero()
    }
}

impl<T: GeoNum> Contains<Point<T>> for Triangle<T> {
    fn contains(&self, point: &Point<T>) -> bool {
        self.contains(&point.0)
    }
}

impl_contains_from_relate!(Triangle<T>, [Line<T>, LineString<T>, Polygon<T>, MultiPoint<T>, MultiLineString<T>, MultiPolygon<T>, GeometryCollection<T>, Rect<T>, Triangle<T>]);
impl_contains_geometry_for!(Triangle<T>);
