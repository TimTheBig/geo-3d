//! A trait to project a coords to a 2d plane, with utils to see which is the right most.\
//! Impls for `Coord`, `Point`, `MultiPoint`, `LineString`, `Line`, `Triangle`, and `Polygon`

use geo_types::{coord, Coord, CoordNum, Line, LineString, MultiPoint, Point, Polygon, Triangle};
use itertools::Itertools;
use super::Vector3DOps;

/// Convert a direction to a unit vector, all fields are 0-1, for use in projection
fn to_unit_vec<T: CoordNum>(dir: Coord<T>) -> Coord<T> {
    // Direction vector of the ray
    let unit_direction = dir / dir.magnitude();
    unit_direction
}

/// Project a coord to a 2d plane
pub trait ProjectToPlane<T: CoordNum>: Sized {
    /// Project to a plane defined by a normal
    // FOR IMPLEMENTERS: this normal must first be converted to a unit vec
    fn proj(&self, plane: Coord<T>) -> Self;

    /// Project to a plane defined by a normal, updating `Self` in place.
    // FOR IMPLEMENTERS: this normal must first be converted to a unit vec
    fn proj_mut(&mut self, plane: Coord<T>);

    /// Project the geometry to the XY plane
    fn project_xy(&self) -> Self {
        self.proj(coord!(T::zero(), T::zero(), T::one()))
    }

    /// Project the geometry to the XZ plane
    fn project_xz(&self) -> Self {
        self.proj(coord!(T::zero(), T::one(), T::zero()))
    }

    /// Project the geometry to the YZ plane
    fn project_yz(&self) -> Self {
        self.proj(coord!(T::one(), T::zero(), T::zero()))
    }
}

/// Project a `Coord` to a plane returning the neerest point on the plane.
///
/// Note: This does not convert plane to a unit vec, the caller must do that
#[inline(always)]
fn proj_coord<T: CoordNum>(c: Coord<T>, plane: Coord<T>) -> Coord<T> {
    debug_assert!(
        plane.x >= T::zero() && plane.y >= T::zero() && plane.z >= T::zero()
        && plane.x <= T::one() && plane.y <= T::one() && plane.z <= T::one(),
        "The plane normal must be a unit vec, all fields are 0-1."
    );

    c - plane * (plane.dot(c))
}

impl<T: CoordNum> ProjectToPlane<T> for Coord<T> {
    fn proj(&self, plane: Coord<T>) -> Self {
        let plane = to_unit_vec(plane);

        proj_coord(*self, plane)
    }

    fn proj_mut(&mut self, plane: Coord<T>) {
        *self = self.proj(plane);
    }
}

impl<T: CoordNum> ProjectToPlane<T> for Point<T> {
    fn proj(&self, plane: Coord<T>) -> Self {
        Point(self.0.proj(plane))
    }

    fn proj_mut(&mut self, plane: Coord<T>) {
        let plane = to_unit_vec(plane);

        self.0 = proj_coord(self.0, plane);
    }
}

impl<T: CoordNum> ProjectToPlane<T> for Line<T> {
    fn proj(&self, plane: Coord<T>) -> Self {
        let plane = to_unit_vec(plane);

        Self {
            start: proj_coord(self.start, plane),
            end: proj_coord(self.end, plane),
        }
    }

    fn proj_mut(&mut self, plane: Coord<T>) {
        let plane = to_unit_vec(plane);

        self.start = proj_coord(self.start, plane);
        self.end = proj_coord(self.end, plane);
    }
}

impl<T: CoordNum> ProjectToPlane<T> for LineString<T> {
    /// Project to a plane defined by a normal.
    ///
    /// Note:
    /// - This allocates a new `Vec` to save an allocation use [`proj_mut`](ProjectToPlane::proj_mut)
    fn proj(&self, plane: Coord<T>) -> Self {
        let plane = to_unit_vec(plane);

        self.0.iter().map(|c| proj_coord(*c, plane)).collect()
    }

    fn proj_mut(&mut self, plane: Coord<T>) {
        let plane = to_unit_vec(plane);

        for i in 0..self.0.len() {
            self[i] = proj_coord(self[i], plane)
        }
    }

    fn project_xy(&self) -> Self {
        self.0.iter().map(|c| coord!(c.x, c.y, T::zero())).collect()
    }

    fn project_xz(&self) -> Self {
        self.0.iter().map(|c| coord!(c.x, T::zero(), c.z)).collect()
    }

    fn project_yz(&self) -> Self {
        self.0.iter().map(|c| coord!(T::zero(), c.y, c.z)).collect()
    }
}

impl<T: CoordNum> ProjectToPlane<T> for MultiPoint<T> {
    /// Project to a plane defined by a normal.
    ///
    /// Note:
    /// - This allocates a new `Vec` to save an allocation use [`proj_mut`](ProjectToPlane::proj_mut)
    fn proj(&self, plane: Coord<T>) -> Self {
        let plane = to_unit_vec(plane);

        self.0.iter().map(|c| proj_coord(c.0, plane)).collect()
    }

    fn proj_mut(&mut self, plane: Coord<T>) {
        let plane = to_unit_vec(plane);

        for i in 0..self.0.len() {
            self[i].0 = proj_coord(self[i].0, plane)
        }
    }

    fn project_xy(&self) -> Self {
        self.0.iter().map(|c| coord!(c.x(), c.y(), T::zero())).collect()
    }

    fn project_xz(&self) -> Self {
        self.0.iter().map(|c| coord!(c.x(), T::zero(), c.z())).collect()
    }

    fn project_yz(&self) -> Self {
        self.0.iter().map(|c| coord!(T::zero(), c.y(), c.z())).collect()
    }
}

impl<T: CoordNum> ProjectToPlane<T> for Triangle<T> {
    fn proj(&self, plane: Coord<T>) -> Self {
        let plane = to_unit_vec(plane);

        Triangle(
            proj_coord(self.0, plane),
            proj_coord(self.1, plane),
            proj_coord(self.2, plane),
        )
    }

    fn proj_mut(&mut self, plane: Coord<T>) {
        let plane = to_unit_vec(plane);

        self.0 = proj_coord(self.0, plane);
        self.1 = proj_coord(self.1, plane);
        self.2 = proj_coord(self.2, plane);
    }
}

impl<T: CoordNum> ProjectToPlane<T> for Polygon<T> {
    /// Project to a plane defined by a normal.
    ///
    /// Note:
    /// - This allocates a new `Vec` to save an allocation use [`proj_mut`](ProjectToPlane::proj_mut)
    fn proj(&self, plane: Coord<T>) -> Self {
        let mut new_poly = self.clone();
        new_poly.proj_mut(plane);

        new_poly
    }

    fn proj_mut(&mut self, plane: Coord<T>) {
        let plane = to_unit_vec(plane);

        set_poly_coord_fields(self, |c| proj_coord(c, plane));
    }

    fn project_xy(&self) -> Self {
        let mut new_poly = self.clone();
        set_poly_coord_fields(
            &mut new_poly,
            |c| coord!(c.x, c.y, T::zero())
        );

        new_poly
    }

    fn project_xz(&self) -> Self {
        let mut new_poly = self.clone();
        set_poly_coord_fields(
            &mut new_poly,
            |c| coord!(c.x, T::zero(), c.z)
        );

        new_poly
    }

    fn project_yz(&self) -> Self {
        let mut new_poly = self.clone();
        set_poly_coord_fields(
            &mut new_poly,
            |c| coord!(T::zero(), c.y, c.z)
        );

        new_poly
    }
}

fn set_poly_coord_fields<T: CoordNum>(poly: &mut Polygon<T>, set_fields: impl Fn(Coord<T>) -> Coord<T>) {
    poly.exterior_mut(|ls| {
        for i in 0..ls.0.len() {
            ls[i] = set_fields(ls[i])
        }
    });

    poly.interiors_mut(|ls_s| {
        for s_i in 0..ls_s.len() {
            for i in 0..ls_s[s_i].0.len() {
                ls_s[s_i][i] = set_fields(ls_s[s_i][i])
            }
        }
    });
}

/// Gets an iterator of the right most `Coord`s
pub trait RightMostIter<T: CoordNum> {
    fn right_most_iter(&self) -> impl Iterator<Item = Coord<T>>;
}

impl<T: CoordNum> RightMostIter<T> for LineString<T> {
    fn right_most_iter(&self) -> impl Iterator<Item = Coord<T>> {
        // todo use optimizations from https://github.com/JernejPuc/convex-hull
        self.0.iter()
            .take(self.0.len() - 1)
            .circular_tuple_windows::<(&Coord<T>, &Coord<T>, &Coord<T>)>()
            .filter_map(|(&next, &point, &prev)| {
                if orientation(&prev, &point, &next) == Orientation::Right {
                    Some(point)
                } else {
                    None
                }
            })
    }
}

impl<T: CoordNum> RightMostIter<T> for Polygon<T> {
    fn right_most_iter(&self) -> impl Iterator<Item = Coord<T>> {
        self.exterior().right_most_iter()
    }
}

// todo check
impl<T: CoordNum> RightMostIter<T> for Triangle<T> {
    fn right_most_iter(&self) -> impl Iterator<Item = Coord<T>> {
        match orientation(&self.0, &self.1, &self.2) {
            Orientation::Left => [self.2, self.1, self.0].into_iter(),
            Orientation::Right => [self.0, self.1, self.2].into_iter(),
            Orientation::Linear => unreachable!(),
        }
    }
}

#[derive(Debug, PartialEq)]
enum Orientation {
    Linear,
    Left,
    Right,
}

/// The 2d Left/Right `Orientation` of **q** relitive to **p** and **r**.
fn orientation<T: CoordNum>(p: &Coord<T>, q: &Coord<T>, r: &Coord<T>) -> Orientation {
    let left_val = (q.x - p.x) * (r.y - p.y);
    let right_val = (q.y - p.y) * (r.x - p.x);

    if left_val == right_val {
        Orientation::Linear
    } else if left_val > right_val {
        Orientation::Left
    } else {
        Orientation::Right
    }
}

#[cfg(test)]
mod test {
    use geo_types::line_string;
    use super::*;

    #[test]
    fn test_orientation() {
        assert_eq!(
            orientation(
                &coord!(-1.0, 0.0, 0.0),
                &coord!(0.0, 0.0, 0.0),
                &coord!(1.0, 0.0, 0.0),
            ),
            Orientation::Linear,
        );

        assert_eq!(
            orientation(
                &coord!(-1.0, 0.1, 0.0),
                &coord!(-1.1, 0.0, 0.0),
                &coord!(1.0, 0.0, 0.0),
            ),
            Orientation::Left,
        );

        assert_eq!(
            orientation(
                &coord!(-1.0, 0.1, 0.0),
                &coord!(1.1, 0.0, 0.0),
                &coord!(1.0, 0.0, 0.0),
            ),
            Orientation::Right,
        )
    }

    #[test]
    fn test_coord_proj() {
        let mut coord = coord!(1.5, 1.5, 2.0);
        coord.proj_mut(coord!(0.7071067812, 0.0, 0.7071067812));
        assert_eq!(
            // 45˚ plane
            coord!(1.5, 1.5, 2.0).proj(coord!(0.7071067812, 0.0, 0.7071067812)),
            coord,
        );

        assert_relative_eq!(
            // 45˚ plane
            coord!(1.5, 1.5, 2.0).proj(coord!(0.7071067812, 0.0, 0.7071067812)),
            coord!(-0.25, 1.5, 0.25),
            epsilon = 7e-16,
        );
        assert_eq!(
            coord!(1.5, 1.5, 2.0).project_xy(),
            coord!(1.5, 1.5, 0.0),
        );
        assert_eq!(
            coord!(1.5, 1.5, 2.0).project_xz(),
            coord!(1.5, 0.0, 2.0),
        );
        assert_eq!(
            coord!(1.5, 1.5, 2.0).project_yz(),
            coord!(0.0, 1.5, 2.0),
        );
    }

    #[test]
    fn test_line_proj() {
        let mut line = Line::new(coord!(1.5, 1.5, 2.0), coord!(9.0, 2.0, 5.5));
        line.proj_mut(coord!(0.7071067812, 0.0, 0.7071067812));
        assert_eq!(
            // 45˚ plane
            Line::new(coord!(1.5, 1.5, 2.0), coord!(9.0, 2.0, 5.5)).proj(coord!(0.7071067812, 0.0, 0.7071067812)),
            line,
        );

        assert_relative_eq!(
            // 45˚ plane
            Line::new(coord!(1.5, 1.5, 2.0), coord!(9.0, 2.0, 5.5)).proj(coord!(0.7071067812, 0.0, 0.7071067812)),
            Line::new(coord!(-0.25, 1.5, 0.25), coord!(1.75, 2.0, -1.75)),
            epsilon = 1e-14,
        );
        assert_eq!(
            Line::new(coord!(1.5, 1.5, 2.0), coord!(9.0, 2.0, 5.5)).project_xy(),
            Line::new(coord!(1.5, 1.5, 0.0), coord!(9.0, 2.0, 0.0)),
        );
        assert_eq!(
            Line::new(coord!(1.5, 1.5, 2.0), coord!(9.0, 2.0, 5.5)).project_xz(),
            Line::new(coord!(1.5, 0.0, 2.0), coord!(9.0, 0.0, 5.5)),
        );
        assert_eq!(
            Line::new(coord!(1.5, 1.5, 2.0), coord!(9.0, 2.0, 5.5)).project_yz(),
            Line::new(coord!(0.0, 1.5, 2.0), coord!(0.0, 2.0, 5.5)),
        );
    }

    #[test]
    fn test_linestring_proj() {
        let mut ls = line_string![
            coord!(0.0, 0.0, 0.0)
            coord!(
            coord!(
            coord!(
            coord!(
        ];
        ls.proj_mut(coord!(0.7071067812, 0.0, 0.7071067812));
        assert_eq!(
            // 45˚ plane
            .proj(coord!(0.7071067812, 0.0, 0.7071067812)),
            ls,
        );

        assert_relative_eq!(
            // 45˚ plane
            .proj(coord!(0.7071067812, 0.0, 0.7071067812)),
            epsilon = 1e-17,
        );
        assert_eq!(
            .project_xy(),
        );
        assert_eq!(
            .project_xz(),
        );
        assert_eq!(
            .project_yz(),
        );
    }
}
