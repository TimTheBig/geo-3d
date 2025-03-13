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

/// This does not convert plane to a unit vec
#[inline(always)]
fn proj_coord<T: CoordNum>(c: Coord<T>, plane: Coord<T>) -> Coord<T> {
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

        self.exterior_mut(|ls| {
            for i in 0..ls.0.len() {
                ls[i] = proj_coord(ls[i], plane)
            }
        });

        self.interiors_mut(|ls_s| {
            for s_i in 0..ls_s.len() {
                for i in 0..ls_s[s_i].0.len() {
                    ls_s[s_i][i] = proj_coord(ls_s[s_i][i], plane)
                }
            }
        });
    }
}

/// Gets an iterator of the right most `Coord`s
pub trait RightMostIter<T: CoordNum> {
    fn right_most_iter(&self) -> impl Iterator<Item = Coord<T>>;
}

impl<T: CoordNum> RightMostIter<T> for LineString<T> {
    fn right_most_iter(&self) -> impl Iterator<Item = Coord<T>> {
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

/// The `Orientation` of **q** relitive to **p** and **r**.
/// 2d Left/Right
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
    use geo_types::coord;
    use super::{orientation, Orientation};

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
}
