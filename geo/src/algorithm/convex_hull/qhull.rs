use geo_types::{coord, Coord, CoordNum, LineString};
use glam::DVec3;
pub use quickhull::{ErrorKind, ConvexHull};
use crate::GeoNum;
use super::trivial_hull;

/// A container to make using [quickhull](https://github.com/TimTheBig/quickhull) with geo simpler
/// 
/// # Example
/// ```
/// # use geo_3d::{coord, Coord, convex_hull::ConvexQHull};
/// # use geo_3d::{polygon, Centroid};
/// let qhull = ConvexQHull::try_new(
///     &[
///     coord!(5.0, 5.0, 6.0),
///     coord!(2.0, 3.0, 4.0),
///     coord!(5.0, 7.0, 4.0),
///     // this one is inside
///     coord!(4.3, 4.1, 5.0),
///     coord!(3.0, 4.0, 3.0),
///     coord!(6.0, 4.5, 5.0),
///     ]).expect("this is valid");
///
///     assert_eq!(
///         qhull.points::<f64>(),
///         vec![coord!(5.0, 5.0, 6.0), coord!(6.0, 5.0, 5.0), coord!(2.0, 3.0, 4.0)]
///     )
/// ```
pub struct ConvexQHull(ConvexHull);

impl ConvexQHull {
    /// Attempts to compute a [`ConvexQHull`] for the given set of points.
    pub fn try_new<T: CoordNum>(points: &[Coord<T>]) -> Result<Self, ErrorKind> {
        match ConvexHull::try_new(
            &(points.iter().map(glam_vec3_from_geo)).collect::<Vec<_>>(),
            None
        ) {
            Ok(ch) => Ok(ConvexQHull(ch)),
            Err(e) => Err(e),
        }
    }

    /// Gets the points of the convex hull.
    pub fn points<T: CoordNum>(&self) -> Vec<Coord<f64>> {
        self.0.points.iter().map(|v3| {
            coord!{
                x: v3.x,
                y: v3.y,
                z: v3.z,
            }
        }).collect::<Vec<_>>()
    }

    /// Gets the points of the convex hull without collecting an iterator.
    pub fn points_iter<'a, T: CoordNum>(&self) -> impl Iterator<Item = Coord<f64>> + '_ {
        self.0.points.iter().map(|v3| coord!(v3.x, v3.y, v3.z))
    }

    /// Adds the given points to the point set, attempting to update the convex hull.
    pub fn add_points<T: CoordNum + Into<f64>>(&mut self, points: &[Coord<T>]) -> Result<(), ErrorKind> {
        self.0.add_iter_points(&mut points.iter().map(glam_vec3_from_geo))
    }
}

/// Convert `Coord` to `glam::DVec3`
fn glam_vec3_from_geo<T: CoordNum>(coord: &Coord<T>) -> DVec3 {
    DVec3::new(
        coord.x.to_f64().expect("CoordNum is float"),
        coord.y.to_f64().expect("CoordNum is float"),
        coord.z.to_f64().expect("CoordNum is float")
    )
}

/// Return the convex_hull of point set `points` as a `Linestring`
pub fn quick_hull<T: CoordNum + GeoNum + From<f64>>(points: &mut [Coord<T>]) -> Result<LineString<T>, ErrorKind> {
    // Can't build a hull from fewer than four points
    if points.len() < 4 {
        return Ok(trivial_hull(points, false));
    }

    match ConvexQHull::try_new(points) {
        Ok(ps) => Ok(LineString(ps.points_iter::<T>().map(|c| { coord!(c.x.into(), c.y.into(), c.z.into()) }).collect())),
        Err(e) => Err(e),
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::IsConvex;

    #[test]
    fn sphere_test() {
        let (_v, _i) = ConvexQHull::try_new(&geo_test_fixtures::sphere::<f64>().0)
            .unwrap()
            .0
            .vertices_indices();
    }

    #[test]
    fn quick_hull_test1() {
        let mut v = vec![
            coord! { x: 0.0, y: 0.0, z: 0.0 },
            coord! { x: 4.0, y: 0.0, z: 4.0 },
            coord! { x: 4.0, y: 1.0, z: 4.0 },
            coord! { x: 1.0, y: 1.0, z: 1.0 },
            coord! { x: 1.0, y: 4.0, z: 1.0 },
            coord! { x: 0.0, y: 4.0, z: 0.0 },
            coord! { x: 0.0, y: 0.0, z: 0.0 },
        ];
        let res = quick_hull(&mut v);
        assert!(res.unwrap().is_strictly_ccw_convex());
    }

    #[test]
    fn quick_hull_test2() {
        let mut v = vec![
            coord! { x: 0., y: 10., z: 0. },
            coord! { x: 1., y: 1., z: 1. },
            coord! { x: 10., y: 0., z: 10. },
            coord! { x: 1., y: -1., z: 1. },
            coord! { x: 0., y: -10., z: 0. },
            coord! { x: -1., y: -1., z: -1. },
            coord! { x: -10., y: 0., z: -10. },
            coord! { x: -1., y: 1., z: -1. },
            coord! { x: 0., y: 10., z: 0. },
        ];
        let correct = vec![
            coord! { x: 0., y: -10., z: 0. },
            coord! { x: 10., y: 0., z: 10. },
            coord! { x: 0., y: 10., z: 0. },
            coord! { x: -10., y: 0., z: -10. },
            coord! { x: 0., y: -10., z: 0. },
        ];
        let res = quick_hull(&mut v);
        assert_eq!(res.unwrap().0, correct);
    }

    #[test]
    // test whether output is ccw
    fn quick_hull_test_ccw() {
        let mut initial = [
            coord!(1.0, 0.0, 1.0),
            coord!(2.0, 1.0, 2.0),
            coord!(1.75, 1.1, 1.2),
            coord!(1.0, 2.0, 1.0),
            coord!(0.0, 1.0, 0.0),
            coord!(1.0, 0.0, 1.0),
        ].to_vec();
        let correct = [
            coord!(1.0, 0.0, 1.0),
            coord!(2.0, 1.0, 2.0),
            coord!(1.0, 2.0, 1.0),
            coord!(0.0, 1.0, 0.0),
            coord!(1.0, 0.0, 1.0)
        ].to_vec();
        let res = quick_hull(&mut initial);
        assert_eq!(res.unwrap().0, correct);
    }

    #[test]
    fn quick_hull_test_ccw_maintain() {
        // initial input begins at min y, is oriented ccw
        let mut initial = [
            coord!(0., 0., 0.),
            coord!(2., 0., 2.),
            coord!(2.5, 1.75, 2.5),
            coord!(2.3, 1.7, 2.3),
            coord!(1.75, 2.5, 1.75),
            coord!(1.3, 2., 1.3),
            coord!(0., 2., 0.),
            coord!(0., 0., 0.),
        ].to_vec();
        let res = quick_hull(&mut initial);
        assert!(res.unwrap().is_strictly_ccw_convex());
    }

    #[test]
    fn quick_hull_test_complex() {
        let mut coords = geo_test_fixtures::poly1::<f64>().0;
        let correct = geo_test_fixtures::poly1_hull::<f64>().0;
        let res = quick_hull(&mut coords);
        assert_eq!(res.unwrap().0, correct);
    }

    #[test]
    fn quick_hull_test_complex_2() {
        let mut coords = geo_test_fixtures::poly2::<f64>().0;
        let correct = geo_test_fixtures::poly2_hull::<f64>().0;
        let res = quick_hull(&mut coords);
        assert_eq!(res.unwrap().0, correct);
    }

    #[test]
    fn quick_hull_test_collinear() {
        // Initial input begins at min x, but not min y
        // There are three points with same x.
        // Output should not contain the middle point.
        let mut initial = [
            coord!(-1., 0., -1.),
            coord!(-1., -1., -1.),
            coord!(-1., 1., -1.),
            coord!(0., 0., 0.),
            coord!(0., -1., 0.),
            coord!(0., 1., 0.),
            coord!(1., 0., 1.),
            coord!(1., -1., 1.),
            coord!(1., 1., 1.),
        ].to_vec();
        let res = quick_hull(&mut initial);
        assert!(res.unwrap().is_strictly_ccw_convex());
    }
}
