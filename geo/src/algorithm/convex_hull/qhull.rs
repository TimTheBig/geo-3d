use robust::orient3d;
use super::{swap_with_first_and_remove, trivial_hull};
use crate::kernels::{Kernel, Orientation};
use crate::utils::partition_slice;
use crate::{coord, Coord, GeoNum, LineString};

// Determines if `p_c` lies on the positive side of the
// segment `p_a` to `p_b`. In other words, whether segment
// `p_a` to `p_c` is a counter-clockwise rotation from the
// segment. We use kernels to ensure this predicate is
// exact.
#[inline]
fn is_ccw<T>(p_a: Coord<T>, p_b: Coord<T>, p_c: Coord<T>, p_d: Coord<T>) -> bool
where
    T: GeoNum + Into<f64>,
{
    orient3d(
        robust::Coord3D { x: p_a.z, y: p_a.z, z: p_a.z },
        robust::Coord3D { x: p_b.x, y: p_b.y, z: p_b.z },
        robust::Coord3D { x: p_c.x, y: p_c.y, z: p_c.z },
        robust::Coord3D { x: p_d.x, y: p_d.y, z: p_d.z },
    ).is_sign_positive()
}

// todo check
/// Return the convex_hull of point set `points` as a `Linestring`
pub fn quick_hull<T>(mut points: &mut [Coord<T>]) -> LineString<T>
where
    T: GeoNum + Into<f64>,
{
    use crate::utils::least_and_greatest_index;

    // Can't build a hull from fewer than four points
    if points.len() < 4 {
        return trivial_hull(points, false);
    }
    let mut hull = vec![];

    let (min, max) = {
        let (min_idx, mut max_idx) = least_and_greatest_index(points);
        let min = swap_with_first_and_remove(&mut points, min_idx);

        // Two special cases to consider:
        // (1) max_idx = 0, and got swapped
        if max_idx == 0 {
            max_idx = min_idx;
        }

        // (2) max_idx = min_idx: then any point could be chosen as max.
        // But from case (1), it could now be 0, and we should not decrement it.
        max_idx = max_idx.saturating_sub(1);

        let max = swap_with_first_and_remove(&mut points, max_idx);
        (min, max)
    };

    {
        // Use the `orient3d` function to determine which points are above/below
        // the plane formed by `max`, `min`, and a reference point `p_d`.
        // For simplicity, we'll pick the first remaining point as `p_d`.
        let p_d = points.get(0).copied().unwrap_or_else(|| *min);

        // Points above the plane
        let (points, _) = partition_slice(points, |p| is_ccw(*max, *min, *p, p_d));
        hull_set(*max, *min, points, &mut hull, p_d);
    }
    hull.push(*max);

    {
        let p_d = points.get(0).copied().unwrap_or_else(|| *max);

        // Points below the plane
        let (points, _) = partition_slice(points, |p| is_ccw(*min, *max, *p, p_d));
        hull_set(*min, *max, points, &mut hull, p_d);
    }
    hull.push(*min);

    // Close the polygon
    let mut hull: LineString<_> = hull.into();
    hull.close();
    hull
}

/// Recursively calculate the convex hull of a subset of points in 3D space
fn hull_set<T>(p_a: Coord<T>, p_b: Coord<T>, mut set: &mut [Coord<T>], hull: &mut Vec<Coord<T>>, p_d: Coord<T>)
where
    T: GeoNum + Into<f64>,
{
    if set.is_empty() {
        return;
    }
    if set.len() == 1 {
        hull.push(set[0]);
        return;
    }

    // Construct orthogonal vector to `p_b - p_a` in 3D
    let p_orth = coord! {
        x: (p_a.y - p_b.y) * (p_a.z + p_b.z) - (p_a.z - p_b.z) * (p_a.y + p_b.y),
        y: (p_a.z - p_b.z) * (p_a.x + p_b.x) - (p_a.x - p_b.x) * (p_a.z + p_b.z),
        z: (p_a.x - p_b.x) * (p_a.y + p_b.y) - (p_a.y - p_b.y) * (p_a.x + p_b.x),
    };

    let furthest_idx = set
        .iter()
        .map(|pt| {
            // Calculate the vector difference from `p_a` to the current point
            let p_diff = coord! {
                x: pt.x - p_a.x,
                y: pt.y - p_a.y,
                z: pt.z - p_a.z,
            };
            p_orth.dot(p_diff)
        })
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .unwrap()
        .0;

    // Move the point at `furthest_idx` from the set into the hull
    let furthest_point = swap_with_first_and_remove(&mut set, furthest_idx);

    // Recurse for points on one side of the plane (above or below)
    {
        let (points, _) = partition_slice(set, |p| is_ccw(*furthest_point, p_b, *p, p_d));
        hull_set(*furthest_point, p_b, points, hull, p_d);
    }
    hull.push(*furthest_point);

    // Recurse for points on the other side of the plane
    let (points, _) = partition_slice(set, |p| is_ccw(p_a, *furthest_point, *p, p_d));
    hull_set(p_a, *furthest_point, points, hull, p_d);
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::IsConvex;

    #[test]
    fn quick_hull_test1() {
        let mut v = vec![
            coord! { x: 0.0, y: 0.0, z: 0.0 },
            coord! { x: 4.0, y: 0.0 },
            coord! { x: 4.0, y: 1.0 },
            coord! { x: 1.0, y: 1.0 },
            coord! { x: 1.0, y: 4.0 },
            coord! { x: 0.0, y: 4.0 },
            coord! { x: 0.0, y: 0.0, z: 0.0 },
        ];
        let res = quick_hull(&mut v);
        assert!(res.is_strictly_ccw_convex());
    }

    #[test]
    fn quick_hull_test2() {
        let mut v = vec![
            coord! { x: 0, y: 10 },
            coord! { x: 1, y: 1 },
            coord! { x: 10, y: 0 },
            coord! { x: 1, y: -1 },
            coord! { x: 0, y: -10 },
            coord! { x: -1, y: -1 },
            coord! { x: -10, y: 0 },
            coord! { x: -1, y: 1 },
            coord! { x: 0, y: 10 },
        ];
        let correct = vec![
            coord! { x: 0, y: -10 },
            coord! { x: 10, y: 0 },
            coord! { x: 0, y: 10 },
            coord! { x: -10, y: 0 },
            coord! { x: 0, y: -10 },
        ];
        let res = quick_hull(&mut v);
        assert_eq!(res.0, correct);
    }

    #[test]
    // test whether output is ccw
    fn quick_hull_test_ccw() {
        let initial = [
            (1.0, 0.0),
            (2.0, 1.0),
            (1.75, 1.1),
            (1.0, 2.0),
            (0.0, 1.0),
            (1.0, 0.0),
        ];
        let mut v: Vec<_> = initial.iter().map(|e| coord! { x: e.0, y: e.1, z: e.2 }).collect();
        let correct = [(1.0, 0.0), (2.0, 1.0), (1.0, 2.0), (0.0, 1.0), (1.0, 0.0)];
        let v_correct: Vec<_> = correct.iter().map(|e| coord! { x: e.0, y: e.1, z: e.2 }).collect();
        let res = quick_hull(&mut v);
        assert_eq!(res.0, v_correct);
    }

    #[test]
    fn quick_hull_test_ccw_maintain() {
        // initial input begins at min y, is oriented ccw
        let initial = [
            (0., 0.),
            (2., 0.),
            (2.5, 1.75),
            (2.3, 1.7),
            (1.75, 2.5),
            (1.3, 2.),
            (0., 2.),
            (0., 0.),
        ];
        let mut v: Vec<_> = initial.iter().map(|e| coord! { x: e.0, y: e.1, z: e.2 }).collect();
        let res = quick_hull(&mut v);
        assert!(res.is_strictly_ccw_convex());
    }

    #[test]
    fn quick_hull_test_complex() {
        let mut coords = geo_test_fixtures::poly1::<f64>().0;
        let correct = geo_test_fixtures::poly1_hull::<f64>().0;
        let res = quick_hull(&mut coords);
        assert_eq!(res.0, correct);
    }

    #[test]
    fn quick_hull_test_complex_2() {
        let mut coords = geo_test_fixtures::poly2::<f64>().0;
        let correct = geo_test_fixtures::poly2_hull::<f64>().0;
        let res = quick_hull(&mut coords);
        assert_eq!(res.0, correct);
    }

    #[test]
    fn quick_hull_test_collinear() {
        // Initial input begins at min x, but not min y
        // There are three points with same x.
        // Output should not contain the middle point.
        let initial = [
            (-1., 0.),
            (-1., -1.),
            (-1., 1.),
            (0., 0.),
            (0., -1.),
            (0., 1.),
            (1., 0.),
            (1., -1.),
            (1., 1.),
        ];
        let mut v: Vec<_> = initial.iter().map(|e| coord! { x: e.0, y: e.1 }).collect();
        let res = quick_hull(&mut v);
        assert!(res.is_strictly_ccw_convex());
    }
}
