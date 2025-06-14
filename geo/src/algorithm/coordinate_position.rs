use geo_types::{coord, CoordNum};
use std::cmp::Ordering;
use std::ops::AddAssign;

use super::Triangulate;
use crate::geometry::*;
use crate::intersects::{point_in_rect, value_in_between};
use crate::kernels::*;
use crate::{BoundingRect, HasDimensions, Intersects};
use crate::{GeoNum, GeometryCow};
use crate::contains::triangle::barycentric;

/// The position of a `Coord` relative to a `Geometry`
#[derive(PartialEq, Eq, Clone, Copy, Debug)]
pub enum CoordPos {
    OnBoundary,
    Inside,
    Outside,
}

/// Determine whether a `Coord` lies inside, outside, or on the boundary of a geometry.
///
/// # Examples
///
/// ```
/// use geo_3d::{polygon, coord};
/// use geo_3d::coordinate_position::{CoordinatePosition, CoordPos};
///
/// let square_poly = polygon![(x: 0.0, y: 0.0, z: 0.0), (x: 2.0, y: 0.0, z: 2.0), (x: 2.0, y: 2.0, z: 2.0), (x: 0.0, y: 2.0, z: 0.0), (x: 0.0, y: 0.0, z: 0.0)];
///
/// let inside_coord = coord! { x: 1.0, y: 1.0, z: 1.0 };
/// assert_eq!(square_poly.coordinate_position(&inside_coord), CoordPos::Inside);
///
/// let boundary_coord = coord! { x: 0.0, y: 1.0, z: 0.0 };
/// assert_eq!(square_poly.coordinate_position(&boundary_coord), CoordPos::OnBoundary);
///
/// let outside_coord = coord! { x: 5.0, y: 5.0, z: 5.0 };
/// assert_eq!(square_poly.coordinate_position(&outside_coord), CoordPos::Outside);
/// ```
pub trait CoordinatePosition<T: GeoNum> {
    fn coordinate_position(&self, coord: &Coord<T>) -> CoordPos {
        let mut is_inside = false;
        let mut boundary_count = 0;

        self.calculate_coordinate_position(coord, &mut is_inside, &mut boundary_count);

        // “The boundary of an arbitrary collection of geometries whose interiors are disjoint
        // consists of geometries drawn from the boundaries of the element geometries by
        // application of the ‘mod 2’ union rule”
        //
        // ― OpenGIS Simple Feature Access § 6.1.15.1
        if boundary_count % 2 == 1 {
            CoordPos::OnBoundary
        } else if is_inside {
            CoordPos::Inside
        } else {
            CoordPos::Outside
        }
    }

    // impls of this trait must:
    //  1. set `is_inside = true` if `coord` is contained within the Interior of any component.
    //  2. increment `boundary_count` for each component whose Boundary contains `coord`.
    fn calculate_coordinate_position(
        &self,
        coord: &Coord<T>,
        is_inside: &mut bool,
        boundary_count: &mut usize,
    );
}

impl<T: GeoNum> CoordinatePosition<T> for Coord<T> {
    fn calculate_coordinate_position(
        &self,
        coord: &Coord<T>,
        is_inside: &mut bool,
        _boundary_count: &mut usize,
    ) {
        if self == coord {
            *is_inside = true;
        }
    }
}

impl<T: GeoNum> CoordinatePosition<T> for Point<T> {
    fn calculate_coordinate_position(
        &self,
        coord: &Coord<T>,
        is_inside: &mut bool,
        _boundary_count: &mut usize,
    ) {
        if &self.0 == coord {
            *is_inside = true;
        }
    }
}

impl<T: GeoNum> CoordinatePosition<T> for Line<T> {
    fn calculate_coordinate_position(
        &self,
        coord: &Coord<T>,
        is_inside: &mut bool,
        boundary_count: &mut usize,
    ) {
        // degenerate line is a point
        if self.start == self.end {
            self.start
                .calculate_coordinate_position(coord, is_inside, boundary_count);
            return;
        }

        if coord == &self.start || coord == &self.end {
            *boundary_count += 1;
        } else if self.intersects(coord) {
            *is_inside = true;
        }
    }
}

impl<T: GeoNum> CoordinatePosition<T> for LineString<T> {
    fn calculate_coordinate_position(
        &self,
        coord: &Coord<T>,
        is_inside: &mut bool,
        boundary_count: &mut usize,
    ) {
        if self.0.len() < 2 {
            debug_assert!(false, "invalid line string with less than 2 coords");
            return;
        }

        if self.0.len() == 2 {
            // line string with two coords is just a line
            Line::new(self.0[0], self.0[1]).calculate_coordinate_position(
                coord,
                is_inside,
                boundary_count,
            );
            return;
        }

        // optimization: return early if there's no chance of an intersection
        // since self.0 is non-empty, safe to `unwrap`
        if !self.bounding_rect().unwrap().intersects(coord) {
            return;
        }

        // A closed linestring has no boundary, per SFS
        if !self.is_closed() {
            // since self.0 is non-empty, safe to `unwrap`
            if coord == self.0.first().unwrap() || coord == self.0.last().unwrap() {
                *boundary_count += 1;
                return;
            }
        }

        if self.intersects(coord) {
            // We've already checked for "Boundary" condition, so if there's an intersection at
            // this point, coord must be on the interior
            *is_inside = true
        }
    }
}

impl<T: GeoNum> CoordinatePosition<T> for Triangle<T> {
    fn calculate_coordinate_position(
        &self,
        coord: &Coord<T>,
        is_inside: &mut bool,
        boundary_count: &mut usize,
    ) {
        *is_inside = {
            let (u, v, w) = barycentric(*coord, self);
            u >= T::zero() && v >= T::zero() && w >= T::zero()
        };
        self.to_lines()
            .iter()
            .for_each(|l| {
                if *is_inside
                    && point_in_rect(*coord, l.start, l.end)
                    && coord.x != l.end.x
                {
                    *boundary_count += 1;
                }
            });
    }
}

impl<T: GeoNum> CoordinatePosition<T> for Rect<T> {
    fn calculate_coordinate_position(
        &self,
        coord: &Coord<T>,
        is_inside: &mut bool,
        boundary_count: &mut usize,
    ) {
        let mut boundary = false;

        let min = self.min();

        match coord.x.partial_cmp(&min.x).unwrap() {
            Ordering::Less => return,
            Ordering::Equal => boundary = true,
            Ordering::Greater => {}
        }
        match coord.y.partial_cmp(&min.y).unwrap() {
            Ordering::Less => return,
            Ordering::Equal => boundary = true,
            Ordering::Greater => {}
        }
        match coord.z.partial_cmp(&min.z).unwrap() {
            Ordering::Less => return,
            Ordering::Equal => boundary = true,
            Ordering::Greater => {}
        }

        let max = self.max();

        match max.x.partial_cmp(&coord.x).unwrap() {
            Ordering::Less => return,
            Ordering::Equal => boundary = true,
            Ordering::Greater => {}
        }
        match max.y.partial_cmp(&coord.y).unwrap() {
            Ordering::Less => return,
            Ordering::Equal => boundary = true,
            Ordering::Greater => {}
        }
        match max.z.partial_cmp(&coord.z).unwrap() {
            Ordering::Less => return,
            Ordering::Equal => boundary = true,
            Ordering::Greater => {}
        }

        if boundary {
            *boundary_count += 1;
        } else {
            *is_inside = true;
        }
    }
}

impl<T: GeoNum> CoordinatePosition<T> for MultiPoint<T> {
    fn calculate_coordinate_position(
        &self,
        coord: &Coord<T>,
        is_inside: &mut bool,
        _boundary_count: &mut usize,
    ) {
        if self.0.iter().any(|p| &p.0 == coord) {
            *is_inside = true;
        }
    }
}

impl<T: GeoNum + Default> CoordinatePosition<T> for Polygon<T> {
    fn calculate_coordinate_position(
        &self,
        coord: &Coord<T>,
        is_inside: &mut bool,
        boundary_count: &mut usize,
    ) {
        if self.is_empty() {
            return;
        }

        match inside_poly(self, coord.clone(), boundary_count) {
            CoordPos::OnBoundary => {
                boundary_count.add_assign(1);
            }
            CoordPos::Inside => *is_inside = true,
            CoordPos::Outside => {}
        }
    }
}

impl<T: GeoNum> CoordinatePosition<T> for MultiLineString<T> {
    fn calculate_coordinate_position(
        &self,
        coord: &Coord<T>,
        is_inside: &mut bool,
        boundary_count: &mut usize,
    ) {
        for line_string in &self.0 {
            line_string.calculate_coordinate_position(coord, is_inside, boundary_count);
        }
    }
}

impl<T: GeoNum + Default> CoordinatePosition<T> for MultiPolygon<T> {
    fn calculate_coordinate_position(
        &self,
        coord: &Coord<T>,
        is_inside: &mut bool,
        boundary_count: &mut usize,
    ) {
        for polygon in &self.0 {
            polygon.calculate_coordinate_position(coord, is_inside, boundary_count);
        }
    }
}

impl<T: GeoNum + Default> CoordinatePosition<T> for GeometryCollection<T> {
    fn calculate_coordinate_position(
        &self,
        coord: &Coord<T>,
        is_inside: &mut bool,
        boundary_count: &mut usize,
    ) {
        for geometry in self {
            geometry.calculate_coordinate_position(coord, is_inside, boundary_count);
        }
    }
}

impl<T: GeoNum + Default> CoordinatePosition<T> for Geometry<T> {
    crate::geometry_delegate_impl! {
        fn calculate_coordinate_position(
            &self,
            coord: &Coord<T>,
            is_inside: &mut bool,
            boundary_count: &mut usize) -> ();
    }
}

impl<T: GeoNum + Default> CoordinatePosition<T> for GeometryCow<'_, T> {
    crate::geometry_cow_delegate_impl! {
        fn calculate_coordinate_position(
            &self,
            coord: &Coord<T>,
            is_inside: &mut bool,
            boundary_count: &mut usize) -> ();
    }
}

// todo add `all_inside_poly` function to save on trianglulation
fn inside_poly<T: CoordNum + Default>(poly: &Polygon<T>, coord: Coord<T>, boundary_count: &mut usize) -> CoordPos {
    assert!(poly.exterior().is_closed());

    for linestring in poly.rings() {
        // LineString without points
        if linestring.0.is_empty() {
            return CoordPos::Outside;
        }
    }

    let segment = Line::new(
        coord,
        coord!(MAX),
    );

    // todo it should not be on boundary if it's on inner triangle boundary
    for triangle in poly.triangles() {
        if line_tri_intersect(segment, triangle) {
            boundary_count.add_assign(1);
        }
    }

    // is odd
    if *boundary_count % 2 != 0 {
        CoordPos::Inside
    } else {
        CoordPos::Outside
    }
}

/// Checks if a ray(`Line`) intersects a `Triangle`
fn line_tri_intersect<T: CoordNum>(segment: Line<T>, tri: Triangle<T>) -> bool {
    use super::kernels;

    let s1 = kernels::robust::RobustKernel::orient3d(segment.start, tri.0, tri.1, tri.2);
    let s2 = kernels::robust::RobustKernel::orient3d(segment.end, tri.0, tri.1, tri.2);
    // Test whether the two extermities of the segment
    // are on the same side of the supporting plane of the triangle
    if s1 == s2 {
        return false;
    }

    // Now we know that the segment 'straddles' the supporting plane.
    // We need to test whether the three tetrahedra formed
    // by the segment and the three edges of the triangle have the same orientation
    let s3 = kernels::robust::RobustKernel::orient3d(segment.start, segment.end, tri.0, tri.1);
    let s4 = kernels::robust::RobustKernel::orient3d(segment.start, segment.end, tri.1, tri.2);
    let s5 = kernels::robust::RobustKernel::orient3d(segment.start, segment.end, tri.2, tri.0);

    s3 == s4 && s4 == s5
}

/// Calculate the position of a `Coord` relative to a closed `LineString`.
#[deprecated(since = "0.30.0", note = "This is 2d, use `CoordinatePosition::coordinate_position` instead")]
pub(crate) fn coord_pos_relative_to_ring<T: GeoNum>(coord: Coord<T>, linestring: &LineString<T>) -> CoordPos {
    assert!(linestring.is_closed(), "The ring(LineString) must be closed.");

    // LineString without points
    if linestring.0.is_empty() {
        return CoordPos::Outside;
    }
    if linestring.0.len() == 1 {
        // If LineString has one point, it will not generate
        // any lines.  So, we handle this edge case separately.
        return if coord == linestring.0[0] {
            CoordPos::OnBoundary
        } else {
            CoordPos::Outside
        };
    }

    // Use winding number algorithm with on boundary short-cicuit
    // See: https://en.wikipedia.org/wiki/Point_in_polygon#Winding_number_algorithm
    let mut winding_number = 0;
    for line in linestring.lines() {
        // Edge Crossing Rules:
        //   1. an upward edge includes its starting endpoint, and excludes its final endpoint;
        //   2. a downward edge excludes its starting endpoint, and includes its final endpoint;
        //   3. horizontal edges are excluded
        //   4. the edge-ray intersection point must be strictly right of the coord.
        if line.start.y <= coord.y {
            if line.end.y >= coord.y {
                let o = T::Ker::orient2d(line.start, line.end, coord);
                if o == Orientation::CounterClockwise && line.end.y != coord.y {
                    winding_number += 1
                } else if o == Orientation::Collinear
                    && value_in_between(coord.x, line.start.x, line.end.x)
                {
                    return CoordPos::OnBoundary;
                }
            };
        } else if line.end.y <= coord.y {
            let o = T::Ker::orient2d(line.start, line.end, coord);
            if o == Orientation::Clockwise {
                winding_number -= 1
            } else if o == Orientation::Collinear
                && value_in_between(coord.x, line.start.x, line.end.x)
            {
                return CoordPos::OnBoundary;
            }
        }
    }
    if winding_number == 0 {
        CoordPos::Outside
    } else {
        CoordPos::Inside
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{line_string, point, polygon};
    use geo_types::coord;

    #[test]
    fn test_empty_poly() {
        let square_poly: Polygon<f64> = Polygon::new(LineString::new(vec![]), vec![]);
        assert_eq!(
            square_poly.coordinate_position(&Coord::zero()),
            CoordPos::Outside
        );
    }

    #[test]
    fn test_simple_poly() {
        let square_poly = polygon![
            (x: 0.0, y: 0.0, z: 0.0), (x: 2.0, y: 0.0, z: 2.0), (x: 2.0, y: 2.0, z: 2.0), (x: 0.0, y: 2.0, z: 0.0), (x: 0.0, y: 0.0, z: 0.0),
        ];

        let inside_coord = coord! { x: 1.0, y: 1.0, z: 1.0 };
        assert_eq!(
            square_poly.coordinate_position(&inside_coord),
            CoordPos::Inside
        );

        let vertex_coord = coord! { x: 0.0, y: 0.0, z: 0.0 };
        assert_eq!(
            square_poly.coordinate_position(&vertex_coord),
            CoordPos::OnBoundary
        );

        let boundary_coord = coord! { x: 0.0, y: 1.0, z: 0.0 };
        assert_eq!(
            square_poly.coordinate_position(&boundary_coord),
            CoordPos::OnBoundary
        );

        let outside_coord = coord! { x: 5.0, y: 5.0, z: 5.0 };
        assert_eq!(
            square_poly.coordinate_position(&outside_coord),
            CoordPos::Outside
        );
    }

    #[test]
    fn test_poly_interior() {
        let poly = polygon![
            exterior: [
                (x: 11., y: 11., z: 11.),
                (x: 20., y: 11., z: 20.),
                (x: 20., y: 20., z: 20.),
                (x: 11., y: 20., z: 11.),
                (x: 11., y: 11., z: 11.),
            ],
            interiors: [
                [
                    (x: 13., y: 13., z: 13.),
                    (x: 13., y: 17., z: 13.),
                    (x: 17., y: 17., z: 17.),
                    (x: 17., y: 13., z: 17.),
                    (x: 13., y: 13., z: 13.),
                ]
            ],
        ];

        let inside_hole = coord! { x: 14.0, y: 14.0, z: 14.0 };
        assert_eq!(poly.coordinate_position(&inside_hole), CoordPos::Outside);

        let outside_poly = coord! { x: 30.0, y: 30.0, z: 30.0 };
        assert_eq!(poly.coordinate_position(&outside_poly), CoordPos::Outside);

        let on_outside_border = coord! { x: 20.0, y: 15.0, z: 10.0 };
        assert_eq!(
            poly.coordinate_position(&on_outside_border),
            CoordPos::OnBoundary
        );

        let on_inside_border = coord! { x: 13.0, y: 15.0, z: 17.0 };
        assert_eq!(
            poly.coordinate_position(&on_inside_border),
            CoordPos::OnBoundary
        );

        let inside_coord = coord! { x: 12.0, y: 12.0, z: 12.0 };
        assert_eq!(poly.coordinate_position(&inside_coord), CoordPos::Inside);
    }

    #[test]
    fn test_simple_line() {
        let line = Line::new(
            coord![x: 0.0, y: 0.0, z: 0.0],
            coord![x: 10.0, y: 10.0, z: 10.0],
        );

        let start = coord! { x: 0.0, y: 0.0, z: 0.0 };
        assert_eq!(line.coordinate_position(&start), CoordPos::OnBoundary);

        let end = coord! { x: 10.0, y: 10.0, z: 10.0 };
        assert_eq!(line.coordinate_position(&end), CoordPos::OnBoundary);

        let interior = coord! { x: 5.0, y: 5.0, z: 5.0 };
        assert_eq!(line.coordinate_position(&interior), CoordPos::Inside);

        let outside = coord! { x: 6.0, y: 5.0, z: 6.0 };
        assert_eq!(line.coordinate_position(&outside), CoordPos::Outside);
    }

    #[test]
    fn test_degenerate_line() {
        let line = Line::new(
            coord![x: 0.0, y: 0.0, z: 0.0],
            coord![x: 0.0, y: 0.0, z: 0.0],
        );

        let start = coord! { x: 0.0, y: 0.0, z: 0.0 };
        assert_eq!(line.coordinate_position(&start), CoordPos::Inside);

        let outside = coord! { x: 10.0, y: 10.0, z: 10.0 };
        assert_eq!(line.coordinate_position(&outside), CoordPos::Outside);
    }

    #[test]
    fn test_point() {
        let p1 = point![x: 2.0, y: 0.0, z: 0.5];

        let c1 = coord! { x: 2.0, y: 0.0, z: 0.5 };
        let c2 = coord! { x: 3.0, y: 3.0, z: 0.0 };

        assert_eq!(p1.coordinate_position(&c1), CoordPos::Inside);
        assert_eq!(p1.coordinate_position(&c2), CoordPos::Outside);

        assert_eq!(c1.coordinate_position(&c1), CoordPos::Inside);
        assert_eq!(c1.coordinate_position(&c2), CoordPos::Outside);
    }

    #[test]
    fn test_simple_line_string() {
        let line_string = line_string![(x: 0.0, y: 0.0, z: 0.0), (x: 1.0, y: 1.0, z: 1.0), (x: 2.0, y: 0.0, z: 2.0), (x: 3.0, y: 0.0, z: 3.0)];

        let start = Coord::zero();
        assert_eq!(
            line_string.coordinate_position(&start),
            CoordPos::OnBoundary
        );

        let midpoint = coord! { x: 0.5, y: 0.5, z: 0.5 };
        assert_eq!(line_string.coordinate_position(&midpoint), CoordPos::Inside);

        let vertex = coord! { x: 2.0, y: 0.0, z: 2.0 };
        assert_eq!(line_string.coordinate_position(&vertex), CoordPos::Inside);

        let end = coord! { x: 3.0, y: 0.0, z: 3.0 };
        assert_eq!(line_string.coordinate_position(&end), CoordPos::OnBoundary);

        let outside = coord! { x: 3.0, y: 1.0, z: 0.6 };
        assert_eq!(line_string.coordinate_position(&outside), CoordPos::Outside);
    }

    #[test]
    fn test_degenerate_line_strings() {
        let line_string = line_string![(x: 0.0, y: 0.0, z: 0.0), (x: 0.0, y: 0.0, z: 0.0)];

        let start = Coord::zero();
        assert_eq!(line_string.coordinate_position(&start), CoordPos::Inside);

        let line_string = line_string![(x: 0.0, y: 0.0, z: 0.0), (x: 2.0, y: 0.0, z: 2.0)];

        let start = Coord::zero();
        assert_eq!(
            line_string.coordinate_position(&start),
            CoordPos::OnBoundary
        );
    }

    #[test]
    fn test_closed_line_string() {
        let line_string = line_string![
            (x: 0.0, y: 0.0, z: 0.0), (x: 1.0, y: 1.0, z: 1.0), (x: 2.0, y: 0.0, z: 2.0), (x: 3.0, y: 2.0, z: 3.0), (x: 0.0, y: 2.0, z: 0.0), (x: 0.0, y: 0.0, z: 0.0),
        ];

        // sanity check
        assert!(line_string.is_closed());

        // closed line strings have no boundary
        let start = Coord::zero();
        assert_eq!(line_string.coordinate_position(&start), CoordPos::Inside);

        let midpoint = coord! { x: 0.5, y: 0.5, z: 0.5 };
        assert_eq!(line_string.coordinate_position(&midpoint), CoordPos::Inside);

        let outside = coord! { x: 3.0, y: 1.0, z: 3.0 };
        assert_eq!(line_string.coordinate_position(&outside), CoordPos::Outside);
    }

    #[test]
    fn test_boundary_rule() {
        let multi_line_string = MultiLineString::new(vec![
            // first two lines have same start point but different end point
            line_string![(x: 0.0, y: 0.0, z: 0.0), (x: 1.0, y: 1.0, z: 1.0)],
            line_string![(x: 0.0, y: 0.0, z: 0.0), (x: -1.0, y: -1.0, z: -1.0)],
            // third line has its own start point, but it's end touches the middle of first line
            line_string![(x: 0.0, y: 1.0, z: 0.0), (x: 0.5, y: 0.5, z: 0.5)],
            // fourth and fifth have independent start points, but both end at the middle of the
            // second line
            line_string![(x: 0.0, y: -1.0, z: 0.0), (x: -0.5, y: -0.5, z: -0.5)],
            line_string![(x: 0.0, y: -2.0, z: 0.0), (x: -0.5, y: -0.5, z: -0.5)],
        ]);

        let outside_of_all = coord! { x: 123.0, y: 123.0, z: 123.0 };
        assert_eq!(
            multi_line_string.coordinate_position(&outside_of_all),
            CoordPos::Outside
        );

        let end_of_one_line = coord! { x: -1.0, y: -1.0, z: -1.0 };
        assert_eq!(
            multi_line_string.coordinate_position(&end_of_one_line),
            CoordPos::OnBoundary
        );

        // in boundary of first and second, so considered *not* in the boundary by mod 2 rule
        let shared_start = Coord::zero();
        assert_eq!(
            multi_line_string.coordinate_position(&shared_start),
            CoordPos::Outside
        );

        // *in* the first line, on the boundary of the third line
        let one_end_plus_midpoint = coord! { x: 0.5, y: 0.5, z: 0.5 };
        assert_eq!(
            multi_line_string.coordinate_position(&one_end_plus_midpoint),
            CoordPos::OnBoundary
        );

        // *in* the first line, on the *boundary* of the fourth and fifth line
        let two_ends_plus_midpoint = coord! { x: -0.5, y: -0.5, z: -0.5 };
        assert_eq!(
            multi_line_string.coordinate_position(&two_ends_plus_midpoint),
            CoordPos::Inside
        );
    }

    #[test]
    fn test_rect() {
        let rect = Rect::new((0.0, 0.0, 0.0), (10.0, 10.0, 10.0));
        assert_eq!(
            rect.coordinate_position(&coord! { x: 5.0, y: 5.0, z: 5.0 }),
            CoordPos::Inside
        );
        assert_eq!(
            rect.coordinate_position(&coord! { x: 0.0, y: 5.0, z: 0.0 }),
            CoordPos::OnBoundary
        );
        assert_eq!(
            rect.coordinate_position(&coord! { x: 15.0, y: 15.0, z: 15.0 }),
            CoordPos::Outside
        );
    }

    #[test]
    fn test_triangle() {
        let triangle = Triangle::new(
            (0.0, 0.0, 0.0).into(),
            (5.0, 10.0, 5.0).into(),
            (10.0, 0.0, 10.0).into(),
        );
        assert_eq!(
            triangle.coordinate_position(&coord! { x: 5.0, y: 5.0, z: 4.9 }),
            CoordPos::Inside
        );
        assert_eq!(
            triangle.coordinate_position(&coord! { x: 5.0, y: 5.0, z: 5.0 }),
            CoordPos::OnBoundary
        );
        assert_eq!(
            triangle.coordinate_position(&coord! { x: 2.5, y: 5.0, z: 2.5 }),
            CoordPos::OnBoundary
        );
        assert_eq!(
            triangle.coordinate_position(&coord! { x: 2.49, y: 5.0, z: 2.5 }),
            CoordPos::Outside
        );
    }

    #[test]
    fn test_collection() {
        let triangle = Triangle::new(
            (0.0, 0.0, 0.0).into(),
            (5.0, 10.0, 15.0).into(),
            (10.0, 0.0, 10.0).into(),
        );
        let rect = Rect::new((0.0, 0.0, 0.0), (10.0, 10.0, 10.0));
        let collection = GeometryCollection::new(vec![triangle.into(), rect.into()]);

        //  outside of both
        assert_eq!(
            collection.coordinate_position(&coord! { x: 15.0, y: 15.0, z: 15.0 }),
            CoordPos::Outside
        );

        // inside both
        assert_eq!(
            collection.coordinate_position(&coord! { x: 5.0, y: 5.0, z: 5.0 }),
            CoordPos::Inside
        );

        // inside one, boundary of other
        assert_eq!(
            collection.coordinate_position(&coord! { x: 2.5, y: 5.0, z: 2.5 }),
            CoordPos::OnBoundary
        );

        //  boundary of both
        assert_eq!(
            collection.coordinate_position(&coord! { x: 5.0, y: 10.0, z: 15.0 }),
            CoordPos::Outside
        );
    }
}
