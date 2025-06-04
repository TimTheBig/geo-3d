use super::{impl_contains_from_relate, impl_contains_geometry_for, Contains, ContainsXY};
use crate::coordinate_position::CoordPos;
use crate::geometry::*;
use crate::{GeoNum, CoordNum};
use crate::{HasDimensions, Relate};

// ┌─────────────────────────────┐
// │      2D Implementations     │
// └─────────────────────────────┘

impl<T: GeoNum> ContainsXY<Coord<T>> for Polygon<T> {
    /// Checks whether a `Coord` is inside a `Polygon`
    fn contains_2d(&self, coord: &Coord<T>) -> bool {
        xy_position_in_poly(self, coord) == CoordPos::Inside
    }
}

impl<T: GeoNum> ContainsXY<Point<T>> for Polygon<T> {
    /// Checks whether a `Point` is inside a `Polygon`
    fn contains_2d(&self, point: &Point<T>) -> bool {
        self.contains_2d(&point.0)
    }
}

/// Returns whether the `Coord` is CoordPos::OnBoundary, CoordPos::Inside, or CoordPos::Outside.
fn xy_position_in_poly<T: GeoNum>(poly: &Polygon<T>, coord: &Coord<T>) -> CoordPos {
    let mut is_inside = false;
    let mut boundary_count = 0;

    calculate_coord_poly_position_xy(poly, coord, &mut is_inside, &mut boundary_count);

    // “The boundary of an arbitrary collection of geometries whose interiors are disjoint
    // consists of geometries drawn from the boundaries of the element geometries by
    // application of the ‘mod 2’ union rule”
    if boundary_count % 2 == 1 {
        CoordPos::OnBoundary
    } else if is_inside {
        CoordPos::Inside
    } else {
        CoordPos::Outside
    }
}

/// Calculate a `Coord`s 2d possiton in a `Polygon`
#[expect(deprecated)]
fn calculate_coord_poly_position_xy<T: GeoNum>(
    poly: &Polygon<T>,
    coord: &Coord<T>,
    is_inside: &mut bool,
    boundary_count: &mut usize,
) {
    use crate::coordinate_position::coord_pos_relative_to_ring;

    if poly.is_empty() {
        return;
    }

    match coord_pos_relative_to_ring(*coord, poly.exterior()) {
        CoordPos::Outside => {}
        CoordPos::OnBoundary => {
            *boundary_count += 1;
        }
        CoordPos::Inside => {
            for hole in poly.interiors() {
                match coord_pos_relative_to_ring(*coord, hole) {
                    CoordPos::Outside => {}
                    CoordPos::OnBoundary => {
                        *boundary_count += 1;
                        return;
                    }
                    CoordPos::Inside => {
                        return;
                    }
                }
            }
            // the coord is *outside* the interior holes, so it's *inside* the polygon
            *is_inside = true;
        }
    }
}

impl<T: GeoNum> ContainsXY<Coord<T>> for MultiPolygon<T> {
    /// Checks whether any `Polygon` contains `Coord`
    fn contains_2d(&self, coord: &Coord<T>) -> bool {
        self.iter().any(|poly| poly.contains_2d(coord))
    }
}

impl<T: GeoNum> ContainsXY<Point<T>> for MultiPolygon<T> {
    /// Checks whether any `Polygon` contains `Point`
    fn contains_2d(&self, p: &Point<T>) -> bool {
        self.contains_2d(&p.0)
    }
}

impl<T: GeoNum> ContainsXY<MultiPoint<T>> for MultiPolygon<T> {
    /// Checks if all points of `rhs` are contained in `self`.
    ///
    /// ## Cost
    /// This is very expensive as it has to check `MultiPolygon::contains_2d::<Point>`,
    /// an expensive function in it's own right, `n` times.
    fn contains_2d(&self, rhs: &MultiPoint<T>) -> bool {
        if self.is_empty() || rhs.is_empty() {
            return false;
        }
        rhs.iter().all(|point| self.contains_2d(point))
    }
}

// ┌──────────────────────────────────┐
// │   Implementations for Polygon    │
// └──────────────────────────────────┘

impl<T: GeoNum> Contains<Coord<T>> for Polygon<T> {
    fn contains(&self, coord: &Coord<T>) -> bool {
        use crate::coordinate_position::CoordinatePosition;

        self.coordinate_position(coord) == CoordPos::Inside
    }
}

impl<T: GeoNum> Contains<Point<T>> for Polygon<T> {
    fn contains(&self, p: &Point<T>) -> bool {
        self.contains(&p.0)
    }
}

impl_contains_from_relate!(Polygon<T>, [Line<T>, LineString<T>, Polygon<T>, MultiPoint<T>, MultiLineString<T>, MultiPolygon<T>, GeometryCollection<T>, Rect<T>, Triangle<T>]);
impl_contains_geometry_for!(Polygon<T>);

// ┌──────────────────────────────────┐
// │ Implementations for MultiPolygon │
// └──────────────────────────────────┘

impl<T: GeoNum> Contains<Coord<T>> for MultiPolygon<T> {
    fn contains(&self, coord: &Coord<T>) -> bool {
        self.iter().any(|poly| poly.contains(coord))
    }
}

impl<T: GeoNum> Contains<Point<T>> for MultiPolygon<T> {
    fn contains(&self, p: &Point<T>) -> bool {
        self.contains(&p.0)
    }
}

impl<T: GeoNum> Contains<MultiPoint<T>> for MultiPolygon<T> {
    fn contains(&self, rhs: &MultiPoint<T>) -> bool {
        if self.is_empty() || rhs.is_empty() {
            return false;
        }
        // rhs.iter().all(|point| self.contains(point))
        multipolygon_contains_many(self, &rhs.0)
    }
}

impl<F: GeoNum> Contains<Line<F>> for MultiPolygon<F> {
    fn contains(&self, rhs: &Line<F>) -> bool {
        rhs.relate(self).is_within()
    }
}

impl<F: GeoNum> Contains<LineString<F>> for MultiPolygon<F> {
    fn contains(&self, rhs: &LineString<F>) -> bool {
        rhs.relate(self).is_within()
    }
}

impl<F: GeoNum> Contains<MultiLineString<F>> for MultiPolygon<F> {
    fn contains(&self, rhs: &MultiLineString<F>) -> bool {
        rhs.relate(self).is_within()
    }
}

impl<F: GeoNum> Contains<Polygon<F>> for MultiPolygon<F> {
    fn contains(&self, rhs: &Polygon<F>) -> bool {
        rhs.relate(self).is_within()
    }
}

impl<F: GeoNum> Contains<MultiPolygon<F>> for MultiPolygon<F> {
    fn contains(&self, rhs: &MultiPolygon<F>) -> bool {
        rhs.relate(self).is_within()
    }
}

impl<F: GeoNum> Contains<GeometryCollection<F>> for MultiPolygon<F> {
    fn contains(&self, rhs: &GeometryCollection<F>) -> bool {
        rhs.relate(self).is_within()
    }
}

impl<F: GeoNum> Contains<Rect<F>> for MultiPolygon<F> {
    fn contains(&self, rhs: &Rect<F>) -> bool {
        rhs.relate(self).is_within()
    }
}

impl<F: GeoNum> Contains<Triangle<F>> for MultiPolygon<F> {
    fn contains(&self, rhs: &Triangle<F>) -> bool {
        rhs.relate(self).is_within()
    }
}

fn multipolygon_contains_many<T: CoordNum, C: Into<Coord<T>> + Copy>(mpoly: &MultiPolygon<T>, coords: &[C]) -> bool {
    mpoly.iter().any(|poly| polygon_contains_many(poly, coords))
}

/// Function to save on trianglulation
fn polygon_contains_many<T: CoordNum, C: Into<Coord<T>> + Copy>(poly: &Polygon<T>, coords: &[C]) -> bool {
    use crate::Triangulate;
    // todo bounding box skip optamizeation
    // use super::BoundingRect;

    assert!(poly.exterior().is_closed());

    if poly.coords_count() < 3 {
        return false;
    }

    let tri_poly = poly.triangles();

    fn inside_tri_poly<T: CoordNum + Default>(poly: &[Triangle<T>], coord: Coord<T>) -> bool {
        let mut boundary_count: usize = 0;

        let segment = Line::new(
            coord,
            crate::coord!(MAX),
        );

        // todo it should not be on boundary if it's on inner triangle boundary
        for triangle in poly {
            if line_tri_intersect(segment, triangle) {
                boundary_count += 1;
            }
        }

        // is odd
        if boundary_count % 2 != 0 {
            // Inside
            true
        } else {
            // Outside
            false
        }
    }

    /// Checks if a ray(`Line`) intersects a `Triangle`
    fn line_tri_intersect<T: CoordNum>(segment: Line<T>, tri: &Triangle<T>) -> bool {
        use crate::kernels::{robust::RobustKernel, Kernel};

        let s1 = RobustKernel::orient3d(segment.start, tri.0, tri.1, tri.2);
        let s2 = RobustKernel::orient3d(segment.end, tri.0, tri.1, tri.2);
        // Test whether the two extermities of the segment
        // are on the same side of the supporting plane of the triangle
        if s1 == s2 {
            return false;
        }

        // Now we know that the segment 'straddles' the supporting plane.
        // We need to test whether the three tetrahedra formed
        // by the segment and the three edges of the triangle have the same orientation
        let s3 = RobustKernel::orient3d(segment.start, segment.end, tri.0, tri.1);
        let s4 = RobustKernel::orient3d(segment.start, segment.end, tri.1, tri.2);
        let s5 = RobustKernel::orient3d(segment.start, segment.end, tri.2, tri.0);

        s3 == s4 && s4 == s5
    }

    coords.iter().all(|c| inside_tri_poly(&tri_poly, Into::<Coord<T>>::into(*c)))
}
