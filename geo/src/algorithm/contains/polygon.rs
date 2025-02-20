use super::{impl_contains_from_relate, impl_contains_geometry_for, Contains, Contains2D};
use crate::coordinate_position::{coord_pos_relative_to_ring, CoordPos};
use crate::geometry::*;
use crate::{GeoFloat, GeoNum};
use crate::{HasDimensions, Relate};

// ┌─────────────────────────────┐
// │      2D Implementations     │
// └─────────────────────────────┘

impl<T> Contains2D<Coord<T>> for Polygon<T>
where
    T: GeoNum,
{
    fn contains_2d(&self, coord: &Coord<T>) -> bool {
        use crate::coordinate_position::CoordPos;

        coord_position_in_poly(self, coord) == CoordPos::Inside
    }
}

fn coord_position_in_poly<T: GeoNum>(poly: &Polygon<T>, coord: &Coord<T>) -> CoordPos {
    let mut is_inside = false;
    let mut boundary_count = 0;

    calculate_coord_poly_position(poly, coord, &mut is_inside, &mut boundary_count);

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

fn calculate_coord_poly_position<T: GeoNum>(
    poly: &Polygon<T>,
    coord: &Coord<T>,
    is_inside: &mut bool,
    boundary_count: &mut usize,
) {
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

impl<T: GeoNum> Contains2D<Coord<T>> for MultiPolygon<T> {
    fn contains_2d(&self, coord: &Coord<T>) -> bool {
        self.iter().any(|poly| poly.contains_2d(coord))
    }
}

impl<T: GeoNum> Contains2D<Point<T>> for MultiPolygon<T> {
    fn contains_2d(&self, p: &Point<T>) -> bool {
        self.contains_2d(&p.0)
    }
}

impl<T: GeoNum> Contains2D<MultiPoint<T>> for MultiPolygon<T> {
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
        use crate::coordinate_position::{CoordPos, CoordinatePosition};

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
        rhs.iter().all(|point| self.contains(point))
    }
}

impl<F: GeoFloat> Contains<Line<F>> for MultiPolygon<F> {
    fn contains(&self, rhs: &Line<F>) -> bool {
        rhs.relate(self).is_within()
    }
}

impl<F: GeoFloat> Contains<LineString<F>> for MultiPolygon<F> {
    fn contains(&self, rhs: &LineString<F>) -> bool {
        rhs.relate(self).is_within()
    }
}

impl<F: GeoFloat> Contains<MultiLineString<F>> for MultiPolygon<F> {
    fn contains(&self, rhs: &MultiLineString<F>) -> bool {
        rhs.relate(self).is_within()
    }
}

impl<F: GeoFloat> Contains<Polygon<F>> for MultiPolygon<F> {
    fn contains(&self, rhs: &Polygon<F>) -> bool {
        rhs.relate(self).is_within()
    }
}

impl<F: GeoFloat> Contains<MultiPolygon<F>> for MultiPolygon<F> {
    fn contains(&self, rhs: &MultiPolygon<F>) -> bool {
        rhs.relate(self).is_within()
    }
}

impl<F: GeoFloat> Contains<GeometryCollection<F>> for MultiPolygon<F> {
    fn contains(&self, rhs: &GeometryCollection<F>) -> bool {
        rhs.relate(self).is_within()
    }
}

impl<F: GeoFloat> Contains<Rect<F>> for MultiPolygon<F> {
    fn contains(&self, rhs: &Rect<F>) -> bool {
        rhs.relate(self).is_within()
    }
}

impl<F: GeoFloat> Contains<Triangle<F>> for MultiPolygon<F> {
    fn contains(&self, rhs: &Triangle<F>) -> bool {
        rhs.relate(self).is_within()
    }
}
