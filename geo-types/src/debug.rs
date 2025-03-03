use core::fmt::{Debug, Formatter};

use crate::geometry::*;
use crate::CoordNum;

impl<T: CoordNum> Debug for Coord<T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> core::fmt::Result {
        write!(f, "COORD({x:?} {y:?} {z:?})", x = self.x, y = self.y, z = self.z)
    }
}

impl<T: CoordNum> Debug for Point<T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> core::fmt::Result {
        write!(f, "POINT({x:?} {y:?} {z:?})", x = self.x(), y = self.y(), z = self.z())
    }
}

impl<T: CoordNum> Debug for Line<T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> core::fmt::Result {
        write!(f, "LINE")?;
        write_coord_seq(f, [self.start, self.end].iter())
    }
}

impl<T: CoordNum> Debug for LineString<T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> core::fmt::Result {
        write!(f, "LINESTRING")?;
        if self.0.is_empty() {
            write!(f, " ")?;
        }
        write_coord_seq(f, self.0.iter())?;
        Ok(())
    }
}

impl<T: CoordNum> Debug for Polygon<T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> core::fmt::Result {
        write!(f, "POLYGON")?;
        if self.exterior().0.is_empty() && self.interiors().is_empty() {
            write!(f, " ")?;
        }
        write_polygon_inner(f, self)
    }
}

impl<T: CoordNum> Debug for MultiPoint<T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> core::fmt::Result {
        write!(f, "MULTIPOINT")?;
        if self.0.is_empty() {
            write!(f, " ")?;
        }
        write_coord_seq(f, self.0.iter().map(|p| &p.0))
    }
}

impl<T: CoordNum> Debug for MultiLineString<T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> core::fmt::Result {
        write!(f, "MULTILINESTRING")?;
        let mut line_strings = self.0.iter();
        let Some(first) = line_strings.next() else {
            return write!(f, " EMPTY");
        };
        write!(f, "(")?;
        write_coord_seq(f, first.0.iter())?;
        for line_string in line_strings {
            write!(f, ",")?;
            write_coord_seq(f, line_string.0.iter())?;
        }
        write!(f, ")")
    }
}
impl<T: CoordNum> Debug for MultiPolygon<T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> core::fmt::Result {
        write!(f, "MULTIPOLYGON")?;
        let mut polygons = self.0.iter();
        let Some(first) = polygons.next() else {
            return write!(f, " EMPTY");
        };
        write!(f, "(")?;
        write_polygon_inner(f, first)?;
        for polygon in polygons {
            write!(f, ",")?;
            write_polygon_inner(f, polygon)?;
        }
        write!(f, ")")
    }
}

impl<T: CoordNum> Debug for Rect<T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> core::fmt::Result {
        write!(f, "RECT")?;
        write_coord_seq(f, [self.min(), self.max()].iter())
    }
}

impl<T: CoordNum> Debug for Triangle<T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> core::fmt::Result {
        write!(f, "TRIANGLE")?;
        write_coord_seq(f, [self.0, self.1, self.2].iter())
    }
}

impl<T: CoordNum> Debug for Geometry<T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> core::fmt::Result {
        match self {
            Geometry::Point(inner) => inner.fmt(f),
            Geometry::Line(inner) => inner.fmt(f),
            Geometry::LineString(inner) => inner.fmt(f),
            Geometry::Polygon(inner) => inner.fmt(f),
            Geometry::MultiPoint(inner) => inner.fmt(f),
            Geometry::MultiLineString(inner) => inner.fmt(f),
            Geometry::MultiPolygon(inner) => inner.fmt(f),
            Geometry::GeometryCollection(inner) => inner.fmt(f),
            Geometry::Rect(inner) => inner.fmt(f),
            Geometry::Triangle(inner) => inner.fmt(f),
        }
    }
}

impl<T: CoordNum> Debug for GeometryCollection<T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> core::fmt::Result {
        write!(f, "GEOMETRYCOLLECTION")?;
        let mut geometries = self.0.iter();
        let Some(first) = geometries.next() else {
            return write!(f, " EMPTY");
        };
        write!(f, "({first:?}")?;
        for geometry in geometries {
            write!(f, ", {geometry:?}")?;
        }
        write!(f, ")")
    }
}

fn write_coord_seq<'a, T: CoordNum + 'a>(
    f: &mut Formatter<'_>,
    mut coords: impl Iterator<Item = &'a Coord<T>>,
) -> core::fmt::Result {
    let Some(coord) = coords.next() else {
        write!(f, "EMPTY")?;
        return Ok(());
    };
    write!(f, "({x:?} {y:?} {z:?}", x = coord.x, y = coord.y, z = coord.z)?;
    for coord in coords {
        write!(f, ",{x:?} {y:?} {z:?}", x = coord.x, y = coord.y, z = coord.z)?;
    }
    write!(f, ")")
}

fn write_polygon_inner<T: CoordNum>(
    f: &mut Formatter<'_>,
    polygon: &Polygon<T>,
) -> core::fmt::Result {
    if polygon.exterior().0.is_empty() {
        let mut interiors = polygon.interiors().iter();
        let Some(interior) = interiors.next() else {
            write!(f, "EMPTY")?;
            return Ok(());
        };

        // Invalid polygon - having interiors but no exterior!
        // Still, we should try to print something meaningful.
        write!(f, "(EMPTY,")?;
        write_coord_seq(f, interior.0.iter())?;
        for interior in interiors {
            write!(f, ",")?;
            write_coord_seq(f, interior.0.iter())?;
        }
        write!(f, ")")?;
    } else {
        write!(f, "(")?;
        write_coord_seq(f, polygon.exterior().0.iter())?;
        for interior in polygon.interiors().iter() {
            write!(f, ",")?;
            write_coord_seq(f, interior.0.iter())?;
        }
        write!(f, ")")?;
    }
    Ok(())
}

#[cfg(feature = "std")]
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coord() {
        let coord = Coord { x: 1.0, y: 2.0, z: 3.0 };
        assert_eq!("COORD(1.0 2.0 3.0)", format!("{coord:?}"));
    }
    #[test]
    fn point() {
        let point = Point::new(1.0, 2.0, 3.0);
        assert_eq!("POINT(1.0 2.0 3.0)", format!("{point:?}"));
    }
    #[test]
    fn line() {
        let line_string = Line::new((1., 2., 3.), (3., 4., 5.));
        assert_eq!("LINE(1.0 2.0 3.0,3.0 4.0 5.0)", format!("{line_string:?}"));
    }
    #[test]
    fn line_string() {
        let line_string = LineString::new(vec![(1., 2., 3.).into(), (3., 4., 5.).into()]);
        assert_eq!("LINESTRING(1.0 2.0 3.0,3.0 4.0 5.0)", format!("{line_string:?}"));
    }
    #[test]
    fn line_string_with_single_element() {
        let line_string = LineString::new(vec![(1., 2., 3.).into()]);
        assert_eq!("LINESTRING(1.0 2.0 3.0)", format!("{line_string:?}"));
    }
    #[test]
    fn empty_line_string() {
        let line_string = LineString::<f32>::new(vec![]);
        assert_eq!("LINESTRING EMPTY", format!("{line_string:?}"));
    }
    #[test]
    fn polygon_no_holes() {
        let polygon = wkt!(POLYGON((1. 2. 3.,3. 4. 5.,5. 6. 7.)));
        assert_eq!("POLYGON((1.0 2.0 3.0,3.0 4.0 5.0,5.0 6.0 7.0,1.0 2.0 3.0))", format!("{polygon:?}"));
    }
    #[test]
    fn polygon_with_hole() {
        let polygon = wkt!(POLYGON(
            (1. 1. 1.,10. 1. 10.,10. 10. 10.,1. 10. 1.,1. 1. 1.),
            (3. 3. 3.,7. 3. 7.,7. 7. 7.,3. 7. 3.,3. 3. 3.)
        ));
        assert_eq!(
            "POLYGON((1.0 1.0 1.0,10.0 1.0 10.0,10.0 10.0 10.0,1.0 10.0 1.0,1.0 1.0 1.0),(3.0 3.0 3.0,7.0 3.0 7.0,7.0 7.0 7.0,3.0 7.0 3.0,3.0 3.0 3.0))",
            format!("{polygon:?}")
        );
    }
    #[test]
    fn polygon_with_multiple_holes() {
        let polygon = wkt!(POLYGON(
            (0.0 0.0 0.0,10.0 0.0 10.0,10.0 10.0 10.0,0.0 10.0 0.0,0.0 0.0 0.0),
            (2.0 2.0 2.0,4.0 2.0 4.0,4.0 4.0 4.0,2.0 4.0 2.0,2.0 2.0 2.0),
            (6.0 6.0 6.0,8.0 6.0 8.0,8.0 8.0 8.0,6.0 8.0 6.0,6.0 6.0 6.0)
        ));
        assert_eq!(
            "POLYGON((0.0 0.0 0.0,10.0 0.0 10.0,10.0 10.0 10.0,0.0 10.0 0.0,0.0 0.0 0.0),(2.0 2.0 2.0,4.0 2.0 4.0,4.0 4.0 4.0,2.0 4.0 2.0,2.0 2.0 2.0),(6.0 6.0 6.0,8.0 6.0 8.0,8.0 8.0 8.0,6.0 8.0 6.0,6.0 6.0 6.0))",
            format!("{polygon:?}")
        );
    }
    #[test]
    fn invalid_polygon_interior_but_no_exterior() {
        // Not a valid polygon, but we should still have reasonable debug output - note this is *not* valid WKT
        let interior = LineString::new(vec![(1.0, 2.0, 3.0).into()]);
        let polygon = Polygon::new(LineString::new(vec![]), vec![interior]);
        assert_eq!("POLYGON(EMPTY,(1.0 2.0 3.0))", format!("{polygon:?}"));
    }
    #[test]
    fn empty_polygon() {
        let polygon: Polygon = wkt!(POLYGON EMPTY);
        assert_eq!("POLYGON EMPTY", format!("{polygon:?}"));
    }
    #[test]
    fn multi_point_empty() {
        let multi_point: MultiPoint = wkt!(MULTIPOINT EMPTY);
        assert_eq!("MULTIPOINT EMPTY", format!("{multi_point:?}"));
    }
    #[test]
    fn multi_point_one_point() {
        let multi_point = wkt!(MULTIPOINT(1.0 2.0 3.0));
        assert_eq!("MULTIPOINT(1.0 2.0 3.0)", format!("{multi_point:?}"));
    }
    #[test]
    fn multi_point_three_points() {
        let multi_point = wkt!(MULTIPOINT(1.0 2.0 3.0,3.0 4.0 5.0,5.0 6.0 7.0));
        assert_eq!("MULTIPOINT(1.0 2.0 3.0,3.0 4.0 5.0,5.0 6.0 7.0)", format!("{multi_point:?}"));
    }
    #[test]
    fn multilinestring_empty() {
        let multi_line_string: MultiLineString = wkt!(MULTILINESTRING EMPTY);
        assert_eq!("MULTILINESTRING EMPTY", format!("{multi_line_string:?}"));
    }

    #[test]
    fn multi_line_string_one_line() {
        let multi_line_string = wkt!(MULTILINESTRING((1.0 2.0 3.0, 3.0 4.0 5.0, 5.0 6.0 7.0)));
        assert_eq!(
            "MULTILINESTRING((1.0 2.0 3.0,3.0 4.0 5.0,5.0 6.0 7.0))",
            format!("{multi_line_string:?}")
        );
    }

    #[test]
    fn multi_line_string_multiple_lines() {
        let multi_line_string = wkt!(MULTILINESTRING(
            (1.0 2.0 2.5, 3.0 4.0 4.5, 5.0 6.0 6.5),
            (7.0 8.0 9.5, 9.0 10.0 10.5, 11.0 12.0 12.5)
        ));
        assert_eq!(
            "MULTILINESTRING((1.0 2.0 2.5,3.0 4.0 4.5,5.0 6.0 6.5),(7.0 8.0 8.5,9.0 10.0 10.5,11.0 12.0 12.5))",
            format!("{multi_line_string:?}")
        );
    }

    #[test]
    fn multi_line_string_multiple_lines_with_empty() {
        let multi_line_string = wkt!(MULTILINESTRING(
            (1.0 2.0 2.5, 3.0 4.0 4.5, 5.0 6.0 6.5),
            EMPTY,
            (7.0 8.0 8.5, 9.0 10.0 10.5, 11.0 12.0 12.5)
        ));
        assert_eq!(
            "MULTILINESTRING((1.0 2.0 2.5,3.0 4.0 4.5,5.0 6.0 6.5),EMPTY,(7.0 8.0 8.5,9.0 10.0 10.5,11.0 12.0 12.5))",
            format!("{multi_line_string:?}")
        );
    }
    #[test]
    fn multi_polygon_empty() {
        let multi_polygon: MultiPolygon = wkt!(MULTIPOLYGON EMPTY);
        assert_eq!("MULTIPOLYGON EMPTY", format!("{multi_polygon:?}"));
    }

    #[test]
    fn multi_polygon_one_polygon() {
        let multi_polygon = wkt!(MULTIPOLYGON(
            ((1.0 2.0 2.5, 3.0 4.0 4.5, 5.0 6.0 6.5, 1.0 2.0 2.5))
        ));
        assert_eq!(
            "MULTIPOLYGON(((1.0 2.0 2.5,3.0 4.0 4.5,5.0 6.0 6.5,1.0 2.0 2.5)))",
            format!("{multi_polygon:?}")
        );
    }

    #[test]
    fn multi_polygon_multiple_polygons() {
        let multi_polygon = wkt!(MULTIPOLYGON(
            ((1.0 2.0 2.5, 3.0 4.0 4.5, 5.0 6.0 6.5, 1.0 2.0 2.5)),
            ((7.0 8.0 8.5, 9.0 10.0 10.5, 11.0 12.0 12.5, 7.0 8.0 8.5))
        ));
        assert_eq!(
            "MULTIPOLYGON(((1.0 2.0 2.5, 3.0 4.0 4.5, 5.0 6.0 6.5, 1.0 2.0 2.5)),((7.0 8.0 8.5, 9.0 10.0 10.5, 11.0 12.0 12.5, 7.0 8.0 8.5)))",
            format!("{multi_polygon:?}")
        );
    }

    #[test]
    fn multi_polygon_with_holes() {
        let multi_polygon = wkt!(MULTIPOLYGON(
            (
                (1.0 1.0 1.0, 10.0 1.0 10.0, 10.0 10.0 10.0, 1.0 10.0 1.0, 1.0 1.0 1.0)
            ),
            (
                (20. 20. 20., 30. 20. 30., 30. 30. 30., 20. 30. 20., 20. 20. 20.),
                (22. 22. 22., 28. 22. 28., 28. 28. 28., 22. 28. 22., 22. 22. 22.)
            )
        ));
        assert_eq!(
            "MULTIPOLYGON(((1.0 1.0 1.0, 10.0 1.0 10.0, 10.0 10.0 10.0, 1.0 10.0 1.0, 1.0 1.0 1.0)),((20.0 20.0 20.0, 30.0 20.0 30.0, 30.0 30.0 30.0, 20.0 30.0 20.0, 20.0 20.0 20.0),(22.0 22.0 22.0, 28.0 22.0 28.0, 28.0 28.0 28.0, 22.0 28.0 22.0, 22.0 22.0 22.0)))",
            format!("{multi_polygon:?}")
        );
    }
    #[test]
    fn multi_polygon_with_holes_and_empty_polygon() {
        let multi_polygon = wkt!(MULTIPOLYGON(
            (
                (1.0 1.0 1.0, 10.0 1.0 10.0, 10.0 10.0 10.0, 1.0 10.0 1.0, 1.0 1.0 1.0)
            ),
            EMPTY,
            (
                (20.0 20.0 20.0, 30.0 20.0 30.0, 30.0 30.0 30.0, 20.0 30.0 20.0, 20.0 20.0 20.0),
                (22.0 22.0 22.0, 28.0 22.0 28.0, 28.0 28.0 28.0, 22.0 28.0 22.0, 22.0 22.0 22.0)
            )
        ));
        assert_eq!(
            "MULTIPOLYGON(((1.0 1.0 1.0, 10.0 1.0 10.0, 10.0 10.0 10.0, 1.0 10.0 1.0, 1.0 1.0 1.0)),EMPTY,((20.0 20.0 20.0, 30.0 20.0 30.0, 30.0 30.0 30.0, 20.0 30.0 20.0, 20.0 20.0 20.0),(22.0 22.0 22.0, 28.0 22.0 28.0, 28.0 28.0 28.0, 22.0 28.0 22.0, 22.0 22.0 22.0)))",
            format!("{multi_polygon:?}")
        );
    }
    #[test]
    fn rect() {
        let rect = Rect::new((1., 2., 1.), (3., 4., 3.));
        assert_eq!("RECT(1.0 2.0 1.0,3.0 4.0 3.0)", format!("{rect:?}"));

        let rect = Rect::new((3.0, 4.0, 3.0), (1.0, 2.0, 3.0));
        // output is always (min, max)
        assert_eq!("RECT(1.0 2.0 3.0,3.0 4.0 3.0)", format!("{rect:?}"));
    }
    #[test]
    fn triangle() {
        let rect = Triangle::new((1.0, 2.0, 3.0).into(), (3.0, 4.0, 5.0).into(), (5.0, 6.0, 7.0).into());
        assert_eq!("TRIANGLE(1.0 2.0 3.0,3.0 4.0 5.0,5.0 6.0 7.0)", format!("{rect:?}"));
    }

    #[test]
    fn geometry() {
        let rect = Geometry::Triangle(Triangle::new((1.0, 2.0, 3.0).into(), (3.0, 4.0, 5.0).into(), (5.0, 6.0, 7.0).into()));
        assert_eq!("TRIANGLE(1.0 2.0 3.0,3.0 4.0 5.0,5.0 6.0 7.0)", format!("{rect:?}"));
    }

    #[test]
    fn geometry_collection() {
        let rect = Geometry::Triangle(Triangle::new((1.0, 2.0, 3.0).into(), (3.0, 4.0, 5.0).into(), (5.0, 6.0, 7.0).into()));
        assert_eq!("TRIANGLE(1.0 2.0 3.0,3.0 4.0 5.0,5.0 6.0 7.0)", format!("{rect:?}"));
    }

    #[test]
    fn empty_geometry_collection() {
        let geometry_collection: GeometryCollection = GeometryCollection::default();
        assert_eq!(
            "GEOMETRYCOLLECTION EMPTY",
            format!("{geometry_collection:?}")
        );
    }

    #[test]
    fn geometry_collection_with_mixed_geometries() {
        let geometry_collection: GeometryCollection<f32> = GeometryCollection::from(vec![
            Geometry::Point(Point::new(1.0, 2.0, 3.0)),
            Geometry::Line(Line::new((1.0, 2.0, 3.0), (3.0, 4.0, 5.0))),
            Geometry::Polygon(Polygon::new(
                LineString::from(vec![(0.0, 0.0, 0.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 0.0, 0.0)]),
                vec![],
            )),
        ]);

        assert_eq!(
            "GEOMETRYCOLLECTION(POINT(1.0 2.0 3.0),LINE(1.0 2.0 3.0,3.0 4.0 5.0),POLYGON((0.0 0.0 0.0,1.0 0.0 1.0,1.0 1.0 1.0,0.0 0.0 0.0)))",
            format!("{geometry_collection:?}")
        );
    }

    #[test]
    fn nested_geometry_collection() {
        let inner_collection: GeometryCollection<f32> = GeometryCollection::from(vec![
            Geometry::Point(Point::new(5.0, 6.0, 7.0)),
            Geometry::LineString(LineString::from(vec![(1.0, 2.0, 3.0), (3.0, 4.0, 5.0)])),
        ]);

        let outer_collection: GeometryCollection<f32> = GeometryCollection::from(vec![
            Geometry::Point(Point::new(1.0, 2.0, 3.0)),
            Geometry::GeometryCollection(inner_collection),
        ]);

        assert_eq!(
            "GEOMETRYCOLLECTION(POINT(1.0 2.0 3.0),GEOMETRYCOLLECTION(POINT(5.0 6.0 7.0),LINESTRING(1.0 2.0 3.0,3.0 4.0 5.0)))",
            format!("{outer_collection:?}")
        );
    }

    #[test]
    fn geometry_collection_with_no_coordinates() {
        let geometry_collection: GeometryCollection<f64> = GeometryCollection::from(vec![
            Geometry::Point(Point::new(0.0, 0.0, 0.0)),
            Geometry::Polygon(Polygon::new(LineString::new(vec![]), vec![])),
        ]);

        assert_eq!(
            "GEOMETRYCOLLECTION(POINT(0.0 0.0 0.0),POLYGON EMPTY)",
            format!("{geometry_collection:?}")
        );
    }
}
