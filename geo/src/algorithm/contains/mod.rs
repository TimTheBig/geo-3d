/// Checks if `rhs` is completely contained within `self`.
/// More formally, the interior of `rhs` has non-empty
/// (set-theoretic) intersection but neither the interior,
/// nor the boundary of `rhs` intersects the exterior of
/// `self`.
///
/// # Examples
///
/// ```
/// use geo::Contains;
/// use geo::{line_string, point, Polygon};
///
/// let line_string = line_string![
///     (x: 0., y: 0., z: 0.),
///     (x: 2., y: 0., z: -2.),
///     (x: 2., y: 2., z: 2.),
///     (x: 0., y: 2., z: 0.),
///     (x: 0., y: 0., z: 0.),
/// ];
///
/// let polygon = Polygon::new(line_string.clone(), vec![]);
///
/// // Point in Point
/// assert!(point!(x: 2., y: 0., z: -2.).contains(&point!(x: 2., y: 0., z: -2.)));
///
/// // Point in Linestring
/// assert!(line_string.contains(&point!(x: 2., y: 0., z: -2.)));
///
/// // Point in Polygon
/// assert!(polygon.contains(&point!(x: 1., y: 1., z: 1.)));
/// ```
pub trait Contains<Rhs = Self> {
    fn contains(&self, rhs: &Rhs) -> bool;
}

/// Checks if `rhs` (x and y) is completely contained within `self`.
/// More formally, the interior of `rhs` has non-empty
/// (set-theoretic) intersection but neither the interior,
/// nor the boundary of `rhs` intersects the exterior of
/// `self`. In other words, the [DE-9IM] intersection matrix
/// of `(rhs, self)` is `T*F**F***`.
///
/// [DE-9IM]: https://en.wikipedia.org/wiki/DE-9IM
///
/// # Examples
///
/// ```
/// use geo::Contains2D;
/// use geo::{line_string, point, Polygon};
///
/// let line_string = line_string![
///     (x: 0., y: 0., z: 0.),
///     (x: 2., y: 0., z: -2.),
///     (x: 2., y: 2., z: 2.),
///     (x: 0., y: 2., z: 0.),
///     (x: 0., y: 0., z: -80000.),
/// ];
///
/// let polygon = Polygon::new(line_string.clone(), vec![]);
///
/// // Point in Point
/// assert!(point!(x: 2., y: 0., z: -2.).contains_2d(&point!(x: 2., y: 0., z: -2.)));
///
/// // Point in Linestring
/// assert!(line_string.contains_2d(&point!(x: 2., y: 0., z: -2.)));
///
/// // Point in Polygon
/// assert!(polygon.contains_2d(&point!(x: 1., y: 1., z: 1.)));
/// ```
pub trait ContainsXY<Rhs = Self> {
    fn contains_2d(&self, rhs: &Rhs) -> bool;
}

mod geometry;
mod geometry_collection;
mod line;
mod line_string;
mod point;
mod polygon;
mod rect;
mod triangle;

macro_rules! impl_contains_from_relate {
    ($for:ty,  [$($target:ty),*]) => {
        $(
            impl<T> Contains<$target> for $for
            where
                T: GeoNum
            {
                fn contains(&self, target: &$target) -> bool {
                    use $crate::algorithm::Relate;
                    self.relate(target).is_contains()
                }
            }
        )*
    };
}
pub(crate) use impl_contains_from_relate;

macro_rules! impl_contains_geometry_for {
    ($geom_type: ty) => {
        impl<T> Contains<Geometry<T>> for $geom_type
        where
            T: GeoNum + std::iter::Sum,
        {
            fn contains(&self, geometry: &Geometry<T>) -> bool {
                match geometry {
                    Geometry::Point(g) => self.contains(g),
                    Geometry::Line(g) => self.contains(g),
                    Geometry::LineString(g) => self.contains(g),
                    Geometry::Polygon(g) => self.contains(g),
                    Geometry::MultiPoint(g) => self.contains(g),
                    Geometry::MultiLineString(g) => self.contains(g),
                    Geometry::MultiPolygon(g) => self.contains(g),
                    Geometry::GeometryCollection(g) => self.contains(g),
                    Geometry::Rect(g) => self.contains(g),
                    Geometry::Triangle(g) => self.contains(g),
                }
            }
        }
    };
}
pub(crate) use impl_contains_geometry_for;

// ┌───────┐
// │ Tests │
// └───────┘

#[cfg(test)]
mod test {
    use crate::line_string;
    use crate::Contains;
    use crate::Relate;
    use crate::{coord, Coord, Line, LineString, MultiPolygon, Point, Polygon, Rect, Triangle};

    #[test]
    // see https://github.com/georust/geo/issues/452
    fn linestring_contains_point() {
        let line_string = LineString::from(vec![(0., 0., 0.), (3., 3., 3.)]);
        let point_on_line = Point::new(1., 1., 1.);
        assert!(line_string.contains(&point_on_line));
    }

    #[test]
    // V doesn't contain rect because two of its edges intersect with V's exterior boundary
    fn polygon_does_not_contain_polygon() {
        let v = Polygon::new(
            vec![
                (150., 350., 150.),
                (100., 350., 100.),
                (210., 160., 210.),
                (290., 350., 290.),
                (250., 350., 250.),
                (200., 250., 200.),
                (150., 350., 150.),
            ]
            .into(),
            vec![],
        );
        let rect = Polygon::new(
            vec![
                (250., 310., 250.),
                (150., 310., 150.),
                (150., 280., 150.),
                (250., 280., 250.),
                (250., 310., 250.),
            ]
            .into(),
            vec![],
        );
        assert!(!v.contains(&rect));
    }
    #[test]
    // V contains rect because all its vertices are contained, and none of its edges intersect with V's boundaries
    fn polygon_contains_polygon() {
        let v = Polygon::new(
            vec![
                (150., 350., 150.),
                (100., 350., 100.),
                (210., 160., 210.),
                (290., 350., 290.),
                (250., 350., 250.),
                (200., 250., 200.),
                (150., 350., 150.),
            ]
            .into(),
            vec![],
        );
        let rect = Polygon::new(
            vec![
                (185., 237., 185.),
                (220., 237., 220.),
                (220., 220., 220.),
                (185., 220., 185.),
                (185., 237., 185.),
            ]
            .into(),
            vec![],
        );
        assert!(v.contains(&rect));
    }
    #[test]
    // LineString is fully contained
    fn linestring_fully_contained_in_polygon() {
        let poly = Polygon::new(
            LineString::from(vec![(0., 0., 0.), (5., 0., 5.), (5., 6., 5.), (0., 6., 0.), (0., 0., 0.)]),
            vec![],
        );
        let ls = LineString::from(vec![(3.0, 0.5, 3.0), (3.0, 3.5, 3.0)]);
        assert!(poly.contains(&ls));
    }
    /// Tests: Point in LineString
    #[test]
    fn empty_linestring_test() {
        let linestring = LineString::new(Vec::new());
        assert!(!linestring.contains(&Point::new(2., 1., 3.)));
    }
    #[test]
    fn linestring_point_is_vertex_test() {
        let linestring = LineString::from(vec![(0., 0., 0.), (2., 0., 2.), (2., 2., 2.)]);
        // Note: the end points of a linestring are not
        // considered to be "contained"
        assert!(linestring.contains(&Point::new(2., 0., 2.)));
        assert!(!linestring.contains(&Point::new(0., 0., 0.)));
        assert!(!linestring.contains(&Point::new(2., 2., 2.)));
    }
    #[test]
    fn linestring_test() {
        let linestring = LineString::from(vec![(0., 0., 1.), (2., 0., 2.), (2., 2., 1.)]);
        assert!(linestring.contains(&Point::new(1., 0., 1.5)));
    }
    /// Tests: Point in Polygon
    #[test]
    fn empty_polygon_test() {
        let linestring = LineString::new(Vec::new());
        let poly = Polygon::new(linestring, Vec::new());
        assert!(!poly.contains(&Point::new(2., 1., 3.)));
    }
    #[test]
    fn polygon_with_one_point_test() {
        let linestring = LineString::from(vec![(2., 1., 0.4)]);
        let poly = Polygon::new(linestring, Vec::new());
        assert!(!poly.contains(&Point::new(3., 1., 0.4)));
    }
    #[test]
    fn polygon_with_one_point_is_vertex_test() {
        let linestring = LineString::from(vec![(2., 1., 2.2)]);
        let poly = Polygon::new(linestring, Vec::new());
        assert!(!poly.contains(&Point::new(2., 1., 2.2)));
    }
    #[test]
    fn polygon_with_point_on_boundary_test() {
        let linestring = LineString::from(vec![(0., 0., 0.), (2., 0., 2.), (2., 2., 2.), (0., 2., 0.), (0., 0., 0.)]);
        let poly = Polygon::new(linestring, Vec::new());
        assert!(!poly.contains(&Point::new(1., 0., 1.)));
        assert!(!poly.contains(&Point::new(2., 1., 2.)));
        assert!(!poly.contains(&Point::new(1., 2., 1.)));
        assert!(!poly.contains(&Point::new(0., 1., 0.)));
    }
    #[test]
    fn point_in_polygon_test() {
        let linestring = LineString::from(vec![(0., 0., 1.), (2., 0., 0.5), (2., 2., 1.), (0., 2., 1.), (0., 0., 1.)]);
        let poly = Polygon::new(linestring, Vec::new());
        assert!(poly.contains(&Point::new(1., 1., 1.)));
    }
    #[test]
    fn point_in_polygon_with_ray_passing_through_a_vertex_test() {
        let linestring = LineString::from(vec![(1., 0., 1.), (0., 1., 0.), (-1., 0., -1.), (0., -1., 0.)]);
        let poly = Polygon::new(linestring, Vec::new());
        assert!(poly.contains(&Point::new(0., 0., 0.)));
    }
    #[test]
    fn point_in_polygon_with_ray_passing_through_a_vertex_and_not_crossing() {
        let linestring = LineString::from(vec![
            (0., 0., 0.),
            (2., 0., 2.),
            (3., 1., 3.),
            (4., 0., 4.),
            (4., 2., 4.),
            (0., 2., 0.),
            (0., 0., 0.),
        ]);
        let poly = Polygon::new(linestring, Vec::new());
        assert!(poly.contains(&Point::new(1., 1., 1.)));
    }
    #[test]
    fn point_out_polygon_test() {
        let linestring = LineString::from(vec![(0., 0., 0.), (2., 0., 2.), (2., 2., 2.), (0., 2., 0.), (0., 0., 0.)]);
        let poly = Polygon::new(linestring, Vec::new());
        assert!(!poly.contains(&Point::new(2.1, 1., 2.)));
        assert!(!poly.contains(&Point::new(1., 2.1, 2.)));
        assert!(!poly.contains(&Point::new(2.1, 2.1, 2.)));
        assert!(!poly.contains(&Point::new(2., 2., 2.1)));
    }
    #[test]
    fn point_polygon_with_inner_test() {
        let linestring = LineString::from(vec![(0., 0., 0.), (2., 0., 2.), (2., 2., 2.), (0., 2., 0.), (0., 0., 0.)]);
        let inner_linestring = LineString::from(vec![
            [0.5, 0.5, 0.5],
            [1.5, 0.5, 1.5],
            [1.5, 1.5, 1.5],
            [0.0, 1.5, 0.0],
            [0.0, 0.0, 0.0],
        ]);
        let poly = Polygon::new(linestring, vec![inner_linestring]);
        assert!(!poly.contains(&Point::new(0.25, 0.25, 0.25)));
        assert!(!poly.contains(&Point::new(1., 1., 1.)));
        assert!(!poly.contains(&Point::new(1.5, 1.5, 1.5)));
        assert!(!poly.contains(&Point::new(1.5, 1., 1.5)));
    }

    /// Tests: Point in MultiPolygon
    #[test]
    fn empty_multipolygon_test() {
        let multipoly = MultiPolygon::new(Vec::new());
        assert!(!multipoly.contains(&Point::new(2., 1., 0.)));
    }
    #[test]
    fn empty_multipolygon_two_polygons_test() {
        let poly1 = Polygon::new(
            LineString::from(vec![(0., 0., 0.), (1., 0., 1.), (1., 1., 1.), (0., 1., 0.), (0., 0., 0.)]),
            Vec::new(),
        );
        let poly2 = Polygon::new(
            LineString::from(vec![(2., 0., 2.), (3., 0., 3.), (3., 1., 3.), (2., 1., 2.), (2., 0., 2.)]),
            Vec::new(),
        );
        let multipoly = MultiPolygon::new(vec![poly1, poly2]);
        assert!(multipoly.contains(&Point::new(0.5, 0.5, 0.5)));
        assert!(multipoly.contains(&Point::new(2.5, 0.5, 2.5)));
        assert!(!multipoly.contains(&Point::new(1.5, 0.5, 1.5)));
    }
    #[test]
    fn empty_multipolygon_two_polygons_and_inner_test() {
        let poly1 = Polygon::new(
            LineString::from(vec![(0., 0., 0.), (5., 0., 5.), (5., 6., 5.), (0., 6., 0.), (0., 0., 0.)]),
            vec![LineString::from(vec![
                (1., 1., 1.),
                (4., 1., 4.),
                (4., 4., 4.),
                (1., 1., 1.),
            ])],
        );
        let poly2 = Polygon::new(
            LineString::from(vec![(9., 0., 9.), (14., 0., 14.), (14., 4., 14.), (9., 4., 9.), (9., 0., 9.)]),
            Vec::new(),
        );

        let multipoly = MultiPolygon::new(vec![poly1, poly2]);
        assert!(multipoly.contains(&Point::new(3., 5., 10.)));
        assert!(multipoly.contains(&Point::new(12., 2., 12.)));
        assert!(!multipoly.contains(&Point::new(3., 2., 9.)));
        assert!(!multipoly.contains(&Point::new(7., 2., 11.)));
        assert!(!multipoly.contains(&Point::new(3., 5., 7.)));
    }
    /// Tests: LineString in Polygon
    #[test]
    fn linestring_in_polygon_with_linestring_is_boundary_test() {
        let linestring = LineString::from(vec![(0., 0., 0.), (2., 0., 2.), (2., 2., 2.), (0., 2., 0.), (0., 0., 0.)]);
        let poly = Polygon::new(linestring.clone(), Vec::new());
        assert!(!poly.contains(&linestring));
        assert!(!poly.contains(&LineString::from(vec![(0., 0., 0.), (2., 0., 2.)])));
        assert!(!poly.contains(&LineString::from(vec![(2., 0., 2.), (2., 2., 2.)])));
        assert!(!poly.contains(&LineString::from(vec![(0., 2., 0.), (0., 0., 0.)])));
        assert!(!poly.contains(&LineString::from(vec![(2., 2., 3.), (0., 0., 0.)])));
    }
    #[test]
    fn linestring_outside_polygon_test() {
        let linestring = LineString::from(vec![(0., 0., 0.), (2., 0., 2.), (2., 2., 2.), (0., 2., 0.), (0., 0., 0.)]);
        let poly = Polygon::new(linestring, Vec::new());
        assert!(!poly.contains(&LineString::from(vec![(1., 1., 1.), (3., 0., 2.)])));
        assert!(!poly.contains(&LineString::from(vec![(1., 1., 1.), (2., 0., 3.)])));
        assert!(!poly.contains(&LineString::from(vec![(3., 0., 3.), (5., 2., 2.)])));
        assert!(!poly.contains(&LineString::from(vec![(3., 0., 3.), (2., 2., 5.)])));
    }
    #[test]
    fn linestring_in_inner_polygon_test() {
        let poly = Polygon::new(
            LineString::from(vec![(0., 0., 0.), (5., 0., 5.), (5., 6., 5.), (0., 6., 0.), (0., 0., 0.)]),
            vec![LineString::from(vec![
                (1., 1., 1.),
                (4., 1., 4.),
                (4., 4., 4.),
                (1., 4., 1.),
                (1., 1. ,1.),
            ])],
        );
        assert!(!poly.contains(&LineString::from(vec![(2., 2., 5.), (3., 3., 3.)])));
        assert!(!poly.contains(&LineString::from(vec![(2., 2., 2.), (2., 5., 2.)])));
        assert!(!poly.contains(&LineString::from(vec![(3., 0.5, 3.), (3., 5., 5.)])));
        assert!(!poly.contains(&LineString::from(vec![(4.5, 1.5, 6.), (3., 5., 7.)])));
    }
    #[test]
    fn bounding_rect_in_inner_bounding_rect_test() {
        let bounding_rect_xl =
            Rect::new(coord! { x: -100., y: -200., z: -300. }, coord! { x: 100., y: 200., z: 300. });
        let bounding_rect_sm = Rect::new(coord! { x: -10., y: -20., z: -30. }, coord! { x: 10., y: 20., z: 30. });
        assert!(bounding_rect_xl.contains(&bounding_rect_sm));
        assert!(!bounding_rect_sm.contains(&bounding_rect_xl));
    }
    #[test]
    fn point_in_line_test() {
        let c = |x, y, z| coord! { x: x, y: y, z: z };
        let p0 = c(2., 4., 6.);
        // vertical line
        let line1 = Line::new(c(2., 0., 2.), c(2., 5., 2.));
        // point on line, but outside line segment
        let line2 = Line::new(c(0., 6., 0.), c(1.5, 4.5, 1.5));
        // point on line
        let line3 = Line::new(c(0., 6., 0.), c(3., 3., 3.));
        assert!(line1.contains(&Point::from(p0)));
        assert!(!line2.contains(&Point::from(p0)));
        assert!(line3.contains(&Point::from(p0)));
    }
    #[test]
    fn line_in_line_test() {
        let c = |x, y, z| coord! { x: x, y: y, z: z };
        let line0 = Line::new(c(0., 1., 0.), c(3., 4., 3.));
        // first point on line0, second not
        let line1 = Line::new(c(1., 2., 1.), c(2., 2., 2.));
        // co-linear, but extends past the end of line0
        let line2 = Line::new(c(1., 2., 1.), c(4., 5., 4.));
        // contained in line0
        let line3 = Line::new(c(1., 2., 1.), c(3., 4., 3.));
        assert!(!line0.contains(&line1));
        assert!(!line0.contains(&line2));
        assert!(line0.contains(&line3));
    }
    #[test]
    fn linestring_in_line_test() {
        let line = Line::from([(0., 10., 0.), (30., 40., 30.)]);
        // linestring0 in line
        let linestring0 = LineString::from(vec![(1., 11., 1.), (10., 20., 10.), (15., 25., 15.)]);
        // linestring1 starts and ends in line, but wanders in the middle
        let linestring1 = LineString::from(vec![(1., 11., 1.), (20., 20., 20.), (15., 25., 15.)]);
        // linestring2 is co-linear, but extends beyond line
        let linestring2 = LineString::from(vec![(1., 11., 1.), (10., 20., 30.), (40., 50., 60.)]);
        // no part of linestring3 is contained in line
        let linestring3 = LineString::from(vec![(11., 11., 11.), (20., 20., 20.), (25., 25., 25.)]);
        // a linestring with singleton interior on the boundary of the line
        let linestring4 = LineString::from(vec![(0., 10., 0.), (0., 10., 20.), (0., 10., 0.)]);
        // a linestring with singleton interior that is contained in the line
        let linestring5 = LineString::from(vec![(1., 11., 111.), (1., 11., 1.), (1., 11., 111.)]);
        assert!(line.contains(&linestring0));
        assert!(!line.contains(&linestring1));
        assert!(!line.contains(&linestring2));
        assert!(!line.contains(&linestring3));
        assert!(!line.contains(&linestring4));
        assert!(line.contains(&linestring5));
    }
    #[test]
    fn line_in_polygon_test() {
        let line = Line::new(coord!(0.0, 10.0, 0.0), coord!(30.0, 40.0, 50.0));
        let linestring0 = line_string![
            coord!(-10.0, 0.0, -10.0),
            coord!(50.0, 0.0, 50.0),
            coord!(50.0, 50.0, 50.0),
            coord!(0.0, 50.0, 0.0),
            coord!(-10.0, 0.0, -10.0)
        ];
        let poly0 = Polygon::new(linestring0, Vec::new());
        let linestring1 = line_string![
            coord!(0.0, 0.0, 0.0),
            coord!(0.0, 20.0, 0.0),
            coord!(20.0, 20.0, 20.0),
            coord!(20.0, 0.0, 20.0),
            coord!(0.0, 0.0, 0.0)
        ];
        let poly1 = Polygon::new(linestring1, Vec::new());
        assert!(poly0.contains(&line));
        assert!(!poly1.contains(&line));
    }
    #[test]
    fn line_in_polygon_edgecases_test() {
        // Some DE-9IM edge cases for checking line is
        // inside polygon The end points of the line can be
        // on the boundary of the polygon.
        // A non-convex polygon
        let linestring0 = line_string![
            coord!(0.0, 0.0, 0.0),
            coord!(1.0, 1.0, 1.0),
            coord!(1.0, -1.0, 1.0),
            coord!(-1.0, -1.0, -1.0),
            coord!(-1.0, 1.0, -1.0)
        ];
        let poly = Polygon::new(linestring0, Vec::new());

        assert!(poly.contains(&Line::new(coord!(0.0, 0.0, 0.0), coord!(1.0, -1.0, 1.0))));
        assert!(poly.contains(&Line::new(coord!(-1.0, 1.0, -1.0), coord!(1.0, -1.0, 1.0))));
        assert!(!poly.contains(&Line::new(coord!(-1.0, 1.0, -1.0), coord!(1.0, 1.0, 1.0))));
    }
    #[test]
    fn line_in_linestring_edgecases() {
        use crate::line_string;
        let mut ls = line_string![coord!(0., 0., 0.), coord!(1., 0., 1.), coord!(0., 1., 0.), coord!(-1., 0., -1.)];
        assert!(!ls.contains(&Line::from([(0., 0., 0.), (0., 0., 0.)])));
        ls.close();
        assert!(ls.contains(&Line::from([(0., 0., 0.), (0., 0., 0.)])));
        assert!(ls.contains(&Line::from([(-1., 0., -1.), (1., 0., 1.)])));
    }
    #[test]
    fn line_in_linestring_test() {
        let line0 = Line::from([(1., 1., 1.), (2., 2., 2.)]);
        // line0 is completely contained in the second segment
        let linestring0 = LineString::from(vec![(0., 0.5, 0.), (0.5, 0.5, 0.5), (3., 3., 3.)]);
        // line0 is contained in the last three segments
        let linestring1 = LineString::from(vec![
            (0., 0.5, 0.),
            (0.5, 0.5, 0.5),
            (1.2, 1.2, 1.2),
            (1.5, 1.5, 1.5),
            (3., 3., 3.),
        ]);
        // line0 endpoints are contained in the linestring, but the fourth point is off the line
        let linestring2 = LineString::from(vec![
            (0., 0.5, 0.),
            (0.5, 0.5, 0.5),
            (1.2, 1.2, 1.2),
            (1.5, 0., 1.5),
            (2., 2., 2.),
            (3., 3., 3.),
        ]);
        assert!(linestring0.contains(&line0));
        assert!(linestring1.contains(&line0));
        assert!(!linestring2.contains(&line0));
    }

    #[test]
    fn f32_bounding_rects() {
        let p: Point<f32> = Point::new(10., 20., 30.);
        let bounding_rect: Rect<f32> = Rect::new(Coord::zero(), coord! { x: 100., y: 100., z: 100. });
        assert!(bounding_rect.contains(&p));
        assert!(!bounding_rect.contains(&Point::new(-10., -10., -10.)));

        let smaller_bounding_rect: Rect<f32> =
            Rect::new(coord! { x: 10., y: 10., z: 10. }, coord! { x: 20., y: 20., z: 20. });
        assert!(bounding_rect.contains(&smaller_bounding_rect));
    }

    #[test]
    fn triangle_not_contains_point_on_edge() {
        let t = Triangle::from([(0.0, 0.0, 0.0), (2.0, 0.0, 2.0), (2.0, 2.0, 2.0)]);
        let p = Point::new(1.0, 0.0, -1.0);
        assert!(!t.contains(&p));
    }

    #[test]
    fn triangle_not_contains_point_on_vertex() {
        let t = Triangle::from([(0.0, 0.0, 0.0), (2.0, 0.0, -2.0), (2.0, 2.0, -2.0)]);
        let p = Point::new(2.0, 0.0, -2.0);
        assert!(!t.contains(&p));
    }

    #[test]
    fn triangle_contains_point_inside() {
        let t = Triangle::from([(0.0, 0.0, 0.0), (2.0, 0.0, 2.0), (2.0, 2.0, 2.0)]);
        let p = Point::new(1.0, 0.5, 1.0);
        assert!(t.contains(&p));
    }

    #[test]
    fn triangle_not_contains_point_above() {
        let t = Triangle::from([(0.0, 0.0, 0.0), (2.0, 0.0, 2.0), (2.0, 2.0, 2.0)]);
        let p = Point::new(1.0, 1.5, 2.0);
        assert!(!t.contains(&p));
    }

    #[test]
    fn triangle_not_contains_point_below() {
        let t = Triangle::from([(0.0, 0.0, 0.0), (2.0, 0.0, -2.0), (2.0, 2.0, -2.0)]);
        let p = Point::new(-1.0, 0.5, 1.0);
        assert!(!t.contains(&p));
    }

    #[test]
    fn triangle_contains_neg_point() {
        let t = Triangle::from([(0.0, 0.0, 0.0), (-2.0, 0.0, 2.0), (-2.0, -2.0, 2.0)]);
        let p = Point::new(-1.0, -0.5, 1.5);
        assert!(t.contains(&p));
    }

    #[test]
    // https://github.com/georust/geo/issues/473
    fn triangle_contains_collinear_points() {
        let origin: Coord = (0., 0., 0.).into();
        let tri = Triangle::new(origin, origin, origin);
        let pt: Point = (0., 1.23456, 1.8938202).into();
        assert!(!tri.contains(&pt));
        let pt: Point = (0., 0., 0.).into();
        assert!(!tri.contains(&pt));
        let origin: Coord = (0., 0., 0.).into();
        let tri = Triangle::new((1., 1., 1.).into(), origin, origin);
        let pt: Point = (1., 1., 1.).into();
        assert!(!tri.contains(&pt));
        let pt: Point = (0.5, 0.5, 0.5).into();
        assert!(!tri.contains(&pt));
    }

    #[test]
    fn rect_contains_polygon() {
        let rect = Rect::new(coord! { x: 90., y: 150., z: 190. }, coord! { x: 300., y: 360., z: 390. });
        let poly = Polygon::new(
            line_string![
                (x: 150., y: 350., z: 190.),
                (x: 100., y: 350., z: 190.),
                (x: 150., y: 350., z: 220.),
                (x: 210., y: 160., z: 235.),
                (x: 290., y: 350., z: 195.),
                (x: 250., y: 350., z: 300.),
                (x: 200., y: 250., z: 370.),
                (x: 150., y: 350., z: 190.),
            ],
            vec![],
        );
        assert_eq!(rect.contains(&poly), rect.relate(&poly).is_contains());
    }

    #[test]
    fn rect_contains_touching_polygon() {
        let rect = Rect::new(coord! { x: 90., y: 150., z: 180. }, coord! { x: 300., y: 360., z: 420. });
        let touching_poly = Polygon::new(
            line_string![
                (x: 150., y: 350., z: 340.),
                (x: 90.,  y: 350., z: 360.),
                (x: 210., y: 160., z: 240.),
                (x: 290., y: 350., z: 190.),
                (x: 250., y: 350., z: 400.),
                (x: 200., y: 250., z: 360.),
                (x: 150., y: 350., z: 340.),
            ],
            vec![],
        );
        assert_eq!(
            rect.contains(&touching_poly),
            rect.relate(&touching_poly).is_contains()
        );

        let touching_rect = Rect::new(coord! { x: 90., y: 200., z: 380. }, coord! { x: 200., y: 300., z: 400. });
        assert_eq!(
            rect.contains(&touching_rect),
            rect.relate(&touching_rect).is_contains()
        );
    }

    #[test]
    fn rect_contains_empty_polygon() {
        let rect = Rect::new(coord! { x: 90., y: 150., z: 270. }, coord! { x: 300., y: 360., z: 390. });
        let empty_poly = Polygon::new(line_string![], vec![]);
        assert_eq!(
            rect.contains(&empty_poly),
            rect.relate(&empty_poly).is_contains()
        );
    }

    #[test]
    fn rect_contains_polygon_empty_area() {
        let rect = Rect::new(coord! { x: 90., y: 150., z: 270. }, coord! { x: 300., y: 360., z: 420. });
        let empty_poly = Polygon::new(
            line_string![
                (x: 100., y: 200., z: 300.),
                (x: 100., y: 200., z: 300.),
                (x: 100., y: 200., z: 300.),
                (x: 100., y: 200., z: 300.),
            ],
            vec![],
        );
        assert_eq!(
            rect.contains(&empty_poly),
            rect.relate(&empty_poly).is_contains()
        );
    }

    #[test]
    fn rect_contains_rect_polygon() {
        let rect = Rect::new(coord! { x: 90., y: 150., z: 210. }, coord! { x: 300., y: 360., z: 410. });
        let rect_poly = rect.to_polygon();
        assert_eq!(
            rect.contains(&rect_poly),
            rect.relate(&rect_poly).is_contains()
        );
    }

    #[test]
    fn rect_contains_polygon_in_boundary() {
        let rect = Rect::new(coord! { x: 90. , y: 150., z: 90. }, coord! { x: 300., y: 360., z: 300. });
        let poly_one_border =
            Rect::new(coord! { x: 90. , y: 150., z: 90. }, coord! { x: 90., y: 360., z: 90. }).to_polygon();
        assert_eq!(
            rect.contains(&poly_one_border),
            rect.relate(&poly_one_border).is_contains()
        );

        let poly_two_borders = Polygon::new(
            line_string![
                (x: 90., y: 150., z: 6.),
                (x: 300., y: 150., z: 6.),
                (x: 90., y: 150., z: 6.),
                (x: 90., y: 360., z: 12.),
                (x: 90., y: 150., z: 6.),
            ],
            vec![],
        );
        assert_eq!(
            rect.contains(&poly_two_borders),
            rect.relate(&poly_two_borders).is_contains()
        );

        let poly_two_borders_triangle = Polygon::new(
            line_string![
                (x: 90., y: 150., z: 0.),
                (x: 300., y: 150., z: 0.),
                (x: 90., y: 360., z: 0.),
                (x: 90., y: 150., z: 0.),
            ],
            vec![],
        );
        assert_eq!(
            rect.contains(&poly_two_borders_triangle),
            rect.relate(&poly_two_borders_triangle).is_contains()
        );
    }

    #[test]
    fn rect_contains_polygon_in_boundary_with_hole() {
        let rect = Rect::new(coord! { x: 90. , y: 150., z: 90. }, coord! { x: 300., y: 360., z: 300. });
        let poly_two_borders_triangle_with_hole = Polygon::new(
            line_string![
                (x: 90., y: 150., z: 6.),
                (x: 300., y: 150., z: 6.),
                (x: 90., y: 360., z: 12.),
                (x: 90., y: 150., z: 6.),
            ],
            vec![line_string![
                (x: 90., y: 150., z: 6.),
                (x: 300., y: 150., z: 6.),
                (x: 90., y: 360., z: 12.),
                (x: 90., y: 150., z: 6.),
            ]],
        );
        assert_eq!(
            rect.contains(&poly_two_borders_triangle_with_hole),
            rect.relate(&poly_two_borders_triangle_with_hole)
                .is_contains()
        );
    }

    #[test]
    fn rect_empty_contains_polygon() {
        let rect = Rect::new(coord! { x: 90. , y: 150., z: 180. }, coord! { x: 90., y: 150., z: 180. });
        let poly = Polygon::new(
            line_string![
                (x: 150., y: 350., z: 550.),
                (x: 100., y: 350., z: 500.),
                (x: 210., y: 160., z: 360.),
                (x: 290., y: 350., z: 176.),
                (x: 250., y: 350., z: 450.),
                (x: 200., y: 250., z: 300.),
                (x: 150., y: 350., z: 550.),
            ],
            vec![],
        );
        assert_eq!(rect.contains(&poly), rect.relate(&poly).is_contains());

        let rect_poly = rect.to_polygon();
        assert_eq!(
            rect.contains(&rect_poly),
            rect.relate(&rect_poly).is_contains()
        );
    }

    #[test]
    fn rect_contains_point() {
        let rect = Rect::new(coord! { x: 90., y: 150., z: 230. }, coord! { x: 300., y: 360., z: 480. });

        let point1 = Point::new(100., 200., 300.);
        assert_eq!(rect.contains(&point1), rect.relate(&point1).is_contains());

        let point2 = Point::new(90., 200., 390.);
        assert_eq!(rect.contains(&point2), rect.relate(&point2).is_contains());
    }
}
