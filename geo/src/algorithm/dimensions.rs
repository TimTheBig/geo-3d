use crate::{geometry::*, Vector3DOps};
use crate::Orientation::Collinear;
use crate::{CoordNum, GeoNum, GeometryCow};

/// Geometries can have 0, 1, or two dimensions. Or, in the case of an [`empty`](#is_empty)
/// geometry, a special `Empty` dimensionality.
///
/// # Examples
///
/// ```
/// use geo_types::{Point, Rect, line_string};
/// use geo_3d::dimensions::{HasDimensions, Dimensions};
///
/// let point = Point::new(0.0, 5.0, 0.0);
/// let line_string = line_string![(x: 0.0, y: 0.0, z: 0.0), (x: 5.0, y: 5.0, z: 5.0), (x: 0.0, y: 5.0, z: 0.0)];
/// let rect = Rect::new((0.0, 0.0, 0.0), (10.0, 10.0, 10.0));
/// assert_eq!(Dimensions::ZeroDimensional, point.dimensions());
/// assert_eq!(Dimensions::OneDimensional, line_string.dimensions());
/// assert_eq!(Dimensions::TwoDimensional, rect.dimensions());
///
/// assert!(point.dimensions() < line_string.dimensions());
/// assert!(rect.dimensions() > line_string.dimensions());
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Ord, PartialOrd)]
pub enum Dimensions {
    /// Some geometries, like a `MultiPoint` or `GeometryCollection` may have no elements - thus no
    /// dimensions. Note that this is distinct from being `ZeroDimensional`, like a `Point`.
    Empty,
    /// Dimension of a point
    ZeroDimensional,
    /// Dimension of a line or curve
    OneDimensional,
    /// Dimension of a surface
    TwoDimensional,
    /// Dimension of a geometry
    ThreeDimensional,
}

/// Operate on the dimensionality of geometries.
pub trait HasDimensions {
    /// Some geometries, like a `MultiPoint`, can have zero coordinates - we call these `empty`.
    ///
    /// Types like `Point` and `Rect`, which have at least one coordinate by construction, can
    /// never be considered empty.
    /// ```
    /// use geo_types::{Point, coord, LineString};
    /// use geo_3d::HasDimensions;
    ///
    /// let line_string = LineString::new(vec![
    ///     coord! { x: 0., y: 0. },
    ///     coord! { x: 10., y: 0. },
    /// ]);
    /// assert!(!line_string.is_empty());
    ///
    /// let empty_line_string: LineString = LineString::new(vec![]);
    /// assert!(empty_line_string.is_empty());
    ///
    /// let point = Point::new(0.0, 0.0, 0.0);
    /// assert!(!point.is_empty());
    /// ```
    fn is_empty(&self) -> bool;

    /// The dimensions of some geometries are fixed, e.g. a Point always has 0 dimensions. However
    /// for others, the dimensionality depends on the specific geometry instance - for example
    /// typical `Rect`s are 3-dimensional, but it's possible to create degenerate `Rect`s which
    /// have either 2, 1 or 0 dimensions.
    ///
    /// ## Examples
    ///
    /// ```
    /// use geo_types::{GeometryCollection, Rect, Point};
    /// use geo_3d::dimensions::{Dimensions, HasDimensions};
    ///
    /// // normal rectangle
    /// let rect = Rect::new((0.0, 0.0, 0.0), (10.0, 10.0, 10.0));
    /// assert_eq!(Dimensions::TwoDimensional, rect.dimensions());
    ///
    /// // "rectangle" with zero height degenerates to a line
    /// let degenerate_line_rect = Rect::new((0.0, 10.0, 0.0), (10.0, 10.0, 10.0));
    /// assert_eq!(Dimensions::OneDimensional, degenerate_line_rect.dimensions());
    ///
    /// // "rectangle" with zero height and zero width degenerates to a point
    /// let degenerate_point_rect = Rect::new((10.0, 10.0, 10.0), (10.0, 10.0, 10.0));
    /// assert_eq!(Dimensions::ZeroDimensional, degenerate_point_rect.dimensions());
    ///
    /// // collections inherit the greatest dimensionality of their elements
    /// let geometry_collection = GeometryCollection::new(vec![degenerate_line_rect.into(), degenerate_point_rect.into()]);
    /// assert_eq!(Dimensions::OneDimensional, geometry_collection.dimensions());
    ///
    /// let point = Point::new(10.0, 10.0, 10.0);
    /// assert_eq!(Dimensions::ZeroDimensional, point.dimensions());
    ///
    /// // An `Empty` dimensionality is distinct from, and less than, being 0-dimensional
    /// let empty_collection = GeometryCollection::<f32>::new(vec![]);
    /// assert_eq!(Dimensions::Empty, empty_collection.dimensions());
    /// assert!(empty_collection.dimensions() < point.dimensions());
    /// ```
    fn dimensions(&self) -> Dimensions;

    /// The dimensions of the `Geometry`'s boundary, as used by OGC-SFA.
    ///
    /// ## Examples
    ///
    /// ```
    /// # use geo_types::{GeometryCollection, Rect, Point};
    /// use geo_3d::dimensions::{Dimensions, HasDimensions};
    ///
    /// // a point has no boundary
    /// let point = Point::new(10.0, 10.0, 10.0);
    /// assert_eq!(Dimensions::Empty, point.boundary_dimensions());
    ///
    /// // a typical rectangle has a *line* (one dimensional) boundary
    /// let rect = Rect::new((0.0, 0.0, 0.0), (10.0, 10.0, 10.0));
    /// assert_eq!(Dimensions::OneDimensional, rect.boundary_dimensions());
    ///
    /// // a "rectangle" with zero height degenerates to a line, whose boundary is two points
    /// let degenerate_line_rect = Rect::new((0.0, 10.0, 0.0), (10.0, 10.0, 10.0));
    /// assert_eq!(Dimensions::ZeroDimensional, degenerate_line_rect.boundary_dimensions());
    ///
    /// // a "rectangle" with zero height and zero width degenerates to a point,
    /// // and points have no boundary
    /// let degenerate_point_rect = Rect::new((10.0, 10.0, 10.0), (10.0, 10.0, 10.0));
    /// assert_eq!(Dimensions::Empty, degenerate_point_rect.boundary_dimensions());
    ///
    /// // collections inherit the greatest dimensionality of their elements
    /// let geometry_collection = GeometryCollection::new(vec![degenerate_line_rect.into(), degenerate_point_rect.into()]);
    /// assert_eq!(Dimensions::ZeroDimensional, geometry_collection.boundary_dimensions());
    ///
    /// let geometry_collection = GeometryCollection::<f32>::new(vec![]);
    /// assert_eq!(Dimensions::Empty, geometry_collection.boundary_dimensions());
    /// ```
    fn boundary_dimensions(&self) -> Dimensions;
}

impl<C: GeoNum> HasDimensions for Geometry<C> {
    crate::geometry_delegate_impl! {
        fn is_empty(&self) -> bool;
        fn dimensions(&self) -> Dimensions;
        fn boundary_dimensions(&self) -> Dimensions;
    }
}

impl<C: GeoNum> HasDimensions for GeometryCow<'_, C> {
    crate::geometry_cow_delegate_impl! {
        fn is_empty(&self) -> bool;
        fn dimensions(&self) -> Dimensions;
        fn boundary_dimensions(&self) -> Dimensions;
    }
}

impl<C: CoordNum> HasDimensions for Point<C> {
    fn is_empty(&self) -> bool {
        false
    }

    fn dimensions(&self) -> Dimensions {
        Dimensions::ZeroDimensional
    }

    fn boundary_dimensions(&self) -> Dimensions {
        Dimensions::Empty
    }
}

impl<C: CoordNum> HasDimensions for Line<C> {
    fn is_empty(&self) -> bool {
        false
    }

    fn dimensions(&self) -> Dimensions {
        if self.start == self.end {
            // degenerate line is a point
            Dimensions::ZeroDimensional
        } else if self.start.x != self.end.x
            && self.start.y != self.end.y
            && self.start.z != self.end.z
        {
            Dimensions::ThreeDimensional
        } else if (self.start.x == self.end.x && self.start.y == self.end.y)
            || (self.start.y == self.end.y && self.start.z == self.end.z)
            || (self.start.z == self.end.z && self.start.x == self.end.x)
        {
            Dimensions::TwoDimensional
        } else {
            Dimensions::OneDimensional
        }
    }

    fn boundary_dimensions(&self) -> Dimensions {
        if self.start == self.end {
            // degenerate line is a point, which has no boundary
            Dimensions::Empty
        } else {
            Dimensions::ZeroDimensional
        }
    }
}

impl<C: CoordNum> HasDimensions for LineString<C> {
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn dimensions(&self) -> Dimensions {
        if self.0.is_empty() {
            return Dimensions::Empty;
        }

        let first = self.0[0];
        if self.0.iter().any(|&coord| first != coord) {
            Dimensions::OneDimensional
        } else {
            // all coords are the same - i.e. a point
            Dimensions::ZeroDimensional
        }
    }

    /// ```
    /// # use geo_types::line_string;
    /// use geo_3d::dimensions::{HasDimensions, Dimensions};
    ///
    /// let ls = line_string![(x: 0.,  y: 0., z: 0.), (x: 0., y: 1., z: 0.), (x: 1., y: 1., z: 1.)];
    /// assert_eq!(Dimensions::ZeroDimensional, ls.boundary_dimensions());
    ///
    /// let ls = line_string![(x: 0.,  y: 0., z:0.), (x: 0., y: 1., z: 0.), (x: 1., y: 1., z:1.), (x: 0., y: 0., z: 0.)];
    /// assert_eq!(Dimensions::Empty, ls.boundary_dimensions());
    ///```
    fn boundary_dimensions(&self) -> Dimensions {
        if self.is_closed() {
            return Dimensions::Empty;
        }

        match self.dimensions() {
            Dimensions::Empty | Dimensions::ZeroDimensional => Dimensions::Empty,
            Dimensions::OneDimensional => Dimensions::ZeroDimensional,
            Dimensions::TwoDimensional => unreachable!("line_string cannot be 2 dimensional"),
            Dimensions::ThreeDimensional => unreachable!("line_string cannot be 3 dimensional"),
        }
    }
}

impl<C: CoordNum> HasDimensions for Polygon<C> {
    fn is_empty(&self) -> bool {
        self.exterior().is_empty()
    }

    fn dimensions(&self) -> Dimensions {
        use crate::CoordsIter;

        let mut coords = self.exterior_coords_iter();

        let Some(first) = coords.next() else {
            // No coordinates - the polygon is empty
            return Dimensions::Empty;
        };

        let Some(second) = coords.find(|next| *next != first) else {
            // All coordinates in the polygon are the same point
            return Dimensions::ZeroDimensional;
        };

        let Some(_third) = coords.find(|next| *next != first && *next != second) else {
            // There are only two distinct coordinates in the Polygon - it's collapsed to a line
            return Dimensions::OneDimensional;
        };

        // Check if all points lie on the same plane or axis
        let has_varying_x = self.rings().flatten().any(|c| c.x != first.x);

        let has_varying_y = self.rings().flatten().any(|c| c.y != first.y);

        let has_varying_z = self.rings().flatten().any(|c| c.z != first.z);

        if has_varying_x && has_varying_y && has_varying_z {
            Dimensions::ThreeDimensional
        } else {
            Dimensions::TwoDimensional
        }
    }

    fn boundary_dimensions(&self) -> Dimensions {
        match self.dimensions() {
            Dimensions::Empty | Dimensions::ZeroDimensional => Dimensions::Empty,
            Dimensions::OneDimensional => Dimensions::ZeroDimensional,
            Dimensions::TwoDimensional => Dimensions::OneDimensional,
            Dimensions::ThreeDimensional => Dimensions::ThreeDimensional,
        }
    }
}

impl<C: CoordNum> HasDimensions for MultiPoint<C> {
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn dimensions(&self) -> Dimensions {
        if self.0.is_empty() {
            return Dimensions::Empty;
        }

        Dimensions::ZeroDimensional
    }

    fn boundary_dimensions(&self) -> Dimensions {
        Dimensions::Empty
    }
}

impl<C: CoordNum> HasDimensions for MultiLineString<C> {
    fn is_empty(&self) -> bool {
        self.iter().all(LineString::is_empty)
    }

    fn dimensions(&self) -> Dimensions {
        let mut max = Dimensions::Empty;
        for line in &self.0 {
            match line.dimensions() {
                Dimensions::Empty => {}
                Dimensions::ZeroDimensional => max = Dimensions::ZeroDimensional,
                Dimensions::OneDimensional => {
                    // return early since we know multi line string dimensionality cannot exceed
                    // 1-d
                    return Dimensions::OneDimensional;
                }
                Dimensions::TwoDimensional => unreachable!("MultiLineString cannot be 2d"),
                Dimensions::ThreeDimensional => unreachable!("MultiLineString cannot be 3d"),
            }
        }
        max
    }

    fn boundary_dimensions(&self) -> Dimensions {
        if self.is_closed() {
            return Dimensions::Empty;
        }

        match self.dimensions() {
            Dimensions::Empty | Dimensions::ZeroDimensional => Dimensions::Empty,
            Dimensions::OneDimensional => Dimensions::ZeroDimensional,
            Dimensions::TwoDimensional => unreachable!("line_string cannot be 2 dimensional"),
            Dimensions::ThreeDimensional => unreachable!("MultiLineString cannot be 3d"),
        }
    }
}

impl<C: CoordNum> HasDimensions for MultiPolygon<C> {
    fn is_empty(&self) -> bool {
        self.iter().all(Polygon::is_empty)
    }

    fn dimensions(&self) -> Dimensions {
        let mut max = Dimensions::Empty;
        for geom in self {
            let dimensions = geom.dimensions();
            if dimensions == Dimensions::ThreeDimensional {
                // short-circuit since we know none can be larger
                return Dimensions::ThreeDimensional;
            }
            max = max.max(dimensions)
        }

        max
    }

    fn boundary_dimensions(&self) -> Dimensions {
        match self.dimensions() {
            Dimensions::Empty | Dimensions::ZeroDimensional => Dimensions::Empty,
            Dimensions::OneDimensional => Dimensions::ZeroDimensional,
            Dimensions::TwoDimensional => Dimensions::OneDimensional,
            Dimensions::ThreeDimensional => Dimensions::ThreeDimensional,
        }
    }
}

impl<C: GeoNum> HasDimensions for GeometryCollection<C> {
    fn is_empty(&self) -> bool {
        if self.0.is_empty() {
            true
        } else {
            self.iter().all(Geometry::is_empty)
        }
    }

    fn dimensions(&self) -> Dimensions {
        let mut max = Dimensions::Empty;
        for geom in self {
            let dimensions = geom.dimensions();
            if dimensions == Dimensions::ThreeDimensional {
                // short-circuit since we know none can be larger
                return Dimensions::ThreeDimensional;
            }
            max = max.max(dimensions)
        }

        max
    }

    fn boundary_dimensions(&self) -> Dimensions {
        let mut max = Dimensions::Empty;
        for geom in self {
            let d = geom.boundary_dimensions();

            if d == Dimensions::OneDimensional {
                return Dimensions::OneDimensional;
            }

            max = max.max(d);
        }
        max
    }
}

impl<C: CoordNum> HasDimensions for Rect<C> {
    fn is_empty(&self) -> bool {
        false
    }

    fn dimensions(&self) -> Dimensions {
        if self.min() == self.max() {
            // degenerate rectangle is a point
            Dimensions::ZeroDimensional
        } else if (self.min().x == self.max().x && self.min().y == self.max().y)
            || (self.min().x == self.max().x && self.min().z == self.max().z)
            || (self.min().y == self.max().y && self.min().z == self.max().z)
        {
            // degenerate rectangle is a line
            Dimensions::OneDimensional
        } else if self.min().x == self.max().x
            || self.min().y == self.max().y
            || self.min().z == self.max().z
        {
            // only two dimensions are used
            Dimensions::TwoDimensional
        } else {
            Dimensions::ThreeDimensional
        }
    }

    fn boundary_dimensions(&self) -> Dimensions {
        match self.dimensions() {
            Dimensions::Empty => {
                unreachable!("even a degenerate rect should be at least 0-Dimensional")
            }
            Dimensions::ZeroDimensional => Dimensions::Empty,
            Dimensions::OneDimensional => Dimensions::ZeroDimensional,
            Dimensions::TwoDimensional => Dimensions::OneDimensional,
            Dimensions::ThreeDimensional => Dimensions::ThreeDimensional,
        }
    }
}

impl<C: GeoNum> HasDimensions for Triangle<C> {
    fn is_empty(&self) -> bool {
        false
    }

    fn dimensions(&self) -> Dimensions {
        use crate::Kernel;

        // Check if all points are identical
        if self.0 == self.1 && self.1 == self.2 {
            // Triangle collapses to a single point
            return Dimensions::ZeroDimensional;
        }

        let cross_product = (self.1 - self.0).cross(self.2 - self.0);

        // Check if the triangle is degenerate and lies on a single line
        if C::Ker::orient2d(self.0, self.1, self.2) == Collinear {
            // Check if the third point is collinear in 3D space (all points lie on the same line)
            if cross_product.magnitude_squared() == C::zero() {
                return Dimensions::OneDimensional; // Degenerate triangle is a line
            }
        }

        // All points lie on the same plane
        Dimensions::TwoDimensional
    }

    fn boundary_dimensions(&self) -> Dimensions {
        match self.dimensions() {
            Dimensions::Empty => unreachable!("even a degenerate triangle should be at least 0-dimensional"),
            Dimensions::ZeroDimensional => Dimensions::Empty,
            Dimensions::OneDimensional => Dimensions::ZeroDimensional,
            Dimensions::TwoDimensional => Dimensions::OneDimensional,
            Dimensions::ThreeDimensional => unreachable!("The points of a triangle always lie on the same plane"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const ONE: Coord = Coord { x: 1.0, y: 1.0, z: 1.0 };
    use crate::wkt;

    #[test]
    fn point() {
        assert_eq!(
            Dimensions::ZeroDimensional,
            wkt!(POINT(1.0 1.0 1.0)).dimensions()
        );
    }

    #[test]
    fn line_string() {
        assert_eq!(
            Dimensions::OneDimensional,
            wkt!(LINESTRING(1.0 1.0 1.0, 2.0 2.0 2.0, 3.0 3.0 3.0)).dimensions()
        );
    }

    #[test]
    fn polygon() {
        assert_eq!(
            Dimensions::TwoDimensional,
            wkt!(POLYGON((1.0 1.0 0.0, 2.0 2.0 0.0, 3.0 3.0 0.0, 1.0 1.0 0.0))).dimensions()
        );
        assert_eq!(
            Dimensions::ThreeDimensional,
            wkt!(POLYGON((1.0 1.0 1.0, 2.0 2.0 2.0, 3.0 3.0 3.0, 1.0 1.0 1.0))).dimensions()
        );
    }

    #[test]
    fn multi_point() {
        assert_eq!(
            Dimensions::ZeroDimensional,
            wkt!(MULTIPOINT(1.0 1.0 1.0)).dimensions()
        );
    }

    #[test]
    fn multi_line_string() {
        assert_eq!(
            Dimensions::OneDimensional,
            wkt!(MULTILINESTRING((1.0 1.0 1.0, 2.0 2.0 2.0, 3.0 3.0 3.0))).dimensions()
        );
    }

    #[test]
    fn multi_polygon() {
        assert_eq!(
            Dimensions::TwoDimensional,
            wkt!(MULTIPOLYGON(((1.0 0.0 1.0, 2.0 0.0 2.0, 3.0 0.0 3.0, 1.0 0.0 1.0)))).dimensions()
        );
        assert_eq!(
            Dimensions::ThreeDimensional,
            wkt!(MULTIPOLYGON(((1.0 1.0 1.0, 2.0 2.0 2.0, 3.0 3.0 3.0, 1.0 1.0 1.0)))).dimensions()
        );
    }

    mod empty {
        use super::*;
        #[test]
        fn empty_line_string() {
            assert_eq!(
                Dimensions::Empty,
                (wkt!(LINESTRING EMPTY) as LineString<f64>).dimensions()
            );
        }

        #[test]
        fn empty_polygon() {
            assert_eq!(
                Dimensions::Empty,
                (wkt!(POLYGON EMPTY) as Polygon<f64>).dimensions()
            );
        }

        #[test]
        fn empty_multi_point() {
            assert_eq!(
                Dimensions::Empty,
                (wkt!(MULTIPOINT EMPTY) as MultiPoint<f64>).dimensions()
            );
        }

        #[test]
        fn empty_multi_line_string() {
            assert_eq!(
                Dimensions::Empty,
                (wkt!(MULTILINESTRING EMPTY) as MultiLineString<f64>).dimensions()
            );
        }

        #[test]
        fn multi_line_string_with_empty_line_string() {
            let empty_line_string = wkt!(LINESTRING EMPTY) as LineString<f64>;
            let multi_line_string = MultiLineString::new(vec![empty_line_string]);
            assert_eq!(Dimensions::Empty, multi_line_string.dimensions());
        }

        #[test]
        fn empty_multi_polygon() {
            assert_eq!(
                Dimensions::Empty,
                (wkt!(MULTIPOLYGON EMPTY) as MultiPolygon<f64>).dimensions()
            );
        }

        #[test]
        fn multi_polygon_with_empty_polygon() {
            let empty_polygon = (wkt!(POLYGON EMPTY) as Polygon<f64>);
            let multi_polygon = MultiPolygon::new(vec![empty_polygon]);
            assert_eq!(Dimensions::Empty, multi_polygon.dimensions());
        }
    }

    mod dimensional_collapse {
        use super::*;

        #[test]
        fn line_collapsed_to_point() {
            assert_eq!(
                Dimensions::ZeroDimensional,
                Line::new(ONE, ONE).dimensions()
            );
        }

        #[test]
        fn line_string_collapsed_to_point() {
            assert_eq!(
                Dimensions::ZeroDimensional,
                wkt!(LINESTRING(1.0 1.0 1.0)).dimensions()
            );
            assert_eq!(
                Dimensions::ZeroDimensional,
                wkt!(LINESTRING(1.0 1.0 1.0, 1.0 1.0 1.0)).dimensions()
            );
        }

        #[test]
        fn polygon_collapsed_to_point() {
            assert_eq!(
                Dimensions::ZeroDimensional,
                wkt!(POLYGON((1.0 1.0 1.0))).dimensions()
            );
            assert_eq!(
                Dimensions::ZeroDimensional,
                wkt!(POLYGON((1.0 1.0 1.0, 1.0 1.0 1.0))).dimensions()
            );
        }

        #[test]
        fn polygon_collapsed_to_line() {
            assert_eq!(
                Dimensions::OneDimensional,
                wkt!(POLYGON((1.0 1.0 1.0, 2.0 2.0 2.0))).dimensions()
            );
        }

        #[test]
        fn multi_line_string_with_line_string_collapsed_to_point() {
            assert_eq!(
                Dimensions::ZeroDimensional,
                wkt!(MULTILINESTRING((1.0 1.0 1.0))).dimensions()
            );
            assert_eq!(
                Dimensions::ZeroDimensional,
                wkt!(MULTILINESTRING((1.0 1.0 1.0, 1.0 1.0 1.0))).dimensions()
            );
            assert_eq!(
                Dimensions::ZeroDimensional,
                wkt!(MULTILINESTRING((1.0 1.0 1.0), (1.0 1.0 1.0))).dimensions()
            );
        }

        #[test]
        fn multi_polygon_with_polygon_collapsed_to_point() {
            assert_eq!(
                Dimensions::ZeroDimensional,
                wkt!(MULTIPOLYGON(((1.0 1.0 1.0)))).dimensions()
            );
            assert_eq!(
                Dimensions::ZeroDimensional,
                wkt!(MULTIPOLYGON(((1.0 1.0 1.0, 1.0 1.0 1.0)))).dimensions()
            );
        }

        #[test]
        fn multi_polygon_with_polygon_collapsed_to_line() {
            assert_eq!(
                Dimensions::OneDimensional,
                wkt!(MULTIPOLYGON(((1.0 1.0 1.0, 2.0 2.0 2.0)))).dimensions()
            );
        }
    }
}
