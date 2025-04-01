use std::cmp::Ordering;

use crate::area::{get_linestring_area, Area};
use crate::dimensions::{Dimensions, Dimensions::*, HasDimensions};
use crate::geometry::*;
use crate::line_measures::Length;
use crate::GeoNum;

/// Calculation of the centroid.
/// The centroid is the arithmetic mean position of all points in the shape.
/// Informally, it is the point at which a cutout of the shape could be perfectly
/// balanced on the tip of a pin.
/// The geometric centroid of a convex object always lies in the object.
/// A non-convex object might have a centroid that _is outside the object itself_.
///
/// # Examples
///
/// ```
/// use geo_3d::Centroid;
/// use geo_3d::{point, polygon};
///
/// // rhombus shaped polygon
/// let polygon = polygon![
///     (x: -2., y: 1., z: -2.),
///     (x: 1., y: 3., z: 1.),
///     (x: 4., y: 1., z: 4.),
///     (x: 1., y: -1., z: -1.),
///     (x: -2., y: 1., z: -2.),
/// ];
///
/// assert_eq!(
///     Some(point!(x: 1., y: 1., z: 1.)),
///     polygon.centroid(),
/// );
/// ```
pub trait Centroid {
    type Output;

    /// See: <https://en.wikipedia.org/wiki/Centroid>
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_3d::Centroid;
    /// use geo_3d::{line_string, point};
    ///
    /// let line_string = line_string![
    ///     (x: 40.02f64, y: 116.34),
    ///     (x: 40.02f64, y: 118.23),
    /// ];
    ///
    /// assert_eq!(
    ///     Some(point!(x: 40.02, y: 117.285, z: 40.02)),
    ///     line_string.centroid(),
    /// );
    /// ```
    fn centroid(&self) -> Self::Output;
}

impl<T> Centroid for Line<T>
where
    T: GeoNum,
{
    type Output = Point<T>;

    /// The Centroid of a [`Line`] is its middle point
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_3d::Centroid;
    /// use geo_3d::{Line, point};
    ///
    /// let line = Line::new(
    ///     point!(x: 1.0f64, y: 3.0, z: 4.0),
    ///     point!(x: 2.0f64, y: 4.0, z: -4.0),
    /// );
    ///
    /// assert_eq!(
    ///     point!(x: 1.5, y: 3.5, z: 0.0),
    ///     line.centroid(),
    /// );
    /// ```
    fn centroid(&self) -> Self::Output {
        let two = T::one() + T::one();
        (self.start_point() + self.end_point()) / two
    }
}

impl<T> Centroid for LineString<T>
where
    T: GeoNum,
{
    type Output = Option<Point<T>>;

    // The Centroid of a [`LineString`] is the mean of the middle of the segment
    // weighted by the length of the segments.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_3d::Centroid;
    /// use geo_3d::{line_string, point};
    ///
    /// let line_string = line_string![
    ///   (x: 1.0f32, y: 1.0, z: 1.0),
    ///   (x: 2.0, y: 2.0, z: 2.0),
    ///   (x: 4.0, y: 4.0, z: 4.0)
    ///   ];
    ///
    /// assert_eq!(
    ///     // (1.0 * (1.5, 1.5, 1.5) + 2.0 * (3.0, 3.0, 3.0)) / 3.0
    ///     Some(point!(x: 2.5, y: 2.5, z: 2.5)),
    ///     line_string.centroid(),
    /// );
    /// ```
    fn centroid(&self) -> Self::Output {
        let mut operation = CentroidOperation::new();
        operation.add_line_string(self);
        operation.centroid()
    }
}

impl<T> Centroid for MultiLineString<T>
where
    T: GeoNum,
{
    type Output = Option<Point<T>>;

    /// The Centroid of a [`MultiLineString`] is the mean of the centroids of all the constituent linestrings,
    /// weighted by the length of each linestring
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_3d::Centroid;
    /// use geo_3d::{MultiLineString, line_string, point};
    ///
    /// let multi_line_string = MultiLineString::new(vec![
    ///     // centroid: (2.5, 2.5, 2.5)
    ///     line_string![(x: 1.0f32, y: 1.0, z: 1.0), (x: 2.0, y: 2.0, z: 2.0), (x: 4.0, y: 4.0, z: 4.0)],
    ///     // centroid: (4.0, 4.0, 4.0)
    ///     line_string![(x: 1.0, y: 1.0, z: 1.0), (x: 3.0, y: 3.0, z: 3.0), (x: 7.0, y: 7.0, z: 7.0)],
    /// ]);
    ///
    /// assert_eq!(
    ///     // ( 3.0 * (2.5, 2.5, 2.5) + 6.0 * (4.0, 4.0, 4.0) ) / 9.0
    ///     Some(point!(x: 3.5, y: 3.5, z: 3.5)),
    ///     multi_line_string.centroid(),
    /// );
    /// ```
    fn centroid(&self) -> Self::Output {
        let mut operation = CentroidOperation::new();
        operation.add_multi_line_string(self);
        operation.centroid()
    }
}

impl<T> Centroid for Polygon<T>
where
    T: GeoNum,
{
    type Output = Option<Point<T>>;

    /// The Centroid of a [`Polygon`] is the mean of its points
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_3d::Centroid;
    /// use geo_3d::{polygon, point};
    ///
    /// let polygon = polygon![
    ///     (x: 0.0f32, y: 0.0, z: 0.0),
    ///     (x: 2.0, y: 0.0, z: ),
    ///     (x: 2.0, y: 1.0, z: ),
    ///     (x: 0.0, y: 1.0, z: ),
    /// ];
    ///
    /// assert_eq!(
    ///     Some(point!(x: 1.0, y: 0.5, z: )),
    ///     polygon.centroid(),
    /// );
    /// ```
    fn centroid(&self) -> Self::Output {
        let mut operation = CentroidOperation::new();
        operation.add_polygon(self);
        operation.centroid()
    }
}

impl<T> Centroid for MultiPolygon<T>
where
    T: GeoNum,
{
    type Output = Option<Point<T>>;

    /// The Centroid of a [`MultiPolygon`] is the mean of the centroids of its polygons, weighted
    /// by the area of the polygons
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_3d::Centroid;
    /// use geo_3d::{MultiPolygon, polygon, point};
    ///
    /// let multi_polygon = MultiPolygon::new(vec![
    ///   // centroid (1.0, 0.5, z: 1.0)
    ///   polygon![
    ///     (x: 0.0f32, y: 0.0, z: 0.0),
    ///     (x: 2.0, y: 0.0, z: ),
    ///     (x: 2.0, y: 1.0, z: ),
    ///     (x: 0.0, y: 1.0, z: ),
    ///   ],
    ///   // centroid (-0.5, 0.0, )
    ///   polygon![
    ///     (x: 1.0, y: 1.0, z: 1.0),
    ///     (x: -2.0, y: 1.0, z: ),
    ///     (x: -2.0, y: -1.0, z: ),
    ///     (x: 1.0, y: -1.0, z: ),
    ///   ]
    /// ]);
    ///
    /// assert_eq!(
    ///     // ( 2.0 * (1.0, 0.5, ) + 6.0 * (-0.5, 0.0, ) ) / 8.0
    ///     Some(point!(x: -0.125, y: 0.125, z: -0.125)),
    ///     multi_polygon.centroid(),
    /// );
    /// ```
    fn centroid(&self) -> Self::Output {
        let mut operation = CentroidOperation::new();
        operation.add_multi_polygon(self);
        operation.centroid()
    }
}

impl<T> Centroid for Rect<T>
where
    T: GeoNum,
{
    type Output = Point<T>;

    /// The Centroid of a [`Rect`] is the mean of its [`Point`]s
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_3d::Centroid;
    /// use geo_3d::{Rect, point};
    ///
    /// let rect = Rect::new(
    ///   point!(x: 0.0f32, y: 0.0, z: 0.0),
    ///   point!(x: 1.0, y: 1.0, z: 1.0),
    /// );
    ///
    /// assert_eq!(
    ///     point!(x: 0.5, y: 0.5, z: 0.5),
    ///     rect.centroid(),
    /// );
    /// ```
    fn centroid(&self) -> Self::Output {
        self.center().into()
    }
}

impl<T> Centroid for Triangle<T>
where
    T: GeoNum,
{
    type Output = Point<T>;

    /// The Centroid of a [`Triangle`] is the mean of its [`Point`]s
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_3d::Centroid;
    /// use geo_3d::{Triangle, coord, point};
    ///
    /// let triangle = Triangle::new(
    ///   coord!(x: 0.0f32, y: -1.0, z: 0.0),
    ///   coord!(x: 3.0, y: 0.0, z: 3.0),
    ///   coord!(x: 0.0, y: 1.0, z: 0.0),
    /// );
    ///
    /// assert_eq!(
    ///     point!(x: 1.0, y: 0.0, z: 1.0),
    ///     triangle.centroid(),
    /// );
    /// ```
    fn centroid(&self) -> Self::Output {
        let mut operation = CentroidOperation::new();
        operation.add_triangle(self);
        operation
            .centroid()
            .expect("triangle cannot have an empty centroid")
    }
}

impl<T> Centroid for Point<T>
where
    T: GeoNum,
{
    type Output = Point<T>;

    /// The Centroid of a [`Point`] is the point itself
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_3d::Centroid;
    /// use geo_3d::point;
    ///
    /// let point = point!(x: 1.0f32, y: 2.0, z: 3.0);
    ///
    /// assert_eq!(
    ///     point!(x: 1.0f32, y: 2.0, z: 3.0),
    ///     point.centroid(),
    /// );
    /// ```
    fn centroid(&self) -> Self::Output {
        *self
    }
}

impl<T> Centroid for MultiPoint<T>
where
    T: GeoNum,
{
    type Output = Option<Point<T>>;

    /// The Centroid of a [`MultiPoint`] is the mean of all [`Point`]s
    ///
    /// # Example
    ///
    /// ```
    /// use geo_3d::Centroid;
    /// use geo_3d::{MultiPoint, Point};
    ///
    /// let empty: Vec<Point> = Vec::new();
    /// let empty_multi_points: MultiPoint<_> = empty.into();
    /// assert_eq!(empty_multi_points.centroid(), None);
    ///
    /// let points: MultiPoint<_> = vec![(5., 1., 5.), (1., 3., 1.), (3., 2., 3.)].into();
    /// assert_eq!(points.centroid(), Some(Point::new(3., 2., 3.)));
    /// ```
    fn centroid(&self) -> Self::Output {
        let mut operation = CentroidOperation::new();
        operation.add_multi_point(self);
        operation.centroid()
    }
}

impl<T> Centroid for Geometry<T>
where
    T: GeoNum,
{
    type Output = Option<Point<T>>;

    crate::geometry_delegate_impl! {
        /// The Centroid of a [`Geometry`] is the centroid of its enum variant
        ///
        /// # Examples
        ///
        /// ```
        /// use geo_3d::Centroid;
        /// use geo_3d::{Geometry, Rect, point};
        ///
        /// let rect = Rect::new(
        ///   point!(x: 0.0f32, y: 0.0, z: 0.0),
        ///   point!(x: 1.0, y: 1.0, z: 1.0),
        /// );
        /// let geometry = Geometry::from(rect.clone());
        ///
        /// assert_eq!(
        ///     Some(rect.centroid()),
        ///     geometry.centroid(),
        /// );
        ///
        /// assert_eq!(
        ///     Some(point!(x: 0.5, y: 0.5, z: 0.5)),
        ///     geometry.centroid(),
        /// );
        /// ```
        fn centroid(&self) -> Self::Output;
    }
}

impl<T> Centroid for GeometryCollection<T>
where
    T: GeoNum,
{
    type Output = Option<Point<T>>;

    /// The Centroid of a [`GeometryCollection`] is the mean of the centroids of elements, weighted
    /// by the area of its elements.
    ///
    /// Note that this means, that elements which have no area are not considered when calculating
    /// the centroid.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_3d::Centroid;
    /// use geo_3d::{Geometry, GeometryCollection, Rect, Triangle, point, coord};
    ///
    /// let rect_geometry = Geometry::from(Rect::new(
    ///   point!(x: 0.0f32, y: 0.0, z: 0.0),
    ///   point!(x: 1.0, y: 1.0, z: 1.0),
    /// ));
    ///
    /// let triangle_geometry = Geometry::from(Triangle::new(
    ///     coord!(x: 0.0f32, y: -1.0),
    ///     coord!(x: 3.0, y: 0.0, z: -2.0),
    ///     coord!(x: 0.0, y: 1.0, z: 2.0),
    /// ));
    ///
    /// let point_geometry = Geometry::from(
    ///   point!(x: 12351.0, y: 129815.0, z: 0.23597)
    /// );
    ///
    /// let geometry_collection = GeometryCollection::new(
    ///   vec![
    ///     rect_geometry,
    ///     triangle_geometry,
    ///     point_geometry
    ///   ]
    /// );
    ///
    /// assert_eq!(
    ///     Some(point!(x: 0.875, y: 0.125, z: 0.5)),
    ///     geometry_collection.centroid(),
    /// );
    /// ```
    fn centroid(&self) -> Self::Output {
        let mut operation = CentroidOperation::new();
        operation.add_geometry_collection(self);
        operation.centroid()
    }
}

struct CentroidOperation<T: GeoNum>(Option<WeightedCentroid<T>>);
impl<T: GeoNum> CentroidOperation<T> {
    const fn new() -> Self {
        CentroidOperation(None)
    }

    fn centroid(&self) -> Option<Point<T>> {
        self.0.as_ref().map(|weighted_centroid| {
            Point::from(weighted_centroid.accumulated / weighted_centroid.weight)
        })
    }

    fn centroid_dimensions(&self) -> Dimensions {
        self.0
            .as_ref()
            .map(|weighted_centroid| weighted_centroid.dimensions)
            .unwrap_or(Empty)
    }

    fn add_coord(&mut self, coord: Coord<T>) {
        self.add_centroid(ZeroDimensional, coord, T::one());
    }

    fn add_line(&mut self, line: &Line<T>) {
        match line.dimensions() {
            ZeroDimensional => self.add_coord(line.start),
            OneDimensional => self.add_centroid(
                OneDimensional,
                line.centroid().0,
                line.length(),
            ),
            TwoDimensional => self.add_centroid(
                TwoDimensional,
                line.centroid().0,
                line.length(),
            ),
            ThreeDimensional => self.add_centroid(
                ThreeDimensional,
                line.centroid().0,
                line.length(),
            ),
            Empty => unreachable!("Line have a dimension"),
        }
    }

    fn add_line_string(&mut self, line_string: &LineString<T>) {
        if self.centroid_dimensions() > OneDimensional {
            return;
        }

        if line_string.0.len() == 1 {
            self.add_coord(line_string.0[0]);
            return;
        }

        for line in line_string.lines() {
            self.add_line(&line);
        }
    }

    fn add_multi_line_string(&mut self, multi_line_string: &MultiLineString<T>) {
        if self.centroid_dimensions() > OneDimensional {
            return;
        }

        for element in &multi_line_string.0 {
            self.add_line_string(element);
        }
    }

    fn add_polygon(&mut self, polygon: &Polygon<T>) {
        // Polygons which are completely covered by their interior rings have zero area, and
        // represent a unique degeneracy into a line_string which cannot be handled by accumulating
        // directly into `self`. Instead, we perform a sub-operation, inspect the result, and only
        // then incorporate the result into `self.

        let mut exterior_operation = CentroidOperation::new();
        exterior_operation.add_ring(polygon.exterior());

        let mut interior_operation = CentroidOperation::new();
        for interior in polygon.interiors() {
            interior_operation.add_ring(interior);
        }

        if let Some(exterior_weighted_centroid) = exterior_operation.0 {
            let mut poly_weighted_centroid = exterior_weighted_centroid;
            if let Some(interior_weighted_centroid) = interior_operation.0 {
                poly_weighted_centroid.sub_assign(interior_weighted_centroid);
                if poly_weighted_centroid.weight.is_zero() {
                    // A polygon with no area `interiors` completely covers `exterior`, degenerating to a linestring
                    self.add_line_string(polygon.exterior());
                    return;
                }
            }
            self.add_weighted_centroid(poly_weighted_centroid);
        }
    }

    fn add_multi_point(&mut self, multi_point: &MultiPoint<T>) {
        if self.centroid_dimensions() > ZeroDimensional {
            return;
        }

        for element in &multi_point.0 {
            self.add_coord(element.0);
        }
    }

    fn add_multi_polygon(&mut self, multi_polygon: &MultiPolygon<T>) {
        for element in &multi_polygon.0 {
            self.add_polygon(element);
        }
    }

    fn add_geometry_collection(&mut self, geometry_collection: &GeometryCollection<T>) {
        for element in &geometry_collection.0 {
            self.add_geometry(element);
        }
    }

    fn add_rect(&mut self, rect: &Rect<T>) {
        match rect.dimensions() {
            ZeroDimensional => self.add_coord(rect.min()),
            OneDimensional => {
                // Degenerate rect is a line, treat it the same way we treat flat polygons
                self.add_line(&Line::new(rect.min(), rect.min()));
                self.add_line(&Line::new(rect.min(), rect.max()));
                self.add_line(&Line::new(rect.max(), rect.max()));
                self.add_line(&Line::new(rect.max(), rect.min()));
            }
            TwoDimensional => {
                self.add_centroid(TwoDimensional, rect.centroid().0, rect.unsigned_area())
            }
            ThreeDimensional => {
                self.add_centroid(ThreeDimensional, rect.centroid().0, rect.unsigned_area())
            }
            Empty => unreachable!("Rect dimensions cannot be empty"),
        }
    }

    fn add_triangle(&mut self, triangle: &Triangle<T>) {
        match triangle.dimensions() {
            ZeroDimensional => self.add_coord(triangle.0),
            OneDimensional => {
                // Degenerate triangle is a line, treat it the same way we treat flat
                // polygons
                let l0_1 = Line::new(triangle.0, triangle.1);
                let l1_2 = Line::new(triangle.1, triangle.2);
                let l2_0 = Line::new(triangle.2, triangle.0);
                self.add_line(&l0_1);
                self.add_line(&l1_2);
                self.add_line(&l2_0);
            }
            TwoDimensional => {
                let centroid = (triangle.0 + triangle.1 + triangle.2) / T::from(3).unwrap();
                self.add_centroid(TwoDimensional, centroid, triangle.unsigned_area());
            }
            ThreeDimensional => {
                let centroid = (triangle.0 + triangle.1 + triangle.2) / T::from(3).unwrap();
                self.add_centroid(ThreeDimensional, centroid, triangle.unsigned_area());
            }
            Empty => unreachable!("Tri dimensions cannot be empty"),
        }
    }

    fn add_geometry(&mut self, geometry: &Geometry<T>) {
        match geometry {
            Geometry::Point(g) => self.add_coord(g.0),
            Geometry::Line(g) => self.add_line(g),
            Geometry::LineString(g) => self.add_line_string(g),
            Geometry::Polygon(g) => self.add_polygon(g),
            Geometry::MultiPoint(g) => self.add_multi_point(g),
            Geometry::MultiLineString(g) => self.add_multi_line_string(g),
            Geometry::MultiPolygon(g) => self.add_multi_polygon(g),
            Geometry::GeometryCollection(g) => self.add_geometry_collection(g),
            Geometry::Rect(g) => self.add_rect(g),
            Geometry::Triangle(g) => self.add_triangle(g),
        }
    }

    fn add_ring(&mut self, ring: &LineString<T>) {
        debug_assert!(ring.is_closed());

        let area = get_linestring_area(ring);
        if area == T::zero() {
            match ring.dimensions() {
                // empty ring doesn't contribute to centroid
                Empty => {}
                // degenerate ring is a point
                ZeroDimensional => self.add_coord(ring[0]),
                // zero-area ring is a line string
                _ => self.add_line_string(ring),
            }
            return;
        }

        // Since area is non-zero, we know the ring has at least one point
        let shift = ring.0[0];

        let accumulated_coord = ring.lines().fold(Coord::zero(), |accum, line| {
            use crate::MapCoords;
            let line = line.map_coords(|c| c - shift);
            let tmp = line.determinant();
            accum + (line.end + line.start) * tmp
        });
        let six = T::from(6).unwrap();
        let centroid = accumulated_coord / (six * area) + shift;
        let weight = area.abs();
        self.add_centroid(TwoDimensional, centroid, weight);
    }

    fn add_centroid(&mut self, dimensions: Dimensions, centroid: Coord<T>, weight: T) {
        let weighted_centroid = WeightedCentroid {
            dimensions,
            weight,
            accumulated: centroid * weight,
        };
        self.add_weighted_centroid(weighted_centroid);
    }

    fn add_weighted_centroid(&mut self, other: WeightedCentroid<T>) {
        match self.0.as_mut() {
            Some(centroid) => centroid.add_assign(other),
            None => self.0 = Some(other),
        }
    }
}

// Aggregated state for accumulating the centroid of a geometry or collection of geometries.
struct WeightedCentroid<T: GeoNum> {
    weight: T,
    accumulated: Coord<T>,
    /// Collections of Geometries can have different dimensionality. Centroids must be considered
    /// separately by dimensionality.
    ///
    /// e.g. If I have several Points, adding a new `Point` will affect their centroid.
    ///
    /// However, because a Point is zero dimensional, it is infinitely small when compared to
    /// any 2-D Polygon. Thus a Point will not affect the centroid of any GeometryCollection
    /// containing a 2-D Polygon.
    ///
    /// So, when accumulating a centroid, we must track the dimensionality of the centroid
    dimensions: Dimensions,
}

impl<T: GeoNum> WeightedCentroid<T> {
    fn add_assign(&mut self, b: WeightedCentroid<T>) {
        match self.dimensions.cmp(&b.dimensions) {
            Ordering::Less => *self = b,
            Ordering::Greater => {}
            Ordering::Equal => {
                self.accumulated = self.accumulated + b.accumulated;
                self.weight = self.weight + b.weight;
            }
        }
    }

    fn sub_assign(&mut self, b: WeightedCentroid<T>) {
        match self.dimensions.cmp(&b.dimensions) {
            Ordering::Less => *self = b,
            Ordering::Greater => {}
            Ordering::Equal => {
                self.accumulated = self.accumulated - b.accumulated;
                self.weight = self.weight - b.weight;
            }
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{coord, line_string, point, polygon, wkt};

    // --- Centroid of LineString ---

    #[test]
    fn empty_linestring_test() {
        let linestring: LineString<f32> = line_string![];
        let centroid = linestring.centroid();
        assert!(centroid.is_none());
    }

    #[test]
    fn linestring_one_point_test() {
        let coord = coord! {
            x: 40.02f64,
            y: 116.34,
            z: 10.0,
        };
        let linestring = line_string![coord];
        let centroid = linestring.centroid();
        assert_eq!(centroid, Some(Point::from(coord)));
    }

    #[test]
    fn linestring_test() {
        let linestring = line_string![
            (x: 1., y: 1., z: 1.),
            (x: 7., y: 1., z: -7.),
            (x: 8., y: 1., z: -8.),
            (x: 9., y: 1., z: -9.),
            (x: 10., y: 1., z: -10.),
            (x: 11., y: 1., z: -11.)
        ];
        // The weighted–average of the segment midpoints gives (6,1,–6)
        assert_eq!(linestring.centroid(), Some(point!(x: 6., y: 1., z: -6.)));
    }

    #[test]
    fn linestring_with_repeated_point_test() {
        let l1 = LineString::from(vec![
            point!(1., 1., 1.),
            point!(1., 1., 1.),
            point!(1., 1., 1.),
        ]);
        assert_eq!(l1.centroid(), Some(point!(1., 1., 1.)));

        let l2 = LineString::from(vec![
            point!(2., 2., 2.),
            point!(2., 2., 2.),
            point!(2., 2., 2.),
        ]);
        let mls = MultiLineString::new(vec![l1, l2]);
        // The two linestrings have lengths zero so the centroid is the average of (1,1,1) and (2,2,2)
        assert_eq!(mls.centroid(), Some(point!(1.5, 1.5, 1.5)));
    }

    // --- Centroid of MultiLineString ---

    #[test]
    fn empty_multilinestring_test() {
        let mls: MultiLineString = MultiLineString::new(vec![]);
        let centroid = mls.centroid();
        assert!(centroid.is_none());
    }

    #[test]
    fn multilinestring_with_empty_line_test() {
        let mls: MultiLineString = MultiLineString::new(vec![line_string![]]);
        let centroid = mls.centroid();
        assert!(centroid.is_none());
    }

    #[test]
    fn multilinestring_length_0_test() {
        let coord = coord! {
            x: 40.02f64,
            y: 116.34,
            z: 15.0,
        };
        let mls: MultiLineString<f64> = MultiLineString::new(vec![
            line_string![coord],
            line_string![coord],
            line_string![coord],
        ]);
        assert_relative_eq!(mls.centroid().unwrap(), Point::from(coord));
    }

    #[test]
    fn multilinestring_one_line_test() {
        let linestring = line_string![
            (x: 1., y: 1., z: 1.),
            (x: 7., y: 1., z: 1.),
            (x: 8., y: 1., z: 8.),
            (x: 9., y: 1., z: 8.),
            (x: 10., y: 1., z: 1.),
            (x: 11., y: 1., z: 1.)
        ];
        let mls: MultiLineString<f64> = MultiLineString::new(vec![linestring]);

        assert_relative_eq!(mls.centroid().unwrap(), point!(x: 6., y: 1., z: 1.));
    }

    #[test]
    fn multilinestring_test() {
        let mls = wkt! {
            MULTILINESTRING(
                (0.0 0.0 0.0, 1.0 10.0 100.0),
                (1.0 10.0 100.0, 2.0 0.0 -2.0, 3.0 1.0 -3.0),
                (-12.0 -100.0 -12.0, 7.0 8.0 9.0)
            )
        };
        // (Note: The numbers below are computed by the centroid algorithm over all segments.)
        // For example, the first linestring contributes a midpoint of (0.5,5,50) weighted by its length, etc.
        // The final weighted–centroid (x,y,z) is approximately:
        assert_relative_eq!(
            mls.centroid().unwrap(),
            point![x: -0.2238, y: -13.024, z: 31.208],
            max_relative = 1e-3
        );
    }

    // --- Centroid of Polygon ---

    #[test]
    fn empty_polygon_test() {
        let poly: Polygon<f64> = polygon![];
        assert!(poly.centroid().is_none());
    }

    #[test]
    fn polygon_one_point_test() {
        // A polygon defined by one point
        let p = point![ x: 2., y: 1., z: 5. ];
        let poly = polygon![p.0];
        assert_relative_eq!(poly.centroid().unwrap(), p);
    }

    #[test]
    fn centroid_polygon_numerical_stability() {
        let polygon = {
            use std::f64::consts::PI;
            const NUM_VERTICES: usize = 10;
            const ANGLE_INC: f64 = 2. * PI / NUM_VERTICES as f64;

            // For variety we set z = cos(angle) + sin(angle)
            Polygon::new(
                (0..NUM_VERTICES)
                    .map(|i| {
                        let angle = i as f64 * ANGLE_INC;
                        coord! {
                            x: angle.cos(),
                            y: angle.sin(),
                            z: angle.cos() + angle.sin(),
                        }
                    })
                    .collect::<Vec<_>>()
                    .into(),
                vec![],
            )
        };

        let centroid = polygon.centroid().unwrap();

        let shift = coord! { x: 1.5e8, y: 1.5e8, z: 1.5e8 };

        use crate::map_coords::MapCoords;
        let polygon = polygon.map_coords(|c| c + shift);
        let new_centroid = polygon.centroid().unwrap().map_coords(|c| c - shift);
        debug!("centroid {:?}", centroid.0);
        debug!("new_centroid {:?}", new_centroid.0);
        assert_relative_eq!(centroid.0.x, new_centroid.0.x, max_relative = 0.0001);
        assert_relative_eq!(centroid.0.y, new_centroid.0.y, max_relative = 0.0001);
        assert_relative_eq!(centroid.0.z, new_centroid.0.z, max_relative = 0.0001);
    }

    #[test]
    fn polygon_test() {
        // Here we define a quadrilateral whose vertices are chosen so that its centroid is (1,1,1).
        // We choose:
        //   A = (0,0,0), B = (2,0,2), C = (2,2,2), D = (0,2,0) and then A again.
        let poly = polygon![
            (x: 0., y: 0., z: 0.),
            (x: 2., y: 0., z: 2.),
            (x: 2., y: 2., z: 2.),
            (x: 0., y: 2., z: 0.),
            (x: 0., y: 0., z: 0.)
        ];
        assert_relative_eq!(poly.centroid().unwrap(), point![x: 1., y: 1., z: 1.]);
    }

    #[test]
    fn polygon_hole_test() {
        // Here we use a WKT string for a polygon with two interior rings.
        // To keep things simple, we supply a constant z–value (z = 1) for every coordinate so that
        // the expected centroid is the same in all three dimensions.
        let p1 = wkt! { POLYGON(
            (5.0 1.0 1.0, 4.0 2.0 1.0, 4.0 3.0 1.0, 5.0 4.0 1.0, 6.0 4.0 1.0, 7.0 3.0 1.0, 7.0 2.0 1.0, 6.0 1.0 1.0, 5.0 1.0 1.0),
            (5.0 1.3 1.0, 5.5 2.0 1.0, 6.0 1.3 1.0, 5.0 1.3 1.0),
            (5.0 2.3 1.0, 5.5 3.0 1.0, 6.0 2.3 1.0, 5.0 2.3 1.0)
        ) };

        let centroid = p1.centroid().unwrap();
        assert_relative_eq!(centroid, point!(x: 5.5, y: 2.5518518518518523, z: 1.0), max_relative = 1e-6);
    }

    #[test]
    fn flat_polygon_test() {
        // A degenerate (flat) polygon defined by three collinear points.
        let poly = wkt! { POLYGON((0. 1. 2., 1. 1. 1., 0. 1. 2.)) };
        assert_eq!(poly.centroid(), Some(point!(0.5, 1., 1.5)));
    }

    #[test]
    fn multi_poly_with_flat_polygon_test() {
        let multipoly = wkt! { MULTIPOLYGON(((0. 0. 0., 1. 0. -1., 0. 0. 0.))) };

        assert_eq!(multipoly.centroid(), Some(point!(0.5, 0., -0.5)));
    }

    #[test]
    fn multi_poly_with_multiple_flat_polygon_test() {
        let multipoly = wkt! { MULTIPOLYGON(
            ((1. 1. 1., 1. 3. 5., 1. 1. 1.)),
            ((2. 2. 2., 6. 2. 6., 2. 2. 2.))
        )};
        // Here the expected overall centroid (averaging the two polygons)
        assert_eq!(multipoly.centroid(), Some(point!(3., 2., 1.)));
    }

    #[test]
    fn multi_poly_with_only_points_test() {
        let p1 = wkt! { POLYGON((1. 1. 1., 1. 1. 1., 1. 1. 1.)) };
        assert_eq!(p1.centroid(), Some(point!(1., 1., 1.)));

        let multipoly = wkt! { MULTIPOLYGON(
            ((1. 1. 1.,1. 1. 1.,1. 1. 1.)),
            ((2. 2. 2.,2. 2. 2.,2. 2. 2.))
        ) };
        // The weighted average of (1,1,1) and (2,2,2) is (1.5,1.5,1.5)
        assert_eq!(multipoly.centroid(), Some(point!(1.5, 1.5, 1.5)));
    }

    #[test]
    fn multi_poly_with_one_ring_and_one_real_poly() {
        // When a multipolygon is composed of a "normal" polygon (non–zero area) and a degenerate ring,
        // the centroid of the multipolygon is that of the normal polygon.
        let normal = Polygon::new(
            LineString::from(vec![
                point!(1., 1., 1.),
                point!(1., 3., -1.),
                point!(3., 1., -3.),
                point!(1., 1., 1.),
            ]),
            vec![],
        );
        let flat = Polygon::new(
            LineString::from(vec![
                point!(2., 2., 2.),
                point!(6., 2., -6.),
                point!(2., 2., 2.),
            ]),
            vec![],
        );
        let multipoly = MultiPolygon::new(vec![normal.clone(), flat]);
        assert_eq!(multipoly.centroid(), normal.centroid());
    }

    #[test]
    fn polygon_flat_interior_test() {
        let poly = Polygon::new(
            LineString::from(vec![
                point!(0., 0., 0.),
                point!(0., 1., 0.),
                point!(1., 1., 1.),
                point!(1., 0., 1.),
                point!(0., 0., 0.)
            ]),
            vec![
                LineString::from(vec![
                    point!(0., 0., 0.),
                    point!(0., 1., 1.),
                    point!(0., 0., 0.),
                ])
            ],
        );
        assert_eq!(poly.centroid(), Some(point!(0.5, 0.5, 0.5)));
    }

    #[test]
    fn empty_interior_polygon_test() {
        let poly = Polygon::new(
            LineString::from(vec![
                point!(0., 0., 0.),
                point!(0., 1., 1.),
                point!(1., 1., 1.),
                point!(1., 0., 0.),
                point!(0., 0., 0.)
            ]),
            vec![LineString::new(vec![])],
        );
        assert_eq!(poly.centroid(), Some(point!(0.5, 0.5, 0.5)));
    }

    #[test]
    fn polygon_ring_test() {
        let square = LineString::from(vec![
            point!(0., 0., 0.),
            point!(0., 1., 1.),
            point!(1., 1., 1.),
            point!(1., 0., 0.),
            point!(0., 0., 0.)
        ]);
        let poly = Polygon::new(square.clone(), vec![square]);
        assert_eq!(poly.centroid(), Some(point!(0.5, 0.5, 0.5)));
    }

    #[test]
    fn polygon_cell_test() {
        // A polygon with interior rings (holes) that cause it to have zero net area.
        let square = LineString::from(vec![
            point!(0., 0., 0.),
            point!(0., 2., 2.),
            point!(2., 2., 2.),
            point!(2., 0., 0.),
            point!(0., 0., 0.)
        ]);
        let bottom = LineString::from(vec![
            point!(0., 0., 0.),
            point!(2., 0., 0.),
            point!(2., 1., 2.),
            point!(0., 1., 1.),
            point!(0., 0., 0.)
        ]);
        let top = LineString::from(vec![
            point!(0., 1., 1.),
            point!(2., 1., 1.),
            point!(2., 2., 2.),
            point!(0., 2., 2.),
            point!(0., 1., 1.)
        ]);
        let poly = Polygon::new(square, vec![top, bottom]);
        assert_eq!(poly.centroid(), Some(point!(1., 1., 1.)));
    }

    // --- Centroid of MultiPolygon ---

    #[test]
    fn empty_multipolygon_polygon_test() {
        assert!(MultiPolygon::<f64>::new(Vec::new()).centroid().is_none());
    }

    #[test]
    fn multipolygon_one_polygon_test() {
        let linestring = LineString::from(vec![
            point!(0., 0., 0.),
            point!(2., 0., 2.),
            point!(2., 2., 2.),
            point!(0., 2., 0.),
            point!(0., 0., 0.)
        ]);
        let poly = Polygon::new(linestring, Vec::new());
        assert_eq!(MultiPolygon::new(vec![poly]).centroid(), Some(point!(1., 1., 1.)));
    }

    #[test]
    fn multipolygon_two_polygons_test() {
        let linestring = LineString::from(vec![
            point!(2., 1., -2.),
            point!(5., 1., -5.),
            point!(5., 3., -5.),
            point!(2., 3., 4.),
            point!(2., 1., -2.)
        ]);
        let poly1 = Polygon::new(linestring, Vec::new());
        // For the second polygon, fill in missing z–values with 0.
        let linestring = LineString::from(vec![
            point!(7., 1., 7.),
            point!(8., 1., 7.),
            point!(8., 2., 7.),
            point!(7., 2., 7.),
            point!(7., 1., 7.)
        ]);
        let poly2 = Polygon::new(linestring, Vec::new());
        let centroid = MultiPolygon::new(vec![poly1, poly2]).centroid().unwrap();

        assert_relative_eq!(
            centroid,
            point![x: 3.477235635767532, y: 2.1334483508915114, z: -1.25808368107084],
            max_relative = 1e-6
        );
    }

    #[test]
    fn multipolygon_two_polygons_of_opposite_clockwise_test() {
        let linestring = LineString::from(vec![
            (0., 0., 0.),
            (2., 0., 0.),
            (2., 2., 2.),
            (0., 2., 0.),
            (0., 0., 0.)
        ]);
        let poly1 = Polygon::new(linestring, Vec::new());
        let linestring = LineString::from(vec![
            (0., 0., 0.),
            (-2., 0., 0.),
            (-2., 2., 0.),
            (0., 2., 0.),
            (0., 0., 0.)
        ]);
        let poly2 = Polygon::new(linestring, Vec::new());
        assert_relative_eq!(
            MultiPolygon::new(vec![poly1, poly2]).centroid().unwrap(),
            point![x: 0., y: 1., z: 2.]
        );
    }

    #[test]
    fn bounding_rect_test() {
        let bounding_rect = Rect::new(
            coord! { x: 0., y: 50., z: 0. },
            coord! { x: 4., y: 100., z: 0. }
        );
        let point = point![x: 2., y: 75., z: 0.];
        assert_eq!(point, bounding_rect.centroid());
    }

    #[test]
    fn line_test() {
        let line1 = Line::new(
            coord!(0., 1., 2.),
            coord!(1., 3., 6.)
        );
        // The centroid (midpoint) is (0.5,2,4)
        assert_eq!(line1.centroid(), point![x: 0.5, y: 2., z: 4.]);
    }

    #[test]
    fn collection_weighting() {
        // Here we define a MultiPoint and a GeometryCollection so that we check weighting across dimensions.
        let p0 = point!(x: 0.0, y: 0.0, z: 0.0);
        let p1 = point!(x: 2.0, y: 0.0, z: 0.0);
        let p2 = point!(x: 2.0, y: 2.0, z: 2.0);
        let p3 = point!(x: 0.0, y: 2.0, z: 2.0);

        let multi_point = MultiPoint::new(vec![p0, p1, p2, p3]);
        assert_eq!(multi_point.centroid().unwrap(), point!(x: 1.0, y: 1.0, z: 1.0));

        let collection =
            GeometryCollection::new(vec![MultiPoint::new(vec![p1, p2, p3]).into(), p0.into()]);
        assert_eq!(collection.centroid().unwrap(), point!(x: 1.0, y: 1.0, z: 1.0));
    }

    #[test]
    fn triangles() {
        // Boring triangle:
        assert_abs_diff_eq!(
            Triangle::new(
                coord!(0., 0., 0.),
                coord!(3., 0., -3.),
                coord!(1.5, 3., -1.5)
            ).centroid(),
            point!(x: 1.5, y: 1.0, z: -1.5),
            epsilon = 0.000000000000001,
        );

        // Flat triangle:
        assert_eq!(
            Triangle::new(
                coord!(0., 0., 0.),
                coord!(3., 0., -3.),
                coord!(1., 0., -1.)
            ).centroid(),
            point!(x: 1.5, y: 0.0, z: -1.5)
        );

        // Flat triangle that’s not axis–aligned:
        assert_eq!(
            Triangle::new(
                coord!(0., 0., 0.),
                coord!(3., 3., 3.),
                coord!(1., 1., 1.)
            ).centroid(),
            point!(x: 1.5, y: 1.5, z: 1.5)
        );

        // Triangle with some repeated points:
        assert_eq!(
            Triangle::new(
                coord!(0., 0., 0.),
                coord!(0., 0., 0.),
                coord!(1., 0., -1.)
            ).centroid(),
            point!(x: 0.5, y: 0.0, z: -0.5)
        );

        // Triangle with all repeated points:
        assert_eq!(
            Triangle::new(
                coord!(0.0, 0.5, 1.0),
                coord!(0.0, 0.5, 1.0),
                coord!(0.0, 0.5, 1.0)
            ).centroid(),
            point!(x: 0., y: 0.5, z: 1.0)
        );
    }

    #[test]
    fn degenerate_triangle_like_ring() {
        let triangle = Triangle::new(
            coord!(0., 0., 0.),
            coord!(1., 1., 1.),
            coord!(2., 2., 2.)
        );
        let poly: Polygon<_> = triangle.into();

        let line = Line::new(
            coord!(0., 1., 2.),
            coord!(1., 3., 5.)
        );

        let g1 = GeometryCollection::new(vec![triangle.into(), line.into()]);
        let g2 = GeometryCollection::new(vec![poly.into(), line.into()]);
        // Although the two collections use different approaches, their centroids should agree.
        assert_eq!(g1.centroid(), g2.centroid());
    }

    #[test]
    fn degenerate_rect_like_ring() {
        let rect = Rect::new(
            coord!(1., 0., 3.),
            coord!(0., 4., 0.)
        );
        let poly: Polygon<_> = rect.into();

        let line = Line::new(
            coord!(0., -1., 2.),
            coord!(1., 3., 5.)
        );

        let g1 = GeometryCollection::new(vec![rect.into(), line.into()]);
        let g2 = GeometryCollection::new(vec![poly.into(), line.into()]);
        assert_eq!(g1.centroid(), g2.centroid());
    }

    #[test]
    fn rectangles() {
        // Boring rect:
        assert_eq!(
            Rect::new(
                coord!(0., 0., 0.),
                coord!(4., 4., 4.)
            ).centroid(),
            point!(x: 2.0, y: 2.0, z: 2.0)
        );

        // Flat rect:
        assert_eq!(
            Rect::new(
                coord!(0., 0., 0.),
                coord!(4., 0., 8.)
            ).centroid(),
            point!(x: 2.0, y: 0.0, z: 4.0)
        );

        // Rect with all repeated points:
        assert_eq!(
            Rect::new(
                coord!(4., 4., 4.),
                coord!(4., 4., 4.)
            ).centroid(),
            point!(x: 4., y: 4., z: 4.)
        );

        // Collection with rect:
        let mut collection = GeometryCollection::new(vec![
            point!(0., 0., 0.).into(),
            point!(6., 0., -6.).into(),
            point!(6., 6., 6.).into(),
        ]);
        // sanity check
        assert_eq!(collection.centroid().unwrap(), point!(x: 4., y: 2., z: 0.));

        // 0-d rect treated like point
        collection.0.push(Rect::new(coord!(0., 6., 0.), coord!(0., 6., 0.)).into());
        assert_eq!(collection.centroid().unwrap(), point!(x: 3., y: 3., z: 0.));

        // 1-d rect treated like line. Since a line has higher dimensions than the rest of the
        // collection, its centroid clobbers everything else in the collection.
        collection.0.push(Rect::new(coord!(0., 0., 0.), coord!(1., 2., 4.)).into());
        assert_eq!(collection.centroid().unwrap(), point!(x: 0.5, y: 1., z: 2.));

        // 2-d has higher dimensions than the rest of the collection, so its centroid clobbers
        // everything else in the collection.
        collection
            .0
            .push(Rect::new(coord!(10., 10., 0.), coord!(11., 11., 0.)).into());
        assert_eq!(collection.centroid().unwrap(), point!(x: 10.5, y: 10.5, z: 0.0));

        // 3-d has higher dimensions than the rest of the collection, so its centroid clobbers
        // everything else in the collection.
        collection
            .0
            .push(Rect::new(coord!(10., 10., 10.), coord!(11., 11., 11.)).into());
        assert_eq!(collection.centroid().unwrap(), point!(x: 10.5, y: 10.5, z: 10.5));
    }
}
