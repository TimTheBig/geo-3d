use crate::{coord, CoordNum, CoordsIter, Polygon, Triangle};
use d_delaunay::delaunay_core::{Point as DDPoint, Tds, Vertex};
use geo_types::{Coord, MultiPoint};
use num_traits::NumCast;

/// Triangulate polygons using an [Delaunay triangulation](https://en.wikipedia.org/wiki/Delaunay_triangulation).
pub trait TriangulateDelaunay<T: CoordNum> {
    /// Collect a Triangulation of points.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_3d::{coord, polygon, Triangle, TriangulateDelaunay};
    ///
    /// let square_polygon = polygon![
    ///     (x: 0., y: 0., z: 0.), // SW
    ///     (x: 10., y: 0., z: 10.), // SE
    ///     (x: 10., y: 10., z: 10.), // NE
    ///     (x: 0., y: 10., z: 0.), // NW
    ///     (x: 0., y: 0., z: 0.), // SW
    /// ];
    ///
    /// let triangles = square_polygon.delaunay_triangles();
    ///
    /// assert_eq!(
    ///     vec![
    ///         Triangle(
    ///             coord! { x: 0., y: 10. }, // NW
    ///             coord! { x: 10., y: 10. }, // NE
    ///             coord! { x: 10., y: 0. }, // SE
    ///         ),
    ///         Triangle(
    ///             coord! { x: 10., y: 0. }, // SE
    ///             coord! { x: 0., y: 0. }, // SW
    ///             coord! { x: 0., y: 10. }, // NW
    ///         ),
    ///     ],
    ///     triangles,
    /// );
    /// ```
    fn delaunay_triangles(&self) -> Vec<Triangle<T>> {
        self.delaunay_triangles_iter().collect()
    }

    /// Triangulate a set of points.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_3d::{coord, polygon, Triangle, TriangulateDelaunay};
    ///
    /// let square_polygon = polygon![
    ///     (x: 0., y: 0., z: 0.), // SW
    ///     (x: 10., y: 0., z: 10.), // SE
    ///     (x: 10., y: 10., z: 10.), // NE
    ///     (x: 0., y: 10., z: 0.), // NW
    ///     (x: 0., y: 0., z: 0.), // SW
    /// ];
    ///
    /// let mut triangles_iter = square_polygon.delaunay_triangles_iter();
    ///
    /// assert_eq!(
    ///     Some(Triangle(
    ///             coord! { x: 0., y: 10. }, // NW
    ///             coord! { x: 10., y: 10. }, // NE
    ///             coord! { x: 10., y: 0. }, // SE
    ///     )),
    ///     triangles_iter.next(),
    /// );
    ///
    /// assert_eq!(
    ///     Some(Triangle(
    ///         coord! { x: 10., y: 0. }, // SE
    ///         coord! { x: 0., y: 0. }, // SW
    ///         coord! { x: 0., y: 10. }, // NW
    ///     )),
    ///     triangles_iter.next(),
    /// );
    ///
    /// assert!(triangles_iter.next().is_none());
    /// ```
    fn delaunay_triangles_iter(&self) -> impl Iterator<Item = Triangle<T>>;
}

impl<T: CoordNum + Default> TriangulateDelaunay<T> for Polygon<T> {
    fn delaunay_triangles_iter(&self) -> impl Iterator<Item = Triangle<T>> {
        let cells = Tds::<f64, usize, usize, 3>::new(polygon_to_delaunay_input(self))
            .bowyer_watson()
            .unwrap()
            .cells;

        cells.into_iter().map(|cell| {
            Triangle::new(
                vertex_to_coord::<T>(cell.1.vertices[0]),
                vertex_to_coord::<T>(cell.1.vertices[1]),
                vertex_to_coord::<T>(cell.1.vertices[2]),
            )
        })
    }
}

impl<T: CoordNum + Default> TriangulateDelaunay<T> for MultiPoint<T> {
    fn delaunay_triangles_iter(&self) -> impl Iterator<Item = Triangle<T>> {
        let cells = Tds::<f64, usize, usize, 3>::new(
            self.0.iter().map(|c| {
                DDPoint::<f64, 3>::new(
                    [NumCast::from(c.x()).unwrap(), NumCast::from(c.y()).unwrap(), NumCast::from(c.z()).unwrap()]
                )
            }).collect()
        ).bowyer_watson()
        .unwrap().cells;

        cells.into_iter().map(|cell| {
            Triangle::new(
                vertex_to_coord::<T>(cell.1.vertices[0]),
                vertex_to_coord::<T>(cell.1.vertices[1]),
                vertex_to_coord::<T>(cell.1.vertices[2]),
            )
        })
    }
}

fn vertex_to_coord<T: CoordNum + Default>(vertex: Vertex<f64, usize, 3>) -> Coord<T> {
    coord! {
        x: NumCast::from(vertex.point.coords[0]).unwrap(),
        y: NumCast::from(vertex.point.coords[1]).unwrap(),
        z: NumCast::from(vertex.point.coords[2]).unwrap(),
    }
}

fn polygon_to_delaunay_input<T: CoordNum + Default>(polygon: &Polygon<T>) -> Vec<DDPoint<f64, 3>> {
    let mut vertices = Vec::with_capacity(polygon.coords_count());
    // todo check if this is needed
    // let mut interior_indexes = Vec::with_capacity(polygon.interiors().len());
    debug_assert!(polygon.exterior().0.len() >= 4);

    flat_line_string_ddpoints(polygon.exterior(), &mut vertices);

    for interior in polygon.interiors() {
        debug_assert!(interior.0.len() >= 4);
        // interior_indexes.push(vertices.len() / 2);
        flat_line_string_ddpoints(interior, &mut vertices);
    }

    vertices
}

fn flat_line_string_ddpoints<T: CoordNum + Default>(
    line_string: &crate::LineString<T>,
    vertices: &mut Vec<DDPoint<f64, 3>>,
) {
    for coord in &line_string.0 {
        vertices.push(DDPoint::<f64, 3>::new([
            NumCast::from(coord.x).unwrap(),
            NumCast::from(coord.y).unwrap(),
            NumCast::from(coord.z).unwrap(),
        ]))
    }
}

#[cfg(test)]
mod test {
    use super::TriangulateDelaunay;
    use crate::{coord, polygon, Triangle};

    #[test]
    fn test_triangle() {
        let triangle_polygon = polygon![
            (x: 0., y: 0., z: 0.),
            (x: 10., y: 0., z: 0.),
            (x: 10., y: 10., z: 10.),
            (x: 0., y: 0., z: 0.),
        ];

        let triangles = triangle_polygon.delaunay_triangles();

        assert_eq!(
            &[Triangle(
                coord! { x: 10.0, y: 0.0, z: 0.0 },
                coord! { x: 0.0, y: 0.0, z: 0.0 },
                coord! { x: 10.0, y: 10.0, z: 10.0 },
            ),][..],
            triangles,
        );
    }

    #[test]
    fn test_square() {
        let square_polygon = polygon![
            (x: 0., y: 0., z: 0.),
            (x: 10., y: 0., z: 10.),
            (x: 10., y: 10., z: 10.),
            (x: 0., y: 10., z: 0.),
            (x: 0., y: 0., z: 0.),
        ];

        let mut triangles = square_polygon.delaunay_triangles();
        triangles.sort_by(|t1, t2| t1.1.x.partial_cmp(&t2.2.x).unwrap());

        assert_eq!(
            &[
                Triangle(
                    coord! { x: 10.0, y: 0.0, z: 10.0 },
                    coord! { x: 0.0, y: 0.0, z: 0.0 },
                    coord! { x: 0.0, y: 10.0, z: 0.0 },
                ),
                Triangle(
                    coord! { x: 0.0, y: 10.0, z: 0.0 },
                    coord! { x: 10.0, y: 10.0, z: 10.0 },
                    coord! { x: 10.0, y: 0.0, z: 10.0 },
                ),
            ][..],
            triangles,
        );
    }
}
