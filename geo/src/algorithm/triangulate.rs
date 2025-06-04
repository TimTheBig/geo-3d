use crate::{coord, CoordNum, Polygon, Triangle};
// use geo_types::{Coord, MultiPoint, MultiPolygon};
use geo_types::{Coord, MultiPolygon};
use num_traits::NumCast;

use earclip_rs::earclip;

/// Triangulate polygons using  triangulation](https://en.wikipedia.org/wiki/Delaunay_triangulation).
pub trait Triangulate<T: CoordNum> {
    /// Collect a Triangulation of points.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_3d::{coord, polygon, Triangle, Triangulate};
    ///
    /// let square_polygon = polygon![
    ///     (x: 0., y: 0., z: 0.), // SW
    ///     (x: 10., y: 0., z: 10.), // SE
    ///     (x: 10., y: 10., z: 10.), // NE
    ///     (x: 0., y: 10., z: 0.), // NW
    ///     (x: 0., y: 0., z: 0.), // SW
    /// ];
    ///
    /// let triangles = square_polygon.triangles();
    ///
    /// assert_eq!(
    ///     vec![
    ///         Triangle(
    ///             coord! { x: 0., y: 10., z: 0. }, // NW
    ///             coord! { x: 10., y: 10., z: 10. }, // NE
    ///             coord! { x: 10., y: 0., z: 10. }, // SE
    ///         ),
    ///         Triangle(
    ///             coord! { x: 10., y: 0., z: 10. }, // SE
    ///             coord! { x: 0., y: 0., z: 0. }, // SW
    ///             coord! { x: 0., y: 10., z: 0. }, // NW
    ///         ),
    ///     ],
    ///     triangles,
    /// );
    /// ```
    fn triangles(&self) -> Vec<Triangle<T>>;
}

// todo check if duplecate (closing) last coord is needed
impl<T: CoordNum> Triangulate<T> for Polygon<T> {
    /// Returns an empty vec if polygon has less then three coords
    fn triangles(&self) -> Vec<Triangle<T>> {
        if self.rings().flatten().count() < 3 {
            // empty or 1/2 coords
            return Vec::new();
        }
        if self.exterior().0.len() == 3
        && (self.interiors().is_empty() || self.interiors().iter().all(|ls| ls.0.is_empty())) {
                return vec![Triangle::new(self.exterior()[0], self.exterior()[1], self.exterior()[2])];
        }

        let mut polygon = Vec::with_capacity(self.interiors().len() + 1);
        polygon.push(self.exterior().0.clone());
        self.interiors().iter().for_each(|ls| polygon.push(ls.0.clone()));

        let (vertices, indices) = earclip(
            &polygon,
            None, None
        );

        let coords: Vec<Coord<T>> = vertices.chunks_exact(3)
            .map(|chunk| {
                coord!(
                    x: NumCast::from(chunk[0]).unwrap_or(T::nan()),
                    y: NumCast::from(chunk[1]).unwrap_or(T::nan()),
                    z: NumCast::from(chunk[2]).unwrap_or(T::nan()),
                )
            }).collect();
        let tris = indices.chunks_exact(3)
            .map(|idxs| {
                // todo check if idx is start or end of coord
                // flip winding
                Triangle::new(
                    coords[idxs[2]],
                    coords[idxs[1]],
                    coords[idxs[0]]
                )
            }).rev()
            .collect();

        tris
    }
}

// todo impl for rect use `.to_points()`

impl<T: CoordNum> Triangulate<T> for MultiPolygon<T> {
    fn triangles(&self) -> Vec<Triangle<T>> {
        self.0.iter().map(|poly| poly.triangles()).flatten().collect()
    }
}

// impl<T: CoordNum + Default> Triangulate<T> for Polygon<T> {
//     fn triangles_iter(&self) -> impl Iterator<Item = Triangle<T>> {
//         if self.is_empty() {
//             return std::iter::empty::<Triangle<T>>();
//         }

//         let bounding_rect = self.bounding_rect().expect("is_empty check above");
//         let mut container = voro_rs::container::ContainerStd::new(
//             coord!(bounding_rect.min().x, bounding_rect.min().y, bounding_rect.min().z),
//             coord!(bounding_rect.max().x, bounding_rect.max().y, bounding_rect.max().z),
//             [2, 2, 2],
//             [false, false, false],
//         );
//         for (i, coord) in
//         container.put(n, xyz, r);
//     }
// }

// impl<T: CoordNum + Default> Triangulate<T> for MultiPoint<T> {
//     fn triangles_iter(&self) -> impl Iterator<Item = Triangle<T>> {
//         let cells = Tds::<f64, usize, usize, 3>::new(
//             self.0.iter().map(|c| {
//                 DDPoint::<f64, 3>::new(
//                     [NumCast::from(c.x()).unwrap(), NumCast::from(c.y()).unwrap(), NumCast::from(c.z()).unwrap()]
//                 )
//             }).collect()
//         ).bowyer_watson()
//         .unwrap().cells;

//         cells.into_iter().map(|cell| {
//             Triangle::new(
//                 vertex_to_coord::<T>(cell.1.vertices[0]),
//                 vertex_to_coord::<T>(cell.1.vertices[1]),
//                 vertex_to_coord::<T>(cell.1.vertices[2]),
//             )
//         })
//     }
// }

// fn vertex_to_coord<T: CoordNum + Default>(vertex: Vertex<f64, usize, 3>) -> Coord<T> {
//     coord! {
//         x: NumCast::from(vertex.point.coords[0]).unwrap(),
//         y: NumCast::from(vertex.point.coords[1]).unwrap(),
//         z: NumCast::from(vertex.point.coords[2]).unwrap(),
//     }
// }

// fn polygon_to_delaunay_input<T: CoordNum + Default>(polygon: &Polygon<T>) -> Vec<DDPoint<f64, 3>> {
//     let mut vertices = Vec::with_capacity(polygon.coords_count());
//     // todo check if this is needed
//     // let mut interior_indexes = Vec::with_capacity(polygon.interiors().len());
//     debug_assert!(polygon.exterior().0.len() >= 4);

//     flat_line_string_ddpoints(polygon.exterior(), &mut vertices);

//     for interior in polygon.interiors() {
//         debug_assert!(interior.0.len() >= 4);
//         // interior_indexes.push(vertices.len() / 2);
//         flat_line_string_ddpoints(interior, &mut vertices);
//     }

//     vertices
// }

// fn flat_line_string_ddpoints<T: CoordNum + Default>(
//     line_string: &crate::LineString<T>,
//     vertices: &mut Vec<DDPoint<f64, 3>>,
// ) {
//     for coord in &line_string.0 {
//         vertices.push(DDPoint::<f64, 3>::new([
//             NumCast::from(coord.x).unwrap(),
//             NumCast::from(coord.y).unwrap(),
//             NumCast::from(coord.z).unwrap(),
//         ]));
//     }
// }

#[cfg(test)]
mod test {
    use super::Triangulate;
    use crate::{coord, polygon, Triangle};

    #[test]
    fn test_triangle() {
        let triangle_polygon = polygon![
            (x: 0., y: 0., z: 0.),
            (x: 10., y: 0., z: 0.),
            (x: 10., y: 10., z: 10.),
            (x: 0., y: 0., z: 0.),
        ];

        let triangles = triangle_polygon.triangles();

        assert_eq!(
            &[Triangle(
                coord!(x: 10.0, y: 0.0, z: 0.0),
                coord!(x: 0.0, y: 0.0, z: 0.0),
                coord!(x: 10.0, y: 10.0, z: 10.0),
            )][..],
            &triangles,
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

        let mut triangles = square_polygon.triangles();
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
