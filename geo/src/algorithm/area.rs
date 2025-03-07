use super::TriangulateDelaunay;
use crate::geometry::*;
use crate::CoordNum;
use std::iter::Sum;

pub(crate) fn twice_signed_ring_area<T>(linestring: &LineString<T>) -> T
where
    T: CoordNum,
{
    // LineString with less than 3 points is empty, or a
    // single point, or is not closed.
    if linestring.0.len() < 3 {
        return T::zero();
    }

    // Above test ensures the vector has at least 2 elements.
    // We check if linestring is closed, and return 0 otherwise.
    if linestring.0.first().unwrap() != linestring.0.last().unwrap() {
        return T::zero();
    }

    // Use a reasonable shift for the line-string coords
    // to avoid numerical-errors when summing the
    // determinants.
    //
    // Note: we can't use the `Centroid` trait as it
    // requires `T: Float` and in fact computes area in the
    // implementation. Another option is to use the average
    // of the coordinates, but it is not fool-proof to
    // divide by the length of the linestring (eg. a long
    // line-string with T = u8)
    let shift = linestring.0[0];

    let mut tmp = T::zero();
    for line in linestring.lines() {
        use crate::MapCoords;
        let line = line.map_coords(|c| c - shift);
        tmp = tmp + line.determinant();
    }

    tmp
}

/// Signed and unsigned planar area of a geometry.
///
/// # Examples
///
/// ```
/// use geo::polygon;
/// use geo::Area;
///
/// let mut polygon = polygon![
///     (x: 0., y: 0., z: 0.),
///     (x: 5., y: 0., z: 5.),
///     (x: 5., y: 6., z: 5.),
///     (x: 0., y: 6., z: 0.),
///     (x: 0., y: 0., z: 0.),
/// ];
///
/// assert_eq!(polygon.signed_area(), 30.);
/// assert_eq!(polygon.unsigned_area(), 30.);
///
/// polygon.exterior_mut(|line_string| {
///     line_string.0.reverse();
/// });
///
/// assert_eq!(polygon.signed_area(), -30.);
/// assert_eq!(polygon.unsigned_area(), 30.);
/// ```
pub trait Area<T>
where
    T: CoordNum,
{
    /// Calculate the signed area of a geometry.
    fn signed_area(&self) -> T;

    /// Calculate the unsigned area of a geometry.
    fn unsigned_area(&self) -> T;
}

// Calculation of simple (no interior holes) Polygon area
pub(crate) fn get_linestring_area<T>(linestring: &LineString<T>) -> T
where
    T: CoordNum,
{
    twice_signed_ring_area(linestring) / (T::one() + T::one())
}

impl<T> Area<T> for Point<T>
where
    T: CoordNum,
{
    fn signed_area(&self) -> T {
        T::zero()
    }

    fn unsigned_area(&self) -> T {
        T::zero()
    }
}

impl<T> Area<T> for LineString<T>
where
    T: CoordNum,
{
    fn signed_area(&self) -> T {
        T::zero()
    }

    fn unsigned_area(&self) -> T {
        T::zero()
    }
}

impl<T> Area<T> for Line<T>
where
    T: CoordNum,
{
    fn signed_area(&self) -> T {
        T::zero()
    }

    fn unsigned_area(&self) -> T {
        T::zero()
    }
}

/// **Note.** The implementation handles polygons whose
/// holes do not all have the same orientation. The sign of
/// the output is the same as that of the exterior shell.
impl<T> Area<T> for Polygon<T>
where
    T: CoordNum + Sum<T>,
{
    fn signed_area(&self) -> T {
        self.delaunay_triangles_iter().map(|tri| tri.signed_area()).sum()
    }

    fn unsigned_area(&self) -> T {
        self.delaunay_triangles_iter().map(|tri| tri.unsigned_area()).sum::<T>().abs()
    }
}

impl<T> Area<T> for MultiPoint<T>
where
    T: CoordNum,
{
    fn signed_area(&self) -> T {
        T::zero()
    }

    fn unsigned_area(&self) -> T {
        T::zero()
    }
}

impl<T> Area<T> for MultiLineString<T>
where
    T: CoordNum,
{
    fn signed_area(&self) -> T {
        T::zero()
    }

    fn unsigned_area(&self) -> T {
        T::zero()
    }
}

/// **Note.** The implementation is a straight-forward
/// summation of the signed areas of the individual
/// polygons. In particular, `unsigned_area` is not
/// necessarily the sum of the `unsigned_area` of the
/// constituent polygons unless they are all oriented the
/// same.
impl<T> Area<T> for MultiPolygon<T>
where
    T: CoordNum + Sum<T>,
{
    fn signed_area(&self) -> T {
        self.0
            .iter()
            .fold(T::zero(), |total, next| total + next.signed_area())
    }

    fn unsigned_area(&self) -> T {
        self.0
            .iter()
            .fold(T::zero(), |total, next| total + next.signed_area().abs())
    }
}

/// Because a `Rect` has no winding order, the area will always be positive.
impl<T: CoordNum> Area<T> for Rect<T> {
    /// Calculate the surface area of a `Rect`.\
    /// Because a `Rect` has no winding order, the area will always be positive.
    fn signed_area(&self) -> T {
        self.unsigned_area()
    }

    /// Calculate the unsigned surface area of a `Rect`.
    fn unsigned_area(&self) -> T {
        let dx = self.width();
        let dy = self.depth();
        let dz = self.height();

        // Surface area formula: 2(lw + lh + wh)
        (dx * dy + dx * dz + dy * dz) * (T::one() + T::one())
    }
}

impl<T: CoordNum> Area<T> for Triangle<T> {
    fn signed_area(&self) -> T {
        self.to_lines()
            .iter()
            .fold(T::zero(), |total, line| total + line.determinant())
            / (T::one() + T::one())
    }

    fn unsigned_area(&self) -> T {
        self.signed_area().abs()
    }
}

impl<T> Area<T> for Geometry<T>
where
    T: CoordNum + Sum<T>,
{
    crate::geometry_delegate_impl! {
        fn signed_area(&self) -> T;
        fn unsigned_area(&self) -> T;
    }
}

impl<T> Area<T> for GeometryCollection<T>
where
    T: CoordNum + Sum<T>,
{
    fn signed_area(&self) -> T {
        self.0
            .iter()
            .map(|g| g.signed_area())
            .fold(T::zero(), |acc, next| acc + next)
    }

    fn unsigned_area(&self) -> T {
        self.0
            .iter()
            .map(|g| g.unsigned_area())
            .fold(T::zero(), |acc, next| acc + next)
    }
}

#[cfg(test)]
mod test {
    use crate::Area;
    use crate::{coord, polygon, wkt, Line, MultiPolygon, Polygon, Rect, Triangle};

    // Area of the polygon
    #[test]
    fn area_empty_polygon_test() {
        let poly: Polygon<f32> = polygon![];
        assert_relative_eq!(poly.signed_area(), 0.);
    }

    #[test]
    fn area_one_point_polygon_test() {
        let poly = wkt! { POLYGON((1. 0. -1.)) };
        assert_relative_eq!(poly.signed_area(), 0.);
    }
    #[test]
    fn area_polygon_test() {
        let polygon = wkt! { POLYGON((0. 0. 0.,5. 0. 5.,5. 6. 5.,0. 6. 0.,0. 0. 0.)) };
        assert_relative_eq!(polygon.signed_area(), 30.);
    }
    #[test]
    fn area_polygon_numerical_stability() {
        let polygon = {
            use std::f64::consts::PI;
            const NUM_VERTICES: usize = 10;
            const ANGLE_INC: f64 = 2. * PI / NUM_VERTICES as f64;

            Polygon::new(
                (0..NUM_VERTICES)
                    .map(|i| {
                        let angle = i as f64 * ANGLE_INC;
                        coord! {
                            x: angle.cos(),
                            y: angle.sin(),
                            z: angle.cos(),
                        }
                    })
                    .collect::<Vec<_>>()
                    .into(),
                vec![],
            )
        };

        let area = polygon.signed_area();

        let shift = coord! { x: 1.5e8, y: 1.5e8, z: 1.5e8 };

        use crate::map_coords::MapCoords;
        let polygon = polygon.map_coords(|c| c + shift);

        let new_area = polygon.signed_area();
        let err = (area - new_area).abs() / area;

        assert!(err < 1e-2);
    }
    #[test]
    fn rectangle_test() {
        let rect1: Rect<f32> = Rect::new(coord! { x: 10., y: 30., z: 50. }, coord! { x: 20., y: 40., z: 60. });
        assert_eq!(rect1.signed_area(), 600.);

        let rect2: Rect = Rect::new(coord! { x: 10., y: 30., z: 0. }, coord! { x: 20., y: 40., z: 1. });
        assert_eq!(rect2.signed_area(), 240.);

        assert_eq!(rect2.signed_area(), rect2.unsigned_area());
    }
    #[test]
    fn area_polygon_inner_test() {
        let poly = polygon![
            exterior: [
                (x: 0., y: 0., z: 0.),
                (x: 10., y: 0., z: 10.),
                (x: 10., y: 10., z: 10.),
                (x: 0., y: 10., z: 0.),
                (x: 0., y: 0., z: 0.)
            ],
            interiors: [
                [
                    (x: 1., y: 1., z: 1.),
                    (x: 2., y: 1., z: 2.),
                    (x: 2., y: 2., z: 2.),
                    (x: 1., y: 2., z: 1.),
                    (x: 1., y: 1., z: 1.),
                ],
                [
                    (x: 5., y: 5., z: 5.),
                    (x: 6., y: 5., z: 6.),
                    (x: 6., y: 6., z: 6.),
                    (x: 5., y: 6., z: 5.),
                    (x: 5., y: 5., z: 5.)
                ],
            ],
        ];
        assert_relative_eq!(poly.signed_area(), 98.);
    }
    #[test]
    fn area_multipolygon_test() {
        let poly0 = polygon![
            (x: 0., y: 0., z: 0.),
            (x: 10., y: 0., z: -10.),
            (x: 10., y: 10., z: 10.),
            (x: 0., y: 10., z: 0.),
            (x: 0., y: 0., z: 0.)
        ];
        let poly1 = polygon![
            (x: 1., y: 1., z: 1.),
            (x: 2., y: 1., z: 2.),
            (x: 2., y: 2., z: 2.),
            (x: 1., y: 2., z: 1.),
            (x: 1., y: 1., z: 1.)
        ];
        let poly2 = polygon![
            (x: 5., y: 5., z: 5.),
            (x: 6., y: 5., z: 6.),
            (x: 6., y: 6., z: 6.),
            (x: 5., y: 6., z: 5.),
            (x: 5., y: 5., z: 5.)
        ];
        let mpoly = MultiPolygon::new(vec![poly0, poly1, poly2]);
        assert_relative_eq!(mpoly.signed_area(), 102.);
        assert_relative_eq!(mpoly.signed_area(), 102.);
    }
    #[test]
    fn area_line_test() {
        let line1 = Line::new(coord! { x: 0.0, y: 0.0, z: 0.0 }, coord! { x: 1.0, y: 1.0, z: 1.0 });
        assert_relative_eq!(line1.signed_area(), 0.);
    }

    #[test]
    fn area_triangle_test() {
        // One on each axis
        let triangle = Triangle::new(
            coord! { x: 0.0, y: 0.0, z: 1.0 },
            coord! { x: 1.0, y: 0.0, z: 0.0 },
            coord! { x: 0.0, y: 1.0, z: 0.0 },
        );
        assert_relative_eq!(triangle.signed_area(), 1.5);
        assert_eq!(triangle.signed_area(), triangle.unsigned_area());

        let triangle = Triangle::new(
            coord! { x: 0.0, y: 0.0, z: 0.0 },
            coord! { x: 1.0, y: 0.0, z: 1.0 },
            coord! { x: 0.0, y: 1.0, z: 0.0 },
        );
        assert_relative_eq!(triangle.signed_area(), 0.707, epsilon = 0.00011);

        let triangle = Triangle::new(
            coord! { x: 0.0, y: 0.0, z: -1.0 },
            coord! { x: 0.0, y: -1.0, z: 0.0 },
            coord! { x: 1.0, y: 0.0, z: -1.0 },
        );
        assert_relative_eq!(triangle.signed_area(), -0.5);
    }

    #[test]
    fn area_multi_polygon_area_reversed() {
        let polygon_cw: Polygon<f32> = polygon![
            coord! { x: 0.0, y: 0.0, z: 0.0 },
            coord! { x: 0.0, y: 1.0, z: 0.0 },
            coord! { x: 1.0, y: 1.0, z: 1.0 },
            coord! { x: 1.0, y: 0.0, z: 1.0 },
            coord! { x: 0.0, y: 0.0, z: 0.0 },
        ];
        let polygon_ccw: Polygon<f32> = polygon![
            coord! { x: 0.0, y: 0.0, z: 0.0 },
            coord! { x: 1.0, y: 0.0, z: 1.0 },
            coord! { x: 1.0, y: 1.0, z: 1.0 },
            coord! { x: 0.0, y: 1.0, z: 0.0 },
            coord! { x: 0.0, y: 0.0, z: 0.0 },
        ];
        let polygon_area = polygon_cw.unsigned_area();

        let multi_polygon = MultiPolygon::new(vec![polygon_cw, polygon_ccw]);

        assert_eq!(polygon_area * 2., multi_polygon.unsigned_area());
    }

    #[test]
    fn area_north_america_cutout() {
        let poly = polygon![
            exterior: [
                (x: -102.902861858977, y: 31.6943450891131, z: 0.0),
                (x: -102.917375513247, y: 31.6990175356827, z: 0.0),
                (x: -102.917887344527, y: 31.7044889522597, z: 0.0),
                (x: -102.938892711173, y: 31.7032871894594, z: 0.0),
                (x: -102.939919687305, y: 31.7142296141915, z: 0.0),
                (x: -102.946922353444, y: 31.713828170995, z: 0.0),
                (x: -102.954642979004, y: 31.7210594956594, z: 0.0),
                (x: -102.960927457803, y: 31.7130240707676, z: 0.0),
                (x: -102.967929895872, y: 31.7126214137469, z: 0.0),
                (x: -102.966383373178, y: 31.6962079209847, z: 0.0),
                (x: -102.973384192133, y: 31.6958049292994, z: 0.0),
                (x: -102.97390013779, y: 31.701276160078, z: 0.0),
                (x: -102.980901394769, y: 31.7008727405409, z: 0.0),
                (x: -102.987902575456, y: 31.7004689164622, z: 0.0),
                (x: -102.986878877087, y: 31.7127206248263, z: 0.0),
                (x: -102.976474089689, y: 31.7054378797983, z: 0.0),
                (x: -102.975448432121, y: 31.7176893134691, z: 0.0),
                (x: -102.96619351228, y: 31.7237224912303, z: 0.0),
                (x: -102.976481009643, y: 31.7286309669534, z: 0.0),
                (x: -102.976997412845, y: 31.7341016591658, z: 0.0),
                (x: -102.978030448215, y: 31.7450427747035, z: 0.0),
                (x: -102.985035821671, y: 31.7446391683265, z: 0.0),
                (x: -102.985552968771, y: 31.7501095683386, z: 0.0),
                (x: -102.992558780682, y: 31.7497055338313, z: 0.0),
                (x: -102.993594334215, y: 31.7606460184322, z: 0.0),
                (x: -102.973746840657, y: 31.7546100958509, z: 0.0),
                (x: -102.966082339116, y: 31.767730116605, z: 0.0),
                (x: -102.959074676589, y: 31.768132602064, z: 0.0),
                (x: -102.95206693787, y: 31.7685346826851, z: 0.0),
                (x: -102.953096767614, y: 31.7794749110023, z: 0.0),
                (x: -102.953611796704, y: 31.7849448911322, z: 0.0),
                (x: -102.952629078076, y: 31.7996518517642, z: 0.0),
                (x: -102.948661251495, y: 31.8072257578725, z: 0.0),
                (x: -102.934638176282, y: 31.8080282207231, z: 0.0),
                (x: -102.927626524626, y: 31.8084288446215, z: 0.0),
                (x: -102.927113253813, y: 31.8029591283411, z: 0.0),
                (x: -102.920102042027, y: 31.8033593239799, z: 0.0),
                (x: -102.919076759513, y: 31.792419577395, z: 0.0),
                (x: -102.912066503301, y: 31.7928193216213, z: 0.0),
                (x: -102.911554491357, y: 31.7873492912889, z: 0.0),
                (x: -102.904544675025, y: 31.7877486073783, z: 0.0),
                (x: -102.904033254331, y: 31.7822784646103, z: 0.0),
                (x: -102.903521909259, y: 31.7768082325431, z: 0.0),
                (x: -102.895800463718, y: 31.7695748336589, z: 0.0),
                (x: -102.889504111843, y: 31.7776055573633, z: 0.0),
                (x: -102.882495099915, y: 31.7780036124077, z: 0.0),
                (x: -102.868476849997, y: 31.7787985077398, z: 0.0),
                (x: -102.866950998738, y: 31.7623869292283, z: 0.0),
                (x: -102.873958615171, y: 31.7619897531194, z: 0.0),
                (x: -102.87888647278, y: 31.7688910039026, z: 0.0),
                (x: -102.879947237315, y: 31.750650764952, z: 0.0),
                (x: -102.886953672823, y: 31.750252825268, z: 0.0),
                (x: -102.89396003296, y: 31.7498544807869, z: 0.0),
                (x: -102.892939355062, y: 31.7389128078806, z: 0.0),
                (x: -102.913954892669, y: 31.7377154844276, z: 0.0),
                (x: -102.913443122277, y: 31.7322445829725, z: 0.0),
                (x: -102.912931427507, y: 31.7267735918962, z: 0.0),
                (x: -102.911908264767, y: 31.7158313407426, z: 0.0),
                (x: -102.904905220014, y: 31.7162307607961, z: 0.0),
                (x: -102.904394266551, y: 31.7107594775392, z: 0.0),
                (x: -102.903372586049, y: 31.6998166417321, z: 0.0),
                (x: -102.902861858977, y: 31.6943450891131, z: 0.0),
            ],
            interiors: [
                [
                    (x: -102.916514879554, y: 31.7650686485918, z: 0.0),
                    (x: -102.921022256876, y: 31.7770831833398, z: 0.0),
                    (x: -102.93367363719, y: 31.771184865332, z: 0.0),
                    (x: -102.916514879554, y: 31.7650686485918, z: 0.0),
                ],
                [
                    (x: -102.935483140202, y: 31.7419852607081, z: 0.0),
                    (x: -102.932452314332, y: 31.7328567234689, z: 0.0),
                    (x: -102.918345099146, y: 31.7326099897391, z: 0.0),
                    (x: -102.925566322952, y: 31.7552505533503, z: 0.0),
                    (x: -102.928990700436, y: 31.747856686604, z: 0.0),
                    (x: -102.935996606762, y: 31.7474559134477, z: 0.0),
                    (x: -102.939021176592, y: 31.7539885279379, z: 0.0),
                    (x: -102.944714388971, y: 31.7488395547293, z: 0.0),
                    (x: -102.935996606762, y: 31.7474559134477, z: 0.0),
                    (x: -102.935483140202, y: 31.7419852607081, z: 0.0),
                ],
                [
                    (x: -102.956498858767, y: 31.7407805824758, z: 0.0),
                    (x: -102.960959476367, y: 31.7475080456347, z: 0.0),
                    (x: -102.972817445204, y: 31.742072061889, z: 0.0),
                    (x: -102.956498858767, y: 31.7407805824758, z: 0.0),
                ]
            ],
        ];
        // Value from shapely
        assert_relative_eq!(
            poly.unsigned_area(),
            0.006547948219252177,
            max_relative = 0.0001
        );
    }
}
