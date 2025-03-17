use crate::algorithm::{AffineOps, AffineTransform, BoundingRect, Centroid};
use crate::geometry::*;
use crate::CoordNum;

/// Rotate a geometry around a point by an angle, in degrees.
///
/// Positive angles are counter-clockwise, and negative angles are clockwise rotations.
///
/// ## Performance
///
/// If you will be performing multiple transformations, like [`Scale`](crate::Scale),
/// [`Skew`](crate::Skew), [`Translate`](crate::Translate) or [`Rotate`], it is more
/// efficient to compose the transformations and apply them as a single operation using the
/// [`AffineOps`] trait.
pub trait Rotate<T: CoordNum> {
    /// Rotate a geometry around its [centroid](Centroid) by an angle, in degrees
    ///
    /// Positive angles are counter-clockwise, and negative angles are clockwise rotations.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo::Rotate;
    /// use geo::line_string;
    /// use approx::assert_relative_eq;
    ///
    /// let line_string = line_string![
    ///     (x: 0.0, y: 0.0),
    ///     (x: 5.0, y: 5.0),
    ///     (x: 10.0, y: 10.0),
    /// ];
    ///
    /// let rotated = line_string.rotate_around_centroid(-45.0);
    ///
    /// let expected = line_string![
    ///     (x: -2.071067811865475, y: 5.0),
    ///     (x: 5.0, y: 5.0),
    ///     (x: 12.071067811865476, y: 5.0),
    /// ];
    ///
    /// assert_relative_eq!(expected, rotated);
    /// ```
    #[must_use]
    fn rotate_around_centroid(&self, degrees: T) -> Self;

    /// Mutable version of [`Self::rotate_around_centroid`]
    fn rotate_around_centroid_mut(&mut self, degrees: T);

    /// Rotate a geometry around the center of its [bounding box](BoundingRect) by an angle, in
    /// degrees.
    ///
    /// Positive angles are counter-clockwise, and negative angles are clockwise rotations.
    ///
    #[must_use]
    fn rotate_around_center(&self, degrees: T) -> Self;

    /// Mutable version of [`Self::rotate_around_center`]
    fn rotate_around_center_mut(&mut self, degrees: T);

    /// Rotate a Geometry around an arbitrary point by an angle, given in degrees
    ///
    /// Positive angles are counter-clockwise, and negative angles are clockwise rotations.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo::Rotate;
    /// use geo::{line_string, point};
    ///
    /// let ls = line_string![
    ///     (x: 0.0, y: 0.0),
    ///     (x: 5.0, y: 5.0),
    ///     (x: 10.0, y: 10.0)
    /// ];
    ///
    /// let rotated = ls.rotate_around_point(
    ///     -45.0,
    ///     point!(x: 10.0, y: 0.0),
    /// );
    ///
    /// approx::assert_relative_eq!(rotated, line_string![
    ///     (x: 2.9289321881345245, y: 7.071067811865475),
    ///     (x: 10.0, y: 7.0710678118654755),
    ///     (x: 17.071067811865476, y: 7.0710678118654755)
    /// ]);
    /// ```
    #[must_use]
    fn rotate_around_point(&self, degrees: T, point: Point<T>) -> Self;

    /// Mutable version of [`Self::rotate_around_point`]
    fn rotate_around_point_mut(&mut self, degrees: T, point: Point<T>);
}

impl<G, IP, IR, T> Rotate<T> for G
where
    T: CoordNum,
    IP: Into<Option<Point<T>>>,
    IR: Into<Option<Rect<T>>>,
    G: Clone + Centroid<Output = IP> + BoundingRect<T, Output = IR> + AffineOps<T>,
{
    fn rotate_around_centroid(&self, degrees: T) -> Self {
        let point = match self.centroid().into() {
            Some(coord) => coord,
            // geometry was empty, so there's nothing to rotate
            None => return self.clone(),
        };
        Rotate::rotate_around_point(self, degrees, point)
    }

    fn rotate_around_centroid_mut(&mut self, degrees: T) {
        let point = match self.centroid().into() {
            Some(coord) => coord,
            // geometry was empty, so there's nothing to rotate
            None => return,
        };
        self.rotate_around_point_mut(degrees, point)
    }

    fn rotate_around_center(&self, degrees: T) -> Self {
        let point = match self.bounding_rect().into() {
            Some(rect) => Point(rect.center()),
            // geometry was empty, so there's nothing to rotate
            None => return self.clone(),
        };
        Rotate::rotate_around_point(self, degrees, point)
    }

    fn rotate_around_center_mut(&mut self, degrees: T) {
        let point = match self.bounding_rect().into() {
            Some(rect) => Point(rect.center()),
            // geometry was empty, so there's nothing to rotate
            None => return,
        };
        self.rotate_around_point_mut(degrees, point)
    }

    fn rotate_around_point(&self, degrees: T, point: Point<T>) -> Self {
        let transform = AffineTransform::rotate_xy(degrees, point);
        eprintln!("{}", transform);
        self.affine_transform(&transform)
    }

    fn rotate_around_point_mut(&mut self, degrees: T, point: Point<T>) {
        let transform = AffineTransform::rotate_xy(degrees, point);
        self.affine_transform_mut(&transform)
    }
}

#[cfg(test)]
mod test {
    use crate::algorithm::Rotate;
    use crate::geometry::*;
    use crate::{line_string, point, polygon};
    use approx::assert_relative_eq;
    use geo_types::coord;

    #[test]
    fn test_rotate_around_point() {
        let p = point!(x: 1.0, y: 5.0, z: 2.0);
        // rotating a point about itself leaves it unchanged (including its z)
        let rotated = p.rotate_around_centroid(30.0);
        assert_eq!(rotated, Point::new(1.0, 5.0, 2.0));
    }

    #[test]
    fn test_rotate_points() {
        let pt = point!(x: 1.0, y: 5.0, z: 2.0);
        let rotated_center = pt.rotate_around_center(30.);
        let rotated_centroid = pt.rotate_around_centroid(30.);
        // A point rotated about any point (including itself) should equal itself.
        assert_eq!(pt, rotated_center);
        assert_eq!(pt, rotated_centroid);
    }

    #[test]
    fn test_rotate_multipoints() {
        let multi_points = MultiPoint::new(vec![
            point!(x: 0., y: 0., z: 0.),
            point!(x: 1., y: 1., z: 1.),
            point!(x: 2., y: 1., z: -1.),
        ]);

        // When rotating about the centroid the rotation is applied only to x and y,
        // so the z values remain unchanged.
        let expected_for_centroid = MultiPoint::new(vec![
            point!(x: 0.7642977396044841, y: -0.5118446353109125, z: 0.0),
            point!(x: 0.7642977396044842, y: 0.9023689270621824, z: 1.0),
            point!(x: 1.471404520791032,  y: 1.60947570824873,  z: -1.0),
        ]);
        assert_relative_eq!(
            multi_points.rotate_around_centroid(45.),
            expected_for_centroid
        );

        let expected_for_center = MultiPoint::new(vec![
            point!(x: 0.6464466094067262, y: -0.5606601717798212, z: 0.0),
            point!(x: 0.6464466094067263, y: 0.8535533905932737,  z: 1.0),
            point!(x: 1.353553390593274,  y: 1.560660171779821,   z: -1.0),
        ]);
        assert_relative_eq!(
            multi_points.rotate_around_center(45.),
            expected_for_center,
            epsilon = 1e-15
        );
    }

    #[test]
    fn test_rotate_linestring() {
        let linestring = line_string![
            (x: 0.0, y: 0.0, z: 0.0),
            (x: 5.0, y: 5.0, z: 1.0),
            (x: 5.0, y: 10.0, z: 2.0)
        ];

        // For a rotation about the centroid the x and y values change but the z values are preserved.
        let rotated_around_centroid = linestring.rotate_around_centroid(-45.0);
        assert_relative_eq!(
            rotated_around_centroid,
            line_string![
                (x: -2.196699141100894, y: 3.838834764831844, z: 0.0),
                (x: 4.874368670764582,  y: 3.838834764831844, z: 1.0),
                (x: 8.40990257669732,   y: 7.374368670764582, z: 2.0)
            ]
        );

        let rotated_around_center = linestring.rotate_around_center(-45.0);
        assert_relative_eq!(
            rotated_around_center,
            line_string![
                (x: -2.803300858899106, y: 3.232233047033631, z: 0.0),
                (x: 4.267766952966369,  y: 3.232233047033632, z: 1.0),
                (x: 7.803300858899107,  y: 6.767766952966369, z: 2.0)
            ],
            epsilon = 1e-12
        );
    }

    #[test]
    fn test_rotate_polygon() {
        // The polygon’s vertices now include a z value.
        let poly1 = polygon![
            (x: 5., y: 1., z: 10.),
            (x: 4., y: 2., z: 11.),
            (x: 4., y: 3., z: 12.),
            (x: 5., y: 4., z: 13.),
            (x: 6., y: 4., z: 14.),
            (x: 7., y: 3., z: 15.),
            (x: 7., y: 2., z: 16.),
            (x: 6., y: 1., z: 17.),
            (x: 5., y: 1., z: 10.)
        ];
        let rotated = poly1.rotate_around_centroid(-15.0);
        let correct = polygon![
            (x: 4.6288085192016855, y: 1.1805207831176578, z: 10.0),
            (x: 3.921701738015137,  y: 2.405265654509247,  z: 11.0),
            (x: 4.180520783117659,  y: 3.3711914807983154, z: 12.0),
            (x: 5.405265654509247,  y: 4.0782982619848624, z: 13.0),
            (x: 6.371191480798316,  y: 3.819479216882342,  z: 14.0),
            (x: 7.0782982619848624, y: 2.594734345490753,  z: 15.0),
            (x: 6.819479216882343,  y: 1.6288085192016848, z: 16.0),
            (x: 5.594734345490753,  y: 0.9217017380151372, z: 17.0),
            (x: 4.6288085192016855, y: 1.1805207831176578, z: 10.0)
        ];
        assert_eq!(rotated, correct);
    }

    #[test]
    fn test_rotate_polygon_holes() {
        let poly1 = polygon![
            exterior: [
                (x: 5.0, y: 1.0,  z: 10.0),
                (x: 4.0, y: 2.0,  z: 11.0),
                (x: 4.0, y: 3.0,  z: 12.0),
                (x: 5.0, y: 4.0,  z: 13.0),
                (x: 6.0, y: 4.0,  z: 14.0),
                (x: 7.0, y: 3.0,  z: 15.0),
                (x: 7.0, y: 2.0,  z: 16.0),
                (x: 6.0, y: 1.0,  z: 17.0),
                (x: 5.0, y: 1.0,  z: 10.0)
            ],
            interiors: [
                [
                    (x: 5.0, y: 1.3, z: 20.0),
                    (x: 5.5, y: 2.0, z: 21.0),
                    (x: 6.0, y: 1.3, z: 22.0),
                    (x: 5.0, y: 1.3, z: 20.0),
                ],
                [
                    (x: 5.0, y: 2.3, z: 30.0),
                    (x: 5.5, y: 3.0, z: 31.0),
                    (x: 6.0, y: 2.3, z: 32.0),
                    (x: 5.0, y: 2.3, z: 30.0),
                ],
            ],
        ];

        // Expected result when rotating around the center (of the bounding box).
        let center_expected = polygon![
            exterior: [
                (x: 4.628808519201685, y: 1.180520783117658, z: 10.0),
                (x: 3.921701738015137, y: 2.405265654509247, z: 11.0),
                (x: 4.180520783117659, y: 3.371191480798315, z: 12.0),
                (x: 5.405265654509247, y: 4.078298261984862, z: 13.0),
                (x: 6.371191480798316, y: 3.819479216882342, z: 14.0),
                (x: 7.078298261984862, y: 2.594734345490753, z: 15.0),
                (x: 6.819479216882343, y: 1.628808519201685, z: 16.0),
                (x: 5.594734345490753, y: 0.9217017380151372, z: 17.0),
                (x: 4.628808519201685, y: 1.180520783117658, z: 10.0),
            ],
            interiors: [
                [
                    (x: 4.706454232732442, y: 1.470298531004379, z: 20.0),
                    (x: 5.37059047744874,  y: 2.017037086855466, z: 21.0),
                    (x: 5.67238005902151,  y: 1.211479485901858, z: 22.0),
                    (x: 4.706454232732442, y: 1.470298531004379, z: 20.0),
                ],
                [
                    (x: 4.965273277834962, y: 2.436224357293447, z: 30.0),
                    (x: 5.62940952255126,  y: 2.982962913144534, z: 31.0),
                    (x: 5.931199104124032, y: 2.177405312190926, z: 32.0),
                    (x: 4.965273277834962, y: 2.436224357293447, z: 30.0),
                ],
            ],
        ];

        let rotated_around_center = poly1.rotate_around_center(-15.);
        assert_relative_eq!(rotated_around_center, center_expected, epsilon = 1e-12);

        // Expected result when rotating around the centroid.
        let centroid_expected = polygon![
            exterior: [
                (x: 4.615388272418591, y: 1.182287592124891, z: 10.0),
                (x: 3.908281491232044, y: 2.40703246351648,  z: 11.0),
                (x: 4.167100536334565, y: 3.372958289805549, z: 12.0),
                (x: 5.391845407726153, y: 4.080065070992097, z: 13.0),
                (x: 6.357771234015222, y: 3.821246025889576, z: 14.0),
                (x: 7.064878015201769, y: 2.596501154497987, z: 15.0),
                (x: 6.806058970099248, y: 1.630575328208918, z: 16.0),
                (x: 5.58131409870766,  y: 0.9234685470223708, z: 17.0),
                (x: 4.615388272418591, y: 1.182287592124891, z: 10.0),
            ],
            interiors: [
                [
                    (x: 4.693033985949348, y: 1.472065340011612, z: 20.0),
                    (x: 5.357170230665646, y: 2.0188038958627,   z: 21.0),
                    (x: 5.658959812238415, y: 1.213246294909091, z: 22.0),
                    (x: 4.693033985949348, y: 1.472065340011612, z: 20.0),
                ],
                [
                    (x: 4.951853031051868, y: 2.43799116630068,  z: 30.0),
                    (x: 5.615989275768166, y: 2.984729722151768, z: 31.0),
                    (x: 5.917778857340937, y: 2.179172121198159, z: 32.0),
                    (x: 4.951853031051868, y: 2.43799116630068,  z: 30.0),
                ],
            ],
        ];
        let rotated_around_centroid = poly1.rotate_around_centroid(-15.);
        assert_relative_eq!(rotated_around_centroid, centroid_expected, epsilon = 1e-12);
    }

    #[test]
    fn test_rotate_around_point_arbitrary() {
        let p = Point::new(5.0, 10.0, 3.0);
        let rotated = p.rotate_around_point(-45., Point::new(10.0, 34.0, 3.0));
        // The rotation is computed in the x–y plane so the z value (3.0) is preserved.
        assert_relative_eq!(
            rotated,
            Point::new(-10.506096654409877, 20.564971157455595, 3.0),
            epsilon = 1e-14
        );
    }

    #[test]
    fn test_rotate_line() {
        let line0 = Line::from([(0., 0., 0.), (0., 2., 0.)]);
        let line1 = Line::from([(1., 1., 0.), (-1., 1., 0.)]);
        assert_relative_eq!(line0.rotate_around_centroid(90.0), line1);
        assert_relative_eq!(line0.rotate_around_center(90.0), line1);
    }

    #[test]
    fn test_rotate_multi_line_string() {
        let ls1 = line_string![
            (x: 0., y: 0., z: 0.),
            (x: 1., y: 1., z: 1.),
            (x: 4., y: 1., z: 2.)
        ];
        let ls2 = line_string![
            (x: 10., y: 10., z: 3.),
            (x: 20., y: 20., z: 4.),
            (x: 40., y: 20., z: 5.)
        ];
        let multi_line_string: MultiLineString = MultiLineString::new(vec![ls1, ls2]);

        let expected_around_centroid = MultiLineString::new(vec![
            line_string![
                (x: -5.062519283392216, y: 19.72288595632566, z: 0.0),
                (x: -3.648305721019121, y: 19.72288595632566, z: 1.0),
                (x: -1.526985377459479, y: 17.60156561276602, z: 2.0)
            ],
            line_string![
                (x: 9.079616340338735, y: 19.72288595632566, z: 3.0),
                (x: 23.22175196406969, y: 19.72288595632566, z: 4.0),
                (x: 37.36388758780063, y: 5.580750332594715, z: 5.0)
            ],
        ]);
        assert_relative_eq!(
            multi_line_string.rotate_around_centroid(-45.),
            expected_around_centroid,
            epsilon = 1e-12
        );

        let expected_around_center = MultiLineString::new(vec![
            line_string![
                (x: -1.213203435596426, y: 17.07106781186548, z: 0.0),
                (x: 0.2010101267766693, y: 17.07106781186548, z: 1.0),
                (x: 2.322330470336312, y: 14.94974746830583, z: 2.0),
            ],
            line_string![
                (x: 12.92893218813452, y: 17.07106781186548, z: 3.0),
                (x: 27.07106781186548, y: 17.07106781186548, z: 4.0),
                (x: 41.21320343559643, y: 2.928932188134528, z: 5.0),
            ],
        ]);
        assert_relative_eq!(
            multi_line_string.rotate_around_center(-45.),
            expected_around_center,
            epsilon = 1e-12
        );
    }

    #[test]
    fn test_rotate_line_around_point() {
        let line0 = Line::new(coord!(0., 0., 0.), coord!(0., 2., 0.));
        let line1 = Line::new(coord!(0., 0., 0.), coord!(-2., 0., 0.));
        assert_relative_eq!(
            line0.rotate_around_point(90., Point::new(0., 0., 0.)),
            line1
        );
    }

    #[test]
    fn test_rotate_multipolygon_around_centroid() {
        let multipolygon: MultiPolygon = vec![
            polygon![
                (x: 0.,  y: 0.,  z: 1.),
                (x: 10., y: 0.,  z: 2.),
                (x: 10., y: 10., z: 3.),
                (x: 0.,  y: 10., z: 4.),
                (x: 0.,  y: 0.,  z: 1.)
            ],
            polygon![
                (x: 0.,  y: 0.,  z: 5.),
                (x: -10., y: 0.,  z: 6.),
                (x: -10., y: -10.,z: 7.),
                (x: 0.,  y: -10.,z: 8.),
                (x: 0.,  y: 0.,  z: 5.)
            ],
        ]
        .into();

        let expected_centroid: MultiPolygon = vec![
            polygon![
                (x: 0., y: 0.,  z: 1.),
                (x: 7.0710678118654755, y: 7.071067811865475, z: 2.),
                (x: 0., y: 14.142135623730951, z: 3.),
                (x: -7.071067811865475, y: 7.0710678118654755, z: 4.),
                (x: 0., y: 0.,  z: 1.)
            ],
            polygon![
                (x: 0., y: 0.,  z: 5.),
                (x: -7.0710678118654755, y: -7.071067811865475, z: 6.),
                (x: 0., y: -14.142135623730951, z: 7.),
                (x: 7.071067811865475, y: -7.0710678118654755, z: 8.),
                (x: 0., y: 0.,  z: 5.)
            ],
        ]
        .into();

        assert_relative_eq!(
            multipolygon.rotate_around_centroid(45.),
            expected_centroid,
            epsilon = 1e-12
        );
    }

    #[test]
    fn test_rotate_multipolygons() {
        let multipolygon: MultiPolygon = vec![
            polygon![
                (x: 1., y: 1.,  z: 1.),
                (x: 2., y: 1.,  z: 2.),
                (x: 2., y: 10., z: 3.),
                (x: 1., y: 10., z: 4.),
                (x: 1., y: 1.,  z: 1.),
            ],
            polygon![
                (x: 10., y: 1.,  z: 5.),
                (x: 12., y: 1.,  z: 6.),
                (x: 12., y: 12., z: 7.),
                (x: 10., y: 12., z: 8.),
                (x: 10., y: 1.,  z: 5.),
            ],
        ]
        .into();

        let expected_center: MultiPolygon = vec![
            polygon![
                (x: -0.2360967926537398, y: 2.610912703473988, z: 1.),
                (x: 0.7298290336353284,  y: 2.352093658371467, z: 2.),
                (x: 3.059200439558015,   y: 11.04542609497308,  z: 3.),
                (x: 2.093274613268947,   y: 11.3042451400756,   z: 4.),
                (x: -0.2360967926537398, y: 2.610912703473988, z: 1.),
            ],
            polygon![
                (x: 8.457235643947875, y: 0.2815412975513012, z: 5.),
                (x: 10.38908729652601, y: -0.2360967926537403, z: 6.),
                (x: 13.23609679265374, y: 10.38908729652601,  z: 7.),
                (x: 11.3042451400756,  y: 10.90672538673105,   z: 8.),
                (x: 8.457235643947875, y: 0.2815412975513012, z: 5.),
            ],
        ]
        .into();

        let expected_centroid: MultiPolygon = vec![
            polygon![
                (x: -0.1016007672888048, y: 3.05186627999456, z: 1.),
                (x: 0.8643250590002634,  y: 2.793047234892039, z: 2.),
                (x: 3.19369646492295,    y: 11.48637967149365,  z: 3.),
                (x: 2.227770638633882,   y: 11.74519871659617,  z: 4.),
                (x: -0.1016007672888048, y: 3.05186627999456, z: 1.),
            ],
            polygon![
                (x: 8.59173166931281, y: 0.7224948740718733, z: 5.),
                (x: 10.52358332189095, y: 0.2048567838668318, z: 6.),
                (x: 13.37059281801868, y: 10.83004087304658,  z: 7.),
                (x: 11.43874116544054, y: 11.34767896325162,  z: 8.),
                (x: 8.59173166931281, y: 0.7224948740718733, z: 5.),
            ],
        ]
        .into();

        assert_relative_eq!(
            multipolygon.rotate_around_center(-15.),
            expected_center,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            multipolygon.rotate_around_centroid(-15.),
            expected_centroid,
            epsilon = 1e-12
        );
    }

    #[test]
    fn test_rotate_empty_geometries_error_gracefully() {
        // line string
        let empty_linestring: LineString = line_string![];
        let rotated_empty_linestring = empty_linestring.rotate_around_centroid(90.);
        assert_eq!(empty_linestring, rotated_empty_linestring);

        // multi line string
        let empty_multilinestring: MultiLineString = MultiLineString::new(vec![]);
        let rotated_empty_multilinestring = empty_multilinestring.rotate_around_centroid(90.);
        assert_eq!(empty_multilinestring, rotated_empty_multilinestring);

        // polygon
        let empty_polygon: Polygon<f64> = polygon![];
        let rotated_empty_polygon = empty_polygon.rotate_around_centroid(90.);
        assert_eq!(empty_polygon, rotated_empty_polygon);

        // multi polygon
        let empty_multipolygon: MultiPolygon = Vec::<Polygon<f64>>::new().into();
        let rotated_empty_multipolygon = empty_multipolygon.rotate_around_centroid(90.);
        assert_eq!(empty_multipolygon, rotated_empty_multipolygon);
    }
}
