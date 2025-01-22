/// Creates a [`Point`] from the given coordinates.
///
/// ```no_run
/// point! { x: <number>, y: <number>, z: <number> }
/// point!(<coordinate>)
/// ```
///
/// # Examples
///
/// Creating a [`Point`], supplying x/y/z values:
///
/// ```
/// use geo_types::{point, coord};
///
/// let p = point! { x: 181.2, y: 51.79, z: 72.52 };
///
/// assert_eq!(p.x(), 181.2);
/// assert_eq!(p.y(), 51.79);
/// assert_eq!(p.z(), 72.52);
///
/// let p = point!(coord! { x: 181.2, y: 51.79, z: 91.04 });
///
/// assert_eq!(p.x(), 181.2);
/// assert_eq!(p.y(), 51.79);
/// assert_eq!(p.z(), 91.04);
/// ```
///
/// [`Point`]: ./struct.Point.html
#[macro_export]
macro_rules! point {
    ( $($tag:tt : $val:expr),* $(,)? ) => {
        $crate::point! ( $crate::coord!($( $tag: $val , )*) )
    };
    ($x:expr, $y:expr, $z:expr) => {
        $crate::coord!($x, $y, $z)
    };
    ( $coord:expr $(,)? ) => {
        $crate::Point::from($coord)
    };
}

/// Creates a [`Coord`] from the given scalars.
///
/// ```no_run
/// coord! { x: <number>, y: <number>, z: <number> }
/// ```
///
/// # Examples
///
/// Creating a [`Coord`], supplying x/y/z values:
///
/// ```
/// use geo_types::coord;
///
/// let c = coord! { x: 181.2, y: 51.79, z: 82.916 };
/// let c1 = coord!(181.2, 51.79, 82.916);
///
/// assert_eq!(c, geo_types::coord! { x: 181.2, y: 51.79, z: 82.916 });
/// assert_eq!(c1, geo_types::coord! { x: 181.2, y: 51.79, z: 82.916 });
/// ```
///
/// [`Coord`]: ./struct.Coord.html
#[macro_export]
macro_rules! coord {
    (x: $x:expr, y: $y:expr $(,)* ) => {
        compile_error!("Coord is 3d, it has x, y, and z fields")
    };
    (x: $x:expr, y: $y:expr, z: $z:expr $(,)* ) => {
        $crate::Coord { x: $x, y: $y, z: $z }
    };
    ($x:expr, $y:expr) => {
        compile_error!("Coord is 3d, it has x, y, and z fields")
    };
    ($x:expr, $y:expr, $z:expr) => {
        $crate::Coord { x: $x, y: $y, z: $z }
    };
}

/// Creates a [`LineString`] containing the given coordinates.
///
/// ```no_run
/// line_string![Coord OR (x: <number>, y: <number>, z: <number>), …]
/// ```
///
/// # Examples
///
/// Creating a [`LineString`], supplying x/y/z values:
///
/// ```
/// use geo_types::line_string;
///
/// let ls = line_string![
///     (x: -21.95156, y: 64.1446, z: 82.074),
///     (x: -21.951, y: 64.14479, z: 82.07814),
///     (x: -21.95044, y: 64.14527, z: 82.941),
///     (x: -21.951445, y: 64.145508, z: 82.07424),
/// ];
///
/// assert_eq!(ls[1], geo_types::coord! {
///     x: -21.951,
///     y: 64.14479,
///     z: 82.07814,
/// });
/// ```
///
/// Creating a [`LineString`], supplying [`Coord`]s:
///
/// ```
/// use geo_types::line_string;
///
/// let coord1 = geo_types::coord! {
///     x: -21.95156,
///     y: 64.1446,
///     z: 82.074,
/// };
/// let coord2 = geo_types::coord! {
///     x: -21.951,
///     y: 64.14479,
///     z: 82.07814,
/// };
/// let coord3 = geo_types::coord! {
///     x: -21.95044,
///     y: 64.14527,
///     z: 82.941,
/// };
/// let coord4 = geo_types::coord! {
///     x: -21.951445,
///     y: 64.145508,
///     z: 82.07424,
/// };
///
/// let ls = line_string![coord1, coord2, coord3, coord4];
///
/// assert_eq!(
///     ls[1],
///     geo_types::coord! {
///         x: -21.951,
///         y: 64.14479,
///         z: 82.07814,
///     }
/// );
/// ```
///
/// [`Coord`]: ./struct.Coord.html
/// [`LineString`]: ./line_string/struct.LineString.html
#[macro_export]
macro_rules! line_string {
    () => { $crate::LineString::new($crate::_alloc::vec![]) };
    (
        $(( $($tag:tt : $val:expr),* $(,)? )),*
        $(,)?
    ) => {
        line_string![
            $(
                $crate::coord! { $( $tag: $val , )* },
            )*
        ]
    };
    (
        $($coord:expr),*
        $(,)?
    ) => {
        $crate::LineString::new(
            $crate::_alloc::vec![
                $($coord),*
            ]
        )
    };
}

/// Creates a [`Polygon`] containing the given coordinates.
///
/// ```no_run
/// polygon![Coord OR (x: <number>, y: <number>, z: <number>), …]
///
/// // or
///
/// polygon!(
///     exterior: [Coord OR (x: <number>, y: <number>, z: <number>), …],
///     interiors: [
///         [Coord OR (x: <number>, y: <number>, z: <number>), …],
///         …
///     ],
/// )
/// ```
///
/// # Examples
///
/// Creating a [`Polygon`] without interior rings, supplying x/y/z values:
///
/// ```
/// use geo_types::polygon;
///
/// let poly = polygon![
///     (x: -111., y: 45., z: 8.4),
///     (x: -111., y: 41., z: 8.8),
///     (x: -104., y: 41., z: 8.3),
///     (x: -104., y: 45., z: 8.9),
/// ];
///
/// assert_eq!(
///     poly.exterior()[1],
///     geo_types::coord! { x: -111., y: 41., z: 8.8 },
/// );
/// ```
///
/// Creating a [`Polygon`], supplying x/y/z values:
///
/// ```
/// use geo_types::polygon;
///
/// let poly = polygon!(
///     exterior: [
///         (x: -111., y: 45., z: 8.6),
///         (x: -111., y: 41., z: 8.3),
///         (x: -104., y: 41., z: 8.1),
///         (x: -104., y: 45., z: 8.0),
///     ],
///     interiors: [
///         [
///             (x: -110., y: 44., z: 8.),
///             (x: -110., y: 42., z: 8.),
///             (x: -105., y: 42., z: 8.),
///             (x: -105., y: 44., z: 8.),
///         ],
///     ],
/// );
///
/// assert_eq!(
///     poly.exterior()[1],
///     geo_types::coord! { x: -111., y: 41., z: 8.3 },
/// );
/// ```
///
/// [`Coord`]: ./struct.Coord.html
/// [`Polygon`]: ./struct.Polygon.html
#[macro_export]
macro_rules! polygon {
    () => { $crate::Polygon::new($crate::line_string![], $crate::_alloc::vec![]) };
    (
        exterior: [
            $(( $($exterior_tag:tt : $exterior_val:expr),* $(,)? )),*
            $(,)?
        ],
        interiors: [
            $([
                $(( $($interior_tag:tt : $interior_val:expr),* $(,)? )),*
                $(,)?
            ]),*
            $(,)?
        ]
        $(,)?
    ) => {
        polygon!(
            exterior: [
                $(
                    $crate::coord! { $( $exterior_tag: $exterior_val , )* },
                )*
            ],
            interiors: [
                $([
                    $($crate::coord! { $( $interior_tag: $interior_val , )* }),*
                ]),*
            ],
        )
    };
    (
        exterior: [
            $($exterior_coord:expr),*
            $(,)?
        ],
        interiors: [
            $([
                $($interior_coord:expr),*
                $(,)?
            ]),*
            $(,)?
        ]
        $(,)?
    ) => {
        $crate::Polygon::new(
            $crate::line_string![
                $($exterior_coord), *
            ],
            $crate::_alloc::vec![
                $(
                    $crate::line_string![$($interior_coord),*]
                ), *
            ]
        )
    };
    (
        $(( $($tag:tt : $val:expr),* $(,)? )),*
        $(,)?
    ) => {
        polygon![
            $($crate::coord! { $( $tag: $val , )* }),*
        ]
    };
    (
        $($coord:expr),*
        $(,)?
    ) => {
        $crate::Polygon::new(
            $crate::line_string![$($coord,)*],
            $crate::_alloc::vec![],
        )
    };
}

#[cfg(test)]
mod test {
    #[test]
    fn test_point() {
        let p = point! { x: 1.2, y: 3.4, z: 7.4 };
        assert_eq!(p.x(), 1.2);
        assert_eq!(p.y(), 3.4);
        assert_eq!(p.z(), 7.4);

        let p = point! {
            x: 1.2,
            y: 3.4,
            z: 7.4,
        };
        assert_eq!(p.x(), 1.2);
        assert_eq!(p.y(), 3.4);
        assert_eq!(p.z(), 7.4);

        let p = point!(coord! { x: 1.2, y: 3.2, z: 7.42 });
        assert_eq!(p.x(), 1.2);
        assert_eq!(p.y(), 3.2);
        assert_eq!(p.z(), 7.42);

        let p = point!(coord! { x: 1.2, y: 3.4, z: 7.4 },);
        assert_eq!(p.x(), 1.2);
        assert_eq!(p.y(), 3.4);
        assert_eq!(p.z(), 7.4);
    }

    #[test]
    fn test_line() {
        let ls = line_string![(x: -1.2f32, y: 3.4f32, z: -7.5f32)];
        assert_eq!(ls[0], coord! { x: -1.2, y: 3.4, z: -7.5 });

        let ls = line_string![
            (x: -1.2f32, y: 3.4f32, z: 7.4f32),
        ];
        assert_eq!(ls[0], coord! { x: -1.2, y: 3.4, z: 7.4 });

        let ls = line_string![(
            x: -1.2f32,
            y: 3.4f32,
            z: 7.1f32,
        )];
        assert_eq!(ls[0], coord! { x: -1.2, y: 3.4, z: 7.1 });

        let ls = line_string![
            (x: -1.2f32, y: 3.4f32, z: 4.3f32),
            (x: -5.6, y: 7.8, z: -9.6),
        ];
        assert_eq!(ls[0], coord! { x: -1.2, y: 3.4, z: 4.3 });
        assert_eq!(ls[1], coord! { x: -5.6, y: 7.8, z: -9.6 });
    }

    #[test]
    fn test_polygon() {
        let p = polygon!(
            exterior: [(x: 1, y: 2, z: 3)],
            interiors: [[(x: 3, y: 4, z: 5)]]
        );
        assert_eq!(p.exterior()[0], coord! { x: 1, y: 2, z: 3 });
        assert_eq!(p.interiors()[0][0], coord! { x: 3, y: 4, z: 5 });

        let p = polygon!(
            exterior: [(x: 1, y: 2, z: 3)],
            interiors: [[(x: 3, y: 4, z: 5)]],
        );
        assert_eq!(p.exterior()[0], coord! { x: 1, y: 2, z: 3 });
        assert_eq!(p.interiors()[0][0], coord! { x: 3, y: 4, z: 5 });

        let p = polygon!(
            exterior: [(x: 1, y: 2, z: 3, )],
            interiors: [[(x: 3, y: 4, z: 5, )]],
        );
        assert_eq!(p.exterior()[0], coord! { x: 1, y: 2, z: 3 });
        assert_eq!(p.interiors()[0][0], coord! { x: 3, y: 4, z: 5 });

        let p = polygon!(
            exterior: [(x: 1, y: 2, z: 3, ), ],
            interiors: [[(x: 3, y: 4, z: 5, ), ]],
        );
        assert_eq!(p.exterior()[0], coord! { x: 1, y: 2, z: 3 });
        assert_eq!(p.interiors()[0][0], coord! { x: 3, y: 4, z: 5 });

        let p = polygon!(
            exterior: [(x: 1, y: 2, z: 3, ), ],
            interiors: [[(x: 3, y: 4, z: 5, ), ], ],
        );
        assert_eq!(p.exterior()[0], coord! { x: 1, y: 2, z: 3 });
        assert_eq!(p.interiors()[0][0], coord! { x: 3, y: 4, z: 5 });
    }
}
