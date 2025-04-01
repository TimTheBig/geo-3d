use crate::{Coord, CoordNum, MapCoords};

/// Convert (infalliby) the type of a geometry’s coordinate value.
///
/// # Examples
///
/// ```rust
/// use geo_3d::{Convert, LineString, line_string};
///
/// let line_string_32: LineString<f32> = line_string![
///     (x: 5., y: 10., z: 15.),
///     (x: 3., y: 1., z: -3),
///     (x: 8., y: 9., z: 10.),
/// ];
///
/// let line_string_64: LineString<f64> = line_string_32.convert();
/// ```
///
pub trait Convert<T, U> {
    type Output;

    fn convert(&self) -> Self::Output;
}
impl<G, T: CoordNum, U: CoordNum> Convert<T, U> for G
where
    G: MapCoords<T, U>,
    U: From<T>,
{
    type Output = <Self as MapCoords<T, U>>::Output;

    fn convert(&self) -> Self::Output {
        self.map_coords(|Coord { x, y, z }| Coord {
            x: x.into(),
            y: y.into(),
            z: z.into(),
        })
    }
}

/// Convert (fallibly) the type of a geometry’s coordinate value.
///
/// # Examples
///
/// ```rust
/// use geo_3d::{TryConvert, LineString, line_string};
///
/// let line_string_64: LineString<f64> = line_string![
///     (x: 5., y: 10., z: 15.),
///     (x: 3., y: 1., z: -3),
///     (x: 8., y: 9., z: 10.),
/// ];
///
/// let line_string_32: Result<LineString<f32>, _> = line_string_64.try_convert();
/// ```
///
pub trait TryConvert<T, U> {
    type Output;

    fn try_convert(&self) -> Self::Output;
}
impl<G, T: CoordNum, U: CoordNum> TryConvert<T, U> for G
where
    G: MapCoords<T, U>,
    U: TryFrom<T>,
{
    type Output = Result<<Self as MapCoords<T, U>>::Output, <U as TryFrom<T>>::Error>;

    fn try_convert(&self) -> Self::Output {
        self.try_map_coords(|Coord { x, y, z }| {
            Ok(Coord {
                x: x.try_into()?,
                y: y.try_into()?,
                z: z.try_into()?,
            })
        })
    }
}
