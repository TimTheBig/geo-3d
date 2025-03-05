use std::iter::Sum;

use geo_types::{CoordNum, LineString, MultiPoint, MultiPolygon, Point, Polygon, Rect};

/// 3D volume of a geometry.
///
/// # Examples
///
/// ```
/// use geo::polygon;
/// use geo::Volume;
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
pub trait Volume<T>
where
    T: CoordNum,
{
    fn volume(&self) -> T;
}

impl<T: CoordNum> Volume<T> for Point<T> {
    fn volume(&self) -> T {
        T::zero()
    }
}

impl<T: CoordNum> Volume<T> for MultiPoint<T> {
    fn volume(&self) -> T {
        T::zero()
    }
}

impl<T: CoordNum> Volume<T> for Rect<T> {
    fn volume(&self) -> T {
        let diff = self.max() - self.min();

        // L * W * H
        diff.x * diff.y * diff.z
    }
}

/// 3D volume of a geometry.
///
/// # Examples
///
/// ```
/// use geo::polygon;
/// use geo::Volume;
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
pub trait TryVolume<T>
where
    T: CoordNum,
{
    /// None means the volume could not be calculated
    fn try_volume(&self) -> Option<T>;
}

impl<T: CoordNum> TryVolume<T> for Point<T> {
    fn try_volume(&self) -> Option<T> {
        None
    }
}

impl<T: CoordNum> TryVolume<T> for MultiPoint<T> {
    fn try_volume(&self) -> Option<T> {
        None
    }
}

impl<T: CoordNum> TryVolume<T> for Rect<T> {
    fn try_volume(&self) -> Option<T> {
        Some(self.volume())
    }
}

impl<T: CoordNum> TryVolume<T> for LineString<T> {
    fn try_volume(&self) -> Option<T> {
        if !self.is_closed() || !self.is_enclosed() {
            return None;
        }

        todo!("not done yet")
    }
}

impl<T: CoordNum + Sum> TryVolume<T> for Polygon<T> {
    fn try_volume(&self) -> Option<T> {
        Some(
            self.exterior().try_volume()?
            -
            self.interiors().iter()
            .map(|ls| ls.try_volume()).try_fold(T::zero(), |total, vol| -> Option<T> { Some(total + vol?) })?
        )
    }
}

impl<T: CoordNum + Sum> TryVolume<T> for MultiPolygon<T> {
    fn try_volume(&self) -> Option<T> {
        // volume of all `Polygon`
        self.iter().map(|ls| ls.try_volume()).sum()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::ConvexHull;
    use geo_types::coord;

    #[test]
    fn rect_hull_volume() {
        // inssrue that polygon volume is consistent with rect volume
        assert_eq!(
            Rect::new(coord!(-10.0, -10.0, -10.0), coord!(10.0, 10.0, 10.0)).try_volume().unwrap(),
            Rect::new(coord!(-10.0, -10.0, -10.0), coord!(10.0, 10.0, 10.0)).to_polygon().try_volume().unwrap(),
        );
        // this tests to make sure the volume of a convex_hull is smaller
        assert!(
            Rect::new(coord!(-10.0, -10.0, -10.0), coord!(10.0, 10.0, 10.0)).to_polygon().convex_hull().unwrap().try_volume().unwrap()
            <
            Rect::new(coord!(-10.0, -10.0, -10.0), coord!(10.0, 10.0, 10.0)).to_polygon().try_volume().unwrap()
        )
    }
}
