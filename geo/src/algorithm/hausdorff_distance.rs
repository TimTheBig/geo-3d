use crate::algorithm::Distance;
use crate::CoordsIter;
use crate::GeoNum;
use geo_types::{Coord, Point};
use num_traits::Bounded;

/// Determine the distance between two geometries using the [Hausdorff distance formula].
///
/// Hausdorff distance is used to compare two point sets. It measures the maximum euclidean
/// distance of a point in one set to the nearest point in another set. Hausdorff distance
/// is often used to measure the amount of mismatch between two sets.
///
/// [Hausdorff distance formula]: https://en.wikipedia.org/wiki/Hausdorff_distance
pub trait HausdorffDistance<T: GeoNum> {
    fn hausdorff_distance<Rhs>(&self, rhs: &Rhs) -> T
    where
        Rhs: CoordsIter<Scalar = T>;
}

impl<T, G> HausdorffDistance<T> for G
where
    T: GeoNum,
    G: CoordsIter<Scalar = T>,
{
    fn hausdorff_distance<Rhs>(&self, rhs: &Rhs) -> T
    where
        Rhs: CoordsIter<Scalar = T>,
    {
        // calculate from A -> B
        let hd1 = self
            .coords_iter()
            .map(|c| {
                rhs.coords_iter()
                    .map(|c2| c.distance(c2))
                    .fold(<T as Bounded>::max_value(), |accum, val| accum.min(val))
            })
            .fold(<T as Bounded>::min_value(), |accum, val| accum.max(val));

        // Calculate from B -> A
        let hd2 = rhs
            .coords_iter()
            .map(|c| {
                self.coords_iter()
                    .map(|c2| c.distance(c2))
                    .fold(<T as Bounded>::max_value(), |accum, val| accum.min(val))
            })
            .fold(<T as Bounded>::min_value(), |accum, val| accum.max(val));

        // The max of the two
        hd1.max(hd2)
    }
}

// ┌───────────────────────────┐
// │ Implementations for Coord │
// └───────────────────────────┘

impl<T> HausdorffDistance<T> for Coord<T>
where
    T: GeoNum,
{
    fn hausdorff_distance<Rhs>(&self, rhs: &Rhs) -> T
    where
        Rhs: CoordsIter<Scalar = T>,
    {
        Point::from(*self).hausdorff_distance(rhs)
    }
}

#[cfg(test)]
mod test {
    use crate::HausdorffDistance;
    use crate::{line_string, polygon, MultiPoint, MultiPolygon};

    #[test]
    fn hd_mpnt_mpnt() {
        let p1: MultiPoint<_> = vec![(0., 0., 0.), (1., 2., 3.)].into();
        let p2: MultiPoint<_> = vec![(2., 3., 4.), (1., 2., 3.)].into();
        assert_relative_eq!(p1.hausdorff_distance(&p2), 3.741657, epsilon = 1e-6);
    }

    #[test]
    fn hd_mpnt_poly() {
        let p1: MultiPoint<_> = vec![(0., 0., 0.), (1., 2., 3.)].into();
        let poly = polygon![
            (x: 1., y: -3.1, z: 1.), (x: 3.7, y: 2.7, z: 3.7),
            (x: 0.9, y: 7.6, z: 0.9), (x: -4.8, y: 6.7, z: -4.8),
            (x: -7.5, y: 0.9, z: -7.5), (x: -4.7, y: -4., z: -4.7),
            (x: 1., y: -3.1, z: 1.)
        ];

        assert_relative_eq!(p1.hausdorff_distance(&poly), 10.644716, epsilon = 1e-6)
    }

    #[test]
    fn hd_mpnt_lns() {
        let p1: MultiPoint<_> = vec![(0., 0., 0.), (1., 2., 3.)].into();
        let lns = line_string![
            (x: 1., y: -3.1, z: 1.), (x: 3.7, y: 2.7, z: 3.7),
            (x: 0.9, y: 7.6, z: 0.9), (x: -4.8, y: 6.7, z: -4.8),
            (x: -7.5, y: 0.9, z: -7.5), (x: -4.7, y: -4., z: -4.7),
            (x: 1., y: -3.1, z: 1.)
        ];

        assert_relative_eq!(p1.hausdorff_distance(&lns), 10.64471, epsilon = 1e-5)
    }

    #[test]
    fn hd_mpnt_mply() {
        let p1: MultiPoint<_> = vec![(0., 0., 0.), (1., 2., 3.)].into();
        let multi_polygon = MultiPolygon::new(vec![
            polygon![
                (x: 0.0f32, y: 0.0, z: 0.0),
                (x: 2.0, y: 0.0, z: 2.0),
                (x: 2.0, y: 1.0, z: 2.0),
                (x: 0.0, y: 1.0, z: 0.0),
            ],
            polygon![
                (x: 1.0, y: 1.0, z: 1.0),
                (x: -2.0, y: 1.0, z: -2.0),
                (x: -2.0, y: -1.0, z: -2.0),
                (x: 1.0, y: -1.0, z: 1.0),
            ],
        ]);

        assert_eq!(p1.hausdorff_distance(&multi_polygon), 3.0)
    }
}
