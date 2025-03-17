use geo_types::Coord;
use geo_types::CoordNum;

use crate::{MapCoords, MapCoordsInPlace};

pub trait ToRadians<T: CoordNum>:
    Sized + MapCoords<T, T, Output = Self> + MapCoordsInPlace<T>
{
    fn to_radians(&self) -> Self {
        self.map_coords(|Coord { x, y, z }| Coord {
            x: x.to_radians(),
            y: y.to_radians(),
            z: z.to_radians(),
        })
    }

    fn to_radians_in_place(&mut self) {
        self.map_coords_in_place(|Coord { x, y, z }| Coord {
            x: x.to_radians(),
            y: y.to_radians(),
            z: z.to_radians(),
        })
    }
}
impl<T: CoordNum, G: MapCoords<T, T, Output = Self> + MapCoordsInPlace<T>> ToRadians<T> for G {}

pub trait ToDegrees<T: CoordNum>:
    Sized + MapCoords<T, T, Output = Self> + MapCoordsInPlace<T>
{
    fn to_degrees(&self) -> Self {
        self.map_coords(|Coord { x, y, z }| Coord {
            x: x.to_degrees(),
            y: y.to_degrees(),
            z: z.to_degrees(),
        })
    }

    fn to_degrees_in_place(&mut self) {
        self.map_coords_in_place(|Coord { x, y, z }| Coord {
            x: x.to_degrees(),
            y: y.to_degrees(),
            z: z.to_degrees(),
        })
    }
}
impl<T: CoordNum, G: MapCoords<T, T, Output = Self> + MapCoordsInPlace<T>> ToDegrees<T> for G {}

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;
    use approx::assert_relative_eq;
    use geo_types::{coord, Line};
    use super::*;

    const fn line_degrees_mock() -> Line {
        Line::new(coord!(90.0, 180., 90.0), coord!(0., -90., 0.))
    }

    fn line_radians_mock() -> Line {
        Line::new(coord!(PI / 2., PI, PI / 2.), coord!(0., -PI / 2., 0.))
    }

    #[test]
    fn converts_to_radians() {
        assert_relative_eq!(line_radians_mock(), line_degrees_mock().to_radians())
    }

    #[test]
    fn converts_to_radians_in_place() {
        let mut line = line_degrees_mock();
        line.to_radians_in_place();
        assert_relative_eq!(line_radians_mock(), line)
    }

    #[test]
    fn converts_to_degrees() {
        assert_relative_eq!(line_degrees_mock(), line_radians_mock().to_degrees())
    }

    #[test]
    fn converts_to_degrees_in_place() {
        let mut line = line_radians_mock();
        line.to_degrees_in_place();
        assert_relative_eq!(line_degrees_mock(), line)
    }
}
