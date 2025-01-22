use geo_types::{CoordNum, Point};

/// Calculate the bearing between two points.
pub trait Bearing<F: CoordNum> {
    /// Calculate the bearing from `origin` to `destination` in degrees.
    ///
    /// See [specific implementations](#implementors) for details.
    ///
    /// # Units
    /// - `origin`, `destination`: Point where the units of x/y depend on the [trait implementation](#implementors).
    /// - returns: degrees, where: North: 0°, East: 90°, South: 180°, West: 270°
    fn bearing(origin: Point<F>, destination: Point<F>) -> F;
}
