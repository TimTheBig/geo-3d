use super::Distance;
use crate::{CoordNum, Line, LineString, MultiLineString, Point};

/// Calculate the length of a `Line`, `LineString`, or `MultiLineString` in a given [metric space](crate::algorithm::line_measures::metric_spaces).
///
/// # Examples
/// ```
/// use geo::algorithm::line_measures::{Length, Euclidean};
///
/// let line_string = geo::wkt!(LINESTRING(
///     0.0 0.0 0.0,
///     3.0 4.0 5.0,
///     3.0 5.0 7.0
/// ));
/// assert_eq!(line_string.length::<Euclidean>(), 6.);
/// ```
pub trait Length<F: CoordNum> {
    fn length<Euclidean: Distance<F, Point<F>, Point<F>>>(&self) -> F;
}

impl<F: CoordNum> Length<F> for Line<F> {
    fn length<Euclidean: Distance<F, Point<F>, Point<F>>>(&self) -> F {
        Euclidean::distance(self.start_point(), self.end_point())
    }
}

impl<F: CoordNum> Length<F> for LineString<F> {
    fn length<Euclidean: Distance<F, Point<F>, Point<F>>>(&self) -> F {
        let mut length = F::zero();
        for line in self.lines() {
            length = length + line.length::<Euclidean>();
        }
        length
    }
}

impl<F: CoordNum> Length<F> for MultiLineString<F> {
    fn length<Euclidean: Distance<F, Point<F>, Point<F>>>(&self) -> F {
        let mut length = F::zero();
        for line in self {
            length = length + line.length::<Euclidean>();
        }
        length
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{coord, Euclidean};

    #[test]
    fn lines() {
        // london to paris
        let line = Line::new(
            coord!(x: -0.1278f64, y: 51.5074),
            coord!(x: 2.3522, y: 48.8566),
        );

        // computing Euclidean length of an unprojected (lng/lat) line gives a nonsense answer
        assert_eq!(
            4., // nonsense!
            line.length::<Euclidean>().round()
        );
        // london to paris in EPSG:3035
        let projected_line = Line::new(
            coord!(x: 3620451.74f64, y: 3203901.44),
            coord!(x: 3760771.86, y: 2889484.80),
        );
        assert_eq!(344_307., projected_line.length::<Euclidean>().round());
    }

    #[test]
    fn line_strings() {
        let line_string = LineString::new(vec![
            coord!(x: -58.3816f64, y: -34.6037, z: 0.0), // Buenos Aires, Argentina
            coord!(x: -77.0428, y: -12.0464, z: 0.0),    // Lima, Peru
            coord!(x: -47.9292, y: -15.7801, z: 0.0),    // Brasília, Brazil
        ]);

        // computing Euclidean length of an unprojected (lng/lat) gives a nonsense answer
        assert_eq!(
            59., // nonsense!
            line_string.length::<Euclidean>().round()
        );
        // EPSG:102033
        let projected_line_string = LineString::from(vec![
            coord!(x: 143042.46f64, y: -1932485.45, z: 0.0), // Buenos Aires, Argentina
            coord!(x: -1797084.08, y: 583528.84, z: 0.0),    // Lima, Peru
            coord!(x: 1240052.27, y: 207169.12, z: 0.0),     // Brasília, Brazil
        ]);
        assert_eq!(
            6_237_538.,
            projected_line_string.length::<Euclidean>().round()
        );
    }

    #[test]
    fn line_strings_3d() {
        let line_string = LineString::new(vec![
            coord!(x: 1.0, y: 1.0, z: 1.0),
            coord!(x: 2.0, y: 2.0, z: 2.0),
            coord!(x: 3.0, y: 3.0, z: 3.0),
        ]);

        assert_eq!(
            2.0,
            line_string.length::<Euclidean>()
        );

        let line_string = LineString::new(vec![
            coord!(x: 1.0, y: -1.0, z: 1.0),
            coord!(x: 2.0, y: -2.0, z: 2.0),
            coord!(x: 3.0, y: -3.0, z: 3.0),
        ]);

        assert_eq!(
            2.0,
            line_string.length::<Euclidean>()
        );
    }
}
