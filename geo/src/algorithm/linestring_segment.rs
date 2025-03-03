use crate::algorithm::{Densify, Length, LineInterpolatePoint, LinesIter};
use crate::geometry::{Coord, LineString, MultiLineString};

/// Segments a LineString into `segment_count` equal length LineStrings as a MultiLineString
/// using Euclidean distance calculations.
///
/// `None` will be returned when `segment_count` is equal to 0 or when a point
/// cannot be interpolated on a `Line` segment.
///
/// # Examples
/// ```
/// use geo_3d::{LineString, MultiLineString, LineStringSegmentize};
/// // Create a simple line string
/// let lns: LineString<f64> = vec![[0.0, 0.0, 0.0], [1.0, 2.0, 3.0], [3.0, 6.0, 3.0]].into();
/// // Segment it into 6 LineStrings inside of a MultiLineString
/// let segmentized = lns.line_segmentize(6).unwrap();
/// // Compare the number of elements
/// assert_eq!(6, segmentized.0.len());
///```
pub trait LineStringSegmentize {
    fn line_segmentize(&self, segment_count: usize) -> Option<MultiLineString>;
}

macro_rules! implement_segmentize {
    ($trait_name:ident, $method_name:ident, $metric_space:ty) => {
        impl $trait_name for LineString {
            fn $method_name(&self, n: usize) -> Option<MultiLineString> {
                if (n == usize::MIN) || (n == usize::MAX) {
                    return None;
                } else if n == 1 {
                    let mlns = MultiLineString::from(self.clone());
                    return Some(mlns);
                }

                let mut res_coords: Vec<Vec<Coord>> = Vec::with_capacity(n);
                let total_length = self.length();
                let mut cum_length = 0_f64;
                let segment_prop = (1_f64) / (n as f64);
                let segment_length = total_length * segment_prop;
                let densified = self.densify(segment_length - f64::EPSILON);

                if densified.lines().count() == n {
                    let linestrings = densified
                        .lines()
                        .map(LineString::from)
                        .collect::<Vec<LineString>>();
                    return Some(MultiLineString::new(linestrings));
                }

                let n_lines = densified.lines().count();
                let lns = densified.lines_iter();
                let mut ln_vec: Vec<Coord> = Vec::new();

                for (i, segment) in lns.enumerate() {
                    if i == 0 {
                        ln_vec.push(segment.start)
                    }

                    let length = segment.length();
                    cum_length += length;

                    if (cum_length >= segment_length) && (i != (n_lines - 1)) {
                        let remainder = cum_length - segment_length;
                        let endpoint =
                            segment.line_interpolate_point((length - remainder) / length)?;

                        ln_vec.push(endpoint.into());
                        let to_push = ln_vec.drain(..);
                        res_coords.push(to_push.collect::<Vec<Coord>>());

                        if i != n_lines {
                            ln_vec.push(endpoint.into());
                        }
                        cum_length = remainder;
                    }
                    ln_vec.push(segment.end);
                }

                res_coords.push(ln_vec);
                let res_lines = res_coords
                    .into_iter()
                    .map(LineString::new)
                    .collect::<Vec<LineString>>();
                Some(MultiLineString::new(res_lines))
            }
        }
    };
}

implement_segmentize!(LineStringSegmentize, line_segmentize, Euclidean);

#[cfg(test)]
mod test {
    use approx::RelativeEq;
    use super::*;
    use crate::LineString;

    #[test]
    fn n_elems_bug() {
        // Test for an edge case that seems to fail:
        // https://github.com/georust/geo/issues/1075
        // https://github.com/JosiahParry/rsgeo/issues/28

        let linestring: LineString = vec![
            [324957.69921197, 673670.123131518, 324957.69921197],
            [324957.873557727, 673680.139281405, 324957.873557727],
            [324959.863123514, 673686.784106964, 324959.863123514],
            [324961.852683597, 673693.428933452, 324961.852683597],
            [324963.822867622, 673698.960855279, 324963.822867622],
            [324969.636546456, 673709.992098018, 324969.636546456],
            [324976.718443977, 673722.114520549, 324976.718443977],
            [324996.443964294, 673742.922904206, 324996.443964294],
        ]
        .into();
        let segments = linestring.line_segmentize(2).unwrap();
        assert_eq!(segments.0.len(), 2);
        let segments = linestring.line_segmentize(3).unwrap();
        assert_eq!(segments.0.len(), 3);
        let segments = linestring.line_segmentize(4).unwrap();
        assert_eq!(segments.0.len(), 4);

        assert_eq!(
            segments.length(),
            linestring.length()
        );
    }

    #[test]
    fn long_end_segment() {
        let linestring: LineString = vec![
            [325581.792390628, 674398.495901267, 325581.792390628],
            [325585.576868499, 674400.657039341, 325585.576868499],
            [325589.966469742, 674401.694493658, 325589.966469742],
            [325593.750940609, 674403.855638851, 325593.750940609],
            [325599.389217394, 674404.871546368, 325599.389217394],
            [325604.422360924, 674407.011146146, 325604.422360924],
            [325665.309662534, 674424.885671739, 325665.309662534],
        ]
        .into();

        let segments = linestring.line_segmentize(5).unwrap();
        assert_eq!(segments.0.len(), 5);
        assert_relative_eq!(
            linestring.length(),
            segments.length(),
            epsilon = f64::EPSILON
        );
    }

    #[test]
    fn two_coords() {
        let linestring: LineString = vec![[0.0, 0.0, 0.0], [0.0, 1.0, 2.0]].into();

        let segments = linestring.line_segmentize(5).unwrap();
        assert_eq!(segments.0.len(), 5);
        assert_relative_eq!(
            linestring.length(),
            segments.length(),
            epsilon = f64::EPSILON
        );
    }

    #[test]
    fn long_middle_segments() {
        let linestring: LineString = vec![
            [325403.816883668, 673966.295402012, 325403.816883668],
            [325410.280933752, 673942.805501254, 325410.280933752],
            [325410.280933752, 673942.805501254, 325410.280933752],
            [325439.782082601, 673951.201057316, 325439.782082601],
            [325439.782082601, 673951.201057316, 325439.782082601],
            [325446.064640793, 673953.318876004, 325446.064640793],
            [325446.064640793, 673953.318876004, 325446.064640793],
            [325466.14184472, 673958.537886844, 325466.14184472],
            [325466.14184472, 673958.537886844, 325466.14184472],
            [325471.799973648, 673960.666539074, 325471.799973648],
            [325471.799973648, 673960.666539074, 325471.799973648],
            [325518.255916084, 673974.335722824, 325518.255916084],
            [325518.255916084, 673974.335722824, 325518.255916084],
            [325517.669972133, 673976.572326305, 325517.669972133],
            [325517.669972133, 673976.572326305, 325517.669972133],
            [325517.084028835, 673978.808929878, 325517.084028835],
            [325517.084028835, 673978.808929878, 325517.084028835],
            [325515.306972763, 673984.405833764, 325515.306972763],
            [325515.306972763, 673984.405833764, 325515.306972763],
            [325513.549152184, 673991.115645844, 325513.549152184],
            [325513.549152184, 673991.115645844, 325513.549152184],
            [325511.772106396, 673996.712551354, 325511.772106396],
        ]
        .into();

        let segments = linestring.line_segmentize(5).unwrap();
        assert_eq!(segments.0.len(), 5);

        assert_relative_eq!(
            linestring.length(),
            segments.length(),
            epsilon = f64::EPSILON
        );
    }

    #[test]
    // that 0 returns None and that usize::MAX returns None
    fn n_is_zero() {
        let linestring: LineString = vec![[-1.0, 0.0, -1.0], [0.5, 1.0, 0.5], [1.0, 2.0, 3.0]].into();
        let segments = linestring.line_segmentize(0);
        assert!(segments.is_none())
    }

    #[test]
    fn n_is_max() {
        let linestring: LineString = vec![[-1.0, 0.0, -2.0], [0.5, 1.0, 0.5], [1.0, 2.0, 3.0]].into();
        let segments = linestring.line_segmentize(usize::MAX);
        assert!(segments.is_none())
    }

    #[test]
    fn n_greater_than_lines() {
        let linestring: LineString = vec![[-1.0, 0.0, -1.0], [0.5, 1.0, 0.5], [1.0, 2.0, 1.0]].into();
        let segments = linestring.line_segmentize(5).unwrap();

        // assert that there are n linestring segments
        assert_eq!(segments.0.len(), 5);

        // assert that the lines are equal length
        let lens = segments
            .into_iter()
            .map(|x| x.length())
            .collect::<Vec<f64>>();

        let first = lens[0];

        assert!(lens
            .iter()
            .all(|x| first.relative_eq(x, f64::EPSILON, 1e-10)))
    }

    #[test]
    // test the cumulative length is the same
    fn cumul_length() {
        let linestring: LineString = vec![[0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [1.0, 2.0, 3.0], [3.0, 3.0, 3.0]].into();
        let segments = linestring.line_segmentize(2).unwrap();

        assert_relative_eq!(
            linestring.length(),
            segments.length(),
            epsilon = f64::EPSILON
        )
    }

    #[test]
    fn n_elems() {
        let linestring: LineString = vec![[0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [1.0, 2.0, 3.0], [3.0, 3.0, 3.0]].into();
        let segments = linestring.line_segmentize(2).unwrap();
        assert_eq!(segments.0.len(), 2)
    }

    #[test]
    fn tiny_distances() {
        // this test is to ensure that at super small distances
        // the number of units is still the specified one.
        let linestring: LineString = vec![
            [-3.19416, 55.95524, -3.19416],
            [-3.19352, 55.95535, -3.19352],
            [-3.19288, 55.95546, -3.19288],
        ]
        .into();

        let n = 8;
        let segments = linestring.line_segmentize(n).unwrap();
        assert_eq!(segments.0.len(), n)
    }
}
