// To implement RStar’s traits in the geo-types crates, we need to access to a
// few geospatial algorithms, which are included in this hidden module. This
// hidden module is public so the geo crate can reuse these algorithms to
// prevent duplication. These functions are _not_ meant for public consumption.

use crate::{Coord, CoordNum, Line, LineString, Point, Rect};

pub fn line_string_bounding_rect<T: CoordNum>(line_string: &LineString<T>) -> Option<Rect<T>> {
    get_bounding_rect(line_string.coords().cloned())
}

pub fn line_bounding_rect<T: CoordNum>(line: Line<T>) -> Rect<T> {
    Rect::new(line.start, line.end)
}

pub fn get_bounding_rect<I, T>(collection: I) -> Option<Rect<T>>
where
    T: CoordNum,
    I: IntoIterator<Item = Coord<T>>,
{
    let mut iter = collection.into_iter();
    if let Some(pnt) = iter.next() {
        let mut xrange = (pnt.x, pnt.x);
        let mut yrange = (pnt.y, pnt.y);
        let mut zrange = (pnt.z, pnt.z);
        for pnt in iter {
            let (px, py, pz) = pnt.x_y_z();
            xrange = get_min_max(px, xrange.0, xrange.1);
            yrange = get_min_max(py, yrange.0, yrange.1);
            zrange = get_min_max(pz, zrange.0, zrange.1);
        }

        return Some(Rect::new(
            coord! {
                x: xrange.0,
                y: yrange.0,
                z: zrange.0,
            },
            coord! {
                x: xrange.1,
                y: yrange.1,
                z: zrange.1,
            },
        ));
    }
    None
}

fn get_min_max<T: PartialOrd>(p: T, min: T, max: T) -> (T, T) {
    if p > max {
        (min, p)
    } else if p < min {
        (p, max)
    } else {
        (min, max)
    }
}

pub fn line_segment_distance<T, C>(point: C, start: C, end: C) -> T
where
    T: CoordNum,
    C: Into<Coord<T>>,
{
    let point = point.into();
    let start = start.into();
    let end = end.into();

    if start == end {
        return line_euclidean_length(Line::new(point, start));
    }

    let dx = end.x - start.x;
    let dy = end.y - start.y;
    let dz = end.z - start.z;

    let d_squared = dx * dx + dy * dy + dz * dz; // Full 3D distance squared

    // Projection parameter (normalized position on the line segment)
    let r = ((point.x - start.x) * dx + (point.y - start.y) * dy + (point.z - start.z) * dz) / d_squared;

    if r <= T::zero() {
        // Closest to the start point
        return line_euclidean_length(Line::new(point, start));
    }
    if r >= T::one() {
        // Closest to the end point
        return line_euclidean_length(Line::new(point, end));
    }

    // Closest point on the segment
    let proj_x = start.x + r * dx;
    let proj_y = start.y + r * dy;
    let proj_z = start.z + r * dz;

    // Distance from the point to the projected point
    let dx_p = point.x - proj_x;
    let dy_p = point.y - proj_y;
    let dz_p = point.z - proj_z;

    (dx_p * dx_p + dy_p * dy_p + dz_p * dz_p).sqrt()
}

pub fn line_euclidean_length<T: CoordNum>(line: Line<T>) -> T {
    line.dx().hypot(line.dy()).hypot(line.dz())
}

pub fn point_line_string_euclidean_distance<T: CoordNum>(p: Point<T>, l: &LineString<T>) -> T {
    // No need to continue if the point is on the LineString, or it's empty
    if line_string_contains_point(l, p) || l.0.is_empty() {
        return T::zero();
    }
    l.lines()
        .map(|line| line_segment_distance(p.0, line.start, line.end))
        .fold(T::max_value(), |accum, val| accum.min(val))
}

pub fn point_line_euclidean_distance<C, T>(p: C, l: Line<T>) -> T
where
    T: CoordNum,
    C: Into<Coord<T>>,
{
    line_segment_distance(p.into(), l.start, l.end)
}

pub fn point_contains_point<T: CoordNum>(p1: Point<T>, p2: Point<T>) -> bool {
    let distance = line_euclidean_length(Line::new(p1.0, p2.0)).to_f64().unwrap();
    approx::relative_eq!(distance, 0.0)
}

pub fn line_string_contains_point<T: CoordNum>(
    line_string: &LineString<T>,
    point: Point<T>,
) -> bool {
    // LineString without points
    if line_string.0.is_empty() {
        return false;
    }
    // LineString with one point equal to the given point
    if line_string.0.len() == 1 {
        return point_contains_point(Point::from(line_string[0]), point);
    }
    // Check if point is a vertex of the LineString
    if line_string.0.contains(&point.0) {
        return true;
    }

    for line in line_string.lines() {
        let tx = if line.dx() == T::zero() {
            None
        } else {
            Some((point.x() - line.start.x) / line.dx())
        };
        let ty = if line.dy() == T::zero() {
            None
        } else {
            Some((point.y() - line.start.y) / line.dy())
        };
        let tz = if line.dz() == T::zero() {
            None
        } else {
            Some((point.z() - line.start.z) / line.dz())
        };

        let contains = match (tx, ty, tz) {
            (None, None, None) => {
                // Degenerate line (line has no length, essentially a single point)
                point.0 == line.start
            }
            (Some(t), None, None) => {
                // Horizontal line
                point.y() == line.start.y && point.z() == line.start.z && T::zero() <= t && t <= T::one()
            }
            (None, Some(t), None) => {
                // Vertical line in Y direction
                point.x() == line.start.x && point.z() == line.start.z && T::zero() <= t && t <= T::one()
            }
            (None, None, Some(t)) => {
                // Vertical line in Z direction
                point.x() == line.start.x && point.y() == line.start.y && T::zero() <= t && t <= T::one()
            }
            (Some(t_x), Some(t_y), None) => {
                // 2D line in XY plane
                (t_x - t_y).abs() <= T::epsilon() && T::zero() <= t_x && t_x <= T::one()
            }
            (Some(t_x), None, Some(t_z)) => {
                // 2D line in XZ plane
                (t_x - t_z).abs() <= T::epsilon() && T::zero() <= t_x && t_x <= T::one()
            }
            (None, Some(t_y), Some(t_z)) => {
                // 2D line in YZ plane
                (t_y - t_z).abs() <= T::epsilon() && T::zero() <= t_y && t_y <= T::one()
            }
            (Some(t_x), Some(t_y), Some(t_z)) => {
                // General 3D line
                (t_x - t_y).abs() <= T::epsilon()
                    && (t_y - t_z).abs() <= T::epsilon()
                    && T::zero() <= t_x
                    && t_x <= T::one()
            }
        };

        if contains {
            return true;
        }
    }
    false
}
