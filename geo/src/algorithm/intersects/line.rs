use super::{point_in_rect, Intersects};
use crate::*;

// todo check
impl<T: CoordNum> Intersects<Coord<T>> for Line<T> {
    fn intersects(&self, rhs: &Coord<T>) -> bool {
        let dist = (*rhs).distance(self);
        dist > T::zero() - T::epsilon() && dist < T::zero() + T::epsilon()
        // point_line_segment_distance(*rhs, self) == T::zero()
    }
}

// /// returns the distance that point is from the line segment
// pub fn point_line_segment_distance<T: CoordNum>(point: Coord<T>, line: &Line<T>) -> T {
//     let (l1, l2) = (line.start, line.end);
//     let dx = line.delta();
//     let m2 = dx.magnitude_squared();
//     // find parameter value of closest point on segment
//     let s12 = ((l2 - point).dot(dx) / m2).clamp(T::zero(), T::one());
//     // and find the distance
//     Euclidean::distance(point, l1 * s12 + l2 * (T::one() - s12))
// }

symmetric_intersects_impl!(Coord<T>, Line<T>);
symmetric_intersects_impl!(Line<T>, Point<T>);

// todo make 3d, with a way to get the point of intersection, try liner systems
impl<T: GeoNum> Intersects<Line<T>> for Line<T> {
    fn intersects(&self, line: &Line<T>) -> bool {
        // Special case: self is equiv. to a point.
        if self.start == self.end {
            return line.intersects(&self.start);
        }

        // Precondition: start and end are distinct.

        // Check if orientation of rhs.{start,end} are different
        // with respect to self.{start,end}.
        let check_1_1 = T::Ker::orient2d(self.start, self.end, line.start);
        let check_1_2 = T::Ker::orient2d(self.start, self.end, line.end);

        if check_1_1 != check_1_2 {
            // Since the checks are different,
            // rhs.{start,end} are distinct, and rhs is not
            // collinear with self. Thus, there is exactly
            // one point on the infinite extensions of rhs,
            // that is collinear with self.

            // By continuity, this point is not on the
            // exterior of rhs. Now, check the same with
            // self, rhs swapped.

            let check_2_1 = T::Ker::orient2d(line.start, line.end, self.start);
            let check_2_2 = T::Ker::orient2d(line.start, line.end, self.end);

            // By similar argument, there is (exactly) one
            // point on self, collinear with rhs. Thus,
            // those two have to be same, and lies (interior
            // or boundary, but not exterior) on both lines.
            check_2_1 != check_2_2
        } else if check_1_1 == Orientation::Collinear {
            // Special case: collinear line segments.

            // Equivalent to 4 point-line intersection
            // checks, but removes the calls to the kernel
            // predicates.
            point_in_rect(line.start, self.start, self.end)
                || point_in_rect(line.end, self.start, self.end)
                || point_in_rect(self.end, line.start, line.end)
                || point_in_rect(self.end, line.start, line.end)
        } else {
            false
        }
    }
}
