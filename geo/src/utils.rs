//! Internal utility functions, types, and data structures.

use geo_types::{Coord, CoordNum};

/// An iterator of `I1` or `I2`
pub enum EitherIter<I1, I2> {
    A(I1),
    B(I2),
}

impl<I1, I2> ExactSizeIterator for EitherIter<I1, I2>
where
    I1: ExactSizeIterator,
    I2: ExactSizeIterator<Item = I1::Item>,
{
    #[inline]
    fn len(&self) -> usize {
        match self {
            EitherIter::A(i1) => i1.len(),
            EitherIter::B(i2) => i2.len(),
        }
    }
}

impl<T, I1, I2> Iterator for EitherIter<I1, I2>
where
    I1: Iterator<Item = T>,
    I2: Iterator<Item = T>,
{
    type Item = T;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        match self {
            EitherIter::A(iter) => iter.next(),
            EitherIter::B(iter) => iter.next(),
        }
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        match self {
            EitherIter::A(iter) => iter.size_hint(),
            EitherIter::B(iter) => iter.size_hint(),
        }
    }
}

// The Rust standard library has `max` for `Ord`, but not for `PartialOrd`
pub fn partial_max<T: PartialOrd>(a: T, b: T) -> T {
    if a > b {
        a
    } else {
        b
    }
}

// The Rust standard library has `min` for `Ord`, but not for `PartialOrd`
pub fn partial_min<T: PartialOrd>(a: T, b: T) -> T {
    if a < b {
        a
    } else {
        b
    }
}

use std::cmp::Ordering;

/// Compare two coordinates lexicographically: first by the
/// x coordinate, and break ties with the y coordinate then z.
/// Expects none of coordinates to be uncomparable (eg. nan)
#[inline]
pub fn lex_cmp<T: CoordNum>(p: &Coord<T>, q: &Coord<T>) -> Ordering {
    p.x.partial_cmp(&q.x).unwrap()
    .then(
        p.y.partial_cmp(&q.y).unwrap()
        .then(p.z.partial_cmp(&q.z).unwrap())
    )
}

/// Compute index of the least point in slice. Comparison is
/// done using [`lex_cmp`].
///
/// Should only be called on a non-empty slice with no `nan`
/// coordinates.
pub fn least_index<T: CoordNum>(pts: &[Coord<T>]) -> usize {
    pts.iter()
        .enumerate()
        .min_by(|(_, p), (_, q)| lex_cmp(p, q))
        .unwrap()
        .0
}

#[cfg(test)]
mod test {
    use super::{partial_max, partial_min};

    #[test]
    fn test_partial_max() {
        assert_eq!(5, partial_max(5, 4));
        assert_eq!(5, partial_max(5, 5));
    }

    #[test]
    fn test_partial_min() {
        assert_eq!(4, partial_min(5, 4));
        assert_eq!(4, partial_min(4, 4));
    }
}
