use std::marker::PhantomData;

#[cfg(feature = "geo-types")]
use geo_types::{Coord, CoordNum};

use crate::Dimensions;

/// A trait for accessing data from a generic Coord.
///
/// Refer to [geo_types::Coord] for information about semantics and validity.
pub trait CoordTrait {
    /// The coordinate type of this geometry
    type T;

    /// Dimensions of the coordinate tuple
    fn dim(&self) -> Dimensions;

    /// Access the n'th (0-based) element of the CoordinateTuple.
    /// Returns `None` if `n >= DIMENSION`.
    ///
    /// See also [`nth_or_panic()`](Self::nth_or_panic) and [`nth_unchecked()`](Self::nth_unchecked).
    ///
    /// # Panics
    ///
    /// This method may panic if [`dim()`](Self::dim) does not correspond to
    /// the actual number of dimensions in this coordinate.
    fn nth(&self, n: usize) -> Option<Self::T> {
        if n < self.dim().size() {
            Some(self.nth_or_panic(n))
        } else {
            None
        }
    }

    /// x component of this coord.
    fn x(&self) -> Self::T;

    /// y component of this coord.
    fn y(&self) -> Self::T;

    /// z component of this coord.
    fn z(&self) -> Self::T;

    /// Returns a tuple that contains the x/horizontal, y/depth & z/vertical component of the coord.
    fn x_y_z(&self) -> (Self::T, Self::T, Self::T) {
        (self.x(), self.y(), self.z())
    }

    /// Returns a tuple that contains the x/horizontal & y/vertical component of the coord.
    fn x_y(&self) -> (Self::T, Self::T) {
        (self.x(), self.y())
    }

    /// Access the n'th (0-based) element of the CoordinateTuple.
    /// May panic if n >= DIMENSION.
    /// See also [`nth()`](Self::nth).
    fn nth_or_panic(&self, n: usize) -> Self::T;

    /// Access the n'th (0-based) element of the CoordinateTuple.
    /// May panic if n >= DIMENSION.
    ///
    /// See also [`nth()`](Self::nth), [`nth_or_panic()`](Self::nth_or_panic).
    ///
    /// You might want to override the default implementation of this method
    /// if you can provide a more efficient implementation.
    ///
    /// # Safety
    ///
    /// Though it may panic, the default implementation actually is safe. However, implementors
    /// are allowed to implement this method themselves with an unsafe implementation. See the
    /// individual implementations for more information on their own Safety considerations.
    unsafe fn nth_unchecked(&self, n: usize) -> Self::T {
        self.nth_or_panic(n)
    }
}

#[cfg(feature = "geo-types")]
impl<T: CoordNum> CoordTrait for Coord<T> {
    type T = T;

    fn nth_or_panic(&self, n: usize) -> Self::T {
        match n {
            0 => self.x(),
            1 => self.y(),
            2 => self.z(),
            _ => panic!("Coord only supports 3 dimensions"),
        }
    }

    fn dim(&self) -> Dimensions {
        Dimensions::Xyz
    }

    fn x(&self) -> Self::T {
        self.x
    }

    fn y(&self) -> Self::T {
        self.y
    }

    fn z(&self) -> Self::T {
        self.z
    }
}

#[cfg(feature = "geo-types")]
impl<T: CoordNum> CoordTrait for &Coord<T> {
    type T = T;

    fn nth_or_panic(&self, n: usize) -> Self::T {
        match n {
            0 => self.x(),
            1 => self.y(),
            2 => self.z(),
            _ => panic!("Coord only supports 3 dimensions"),
        }
    }

    fn dim(&self) -> Dimensions {
        Dimensions::Xyz
    }

    fn x(&self) -> Self::T {
        self.x
    }

    fn y(&self) -> Self::T {
        self.y
    }

    fn z(&self) -> Self::T {
        self.z
    }
}

impl<T: Copy> CoordTrait for (T, T, T) {
    type T = T;

    fn nth_or_panic(&self, n: usize) -> Self::T {
        match n {
            0 => self.x(),
            1 => self.y(),
            2 => self.z(),
            _ => panic!("(T, T, T) only supports 3 dimensions"),
        }
    }

    fn dim(&self) -> Dimensions {
        Dimensions::Xyz
    }

    fn x(&self) -> Self::T {
        self.0
    }

    fn y(&self) -> Self::T {
        self.1
    }

    fn z(&self) -> Self::T {
        self.2
    }
}

/// An empty struct that implements [CoordTrait].
///
/// This can be used as the `CoordType` of the `GeometryTrait` by implementations that don't have a
/// Coord concept
pub struct UnimplementedCoord<T>(PhantomData<T>);

impl<T> CoordTrait for UnimplementedCoord<T> {
    type T = T;

    fn dim(&self) -> Dimensions {
        unimplemented!()
    }

    fn nth_or_panic(&self, _n: usize) -> Self::T {
        unimplemented!()
    }

    fn x(&self) -> Self::T {
        unimplemented!()
    }

    fn y(&self) -> Self::T {
        unimplemented!()
    }

    fn z(&self) -> Self::T {
        unimplemented!()
    }
}
