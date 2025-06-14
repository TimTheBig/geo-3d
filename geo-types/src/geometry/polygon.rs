use crate::{coord, Coord, CoordNum, Line, LineString, Rect, Triangle};
use alloc::vec;
use alloc::vec::Vec;

#[cfg(any(feature = "approx", test))]
use approx::{AbsDiffEq, RelativeEq};

/// A bounded three-dimensional area.
///
/// A `Polygon`’s outer boundary (_exterior ring_) is represented by a
/// [`LineString`]. It may contain zero or more holes (_interior rings_), also
/// represented by `LineString`s.
///
/// A `Polygon` can be created with the [`Polygon::new`] constructor or the [`polygon!`][`crate::polygon!`] macro.
///
/// # Semantics
///
/// The _boundary_ of the polygon is the union of the
/// boundaries of the exterior and interiors. The interior
/// is all the points inside the polygon (not on the
/// boundary).
///
/// The `Polygon` structure guarantees that all exterior and interior rings will
/// be _closed_, such that the first and last `Coord` of each ring has
/// the same value.
///
/// # Validity
///
/// - The exterior and interior rings must be valid
///   `LinearRing`s (see [`LineString`]).
///
/// - No two rings in the boundary may cross, and may
///   intersect at a `Point` only as a tangent. In other
///   words, the rings must be distinct, and for every pair of
///   common points in two of the rings, there must be a
///   neighborhood (a topological open set) around one that
///   does not contain the other point.
///
/// - The closure of the interior of the `Polygon` must
///   equal the `Polygon` itself. For instance, the exterior
///   may not contain a spike.
///
/// - The interior of the polygon must be a connected
///   point-set. That is, any two distinct points in the
///   interior must admit a curve between these two that lies
///   in the interior.
///
/// Refer to section 6.1.11.1 of the OGC-SFA for a formal
/// definition of validity. Besides the closed `LineString`
/// guarantee, the `Polygon` structure does not enforce
/// validity at this time. For example, it is possible to
/// construct a `Polygon` that has:
///
/// - fewer than 3 coordinates per `LineString` ring
/// - interior rings that intersect other interior rings
/// - interior rings that extend beyond the exterior ring
///
/// # `LineString` closing operation
///
/// Some APIs on `Polygon` result in a closing operation on a `LineString`. The
/// operation is as follows:
///
/// If a `LineString`’s first and last `Coord` have different values, a
/// new `Coord` will be appended to the `LineString` with a value equal to
/// the first `Coord`.
///
/// [`LineString`]: line_string/struct.LineString.html
#[derive(Eq, PartialEq, Clone, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Polygon<T: CoordNum = f64> {
    exterior: LineString<T>,
    interiors: Vec<LineString<T>>,
}

impl<T: CoordNum> Polygon<T> {
    /// Create a new `Polygon` with the provided exterior `LineString` ring and
    /// interior `LineString` rings.
    ///
    /// Upon calling `new`, the exterior and interior `LineString` rings [will
    /// be closed].
    ///
    /// [will be closed]: #linestring-closing-operation
    ///
    /// # Examples
    ///
    /// Creating a `Polygon` with no interior rings:
    ///
    /// ```
    /// use geo_3d_types::{LineString, Polygon};
    ///
    /// let polygon = Polygon::new(
    ///     LineString::from(vec![(0., 0., 0.), (1., 1., 1.), (1., 0., 1.), (0., 0., 1.)]),
    ///     vec![],
    /// );
    /// ```
    ///
    /// Creating a `Polygon` with an interior ring:
    ///
    /// ```
    /// use geo_3d_types::{LineString, Polygon};
    ///
    /// let polygon = Polygon::new(
    ///     LineString::from(vec![(0., 0., 0.), (1., 1., 1.), (1., 0., 1.), (0., 0., 1.)]),
    ///     vec![LineString::from(vec![
    ///         (0.1, 0.1, 0.1),
    ///         (0.9, 0.9, 0.9),
    ///         (0.9, 0.1, 0.9),
    ///         (0.1, 0.1, 0.1),
    ///     ])],
    /// );
    /// ```
    ///
    /// If the first and last `Coord`s of the exterior or interior
    /// `LineString`s no longer match, those `LineString`s [will be closed]:
    ///
    /// ```
    /// use geo_3d_types::{coord, LineString, Polygon};
    ///
    /// let mut polygon = Polygon::new(LineString::from(vec![(0., 0., 0.), (1., 1., 1.), (1., 0., 1.), (0., 0., 1.)]), vec![]);
    ///
    /// assert_eq!(
    ///     polygon.exterior(),
    ///     &LineString::from(vec![(0., 0., 0.), (1., 1., 1.), (1., 0., 1.), (0., 0., 1.)])
    /// );
    /// ```
    pub fn new(mut exterior: LineString<T>, mut interiors: Vec<LineString<T>>) -> Self {
        exterior.close();
        for interior in &mut interiors {
            interior.close();
        }
        Self {
            exterior,
            interiors,
        }
    }

    /// Returns the rings of the polygon
    pub fn rings(&self) -> impl Iterator<Item = &LineString<T>> {
        std::iter::once(self.exterior()).chain(self.interiors())
    }

    /// Consume the `Polygon`, returning the exterior `LineString` ring and
    /// a vector of the interior `LineString` rings.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_3d_types::{LineString, Polygon};
    ///
    /// let mut polygon = Polygon::new(
    ///     LineString::from(vec![(0., 0., 0.), (1., 1., 1.), (1., 0., 1.), (0., 0., 1.)]),
    ///     vec![LineString::from(vec![
    ///         (0.1, 0.1, 0.1),
    ///         (0.9, 0.9, 0.9),
    ///         (0.9, 0.1, 0.9),
    ///         (0.1, 0.1, 0.1),
    ///     ])],
    /// );
    ///
    /// let (exterior, interiors) = polygon.into_inner();
    ///
    /// assert_eq!(
    ///     exterior,
    ///     LineString::from(vec![(0., 0., 0.), (1., 1., 1.), (1., 0., 1.), (0., 0., 1.)])
    /// );
    ///
    /// assert_eq!(
    ///     interiors,
    ///     vec![LineString::from(vec![
    ///         (0.1, 0.1, 0.1),
    ///         (0.9, 0.9, 0.9),
    ///         (0.9, 0.1, 0.9),
    ///         (0.1, 0.1, 0.1),
    ///     ])]
    /// );
    /// ```
    pub fn into_inner(self) -> (LineString<T>, Vec<LineString<T>>) {
        (self.exterior, self.interiors)
    }

    /// Return a reference to the exterior `LineString` ring.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_3d_types::{LineString, Polygon};
    ///
    /// let exterior = LineString::from(vec![(0., 0., 0.), (1., 1., 1.), (1., 0., 1.), (0., 0., 1.)]);
    ///
    /// let polygon = Polygon::new(exterior.clone(), vec![]);
    ///
    /// assert_eq!(polygon.exterior(), &exterior);
    /// ```
    pub const fn exterior(&self) -> &LineString<T> {
        &self.exterior
    }

    /// Take the exterior `LineString` ring leving it empty.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_3d_types::{LineString, Polygon};
    ///
    /// let exterior = LineString::from(vec![(0., 0., 0.), (1., 1., 1.), (1., 0., 1.), (0., 0., 1.)]);
    ///
    /// let polygon = Polygon::new(exterior.clone(), vec![]);
    ///
    /// assert_eq!(polygon.exterior(), &exterior);
    /// ```
    pub fn take_exterior(&mut self) -> LineString<T> {
        LineString(core::mem::take(&mut self.exterior.0))
    }

    /// Execute the provided closure `f`, which is provided with a mutable
    /// reference to the exterior `LineString` ring.
    ///
    /// After the closure executes, the exterior `LineString` [will be closed].
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_3d_types::{coord, LineString, Polygon};
    ///
    /// let mut polygon = Polygon::new(
    ///     LineString::from(vec![(0., 0., 0.), (1., 1., 1.), (1., 0., 1.), (0., 0., 1.)]),
    ///     vec![],
    /// );
    ///
    /// polygon.exterior_mut(|exterior| {
    ///     exterior.0[1] = coord! { x: 1., y: 2., z: 1. };
    /// });
    ///
    /// assert_eq!(
    ///     polygon.exterior(),
    ///     &LineString::from(vec![(0., 0., 0.), (1., 2., 1.), (1., 0., 1.), (0., 0., 0.),])
    /// );
    /// ```
    ///
    /// If the first and last `Coord`s of the exterior `LineString` no
    /// longer match, the `LineString` [will be closed]:
    ///
    /// ```
    /// use geo_3d_types::{coord, LineString, Polygon};
    ///
    /// let mut polygon = Polygon::new(
    ///     LineString::from(vec![(0., 0., 0.), (1., 1., 1.), (1., 0., 1.), (0., 0., 1.)]),
    ///     vec![],
    /// );
    ///
    /// polygon.exterior_mut(|exterior| {
    ///     exterior.0[0] = coord! { x: 0., y: 1., z: 0. };
    /// });
    ///
    /// assert_eq!(
    ///     polygon.exterior(),
    ///     &LineString::from(vec![(0., 1., 0.), (1., 1., 1.), (1., 0., 1.), (0., 0., 0.), (0., 1., 0.),])
    /// );
    /// ```
    ///
    /// [will be closed]: #linestring-closing-operation
    pub fn exterior_mut<F>(&mut self, f: F)
    where
        F: FnOnce(&mut LineString<T>),
    {
        f(&mut self.exterior);
        self.exterior.close();
    }

    /// Fallible alternative to [`exterior_mut`](Polygon::exterior_mut).
    pub fn try_exterior_mut<F, E>(&mut self, f: F) -> Result<(), E>
    where
        F: FnOnce(&mut LineString<T>) -> Result<(), E>,
    {
        f(&mut self.exterior)?;
        self.exterior.close();
        Ok(())
    }

    /// Return a slice of the interior `LineString` rings.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_3d_types::{coord, LineString, Polygon};
    ///
    /// let interiors = vec![LineString::from(vec![
    ///     (0.1, 0.1, 0.1),
    ///     (0.9, 0.9, 0.9),
    ///     (0.9, 0.1, 0.9),
    ///     (0.1, 0.1, 0.1),
    /// ])];
    ///
    /// let polygon = Polygon::new(
    ///     LineString::from(vec![(0., 0., 0.), (1., 1., 1.), (1., 0., 1.), (0., 0., 1.)]),
    ///     interiors.clone(),
    /// );
    ///
    /// assert_eq!(interiors, polygon.interiors());
    /// ```
    pub fn interiors(&self) -> &[LineString<T>] {
        &self.interiors
    }

    /// Take the interior `LineString` rings leaving them empty.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_3d_types::{coord, LineString, Polygon};
    ///
    /// let interiors = vec![LineString::from(vec![
    ///     (0.1, 0.1, 0.1),
    ///     (0.9, 0.9, 0.9),
    ///     (0.9, 0.1, 0.9),
    ///     (0.1, 0.1, 0.1),
    /// ])];
    ///
    /// let mut polygon = Polygon::new(
    ///     LineString::from(vec![(0., 0., 0.), (1., 1., 1.), (1., 0., 1.), (0., 0., 1.)]),
    ///     interiors.clone(),
    /// );
    ///
    /// assert_eq!(interiors, polygon.take_interiors());
    /// assert!(polygon.interiors().is_empty());
    /// ```
    pub fn take_interiors(&mut self) -> Vec<LineString<T>> {
        core::mem::take(&mut self.interiors)
    }

    /// Execute the provided closure `f`, which is provided with a mutable
    /// reference to the interior `LineString` rings.
    ///
    /// After the closure executes, each of the interior `LineString`s [will be
    /// closed].
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_3d_types::{coord, LineString, Polygon};
    ///
    /// let mut polygon = Polygon::new(
    ///     LineString::from(vec![(0., 0., 0.), (1., 1., 1.), (1., 0., 1.), (0., 0., 1.)]),
    ///     vec![LineString::from(vec![
    ///         (0.1, 0.1, 0.1),
    ///         (0.9, 0.9, 0.9),
    ///         (0.9, 0.1, 0.9),
    ///         (0.1, 0.1, 0.1),
    ///     ])],
    /// );
    ///
    /// polygon.interiors_mut(|interiors| {
    ///     interiors[0].0[1] = coord! { x: 0.8, y: 0.8, z: 0.8 };
    /// });
    ///
    /// assert_eq!(
    ///     polygon.interiors(),
    ///     &[LineString::from(vec![
    ///         (0.1, 0.1, 0.1),
    ///         (0.8, 0.8, 0.8),
    ///         (0.9, 0.1, 0.9),
    ///         (0.1, 0.1, 0.1),
    ///     ])]
    /// );
    /// ```
    ///
    /// If the first and last `Coord`s of any interior `LineString` no
    /// longer match, those `LineString`s [will be closed]:
    ///
    /// ```
    /// use geo_3d_types::{coord, LineString, Polygon};
    ///
    /// let mut polygon = Polygon::new(
    ///     LineString::from(vec![(0., 0., 0.), (1., 1., 1.), (1., 0., 1.), (0., 0., 1.)]),
    ///     vec![LineString::from(vec![
    ///         (0.1, 0.1, 0.1),
    ///         (0.9, 0.9, 0.9),
    ///         (0.9, 0.1, 0.9),
    ///         (0.1, 0.1, 0.1),
    ///     ])],
    /// );
    ///
    /// polygon.interiors_mut(|interiors| {
    ///     interiors[0].0[0] = coord! { x: 0.1, y: 0.2, z: 0.3 };
    /// });
    ///
    /// assert_eq!(
    ///     polygon.interiors(),
    ///     &[LineString::from(vec![
    ///         (0.1, 0.2, 0.1),
    ///         (0.9, 0.9, 0.9),
    ///         (0.9, 0.1, 0.9),
    ///         (0.1, 0.1, 0.1),
    ///         (0.1, 0.2, 0.1),
    ///     ])]
    /// );
    /// ```
    ///
    /// [will be closed]: #linestring-closing-operation
    pub fn interiors_mut<F>(&mut self, f: F)
    where
        F: FnOnce(&mut [LineString<T>]),
    {
        f(&mut self.interiors);
        for interior in &mut self.interiors {
            interior.close();
        }
    }

    /// Fallible alternative to [`interiors_mut`](Self::interiors_mut).
    pub fn try_interiors_mut<F, E>(&mut self, f: F) -> Result<(), E>
    where
        F: FnOnce(&mut [LineString<T>]) -> Result<(), E>,
    {
        f(&mut self.interiors)?;
        for interior in &mut self.interiors {
            interior.close();
        }
        Ok(())
    }

    /// Add an interior ring to the `Polygon`.
    ///
    /// The new `LineString` interior ring [will be closed]:
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_3d_types::{coord, LineString, Polygon};
    ///
    /// let mut polygon = Polygon::new(
    ///     LineString::from(vec![(0., 0., 0.), (1., 1., 1.), (1., 0., 1.), (0., 0., 1.)]),
    ///     vec![],
    /// );
    ///
    /// assert_eq!(polygon.interiors().len(), 0);
    ///
    /// polygon.interiors_push(vec![(0.1, 0.1, 0.1), (0.9, 0.9, 0.9), (0.9, 0.1, 0.9)]);
    ///
    /// assert_eq!(
    ///     polygon.interiors(),
    ///     &[LineString::from(vec![
    ///         (0.1, 0.1, 0.1),
    ///         (0.9, 0.9, 0.9),
    ///         (0.9, 0.1, 0.9),
    ///         (0.1, 0.1, 0.1),
    ///     ])]
    /// );
    /// ```
    ///
    /// [will be closed]: #linestring-closing-operation
    pub fn interiors_push(&mut self, new_interior: impl Into<LineString<T>>) {
        let mut new_interior = new_interior.into();
        new_interior.close();
        self.interiors.push(new_interior);
    }

    /// Count the total number of rings (interior and exterior) in the polygon
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_3d_types::{coord, LineString, Polygon};
    ///
    /// let polygon = Polygon::new(
    ///     LineString::from(vec![(0., 0., 0.), (1., 1., 1.), (1., 0., 1.), (0., 0., 1.)]),
    ///     vec![],
    /// );
    ///
    /// assert_eq!(polygon.num_rings(), 1);
    ///
    /// let polygon = Polygon::new(
    ///     LineString::from(vec![(0., 0., 0.), (1., 1., 1.), (1., 0., 1.), (0., 0., 1.)]),
    ///     vec![LineString::from(vec![(0.1, 0.1, 0.1), (0.9, 0.9, 0.9), (0.9, 0.1, 0.9)])],
    /// );
    ///
    /// assert_eq!(polygon.num_rings(), 2);
    /// ```
    pub fn num_rings(&self) -> usize {
        self.num_interior_rings() + 1
    }

    /// Count the number of interior rings in the polygon
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_3d_types::{coord, LineString, Polygon};
    ///
    /// let polygon = Polygon::new(
    ///     LineString::from(vec![(0., 0., 0.), (1., 1., 1.), (1., 0., 1.), (0., 0., 1.)]),
    ///     vec![],
    /// );
    ///
    /// assert_eq!(polygon.num_interior_rings(), 0);
    ///
    /// let polygon = Polygon::new(
    ///     LineString::from(vec![(0., 0., 0.), (1., 1., 1.), (1., 0., 1.), (0., 0., 1.)]),
    ///     vec![LineString::from(vec![(0.1, 0.1, 0.1), (0.9, 0.9, 0.9), (0.9, 0.1, 0.9)])],
    /// );
    ///
    /// assert_eq!(polygon.num_interior_rings(), 1);
    /// ```
    pub fn num_interior_rings(&self) -> usize {
        self.interiors.len()
    }

    /// Makes a `Polygon` from a 2D shape, with a copy of the shape at the start and end heights
    pub fn from_2d_with_height(shape: &[(T, T)], start_height: T, end_height: T) -> Self {
        let exterior =
            // start height
            shape.iter().map(|p| coord!(p.0, p.1, start_height))
            // end height
            .chain(shape.iter().map(|p| coord!(p.0, p.1, end_height)))
            .collect();
        Polygon::new(exterior, vec![])
    }

    // burglarized from [`csgrs`](https://github.com/timschmidt/csgrs)
    /// Return an iterator over paired vertices each forming an edge of the polygon.
    pub fn exterior_edges(&self) -> impl Iterator<Item=(&Coord<T>, &Coord<T>)> {
        self.exterior.0.iter().zip(self.exterior.0.iter().cycle().skip(1))
    }

    // burglarized from [`csgrs`](https://github.com/timschmidt/csgrs)
    /// Return an iterator over paired vertices each forming an edge of the polygon.\
    /// Note: This copies the `Coord`s if you only need a reference to the coords try [`Polygon::exterior_edges`].
    ///
    /// # Examples
    ///
    /// ```
    /// # use geo_3d_types::{coord, Polygon, Line, LineString};
    /// #
    /// let polygon = Polygon::new(
    ///     LineString::from(vec![(0., 0., 0.), (1., 1., 1.), (1., 0., 1.), (0., 0., 1.)]),
    ///     vec![],
    /// );
    ///
    /// let expected_edges = [
    ///     Line::new(coord!(0., 0., 0.), coord!(1., 1., 1.)),
    ///     Line::new(coord!(1., 1., 1.), coord!(1., 0., 1.)),
    ///     Line::new(coord!(1., 0., 1.), coord!(0., 0., 1.)),
    ///     Line::new(coord!(0., 0., 1.), coord!(0., 0., 0.)),
    ///     Line::new(coord!(0., 0., 0.), coord!(0., 0., 0.)),
    /// ];
    ///
    /// # assert_eq!(polygon.exterior_lines().collect::<Vec<_>>(), expected_edges);
    /// ```
    /// ```ignore
    ///
    /// assert_eq!(polygon.exterior_lines(), expected_edges.into_iter());
    /// ```
    pub fn exterior_lines(&self) -> impl Iterator<Item=Line<T>> {
        self.exterior.0.iter().zip(self.exterior.0.iter().cycle().skip(1))
            .map(|(start, end)| Line::new(*start, *end))
    }

    /// How many coords a `Polygon` contains in all rings
    /// # Examples
    ///
    /// ```
    /// # use geo_3d_types::{coord, Polygon, Line, LineString};
    /// #
    /// let polygon = Polygon::new(
    ///     LineString::from(vec![(0., 0., 0.), (1., 1., 1.), (1., 0., 1.), (0., 0., 1.), (0., 0., 0.)]),
    /// #    vec![],
    /// );
    ///
    /// assert_eq!(polygon.coords_count(), 5);
    /// ```
    pub fn coords_count(&self) -> usize {
        self.rings().flatten().count()
    }
}

// impl<T: CoordNum> From<Rect<T>> for Polygon<T> {
//     fn from(r: Rect<T>) -> Self {
//         // todo find z min/max order
//         todo!("find z min/max order");
//         Polygon::new(
//             vec![
//                 (r.min().x, r.min().y),
//                 (r.max().x, r.min().y),
//                 (r.max().x, r.max().y),
//                 (r.min().x, r.max().y),
//                 (r.min().x, r.min().y),
//             ]
//             .into(),
//             Vec::new(),
//         )
//     }
// }
// todo check
impl<T: CoordNum> From<Rect<T>> for Polygon<T> {
    fn from(r: Rect<T>) -> Self {
        // Determine the z-order: min z and max z
        let z_min = r.min().z;
        let z_max = r.max().z;

        // Construct the 3D polygon. For a rectangle in 3D space,
        // we define two "levels" (z_min and z_max) and then close the shape.
        let exterior_coords = vec![
            // Bottom face (z_min)
            (r.min().x, r.min().y, z_min),
            (r.max().x, r.min().y, z_min),
            (r.max().x, r.max().y, z_min),
            (r.min().x, r.max().y, z_min),
            (r.min().x, r.min().y, z_min),
            // Top face (z_max)
            (r.min().x, r.min().y, z_max),
            (r.max().x, r.min().y, z_max),
            (r.max().x, r.max().y, z_max),
            (r.min().x, r.max().y, z_max),
            (r.min().x, r.min().y, z_max),
        ];

        // Convert exterior coordinates to polygon
        Polygon::new(exterior_coords.into(), Vec::new())
    }
}

impl<T: CoordNum> From<Triangle<T>> for Polygon<T> {
    fn from(t: Triangle<T>) -> Self {
        Polygon::new(vec![t.0, t.1, t.2, t.0].into(), Vec::new())
    }
}

#[cfg(any(feature = "approx", test))]
impl<T> RelativeEq for Polygon<T>
where
    T: AbsDiffEq<Epsilon = T> + CoordNum + RelativeEq,
{
    #[inline]
    fn default_max_relative() -> Self::Epsilon {
        T::default_max_relative()
    }

    /// Equality assertion within a relative limit.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_3d_types::{Polygon, polygon};
    ///
    /// let a: Polygon<f32> = polygon![(x: 0., y: 0.), (x: 5., y: 0.), (x: 7., y: 9.), (x: 0., y: 0.)];
    /// let b: Polygon<f32> = polygon![(x: 0., y: 0.), (x: 5., y: 0.), (x: 7.01, y: 9.), (x: 0., y: 0.)];
    ///
    /// approx::assert_relative_eq!(a, b, max_relative=0.1);
    /// approx::assert_relative_ne!(a, b, max_relative=0.001);
    /// ```
    ///
    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        if !self
            .exterior
            .relative_eq(&other.exterior, epsilon, max_relative)
        {
            return false;
        }

        if self.interiors.len() != other.interiors.len() {
            return false;
        }
        let mut zipper = self.interiors.iter().zip(other.interiors.iter());
        zipper.all(|(lhs, rhs)| lhs.relative_eq(rhs, epsilon, max_relative))
    }
}

#[cfg(any(feature = "approx", test))]
impl<T: AbsDiffEq<Epsilon = T> + CoordNum> AbsDiffEq for Polygon<T> {
    type Epsilon = T;

    #[inline]
    fn default_epsilon() -> Self::Epsilon {
        T::default_epsilon()
    }

    /// Equality assertion with an absolute limit.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo_3d_types::{Polygon, polygon};
    ///
    /// let a: Polygon<f32> = polygon![(x: 0., y: 0., z: 0.), (x: 5., y: 0., z: 5.), (x: 7., y: 9., z: 5.), (x: 0., y: 0., z: 0.01)];
    /// let b: Polygon<f32> = polygon![(x: 0., y: 0., z: 0.), (x: 5., y: 0., z: 5.), (x: 7.01, y: 9., z: 5.01), (x: 0., y: 0., z: 0.)];
    ///
    /// approx::assert_abs_diff_eq!(a, b, epsilon=0.1);
    /// approx::assert_abs_diff_ne!(a, b, epsilon=0.001);
    /// ```
    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        if !self.exterior.abs_diff_eq(&other.exterior, epsilon) {
            return false;
        }

        if self.interiors.len() != other.interiors.len() {
            return false;
        }
        let mut zipper = self.interiors.iter().zip(other.interiors.iter());
        zipper.all(|(lhs, rhs)| lhs.abs_diff_eq(rhs, epsilon))
    }
}

#[cfg(feature = "rstar")]
impl<T> rstar::RTreeObject for Polygon<T>
    where T: num_traits::Float + rstar::RTreeNum + Default
{
    type Envelope = rstar::AABB<super::Point<T>>;

    fn envelope(&self) -> Self::Envelope {
        self.exterior.envelope()
    }
}
