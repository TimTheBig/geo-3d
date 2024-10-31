mod i_overlay_integration;
#[cfg(test)]
mod tests;

pub use i_overlay_integration::BoolOpsNum;

use crate::geometry::{LineString, MultiLineString, MultiPolygon, Polygon};
use rstar::{ParentNode, RTree, RTreeNode, RTreeObject};

/// Boolean Operations on geometry.
///
/// Boolean operations are set operations on geometries considered as a subset
/// of the 2-D plane. The operations supported are: intersection, union, xor or
/// symmetric difference, and set-difference on pairs of 2-D geometries and
/// clipping a 1-D geometry with self.
///
/// These operations are implemented on [`Polygon`] and the [`MultiPolygon`]
/// geometries.
///
/// # Validity
///
/// Note that the operations are strictly well-defined only on *valid*
/// geometries. However, the implementation generally works well as long as the
/// interiors of polygons are contained within their corresponding exteriors.
///
/// Degenerate 2-d geoms with 0 area are handled, and ignored by the algorithm.
/// In particular, taking `union` with an empty geom should remove degeneracies
/// and fix invalid polygons as long the interior-exterior requirement above is
/// satisfied.
pub trait BooleanOps {
    type Scalar: BoolOpsNum;

    /// The exterior and interior rings of the geometry.
    ///
    /// It doesn't particularly matter which order they are in, as the topology algorithm counts crossings
    /// to determine the interior and exterior of the polygon.
    ///
    /// It is required that the rings are from valid geometries, that the rings not overlap.
    /// In the case of a MultiPolygon, this requires that none of its polygon's interiors may overlap.
    fn rings(&self) -> impl Iterator<Item = &LineString<Self::Scalar>>;

    fn boolean_op(
        &self,
        other: &impl BooleanOps<Scalar = Self::Scalar>,
        op: OpType,
    ) -> MultiPolygon<Self::Scalar> {
        use i_overlay::core::fill_rule::FillRule;
        use i_overlay::core::overlay::ShapeType;
        use i_overlay_integration::{convert, BoolOpsOverlay, BoolOpsOverlayGraph};
        let mut overlay = <Self::Scalar as BoolOpsNum>::OverlayType::new();

        for ring in self.rings() {
            overlay.add_path(convert::ring_to_shape_path(ring), ShapeType::Subject);
        }
        for ring in other.rings() {
            overlay.add_path(convert::ring_to_shape_path(ring), ShapeType::Clip);
        }

        let graph = overlay.into_graph(FillRule::EvenOdd);
        let shapes = graph.extract_shapes(op.into());

        convert::multi_polygon_from_shapes(shapes)
    }

    fn intersection(
        &self,
        other: &impl BooleanOps<Scalar = Self::Scalar>,
    ) -> MultiPolygon<Self::Scalar> {
        self.boolean_op(other, OpType::Intersection)
    }
    fn union(&self, other: &impl BooleanOps<Scalar = Self::Scalar>) -> MultiPolygon<Self::Scalar> {
        self.boolean_op(other, OpType::Union)
    }
    fn xor(&self, other: &impl BooleanOps<Scalar = Self::Scalar>) -> MultiPolygon<Self::Scalar> {
        self.boolean_op(other, OpType::Xor)
    }
    fn difference(
        &self,
        other: &impl BooleanOps<Scalar = Self::Scalar>,
    ) -> MultiPolygon<Self::Scalar> {
        self.boolean_op(other, OpType::Difference)
    }

    /// Clip a 1-D geometry with self.
    ///
    /// Returns the portion of `ls` that lies within `self` (known as the set-theoeretic
    /// intersection) if `invert` is false, and the difference (`ls - self`) otherwise.
    fn clip(
        &self,
        multi_line_string: &MultiLineString<Self::Scalar>,
        invert: bool,
    ) -> MultiLineString<Self::Scalar> {
        use i_overlay::core::fill_rule::FillRule;
        use i_overlay::string::clip::ClipRule;
        use i_overlay_integration::{convert, BoolOpsStringGraph, BoolOpsStringOverlay};

        let mut overlay = <Self::Scalar as BoolOpsNum>::StringOverlayType::new();

        for ring in self.rings() {
            overlay.add_shape_path(convert::ring_to_shape_path(ring));
        }
        for line_string in multi_line_string {
            for line in line_string.lines() {
                let line = [
                    Self::Scalar::to_bops_coord(line.start),
                    Self::Scalar::to_bops_coord(line.end),
                ];
                overlay.add_string_line(line)
            }
        }

        let graph = overlay.into_graph(FillRule::EvenOdd);
        let paths = graph.clip_string_lines(ClipRule {
            invert,
            boundary_included: true,
        });
        convert::multi_line_string_from_paths(paths)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum OpType {
    Intersection,
    Union,
    Difference,
    Xor,
}

/// Returns the [BooleanOps::union] of all contained geometries
///
/// This is more efficient than a manual union, but may result in increased memory use due to the use of an R*-tree
pub trait UnaryUnion {
    type Scalar: BoolOpsNum;

    /// Construct a tree of all the input geometries and progressively union them from the "bottom up"
    ///
    /// This is considerably more efficient than successively adding new [Polygon]s to an ever more complex [Polygon].
    /// The output [MultiPolygon] will contain a single member.
    fn unary_union(&self) -> MultiPolygon<Self::Scalar>;
}

fn bottom_up_fold_reduce<T, S, I, F, R>(
    tree: &RTree<T>,
    mut init: I,
    mut fold: F,
    mut reduce: R,
) -> S
where
    T: RTreeObject,
    I: FnMut() -> S,
    F: FnMut(S, &T) -> S,
    R: FnMut(S, S) -> S,
{
    fn inner<T, S, I, F, R>(parent: &ParentNode<T>, init: &mut I, fold: &mut F, reduce: &mut R) -> S
    where
        T: RTreeObject,
        I: FnMut() -> S,
        F: FnMut(S, &T) -> S,
        R: FnMut(S, S) -> S,
    {
        parent
            .children()
            .iter()
            .fold(init(), |accum, child| match child {
                RTreeNode::Leaf(value) => fold(accum, value),
                RTreeNode::Parent(parent) => {
                    let value = inner(parent, init, fold, reduce);

                    reduce(accum, value)
                }
            })
    }

    inner(tree.root(), &mut init, &mut fold, &mut reduce)
}

impl<T: BoolOpsNum> BooleanOps for Polygon<T> {
    type Scalar = T;

    fn rings(&self) -> impl Iterator<Item = &LineString<Self::Scalar>> {
        std::iter::once(self.exterior()).chain(self.interiors().iter())
    }
}

impl<T: BoolOpsNum> BooleanOps for MultiPolygon<T> {
    type Scalar = T;

    fn rings(&self) -> impl Iterator<Item = &LineString<Self::Scalar>> {
        self.0.iter().flat_map(|p| p.rings())
    }
}

/// Allows the unary union operation to be performed on any container which can produce items of type `Polygon<T>`
impl<T: BoolOpsNum, C> UnaryUnion for C
where
    C: IntoIterator<Item = Polygon<T>> + Clone,
    Polygon<T>: RTreeObject,
{
    type Scalar = T;

    fn unary_union(&self) -> MultiPolygon<Self::Scalar> {
        let init = || MultiPolygon::<T>::new(vec![]);

        let fold = |mut accum: MultiPolygon<T>, poly: &Polygon<T>| -> MultiPolygon<T> {
            accum = accum.union(poly);
            accum
        };

        let reduce = |accum1: MultiPolygon<T>, accum2: MultiPolygon<T>| -> MultiPolygon<T> {
            accum1.union(&accum2)
        };
        let rtree = RTree::bulk_load(self.clone().into_iter().collect());

        bottom_up_fold_reduce(&rtree, init, fold, reduce)
    }
}
