use crate::Coord;
use crate::GeoNum;

#[derive(Debug, Clone)]
pub(crate) struct Segment<F: GeoNum + rstar::RTreeNum> {
    pub edge_idx: usize,
    pub segment_idx: usize,
    pub envelope: rstar::AABB<Coord<F>>,
}

impl<F> Segment<F>
where
    F: GeoNum + rstar::RTreeNum,
{
    pub fn new(edge_idx: usize, segment_idx: usize, p1: Coord<F>, p2: Coord<F>) -> Self {
        Self {
            edge_idx,
            segment_idx,
            envelope: rstar::AABB::from_corners(p1, p2),
        }
    }
}

impl<F> rstar::RTreeObject for Segment<F>
where
    F: GeoNum + rstar::RTreeNum,
{
    type Envelope = rstar::AABB<Coord<F>>;

    fn envelope(&self) -> Self::Envelope {
        self.envelope
    }
}
