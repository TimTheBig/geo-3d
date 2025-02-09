//! Line measurements like [`Bearing`] and [`Distance`] for the [`Euclidean`] metric space.

mod bearing;
pub use bearing::Bearing;

mod destination;
pub use destination::Destination;

mod distance;
pub use distance::Distance;

mod interpolate_point;
pub use interpolate_point::InterpolatePoint;

mod length;
pub use length::Length;

mod densify;
pub use densify::Densify;

pub mod metric_spaces;
pub use metric_spaces::Euclidean;
