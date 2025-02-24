#![no_main]

use geo_3d::Contains;
use libfuzzer_sys::fuzz_target;

fuzz_target!(|tuple: (geo_types::Geometry<f32>, geo_types::Coord<f32>)| {
    let (geometry, point) = tuple;

    geometry.contains(point);
});
