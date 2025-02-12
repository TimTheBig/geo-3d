use geo_3d::Point;
use geo_types::point;

fn main() {
    let p = point! {
        x: 40.02f64,
        y: 116.34,
        z: 450.9,
    };

    let Point(coord) = p;
    println!("Point at ({}, {}, {})", coord.x, coord.y, coord.z);
}
