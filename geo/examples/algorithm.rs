use geo_3d::{line_string, Centroid};

fn main() {
    let linestring = line_string![
        (x: 40.02f64, y: 116.34, z: 1000.0),
        (x: 41.02f64, y: 116.34, z: -1000.1),
    ];
    println!("Centroid {:?}", linestring.centroid());
}
