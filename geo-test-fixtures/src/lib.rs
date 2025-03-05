use std::{fs, iter::FromIterator, path::{Path, PathBuf}, str::FromStr};

use geo_types::{coord, Coord, CoordNum, LineString, MultiPolygon, Point, Polygon};
use num_traits::FloatConst;
use wkt::{Wkt, WktFloat};

pub fn louisiana<T>() -> LineString<T>
where
    T: WktFloat + Default + FromStr,
{
    line_string(Path::new("louisiana.wkt"))
}

pub fn baton_rouge<T>() -> Point<T>
where
    T: WktFloat + Default + FromStr,
{
    let x = T::from(-91.147385).unwrap();
    let y = T::from(30.471165).unwrap();
    let z = T::from(4.571367).unwrap();
    Point::new(x, y, z)
}

pub fn east_baton_rouge<T>() -> Polygon<T>
where
    T: WktFloat + Default + FromStr,
{
    polygon(Path::new("east_baton_rouge.wkt"))
}

pub fn norway_main<T>() -> LineString<T>
where
    T: WktFloat + Default + FromStr,
{
    line_string(Path::new("norway_main.wkt"))
}

pub fn norway_concave_hull<T>() -> LineString<T>
where
    T: WktFloat + Default + FromStr,
{
    line_string(Path::new("norway_concave_hull.wkt"))
}

pub fn norway_convex_hull<T>() -> LineString<T>
where
    T: WktFloat + Default + FromStr,
{
    line_string(Path::new("norway_convex_hull.wkt"))
}

pub fn norway_nonconvex_hull<T>() -> LineString<T>
where
    T: WktFloat + Default + FromStr,
{
    line_string(Path::new("norway_nonconvex_hull.wkt"))
}

pub fn vw_orig<T>() -> LineString<T>
where
    T: WktFloat + Default + FromStr,
{
    line_string(Path::new("vw_orig.wkt"))
}

pub fn vw_simplified<T>() -> LineString<T>
where
    T: WktFloat + Default + FromStr,
{
    line_string(Path::new("vw_simplified.wkt"))
}

pub fn poly1<T>() -> LineString<T>
where
    T: WktFloat + Default + FromStr,
{
    line_string(Path::new("poly1.wkt"))
}

pub fn poly1_hull<T>() -> LineString<T>
where
    T: WktFloat + Default + FromStr,
{
    line_string(Path::new("poly1_hull.wkt"))
}

pub fn poly2<T>() -> LineString<T>
where
    T: WktFloat + Default + FromStr,
{
    line_string(Path::new("poly2.wkt"))
}

pub fn poly2_hull<T>() -> LineString<T>
where
    T: WktFloat + Default + FromStr,
{
    line_string(Path::new("poly2_hull.wkt"))
}

pub fn poly_in_ring<T>() -> LineString<T>
where
    T: WktFloat + Default + FromStr,
{
    line_string(Path::new("poly_in_ring.wkt"))
}

pub fn sphere<T: CoordNum + FloatConst>() -> LineString<T> {
    fn sphere_points<T: CoordNum + FloatConst>(divisions: usize) -> LineString<T> {
        #[inline]
        fn rot_z<T: CoordNum>(point: Coord<T>, angle: T) -> Coord<T> {
            let e1 = angle.cos() * point.x - angle.sin() * point.y;
            let e2 = angle.sin() * point.x + angle.cos() * point.y;
            let e3 = point.z;
            coord!(e1, e2, e3)
        }

        #[inline]
        fn rot_x<T: CoordNum>(point: Coord<T>, angle: T) -> Coord<T> {
            let e1 = point.x;
            let e2 = angle.cos() * point.y - angle.sin() * point.z;
            let e3 = angle.sin() * point.y + angle.cos() * point.z;
            coord!(e1, e2, e3)
        }

        let mut points = Vec::with_capacity(divisions * divisions);
        let unit_y = coord!(T::zero(), T::one(), T::zero());
        for step_x in 0..divisions {
            let angle_x: T = (T::one() + T::one()) * T::PI() * (T::from(step_x).unwrap() / T::from(divisions).unwrap());
            let p = rot_x(unit_y, angle_x);
            for step_z in 0..divisions {
                let angle_z = (T::one() + T::one()) * T::PI() * (T::from(step_z).unwrap() / T::from(divisions).unwrap());
                let p = rot_z(p, angle_z);
                points.push(p);
            }
        }

        points.into()
    }
    sphere_points(104)
}

pub fn ring<T>() -> LineString<T>
where
    T: WktFloat + Default + FromStr,
{
    line_string(Path::new("ring.wkt"))
}

pub fn shell<T>() -> LineString<T>
where
    T: WktFloat + Default + FromStr,
{
    line_string(Path::new("shell.wkt"))
}

// From https://geodata.nationaalgeoregister.nl/kadastralekaart/wfs/v4_0?request=GetFeature&service=WFS&srsName=EPSG:4326&typeName=kadastralekaartv4:perceel&version=2.0.0&outputFormat=json&bbox=165593,480993,166125,481552
pub fn nl_zones<T>() -> MultiPolygon<T>
where
    T: WktFloat + Default + FromStr,
{
    multi_polygon(Path::new("nl_zones.wkt"))
}

// From https://afnemers.ruimtelijkeplannen.nl/afnemers/services?request=GetFeature&service=WFS&srsName=EPSG:4326&typeName=Enkelbestemming&version=2.0.0&bbox=165618,480983,166149,481542";
pub fn nl_plots_wgs84<T>() -> MultiPolygon<T>
where
    T: WktFloat + Default + FromStr,
{
    multi_polygon(Path::new("nl_plots.wkt"))
}

pub fn nl_plots_epsg_28992<T>() -> MultiPolygon<T>
where
    T: WktFloat + Default + FromStr,
{
    // https://epsg.io/28992
    multi_polygon(Path::new("nl_plots_epsg_28992.wkt"))
}

fn line_string<T>(name: &Path) -> LineString<T>
where
    T: WktFloat + Default + FromStr,
{
    let wkt = Wkt::from_str(&file(name)).unwrap();
    match wkt {
        Wkt::LineString(line_string) => wkt_line_string_to_geo(&line_string),
        _ => unreachable!(),
    }
}

pub fn polygon<T>(name: &Path) -> Polygon<T>
where
    T: WktFloat + Default + FromStr,
{
    let wkt = Wkt::from_str(&file(name)).unwrap();
    match wkt {
        Wkt::Polygon(wkt_polygon) => wkt_polygon_to_geo(&wkt_polygon),
        _ => unreachable!(),
    }
}

pub fn multi_polygon<T>(name: &Path) -> MultiPolygon<T>
where
    T: WktFloat + Default + FromStr,
{
    let wkt = Wkt::from_str(&file(name)).unwrap();
    match wkt {
        Wkt::MultiPolygon(multi_polygon) => wkt_multi_polygon_to_geo(&multi_polygon),
        _ => unreachable!(),
    }
}

fn wkt_line_string_to_geo<T>(line_string: &wkt::types::LineString<T>) -> LineString<T>
where
    T: WktFloat + Default + FromStr,
{
    LineString::from_iter(
        line_string
            .0.iter()
            .map(|coord| (coord.x, coord.y, coord.z)),
    )
}

fn wkt_polygon_to_geo<T>(polygon: &wkt::types::Polygon<T>) -> Polygon<T>
where
    T: WktFloat + Default + FromStr,
{
    let exterior: LineString<T> = wkt_line_string_to_geo(&polygon.0[0]);
    let interiors: Vec<LineString<T>> = polygon.0[1..].iter().map(wkt_line_string_to_geo).collect();

    Polygon::new(exterior, interiors)
}

fn wkt_multi_polygon_to_geo<T>(multi_polygon: &wkt::types::MultiPolygon<T>) -> MultiPolygon<T>
where
    T: WktFloat + Default + FromStr,
{
    let polygons: Vec<Polygon<T>> = multi_polygon.0.iter().map(wkt_polygon_to_geo).collect();
    MultiPolygon(polygons)
}

pub(crate) fn file(name: &Path) -> String {
    let mut res = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    res.push("fixtures");
    res.push(name);
    if let None = name.extension() {
        res.set_extension("wkt");
    }

    fs::read_to_string(res).unwrap()
}
