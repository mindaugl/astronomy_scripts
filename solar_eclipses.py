"""
Calculate solar eclipse times - the times when the angle
between moon and sun centers from earth are the smallest.
Differentiate between Total/Hybrid, Partial and Annular eclipses.
Determine magnitude and obscurity.
"""

from numpy import arctan2, pi, sin, byte, arccos, sqrt, array, dot
from skyfield.api import load
from skyfield.constants import AU_KM, C_AUDAY, ERAD
from skyfield.functions import angle_between, length_of
from skyfield.searchlib import find_minima
from skyfield.relativity import add_aberration

SOLAR_ECLIPSES = ["Partial", "Total/Hybrid", "Annular"]


def cross_area(d, r, R):
    return (
        r * r * arccos((d * d + r * r - R * R) / (2 * d * r))
        + R * R * arccos((d * d + R * R - r * r) / (2 * d * R))
        - 1 / 2 * sqrt((-d + r + R) * (d + r - R) * (d - r + R) * (d + r + R))
    )


def line_sphere_intersection(unit_vec, sphere_center, sphere_radius):
    center_sq = dot(sphere_center, sphere_center)
    prod1 = dot(unit_vec, sphere_center)
    delta = prod1 * prod1 - center_sq + sphere_radius * sphere_radius
    if delta < 0:
        s_to_point = unit_vec * prod1
        unit_e_to_point = (s_to_point - sphere_center) / length_of(
            s_to_point - sphere_center
        )
        return sphere_center + unit_e_to_point * sphere_radius
    else:
        d = prod1 - sqrt(delta)
        return unit_vec * d


def magnitude_obscurity(e_to_s, e_to_m, s_radius, m_radius, e_radius):
    s_to_m = e_to_m - e_to_s
    s_to_ecl = line_sphere_intersection(s_to_m / length_of(s_to_m), -e_to_s, e_radius)
    m_to_ecl = s_to_ecl - s_to_m
    scale = length_of(s_to_ecl) / length_of(m_to_ecl)
    s_radius_res = s_radius / scale
    s_to_ecl_res = s_to_ecl / scale
    dist = length_of(m_to_ecl - s_to_ecl_res)
    if dist >= s_radius_res + m_radius:
        return 0.0, 0.0
    if m_radius + dist <= s_radius_res or s_radius_res + dist <= m_radius:
        mag = m_radius / s_radius_res
        obs = 1.0
    else:
        mag = (s_radius_res + m_radius - dist) / (2 * s_radius_res)
        cross_a = cross_area(dist, s_radius_res, m_radius)
        obs = cross_a / (pi * s_radius_res * s_radius_res)
    return mag, obs


def f(t):
    """
    Taken from skyfield.eclipselib.lunar_eclipses()
    """
    jd, fr = t.whole, t.tdb_fraction
    b, velocity = earth_barycenter.compute_and_differentiate(jd, fr)
    e = earth.compute(jd, fr)
    m = moon.compute(jd, fr)
    s = sun.compute(jd, fr)

    earth_to_sun = s - b - e
    earth_to_moon = m - e

    # The aberration routine requires specific units.  (We can leave
    # the `earth_to_moon` vector unconverted because we only need
    # its direction.)  We approximate the Earth’s velocity as being
    # that of the Earth-Moon barycenter.

    earth_to_sun /= AU_KM
    velocity /= AU_KM
    light_travel_time = length_of(earth_to_sun) / C_AUDAY
    add_aberration(earth_to_sun, velocity, light_travel_time)

    return angle_between(earth_to_sun, earth_to_moon)


ts = load.timescale()
eph = load("de421.bsp")
# eph = load("de422.bsp")

start_time = ts.tt(1901, 1, 1)
end_time = ts.tt(2001, 1, 1)

solar_radius_km = 696340.0
moon_radius_km = 1737.1
earth_radius_km = ERAD / 1000

sdict = dict(((s.center, s.target), s.spk_segment) for s in eph.segments)
sun = sdict[0, 10]
earth_barycenter = sdict[0, 3]
earth = sdict[3, 399]
moon = sdict[3, 301]

f.step_days = 5.0
t, y = find_minima(start_time, end_time, f, num=4)

jd, fr = t.whole, t.tdb_fraction
b, velocity = earth_barycenter.compute_and_differentiate(jd, fr)
e = earth.compute(jd, fr)
m = moon.compute(jd, fr)
s = sun.compute(jd, fr)

earth_to_sun = s - b - e
moon_to_earth = e - m
moon_to_sun = s - b - m

earth_to_sun_AU = earth_to_sun / AU_KM
add_aberration(earth_to_sun_AU, velocity / AU_KM, length_of(earth_to_sun_AU) / C_AUDAY)
earth_to_sun = earth_to_sun_AU * AU_KM

print("angle:", angle_between(earth_to_sun, -moon_to_earth) * 180 / pi)
print("y:", y * 180 / pi)

###########

angle_penumbra = arctan2(moon_radius_km + solar_radius_km, length_of(moon_to_sun))
origin_pen_from_sun = (
    length_of(moon_to_sun) * solar_radius_km / (moon_radius_km + solar_radius_km)
)
origin_pen = s - origin_pen_from_sun * moon_to_sun / length_of(moon_to_sun)
pen_to_earth = b + e - origin_pen
angle_earth_from_pen = angle_between(-moon_to_sun, pen_to_earth)

earth_from_pen_cone = sin(angle_earth_from_pen - angle_penumbra) * length_of(
    pen_to_earth
)
is_eclipse = (angle_earth_from_pen < angle_penumbra) + (
    earth_from_pen_cone <= earth_radius_km
)

#################

origin_umbra_from_sun = (
    solar_radius_km * length_of(moon_to_sun) / (solar_radius_km - moon_radius_km)
)
angle_umbra = arctan2(solar_radius_km, origin_umbra_from_sun)
origin_umbra = s - origin_umbra_from_sun * moon_to_sun / length_of(moon_to_sun)
umbra_to_earth = b + e - origin_umbra

angle_earth_from_umbra = angle_between(moon_to_sun, umbra_to_earth)


is_angle_inside_umbra = angle_umbra - angle_earth_from_umbra > 0
is_umbra_inside_earth_radius = length_of(umbra_to_earth) <= earth_radius_km

angle_antumbra = pi - angle_umbra
earth_from_antumbra_cone = sin(angle_antumbra - angle_earth_from_umbra) * length_of(
    umbra_to_earth
)
is_eclipse_annular = (angle_earth_from_umbra - angle_antumbra > 0) + (
    earth_from_antumbra_cone <= earth_radius_km
)
is_eclipse_annular *= [not x for x in is_umbra_inside_earth_radius]

earth_from_umbra_cone = sin(angle_earth_from_umbra - angle_umbra) * length_of(
    umbra_to_earth
)
is_earth_radius_inside_umbra = earth_radius_km >= earth_from_umbra_cone

is_eclipse_total = is_eclipse * (
    is_angle_inside_umbra
    + ((is_earth_radius_inside_umbra) * (angle_earth_from_umbra - angle_umbra < pi / 2))
    + is_umbra_inside_earth_radius
)

times = t[is_eclipse]
code = is_eclipse_total.astype(byte)
code += 2 * (
    is_eclipse_annular.astype(byte)
    - (is_eclipse_annular * is_eclipse_total).astype(byte)
)
code = code[is_eclipse]

counts_by_type = dict(zip(SOLAR_ECLIPSES, [0, 0, 0]))
for i in range(len(times)):
    mag, obs = magnitude_obscurity(
        array([array(x[is_eclipse][i]) for x in earth_to_sun]),
        -array([array(x[is_eclipse][i]) for x in moon_to_earth]),
        solar_radius_km,
        moon_radius_km,
        earth_radius_km,
    )
    print(f"{times[i].tt_strftime()}, {SOLAR_ECLIPSES[code[i]]}, {mag:.5}, {obs:.3}")
    counts_by_type[SOLAR_ECLIPSES[code[i]]] += 1

print("Number of eclipses:", len(times))
print(counts_by_type)

