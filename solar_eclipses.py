"""
Calculate solar eclipse times - the times when the angle
between moon and sun centers from earth are the smallest.
Differentiate between Total/Hybrid, Partial and Annular eclipses.
"""

from numpy import arctan2, pi, sin
from skyfield.api import load
from skyfield.constants import AU_KM, C_AUDAY, ERAD
from skyfield.functions import angle_between, length_of
from skyfield.searchlib import find_minima
from skyfield.relativity import add_aberration


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

start_time = ts.utc(2011, 1, 1)
end_time = ts.utc(2021, 1, 1)

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
b = earth_barycenter.compute(jd, fr)
e = earth.compute(jd, fr)
m = moon.compute(jd, fr)
s = sun.compute(jd, fr)

earth_to_sun = s - b - e
moon_to_earth = e - m
moon_to_sun = s - b - m

###########

angle_penumbra = arctan2(moon_radius_km + solar_radius_km, length_of(moon_to_sun))
origin_pen_from_sun = (
    length_of(moon_to_sun) * solar_radius_km / (moon_radius_km + solar_radius_km)
)
origin_pen = s - origin_pen_from_sun * moon_to_sun / length_of(moon_to_sun)
pen_to_earth = b + e - origin_pen
angle_earth_from_pen = angle_between(-moon_to_sun, pen_to_earth)
angle_earth_radius_pen = arctan2(earth_radius_km, length_of(pen_to_earth))

angle_diffs_pen = angle_penumbra + angle_earth_radius_pen - angle_earth_from_pen
is_eclipse_pen = angle_diffs_pen > 0

#################

origin_umbra_from_sun = (
    solar_radius_km * length_of(moon_to_sun) / (solar_radius_km - moon_radius_km)
)
angle_umbra = arctan2(solar_radius_km, origin_umbra_from_sun)
origin_umbra = s - origin_umbra_from_sun * moon_to_sun / length_of(moon_to_sun)
umbra_to_earth = b + e - origin_umbra

angle_earth_from_umbra = angle_between(moon_to_sun, umbra_to_earth)

distance_from_umbra_cone = sin(angle_earth_from_umbra - angle_umbra) * length_of(
    umbra_to_earth
)
is_angle_inside_umbra = angle_umbra - angle_earth_from_umbra > 0
is_earth_radius_inside_umbra = earth_radius_km >= distance_from_umbra_cone
is_umbra_inside_earth_radius = length_of(umbra_to_earth) <= earth_radius_km

angle_antumbra = pi - angle_umbra
distance_from_antumbra_cone = sin(angle_antumbra - angle_earth_from_umbra) * length_of(
    umbra_to_earth
)
is_eclipse_annular = (angle_earth_from_umbra - angle_antumbra > 0) + (
    distance_from_antumbra_cone <= earth_radius_km
)
is_eclipse_annular *= [not x for x in is_umbra_inside_earth_radius]

is_eclipse_total = is_eclipse_pen * (
    is_angle_inside_umbra
    + ((is_earth_radius_inside_umbra) * (angle_earth_from_umbra - angle_umbra < pi / 2))
    + is_umbra_inside_earth_radius
)

#################

for i in range(len(t)):
    if is_eclipse_total[i]:
        print(t[i].tt_strftime(), "Total/Hybrid")
    elif is_eclipse_annular[i]:
        print(t[i].tt_strftime(), "Annular")
    elif is_eclipse_pen[i]:
        print(t[i].tt_strftime(), "Partial")

