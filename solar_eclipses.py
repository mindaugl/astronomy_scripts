"""
Calculate solar eclipse (partial and full) times - the times when the angle
between moon and sun centers from earth are the smallest.
"""

from numpy import arctan2, pi
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

start_time = ts.utc(2018, 1, 1)
end_time = ts.utc(2022, 1, 1)

solar_radius_km = 696340.0
moon_radius_km = 1737.1

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

angle_penumbra = arctan2(moon_radius_km + solar_radius_km, length_of(moon_to_sun))
origin_pen_from_sun = (
    length_of(moon_to_sun) * solar_radius_km / (moon_radius_km + solar_radius_km)
)
origin_pen = s - origin_pen_from_sun * moon_to_sun / length_of(moon_to_sun)
pen_to_earth = b + e - origin_pen
angle_earth_from_pen = angle_between(-moon_to_sun, pen_to_earth)
angle_earth_radius_pen = arctan2(ERAD / 1000, length_of(pen_to_earth))

angle_diffs = angle_penumbra + angle_earth_radius_pen - angle_earth_from_pen
is_eclipse = angle_diffs > 0

print("solar eclipse times: ", t[is_eclipse].tt_strftime())
print("Solar eclipse angle differences [deg]: ", angle_diffs[is_eclipse] * 180 / pi)
