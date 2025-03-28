"""
Display RA/Dec of specific planet given start date, end date, and number date points
bokeh serve display_radec.py --show
http://localhost:5006/display_radec
"""

from numpy import arccos, arcsin, sin, cos, sqrt, pi, array, linspace, empty
from skyfield.api import load
from bokeh.plotting import figure, show
from bokeh.layouts import row, column
from bokeh.io import curdoc
from bokeh.models import ColumnDataSource, Slider, Select, TextInput, CustomJS, Dropdown, DatePicker

def update_data(att, old, new):
    def date_from_str(s):
        els = s.split('-')
        return ts.tt(int(els[0]), int(els[1]), int(els[2]))

    t_start = date_from_str(date_picker_s.value)
    t_end = date_from_str(date_picker_e.value)
    N = slider_n.value

    times = ts.tt_jd(linspace(t_start.tt, t_end.tt, N))
    try:
        radec_s = eph_s["earth"].at(times).observe(eph_s[planet]).radec()
    except KeyError:
        radec_s = eph_s["earth"].at(times).observe(eph_s[planet + "_barycenter"]).radec()

    ra = radec_s[0].radians * 12 / pi
    dec = radec_s[1].radians * 180 / pi
    colors = [(f"#{255:02x}{int(255 * (times.tt[-1] - t) / (times.tt[-1] - times.tt[0])):02x}"
               f"{255:02x}") for t in times.tt]
    source.data = dict(ra=ra, dec=dec, times_tt=times.tt, colors=colors)


ts = load.timescale()
planets = ['sun', 'mercury', 'moon', 'venus', 'mars', 'jupiter', 'saturn', 'uranus', 'neptune']
planet = 'sun'
eph_s_name = 'de421.bsp'
# eph_s_name = 'de422.bsp'
t_start = ts.tt(2000, 6, 21, 0, 0, 0)
t_end = ts.tt(2001, 6, 21, 0, 0, 0)
N = 100
eph_s = load(eph_s_name)

times = ts.tt_jd(linspace(t_start.tt, t_end.tt, N))
try:
    radec_s = eph_s["earth"].at(times).observe(eph_s[planet]).radec()
except KeyError:
    radec_s = eph_s["earth"].at(times).observe(eph_s[planet + "_barycenter"]).radec()

ra = radec_s[0].radians * 12/pi
dec = radec_s[1].radians * 180/pi
colors = [(f"#{255:02x}{int(255 * (times.tt[-1] - t) / (times.tt[-1] - times.tt[0])):02x}"
           f"{255:02x}") for t in times.tt]
source = ColumnDataSource(data=dict(ra=ra, dec=dec, times_tt=times.tt, colors=colors))
fig1 = figure(title="RA/Dec: " + planet.capitalize(), x_axis_label="RA(hours)",
              y_axis_label="DEC(deg)", width=1000, height=750)
fig1.scatter(x='ra', y='dec', source=source, fill_color='colors', line_color="red", size=4)

date_picker_s = DatePicker(
    title="Select start date",
    value=t_start.tt_strftime('%Y-%m-%d'),
    min_date="1899-07-29",
    max_date="2053-10-09",
)
date_picker_e = DatePicker(
    title="Select end date",
    value=t_end.tt_strftime('%Y-%m-%d'),
    min_date="1899-07-29",
    max_date="2053-10-09",
)
slider_n = Slider(title=f"Date points", value=60, start=5, end=1000, step=5)

date_pickers = [date_picker_s, date_picker_e]
widgets = date_pickers + [slider_n]
for w in widgets:
    w.on_change('value', update_data)

inputs = column(*widgets)
curdoc().add_root(row(inputs, fig1))
curdoc().title = f"RA/Dec"
