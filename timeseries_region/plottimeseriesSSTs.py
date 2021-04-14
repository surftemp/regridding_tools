# -*- coding: utf-8 -*-

#    regridding_tools
#    Copyright (C) 2020  National Centre for Earth Observation (NCEO)
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""
module for extracting SST time series and then generating test plots of the time series

call on the command line using the same options as extract_timeseries.py, with the following additions:

--out_plot_path <path to html or svg file to contain plot>

"""

import datetime
import csv
from .maketimeseriesSSTs import createParser, dispatch

# plot uses visigoth so "pip install visigoth"
from visigoth import Diagram
from visigoth.charts import Line
from visigoth.common import Text
from visigoth.containers import Map, Sequence
from visigoth.map_layers import WMS, Geoplot
from visigoth.map_layers.geoplot import Multipolygon
from visigoth.utils.mapping.projections import Projections

def createTimeSeriesPlot(args):
    reader = csv.reader(open(args.out_path))
    line = 0
    lookup = {}
    data = []
    parsers = {
        "mean temperature kelvin": lambda x:float(x["mean temperature kelvin"]),
        "mean temperature deg C": lambda x: float(x["mean temperature deg C"]),
        "mean temperature standard_error": lambda x: float(x["mean temperature standard_error"]),
        "fraction of sea-ice-covered ocean": lambda x: float(x["fraction of sea-ice-covered ocean"]),
        "date": lambda x:datetime.datetime(int(x["year"]),int(x["month"]),int(x["day"]))
    }
    for row in reader:
        line += 1
        if line < 3:
            pass
        elif line == 3:
            print(str(row))
            for idx in range(len(row)):
                lookup[row[idx]] = idx
        else:
            raw = {}
            datapoint = {}
            for key in lookup:
                raw[key] = row[lookup[key]]
            for key in parsers:
                try:
                    datapoint[key] = parsers[key](raw)
                except:
                    pass
            data.append(datapoint)

    d = Diagram(fill="#EEEEEE",spacing=40)
    d.add(Text("Lon: %0.2f - %0.2f, Lat: %0.2f - %0.2f"%(args.lon_min,args.lon_max,args.lat_min,args.lat_max)))

    lat_range = args.lat_max - args.lat_min
    lon_range = args.lon_max - args.lon_min

    lat_min = max(args.lat_min - lat_range*10,-90)
    lat_max = min(args.lat_max + lat_range*10,90)
    lon_min = max(args.lon_min - lon_range*10,-180)
    lon_max = min(args.lon_max + lon_range*10,180)

    m = Map(zoom_to=2,width=512,boundaries=((lon_min,lat_min),(lon_max,lat_max)),projection=Projections.ESPG_4326)
    wm = Map(width=256, boundaries=((-180,-90),(180,90)), projection=Projections.ESPG_4326)

    lps = []

    lps.append((args.lon_min, args.lat_min))
    lps.append((args.lon_max, args.lat_min))
    lps.append((args.lon_max, args.lat_max))
    lps.append((args.lon_min, args.lat_max))

    poly = Multipolygon([[lps]],stroke_width=0,fill="#FF000080")

    gp = Geoplot(multipolys=[poly])
    m.addLayer(WMS())
    m.addLayer(gp)

    birdview = []
    birdview.append((lon_min, lat_min))
    birdview.append((lon_max, lat_min))
    birdview.append((lon_max, lat_max))
    birdview.append((lon_min, lat_max))

    print(lat_max)

    birdpoly = Multipolygon([[birdview]], stroke_width=0, fill="#FF000080")

    wm.addLayer(WMS())
    birdplot = Geoplot(multipolys=[birdpoly])
    wm.addLayer(birdplot)

    s = Sequence(orientation="horizontal")
    s.add(m).add(wm)
    d.add(s)

    for key in ["mean temperature kelvin","mean temperature deg C","mean temperature standard_error","fraction of sea-ice-covered ocean"]:
        if key in data[0]:
            al = Line(data, height=256, line_width=4, width=1024, y=key, x="date", smoothing=0.3, font_height=10)
            (ax, ay) = al.getAxes()
            ay.setFontHeight(10)
            ax.setFontHeight(10)

            if key == "fraction of sea-ice-covered ocean":
                ay.setMinValue(0)
                ay.setMaxValue(1.0)
            else:
                y_min = min(datapoint[key] for datapoint in data)
                y_max = max(datapoint[key] for datapoint in data)
                if y_max == y_min:
                    ay.setMinValue(y_max - 0.5)
                    ay.setMaxValue(y_max + 0.5)

            ax.labelfn = lambda dt:dt.strftime("%Y-%m-%d")
            d.add(al)

    out_path = args.out_plot_path
    if out_path.endswith(".svg"):
        content=d.draw(format="svg")
    else:
        content=d.draw(format="html", html_title="SST time series")

    open(out_path,"w").write(content)


if __name__ == '__main__':
    parser = createParser()
    parser.add_argument('--out_plot_path', default="/tmp/out.html",
                        help='The path in which to write the output time series plot (should be .svg or .html).')
    args = parser.parse_args()
    if not args.out_path.endswith(".csv"):
        raise Exception("Error, --out_path should specify a comma separated variable (.csv) file")
    dispatch(args)
    createTimeSeriesPlot(args)
    print("Plot written to ",args.out_plot_path)