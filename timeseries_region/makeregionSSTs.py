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
This module extracts regions of data from the L4 C3S/CCI dataset.  The time series:

* aggregates data over a number of different time granularities (daily, N-daily, pentads, dekads, monthly)
* region is defined using a spatial bounding box aligned on 0.05 degree boundaries and up to 5 degrees on each side
* can compute climatology-based anomalies or absolute SSTs
* allows the user to ignore input cells with sea ice fractions above a user specified value (f_max)
* outputs a a standard error based on user specified parameters tau and spatial lambda
* optionally, outputs the sea ice fraction of the ocean area
* outputs to netcdf4 or comma separated variable formats

This data accepts as input climatology and C3S/CCI L4 SST files that have been spatially resliced from the original
one-file-per-day format.  The resliced format stores data in zarr format, with 5 degree spatial and 7 day temporal chunking
and is efficient for the retrieval of time series data for a small area.

This module can be invoked from the command line, use --help to show the options.
This module can also be invoked by importing and calling the makeTimeSeriesSSTs function
"""

import cf
import os.path
import datetime
import numpy as np
import json
import copy
from uuid import uuid4
from utils import TimeSeriesUtils
from extractor import Extractor
from aggregator import Aggregator
from netcdf4formatter import NetCDF4Formatter

# define some defaults
_default_out_path = "/tmp/region.nc"

_default_start_date = datetime.datetime(2018,1,1)
_default_end_date = datetime.datetime(2018,12,31)

_default_in_path = "/group_workspaces/jasmin2/nceo_uor/niall/reslice"
_default_out_path = "timeseries.csv"

_default_f_max = 1.0
_default_tau = 3
_default_spatial_lambda = 1.0

_default_anomalies = False
_default_no_sea_ice_fraction = False

FIELD_PROPERTIES = json.loads(open(os.path.join(os.path.split(__file__)[0],"field_properties_template.json"),"r").read())
GLOBAL_PROPERTIES = json.loads(open(os.path.join(os.path.split(__file__)[0],"global_properties_template.json"),"r").read())

HEADER2a = "Region of %(time_resolution)s %(var_name)s averaged across ocean surface"+ \
    " within latitudes %(south_lat)0.2f and %(north_lat)0.2f and longitudes %(west_lon)0.2f and %(east_lat)0.2f," + \
    " in kelvin and degree Celsius. "
HEADER2b = "Areas where sea-ice concentration exceeded %(fmax_pct)d %% were excluded in the average temperature calculation. "

def makeRegionSSTs(lon_min,lon_max,lat_min,lat_max,time_resolution,
                       start_date=_default_start_date,
                       end_date=_default_end_date,
                       input_path=_default_in_path,out_path=_default_out_path,
                       tau=_default_tau,spatial_lambda=_default_spatial_lambda,f_max=_default_f_max,
                       anomalies=_default_anomalies,no_sea_ice_fraction=_default_no_sea_ice_fraction,
                       comment=""):
    """
    Obtain a time series from SST data.

    :param lon_min:
        The minimum longitude value of the spatial area to aggregate.  Must be aligned on 0.05 degree boundary.

    :param lon_max:
        The maximum longitude value of the spatial area to aggregate.  Must be aligned on 0.05 degree boundary.

    :param lat_min:
        The minimum latitude value of the spatial area to aggregate.  Must be aligned on 0.05 degree boundary.

    :param lat_max:
        The maximum latitude value of the spatial area to aggregate.  Must be aligned on 0.05 degree boundary.

    :param time_resolution:
        The time resolution to aggregate, one of "daily","pentad","dekad","monthly" or "N" where N is a number of days >= 1

    :param input_path:
        Path of the folder containing spatially resliced SST and climatology data providing the input data.

    :param output_path:
        Path to write the time series data.  Must be a filename, data will be written in CSV format if the file extension is .csv

    :param f_max:
        Ignore input cells with a fraction of sea ice higher than this

    :param spatial_lambda:
        Distance in degrees over which input data errors are assumed to be correlated

    :param tau:
        Time in days over which input data errors are assumed to be correlated

    :param anomalies:
        Output anomalies rather than absolute SSTs

    :param no_sea_ice_fraction:
        Do not output the aggregated sea ice fraction

    :param comment:
        An extra comment string to add to the output file metadata
    """

    output_sea_ice = not no_sea_ice_fraction

    # work out the output format from the file extension
    output_csv = False
    if out_path.endswith(".csv"):
        output_csv = True

    # create an extractor to read the relevant part of the input data covering the extraction times and spatial boundaries
    extractor = Extractor(input_path)

    # create an aggregator to aggregate each period in the extracted data
    aggregator = Aggregator(mode="region",lon_min=lon_min, lat_min=lat_min, lon_max=lon_max,
            lat_max=lat_max, f_max=f_max, output_sea_ice=output_sea_ice, output_anomaly=anomalies,
            spatial_lambda=spatial_lambda, tau=tau)

    # create a formatter (either CSV or netcdf4 based) to handle writing the aggregated data to file
    formatter_args = {
        "path":out_path,
        "f_max":f_max, "tau":tau, "spatial_lambda":spatial_lambda,
        "output_anomaly":anomalies,"output_sea_ice":output_sea_ice,
        "time_resolution":time_resolution,
        "lon_min":lon_min, "lat_min":lat_min,"lon_max":lon_max,"lat_max":lat_max,
        "comment":comment,
        "header":HEADER2a%({
            "time_resolution":time_resolution,
            "var_name": "sea surface temperature anomaly" if anomalies else "sea surface_temperature",
            "south_lat": lat_min,
            "north_lat": lat_max,
            "west_lon": lon_min,
            "east_lat": lon_max
        })
    }

    formatter = NetCDF4Formatter(**formatter_args)

    # if outputting anomalies, we'll need the climatology.  Extract it and pass it to the aggregator
    if anomalies:
        climatology = extractor.extractClimatology(min_lon=lon_min, min_lat=lat_min, max_lon=lon_max, max_lat=lat_max)
        aggregator.setClimatology(climatology)


    period_duration = end_date.timestamp() - start_date.timestamp()

    # loop over each time period in the required date range...
    for (dates,slice_data) in extractor.generateData(start_dt=start_date, end_dt=end_date, time_resolution=time_resolution,
                                                    min_lon=lon_min, min_lat=lat_min, max_lon=lon_max, max_lat=lat_max):

        # get the first,middle and end date of the period
        (s_dt,mid_dt,e_dt) = dates

        # aggregate this time period...
        (sst_or_anomaly,uncertainty,sea_ice_fraction) = aggregator.aggregate(start_dt=s_dt,end_dt=e_dt,data=slice_data)
        # print("slice:",mid_dt,sst_or_anomaly,uncertainty,sea_ice_fraction)

        # and append it to the output file
        formatter.write(s_dt,mid_dt,e_dt,sst_or_anomaly,uncertainty,sea_ice_fraction)

    formatter.close()

def createParser():
    import argparse
    parser = argparse.ArgumentParser(description='extract SST data time series.')

    parser.add_argument('lon_min', type=float,
                        help='The minimum longitude value in degrees. This must be a multiple ' +
                             'of 0.05 degrees and range between -180.00 and +179.95. It must ' +
                             'also be less than the maximum longitude value.')

    parser.add_argument('lon_max', type=float,
                        help='The maximum longitude value in degrees. This must be a multiple ' +
                             'of 0.05 degrees and range between -179.95 and +180.00. It must ' +
                             'also be greater than the minimum longitude value.')

    parser.add_argument('lat_min', type=float,
                        help='The minimum latitude value in degrees. This must be a multiple ' +
                             'of 0.05 degrees and range between -90.00 and +89.95. It must ' +
                             'also be less than the maximum latitude value.')

    parser.add_argument('lat_max', type=float,
                        help='The maximum latitude value in degrees. This must be a multiple ' +
                             'of 0.05 degrees and range between -89.95 and +90.00. It must ' +
                             'also be greater than the minimum latitude value.')

    parser.add_argument('time_resolution', help="The target time resolution. This can be 'monthly', 'daily'," +
                                                "'10-day' for dekads, '5-day' for pentads or an integer for regular " +
                                                " N day regridding aligned with the start of the year, or daily.")

    parser.add_argument('--start_year', type=int, default=_default_start_date.year,
                        help='The start year of the time series.')

    parser.add_argument('--start_month', type=int, default=_default_start_date.month,
                        help='The start month of the time series.')

    parser.add_argument('--start_day', type=int, default=_default_start_date.day,
                        help='The start day of the time series.')

    parser.add_argument('--end_year', type=int, default=_default_end_date.year,
                        help='The end year of the time series.')

    parser.add_argument('--end_month', type=int, default=_default_end_date.month,
                        help='The end month of the time series.')

    parser.add_argument('--end_day', type=int, default=_default_end_date.day,
                        help='The end day of the time series.')

    parser.add_argument('--in_path', default=_default_in_path,
                        help='Path to the resliced SST C3S/CCI Analysis Level 4 input data.')

    parser.add_argument('--out_path', default=_default_out_path,
                        help='The path in which to write the output time series (output will be written in CSV format if file extension is .csv, otherwise netcdf4).')

    parser.add_argument('--f_max', type=float, default=_default_f_max,
                        help='The fraction of sea ice above which a cell is ignored in calculating the regridded ' +
                             'SST. Default is 1.0. When calculating anomalies, the climatology is always calculated ' +
                             'with an f_max of 1.0 regardless of this value.')

    parser.add_argument('--tau', type=int, default=_default_tau,
                        help='Timescale within which errors are assumed to be fully correlated in days. ' +
                             'Default is 3 days.')

    parser.add_argument('--spatial_lambda', type=float, default=_default_spatial_lambda,
                        help='Spatial scale within which errors are assumed to be fully correlated in degrees. ' +
                             'Default is 1.0 degrees.')

    parser.add_argument('--anomalies', action='store_true', default=_default_anomalies,
                        help='Output anomalies instead of absolute SSTs.')

    parser.add_argument('--no_sea_ice_fraction', action='store_true', default=_default_no_sea_ice_fraction,
                        help='Do not output the sea ice fraction.')

    parser.add_argument('--comment', type=str, default="",
                        help='Supply an extra comment string to be added to the output file.')

    return parser

def dispatch(args):
    # assemble the start and end dates from the input/output year/month/day components
    start_dt = datetime.datetime(args.start_year, args.start_month, args.start_day, 12, 0, 0)
    end_dt = datetime.datetime(args.end_year, args.end_month, args.end_day, 12, 0, 0)

    makeRegionSSTs(lon_min=args.lon_min, lat_min=args.lat_min, lon_max=args.lon_max, lat_max=args.lat_max,
                       time_resolution=args.time_resolution, start_date=start_dt, end_date=end_dt,
                       out_path=args.out_path, input_path=args.in_path,
                       f_max=args.f_max, spatial_lambda=args.spatial_lambda, tau=args.tau,
                       anomalies=args.anomalies, no_sea_ice_fraction=args.no_sea_ice_fraction,
                       comment=args.comment)


if __name__ == '__main__':
    parser = createParser()
    args = parser.parse_args()
    dispatch(args)

