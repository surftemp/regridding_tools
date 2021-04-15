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

class RegionNetCDF4Formatter(object):
    """
    Create a formatter for writing data to a netcdf4 output file
    """

    def __init__(self,path="",time_resolution="daily",f_max=1.0,spatial_lambda=1.0,tau=3,output_anomaly=False,output_sea_ice=False,lon_min=0,lat_min=0,lon_max=5,lat_max=5,comment=""):
        """
        Construct the netcdf4 formatter using options

        :param path: the output path
        :param time_resolution: the time resolution
        :param f_max: sea ice threshold
        :param spatial_lambda: parameter for uncertainty calculation
        :param tau: parameter for uncertainty calculation
        :param output_anomaly: output anomaly rather than absolute SST
        :param output_sea_ice: include sea ice fraction in output
        :param lon_min:  minimum longitude of box, must be aligned on 0.05 degree boundary
        :param lat_min:  minimum latitude of box, must be aligned on 0.05 degree boundary
        :param lon_max:  maximum longitude of box, must be aligned on 0.05 degree boundary
        :param lat_max:  maximum latitude of box, must be aligned on 0.05 degree boundary
        :param comment: an extra string to add to the output file metadata
        """
        self.output_path = path
        self.time_resolution = time_resolution
        self.output_anomaly = output_anomaly
        self.output_sea_ice = output_sea_ice
        self.lon_min = lon_min
        self.lon_max = lon_max
        self.lat_min = lat_min
        self.lat_max = lat_max
        self.comment = comment
        self.f_max = f_max
        self.spatial_lambda = spatial_lambda
        self.tau = tau
        self.datapoints = []
        self.min_dt = None
        self.max_dt = None
        self.uuid = str(uuid4())

    def write(self,start_dt,mid_dt,end_dt,sst_or_anomaly,uncertainty,sea_ice_fraction):
        """
        Write an entry to the output file covering a time period
        :param start_dt: start date of the period
        :param mid_dt: mid date of the period
        :param end_dt: end date of the period
        :param sst_or_anomaly: the absolute SST or anomaly SST value
        :param uncertainty: the estimated uncertainty standard error
        :param sea_ice_fraction: the sea ice fraction
        """
        self.datapoints.append((start_dt,mid_dt,end_dt,sst_or_anomaly,uncertainty,sea_ice_fraction))

    def createField(self,values,name,long_name,units,scalar_int=False):
        """
        Create a cf.Field object with a time axis and metadata
        :param values: the array of values to form the field's data
        :param name: the field short name
        :param long_name: the field long name
        :param units: the field's measurement units
        :param scalar_int: True iff the field contains integers
        :return: a cf.Field object
        """

        properties = copy.deepcopy(FIELD_PROPERTIES)
        field = cf.Field(properties=properties)

        for key in GLOBAL_PROPERTIES:
            value = GLOBAL_PROPERTIES[key]
            field.nc_set_global_attribute(key,value)

        field.nc_set_global_attribute('f_max', self.f_max)
        field.nc_set_global_attribute('lambda', self.spatial_lambda)
        field.nc_set_global_attribute('tau', np.int32(self.tau))

        field.nc_set_variable(name)
        field.set_property("standard_name",long_name)
        field.long_name = long_name

        # Note - domain axes are constructed based on running .creation_commands() on a regridded cf.Field

        # domain axis - time
        c = cf.DomainAxis(size=len(values))
        c.nc_set_dimension('time')
        field.set_construct(c, key='domainaxis0')

        array_width = 0
        array_height = 0
        if not scalar_int:
            array_height = values[0].shape[0]
            array_width = values[0].shape[1]

            # domain_axis - latitude
            c = cf.DomainAxis(size=array_height)
            c.nc_set_dimension('lat')
            field.set_construct(c, key='domainaxis1')

            # domain_axis - longitude
            c = cf.DomainAxis(size=array_width)
            c.nc_set_dimension('lon')
            field.set_construct(c, key='domainaxis2')

        # Set the field data from the passed values
        if not scalar_int:
            values = list(map(lambda arr: arr.reshape(-1,array_height,array_width),values))
            arr = np.concatenate(values,axis=0) # values is a list of np array of the correct shape
        else:
            arr = np.array(values).reshape((-1)) # values is a list of scalar values, convert to an array
        m_arr = np.ma.masked_invalid(arr)
        data = cf.Data(m_arr, units=units, dtype='i4' if scalar_int else "f4")

        if not scalar_int:
            field.set_data(data, axes=('domainaxis0', 'domainaxis1', 'domainaxis2'))
        else:
            field.set_data(data, axes=('domainaxis0'))


        # Note - dimension coordinates below are constructed based on running .creation_commands() on a regridded cf.Field

        # dimension_coordinate - Time
        c = cf.DimensionCoordinate()
        c.set_properties({'long_name': 'reference time of sst field', 'standard_name': 'time', 'axis': 'T',
                          'units': 'seconds since 1981-01-01 00:00:00', 'calendar': 'gregorian', 'comment': ''})
        c.nc_set_variable('time')
        data = cf.Data([mdt for (_,mdt,_) in self.times], units='seconds since 1981-01-01 00:00:00', calendar='gregorian', dtype='i4')
        c.set_data(data)
        b = cf.Bounds()
        b.set_properties({'long_name': 'Time cell boundaries',
                          'comment': 'Contains the start and end times for the time period the data represent.',
                          'units': 'seconds since 1981-01-01 00:00:00', 'calendar': 'gregorian'})
        b.nc_set_variable('time_bnds')
        data = cf.Data([[sdt, edt] for (sdt,_,edt) in self.times], units='seconds since 1981-01-01 00:00:00', calendar='gregorian',
                       dtype='i4')
        b.set_data(data)
        c.set_bounds(b)
        field.set_construct(c, axes=('domainaxis0',), key='dimensioncoordinate0', copy=False)

        if not scalar_int:
            # dimension_coordinate - Latitude
            c = cf.DimensionCoordinate()
            c.set_properties(
                {'standard_name': 'latitude', 'long_name': 'latitude', 'units': 'degrees_north', 'valid_min': -90.0,
                 'valid_max': 90.0, 'axis': 'Y', 'comment': ' Latitude geographical coordinates,WGS84 projection'})
            c.nc_set_variable('lat')
            cell_height = (self.lat_max-self.lat_min)/array_height
            data = cf.Data([np.float32(self.lat_min + i*cell_height + cell_height/2.0) for i in range(0,array_height)], units='degrees_north', dtype='f4')
            c.set_data(data)
            b = cf.Bounds()
            b.set_properties({'units': 'degrees_north'})
            b.nc_set_variable('lat_bnds')
            data = cf.Data([[np.float32(self.lat_min+i*cell_height), np.float32(self.lat_min+(i+1)*cell_height)] for i in range(0,array_height)], units='degrees_north', dtype='f4')
            b.set_data(data)
            c.set_bounds(b)
            field.set_construct(c, axes=('domainaxis1',), key='dimensioncoordinate1', copy=False)

            # dimension_coordinate - Longitude
            c = cf.DimensionCoordinate()
            c.set_properties(
                {'standard_name': 'longitude', 'long_name': 'longitude', 'units': 'degrees_east', 'valid_min': -180.0,
                 'valid_max': 180.0, 'axis': 'X', 'comment': ' Longitude geographical coordinates,WGS84 projection'})
            c.nc_set_variable('lon')
            cell_width = (self.lon_max - self.lon_min) / array_width
            data = cf.Data([np.float32(self.lon_min + i*cell_width + cell_width/2.0) for i in range(0,array_width)], units='degrees_east', dtype='f4')
            c.set_data(data)
            b = cf.Bounds()
            b.set_properties({'units': 'degrees_east'})
            b.nc_set_variable('lon_bnds')
            data = cf.Data([[np.float32(self.lon_min+i*cell_width), np.float32(self.lon_min+(i+1)*cell_width)] for i in range(0,array_width)], units='degrees_east', dtype='f4')
            b.set_data(data)
            c.set_bounds(b)
            field.set_construct(c, axes=('domainaxis2',), key='dimensioncoordinate2', copy=False)


        comment = HEADER2a%({
            "time_resolution":self.time_resolution,
            "var_name": "sea surface temperature anomaly" if self.output_anomaly else "sea surface_temperature",
            "south_lat": self.lat_min,
            "north_lat": self.lat_max,
            "west_lon": self.lon_min,
            "east_lat": self.lon_max
        })

        if self.comment:
            comment += self.comment

        if self.comment:
            field.nc_set_global_attribute("comment",comment)

        return field

    def close(self):
        """Close the formatter and flush all changes to disk"""
        years = []
        months = []
        days = []
        doys = []
        uncertainties = []
        sea_ice_fractions = []
        self.times = []
        ssts_or_anomalies = []
        self.min_dt = self.datapoints[0][0]
        self.max_dt = self.datapoints[-1][0]

        # go through each of the collected data points and append the fields to arrays, ready to package into cf.Field
        for datapoint in self.datapoints:
            (sdt,dt,edt,sst_or_anomaly,uncertainty,sea_ice_fraction) = datapoint
            self.times.append((TimeSeriesUtils.seconds_since_1981(sdt),TimeSeriesUtils.seconds_since_1981(dt),TimeSeriesUtils.seconds_since_1981(edt)))
            years.append(dt.year)
            months.append(dt.month)
            days.append(dt.day)
            doys.append(dt.timetuple().tm_yday)
            ssts_or_anomalies.append(sst_or_anomaly)
            uncertainties.append(uncertainty)
            sea_ice_fractions.append(sea_ice_fraction)

        # create the fields that will be written to file
        year_field = self.createField(years,"calendar_year","calendar year","1",scalar_int=True)
        month_field = self.createField(months, "calendar_month", "calendar month","1",scalar_int=True)
        day_field = self.createField(days, "day_of_month", "day of month","1",scalar_int=True)
        doy_field = self.createField(doys, "day_of_year", "day of year","1",scalar_int=True)
        sst_or_anomaly_field = self.createField(ssts_or_anomalies,
                    "sst_anomaly" if self.output_anomaly else "sst",
                    "sea_water_temperature_anomaly" if self.output_anomaly else "sea_water_temperature","K")
        uncertainty_field = self.createField(uncertainties, "sst_uncertainty", "sea_water_temperature uncertainty","K")
        fields = [year_field,month_field,day_field,doy_field,sst_or_anomaly_field,uncertainty_field]
        if self.output_sea_ice:
            sea_ice_field = self.createField(sea_ice_fractions, "sea_ice_area_fraction", "sea_ice_area_fraction","K")
            fields.append(sea_ice_field)

        # write out to file
        fl = cf.FieldList(fields)
        cf.write(fl,self.output_path,datatype={np.dtype('float64'): np.dtype('float32')}, compress=1)
        fl.close()


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
        "comment":comment
    }

    formatter = RegionNetCDF4Formatter(**formatter_args)

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

