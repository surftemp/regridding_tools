# -*- coding: utf-8 -*-

#    sst-services
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
This module extracts time series data from the L4 C3S/CCI dataset.  The time series:

* aggregates data over a number of different time granularities (daily, N-daily, pentads, dekads, monthly)
* aggregates data over a spatial bounding box aligned on 0.05 degree boundaries and up to 5 degrees on each side
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
import csv
import numpy as np
from math import cos, radians, sqrt, isnan, isinf, nan
from calendar import monthrange
import json
import copy
import sys
from uuid import uuid4
import zarr
import math

# define some defaults
_default_out_path = "/tmp/ts.csv"

_default_start_date = datetime.datetime(2018,1,1)
_default_end_date = datetime.datetime(2018,12,31)

_default_in_path = "/group_workspaces/jasmin2/nceo_uor/niall/reslice"
_default_out_path = "timeseries.csv"

_default_f_max = 1.0
_default_tau = 7
_default_spatial_lambda = 3.0

_default_anomalies = False
_default_no_sea_ice_fraction = False

SST_FIELD_NAMES = ['sea_water_temperature', 'sea_water_temperature standard_error', 'sea_ice_area_fraction']
CLIMATOLOGY_FIELD_NAMES = ['sea_water_temperature', 'sea_ice_area_fraction']

FIELD_PROPERTIES = json.loads(open(os.path.join(os.path.split(__file__)[0],"field_properties_template.json"),"r").read())
GLOBAL_PROPERTIES = json.loads(open(os.path.join(os.path.split(__file__)[0],"global_properties_template.json"),"r").read())

class Progress(object):

    def __init__(self,label):
        self.label = label
        self.last_progress_frac = None

    def report(self,msg,progress_frac):
        if self.last_progress_frac == None or (progress_frac - self.last_progress_frac) >= 0.01:
            self.last_progress_frac = progress_frac
            i = int(100*progress_frac)
            if i > 100:
                i = 100
            si = i // 2
            sys.stdout.write("\r%s %s %-05s %s" % (self.label,msg,str(i)+"%","#"*si))
            sys.stdout.flush()

    def complete(self,msg):
        sys.stdout.write("\n%s %s\n" % (self.label,msg))
        sys.stdout.flush()

    def clear(self):
        sys.stdout.write("\r")


class TimeSeriesUtils(object):
    """
    Implement some handy utility methods to help with time series extraction
    """

    # some metdata uses seconds since 1981
    EPOCH = datetime.datetime(1981, 1, 1)

    @staticmethod
    def seconds_since_1981(dt):
        """Compute the number of seconds since Jan 1st 1982"""
        return int((dt - TimeSeriesUtils.EPOCH).total_seconds())

    @staticmethod
    def createLatitudeWeighting(lat_min, lat_max):
        """Create an array of latitude area weighting values for the given latitude range"""
        height = round((lat_max - lat_min) / 0.05)
        arraydata = []
        for y in range(0, height):
            lat = lat_min + (0.05 * y) + 0.025
            weight = cos(radians(lat))
            arraydata.append(weight)
        # reshape so can be multiplied along the latitude access of an array organised by (time,latitude,longitude)
        return np.array(arraydata).reshape(-1,1)

    @staticmethod
    def kToDegC(k):
        """Convert from kelvin to degrees centigrade"""
        return k - 273.15

    @staticmethod
    def getDaysInYear(year):
        """Get the number of days in a year"""
        d1 = datetime.date(year, 1, 1)
        d2 = datetime.date(year, 12, 31)
        return 1 + (d2 - d1).days

    @staticmethod
    def lastDayInMonth(year, month):
        """Work out how many days there are in a given month"""
        return monthrange(year, month)[1]

class TimeSeriesAggregator(object):
    """
    Aggregate a 3 dimensional box (bounded by time, latitude and longitude) to yield sst, anomaly, uncertainty
    and sea ice fraction values according to the requirements
    """

    def __init__(self,lon_min,lat_min,lon_max,lat_max,f_max=1.0,output_anomaly=False,output_sea_ice=False,spatial_lambda=1.0,tau=3):
        """
        Construct the aggregator

        :param lon_min:  minimum longitude of box, must be aligned on 0.05 degree boundary
        :param lat_min:  minimum latitude of box, must be aligned on 0.05 degree boundary
        :param lon_max:  maximum longitude of box, must be aligned on 0.05 degree boundary
        :param lat_max:  maximum latitude of box, must be aligned on 0.05 degree boundary
        :param f_max: sea ice fraction threshold
        :param output_anomaly: output anomalies rather than absolute SST
        :param output_sea_ice: output sea ice fraction
        :param spatial_lambda: parameter for uncertainty calculation
        :param tau: parameter for uncertainty calculation
        """
        self.lon_min = lon_min
        self.lon_max = lon_max
        self.lat_min = lat_min
        self.lat_max = lat_max
        self.lat_weighting = None # TimeSeriesUtils.createLatitudeWeighting(lat_min, lat_max)
        self.f_max = f_max
        self.output_anomaly = output_anomaly
        self.output_sea_ice = output_sea_ice
        self.spatial_lambda = spatial_lambda
        self.tau = tau
        self.climatology = None


    def aggregate(self,start_dt,end_dt,data):
        """
        Perform aggregation on the data pertaining to a particular time period

        :param start_dt: the start date of the period (inclusive)
        :param end_dt: the end date of the period (inclusive)
        :param fieldset: a cf.FieldSet containing the daily data for the period and space at 0.05 degree resolution

        :return: tuple containing aggregted values (sst_or_anomaly,uncertainty,sea_ice_fraction)
        """
        sea_ice_fraction_index = SST_FIELD_NAMES.index("sea_ice_area_fraction")
        sea_water_temp_index = SST_FIELD_NAMES.index("sea_water_temperature")
        standard_error_index = SST_FIELD_NAMES.index("sea_water_temperature standard_error")

        sea_ice_fractions = np.nan_to_num(data[:,:,:,sea_ice_fraction_index],copy=False)
        sea_water_temp = np.nan_to_num(data[:, :, :, sea_water_temp_index],copy=False)
        standard_error = np.nan_to_num(data[:, :, :, standard_error_index],copy=False)
        land_mask = np.where(sea_water_temp > 0.0,1.0,0.0)

        if self.lat_weighting is None:
            self.lat_weighting = TimeSeriesUtils.createLatitudeWeighting(lat_min=self.lat_min,lat_max=self.lat_max)

        sea_ice_mask = np.where(sea_ice_fractions > self.f_max, 0.0, 1.0)


        mean_sea_water_temp_numerator = (sea_water_temp*self.lat_weighting*sea_ice_mask*land_mask).sum()
        mean_sea_water_temp_denominator = (self.lat_weighting*sea_ice_mask*land_mask).sum()
        if mean_sea_water_temp_denominator == 0:
            mean_sea_water_temp = nan
        else:
            mean_sea_water_temp = mean_sea_water_temp_numerator / mean_sea_water_temp_denominator

        is_leap_year = TimeSeriesUtils.getDaysInYear(start_dt.year) == 366

        if self.output_anomaly:
            # work out the relevant start and end time climatology indices in the range 0..364 (0..365 in leap years)
            start_tindex = start_dt.timetuple().tm_yday-1
            end_tindex = end_dt.timetuple().tm_yday-1

            if is_leap_year:
                climatology = self.leap_climatology[start_tindex:end_tindex + 1, :,:]
            else:
                climatology = self.climatology[start_tindex:end_tindex+1,:,:]
            mean_sea_water_climatology_numerator = (climatology*self.lat_weighting*sea_ice_mask*land_mask).sum()
            mean_sea_water_climatology_denominator = (self.lat_weighting*sea_ice_mask*land_mask).sum()
            if mean_sea_water_climatology_denominator == 0:
                mean_sea_water_climatology = nan
            else:
                mean_sea_water_climatology = mean_sea_water_climatology_numerator / mean_sea_water_climatology_denominator
            sst_or_anomaly = mean_sea_water_temp - mean_sea_water_climatology
        else:
            sst_or_anomaly = mean_sea_water_temp

        sea_ice_fraction_numerator = (self.lat_weighting*sea_ice_fractions*land_mask).sum()
        sea_ice_fraction_denominator = (self.lat_weighting*land_mask).sum()
        if sea_ice_fraction_denominator == 0:
            sea_ice_fraction = 0.0
        else:
            sea_ice_fraction = sea_ice_fraction_numerator / sea_ice_fraction_denominator


        # uncertainty calculation

        standard_squared_error = np.square(standard_error)
        N = (sea_ice_mask*land_mask).sum() # count used ssts
        days_duration = 1+(end_dt - start_dt).days
        mid_lat = (self.lat_max + self.lat_min) / 2
        k = max((self.lat_max - self.lat_min) / self.spatial_lambda, 1.0) \
                * max((self.lon_max - self.lon_min) / (self.spatial_lambda/cos(radians(mid_lat))), 1.0) \
                * max(days_duration / float(self.tau), 1.0)
        n = min(k,N)
        uncertainty_numerator = (standard_squared_error*self.lat_weighting*sea_ice_mask*land_mask).sum()
        uncertainty_demoninator = (self.lat_weighting*sea_ice_mask*land_mask).sum() * n
        if uncertainty_demoninator == 0:
            uncertainty = nan
        else:
            uncertainty = sqrt(uncertainty_numerator/uncertainty_demoninator)

        return (sst_or_anomaly,uncertainty,sea_ice_fraction)

    def setClimatology(self,climatology):
        """
        Assign a climatology (needed before calling aggregate if anomalies need to be calculated)
        :param climatology: a cf.FieldList containing the 365-day climatology and 0.05 degree resolution for the bounded area
        """
        self.climatology = np.nan_to_num(climatology[:,:,:,CLIMATOLOGY_FIELD_NAMES.index("sea_water_temperature")],copy=False)

        # leap climatology adds a extra day.  28th Feb is the day with index 58.
        d59 = self.climatology[58:59,:,:]
        d60 = self.climatology[59:60,:,:]
        leap_day = (d59+d60)/2
        self.leap_climatology = np.append(self.climatology[:59,:,:],np.append(leap_day,self.climatology[59:,:,:],axis=0),axis=0)

HEADER2a = "Time series of %(time_resolution)s %(var_name)s averaged across ocean surface"+ \
    " within latitudes %(south_lat)0.2f and %(north_lat)0.2f and longitudes %(west_lon)0.2f and %(east_lat)0.2f," + \
    " in kelvin and degree Celsius. "
HEADER2b = "Areas where sea-ice concentration exceeded %(fmax_pct)d %% were excluded in the average temperature calculation. "

class TimeSeriesCSVFormatter(object):

    def __init__(self,path="",time_resolution="daily",f_max=1.0,spatial_lambda=1.0,tau=3,output_anomaly=False,output_sea_ice=False,lon_min=0,lat_min=0,lon_max=5,lat_max=5,comment=""):
        """
        Construct the csv formatter using options

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
        :param comment: an extra comment string to add to the output file metadata
        """
        self.output_path = path
        self.time_resolution = time_resolution
        self.output_anomaly = output_anomaly
        self.output_sea_ice = output_sea_ice
        self.f_max = f_max
        self.spatial_lambda = spatial_lambda
        self.tau = tau

        # open a writer to handle the CSV format
        self.outfile = open(self.output_path, "w")
        self.writer = csv.writer(self.outfile)

        # define the formats for writing floating point values to CSV, using appropriate precisions
        def invalid(val):
            return isnan(val) or isinf(val)

        self.format_sst_or_anomaly = lambda val: "" if invalid(val) else "%0.5f"%val
        self.format_sea_ice_fraction = lambda val:"" if invalid(val) else "%0.5f"%val
        self.format_uncertainty = lambda val: "" if invalid(val) else "%0.5f"%val

        # get the three header lines ready and write them to the file
        header1 = "%(license)s. %(project)s, %(creator_processing_institution)s (%(publisher_url)s)"%GLOBAL_PROPERTIES
        header2 = HEADER2a%({
            "time_resolution":self.time_resolution,
            "var_name": "sea surface temperature anomaly" if self.output_anomaly else "sea surface_temperature",
            "south_lat": lat_min,
            "north_lat": lat_max,
            "west_lon": lon_min,
            "east_lat": lon_max
        })

        if self.f_max < 1.0:
            header2 += HEADER2b%({ "fmax_pct":self.f_max*100 })

        if comment:
            header2 += comment

        var_name = "mean temperature"
        if self.output_anomaly:
            var_name += " anomaly"
        column_headers = ["year","month","day",var_name+ " kelvin"]
        if not self.output_anomaly:
            # if producing absolute values as output, include a columns with kelvin centigrade units
            column_headers.append(var_name+ " deg C")
        column_headers.append(var_name+ " uncertainty")
        if self.output_sea_ice:
            column_headers += ["fraction of sea-ice-covered ocean"]

        self.writer.writerow([header1])
        self.writer.writerow([header2])
        self.writer.writerow(column_headers)

    def write(self,s_dt,mid_dt,e_dt,sst_or_anomaly,uncertainty,sea_ice_fraction):
        """
        Write an entry to the output file covering a time period
        :param start_dt: start date of the period
        :param mid_dt: mid date of the period
        :param end_dt: end date of the period
        :param sst_or_anomaly: the absolute SST or anomaly SST value
        :param uncertainty: the estimated uncertainty standard error
        :param sea_ice_fraction: the sea ice fraction
        """
        # anomalies are the difference between temperatures so do not need to be converted from kelvin to centigrade
        sst_or_anomaly_degC = [TimeSeriesUtils.kToDegC(sst_or_anomaly)] if not self.output_anomaly else []

        outrow = [mid_dt.year, mid_dt.month, mid_dt.day,
                  self.format_sst_or_anomaly(sst_or_anomaly)] + sst_or_anomaly_degC + [self.format_uncertainty(uncertainty)]
        if self.output_sea_ice:
            outrow.append(self.format_sea_ice_fraction(sea_ice_fraction))
        self.writer.writerow(outrow)

    def close(self):
        """close the formatter and flush changes to disk"""
        del self.writer
        self.writer = None
        self.outfile.close()

class TimeSeriesNetCDF4Formatter(object):
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

    def createField(self,values,name,long_name,units):
        """
        Create a cf.Field object with a time axis and metadata
        :param values: the array of values to form the field's data
        :param name: the field short name
        :param long_name: the field long name
        :param units: the field's measurement units
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

        # Set the field data from the passed values
        arr = np.array(values).reshape((-1))
        m_arr = np.ma.masked_invalid(arr)
        data = cf.Data(m_arr, units=units, dtype='i4' if isinstance(values[0],int) else "f4")
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
        year_field = self.createField(years,"calendar_year","calendar year","1")
        month_field = self.createField(months, "calendar_month", "calendar month","1")
        day_field = self.createField(days, "day_of_month", "day of month","1")
        doy_field = self.createField(doys, "day_of_year", "day of year","1")
        sst_or_anomaly_field = self.createField(ssts_or_anomalies,
                    "sst_anomaly" if self.output_anomaly else "sst",
                    "sea_water_temperature_anomaly" if self.output_anomaly else "sea_water_temperature","K")
        uncertainty_field = self.createField(uncertainties, "sst_uncertainty", "sea_water_temperature uncertainty","K")
        fields = [year_field,month_field,day_field,doy_field,sst_or_anomaly_field,uncertainty_field]
        if self.output_sea_ice:
            sea_ice_field = self.createField(sea_ice_fractions, "sea_ice_area_fraction", "sea_ice_area_fraction","1")
            fields.append(sea_ice_field)

        # write out to file
        fl = cf.FieldList(fields)
        cf.write(fl,self.output_path)
        fl.close()



class TimeSeriesExtractor(object):

    """
    This class handles the efficient extraction of data from the input (spatially resliced CCI/C3S/Climatology files)
    for the exact time and space bounded regions to be processed, and providing an iterator to lazily return data
    for each discrete time period.
    """

    def __init__(self,base_folder):
        """
        Constructor

        :param base_folder:
            the folder where the input (spatially resliced) dataset is stored.
        """
        self.base_folder = base_folder

        self.sst_field_names = SST_FIELD_NAMES
        self.climatology_field_names = CLIMATOLOGY_FIELD_NAMES


    def getOffsets(self,min_lon,min_lat,max_lon,max_lat):
        """
        Work out a set of files required to span a bounding box, the offsets into the area covered by the files,
        and the width and height of the bounding box in 0.05 degree increments

        NOTE: assumes input data has 0,05 degree resolution!

        :param min_lon:  minimum longitude of box, must be aligned on 0.05 degree boundary
        :param min_lat:  minimum latitude of box, must be aligned on 0.05 degree boundary
        :param max_lon:  maximum longitude of box, must be aligned on 0.05 degree boundary
        :param max_lat:  maximum latitude of box, must be aligned on 0.05 degree boundary

        :return: (start-offset in longitude dimension, start-offset in latitude dimension,width of longitude slice,height of latitude slice)
        """
        height = round((max_lat - min_lat) / 0.05)
        width = round((max_lon - min_lon) / 0.05)
        lat_offset = round((min_lat+90)/0.05)
        lon_offset = round((min_lon+180)/0.05)

        return (lon_offset,lat_offset,width,height)

    @staticmethod
    def createTimePeriods(time_resolution,start_date,end_date):
        """Given a time resolution and an inclusive start and end date, find the contained time periods

        :param time_resolution:  the time resolution as "daily"|"pentad"|"dekad"|"N" where N is an integer number of days
        :param start_date: the datetime of the start day (inclusive).  Time must be set to mid day.
        :param end_date: the datetime of the end day (inclusive).  Time must be set to mid day.

        start_date and end_date must be within the same year

        :return: a list of time periods in the date range expressed as (start,middle,end) tuples

        where start is a datetime set to mid-day on the first day of the period
        where end is a datetime set to mid-day on the last day of the period
        where middle is a datetime representing the middle of the period.  The time may be set to mid-day or midnight.
        """
        periods = []
        dt = start_date
        last_day_in_end_month = TimeSeriesUtils.lastDayInMonth(end_date.year,end_date.month)

        # basic date range checks
        if start_date.hour != 12 or start_date.minute != 0 or start_date.second != 0:
            raise Exception("start date time is not midday")
        if end_date.hour != 12 or end_date.minute != 0 or end_date.second != 0:
            raise Exception("end date time is not midday")
        if start_date.year != end_date.year:
            raise Exception("start date and end date must be within the same year")
        if start_date > end_date:
            raise Exception("start date cannot be later than end date")

        # work out the periods based on the desired time resolution
        if time_resolution == "daily":
            while dt <= end_date:
                periods.append((dt-datetime.timedelta(seconds=12*60*60),dt,dt+datetime.timedelta(seconds=12*60*60)))
                dt += datetime.timedelta(1)
        elif time_resolution == "monthly":
            # check start and end dates are aligned on month
            if start_date.day != 1 or end_date.day != last_day_in_end_month:
                raise Exception("internal error, start/end dates not correctly month aligned")
            while dt <= end_date:
                y = dt.year
                m = dt.month
                dim = monthrange(y,m)[1]
                periods.append((datetime.datetime(y,m,1,12,0,0), datetime.datetime(y,m,15,12,0,0), datetime.datetime(y,m,dim,12,0,0)))
                dt = dt + datetime.timedelta(dim)
        elif time_resolution == "5-day":
            if start_date.day not in [1, 6, 11, 16, 21, 26] or end_date.day not in [5, 10, 15, 20, 25,
                                                                     TimeSeriesUtils.lastDayInMonth(end_date.year,
                                                                                                        end_date.month)]:
                raise Exception("start or end date not correctly aligned for pentads")
            while dt <= end_date:
                len_days = 5
                mid_dt = dt + datetime.timedelta(2)
                if dt.day == 26:
                    # for the last pentad, include all remaining days in the month
                    len_days = monthrange(dt.year, dt.month)[1] - 25
                    if dt.month==2:
                        mid_dt = dt + datetime.timedelta(1) # use 27th as mid point of the last short pentad in february
                periods.append((dt, mid_dt, dt + datetime.timedelta(len_days-1)))
                dt = dt + datetime.timedelta(len_days)
        elif time_resolution == "10-day":
            if start_date.day not in [1,11,21] or end_date.day not in [10,20,TimeSeriesUtils.lastDayInMonth(end_date.year,end_date.month)]:
                raise Exception("start or end date not correctly aligned for dekads")
            while dt <= end_date:
                len_days = 10
                mid_dt = dt + datetime.timedelta(4)
                if dt.day == 21:
                    # for the last dekad, include all remaining days in the month
                    len_days = monthrange(dt.year, dt.month)[1] - 20
                periods.append((dt, mid_dt, dt + datetime.timedelta(len_days-1)))
                dt = dt + datetime.timedelta(len_days)
        else:
            # try to treat time_resolution as an integer number of days
            try:
                time_resolution = int(time_resolution)
            except:
                raise Exception("Invalid time resolution %s"%(time_resolution))
            dt = start_date
            while dt <= end_date:
                len_days = time_resolution
                peek_dt = dt + datetime.timedelta(len_days)
                if peek_dt > end_date:
                    # do not "overflow" the year, make the last period run to end
                    len_days = 1+int((end_date - dt).total_seconds()/(60*60*24))
                mid_dt = dt - datetime.timedelta(0.5) + datetime.timedelta(len_days/2)
                periods.append((dt, mid_dt, dt + datetime.timedelta(len_days-1)))
                dt += datetime.timedelta(len_days)
            # now check if the last period is too short
            (last_start,_,last_end) = periods[-1]

            # truncate the last period to fit into the date range, first work out the length in days
            last_dur = 1+int((end_date - last_start).total_seconds()/(60*60*24))

            if last_dur < time_resolution/2 and len(periods) > 1:
                # then merge last two time periods to yield a period that will be less than N*1.5
                (second_last_start, _, _) = periods[-2]
                periods = periods[:-1]
                last_dur = 1+int((end_date - second_last_start).total_seconds() / (60 * 60 * 24))
                mid = second_last_start-datetime.timedelta(0.5)+datetime.timedelta(last_dur/2)
                periods[-1] = (second_last_start,mid,end_date)

        return periods

    def generateYearData(self,start_date,end_date,time_resolution,min_lon,min_lat,max_lon,max_lat):
        """Generator that yields the time period within a year

        :param start_date: the datetime of the start day (inclusive).  Time must be set to mid day.
        :param end_date: the datetime of the end day (inclusive).  Time must be set to mid day.
        :param time_resolution:  the time resolution as "daily"|"pentad"|"dekad"|"N" where N is an integer number of days
        :param min_lon:  minimum longitude of box, must be aligned on 0.05 degree boundary
        :param min_lat:  minimum latitude of box, must be aligned on 0.05 degree boundary
        :param max_lon:  maximum longitude of box, must be aligned on 0.05 degree boundary
        :param max_lat:  maximum latitude of box, must be aligned on 0.05 degree boundary

        The generator yields ((start_dt,mid_dt,end_dt),cf.FieldList) tuples
        """

        # first work out how to map the spatial box to the input data

        input_path = os.path.join(self.base_folder,"%d.zarr"%(start_date.year))
        z = zarr.open(input_path,mode='r')

        (lon_offset,lat_offset,width,height) = self.getOffsets(min_lon,min_lat,max_lon,max_lat)

        time_periods = TimeSeriesExtractor.createTimePeriods(time_resolution,start_date,end_date)

        for(period_start_dt,period_mid_dt,period_end_dt) in time_periods:
            doy_start = period_start_dt.timetuple().tm_yday - 1
            doy_end = period_end_dt.timetuple().tm_yday - 1
            subf = z[doy_start:doy_end+1,lat_offset:lat_offset+height,lon_offset:lon_offset+width,:]
            yield ((period_start_dt,period_mid_dt,period_end_dt),subf)

    def generateData(self, start_dt, end_dt, time_resolution, min_lon, min_lat, max_lon, max_lat):
        """Generator that lazily yields the time period data for a given time and space range

        :param start_date: the datetime of the start day (inclusive).  Time must be set to mid day.
        :param end_date: the datetime of the end day (inclusive).  Time must be set to mid day.
        :param time_resolution:  the time resolution as "daily"|"pentad"|"dekad"|"N" where N is an integer number of days
        :param min_lon:  minimum longitude of box, must be aligned on 0.05 degree boundary
        :param min_lat:  minimum latitude of box, must be aligned on 0.05 degree boundary
        :param max_lon:  maximum longitude of box, must be aligned on 0.05 degree boundary
        :param max_lat:  maximum latitude of box, must be aligned on 0.05 degree boundary

        The generator yields ((start_dt,mid_dt,end_dt),cf.FieldList) tuples
        """
        year = start_dt.year
        while year <= end_dt.year:
            # go through each year in turn...
            slice_end_dt = datetime.datetime(year,12,31,12,0,0) if year < end_dt.year else end_dt
            slice_start_dt = datetime.datetime(year,1,1,12,0,0) if year > start_dt.year else start_dt
            # yield from that year's generator until exhausted
            yield from self.generateYearData(slice_start_dt,slice_end_dt,time_resolution,min_lon,min_lat,max_lon,max_lat)
            # move to the next year
            year += 1

    def extractClimatology(self, min_lon, min_lat, max_lon, max_lat):
        """Extract the climatology for a given bounding box

        :param min_lon:  minimum longitude of box, must be aligned on 0.05 degree boundary
        :param min_lat:  minimum latitude of box, must be aligned on 0.05 degree boundary
        :param max_lon:  maximum longitude of box, must be aligned on 0.05 degree boundary
        :param max_lat:  maximum latitude of box, must be aligned on 0.05 degree boundary

        :return: a cf.FieldList containing the 365-day climatology covering the bounded area
        """
        input_path = os.path.join(self.base_folder, "climatology.zarr")
        z = zarr.open(input_path, mode='r')

        (lon_offset, lat_offset, width, height) = self.getOffsets(min_lon, min_lat, max_lon, max_lat)

        return z[:, lat_offset:lat_offset + height, lon_offset:lon_offset + width, :]


def makeTimeSeriesSSTs(lon_min,lon_max,lat_min,lat_max,time_resolution,
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
    extractor = TimeSeriesExtractor(input_path)

    # create an aggregator to aggregate each period in the extracted data
    aggregator = TimeSeriesAggregator(lon_min=lon_min, lat_min=lat_min, lon_max=lon_max,
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

    if output_csv:
        formatter = TimeSeriesCSVFormatter(**formatter_args)
    else:
        formatter = TimeSeriesNetCDF4Formatter(**formatter_args)

    # if outputting anomalies, we'll need the climatology.  Extract it and pass it to the aggregator
    if anomalies:
        climatology = extractor.extractClimatology(min_lon=lon_min, min_lat=lat_min, max_lon=lon_max, max_lat=lat_max)
        aggregator.setClimatology(climatology)

    p = Progress("makeTimeSeriesSSTs")
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

        frac = (mid_dt.timestamp() - start_date.timestamp()) / period_duration
        p.report("Processing",frac)

    p.clear()

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

    makeTimeSeriesSSTs(lon_min=args.lon_min, lat_min=args.lat_min, lon_max=args.lon_max, lat_max=args.lat_max,
                       time_resolution=args.time_resolution, start_date=start_dt, end_date=end_dt,
                       out_path=args.out_path, input_path=args.in_path,
                       f_max=args.f_max, spatial_lambda=args.spatial_lambda, tau=args.tau,
                       anomalies=args.anomalies, no_sea_ice_fraction=args.no_sea_ice_fraction,
                       comment=args.comment)


if __name__ == '__main__':
    parser = createParser()
    args = parser.parse_args()
    dispatch(args)

