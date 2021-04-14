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
This module handles the extraction of SST, uncertainty and sea-ice fraction data from L4 SST data saved in zarr format
(see the "timeseries_region.reslicing" module for the code that creates the zarr format data)
"""

import datetime
import zarr
import os.path

from .utils import SST_FIELD_NAMES, CLIMATOLOGY_FIELD_NAMES, createTimePeriods

class Extractor(object):

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

        time_periods = createTimePeriods(time_resolution,start_date,end_date)

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

