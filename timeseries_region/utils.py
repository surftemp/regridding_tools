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
This module contains common utilities used by the timeseries and region regridding code
"""

import datetime
from math import cos, radians, nan, sqrt
from calendar import monthrange
import numpy as np

SST_FIELD_NAMES = ['sea_water_temperature', 'sea_water_temperature standard_error', 'sea_ice_area_fraction']
CLIMATOLOGY_FIELD_NAMES = ['sea_water_temperature', 'sea_ice_area_fraction']

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


def createTimePeriods(time_resolution, start_date, end_date):
    """Given a time resolution and an inclusive start and end date, find the contained time periods

    :param time_resolution:  the time resolution as "pentad"|"dekad"|"N" where N is an integer number of days
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
    last_day_in_end_month = TimeSeriesUtils.lastDayInMonth(end_date.year, end_date.month)

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
    if time_resolution == "monthly":
        # check start and end dates are aligned on month
        if start_date.day != 1 or end_date.day != last_day_in_end_month:
            raise Exception("internal error, start/end dates not correctly month aligned")
        while dt <= end_date:
            y = dt.year
            m = dt.month
            dim = monthrange(y, m)[1]
            periods.append((datetime.datetime(y, m, 1, 12, 0, 0), datetime.datetime(y, m, 15, 12, 0, 0),
                            datetime.datetime(y, m, dim, 12, 0, 0)))
            dt = dt + datetime.timedelta(dim)
    elif time_resolution == "5-day":
        if start_date.day not in [1, 6, 11, 16, 21, 26] or end_date.day not in [5, 10, 15, 20, 25,
                                                                                TimeSeriesUtils.lastDayInMonth(
                                                                                    end_date.year,
                                                                                    end_date.month)]:
            raise Exception("start or end date not correctly aligned for pentads")
        while dt <= end_date:
            len_days = 5
            mid_dt = dt + datetime.timedelta(2)
            if dt.day == 26:
                # for the last pentad, include all remaining days in the month
                len_days = monthrange(dt.year, dt.month)[1] - 25
                if dt.month == 2:
                    mid_dt = dt + datetime.timedelta(1)  # use 27th as mid point of the last short pentad in february
            periods.append((dt, mid_dt, dt + datetime.timedelta(len_days - 1)))
            dt = dt + datetime.timedelta(len_days)
    elif time_resolution == "10-day":
        if start_date.day not in [1, 11, 21] or end_date.day not in [10, 20,
                                                                     TimeSeriesUtils.lastDayInMonth(end_date.year,
                                                                                                    end_date.month)]:
            raise Exception("start or end date not correctly aligned for dekads")
        while dt <= end_date:
            len_days = 10
            mid_dt = dt + datetime.timedelta(4)
            if dt.day == 21:
                # for the last dekad, include all remaining days in the month
                len_days = monthrange(dt.year, dt.month)[1] - 20
            periods.append((dt, mid_dt, dt + datetime.timedelta(len_days - 1)))
            dt = dt + datetime.timedelta(len_days)
    else:
        # try to treat time_resolution as an integer number of days
        try:
            time_resolution = int(time_resolution)
        except:
            raise Exception("Invalid time resolution %s" % (time_resolution))
        dt = start_date
        while dt <= end_date:
            len_days = time_resolution
            peek_dt = dt + datetime.timedelta(len_days)
            if peek_dt > end_date:
                # do not "overflow" the year, make the last period run to end
                len_days = 1 + int((end_date - dt).total_seconds() / (60 * 60 * 24))
            mid_dt = dt - datetime.timedelta(0.5) + datetime.timedelta(len_days / 2)
            periods.append((dt, mid_dt, dt + datetime.timedelta(len_days - 1)))
            dt += datetime.timedelta(len_days)
        # now check if the last period is too short
        (last_start, _, last_end) = periods[-1]

        # truncate the last period to fit into the date range, first work out the length in days
        last_dur = 1 + int((end_date - last_start).total_seconds() / (60 * 60 * 24))

        if last_dur < time_resolution / 2 and len(periods) > 1:
            # then merge last two time periods to yield a period that will be less than N*1.5
            (second_last_start, _, _) = periods[-2]
            periods = periods[:-1]
            last_dur = 1 + int((end_date - second_last_start).total_seconds() / (60 * 60 * 24))
            mid = second_last_start - datetime.timedelta(0.5) + datetime.timedelta(last_dur / 2)
            periods[-1] = (second_last_start, mid, end_date)

    return periods