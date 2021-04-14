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
This module handles the aggregation of SST, uncertainty and sea ice fractions for a bounded area and time period

Input data is provided as a 3 dimensional array organised by (time,lat,lon)

Aggregation for regions aplies aggregation in the time dimension only
Aggregation for timeseries applies aggregation across all dimensions
"""
from math import cos, radians, nan, sqrt
import numpy as np
from .utils import SST_FIELD_NAMES, CLIMATOLOGY_FIELD_NAMES, TimeSeriesUtils

class Aggregator(object):
    """
    Aggregate SST data from a 3 dimensional box (bounded by time, latitude and longitude) to yield sst, anomaly, uncertainty
    and sea ice fraction values according to the requirements

    In "region" mode the yielded values are 2D numpy arrays, in timeseries mode they are scalar values

    Note:
        after creating an Aggregator instance, call setClimatology if
    """

    def __init__(self, mode, lon_min, lat_min, lon_max, lat_max, f_max=1.0, output_anomaly=False, output_sea_ice=False,
                 spatial_lambda=1.0, tau=3):
        """
        Construct the aggregator

        :param mode: "timeseries" or "region": whether to output data for time series (scalar values) or region (arrays)
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
        if mode not in ["timeseries", "region"]:
            raise Exception("invalid value specified for mode (%s) should be either 'timeseries' or 'region'" % (mode))
        self.mode = mode
        self.lon_min = lon_min
        self.lon_max = lon_max
        self.lat_min = lat_min
        self.lat_max = lat_max
        self.lat_weighting = None  # TimeSeriesUtils.createLatitudeWeighting(lat_min, lat_max)
        self.f_max = f_max
        self.output_anomaly = output_anomaly
        self.output_sea_ice = output_sea_ice
        self.spatial_lambda = spatial_lambda
        self.tau = tau
        self.climatology = None

    def aggregate(self, start_dt, end_dt, data):
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

        sea_ice_fractions = np.nan_to_num(data[:, :, :, sea_ice_fraction_index], copy=False)
        sea_water_temp = np.nan_to_num(data[:, :, :, sea_water_temp_index], copy=False)
        standard_error = np.nan_to_num(data[:, :, :, standard_error_index], copy=False)
        land_mask = np.where(sea_water_temp > 0.0, 1.0, 0.0)

        if self.lat_weighting is None:
            self.lat_weighting = TimeSeriesUtils.createLatitudeWeighting(lat_min=self.lat_min, lat_max=self.lat_max)

        sea_ice_mask = np.where(sea_ice_fractions > self.f_max, 0.0, 1.0)

        axis_params = {"axis": 0} if self.mode == "region" else {}

        mean_sea_water_temp_numerator = (sea_water_temp * self.lat_weighting * sea_ice_mask * land_mask).sum(
            **axis_params)
        mean_sea_water_temp_denominator = (self.lat_weighting * sea_ice_mask * land_mask).sum(**axis_params)
        if self.mode == "timeseries" and mean_sea_water_temp_denominator == 0:
            mean_sea_water_temp = nan
        else:
            mean_sea_water_temp = mean_sea_water_temp_numerator / mean_sea_water_temp_denominator

        is_leap_year = TimeSeriesUtils.getDaysInYear(start_dt.year) == 366

        if self.output_anomaly:
            # work out the relevant start and end time climatology indices in the range 0..364 (0..365 in leap years)
            start_tindex = start_dt.timetuple().tm_yday - 1
            end_tindex = end_dt.timetuple().tm_yday - 1

            if is_leap_year:
                climatology = self.leap_climatology[start_tindex:end_tindex + 1, :, :]
            else:
                climatology = self.climatology[start_tindex:end_tindex + 1, :, :]
            mean_sea_water_climatology_numerator = (climatology * self.lat_weighting * sea_ice_mask * land_mask).sum(
                **axis_params)
            mean_sea_water_climatology_denominator = (self.lat_weighting * sea_ice_mask * land_mask).sum(**axis_params)
            if self.mode == "timeseries" and mean_sea_water_climatology_denominator == 0:
                mean_sea_water_climatology = nan
            else:
                mean_sea_water_climatology = mean_sea_water_climatology_numerator / mean_sea_water_climatology_denominator

            sst_or_anomaly = mean_sea_water_temp - mean_sea_water_climatology
        else:
            sst_or_anomaly = mean_sea_water_temp

        sea_ice_fraction_numerator = (self.lat_weighting * sea_ice_fractions * land_mask).sum(**axis_params)
        sea_ice_fraction_denominator = (self.lat_weighting * land_mask).sum(**axis_params)
        sea_ice_fraction = sea_ice_fraction_numerator / sea_ice_fraction_denominator

        # uncertainty calculation

        standard_squared_error = np.square(standard_error)
        N = (sea_ice_mask * land_mask).sum(**axis_params)  # count used ssts for each pixel
        days_duration = 1 + (end_dt - start_dt).days

        if self.mode == "region":
            k = max(days_duration / float(self.tau), 1.0)
            n = np.minimum(k, N)
            uncertainty_numerator = (standard_squared_error * self.lat_weighting * sea_ice_mask * land_mask).sum(
                **axis_params)
            uncertainty_demoninator = (self.lat_weighting * sea_ice_mask * land_mask).sum(**axis_params) * n
            uncertainty = np.sqrt(uncertainty_numerator / uncertainty_demoninator)
        else:
            mid_lat = (self.lat_max + self.lat_min) / 2
            k = max((self.lat_max - self.lat_min) / self.spatial_lambda, 1.0) \
                * max((self.lon_max - self.lon_min) / (self.spatial_lambda / cos(radians(mid_lat))), 1.0) \
                * max(days_duration / float(self.tau), 1.0)
            n = min(k, N)
            uncertainty_numerator = (standard_squared_error * self.lat_weighting * sea_ice_mask * land_mask).sum()
            uncertainty_demoninator = (self.lat_weighting * sea_ice_mask * land_mask).sum() * n
            if uncertainty_demoninator == 0:
                uncertainty = nan
            else:
                uncertainty = sqrt(uncertainty_numerator / uncertainty_demoninator)

        return (sst_or_anomaly, uncertainty, sea_ice_fraction)

    def setClimatology(self, climatology):
        """
        Assign a climatology (needed before calling aggregate if anomalies need to be calculated)
        :param climatology: a cf.FieldList containing the 365-day climatology and 0.05 degree resolution for the bounded area
        """
        self.climatology = np.nan_to_num(climatology[:, :, :, CLIMATOLOGY_FIELD_NAMES.index("sea_water_temperature")],
                                         copy=False)

        # leap climatology adds a extra day.  28th Feb is the day with index 58.
        d59 = self.climatology[58:59, :, :]
        d60 = self.climatology[59:60, :, :]
        leap_day = (d59 + d60) / 2
        self.leap_climatology = np.append(self.climatology[:59, :, :],
                                          np.append(leap_day, self.climatology[59:, :, :], axis=0), axis=0)
