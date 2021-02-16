#!/usr/bin/env python3

import argparse
import os
import zipfile
from calendar import monthrange, isleap
from datetime import date, datetime, timedelta
from uuid import uuid4

import cf
import numpy as np

import regridding_constants
import regridding_utilities

_default_year = regridding_constants.cci_start_year
_default_start_month = 1
_default_end_month = 12
_default_anomalies = False
_default_min_qlevel = 4
_default_out_path = '/Users/charles/Data/l3u_regrid_test'
_default_zip_name = ''


class L3USSTRegridder(regridding_utilities.Regridder):
    """
    Class for regridding L3U SSTs.

    Should be instantiated first, optionally with paths to the input data and then a call made to the resulting object
    to perform the regridding at a particular resolution with other options available.
    """
    def __init__(self, sst_l3u_path=regridding_constants.default_sst_l3u_path,
                 sst_cci_climatology_path=regridding_constants.default_sst_cci_climatology_path, ):
        """
        Initiate an L3U SST regridder optionally specfying the input file paths.

        :param sst_l3u_path:
            Path to the SST L3U input data.

        :param sst_cci_climatology_path:
            Path to the SST CCI Climatology input data.
        """
        self.input_type = regridding_utilities.InputType.SST_L3U

        self.sst_l3u_path = sst_l3u_path
        self.sst_cci_climatology_path = sst_cci_climatology_path

        self.base_resolution = 0.05
        self.min_resolution = 0.25
        self.max_resolution = 5.0
        self.lon_range = 360.0
        self.lat_range = 180.0

        self.allowed_lon_resolutions = regridding_utilities.allowed_lonlat_resolutions(self.base_resolution,
                                                                                       self.max_resolution,
                                                                                       self.lon_range,
                                                                                       self.min_resolution)
        self.allowed_lat_resolutions = regridding_utilities.allowed_lonlat_resolutions(self.base_resolution,
                                                                                       self.max_resolution,
                                                                                       self.lat_range,
                                                                                       self.min_resolution)

    def __call__(self, lon_resolution, lat_resolution, time_resolution, sensor, year=_default_year,
                 start_month=_default_start_month, end_month=_default_end_month, anomalies=_default_anomalies,
                 min_qlevel=_default_min_qlevel, out_path=_default_out_path, zip_name=_default_zip_name):
        """
        Regrid the L3U SSTs to a target resolution.

        :param lon_resolution:
            The target longitude resolution. This must be a multiple of 0.05 degrees and a factor of 360 degrees. It
            must also not be greater than 5 degrees.

        :param lat_resolution:
            The target latitude resolution. This must be a multiple of 0.05 degrees and a factor of 180 degrees. It must
            also not be greater than 5 degrees.

        :param time_resolution:
            The target time resolution. This can be 'annual', 'monthly', '10-day' for dekads.

        :param sensor:
            The name of the sensor for which to perform the regridding.

        :param year:
            The year of the data to be regridded.

        :param start_month:
            The first month of the data to be regridded.

        :param end_month:
            The final month of the data to be regridded (the whole of this month is included).

        :param anomalies:
            Whether to output anomalies instead of absolute SSTs. Default set by _default_anomalies.

        :param min_qlevel:
            The minimum quality level of the L3U cell SST for inclusion. Default set by _default_min_qlevel. If 4 then
            quality levels 4 and 5 are included. Minimum is 3 such that quality levels 3, 4 and 5 anre included.
            Maximum is 5 such that only quality level 5 is included.

        :param out_path:
            The path in which to write the output file of regridded data.

        :param zip_name:
            The name of a zip file to combine all output files into.

        :return:
        """
        # Check the longitude and latitude resolutions are allowed and get the corresponding integer steps
        if lon_resolution not in self.allowed_lon_resolutions:
            raise ValueError('lon_resolution must be a multiple of ' + str(self.base_resolution) + ', a factor of ' +
                             str(self.lon_range) + ' and not greater than ' + str(self.max_resolution) + ' degrees.')
        if lat_resolution not in self.allowed_lat_resolutions:
            raise ValueError('lat_resolution must be a multiple of ' + str(self.base_resolution) + ', a factor of ' +
                             str(self.lat_range) + ' and not greater than ' + str(self.max_resolution) + ' degrees.')
        self.lon_step = self.allowed_lon_resolutions[lon_resolution]
        self.lat_step = self.allowed_lat_resolutions[lat_resolution]
        self.lon_resolution = lon_resolution
        self.lat_resolution = lat_resolution

        # Set the dates and time resolution
        self.year = year
        self.start_month = start_month
        self.start_day = 1
        self.end_month = end_month
        self.end_day = monthrange(year, end_month)[1]
        self.time_resolution = time_resolution

        # Set the sensor
        self.sensor = sensor

        # Set flags for whether to calculate anomalies
        self.anomalies = anomalies

        # Set the minimum quality level
        assert type(min_qlevel) is int, 'Minimum quality level must be an integer.'
        assert 3 <= min_qlevel <= 5, 'Minimum quality level out of range.'
        self.min_qlevel = min_qlevel

        # Check the dates and time resolution
        self.check_dates()

        # Generate a list of lists of filenames
        self.filename_groups = self.create_filename_groups()

        # Generate a list of lists of filenames for the climatology
        self.climatology_file_name_groups = self.create_filename_groups(climatology=True)

        # Set the output paths
        self.out_path = out_path
        self.zip_name = zip_name

        # Set the date created and the uuid
        self.date_created = datetime.now().strftime('%Y-%m-%dT%H:%M:%SZ')
        self.uuid = str(uuid4())

        # Read in the first SST field
        flat_list = [item for sublist in self.filename_groups for item in sublist]
        is_data = False
        for item in flat_list:
            if item:
                fl = cf.read(item[0])
                is_data = True
                break
        if not is_data:
            raise ValueError('No data found.')
        sst = fl.select_by_property(standard_name='sea_surface_skin_temperature')[0]

        # Create longitude and latitude bounds if necessary
        regridding_utilities.create_lonlat_bounds(sst)

        # Make the weights
        self.weights = sst.weights('area', scale=1.0)

        # Close the files
        fl.close()

        # Do the regridding
        self.make_gridded_l3u_ssts()

    def check_dates(self):
        """
        Check the dates and time resolution are allowed.
        """
        if self.start_month < 1 or self.start_month > 12:
            raise ValueError('Start month out of range.')
        if self.end_month < 1 or self.end_month > 12:
            raise ValueError('End month out of range.')
        if date(self.year, self.start_month, self.start_day) > date(self.year, self.end_month, self.end_day):
            raise ValueError('Start date cannot be greater than end date.')
        if self.time_resolution == 'annual':  # Annual resampling
            if self.start_month != 1 or self.end_month != 12:
                raise ValueError('A whole number of years is required for annual regridding.')
        elif self.time_resolution == 'monthly':  # Monthly resampling
            pass
        elif self.time_resolution == '10-day':  # Nominal 10-day resampling
            pass
        else:  # The time resolution is not recognised
            raise ValueError('Time resolution not recognised.')

    def make_gridded_l3u_ssts(self):
        """
        Calculate the regridded L3U SSTs.
        """
        # TODO
        pass


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Regrid L3U SST data for a particular year.')
    parser.add_argument('lon_resolution', type=float, help='The target longitude resolution. This must be a multiple ' +
                                                           'of 0.05 degrees and a factor of 360 degrees. It must ' +
                                                           'also not be greater than 5 degrees.')
    parser.add_argument('lat_resolution', type=float, help='The target latitude resolution. This must be a multiple ' +
                                                           'of 0.05 degrees and a factor of 180 degrees. It must ' +
                                                           'also not be greater than 5 degrees.')
    parser.add_argument('time_resolution', help="The target time resolution. This can be 'annual', 'monthly', " +
                                                "'10-day' for dekads.")
    parser.add_argument('sensor', help='The sensor for which to perform regridding.')
    parser.add_argument('--year', type=int, default=_default_year,
                        help='The year of the data to be regridded. Default is ' + str(_default_year) + '.')
    parser.add_argument('--start_month', type=int, default=_default_start_month,
                        help='The first month of the data to be regridded. Default is ' + str(_default_start_month)
                             + '.')
    parser.add_argument('--end_month', type=int, default=_default_end_month,
                        help='The final month of the data to be regridded. Default is ' + str(_default_end_month) + '.')
    parser.add_argument('--anomalies', action='store_true', default=False,
                        help='Output anomalies instead of absolute SSTs.')
    parser.add_argument('--min_qlevel', default=_default_min_qlevel,
                        help='The minimum quality level of the L3U cell SST for inclusion. Default is '
                             + str(_default_min_qlevel) + '. If 4 then quality levels 4 and 5 are included. Minimum is '
                             + '3 such that quality levels 3, 4 and 5 anre included. Maximum is 5 such that only '
                             + 'quality level 5 is included.')
    parser.add_argument('--sst_l3u_path', default=regridding_constants.default_sst_l3u_path,
                        help='Path to the SST Level 3 Uncollated input data.')
    parser.add_argument('--sst_cci_climatology_path', default=regridding_constants.default_sst_cci_climatology_path,
                        help='Path to the SST CCI Climatology input data.')
    parser.add_argument('--out_path', default=_default_out_path,
                        help='The path in which to write the output file of regridded data.')
    parser.add_argument('--zip_name', default=_default_zip_name,
                        help='Combine all output files into a single zip file with this name.')

    args = parser.parse_args()

    l3u_sst_regridder = L3USSTRegridder(sst_l3u_path=args.sst_l3u_path,
                                        sst_cci_climatology_path=args.sst_cci_climatology_path)

    l3u_sst_regridder(args.lon_resolution, args.lat_resolution, args.time_resolution, args.sensor, year=args.year,
                      start_month=args.start_month, end_month=args.end_month, anomalies=args.anomalies,
                      min_qlevel=args.min_qlevel, out_path=args.out_path, zip_name=args.zip_name)
