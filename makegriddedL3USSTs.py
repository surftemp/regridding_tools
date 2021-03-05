#!/usr/bin/env python3

import argparse
import os
import zipfile
from calendar import monthrange
from datetime import date, datetime
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

        # Read in the first SST L3U field
        flat_list = [item for sublist in self.filename_groups for subsublist in sublist for item in subsublist]
        if not flat_list:
            raise ValueError('No input data found.')
        fl = cf.read(flat_list[0])
        sst = fl.select_by_property(standard_name='sea_water_temperature')[0]

        # Make the weights
        self.weights = sst.weights('area', scale=1.0)

        # Save the original longitude and latitude as those of the climatology need to be replaced with these as
        # they differ in metadata and by rounding errors.
        self.orig_lon = sst.dim('X').copy()
        self.orig_lat = sst.dim('Y').copy()

        # Resample the longitude and latitude
        self.resample_lonlat(sst)

        # Close the file
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
        output_paths = []
        # The main code to perform the regridding
        for filename_group, climatology_filename_group in zip(self.filename_groups, self.climatology_file_name_groups):
            resampled_sst_climatology_data = None
            sst_climatology_denominator = None
            resampled_sst_data = None
            sst_denominator = None
            n_data = None

            # Create a flat list of the input files and check it is not empty
            flat_list = [item for sublist in filename_group for item in sublist]
            if not flat_list:
                continue

            # Read in and regrid the SST
            for filenames, climatology_filename in zip(filename_group, climatology_filename_group):
                # Read in and regrid the climatology
                if type(climatology_filename) is str:
                    gl = cf.read(climatology_filename, aggregate=False)
                    sst_climatology = gl.select_by_property(standard_name='sea_water_temperature')[0]
                else:
                    # The climatology is a tuple of two fields on 29 February, which must be averaged
                    gl = cf.read(climatology_filename[0], aggregate=False)
                    sst_climatology = gl.select_by_property(standard_name='sea_water_temperature')[0]
                    hl = cf.read(climatology_filename[1], aggregate=False)
                    sst_climatology += hl.select_by_property(standard_name='sea_water_temperature')[0]
                    sst_climatology /= 2
                # Replace the latitude and longitude dimension coordinates in the climatology with the ones from the L3U
                # data as they are incompatible because of metadata and rounding differences.
                sst_climatology.replace_construct('X', self.orig_lon)
                sst_climatology.replace_construct('Y', self.orig_lat)
                # Resample the climatology to the lower resolution and aggregate it
                data = self.spatially_resample_data(sst_climatology, weights=True)
                if resampled_sst_climatology_data is None:
                    resampled_sst_climatology_data = data
                else:
                    resampled_sst_climatology_data = regridding_utilities.add_data(resampled_sst_climatology_data, data)
                gl.close()
                if type(climatology_filename) is not str:
                    hl.close()
                # Calculate the denominator for averaging
                data = self.spatially_resample_data(sst_climatology.where(True, 1.0), weights=True)
                if sst_climatology_denominator is None:
                    sst_climatology_denominator = data
                else:
                    sst_climatology_denominator = regridding_utilities.add_data(sst_climatology_denominator, data)

                for filename in filenames:
                    # Read in data
                    fl = cf.read(filename, aggregate=False)

                    # Select the SST and create longitude and latitude bounds if necessary
                    sst = fl.select_by_property(standard_name='sea_water_temperature')[0]
                    regridding_utilities.create_lonlat_bounds(sst)

                    # Select the quality level
                    qlevel = fl.select_by_ncvar('quality_level')[0]

                    # Calculate the regridded SST anomaly
                    sst -= sst_climatology
                    data = self.spatially_resample_data(sst.where(qlevel < self.min_qlevel, cf.masked), weights=True)
                    if resampled_sst_data is None:
                        resampled_sst_data = data
                    else:
                        resampled_sst_data = regridding_utilities.add_data(resampled_sst_data, data)

                    # Calculate the denominator for averaging
                    data = self.spatially_resample_data(sst.where(qlevel < self.min_qlevel, cf.masked, 1.0), weights=True)
                    if sst_denominator is None:
                        sst_denominator = data
                    else:
                        sst_denominator = regridding_utilities.add_data(sst_denominator, data)

                    # Calculate the number of observations used in each target cell
                    f_tmp = sst.where(qlevel < self.min_qlevel, 0, 1)
                    f_tmp.data.filled(0, inplace=True)
                    f_tmp.dtype = np.int32
                    data = self.spatially_resample_data(f_tmp)
                    if n_data is None:
                        n_data = data
                    else:
                        n_data = regridding_utilities.add_data(n_data, data)

                    fl.close()

            # Finalise the climatology regridding calculation
            resampled_sst_climatology_data /= sst_climatology_denominator

            # Average the summed data
            resampled_sst_data /= sst_denominator
            if not self.anomalies:
                resampled_sst_data += resampled_sst_climatology_data.array

            # Create the resampled sst field
            resampled_sst = self.update_field(flat_list, 'sea_water_temperature', resampled_sst_data)
            if self.anomalies:
                resampled_sst.nc_set_variable('sst_anomaly')
                resampled_sst.standard_name = 'sea_water_temperature_anomaly'
                resampled_sst.long_name = 'regridded L3U sea surface temperature anomaly'
            else:
                resampled_sst.nc_set_variable('sst')
                resampled_sst.standard_name = 'sea_water_temperature'
                resampled_sst.long_name = 'regridded L3U sea surface temperature'
            resampled_sst.override_units('kelvin', inplace=True)
            resampled_sst.comment = 'These data were produced by the University of Reading as part of the ESA ' + \
                                    'CCI project.'

            # Create the field n with the number of data points contributing to each cell
            n = self.update_field(flat_list, 'sea_water_temperature', n_data)
            n.nc_set_variable('n')
            n.standard_name = 'number_of_observations'
            n.long_name = 'number of observations'
            n.override_units('1', inplace=True)
            n.comment = 'Number of L3U cells contributing to the average.'
            n.del_property('depth')
            n.del_property('source')

            # Get the date of the resampled SST
            dt = resampled_sst.dimension_coordinate('T').datetime_array[0]

            # Create a field list with all the fields in it
            fl = cf.FieldList()
            fl.append(resampled_sst)
            fl.append(n)
            fl.append(regridding_utilities.create_time_field(n, 'calendar_year', 'calendar year', dt.year))
            fl.append(regridding_utilities.create_time_field(n, 'calendar_month', 'calendar month', dt.month))
            fl.append(regridding_utilities.create_time_field(n, 'day_of_month', 'day of month', dt.day))
            fl.append(regridding_utilities.create_time_field(n, 'day_of_year', 'day of year', dt.dayofyr))

            # Write the data
            output_path = os.path.join(self.out_path, dt.strftime('%Y%m%d') + '_regridded_sst.nc')
            cf.write(fl, output_path, datatype={np.dtype('float64'): np.dtype('float32'),
                                                np.dtype('int64'): np.dtype('int32')},
                     compress=1, least_significant_digit=3)
            output_paths.append(output_path)

        # if specified, zip up all output files and remove the originals
        if self.zip_name:
            zip_path = os.path.join(self.out_path, self.zip_name)
            with zipfile.ZipFile(zip_path, "w") as z:
                for output_path in output_paths:
                    z.write(output_path, os.path.split(output_path)[1])
                    os.unlink(output_path)


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
