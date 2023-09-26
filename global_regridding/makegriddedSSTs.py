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

_default_year = 2000
_default_start_month = 1
_default_start_day = 1
_default_end_month = 12
_default_end_day = 31
_default_anomalies = False
_default_no_sea_ice_fraction = False
_default_f_max = 1.0
_default_tau = 7
_default_spatial_lambda = 3.0
_default_out_path = '/Users/charles/Data/regrid_test'
_default_zip_name = ''


class SSTRegridder(regridding_utilities.Regridder):
    """
    Class for regridding L4 SSTs.

    Should be instantiated first, optionally with paths to the input data and then a call made to the resulting object
    to perform the regridding at a particular resolution with other options available.
    """
    def __init__(self, sst_analysis_l4_path=regridding_constants.default_sst_analysis_l4_path,
                 sst_climatology_path=regridding_constants.default_sst_climatology_path):
        """
        Initiate an SST regridder optionally specfying the input file paths.

        :param sst_cci_climatology_path:
            Path to the SST CCI Climatology input data.

        :param sst_analysis_l4_path:
            Path to the SST Analysis Level 4 input data.
        """
        self.input_type = regridding_utilities.InputType.SST_L4

        self.sst_analysis_l4_path = sst_analysis_l4_path
        self.sst_climatology_path = sst_climatology_path

        self.base_resolution = 0.05
        self.max_resolution = 5.0
        self.lon_range = 360.0
        self.lat_range = 180.0

        self.allowed_lon_resolutions = regridding_utilities.allowed_lonlat_resolutions(self.base_resolution,
                                                                                       self.max_resolution,
                                                                                       self.lon_range)
        self.allowed_lat_resolutions = regridding_utilities.allowed_lonlat_resolutions(self.base_resolution,
                                                                                       self.max_resolution,
                                                                                       self.lat_range)

    def __call__(self, lon_resolution, lat_resolution, time_resolution, year=_default_year,
                 start_month=_default_start_month, start_day=_default_start_day,
                 end_month=_default_end_month, end_day=_default_end_day, anomalies=_default_anomalies,
                 no_sea_ice_fraction=_default_no_sea_ice_fraction, f_max=_default_f_max, tau=_default_tau,
                 spatial_lambda=_default_spatial_lambda, out_path=_default_out_path, zip_name=_default_zip_name):
        """
        Regrid the L4 SSTs to a target resolution.

        :param lon_resolution:
            The target longitude resolution. This must be a multiple of 0.05 degrees and a factor of 360 degrees. It
            must also not be greater than 5 degrees.

        :param lat_resolution:
            The target latitude resolution. This must be a multiple of 0.05 degrees and a factor of 180 degrees. It must
            also not be greater than 5 degrees.

        :param time_resolution:
            The target time resolution. This can be 'annual', 'monthly', '10-day' for dekads, '5-day' for pentads or an
            integer for regular N day regridding aligned with the start of the year.

        :param year:
            The year of the data to be regridded.

        :param start_month:
            The first month of the data to be regridded.

        :param start_day:
            The first day of the month of the data to be regridded.

        :param end_month:
            The final month of the data to be regridded.

        :param end_day:
            The final day of the month of the data to be regridded.

        :param anomalies:
            Whether to output anomalies instead of absolute SSTs. Default set by _default_anomalies.

        :param no_sea_ice_fraction:
            Whether to output the sea ice fraction or not. Default set by _default_no_sea_ice_fraction.

        :param f_max:
            The fraction of sea ice above which a cell is ignored in calculating the regridded SST. Default set by
            _default_f_max. When calculating anomalies, the climatology is always calculated with an f_max of 1.0
            regardless of this value.

        :param tau:
            Timescale within which errors are assumed to be fully correlated in days. Default set by _default_tau.

        :param spatial_lambda:
            Spatial scale within which errors are assumed to be fully correlated in degrees. Default set by
            _default_spatial_lambda.

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
        self.start_day = start_day
        self.end_month = end_month
        self.end_day = end_day
        self.time_resolution = time_resolution

        # Set flags for whether to calculate anomalies or the sea ice fraction
        self.anomalies = anomalies
        self.sea_ice_fraction = not no_sea_ice_fraction

        # Check the dates and time resolution
        self.check_dates()

        # Generate a list of lists of filenames
        self.filename_groups = self.create_filename_groups()

        # Generate a list of lists of filenames for the climatology
        self.climatology_file_name_groups = self.create_filename_groups(climatology=True)

        # Set the threshold for sea ice fraction
        self.f_max = f_max,

        # Set the output paths
        self.out_path = out_path
        self.zip_name = zip_name

        # Set the date created and the uuid
        self.date_created = datetime.now().strftime('%Y-%m-%dT%H:%M:%SZ')
        self.uuid = str(uuid4())

        # Read in a sea ice fraction field
        fl = cf.read(self.filename_groups[0][0])
        sif = fl.select_by_property(standard_name='sea_ice_area_fraction')[0]

        # Create longitude and latitude bounds if necessary
        regridding_utilities.create_lonlat_bounds(sif)

        # Make the weights
        self.weights = sif.weights('area', scale=1.0)

        # Resample the longitude and latitude
        self.resample_lonlat(sif)

        # Set the timescale and spatial scale within which errors are assumed to be fully correlated in days
        self.tau = tau
        self.spatial_lambda = spatial_lambda

        self.calculate_k_xy()

        # Calculate the fraction of ocean in each target cell
        self.sf_data = self.spatially_resample_data(sif.where(True, 1.0), weights=True)
        self.sf_data /= self.spatially_resample_data(self.weights)
        self.sf_data.filled(fill_value=0.0, inplace=True)

        # Close the files
        fl.close()

        # Do the regridding
        self.make_gridded_ssts()

    def calculate_k_xy(self):
        """
        Calculate the spatial component of k for calculating the effective number of assumed-independent observations.
        """
        lon_bounds_array = self.lon.bounds.array
        lat_bounds_array = self.lat.bounds.array
        cos_theta = np.cos(np.radians(self.lat.array)).reshape(-1, 1)
        k_x_denominator = self.spatial_lambda / cos_theta
        k_x = np.maximum((lon_bounds_array[:, 1] - lon_bounds_array[:, 0]) / k_x_denominator, 1.0)
        k_y = np.maximum((lat_bounds_array[:, 1] - lat_bounds_array[:, 0]) / self.spatial_lambda, 1.0)
        self._k_xy = cf.Data(k_x.reshape((1, self.lat.size, self.lon.size)) * k_y.reshape((1, self.lat.size, 1)),
                             units='1')

    def check_dates(self):
        """
        Check the dates and time resolution are allowed.
        """
        if self.start_month < 1 or self.start_month > 12:
            raise ValueError('Start month out of range.')
        if self.end_month < 1 or self.end_month > 12:
            raise ValueError('End month out of range.')
        days_in_start_month = monthrange(self.year, self.start_month)[1]
        days_in_end_month = monthrange(self.year, self.end_month)[1]
        if self.start_day < 1 or self.start_day > days_in_start_month:
            raise ValueError('Start day out of range.')
        if self.end_day < 1 or self.end_day > days_in_end_month:
            raise ValueError('End day out of range.')
        if date(self.year, self.start_month, self.start_day) > date(self.year, self.end_month, self.end_day):
            raise ValueError('Start date cannot be greater than end date.')
        if self.time_resolution == 'annual':  # Annual resampling
            if self.start_month != 1 or self.start_day != 1 or self.end_month != 12 or self.end_day != 31:
                raise ValueError('A whole number of years is required for annual regridding.')
        elif self.time_resolution == 'monthly':  # Monthly resampling
            if self.start_day != 1 or self.end_day != days_in_end_month:
                raise ValueError('A whole number of months is required for monthly regridding.')
        elif self.time_resolution == '10-day':  # Nominal 10-day resampling
            if self.start_day not in [1, 11, 21] or self.end_day not in [10, 20, days_in_end_month]:
                raise ValueError('A whole number of dekads is required for 10-day regridding.')
        elif self.time_resolution == '5-day':  # Nominal 5-day resampling
            if self.start_day not in [1, 6, 11, 16, 21, 26] or \
                    self.end_day not in [5, 10, 15, 20, 25, days_in_end_month]:
                raise ValueError('A whole number of pentads is required for 5-day regridding.')
        elif type(self.time_resolution) is int:  # Regular N-day resampling
            if self.time_resolution < 1:
                raise ValueError('Time resolution out of range.')
            elif self.time_resolution > 1:
                if self.start_month != 1 or self.start_day != 1 or self.end_month != 12 or self.end_day != 31:
                    raise ValueError('A whole number of years is required for regular N-day regridding.')
        else:  # The time resolution is not recognised
            raise ValueError('Time resolution not recognised.')

    def make_gridded_ssts(self):
        """
        Calculate the regridded L4 SSTs.
        """
        output_paths = []
        # The main code to perform the regridding
        for filenames, climatology_filenames in zip(self.filename_groups, self.climatology_file_name_groups):
            resampled_sst_data = None
            sst_denominator = None
            resampled_sst_uncert_data = None
            n = None
            if self.sea_ice_fraction:
                resampled_sif_data = None
                sif_denominator = None
            if self.anomalies:
                resampled_sst_climatology_data = None

            for filename, climatology_filename in zip(filenames, climatology_filenames):
                # Read in data
                print("Reading: "+filename)
                fl = cf.read(filename, aggregate=False,extra='field_ancillary')

                # Select the SST and create longitude and latitude bounds if necessary
                sst = fl.select_by_property(standard_name='sea_water_temperature')[0]
                regridding_utilities.create_lonlat_bounds(sst)

                # Select the sea ice fraction and create longitude and latitude bounds if necessary
                sif = fl.select_by_property(standard_name='sea_ice_area_fraction')[0]

                # OR the sea ice fraction mask with the SST mask.  In years >= c3s_start_year, a few cells in the SST
                # field are masked where equivalent cells in other fields are not.
                # if self.year >= regridding_constants.c3s_start_year:
                #    np.ma.masked_where(np.ma.getmask(sst.array) | np.ma.getmask(sif.array), sif.varray, copy=False)

                regridding_utilities.create_lonlat_bounds(sif)

                # Select the SST uncertainty and create longitude and latitude bounds if necessary
                sst_uncert = fl.select_by_property(standard_name='sea_water_temperature standard_error')[0]

                # OR the SST uncertainty mask with the SST mask.  In years >= c3s_start_year, a few cells in the SST
                # field are masked where equivalent cells in other fields are not
                # if self.year >= regridding_constants.c3s_start_year:
                #    np.ma.masked_where(np.ma.getmask(sst.array) | np.ma.getmask(sst_uncert.array),
                #                       sst_uncert.varray, copy=False)

                regridding_utilities.create_lonlat_bounds(sst_uncert)

                # Calculate the regridded absolute SST
                data = self.spatially_resample_data(sst.where(sif > self.f_max, cf.masked), weights=True)
                if resampled_sst_data is None:
                    resampled_sst_data = data
                else:
                    resampled_sst_data = regridding_utilities.add_data(resampled_sst_data, data)

                # Calculate the denominator for averaging
                data = self.spatially_resample_data(sif.where(cf.le(self.f_max), 1.0, cf.masked), weights=True)
                if sst_denominator is None:
                    sst_denominator = data
                else:
                    sst_denominator = regridding_utilities.add_data(sst_denominator, data)

                # Calculate the regridded SST uncertainty
                data = self.spatially_resample_data(sst_uncert.where(sif > self.f_max, cf.masked) ** 2.0, weights=True)
                if resampled_sst_uncert_data is None:
                    resampled_sst_uncert_data = data
                else:
                    resampled_sst_uncert_data = regridding_utilities.add_data(resampled_sst_uncert_data, data)

                # Calculate the number of observations used in each target cell
                data = self.spatially_resample_data(sif.where(cf.le(self.f_max), 1.0, cf.masked))
                if n is None:
                    n = data
                else:
                    n = regridding_utilities.add_data(n, data)

                if self.sea_ice_fraction:
                    # Calculate the regridded sea ice fraction
                    data = self.spatially_resample_data(sif, weights=True)
                    if resampled_sif_data is None:
                        resampled_sif_data = data
                    else:
                        resampled_sif_data = regridding_utilities.add_data(resampled_sif_data, data)

                    # Calculate the denominator for averaging
                    data = self.spatially_resample_data(sif.where(True, 1.0), True)
                    if sif_denominator is None:
                        sif_denominator = data
                    else:
                        sif_denominator = regridding_utilities.add_data(sif_denominator, data)

                if self.anomalies:
                    # Read in and regrid the climatology
                    if type(climatology_filename) is str:
                        gl = cf.read(climatology_filename, aggregate=False)
                        sst_climatology = gl.select_by_property(standard_name='sea_water_temperature')[0]
                    else:
                        gl = cf.read(climatology_filename[0], aggregate=False)
                        sst_climatology = gl.select_by_property(standard_name='sea_water_temperature')[0]
                        hl = cf.read(climatology_filename[1], aggregate=False)
                        sst_climatology += hl.select_by_property(standard_name='sea_water_temperature')[0]
                        sst_climatology /= 2

                    # OR the climatology mask with the SST mask.  In years >= c3s_start_year, a few cells in the SST
                    # field are masked where equivalent cells in other fields are not
                    # if self.year >= regridding_constants.c3s_start_year:
                    #    np.ma.masked_where(np.ma.getmask(sst.array) | np.ma.getmask(sst_climatology.array),
                    #                       sst_climatology.varray, copy=False)

                    data = self.spatially_resample_data(sst_climatology.where(sif > self.f_max, cf.masked), weights=True)
                    if resampled_sst_climatology_data is None:
                        resampled_sst_climatology_data = data
                    else:
                        resampled_sst_climatology_data = regridding_utilities.add_data(resampled_sst_climatology_data, data)
                    gl.close()
                    if type(climatology_filename) is not str:
                        hl.close()

                fl.close()

            # Average the summed data
            resampled_sst_data /= sst_denominator
            resampled_sst_uncert_data /= sst_denominator
            if self.sea_ice_fraction:
                resampled_sif_data /= sif_denominator
            if self.anomalies:
                resampled_sst_climatology_data /= sst_denominator
                resampled_sst_data -= resampled_sst_climatology_data

            # Finalise calculation of the SST uncertainty
            k = self._k_xy * max(len(filenames) / self.tau, 1.0)
            n = np.minimum(k, n)
            resampled_sst_uncert_data /= n
            resampled_sst_uncert_data **= 0.5

            # Create the resampled sst field
            resampled_sst = self.update_field(filenames, 'sea_water_temperature', resampled_sst_data)
            if self.anomalies:
                resampled_sst.nc_set_variable('sst_anomaly')
                resampled_sst.standard_name = 'sea_water_temperature_anomaly'
                resampled_sst.long_name = 'regridded analysed sea surface temperature anomaly'
            else:
                resampled_sst.nc_set_variable('sst')
                resampled_sst.long_name = 'regridded analysed sea surface temperature'
            resampled_sst.comment = 'These data were produced by the University of Reading as part of the ESA ' + \
                                    'CCI project. Points where the sea ice fraction was above f_max were excluded. ' + \
                                    'Please refer to the global attributes for its value.'

            # Create the resampled sst uncertainty field
            resampled_sst_uncert = self.update_field(filenames, 'sea_water_temperature standard_error',
                                                     resampled_sst_uncert_data)
            resampled_sst_uncert.nc_set_variable('sst_uncertainty')
            resampled_sst_uncert.long_name = 'estimated error standard deviation of regridded analysed sea surface ' + \
                                             'temperature'
            resampled_sst_uncert.comment = 'These data were produced by the University of Reading as part of the ' + \
                                           'ESA CCI project. Points where the sea ice fraction was above f_max ' + \
                                           'were excluded. The approximate model is that errors are fully ' + \
                                           'correlated for timescales shorter than tau days and spacescales ' + \
                                           'shorter than lambda degrees. Please refer to the global attributes for ' + \
                                           'their values.'

            # Create the ocean fraction field
            sf = self.update_field(filenames, 'sea_ice_area_fraction', self.sf_data)
            sf.nc_set_variable('sea_fraction')
            sf.standard_name = 'sea_area_fraction'
            sf.long_name = 'sea area fraction'
            sf.comment = 'Fraction of cell that is covered by sea (including sea ice)'
            del sf.source

            # Create the resampled sif field
            if self.sea_ice_fraction:
                resampled_sif = self.update_field(filenames, 'sea_ice_area_fraction', resampled_sif_data)
                resampled_sif.long_name = 'regridded sea ice area fraction'
                resampled_sif.comment = 'Regridded sea ice area fraction'

            # Get the date of the resampled SST
            dt = resampled_sst.dimension_coordinate('T').datetime_array[0]

            # Create a field list with all the fields in it
            fl = cf.FieldList()
            fl.append(resampled_sst)
            fl.append(resampled_sst_uncert)
            if self.sea_ice_fraction:
                fl.append(resampled_sif)
            fl.append(sf)
            fl.append(regridding_utilities.create_time_field(sf, 'calendar_year', 'calendar year', dt.year))
            fl.append(regridding_utilities.create_time_field(sf, 'calendar_month', 'calendar month', dt.month))
            fl.append(regridding_utilities.create_time_field(sf, 'day_of_month', 'day of month', dt.day))
            fl.append(regridding_utilities.create_time_field(sf, 'day_of_year', 'day of year', dt.dayofyr))

            # Write the data
            output_path = os.path.join(self.out_path, dt.strftime('%Y%m%d') + '_regridded_sst.nc')
            cf.write(fl, output_path, datatype={np.dtype('float64'): np.dtype('float32')}, compress=1,
                     least_significant_digit=3)
            output_paths.append(output_path)

        # if specified, zip up all output files and remove the originals
        if self.zip_name:
            zip_path = os.path.join(self.out_path, self.zip_name)
            with zipfile.ZipFile(zip_path, "w") as z:
                for output_path in output_paths:
                    z.write(output_path, os.path.split(output_path)[1])
                    os.unlink(output_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Regrid L4 SST data for a particular year.')
    parser.add_argument('lon_resolution', type=float, help='The target longitude resolution. This must be a multiple ' +
                                                           'of 0.05 degrees and a factor of 360 degrees. It must ' +
                                                           'also not be greater than 5 degrees.')
    parser.add_argument('lat_resolution', type=float, help='The target latitude resolution. This must be a multiple ' +
                                                           'of 0.05 degrees and a factor of 180 degrees. It must ' +
                                                           'also not be greater than 5 degrees.')
    parser.add_argument('time_resolution', help="The target time resolution. This can be 'annual', 'monthly', " +
                                                "'10-day' for dekads, '5-day' for pentads or an integer for regular " +
                                                " N day regridding aligned with the start of the year.")
    parser.add_argument('--year', type=int, default=_default_year,
                        help='The year of the data to be regridded. Default is ' + str(_default_year) + '.')
    parser.add_argument('--start_month', type=int, default=_default_start_month,
                        help='The first month of the data to be regridded. Default is ' + str(_default_start_month)
                             + '.')
    parser.add_argument('--start_day', type=int, default=_default_start_day,
                        help='The first day of the month of the data to be regridded. Default is '
                             + str(_default_start_day) + '.')
    parser.add_argument('--end_month', type=int, default=_default_end_month,
                        help='The final month of the data to be regridded. Default is ' + str(_default_end_month) + '.')
    parser.add_argument('--end_day', type=int, default=_default_end_day,
                        help='The final day of the month of the data to be regridded. Default is '
                             + str(_default_end_day) + '.')
    parser.add_argument('--anomalies', action='store_true', default=False,
                        help='Output anomalies instead of absolute SSTs.')
    parser.add_argument('--no_sea_ice_fraction', action='store_true', default=False,
                        help='Do not output the sea ice fraction.')
    parser.add_argument('--f_max', type=float, default=_default_f_max,
                        help='The fraction of sea ice above which a cell is ignored in calculating the regridded SST '
                             + 'SST. Default is ' + str(_default_f_max) + '. When calculating anomalies, the '
                             + 'climatology is always calculated with an f_max of 1.0 regardless of this value.')
    parser.add_argument('--tau', type=int, default=_default_tau,
                        help='Timescale within which errors are assumed to be fully correlated in days. '
                             + 'Default is ' + str(_default_tau) + ' days.')
    parser.add_argument('--spatial_lambda', type=float, default=_default_spatial_lambda,
                        help='Spatial scale within which errors are assumed to be fully correlated in degrees. '
                             + 'Default is ' + str(_default_spatial_lambda) + ' degrees.')
    parser.add_argument('--sst_analysis_l4_path', default=regridding_constants.default_sst_analysis_l4_path,
                        help='Path to the SST Analysis Level 4 input data.')
    parser.add_argument('--sst_climatology_path', default=regridding_constants.default_sst_climatology_path,
                        help='Path to the SST CCI Climatology input data.')
    parser.add_argument('--out_path', default=_default_out_path,
                        help='The path in which to write the output file of regridded data.')
    parser.add_argument('--zip_name', default=_default_zip_name,
                        help='Combine all output files into a single zip file with this name.')

    args = parser.parse_args()

    try:
        time_resolution = int(args.time_resolution)
    except ValueError:
        time_resolution = args.time_resolution

    sst_regridder = SSTRegridder(sst_analysis_l4_path=args.sst_analysis_l4_path,
                                 sst_climatology_path=args.sst_climatology_path)

    sst_regridder(args.lon_resolution, args.lat_resolution, time_resolution, year=args.year,
                  start_month=args.start_month, start_day=args.start_day, end_month=args.end_month,
                  end_day=args.end_day, anomalies=args.anomalies, no_sea_ice_fraction=args.no_sea_ice_fraction,
                  f_max=args.f_max, tau=args.tau, spatial_lambda=args.spatial_lambda, out_path=args.out_path,
                  zip_name=args.zip_name)
