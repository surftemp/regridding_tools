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

_default_year = regridding_constants._cci_start_year
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


class SSTRegridder(object):
    """
    SST regridder class. Should be instantiated first, optionally with paths to the input data and then a call made
    to the resulting object to perform the regridding at a particular resolution with other options available.
    """
    def __init__(self, sst_cci_analysis_l4_path=regridding_constants._default_sst_cci_analysis_l4_path,
                 c3s_sst_analysis_l4_path=regridding_constants._default_c3s_sst_analysis_l4_path,
                 sst_cci_climatology_path=regridding_constants._default_sst_cci_climatology_path, ):
        """
        Initiate an SST regridder optionally specfying the input file paths.

        :param sst_cci_climatology_path:
            Path to the SST CCI Climatology input data.

        :param c3s_sst_analysis_l4_path:
            Path to the C3S SST Analysis Level 4 input data.

        :param sst_cci_analysis_l4_path:
            Path to the SST CCI Analysis Level 4 input data.
        """
        self._sst_cci_analysis_l4_path = sst_cci_analysis_l4_path
        self._c3s_sst_analysis_l4_path = c3s_sst_analysis_l4_path
        self._sst_cci_climatology_path = sst_cci_climatology_path

        self._base_resolution = 0.05
        self._max_resolution = 5.0
        self._lon_range = 360.0
        self._lat_range = 180.0

        self._allowed_lon_resolutions = SSTRegridder._allowed_lonlat_resolutions(self._base_resolution,
                                                                                 self._max_resolution,
                                                                                 self._lon_range)
        self._allowed_lat_resolutions = SSTRegridder._allowed_lonlat_resolutions(self._base_resolution,
                                                                                 self._max_resolution,
                                                                                 self._lat_range)

    def __call__(self, lon_resolution, lat_resolution, time_resolution, year=_default_year,
                 start_month=_default_start_month, start_day=_default_start_day,
                 end_month=_default_end_month, end_day=_default_end_day, anomalies=_default_anomalies,
                 no_sea_ice_fraction=_default_no_sea_ice_fraction, f_max=_default_f_max, tau=_default_tau,
                 spatial_lambda=_default_spatial_lambda, out_path=_default_out_path, zip_name=_default_zip_name):
        """
        Regrid the SSTs to a target resolution.

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
        if lon_resolution not in self._allowed_lon_resolutions:
            raise ValueError('lon_resolution must be a multiple of ' + str(self._base_resolution) + ', a factor of ' +
                             str(self._lon_range) + ' and not greater than ' + str(self._max_resolution) + ' degrees.')
        if lat_resolution not in self._allowed_lat_resolutions:
            raise ValueError('lat_resolution must be a multiple of ' + str(self._base_resolution) + ', a factor of ' +
                             str(self._lat_range) + ' and not greater than ' + str(self._max_resolution) + ' degrees.')
        self._lon_step = self._allowed_lon_resolutions[lon_resolution]
        self._lat_step = self._allowed_lat_resolutions[lat_resolution]
        self._lon_resolution = lon_resolution
        self._lat_resolution = lat_resolution

        # Set the dates and time resolution
        self._year = year
        self._start_month = start_month
        self._start_day = start_day
        self._end_month = end_month
        self._end_day = end_day
        self._time_resolution = time_resolution

        # Set flags for whether to calculate anomalies or the sea ice fraction
        self._anomalies = anomalies
        self._sea_ice_fraction = not no_sea_ice_fraction

        # Check the dates and time resolution
        self._check_dates()

        # Generate a list of lists of filenames
        self._filename_groups = self._create_filename_groups()

        # Generate a list of lists of filenames for the climatology
        self._climatology_file_name_groups = self._create_filename_groups(climatology=True)

        # Set the threshold for sea ice fraction
        self._f_max = f_max,

        # Set the output paths
        self._out_path = out_path
        self._zip_name = zip_name

        # Set the date created and the uuid
        self._date_created = datetime.now().strftime('%Y-%m-%dT%H:%M:%SZ')
        self._uuid = str(uuid4())

        # Read in a sea ice fraction field
        fl = cf.read(self._filename_groups[0][0])
        sif = fl.select_by_property(standard_name='sea_ice_area_fraction')[0]

        # Create longitude and latitude bounds if necessary
        SSTRegridder._create_lonlat_bounds(sif)

        # Make the weights
        self._weights = sif.weights('area', scale=1.0)

        # Resample the longitude and latitude
        self._resample_lonlat(sif)

        # Set the timescale and spatial scale within which errors are assumed to be fully correlated in days
        self._tau = tau
        self._spatial_lambda = spatial_lambda

        self._calculate_k_xy()

        # Calculate the fraction of ocean in each target cell
        self._sf_data = self._spatially_resample_data(sif.where(True, 1.0), weights=True)
        self._sf_data /= self._spatially_resample_data(self._weights)
        self._sf_data.filled(fill_value=0.0, inplace=True)

        # Close the files
        fl.close()

        # Do the regridding
        self.make_gridded_ssts()

    @staticmethod
    def _allowed_lonlat_resolutions(base_resolution, max_resolution, lonlat_range):
        """
        Calculate the allowed resolutions for averaging over longitude or latitude given the resolution and range of the
        longitude or latitude of the input grid in degrees. The new resolution must be a multiple of the base resolution
        and a factor of the range.

        :param base_resolution:
            The resolution of the longitude or latitude of the input grid in degrees.

        :param max_resolution:
            The maximum resolution for the longitude or latitude in degrees.

        :param lonlat_range:
            The range of the longitude or latitude of the input grid in degrees.

        :return:
            A list of the allowed longitude or latitude steps for averaging in degrees.
        """
        allowed_resolutions = {}
        size = int(round(lonlat_range / base_resolution))
        for step in range(1, size):
            resolution = round(step * base_resolution, 2)
            if resolution > max_resolution:
                break
            if size % step == 0:
                allowed_resolutions[resolution] = step
        return allowed_resolutions

    @staticmethod
    def _create_lonlat_bounds(f):
        """
        Create longitude and latitude bounds for a field if necessary

        :param f:
            The field for which to create longitude and latitude bounds
        """
        lon = f.dimension_coordinate('X')
        if not lon.has_bounds():
            lon_bounds = lon.create_bounds(cellsize=0.05, min=-180.0, max=180.0)
            lon_bounds.nc_set_variable('lon_bnds')
            lon_bounds.nc_set_dimension('bnds')
            lon.set_bounds(lon_bounds)
        lat = f.dimension_coordinate('Y')
        if not lat.has_bounds():
            lat_bounds = lat.create_bounds(cellsize=0.05, min=-90.0, max=90.0)
            lat_bounds.nc_set_variable('lat_bnds')
            lat_bounds.nc_set_dimension('bnds')
            lat.set_bounds(lat_bounds)

    def _resample_lonlat(self, f):
        """
        Create resampled longitude and latitude cf dimension coordinates.

        :param f:
             cf Field object with longitudes and latitudes.
        """
        # Resample the latitude and longitude arrays and their bounds
        lon = f.dimension_coordinate('X')
        lat = f.dimension_coordinate('Y')

        lon_bounds_array = lon.bounds.array
        lat_bounds_array = lat.bounds.array

        resampled_lon_bounds_array = np.column_stack((lon_bounds_array[0::self._lon_step, 0],
                                                      lon_bounds_array[self._lon_step - 1::self._lon_step, 1]))
        resampled_lat_bounds_array = np.column_stack((lat_bounds_array[0::self._lat_step, 0],
                                                      lat_bounds_array[self._lat_step - 1::self._lat_step, 1]))

        # If the step is even use the value of the bounds of the old cell in the centre of the new cell as the new
        # coordinate values, otherwise use the coordinate values in the centre of the old cell.
        if self._lon_step % 2 == 0:
            resampled_lon_array = lon_bounds_array[self._lon_step // 2::self._lon_step, 0]
        else:
            resampled_lon_array = lon.array[self._lon_step // 2::self._lon_step]
        if self._lat_step % 2 == 0:
            resampled_lat_array = lat_bounds_array[self._lat_step // 2::self._lat_step, 0]
        else:
            resampled_lat_array = lat.array[self._lat_step // 2::self._lat_step]

        # Create new cf dimension coordinates
        lon_units = lon.units
        lat_units = lat.units
        lon_dtype = lon.dtype
        lat_dtype = lat.dtype

        resampled_lon = lon.copy(data=False)
        resampled_lat = lat.copy(data=False)

        resampled_lon_bounds = lon.bounds.copy(data=False)
        resampled_lat_bounds = lat.bounds.copy(data=False)

        resampled_lon_bounds.set_data(cf.Data(resampled_lon_bounds_array, units=lon_units, dtype=lon_dtype))
        resampled_lat_bounds.set_data(cf.Data(resampled_lat_bounds_array, units=lat_units, dtype=lat_dtype))

        resampled_lon.set_data(cf.Data(resampled_lon_array, units=lon_units, dtype=lon_dtype))
        resampled_lon.set_bounds(resampled_lon_bounds)
        resampled_lat.set_data(cf.Data(resampled_lat_array, units=lat_units, dtype=lat_dtype))
        resampled_lat.set_bounds(resampled_lat_bounds)

        self._lon = resampled_lon
        self._lat = resampled_lat

    def _calculate_k_xy(self):
        """
        Calculate the spatial component of k for calculating the effective number of assumed-independent observations.
        """
        lon_bounds_array = self._lon.bounds.array
        lat_bounds_array = self._lat.bounds.array
        cos_theta = np.cos(np.radians(self._lat.array)).reshape(-1, 1)
        k_x_denominator = self._spatial_lambda / cos_theta
        k_x = np.maximum((lon_bounds_array[:, 1] - lon_bounds_array[:, 0]) / k_x_denominator, 1.0)
        k_y = np.maximum((lat_bounds_array[:, 1] - lat_bounds_array[:, 0]) / self._spatial_lambda, 1.0)
        self._k_xy = cf.Data(k_x.reshape((1, self._lat.size, self._lon.size)) * k_y.reshape((1, self._lat.size, 1)),
                             units='1')


    def _spatially_resample_data(self, f, weights=False):
        """
        Resample a field's data spatially.

        :param f:
            The field whose data is to be resampled.

        :return:
            The resampled data as a cf.Data object.
        """
        if weights:
            array = (f * self._weights).array
        else:
            array = f.array

        array = array.reshape((1, self._lat.size, self._lat_step, self._lon.size, self._lon_step))
        array = array.transpose((0, 1, 3, 2, 4))
        array = array.reshape((1, self._lat.size, self._lon.size, self._lat_step * self._lon_step))
        array = np.sum(array, axis=3)

        data = cf.Data(array, units=f.units, fill_value=f.fill_value(default='netCDF'), dtype=f.dtype)
        return data

    @staticmethod
    def add_data(d1, d2):
        """
        Add two Data objects together such that masked points are ignored, but a masked point plus a non-masked point
        equals the value of the non-masked point. The data objects are concatenated along the first dimension, which is
        expected to be the size 1 time dimension.

        :param d1:
            The first operand as a cf.Data object.

        :param d2:
            The second operand as a cf.Data object.

        :return:
            The result of summing the two operands as a cf.Data object.
        """
        # Concatenate the data arrays together along the size 1 time dimension.
        e = cf.Data.concatenate([d1, d2], axis=0)

        # Return the sum along the time dimension.
        return e.sum(axes=0)

    def _update_field(self, filenames, standard_name, resampled_data):
        """
        Create an updated field with new data, longitudes, latitudes and a modified time stamp.

        :param filenames:
            The filenames in the time series.

        :param standard_name:
            The standard name of the field to be updated.

        :param resampled_data:
            The new data.

        :return:
            The new field.
        """
        # Get the first and the last field
        fl = cf.read(filenames[0], aggregate=False)
        first_field = fl.select_by_property(standard_name=standard_name)[0]

        gl = cf.read(filenames[-1], aggregate=False)
        last_field = gl.select_by_property(standard_name=standard_name)[0]

        # Copy the field without the data
        g = first_field.copy(data=False)

        # Replace the longitudes and latitudes with the resampled ones
        lon_key = g.dimension_coordinate('X', key=True)
        lat_key = g.dimension_coordinate('Y', key=True)
        lon_axis = g.domain_axis(lon_key, key=True)
        lat_axis = g.domain_axis(lat_key, key=True)
        g.domain_axes[lon_axis].set_size(self._lon.size)
        g.domain_axes[lat_axis].set_size(self._lat.size)
        g.set_construct(self._lon, axes=(lon_key,))
        g.set_construct(self._lat, axes=(lat_key,))

        # Get the time coordinates
        time = first_field.dimension_coordinate('T').copy()
        time_units = time.units
        time_dtype = time.dtype
        dtarray = time.datetime_array
        if time.has_bounds():
            bounds_dtarray = time.bounds.datetime_array
        else:
            bounds = time.create_bounds(cellsize=cf.D())
            bounds.nc_set_variable('time_bnds')
            bounds.nc_set_dimension('bnds')
            bounds_dtarray = bounds.datetime_array
            time.set_bounds(bounds)
        dt = dtarray[0]
        start_dt = bounds_dtarray[0, 0]
        time2 = last_field.dimension_coordinate('T')
        if time2.has_bounds():
            end_dt = time2.bounds.datetime_array[0, 1]
        else:
            end_dt = time2.create_bounds(cellsize=cf.D()).datetime_array[0, 1]

        # Modify the time stamps
        if self._time_resolution == 'annual':
            dtarray[0] = cf.dt(dt.year, 7, 1, calendar='gregorian')
        elif self._time_resolution == 'monthly':
            dtarray[0] = cf.dt(dt.year, dt.month, 15, 12, calendar='gregorian')
        elif self._time_resolution == '10-day':
            day = dt.day
            if day <= 10:
                new_day = 5
            elif day <= 20:
                new_day = 15
            else:
                new_day = 25
            dtarray[0] = cf.dt(dt.year, dt.month, new_day, 12, calendar='gregorian')
        elif self._time_resolution == '5-day':
            day = dt.day
            if day <= 5:
                new_day = 3
            elif day <= 10:
                new_day = 8
            elif day <= 15:
                new_day = 13
            elif day <= 20:
                new_day = 18
            elif day <= 25:
                new_day = 23
            else:
                if dt.month == 2:
                    new_day = 27
                else:
                    new_day = 28
            dtarray[0] = cf.dt(dt.year, dt.month, new_day, 12, calendar='gregorian')
        else:
            dtarray[0] = start_dt + ((end_dt - start_dt) / 2)
        bounds_dtarray[0, 1] = end_dt

        # Insert the modified time stamps
        time.set_data(cf.Data(dtarray, units=time_units, dtype=time_dtype))
        time.bounds.set_data(cf.Data(bounds_dtarray, units=time_units, dtype=time_dtype))

        # Insert the updated time coordinate
        time_key = g.dimension_coordinate('T', key=True)
        time_axis = g.domain_axis(time_key, key=True)
        g.set_construct(time, axes=(time_key,))

        # Insert the data in the new field
        g.set_data(resampled_data, axes=(time_axis, lat_axis, lon_axis))

        # Remove valid_min and valid_max
        del g.valid_min
        del g.valid_max

        # Update time related metadata
        g.set_property('stop_time', last_field.get_property('stop_time'))
        g.set_property('time_coverage_end', last_field.get_property('time_coverage_end'))
        if self._time_resolution == 'annual':
            time_coverage_duration = 'P1Y'
        elif self._time_resolution == 'monthly':
            time_coverage_duration = 'P1M'
        else:
            time_coverage_duration = 'P' + str(len(filenames)) + 'D'
        g.set_property('time_coverage_duration', time_coverage_duration)
        g.set_property('time_coverage_resolution', time_coverage_duration)

        # Change global attributes
        g.set_property('geospatial_lon_resolution', np.float32(self._lon_resolution))
        g.set_property('geospatial_lat_resolution', np.float32(self._lat_resolution))
        if self._lon_resolution == self._lat_resolution:
            g.set_property('spatial_resolution', str(self._lon_resolution) + ' degree')
        else:
            g.set_property('spatial_resolution', str(self._lon_resolution) + ' x ' + str(self._lat_resolution) +
                           ' degree')
        g.set_property('creator_processing_institution', 'These data were produced by the University of ' +
                       'Reading as part of the ESA CCI project.')
        g.set_property('acknowledgment', 'Regridding service funded by National Centre for Earth ' +
                       'Observation, UK')

        # Tweak metadata as discussed in https://github.com/surftemp/sst-services/issues/17
        g.set_property('summary', 'Regridding of L4 analysis product from the ESA SST CCI project')
        g.set_property('creator_email','c.j.merchant@reading.ac.uk')
        g.set_property('publisher_email', 'c.j.merchant@reading.ac.uk')
        g.set_property('references','Merchant, C.J., Embury, O., Bulgin, C.E., Block, T., Corlett, G.K., Fiedler, E., Good, S.A., Mittaz, J., Rayner, N.A., Berry, D., Eastwood, S., Taylor, M., Tsushima, Y., Waterfall, A., Wilson, R. and Donlon, C. (2019), Satellite-based time-series of sea-surface temperature since 1981 for climate applications. Scientific Data 6, 223, doi:10.1038/s41597-019-0236-x, and Merchant, C. J. and Embury, O. (2020) Adjusting for desert-dust-related biases in a climate data record of sea surface temperature. Remote Sensing, 12 (16). 2554. ISSN 2072-4292 doi:10.3390/rs12162554')
        g.set_property('creator_url','http://climate.esa.int/')
        g.set_property('metadata_link','http://climate.esa.int/')
        g.set_property('title','ESA SST CCI L4 analysis')
        g.set_property('history','Created by makegriddedSSTs.py v1.0 operating on SST CCI L4 analysis products')

        g.del_property('product_specification_version')
        g.nc_set_global_attribute('f_max', self._f_max)
        g.nc_set_global_attribute('lambda', self._spatial_lambda)
        g.nc_set_global_attribute('tau', np.int32(self._tau))
        g.set_property('date_created', self._date_created)
        g.set_property('uuid', self._uuid)
        g.set_property('tracking_id', self._uuid)

        # Close the fields
        fl.close()
        gl.close()

        return g

    @staticmethod
    def _create_time_field(f, variable_name, long_name, data_value):
        """
        Create "convenience" fields for calendar year, calendar month, day of month and day of year.

        :param f:
            A field whose metadata will be copied.

        :param variable_name:
            The name of the netCDF variable.

        :param long_name:
            The long name of the variable.

        :param data_value:
            The value of the variable.

        :return:
            A new cf.Field object.
        """
        # Field's data is copied lazily and deleted before being used
        g = f.copy()
        g.del_data()
        g.del_data_axes()

        # Get the identifiers for the dimension coordinates
        lon_key = g.dimension_coordinate('X', key=True)
        lat_key = g.dimension_coordinate('Y', key=True)

        # Get the identifiers for the corresponding domain axes
        lon_axis = g.domain_axis(lon_key, key=True)
        lat_axis = g.domain_axis(lat_key, key=True)

        # Delete the latitude and longitude coordinates
        g.del_construct(lon_key)
        g.del_construct(lat_key)

        # Delete the corresponding domain axes
        g.del_construct(lon_axis)
        g.del_construct(lat_axis)

        # Insert the data
        time_key = g.dimension_coordinate('T', key=True)
        g.set_data(cf.Data(np.array([data_value]), dtype='int32'), axes=(time_key,))

        # Update metadata
        del g.standard_name
        g.nc_set_variable(variable_name)
        g.long_name = long_name
        del g.comment

        return g

    @staticmethod
    def daterange(start_date, end_date, step=1):
        """
        Generator to yield a range of dates a day at a time. Note that the end date is not included to be consistent
        with the internal range function.

        Code modified from: https://stackoverflow.com/questions/1060279/iterating-through-a-range-of-dates-in-python

        :param start_date:
            The first date in the range as a date object.

        :param end_date:
            The last date in the range as a date object.

        :param step
            The step in days.

        :return:
            Yield a date object.
        """
        for n in range(0, int((end_date - start_date).days), step):
            yield start_date + timedelta(n)

    def _get_filename(self, d, climatology):
        """
        Generate the name of the input file for a particular date.

        :param d:
            Python data object.

        :param climatology:
            Boolean indicating whether to return filenames for the climatology (if True) or the analysed SST (if False).

        :return:
            The filename as a string.
        """
        if climatology:
            if d.month == 2 and d.day == 29:
                # If the day is the 29th of February the files for D59 and D60 will be averaged.
                filename = (os.path.join(self._sst_cci_climatology_path,
                                         'D059-ESACCI-L4_GHRSST-SSTdepth-OSTIA-GLOB_CDR2.1-v02.0-fv01.0.nc'),
                            os.path.join(self._sst_cci_climatology_path,
                                         'D060-ESACCI-L4_GHRSST-SSTdepth-OSTIA-GLOB_CDR2.1-v02.0-fv01.0.nc'))
            else:
                # Get day of year as if the year was a non-leap year.
                day_of_year = date(1982, d.month, d.day).timetuple().tm_yday
                filename = (os.path.join(self._sst_cci_climatology_path, 'D' + str(day_of_year).zfill(3) +
                                         '-ESACCI-L4_GHRSST-SSTdepth-OSTIA-GLOB_CDR2.1-v02.0-fv01.0.nc'))
        else:
            if d.year >= regridding_constants._c3s_start_year:
                filename = os.path.join(self._c3s_sst_analysis_l4_path,
                                        d.strftime('%Y'), d.strftime('%m'), d.strftime('%d'), d.strftime('%Y%m%d') +
                                        '120000-C3S-L4_GHRSST-SSTdepth-OSTIA-GLOB_ICDR2.0-v02.0-fv01.0.nc')
            else:
                filename = os.path.join(self._sst_cci_analysis_l4_path,
                                        d.strftime('%Y'), d.strftime('%m'), d.strftime('%d'), d.strftime('%Y%m%d') +
                                        '120000-ESACCI-L4_GHRSST-SSTdepth-OSTIA-GLOB_CDR2.1-v02.0-fv01.0.nc')
        return filename

    def _get_start_end_dates(self, month):
        """
        Get the start dates and end dates for a particular month of the year.

        :param month:
            The month of the year.

        :return:
            The start date and end date as Python date objects.
        """
        days_in_month = monthrange(self._year, month)[1]
        if month == self._start_month == self._end_month:
            start_date = date(self._year, month, self._start_day)
            end_date = date(self._year, month, self._end_day)
        elif month == self._start_month:
            start_date = date(self._year, month, self._start_day)
            end_date = date(self._year, month, days_in_month)
        elif month == self._end_month:
            start_date = date(self._year, month, 1)
            end_date = date(self._year, month, self._end_day)
        else:
            start_date = date(self._year, month, 1)
            end_date = date(self._year, month, days_in_month)
        return start_date, end_date + timedelta(1)

    def _check_dates(self):
        """
        Check the dates and time resolution are allowed.
        """
        if self._year < 1981:
            raise ValueError('Year out of range.')
        if self._start_month < 1 or self._start_month > 12:
            raise ValueError('Start month out of range.')
        if self._end_month < 1 or self._end_month > 12:
            raise ValueError('End month out of range.')
        days_in_start_month = monthrange(self._year, self._start_month)[1]
        days_in_end_month = monthrange(self._year, self._end_month)[1]
        if self._start_day < 1 or self._start_day > days_in_start_month:
            raise ValueError('Start day out of range.')
        if self._end_day < 1 or self._end_day > days_in_end_month:
            raise ValueError('End day out of range.')
        if date(self._year, self._start_month, self._start_day) > date(self._year, self._end_month, self._end_day):
            raise ValueError('Start date cannot be greater than end date.')
        if self._time_resolution == 'annual':  # Annual resampling
            if self._start_month != 1 or self._start_day != 1 or self._end_month != 12 or self._end_day != 31:
                raise ValueError('A whole number of years is required for annual regridding.')
        elif self._time_resolution == 'monthly':  # Monthly resampling
            if self._start_day != 1 or self._end_day != days_in_end_month:
                raise ValueError('A whole number of months is required for monthly regridding.')
        elif self._time_resolution == '10-day':  # Nominal 10-day resampling
            if self._start_day not in [1, 11, 21] or self._end_day not in [10, 20, days_in_end_month]:
                raise ValueError('A whole number of dekads is required for 10-day regridding.')
        elif self._time_resolution == '5-day':  # Nominal 5-day resampling
            if self._start_day not in [1, 6, 11, 16, 21, 26] or \
                    self._end_day not in [5, 10, 15, 20, 25, days_in_end_month]:
                raise ValueError('A whole number of pentads is required for 5-day regridding.')
        elif type(self._time_resolution) is int:  # Regular N-day resampling
            if self._time_resolution < 1:
                raise ValueError('Time resolution out of range.')
            elif self._time_resolution > 1:
                if self._start_month != 1 or self._start_day != 1 or self._end_month != 12 or self._end_day != 31:
                    raise ValueError('A whole number of years is required for regular N-day regridding.')
        else:  # The time resolution is not recognised
            raise ValueError('Time resolution not recognised.')

    def _create_filename_groups(self, climatology=False):
        """
        Create a list of lists of filenames. Each group of filenames will be averaged together over time. It is assumed
        that the dates and time resolution have already been checked for compatibility.

        :param climatology:
            Boolean indicating whether to return filenames for the climatology (if True) or the analysed SST (if False).

        :return:
            A list of lists of filenames.
        """
        filename_groups = []
        if self._time_resolution == 'annual':
            filenames = []
            start_date = date(self._year, self._start_month, self._start_day)
            end_date = date(self._year, self._end_month, self._end_day) + timedelta(1)
            for d in SSTRegridder.daterange(start_date, end_date):
                filename = self._get_filename(d, climatology)
                filenames.append(filename)
            filename_groups.append(filenames)
        elif self._time_resolution == 'monthly':
            for month in range(self._start_month, self._end_month + 1):
                start_date = date(self._year, month, 1)
                days_in_month = monthrange(self._year, month)[1]
                end_date = date(self._year, month, days_in_month) + timedelta(1)
                filenames = []
                for d in SSTRegridder.daterange(start_date, end_date):
                    filename = self._get_filename(d, climatology)
                    filenames.append(filename)
                filename_groups.append(filenames)
        elif self._time_resolution == '10-day' or self._time_resolution == '5-day':
            step = 10 if self._time_resolution == '10-day' else 5
            break_day = 21 if self._time_resolution == '10-day' else 26
            for month in range(self._start_month, self._end_month + 1):
                days_in_month = monthrange(self._year, month)[1]
                start_date, end_date = self._get_start_end_dates(month)
                for outer_d in SSTRegridder.daterange(start_date, end_date, step=step):
                    if outer_d.day < break_day:
                        inner_drange = SSTRegridder.daterange(outer_d, outer_d + timedelta(step))
                    else:
                        inner_drange = SSTRegridder.daterange(outer_d, date(self._year, month, days_in_month) +
                                                              timedelta(1))
                    filenames = []
                    for inner_d in inner_drange:
                        filename = self._get_filename(inner_d, climatology)
                        filenames.append(filename)
                    filename_groups.append(filenames)
                    if outer_d.day >= break_day:
                        break
        else:
            start_date = date(self._year, self._start_month, self._start_day)
            end_date = date(self._year, self._end_month, self._end_day) + timedelta(1)
            days_in_year = 366 if isleap(self._year) else 365
            finish = False
            for i in range(0, days_in_year, self._time_resolution):
                if i + self._time_resolution >= days_in_year:
                    drange = SSTRegridder.daterange(start_date + timedelta(i), end_date)
                elif i + 3 * self._time_resolution / 2 > days_in_year:
                    drange = SSTRegridder.daterange(start_date + timedelta(i), end_date)
                    finish = True
                else:
                    drange = SSTRegridder.daterange(start_date + timedelta(i),
                                                    start_date + timedelta(i + self._time_resolution))
                filenames = []
                for d in drange:
                    filename = self._get_filename(d, climatology)
                    filenames.append(filename)
                filename_groups.append(filenames)
                if finish:
                    break
        return filename_groups

    def make_gridded_ssts(self):
        """
        Calculate the regridded SSTs.
        """
        output_paths = []
        # The main code to perform the regridding
        for filenames, climatology_filenames in zip(self._filename_groups, self._climatology_file_name_groups):
            resampled_sst_data = None
            sst_denominator = None
            resampled_sst_uncert_data = None
            n = None
            if self._sea_ice_fraction:
                resampled_sif_data = None
                sif_denominator = None
            if self._anomalies:
                resampled_sst_climatology_data = None

            for filename, climatology_filename in zip(filenames, climatology_filenames):
                # Read in data
                fl = cf.read(filename, aggregate=False)

                # Select the SST and create longitude and latitude bounds if necessary
                sst = fl.select_by_property(standard_name='sea_water_temperature')[0]
                SSTRegridder._create_lonlat_bounds(sst)

                # Select the sea ice fraction and create longitude and latitude bounds if necessary
                sif = fl.select_by_property(standard_name='sea_ice_area_fraction')[0]

                # OR the sea ice fraction mask with the SST mask.  In years >= _c3s_start_year, a few cells in the SST
                # field are masked where equivalent cells in other fields are not.
                if self._year >= regridding_constants._c3s_start_year:
                    np.ma.masked_where(np.ma.getmask(sst.array) | np.ma.getmask(sif.array), sif.varray, copy=False)

                SSTRegridder._create_lonlat_bounds(sif)

                # Select the SST uncertainty and create longitude and latitude bounds if necessary
                sst_uncert = fl.select_by_property(standard_name='sea_water_temperature standard_error')[0]

                # OR the SST uncertainty mask with the SST mask.  In years >= _c3s_start_year, a few cells in the SST
                # field are masked where equivalent cells in other fields are not
                if self._year >= regridding_constants._c3s_start_year:
                    np.ma.masked_where(np.ma.getmask(sst.array) | np.ma.getmask(sst_uncert.array),
                                       sst_uncert.varray, copy=False)

                SSTRegridder._create_lonlat_bounds(sst_uncert)

                # Calculate the regridded absolute SST
                data = self._spatially_resample_data(sst.where(sif > self._f_max, cf.masked), weights=True)
                if resampled_sst_data is None:
                    resampled_sst_data = data
                else:
                    resampled_sst_data = SSTRegridder.add_data(resampled_sst_data, data)

                # Calculate the denominator for averaging
                data = self._spatially_resample_data(sif.where(cf.le(self._f_max), 1.0, cf.masked), weights=True)
                if sst_denominator is None:
                    sst_denominator = data
                else:
                    sst_denominator = SSTRegridder.add_data(sst_denominator, data)

                # Calculate the regridded SST uncertainty
                data = self._spatially_resample_data(sst_uncert.where(sif > self._f_max, cf.masked) ** 2.0,
                                                     weights=True)
                if resampled_sst_uncert_data is None:
                    resampled_sst_uncert_data = data
                else:
                    resampled_sst_uncert_data = SSTRegridder.add_data(resampled_sst_uncert_data, data)

                # Calculate the number of observations used in each target cell
                data = self._spatially_resample_data(sif.where(cf.le(self._f_max), 1.0, cf.masked))
                if n is None:
                    n = data
                else:
                    n = SSTRegridder.add_data(n, data)

                if self._sea_ice_fraction:
                    # Calculate the regridded sea ice fraction
                    data = self._spatially_resample_data(sif, weights=True)
                    if resampled_sif_data is None:
                        resampled_sif_data = data
                    else:
                        resampled_sif_data = SSTRegridder.add_data(resampled_sif_data, data)

                    # Calculate the denominator for averaging
                    data = self._spatially_resample_data(sif.where(True, 1.0), True)
                    if sif_denominator is None:
                        sif_denominator = data
                    else:
                        sif_denominator = SSTRegridder.add_data(sif_denominator, data)

                if self._anomalies:
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

                    # OR the climatology mask with the SST mask.  In years >= _c3s_start_year, a few cells in the SST
                    # field are masked where equivalent cells in other fields are not
                    if self._year >= regridding_constants._c3s_start_year:
                        np.ma.masked_where(np.ma.getmask(sst.array) | np.ma.getmask(sst_climatology.array),
                                           sst_climatology.varray, copy=False)

                    data = self._spatially_resample_data(sst_climatology.where(sif > self._f_max, cf.masked),
                                                         weights=True)
                    if resampled_sst_climatology_data is None:
                        resampled_sst_climatology_data = data
                    else:
                        resampled_sst_climatology_data = SSTRegridder.add_data(resampled_sst_climatology_data, data)
                    gl.close()
                    if type(climatology_filename) is not str:
                        hl.close()

                fl.close()

            # Average the summed data
            resampled_sst_data /= sst_denominator
            resampled_sst_uncert_data /= sst_denominator
            if self._sea_ice_fraction:
                resampled_sif_data /= sif_denominator
                resampled_sif_data.filled(fill_value=0.0, inplace=True)
            if self._anomalies:
                resampled_sst_climatology_data /= sst_denominator
                resampled_sst_data -= resampled_sst_climatology_data

            # Finalise calculation of the SST uncertainty
            k = self._k_xy * max(len(filenames) / self._tau, 1.0)
            n = np.minimum(k, n)
            resampled_sst_uncert_data /= n
            resampled_sst_uncert_data **= 0.5

            # Create the resampled sst field
            resampled_sst = self._update_field(filenames, 'sea_water_temperature', resampled_sst_data)
            if self._anomalies:
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
            resampled_sst_uncert = self._update_field(filenames, 'sea_water_temperature standard_error',
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
            sf = self._update_field(filenames, 'sea_ice_area_fraction', self._sf_data)
            sf.nc_set_variable('sea_fraction')
            sf.standard_name = 'sea_area_fraction'
            sf.long_name = 'sea area fraction'
            sf.comment = 'Fraction of cell that is covered by sea (including sea ice)'
            del sf.source

            # Create the resampled sif field
            if self._sea_ice_fraction:
                resampled_sif = self._update_field(filenames, 'sea_ice_area_fraction', resampled_sif_data)
                resampled_sif.long_name = 'regridded sea ice area fraction'
                resampled_sif.comment = 'Regridded sea ice area fraction'

            # Get the date of the resampled SST
            dt = resampled_sst.dimension_coordinate('T').datetime_array[0]

            # Create a field list with all the fields in it
            fl = cf.FieldList()
            fl.append(resampled_sst)
            fl.append(resampled_sst_uncert)
            if self._sea_ice_fraction:
                fl.append(resampled_sif)
            fl.append(sf)
            fl.append(SSTRegridder._create_time_field(sf, 'calendar_year', 'calendar year', dt.year))
            fl.append(SSTRegridder._create_time_field(sf, 'calendar_month', 'calendar month', dt.month))
            fl.append(SSTRegridder._create_time_field(sf, 'day_of_month', 'day of month', dt.day))
            fl.append(SSTRegridder._create_time_field(sf, 'day_of_year', 'day of year', dt.dayofyr))

            # Write the data
            output_path = os.path.join(self._out_path, dt.strftime('%Y%m%d') + '_regridded_sst.nc')
            cf.write(fl, output_path, datatype={np.dtype('float64'): np.dtype('float32')}, compress=1,
                     least_significant_digit=3)
            output_paths.append(output_path)

        # if specified, zip up all output files and remove the originals
        if self._zip_name:
            zip_path = os.path.join(self._out_path, self._zip_name)
            with zipfile.ZipFile(zip_path, "w") as z:
                for output_path in output_paths:
                    z.write(output_path, os.path.split(output_path)[1])
                    os.unlink(output_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Regrid SST data for a particular year.')
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
    parser.add_argument('--sst_cci_analysis_l4_path', default=regridding_constants._default_sst_cci_analysis_l4_path,
                        help='Path to the SST CCI Analysis Level 4 input data.')
    parser.add_argument('--c3s_sst_analysis_l4_path', default=regridding_constants._default_c3s_sst_analysis_l4_path,
                        help='Path to the C3S SST Analysis Level 4 input data.')
    parser.add_argument('--sst_cci_climatology_path', default=regridding_constants._default_sst_cci_climatology_path,
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

    sst_regridder = SSTRegridder(sst_cci_analysis_l4_path=args.sst_cci_analysis_l4_path,
                                 c3s_sst_analysis_l4_path=args.c3s_sst_analysis_l4_path,
                                 sst_cci_climatology_path=args.sst_cci_climatology_path)

    sst_regridder(args.lon_resolution, args.lat_resolution, time_resolution, year=args.year,
                  start_month=args.start_month, start_day=args.start_day, end_month=args.end_month,
                  end_day=args.end_day, anomalies=args.anomalies, no_sea_ice_fraction=args.no_sea_ice_fraction,
                  f_max=args.f_max, tau=args.tau, spatial_lambda=args.spatial_lambda, out_path=args.out_path,
                  zip_name=args.zip_name)
