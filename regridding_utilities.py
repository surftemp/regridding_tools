"""
Common functionality between the regridding tools.
"""

import cf
from calendar import monthrange, isleap
from datetime import date, timedelta
from enum import Enum, auto
import glob
import numpy as np
import os

import regridding_constants


class InputType(Enum):
    SST_L4 = auto()
    SST_L3U = auto()


def allowed_lonlat_resolutions(base_resolution, max_resolution, lonlat_range, min_resolution=None):
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

    :param min_resolution:
        The minimum resolution for the longitude or latitude in degrees. Default is base_resolution (None).

    :return:
        A list of the allowed longitude or latitude steps for averaging in degrees.
    """
    allowed_resolutions = {}
    size = int(round(lonlat_range / base_resolution))
    for step in range(1, size):
        resolution = round(step * base_resolution, 2)
        if resolution > max_resolution:
            break
        if min_resolution is not None and resolution < min_resolution:
            continue
        if size % step == 0:
            allowed_resolutions[resolution] = step
    return allowed_resolutions


def create_lonlat_bounds(f):
    """
    Create longitude and latitude bounds for a field if necessary.

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


def create_time_field(f, variable_name, long_name, data_value):
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


class Regridder(object):
    """
    Base class for regridding high resolution (SST) data to more managable resolutions.
    """
    def resample_lonlat(self, f):
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

        resampled_lon_bounds_array = np.column_stack((lon_bounds_array[0::self.lon_step, 0],
                                                      lon_bounds_array[self.lon_step - 1::self.lon_step, 1]))
        resampled_lat_bounds_array = np.column_stack((lat_bounds_array[0::self.lat_step, 0],
                                                      lat_bounds_array[self.lat_step - 1::self.lat_step, 1]))

        # If the step is even use the value of the bounds of the old cell in the centre of the new cell as the new
        # coordinate values, otherwise use the coordinate values in the centre of the old cell.
        if self.lon_step % 2 == 0:
            resampled_lon_array = lon_bounds_array[self.lon_step // 2::self.lon_step, 0]
        else:
            resampled_lon_array = lon.array[self.lon_step // 2::self.lon_step]
        if self.lat_step % 2 == 0:
            resampled_lat_array = lat_bounds_array[self.lat_step // 2::self.lat_step, 0]
        else:
            resampled_lat_array = lat.array[self.lat_step // 2::self.lat_step]

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

        self.lon = resampled_lon
        self.lat = resampled_lat

    def spatially_resample_data(self, f, weights=False):
        """
        Resample a field's data spatially.

        :param f:
            The field whose data is to be resampled.

        :return:
            The resampled data as a cf.Data object.
        """
        if weights:
            array = (f * self.weights).array
        else:
            array = f.array

        array = array.reshape((1, self.lat.size, self.lat_step, self.lon.size, self.lon_step))
        array = array.transpose((0, 1, 3, 2, 4))
        array = array.reshape((1, self.lat.size, self.lon.size, self.lat_step * self.lon_step))
        array = np.sum(array, axis=3)

        data = cf.Data(array, units=f.units, fill_value=f.fill_value(default='netCDF'), dtype=f.dtype)
        return data

    def get_filename(self, d, climatology):
        """
        Generate the name of the input file for a particular date.

        :param d:
            Python data object.

        :param climatology:
            Boolean indicating whether to return the filename for the climatology.

        :return:
            The filename as a string.
        """
        if climatology:
            if d.month == 2 and d.day == 29:
                # If the day is the 29th of February the files for D59 and D60 will be averaged.
                filename = (os.path.join(self.sst_cci_climatology_path,
                                         'D059-ESACCI-L4_GHRSST-SSTdepth-OSTIA-GLOB_CDR2.1-v02.0-fv01.0.nc'),
                            os.path.join(self.sst_cci_climatology_path,
                                         'D060-ESACCI-L4_GHRSST-SSTdepth-OSTIA-GLOB_CDR2.1-v02.0-fv01.0.nc'))
            else:
                # Get day of year as if the year was a non-leap year.
                day_of_year = date(1982, d.month, d.day).timetuple().tm_yday
                filename = (os.path.join(self.sst_cci_climatology_path, 'D' + str(day_of_year).zfill(3) +
                                         '-ESACCI-L4_GHRSST-SSTdepth-OSTIA-GLOB_CDR2.1-v02.0-fv01.0.nc'))
        elif self.input_type is InputType.SST_L4:
            if d.year >= regridding_constants.c3s_start_year:
                filename = os.path.join(self.c3s_sst_analysis_l4_path,
                                        d.strftime('%Y'), d.strftime('%m'), d.strftime('%d'), d.strftime('%Y%m%d') +
                                        '120000-C3S-L4_GHRSST-SSTdepth-OSTIA-GLOB_ICDR2.0-v02.0-fv01.0.nc')
            else:
                filename = os.path.join(self.sst_cci_analysis_l4_path,
                                        d.strftime('%Y'), d.strftime('%m'), d.strftime('%d'), d.strftime('%Y%m%d') +
                                        '120000-ESACCI-L4_GHRSST-SSTdepth-OSTIA-GLOB_CDR2.1-v02.0-fv01.0.nc')
        elif self.input_type is InputType.SST_L3U:
            filename = sorted(glob.glob(os.path.join(self.sst_l3u_path, self.sensor, d.strftime('%Y'), d.strftime('%m'),
                                                     d.strftime('%d'), '*.nc')))
        else:
            raise ValueError('Input type not recognized.')
        return filename

    def get_start_end_dates(self, month):
        """
        Get the start dates and end dates for a particular month of the year.

        :param month:
            The month of the year.

        :return:
            The start date and end date as Python date objects.
        """
        days_in_month = monthrange(self.year, month)[1]
        if month == self.start_month == self.end_month:
            start_date = date(self.year, month, self.start_day)
            end_date = date(self.year, month, self.end_day)
        elif month == self.start_month:
            start_date = date(self.year, month, self.start_day)
            end_date = date(self.year, month, days_in_month)
        elif month == self.end_month:
            start_date = date(self.year, month, 1)
            end_date = date(self.year, month, self.end_day)
        else:
            start_date = date(self.year, month, 1)
            end_date = date(self.year, month, days_in_month)
        return start_date, end_date + timedelta(1)

    def create_filename_groups(self, climatology=False):
        """
        Create a list of lists of filenames. Each group of filenames will be averaged together over time. It is assumed
        that the dates and time resolution have already been checked for compatibility.

        :param climatology:
            Boolean indicating whether to generate the filename groups for the climatology.

        :return:
            A list of lists of filenames.
        """
        filename_groups = []
        if self.time_resolution == 'annual':
            filenames = []
            start_date = date(self.year, self.start_month, self.start_day)
            end_date = date(self.year, self.end_month, self.end_day) + timedelta(1)
            for d in daterange(start_date, end_date):
                filename = self.get_filename(d, climatology)
                filenames.append(filename)
            filename_groups.append(filenames)
        elif self.time_resolution == 'monthly':
            for month in range(self.start_month, self.end_month + 1):
                start_date = date(self.year, month, 1)
                days_in_month = monthrange(self.year, month)[1]
                end_date = date(self.year, month, days_in_month) + timedelta(1)
                filenames = []
                for d in daterange(start_date, end_date):
                    filename = self.get_filename(d, climatology)
                    filenames.append(filename)
                filename_groups.append(filenames)
        elif self.time_resolution == '10-day' or self.time_resolution == '5-day':
            step = 10 if self.time_resolution == '10-day' else 5
            break_day = 21 if self.time_resolution == '10-day' else 26
            for month in range(self.start_month, self.end_month + 1):
                days_in_month = monthrange(self.year, month)[1]
                start_date, end_date = self.get_start_end_dates(month)
                for outer_d in daterange(start_date, end_date, step=step):
                    if outer_d.day < break_day:
                        inner_drange = daterange(outer_d, outer_d + timedelta(step))
                    else:
                        inner_drange = daterange(outer_d, date(self.year, month, days_in_month) + timedelta(1))
                    filenames = []
                    for inner_d in inner_drange:
                        filename = self.get_filename(inner_d, climatology)
                        filenames.append(filename)
                    filename_groups.append(filenames)
                    if outer_d.day >= break_day:
                        break
        else:
            start_date = date(self.year, self.start_month, self.start_day)
            end_date = date(self.year, self.end_month, self.end_day) + timedelta(1)
            days_in_year = 366 if isleap(self.year) else 365
            finish = False
            for i in range(0, days_in_year, self.time_resolution):
                if i + self.time_resolution >= days_in_year:
                    drange = daterange(start_date + timedelta(i), end_date)
                elif i + 3 * self.time_resolution / 2 > days_in_year:
                    drange = daterange(start_date + timedelta(i), end_date)
                    finish = True
                else:
                    drange = daterange(start_date + timedelta(i), start_date + timedelta(i + self.time_resolution))
                filenames = []
                for d in drange:
                    filename = self.get_filename(d, climatology)
                    filenames.append(filename)
                filename_groups.append(filenames)
                if finish:
                    break
        return filename_groups

    def update_field(self, filenames, standard_name, resampled_data):
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
        g.domain_axes[lon_axis].set_size(self.lon.size)
        g.domain_axes[lat_axis].set_size(self.lat.size)
        g.set_construct(self.lon, axes=(lon_key,))
        g.set_construct(self.lat, axes=(lat_key,))

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
        if self.time_resolution == 'annual':
            dtarray[0] = cf.dt(dt.year, 7, 1, calendar='gregorian')
        elif self.time_resolution == 'monthly':
            dtarray[0] = cf.dt(dt.year, dt.month, 15, 12, calendar='gregorian')
        elif self.time_resolution == '10-day':
            day = dt.day
            if day <= 10:
                new_day = 5
            elif day <= 20:
                new_day = 15
            else:
                new_day = 25
            dtarray[0] = cf.dt(dt.year, dt.month, new_day, 12, calendar='gregorian')
        elif self.time_resolution == '5-day':
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
        if self.time_resolution == 'annual':
            time_coverage_duration = 'P1Y'
        elif self.time_resolution == 'monthly':
            time_coverage_duration = 'P1M'
        else:
            time_coverage_duration = 'P' + str(len(filenames)) + 'D'
        g.set_property('time_coverage_duration', time_coverage_duration)
        g.set_property('time_coverage_resolution', time_coverage_duration)

        # Change global attributes
        g.set_property('geospatial_lon_resolution', np.float32(self.lon_resolution))
        g.set_property('geospatial_lat_resolution', np.float32(self.lat_resolution))
        if self.lon_resolution == self.lat_resolution:
            g.set_property('spatial_resolution', str(self.lon_resolution) + ' degree')
        else:
            g.set_property('spatial_resolution', str(self.lon_resolution) + ' x ' + str(self.lat_resolution) +
                           ' degree')
        g.set_property('creator_processing_institution', 'These data were produced by the University of ' +
                       'Reading as part of the ESA CCI project.')
        if self.input_type is InputType.SST_L4:
            g.set_property('acknowledgment', 'Regridding service funded by National Centre for Earth ' +
                           'Observation, UK')

        # Tweak metadata as discussed in https://github.com/surftemp/sst-services/issues/17
        if self.input_type is InputType.SST_L4:
            g.set_property('summary', 'Regridding of L4 analysis product from the ESA SST CCI project')
        g.set_property('creator_email','c.j.merchant@reading.ac.uk')
        g.set_property('publisher_email', 'c.j.merchant@reading.ac.uk')
        if self.input_type is InputType.SST_L4:
            g.set_property('references','Merchant, C.J., Embury, O., Bulgin, C.E., Block, T., Corlett, G.K., Fiedler, E., Good, S.A., Mittaz, J., Rayner, N.A., Berry, D., Eastwood, S., Taylor, M., Tsushima, Y., Waterfall, A., Wilson, R. and Donlon, C. (2019), Satellite-based time-series of sea-surface temperature since 1981 for climate applications. Scientific Data 6, 223, doi:10.1038/s41597-019-0236-x, and Merchant, C. J. and Embury, O. (2020) Adjusting for desert-dust-related biases in a climate data record of sea surface temperature. Remote Sensing, 12 (16). 2554. ISSN 2072-4292 doi:10.3390/rs12162554')
        g.set_property('creator_url','http://climate.esa.int/')
        g.set_property('metadata_link','http://climate.esa.int/')
        if self.input_type is InputType.SST_L4:
            g.set_property('title','ESA SST CCI L4 analysis')
            g.set_property('history','Created by makegriddedSSTs.py v1.0 operating on SST CCI L4 analysis products')
        elif self.input_type is InputType.SST_L3U:
            g.set_property('history', 'Created by makegriddedL3USSTs.py v1.0 operating on SST CCI L3U products')

        g.del_property('product_specification_version')
        if self.input_type is InputType.SST_L4:
            g.nc_set_global_attribute('f_max', self.f_max)
            g.nc_set_global_attribute('lambda', self.spatial_lambda)
            g.nc_set_global_attribute('tau', np.int32(self.tau))
        g.set_property('date_created', self.date_created)
        g.set_property('uuid', self.uuid)
        g.set_property('tracking_id', self.uuid)

        # Close the fields
        fl.close()
        gl.close()

        return g
