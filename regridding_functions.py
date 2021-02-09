"""
Common functionality between the regridding tool.
"""

import cf
from datetime import timedelta
import numpy as np


def allowed_lonlat_resolutions(base_resolution, max_resolution, lonlat_range):
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


def resample_lonlat(regridder, f):
    """
    Create resampled longitude and latitude cf dimension coordinates.

    :param regridder:
        regridder object with lon_step and lat_step.
    :param f:
        cf Field object with longitudes and latitudes.
    """
    # Resample the latitude and longitude arrays and their bounds
    lon = f.dimension_coordinate('X')
    lat = f.dimension_coordinate('Y')

    lon_bounds_array = lon.bounds.array
    lat_bounds_array = lat.bounds.array

    resampled_lon_bounds_array = np.column_stack((lon_bounds_array[0::regridder.lon_step, 0],
                                                  lon_bounds_array[regridder.lon_step - 1::regridder.lon_step, 1]))
    resampled_lat_bounds_array = np.column_stack((lat_bounds_array[0::regridder.lat_step, 0],
                                                  lat_bounds_array[regridder.lat_step - 1::regridder.lat_step, 1]))

    # If the step is even use the value of the bounds of the old cell in the centre of the new cell as the new
    # coordinate values, otherwise use the coordinate values in the centre of the old cell.
    if regridder.lon_step % 2 == 0:
        resampled_lon_array = lon_bounds_array[regridder.lon_step // 2::regridder.lon_step, 0]
    else:
        resampled_lon_array = lon.array[regridder.lon_step // 2::regridder.lon_step]
    if regridder.lat_step % 2 == 0:
        resampled_lat_array = lat_bounds_array[regridder.lat_step // 2::regridder.lat_step, 0]
    else:
        resampled_lat_array = lat.array[regridder.lat_step // 2::regridder.lat_step]

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

    regridder.lon = resampled_lon
    regridder.lat = resampled_lat


def spatially_resample_data(regridder, f, weights=False):
    """
    Resample a field's data spatially.

    :param regridder:
        regridder object with information for regridding.
    :param f:
        The field whose data is to be resampled.

    :return:
        The resampled data as a cf.Data object.
    """
    if weights:
        array = (f * regridder.weights).array
    else:
        array = f.array

    array = array.reshape((1, regridder.lat.size, regridder.lat_step, regridder.lon.size, regridder.lon_step))
    array = array.transpose((0, 1, 3, 2, 4))
    array = array.reshape((1, regridder.lat.size, regridder.lon.size, regridder.lat_step * regridder.lon_step))
    array = np.sum(array, axis=3)

    data = cf.Data(array, units=f.units, fill_value=f.fill_value(default='netCDF'), dtype=f.dtype)
    return data


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
