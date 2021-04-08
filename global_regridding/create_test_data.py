#!/usr/bin/env python3

import cf
import numpy as np

# Number of latitudes and longitudes
lat_size = 72
lon_size = 144

# Maximum and minimum latitudes and longitudes
min_lat = 0.0
max_lat = 180.0
min_lon = 0.0
max_lon = 360.0

# Length of time (days)
time_size = 1

# Length of day (seconds)
day_length = 24 * 60 * 60


def create_test_data(lon_step=2, lat_step=2):
    # Create artificial SST data with the same value in each block to be averaged
    array = np.empty((time_size, lat_size, lon_size))
    for j in range(lat_size):
        for i in range(lon_size):
            array[:, j, i] = i // lon_step + (j // lat_step) * (lon_size // lon_step)
    data = cf.Data(array, units='kelvin')

    # Create time dimension
    array = np.arange(12 * 60 * 60, time_size * day_length, day_length, dtype='f8')
    bounds_array = np.column_stack((np.arange(0, time_size * day_length, day_length, dtype='f8'),
                                    np.arange(day_length, (time_size + 1) * day_length, day_length, dtype='f8')))
    bounds = cf.Bounds(data=cf.Data(bounds_array, units='seconds since 2017-01-01'))
    bounds.nc_set_variable('time_bnds')
    time = cf.DimensionCoordinate(data=cf.Data(array, units='seconds since 2017-01-01'), bounds=bounds)
    time.standard_name = 'time'
    time.nc_set_variable('time')

    # Create latitude dimension
    step = (max_lat - min_lat) / lat_size
    array = np.arange(min_lat + step / 2.0, max_lat, step)
    bounds_array = np.column_stack((np.arange(min_lat, max_lat, step), np.arange(min_lat + step, max_lat + step, step)))
    bounds = cf.Bounds(data=cf.Data(bounds_array, units='degrees_north'))
    bounds.nc_set_variable('lat_bnds')
    lat = cf.DimensionCoordinate(data=cf.Data(array, units='degrees_north'), bounds=bounds)
    lat.standard_name = 'latitude'
    lat.nc_set_variable('lat')

    # Create longitude dimension
    step = (max_lon - min_lon) / lon_size
    array = np.arange(min_lon + step / 2.0, max_lon, step)
    bounds_array = np.column_stack((np.arange(min_lon, max_lon, step), np.arange(min_lon + step, max_lon + step, step)))
    bounds = cf.Bounds(data=cf.Data(bounds_array, units='degrees_east'))
    bounds.nc_set_variable('lon_bnds')
    lon = cf.DimensionCoordinate(data=cf.Data(array, units='degrees_east'), bounds=bounds)
    lon.standard_name = 'longitude'
    lon.nc_set_variable('lon')

    # Create field
    f = cf.Field()
    f.set_construct(cf.DomainAxis(size=time_size), key='domainaxis0')
    f.set_construct(cf.DomainAxis(size=lat_size), key='domainaxis1')
    f.set_construct(cf.DomainAxis(size=lon_size), key='domainaxis2')
    f.set_data(data, axes=('domainaxis0', 'domainaxis1', 'domainaxis2'))
    f.standard_name = 'sea_water_temperature'
    f.nc_set_variable('analysed_sst')
    f.set_construct(time, axes=('domainaxis0',), key='dimensioncoordinate0', copy=False)
    f.set_construct(lat, axes=('domainaxis1',), key='dimensioncoordinate1', copy=False)
    f.set_construct(lon, axes=('domainaxis2',), key='dimensioncoordinate2', copy=False)

    # Write field to disk
    cf.write(f, 'test_data_' + str(lon_step) + '_' + str(lat_step) + '.nc')


if __name__ == '__main__':
    create_test_data(2, 2)
    create_test_data(3, 3)
