#!/usr/bin/env python3

import unittest

import cf
import numpy as np

from create_test_data import create_test_data

from makegriddedSSTs import SSTRegridder
import regridding_functions

class RegriddingFunctionsTestCase(unittest.TestCase):
    def test_allowed_lonlat_resolutions(self):
        # Set parameters
        base_resolution = 0.05
        max_resolution = 5.0
        lon_range = 360.0
        lat_range = 180.0
        n_lon_resolutions = 30
        n_lat_resolutions = 28

        # Check that the allowed longitude resolutions are as expected
        allowed_lon_resolutions = regridding_functions.allowed_lonlat_resolutions(base_resolution, max_resolution, lon_range)
        for lon_resolution in allowed_lon_resolutions:
            self.assertAlmostEqual(round(lon_resolution / base_resolution) * base_resolution, lon_resolution, 10,
                                   'lon_resolution is not a multiple of ' + str(base_resolution) + '.')
            self.assertAlmostEqual(round(lon_range / lon_resolution) * lon_resolution, lon_range, 10,
                                   'lon_resolution is not a factor of ' + str(lon_range) + '.')
        self.assertTrue(len(allowed_lon_resolutions) == n_lon_resolutions,
                        'Expected ' + str(n_lon_resolutions) + ' allowed longitude resolutions.')

        # Check that the allowed latitude resolutions are as expected
        allowed_lat_resolutions = regridding_functions.allowed_lonlat_resolutions(base_resolution, max_resolution, lat_range)
        for lat_resolution in allowed_lat_resolutions:
            self.assertAlmostEqual(round(lat_resolution / base_resolution) * base_resolution, lat_resolution, 10,
                                   'lat_resolution is not a multiple of ' + str(base_resolution) + '.')
            self.assertAlmostEqual(round(lat_range / lat_resolution) * lat_resolution, lat_range, 10,
                                   'lat_resolution is not a factor of ' + str(lat_range) + '.')
        self.assertTrue(len(allowed_lat_resolutions) == n_lat_resolutions,
                        'Expected ' + str(n_lat_resolutions) + ' allowed latitude resolutions.')

    def test_resample_lonlat(self):
        def check_lonlat():
            # Create resampled longitude and latitude coordinates
            regridding_functions.resample_lonlat(sst_regridder, f)

            # Check that the coordinates' data is as expected
            self.assertTrue((sst_regridder.lon.array == expected_lon).all(),
                            'Longitude coordinates not as expected for lon_step of ' + str(sst_regridder.lon_step))
            self.assertTrue((sst_regridder.lat.array == expected_lat).all(),
                            'Latitude coordinates not as expected for lat_step of ' + str(sst_regridder.lat_step))

            # Check the new coordinates' bounds are contiguous and nonoverlapping
            self.assertTrue(sst_regridder.lon.contiguous(overlap=False),
                            'Bounds not contiguous and nonoverlapping for lon_step of ' + str(sst_regridder.lon_step))
            self.assertTrue(sst_regridder.lat.contiguous(overlap=False),
                            'Bounds not contiguous and nonoverlapping for lat_step of ' + str(sst_regridder.lat_step))

        # Initalise an SST regridder
        sst_regridder = SSTRegridder()

        # Read in a test field
        f = cf.read('test_data_2_2.nc')[0]

        # The integer steps with which to resample the longitude and latitude
        sst_regridder.lon_step = 2
        sst_regridder.lat_step = 2

        # The expected lon and lat coordinates
        expected_lon = np.array(np.arange(2.5, 360.0, 5.0))
        expected_lat = np.array(np.arange(2.5, 180.0, 5.0))

        # Check the lon and the lat are as expected
        check_lonlat()

        # The integer steps with which to resample the longitude and latitude
        sst_regridder.lon_step = 3
        sst_regridder.lat_step = 3

        # The expected lon and lat coordinates
        expected_lon = np.array(np.arange(3.75, 360.0, 7.5))
        expected_lat = np.array(np.arange(3.75, 180.0, 7.5))

        # Check the lon and the lat are as expected
        check_lonlat()

    def test_spatially_resample_data(self):
        def check_spatially_resample_data():
            # Read in the test field
            f = cf.read('test_data_' + str(sst_regridder.lon_step) + '_' + str(sst_regridder.lat_step) + '.nc')[0]

            # Create resampled longitude and latitude coordinates
            regridding_functions.resample_lonlat(sst_regridder, f)

            # Resample the field's data
            g = regridding_functions.spatially_resample_data(sst_regridder, f[0])
            g /= regridding_functions.spatially_resample_data(sst_regridder, f[0].where(True, 1.0))

            # Check the data is as expected
            array = g.array
            for j in range(sst_regridder.lat.size):
                for i in range(sst_regridder.lon.size):
                    self.assertTrue((array[:, j, i] == i + j * sst_regridder.lon.size).all(),
                                    'Data of field is not as expected for lon_step = ' + str(sst_regridder.lon_step) +
                                    ' and lat_step = ' + str(sst_regridder.lat_step) + '.')

        # Initalise an SST regridder
        sst_regridder = SSTRegridder()

        # The integer steps with which to resample the longitude and latitude
        sst_regridder.lon_step = 2
        sst_regridder.lat_step = 2

        # Check spatial resampling
        check_spatially_resample_data()

        # The integer steps with which to resample the longitude and latitude
        sst_regridder.lon_step = 3
        sst_regridder.lat_step = 3

        # Check spatial resampling
        check_spatially_resample_data()


if __name__ == '__main__':
    # Create test data files with the specified longitude and latitude steps
    create_test_data(lon_step=2, lat_step=2)
    create_test_data(lon_step=3, lat_step=3)

    unittest.main()
