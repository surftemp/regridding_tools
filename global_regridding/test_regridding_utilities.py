#!/usr/bin/env python3

import unittest

import cf
from datetime import datetime
import numpy as np
import os

from create_test_data import create_test_data

import regridding_constants
import regridding_utilities

class RegriddingUtilitiesTestCase(unittest.TestCase):
    def test_allowed_lonlat_resolutions(self):
        # Set parameters
        base_resolution = 0.05
        max_resolution = 5.0
        lon_range = 360.0
        lat_range = 180.0
        n_lon_resolutions = 30
        n_lat_resolutions = 28

        # Check that the allowed longitude resolutions are as expected
        allowed_lon_resolutions = regridding_utilities.allowed_lonlat_resolutions(base_resolution, max_resolution, lon_range)
        for lon_resolution in allowed_lon_resolutions:
            self.assertAlmostEqual(round(lon_resolution / base_resolution) * base_resolution, lon_resolution, 10,
                                   'lon_resolution is not a multiple of ' + str(base_resolution) + '.')
            self.assertAlmostEqual(round(lon_range / lon_resolution) * lon_resolution, lon_range, 10,
                                   'lon_resolution is not a factor of ' + str(lon_range) + '.')
        self.assertTrue(len(allowed_lon_resolutions) == n_lon_resolutions,
                        'Expected ' + str(n_lon_resolutions) + ' allowed longitude resolutions.')

        # Check that the allowed latitude resolutions are as expected
        allowed_lat_resolutions = regridding_utilities.allowed_lonlat_resolutions(base_resolution, max_resolution, lat_range)
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
            regridder.resample_lonlat(f)

            # Check that the coordinates' data is as expected
            self.assertTrue((regridder.lon.array == expected_lon).all(),
                            'Longitude coordinates not as expected for lon_step of ' + str(regridder.lon_step))
            self.assertTrue((regridder.lat.array == expected_lat).all(),
                            'Latitude coordinates not as expected for lat_step of ' + str(regridder.lat_step))

            # Check the new coordinates' bounds are contiguous and nonoverlapping
            self.assertTrue(regridder.lon.contiguous(overlap=False),
                            'Bounds not contiguous and nonoverlapping for lon_step of ' + str(regridder.lon_step))
            self.assertTrue(regridder.lat.contiguous(overlap=False),
                            'Bounds not contiguous and nonoverlapping for lat_step of ' + str(regridder.lat_step))

        # Initalise an SST regridder
        regridder = regridding_utilities.Regridder()

        # Read in a test field
        f = cf.read('test_data_2_2.nc')[0]

        # The integer steps with which to resample the longitude and latitude
        regridder.lon_step = 2
        regridder.lat_step = 2

        # The expected lon and lat coordinates
        expected_lon = np.array(np.arange(2.5, 360.0, 5.0))
        expected_lat = np.array(np.arange(2.5, 180.0, 5.0))

        # Check the lon and the lat are as expected
        check_lonlat()

        # The integer steps with which to resample the longitude and latitude
        regridder.lon_step = 3
        regridder.lat_step = 3

        # The expected lon and lat coordinates
        expected_lon = np.array(np.arange(3.75, 360.0, 7.5))
        expected_lat = np.array(np.arange(3.75, 180.0, 7.5))

        # Check the lon and the lat are as expected
        check_lonlat()

    def test_spatially_resample_data(self):
        def check_spatially_resample_data():
            # Read in the test field
            f = cf.read('test_data_' + str(regridder.lon_step) + '_' + str(regridder.lat_step) + '.nc')[0]

            # Create resampled longitude and latitude coordinates
            regridder.resample_lonlat(f)

            # Resample the field's data
            g = regridder.spatially_resample_data(f[0])
            g /= regridder.spatially_resample_data(f[0].where(True, 1.0))

            # Check the data is as expected
            array = g.array
            for j in range(regridder.lat.size):
                for i in range(regridder.lon.size):
                    self.assertTrue((array[:, j, i] == i + j * regridder.lon.size).all(),
                                    'Data of field is not as expected for lon_step = ' + str(regridder.lon_step) +
                                    ' and lat_step = ' + str(regridder.lat_step) + '.')

        # Initalise an SST regridder
        regridder = regridding_utilities.Regridder()

        # The integer steps with which to resample the longitude and latitude
        regridder.lon_step = 2
        regridder.lat_step = 2

        # Check spatial resampling
        check_spatially_resample_data()

        # The integer steps with which to resample the longitude and latitude
        regridder.lon_step = 3
        regridder.lat_step = 3

        # Check spatial resampling
        check_spatially_resample_data()

    def test_create_filename_groups(self):
        def check_filename_groups():
            filename_groups = regridder.create_filename_groups(regridding_utilities.InputType.L4)
            filenames = [filename for filename_group in filename_groups for filename in filename_group]
            previous_date = None
            for filename in filenames:
                filename = os.path.split(filename)[1]
                date = datetime.strptime(filename[:8], '%Y%m%d')
                if previous_date is not None:
                    self.assertEqual((date - previous_date).days, 1, 'Filenames not contiguous for time resolution: ' +
                                     str(regridder.time_resolution))
                previous_date = date
            self.assertEqual(len(filenames), n_days, 'Unexpected number of filenames for time resolution: ' +
                             str(regridder.time_resolution))
            self.assertEqual(len(filename_groups), n_groups,
                             'Unexpected number of filename groups for time resolution: ' +
                             str(regridder.time_resolution))

        # Initalise an SST regridder
        regridder = regridding_utilities.Regridder()

        regridder.sst_cci_analysis_l4_path = regridding_constants.default_sst_cci_analysis_l4_path
        regridder.year = 1982
        regridder.start_month = 1
        regridder.start_day = 1
        regridder.end_month = 12
        regridder.end_day = 31
        n_days = 365

        regridder.time_resolution = 'annual'
        n_groups = 1
        check_filename_groups()

        regridder.time_resolution = 'monthly'
        n_groups = 12
        check_filename_groups()

        regridder.time_resolution = '10-day'
        n_groups = 36
        check_filename_groups()

        regridder.time_resolution = '5-day'
        n_groups = 72
        check_filename_groups()

        regridder.time_resolution = 10
        n_groups = 37
        check_filename_groups()

        regridder.time_resolution = 30
        n_groups = 12
        check_filename_groups()


if __name__ == '__main__':
    # Create test data files with the specified longitude and latitude steps
    create_test_data(lon_step=2, lat_step=2)
    create_test_data(lon_step=3, lat_step=3)

    unittest.main()
