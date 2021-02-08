#!/usr/bin/env python3

import unittest

import cf
from datetime import datetime
import numpy as np
import os

from create_test_data import create_test_data

from makegriddedSSTs import SSTRegridder

class MakeGriddedSSTsTestCase(unittest.TestCase):
    def test_allowed_lonlat_resolutions(self):
        # Set parameters
        base_resolution = 0.05
        max_resolution = 5.0
        lon_range = 360.0
        lat_range = 180.0
        n_lon_resolutions = 30
        n_lat_resolutions = 28

        # Check that the allowed longitude resolutions are as expected
        allowed_lon_resolutions = SSTRegridder._allowed_lonlat_resolutions(base_resolution, max_resolution, lon_range)
        for lon_resolution in allowed_lon_resolutions:
            self.assertAlmostEqual(round(lon_resolution / base_resolution) * base_resolution, lon_resolution, 10,
                                   'lon_resolution is not a multiple of ' + str(base_resolution) + '.')
            self.assertAlmostEqual(round(lon_range / lon_resolution) * lon_resolution, lon_range, 10,
                                   'lon_resolution is not a factor of ' + str(lon_range) + '.')
        self.assertTrue(len(allowed_lon_resolutions) == n_lon_resolutions,
                        'Expected ' + str(n_lon_resolutions) + ' allowed longitude resolutions.')

        # Check that the allowed latitude resolutions are as expected
        allowed_lat_resolutions = SSTRegridder._allowed_lonlat_resolutions(base_resolution, max_resolution, lat_range)
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
            sst_regridder._resample_lonlat(f)

            # Check that the coordinates' data is as expected
            self.assertTrue((sst_regridder._lon.array == expected_lon).all(),
                            'Longitude coordinates not as expected for lon_step of ' + str(sst_regridder._lon_step))
            self.assertTrue((sst_regridder._lat.array == expected_lat).all(),
                            'Latitude coordinates not as expected for lat_step of ' + str(sst_regridder._lat_step))

            # Check the new coordinates' bounds are contiguous and nonoverlapping
            self.assertTrue(sst_regridder._lon.contiguous(overlap=False),
                            'Bounds not contiguous and nonoverlapping for lon_step of ' + str(sst_regridder._lon_step))
            self.assertTrue(sst_regridder._lat.contiguous(overlap=False),
                            'Bounds not contiguous and nonoverlapping for lat_step of ' + str(sst_regridder._lat_step))

        # Initalise an SST regridder
        sst_regridder = SSTRegridder()

        # Read in a test field
        f = cf.read('test_data_2_2.nc')[0]

        # The integer steps with which to resample the longitude and latitude
        sst_regridder._lon_step = 2
        sst_regridder._lat_step = 2

        # The expected lon and lat coordinates
        expected_lon = np.array(np.arange(2.5, 360.0, 5.0))
        expected_lat = np.array(np.arange(2.5, 180.0, 5.0))

        # Check the lon and the lat are as expected
        check_lonlat()

        # The integer steps with which to resample the longitude and latitude
        sst_regridder._lon_step = 3
        sst_regridder._lat_step = 3

        # The expected lon and lat coordinates
        expected_lon = np.array(np.arange(3.75, 360.0, 7.5))
        expected_lat = np.array(np.arange(3.75, 180.0, 7.5))

        # Check the lon and the lat are as expected
        check_lonlat()

    def test_calculate_k_xy(self):
        def check_k_xy():
            # Create resampled longitude and latitude coordinates
            sst_regridder._resample_lonlat(f)

            # Calculate k_xy
            sst_regridder._calculate_k_xy()

            # Calculate the expected value of k_xy at each point
            cos_theta = np.cos(np.radians(sst_regridder._lat.array)).reshape(-1, 1)
            k_x_denominator = sst_regridder._spatial_lambda / cos_theta
            k_x = np.maximum(lon_resolution / k_x_denominator, 1.0)
            k_y = np.maximum(lat_resolution / sst_regridder._spatial_lambda, 1.0)
            expected_k_xy = k_x * k_y

            # Check that the k_xy data is as expected
            self.assertTrue((sst_regridder._k_xy == expected_k_xy).all(),
                            'Data of k_xy not as expected for a lon_step of ' + str(sst_regridder._lon_step) +
                            ' and a lat_step of ' + str(sst_regridder._lat_step))

            # Check the shape of k_xy
            self.assertEqual(sst_regridder._k_xy.shape, (1, sst_regridder._lat.size, sst_regridder._lon.size),
                             'Shape of k_xy not as expected for a lon_step of ' + str(sst_regridder._lon_step) +
                             ' and a lat_step of ' + str(sst_regridder._lat_step))

        # Initalise an SST regridder
        sst_regridder = SSTRegridder()

        # Read in a test field
        f = cf.read('test_data_2_2.nc')[0]

        # The spatial scale within which the observations are assumed to be fully correlated in degrees
        sst_regridder._spatial_lambda = 3.0

        # The integer steps with which to resample the longitude and latitude
        sst_regridder._lon_step = 2
        sst_regridder._lat_step = 2

        # The size of each cell in degrees
        lon_resolution = 5.0
        lat_resolution = 5.0

        # Check the lon and the lat are as expected
        check_k_xy()

        # The integer steps with which to resample the longitude and latitude
        sst_regridder._lon_step = 3
        sst_regridder._lat_step = 3

        # The size of each cell in degrees
        lon_resolution = 7.5
        lat_resolution = 7.5

        # Check the lon and the lat are as expected
        check_k_xy()

        # The integer steps with which to resample the longitude and latitude
        sst_regridder._lon_step = 3
        sst_regridder._lat_step = 2

        # The size of each cell in degrees
        lon_resolution = 7.5
        lat_resolution = 5.0

        # Check the lon and the lat are as expected
        check_k_xy()

    def test_spatially_resample_data(self):
        def check_spatially_resample_data():
            # Read in the test field
            f = cf.read('test_data_' + str(sst_regridder._lon_step) + '_' + str(sst_regridder._lat_step) + '.nc')[0]

            # Create resampled longitude and latitude coordinates
            sst_regridder._resample_lonlat(f)

            # Resample the field's data
            g = sst_regridder._spatially_resample_data(f[0])
            g /= sst_regridder._spatially_resample_data(f[0].where(True, 1.0))

            # Check the data is as expected
            array = g.array
            for j in range(sst_regridder._lat.size):
                for i in range(sst_regridder._lon.size):
                    self.assertTrue((array[:, j, i] == i + j * sst_regridder._lon.size).all(),
                                    'Data of field is not as expected for lon_step = ' + str(sst_regridder._lon_step) +
                                    ' and lat_step = ' + str(sst_regridder._lat_step) + '.')

        # Initalise an SST regridder
        sst_regridder = SSTRegridder()

        # The integer steps with which to resample the longitude and latitude
        sst_regridder._lon_step = 2
        sst_regridder._lat_step = 2

        # Check spatial resampling
        check_spatially_resample_data()

        # The integer steps with which to resample the longitude and latitude
        sst_regridder._lon_step = 3
        sst_regridder._lat_step = 3

        # Check spatial resampling
        check_spatially_resample_data()

    def test_create_filename_groups(self):
        def check_filename_groups():
            filename_groups = sst_regridder._create_filename_groups()
            filenames = [filename for filename_group in filename_groups for filename in filename_group]
            previous_date = None
            for filename in filenames:
                filename = os.path.split(filename)[1]
                date = datetime.strptime(filename[:8], '%Y%m%d')
                if previous_date is not None:
                    self.assertEqual((date - previous_date).days, 1, 'Filenames not contiguous for time resolution: ' +
                                     str(sst_regridder._time_resolution))
                previous_date = date
            self.assertEqual(len(filenames), n_days, 'Unexpected number of filenames for time resolution: ' +
                             str(sst_regridder._time_resolution))
            self.assertEqual(len(filename_groups), n_groups,
                             'Unexpected number of filename groups for time resolution: ' +
                             str(sst_regridder._time_resolution))

        # Initalise an SST regridder
        sst_regridder = SSTRegridder()

        sst_regridder._year = 1982
        sst_regridder._start_month = 1
        sst_regridder._start_day = 1
        sst_regridder._end_month = 12
        sst_regridder._end_day = 31
        n_days = 365

        sst_regridder._time_resolution = 'annual'
        n_groups = 1
        check_filename_groups()

        sst_regridder._time_resolution = 'monthly'
        n_groups = 12
        check_filename_groups()

        sst_regridder._time_resolution = '10-day'
        n_groups = 36
        check_filename_groups()

        sst_regridder._time_resolution = '5-day'
        n_groups = 72
        check_filename_groups()

        sst_regridder._time_resolution = 10
        n_groups = 37
        check_filename_groups()

        sst_regridder._time_resolution = 30
        n_groups = 12
        check_filename_groups()


if __name__ == '__main__':
    # Create test data files with the specified longitude and latitude steps
    create_test_data(lon_step=2, lat_step=2)
    create_test_data(lon_step=3, lat_step=3)

    unittest.main()
