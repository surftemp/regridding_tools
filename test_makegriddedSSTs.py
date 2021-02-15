#!/usr/bin/env python3

import unittest

import cf
from datetime import datetime
import numpy as np
import os

from create_test_data import create_test_data

from makegriddedSSTs import SSTRegridder
import regridding_utilities

class MakeGriddedSSTsTestCase(unittest.TestCase):
    def test_calculate_k_xy(self):
        def check_k_xy():
            # Create resampled longitude and latitude coordinates
            sst_regridder.resample_lonlat(f)

            # Calculate k_xy
            sst_regridder.calculate_k_xy()

            # Calculate the expected value of k_xy at each point
            cos_theta = np.cos(np.radians(sst_regridder.lat.array)).reshape(-1, 1)
            k_x_denominator = sst_regridder.spatial_lambda / cos_theta
            k_x = np.maximum(lon_resolution / k_x_denominator, 1.0)
            k_y = np.maximum(lat_resolution / sst_regridder.spatial_lambda, 1.0)
            expected_k_xy = k_x * k_y

            # Check that the k_xy data is as expected
            self.assertTrue((sst_regridder._k_xy == expected_k_xy).all(),
                            'Data of k_xy not as expected for a lon_step of ' + str(sst_regridder.lon_step) +
                            ' and a lat_step of ' + str(sst_regridder.lat_step))

            # Check the shape of k_xy
            self.assertEqual(sst_regridder._k_xy.shape, (1, sst_regridder.lat.size, sst_regridder.lon.size),
                             'Shape of k_xy not as expected for a lon_step of ' + str(sst_regridder.lon_step) +
                             ' and a lat_step of ' + str(sst_regridder.lat_step))

        # Initalise an SST regridder
        sst_regridder = SSTRegridder()

        # Read in a test field
        f = cf.read('test_data_2_2.nc')[0]

        # The spatial scale within which the observations are assumed to be fully correlated in degrees
        sst_regridder.spatial_lambda = 3.0

        # The integer steps with which to resample the longitude and latitude
        sst_regridder.lon_step = 2
        sst_regridder.lat_step = 2

        # The size of each cell in degrees
        lon_resolution = 5.0
        lat_resolution = 5.0

        # Check the lon and the lat are as expected
        check_k_xy()

        # The integer steps with which to resample the longitude and latitude
        sst_regridder.lon_step = 3
        sst_regridder.lat_step = 3

        # The size of each cell in degrees
        lon_resolution = 7.5
        lat_resolution = 7.5

        # Check the lon and the lat are as expected
        check_k_xy()

        # The integer steps with which to resample the longitude and latitude
        sst_regridder.lon_step = 3
        sst_regridder.lat_step = 2

        # The size of each cell in degrees
        lon_resolution = 7.5
        lat_resolution = 5.0

        # Check the lon and the lat are as expected
        check_k_xy()

    def test_create_filename_groups(self):
        def check_filename_groups():
            filename_groups = sst_regridder.create_filename_groups()
            filenames = [filename for filename_group in filename_groups for filename in filename_group]
            previous_date = None
            for filename in filenames:
                filename = os.path.split(filename)[1]
                date = datetime.strptime(filename[:8], '%Y%m%d')
                if previous_date is not None:
                    self.assertEqual((date - previous_date).days, 1, 'Filenames not contiguous for time resolution: ' +
                                     str(sst_regridder.time_resolution))
                previous_date = date
            self.assertEqual(len(filenames), n_days, 'Unexpected number of filenames for time resolution: ' +
                             str(sst_regridder.time_resolution))
            self.assertEqual(len(filename_groups), n_groups,
                             'Unexpected number of filename groups for time resolution: ' +
                             str(sst_regridder.time_resolution))

        # Initalise an SST regridder
        sst_regridder = SSTRegridder()

        sst_regridder.year = 1982
        sst_regridder.start_month = 1
        sst_regridder.start_day = 1
        sst_regridder.end_month = 12
        sst_regridder.end_day = 31
        n_days = 365

        sst_regridder.time_resolution = 'annual'
        n_groups = 1
        check_filename_groups()

        sst_regridder.time_resolution = 'monthly'
        n_groups = 12
        check_filename_groups()

        sst_regridder.time_resolution = '10-day'
        n_groups = 36
        check_filename_groups()

        sst_regridder.time_resolution = '5-day'
        n_groups = 72
        check_filename_groups()

        sst_regridder.time_resolution = 10
        n_groups = 37
        check_filename_groups()

        sst_regridder.time_resolution = 30
        n_groups = 12
        check_filename_groups()


if __name__ == '__main__':
    # Create test data files with the specified longitude and latitude steps
    create_test_data(lon_step=2, lat_step=2)
    create_test_data(lon_step=3, lat_step=3)

    unittest.main()
