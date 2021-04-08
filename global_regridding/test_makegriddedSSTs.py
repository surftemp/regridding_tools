#!/usr/bin/env python3

import unittest

import cf
import numpy as np

from create_test_data import create_test_data

from makegriddedSSTs import SSTRegridder

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


if __name__ == '__main__':
    # Create test data files with the specified longitude and latitude steps
    create_test_data(lon_step=2, lat_step=2)
    create_test_data(lon_step=3, lat_step=3)

    unittest.main()
