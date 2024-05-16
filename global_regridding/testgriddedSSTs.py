#!/usr/bin/env python3

import argparse

import cf
import numpy as np

import glob

# check SST data files for the mask mismatch issue reported in https://github.com/surftemp/common/issues/4
# makegriddedSSTs.py assumes that the masks are identical for the SST, sea ice fraction and SST uncertainty fields
# if any SST files fail the test, makegriddedSSTs.py may produce incorrect results

class SSTTester:

    def __init__(self):
        pass

    def test(self, input_path):

        filenames = glob.glob(input_path, recursive=True)
        count = 0
        total = len(filenames)
        for filename in filenames:

            fl = cf.read(filename, aggregate=False,extra='field_ancillary')

            sst = fl.select_by_property(standard_name='sea_water_temperature')[0]

            sif = fl.select_by_property(standard_name='sea_ice_area_fraction')[0]

            sst_uncert = fl.select_by_property(standard_name='sea_water_temperature standard_error')[0]

            # check that the SST mask is the same as the SIF mask, print the filename if not
            if not np.array_equal(np.ma.getmask(sst.array), np.ma.getmask(sif.array), True):
                print(filename)

            # check that the SST mask is the same as the uncertainty mask, print the filename if not
            if not np.array_equal(np.ma.getmask(sst.array), np.ma.getmask(sst_uncert.array), True):
                print(filename)

            fl.close()

            count += 1
            if count % 20 == 0:
                print("Processed: "+str(count)+" Fraction complete: "+str(count/total))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Check SST files for mismatched mask issue (https://github.com/surftemp/common/issues/4).')

    parser.add_argument('--path', help='path with wildcards to netcdf4 files to check')

    args = parser.parse_args()

    sstt = SSTTester()

    sstt.test(args.path)


