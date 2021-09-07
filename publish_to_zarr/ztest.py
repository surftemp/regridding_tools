# -*- coding: utf-8 -*-

#    sst-services
#    Copyright (C) 2020  National Centre for Earth Observation (NCEO)
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""
Publish a netcdf4 SST dataset to zarr format, transporting all required fields and metadata
"""

import os.path
from logging import Logger, INFO, DEBUG, StreamHandler, Formatter
import zarr
import s3fs
import sys
import numpy
import datetime

from numcodecs import Zstd
compressor=Zstd(level=3)

# for S3 see https://docs.aws.amazon.com/general/latest/gr/aws-sec-cred-types.html
# access token and secret key need to be in a file ~/.aws/credentials, like:
#
# [default]
# aws_access_key_id=<access-key>
# aws_secret_access_key=<access-secret>

import xarray as xr
import numpy as np

DEFAULT_SST_CCI_ANALYSIS_L4_PATH = "/gws/nopw/j04/esacci_sst/users/niallmcc/CDR_v2/Analysis/L4/v2.1"
DEFAULT_C3S_SST_ANALYSIS_L4_PATH = "/gws/nopw/j04/esacci_sst/users/niallmcc/ICDR_v2/Analysis/L4/v2.0/"
DEFAULT_SST_CCI_CLIMATOLOGY_PATH = "/gws/nopw/j04/esacci_sst/users/niallmcc/CDR_v2/Climatology/L4/v2.1/"

DEFAULT_ZARR_PATH = "s3://surftemp-sst/sst.zarr"
# DEFAULT_ZARR_PATH="/gws/nopw/j04/esacci_sst/users/niallmcc/sst.zarr" # "/group_workspaces/jasmin2/nceo_uor/niall/reslice"

cci_sst_field_names = ['analysed_sst', 'analysed_sst_uncertainty', 'sea_ice_fraction', 'mask']
c3s_sst_field_names = ['analysed_sst', 'analysis_uncertainty', 'sea_ice_fraction', 'mask']
climatology_field_names = ['analysed_sst', 'sea_ice_fraction', 'mask']

c3s_rename = { 'analysis_uncertainty':'analysed_sst_uncertainty' }

logger = Logger("publisher")
logger.setLevel(INFO)
sh = StreamHandler(sys.stdout)
formatter = Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
sh.setFormatter(formatter)
logger.addHandler(sh)

class Tester(object):

    @staticmethod
    def parseTime(s):
        return datetime.datetime.strptime(s, "%Y%m%dT%H%M%SZ")

    def __init__(self,input_cci_folder,input_c3s_folder,climatology_folder):
        self.input_cci_folder = input_cci_folder
        self.input_c3s_folder = input_c3s_folder
        self.climatology_folder = climatology_folder

    def open_store(self,path):
        if path.startswith("s3:"):
            s3 = s3fs.S3FileSystem(anon=True)
            store = s3fs.S3Map(root=path, s3=s3, create=False)
        else:
            store = path
        return store

    def test_climatology(self,path):

        store = self.open_store(path)
        zarr_ds = xr.open_zarr(store, mode="r")

        for day in range(0,365):
            in_path = os.path.join(self.climatology_folder,
                                   "D%03d-ESACCI-L4_GHRSST-SSTdepth-OSTIA-GLOB_CDR2.1-v02.0-fv01.0.nc" % (
                                           1 + day))
            original_ds = xr.open_dataset(in_path)

            test_label = "test climatology day %d" % (day)
            for field_name in climatology_field_names:
                original_arr = original_ds[field_name].data.reshape((3600,7200))
                zarr_arr = zarr_ds[field_name].isel(time=day).data
                numpy.testing.assert_almost_equal(original_arr,zarr_arr)

            logger.info("Test %s completed", test_label)


    def test(self,path,start_dt,end_dt):

        store = self.open_store(path)

        zarr_ds = xr.open_zarr(store,mode="r")

        dt = start_dt
        while dt <= end_dt:
            year = dt.year
            field_names = cci_sst_field_names

            test_label = str(dt)
            if year < 2017:
                in_path = os.path.join(self.input_cci_folder,str(year),"%02d"%(dt.month),"%02d"%(dt.day),
                        "%s-ESACCI-L4_GHRSST-SSTdepth-OSTIA-GLOB_CDR2.1-v02.0-fv01.0.nc"%(dt.strftime("%Y%m%d120000")))
            else:
                in_path = os.path.join(self.input_c3s_folder,str(year),"%02d"%(dt.month),"%02d"%(dt.day),
                                           "%s-C3S-L4_GHRSST-SSTdepth-OSTIA-GLOB_ICDR2.0-v02.0-fv01.0.nc"%(dt.strftime("%Y%m%d120000")))

            original_ds = xr.open_dataset(in_path)

            for field_name in field_names:
                if year >= 2017 and field_name in c3s_rename:
                    zarr_field_name = c3s_rename[field_name]
                else:
                    zarr_field_name = field_name
                original_arr = original_ds[field_name].data.reshape((3600,7200))
                zarr_arr = zarr_ds[zarr_field_name].sel(time=dt).data
                numpy.testing.assert_almost_equal(original_arr,zarr_arr)

            dt = dt + datetime.timedelta(days=1)

            original_ds.close()

            logger.info("Test %s completed",test_label)

        logger.info("Completed comparison testing")

    @staticmethod
    def getDaysInYear(year):
        d1 = datetime.date(year, 1, 1)
        d2 = datetime.date(year, 12, 31)
        return 1 + (d2 - d1).days

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('--input-cci', action='store',
                        dest='input_cci',
                        type=str,
                        help='Specify the location of input CCI files',
                        default=DEFAULT_SST_CCI_ANALYSIS_L4_PATH)

    parser.add_argument('--input-c3s', action='store',
                        dest='input_c3s',
                        type=str,
                        help='Specify the location of input C3S files',
                        default=DEFAULT_C3S_SST_ANALYSIS_L4_PATH)

    parser.add_argument('--input-climatology', action='store',
                        dest='climatology_folder',
                        type=str,
                        help='Specify the location of input climatology files',
                        default=DEFAULT_SST_CCI_CLIMATOLOGY_PATH)

    parser.add_argument('--zarr_path',
                        dest='zarr_path',
                        help='path/URL to the zarr file',
                        default=DEFAULT_ZARR_PATH)

    parser.add_argument('--climatology', action='store_true',
                        dest='climatology',
                        help='Run the climatology')

    parser.add_argument('--start-date', action='store',
                        dest='start_date',
                        type=str,
                        help='Specify the start date for testing SST data, format YYYY-MM-DD',
                        metavar="YEAR",
                        default="2016-01-01")

    parser.add_argument('--end-date', action='store',
                        dest='end_date',
                        type=str,
                        help='Specify the end date for testing SST data, format YYYY-MM-DD',
                        metavar="YEAR",
                        default="2016-01-31")

    parser.add_argument('--verbose', action='store_true',
                        dest='verbose',
                        help='log debug output')

    args = parser.parse_args()
    if args.verbose:
        logger.setLevel(DEBUG)

    tester = Tester(args.input_cci, args.input_c3s, args.climatology_folder)

    if not args.climatology and (not args.start_date or not args.end_date):
        logger.warning("Please specify EITHER (--start-year <YEAR> AND --end-year <YEAR>) OR --climatology.")
        sys.exit(-1)

    start_dt = None if args.climatology else (datetime.datetime.strptime(args.start_date,"%Y-%m-%d") + datetime.timedelta(hours=12))
    end_dt = None if args.climatology else (datetime.datetime.strptime(args.end_date,"%Y-%m-%d") + datetime.timedelta(hours=12))

    tester.test(args.zarr_path,start_dt,end_dt)


