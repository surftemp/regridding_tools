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
import datetime
from logging import Logger, INFO, DEBUG, StreamHandler, Formatter
import zarr
import s3fs
import sys
import numpy

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

DEFAULT_SST_CCI_ANALYSIS_L4_PATH = "/Users/cv922550/Data/sst/data/CDR_v2/Analysis/L4/v2.1"  # "/neodc/esacci/sst/data/CDR_v2/Analysis/L4/v2.1"
DEFAULT_C3S_SST_ANALYSIS_L4_PATH = "/Users/cv922550/Data/sst/data/ICDR_v2/Analysis/L4/v2.0/"  # "/neodc/c3s_sst/data/ICDR_v2/Analysis/L4/v2.0/"
DEFAULT_SST_CCI_CLIMATOLOGY_PATH = "/Users/cv922550/Data/sst/data/CDR_v2/Climatology/L4/v2.1/"  # "/neodc/esacci/sst/data/CDR_v2/Climatology/L4/v2.1/"

DEFAULT_ZARR_PATH="/Users/cv922550/Data/2016.zarr" # "/group_workspaces/jasmin2/nceo_uor/niall/reslice"

cci_sst_field_names = ['analysed_sst', 'analysed_sst_uncertainty', 'sea_ice_fraction', 'mask']
c3s_sst_field_names = ['analysed_sst', 'analysis_uncertainty', 'sea_ice_fraction', 'mask']
climatology_field_names = ['analysed_sst', 'sea_ice_fraction', 'mask']

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

    def test(self,path,start_year,end_year):

        climatology = True if start_year is None else False

        if climatology:
            start_year = 0
            end_year = 0

        if path.startswith("s3:"):
            s3 = s3fs.S3FileSystem(anon=False)
            store = s3fs.S3Map(root=path, s3=s3, create=False)
        else:
            store = path

        zarr_ds = xr.open_zarr(store)

        for year in range(start_year,end_year+1):
            if climatology:
                field_names = climatology_field_names
            else:
                if year < 2017:
                    field_names = cci_sst_field_names
                else:
                    field_names = c3s_sst_field_names

            days = 365 if climatology else self.getDaysInYear(year)

            logger.info("Starting comparison testing for year %s"%(str(year)))

            day = 0

            days_processed = 0

            if year == 1981:
                day = 243 # only data from september 1st

            while day < days:

                if climatology:
                    in_path = os.path.join(self.climatology_folder,"D%03d-ESACCI-L4_GHRSST-SSTdepth-OSTIA-GLOB_CDR2.1-v02.0-fv01.0.nc"%(1+day+chunk_day))
                    test_label = "D%003d"%(day)
                else:
                    dt = datetime.datetime(year,1,1,12,0,0) + datetime.timedelta(days=day)
                    test_label = str(dt)
                    if year < 2017:
                        in_path = os.path.join(self.input_cci_folder,str(year),"%02d"%(dt.month),"%02d"%(dt.day),
                            "%s-ESACCI-L4_GHRSST-SSTdepth-OSTIA-GLOB_CDR2.1-v02.0-fv01.0.nc"%(dt.strftime("%Y%m%d120000")))
                    else:
                        in_path = os.path.join(self.input_c3s_folder,str(year),"%02d"%(dt.month),"%02d"%(dt.day),
                                               "%s-C3S-L4_GHRSST-SSTdepth-OSTIA-GLOB_ICDR2.0-v02.0-fv01.0.nc"%(dt.strftime("%Y%m%d120000")))

                original_ds = xr.open_dataset(in_path)

                for field_name in field_names:
                    original_arr = original_ds[field_name].data.reshape((3600,7200))
                    zarr_arr = zarr_ds[field_name].sel(time=dt).data
                    numpy.testing.assert_almost_equal(original_arr,zarr_arr)
                day += 1

                original_ds.close()

                logger.debug("Test %s completed",test_label)

            logger.info("Completed comparison testing for year %s to %s (processed %d days)"%(str(year),path,days_processed))

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

    parser.add_argument('--start-year', action='store',
                        dest='start_year',
                        type=int,
                        help='Specify the start year for reslicing SST data',
                        metavar="YEAR",
                        default=2016)

    parser.add_argument('--end-year', action='store',
                        dest='end_year',
                        type=int,
                        help='Specify the end year for reslicing SST data',
                        metavar="YEAR",
                        default=2016)

    parser.add_argument('--verbose', action='store_true',
                        dest='verbose',
                        help='log debug output')

    args = parser.parse_args()
    if args.verbose:
        logger.setLevel(DEBUG)

    tester = Tester(args.input_cci, args.input_c3s, args.climatology_folder)

    if not args.climatology and (not args.start_year or not args.end_year):
        logger.warning("Please specify EITHER (--start-year <YEAR> AND --end-year <YEAR>) OR --climatology.")
        sys.exit(-1)

    start_year = None if args.climatology else args.start_year
    end_year = None if args.climatology else args.end_year

    tester.test(args.zarr_path,start_year,end_year)


