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
import sys
from logging import Logger, INFO, DEBUG, StreamHandler, Formatter
import zarr
import s3fs
import sys

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

DEFAULT_SST_CCI_ANALYSIS_L4_PATH = "/neodc/esacci/sst/data/CDR_v2/Analysis/L4/v2.1"
DEFAULT_C3S_SST_ANALYSIS_L4_PATH = "/neodc/c3s_sst/data/ICDR_v2/Analysis/L4/v2.0/"
DEFAULT_SST_CCI_CLIMATOLOGY_PATH = "/neodc/esacci/sst/data/CDR_v2/Climatology/L4/v2.1/"

DEFAULT_OUTPUT_PATH="s3://surftemp-sst/sst.zarr"

cci_sst_field_names = ['analysed_sst', 'analysed_sst_uncertainty', 'sea_ice_fraction', 'mask']
c3s_sst_field_names = ['analysed_sst', 'analysis_uncertainty', 'sea_ice_fraction', 'mask']
climatology_field_names = ['analysed_sst', 'sea_ice_fraction', 'mask']

logger = Logger("publisher")
logger.setLevel(INFO)
sh = StreamHandler(sys.stdout)
formatter = Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
sh.setFormatter(formatter)
logger.addHandler(sh)

class Publisher(object):

    @staticmethod
    def parseTime(s):
        return datetime.datetime.strptime(s, "%Y%m%dT%H%M%SZ")

    def __init__(self,input_cci_folder,input_c3s_folder,climatology_folder,chunk_size):
        self.input_cci_folder = input_cci_folder
        self.input_c3s_folder = input_c3s_folder
        self.climatology_folder = climatology_folder
        self.chunk_size = chunk_size

    def publish(self,path,start_year,end_year,first_days_only=None):

        climatology = True if start_year is None else False

        creating = True

        if climatology:
            start_year = 0
            end_year = 0

        if path.startswith("s3:"):
            s3 = s3fs.S3FileSystem(anon=False)
            store = s3fs.S3Map(root=path, s3=s3, create=False)
        else:
            store = path

        if not climatology:
            try:
                ds = zarr.open(store=store,mode="r")
                print(ds.info)
                print(ds.tree())
                creating = False
                for key in ds.attrs:
                    print(key,ds.attrs[key])
                # get start and end times encoded as YYYYMMDDT000000Z
                start_time = ds.attrs["start_time"]
                end_time = ds.attrs["stop_time"]
                current_end_year = int(end_time[0:4])-1 # the stored value should be Jan 1 00:00:00 of the year following the last added data
                logger.info("Appending to existing zarr store covering (%s - %s)"%(start_time,end_time))
                if current_end_year+1 != start_year:
                    logger.error("Error, please append from year %d"%(current_end_year+1))
                    sys.exit(-1)

            except zarr.errors.PathNotFoundError:
                start_time = "%04d0101T000000Z"%(start_year)
                logger.info("No existing zarr store to append to")

        for year in range(start_year,end_year+1):
            if climatology:
                field_names = climatology_field_names
            else:
                if year < 2017:
                    field_names = cci_sst_field_names
                else:
                    field_names = c3s_sst_field_names

            if first_days_only is not None:
                days = first_days_only
            else:
                days = 365 if climatology else self.getDaysInYear(year)

            chunk_size_time = self.chunk_size[0]
            chunk_size_lat = self.chunk_size[1]
            chunk_size_lon = self.chunk_size[2]

            logger.info("Starting publication for year %s"%(str(year)))

            day = 0

            premature_end_of_data = False
            days_processed = 0

            if year == 1981:
                day = 243 # only data from september 1st

            while day < days:
                chunk_day = 0
                chunk_paths = []
                while (chunk_day < chunk_size_time) and (not premature_end_of_data):
                    if day + chunk_day >= days:
                        break
                    if climatology:
                        dt = datetime.datetime(2018, 1, 1) + datetime.timedelta(day+chunk_day)
                    else:
                        dt = datetime.datetime(year, 1, 1) + datetime.timedelta(day+chunk_day)

                    if climatology:
                        in_path = os.path.join(self.climatology_folder,"D%03d-ESACCI-L4_GHRSST-SSTdepth-OSTIA-GLOB_CDR2.1-v02.0-fv01.0.nc"%(1+day+chunk_day))
                    elif year < 2017:
                        in_path = os.path.join(self.input_cci_folder,str(year),"%02d"%(dt.month),"%02d"%(dt.day),
                            "%s-ESACCI-L4_GHRSST-SSTdepth-OSTIA-GLOB_CDR2.1-v02.0-fv01.0.nc"%(dt.strftime("%Y%m%d120000")))
                    else:
                        in_path = os.path.join(self.input_c3s_folder,str(year),"%02d"%(dt.month),"%02d"%(dt.day),
                                                   "%s-C3S-L4_GHRSST-SSTdepth-OSTIA-GLOB_ICDR2.0-v02.0-fv01.0.nc"%(dt.strftime("%Y%m%d120000")))

                    if not os.path.exists(in_path):
                        msg = "could not locate input file %s"%(in_path)
                        logger.warning(msg)
                        # NOTE - not a problem if the data does not exist for the year
                        premature_end_of_data = True
                        break

                    logger.debug("reading %s"%(in_path))
                    chunk_paths.append(in_path)
                    chunk_day += 1
                    days_processed += 1

                # AWS S3 path
                if path.startswith("s3:"):
                    s3 = s3fs.S3FileSystem(anon=False)
                    store = s3fs.S3Map(root=path, s3=s3, create=True)
                else:
                    store = path

                if len(chunk_paths):
                    if creating:
                        ds = xr.open_dataset(chunk_paths[0]).chunk(
                            {"time": chunk_size_time, "lat": chunk_size_lat, "lon": chunk_size_lon})
                        for name in list(ds.variables):
                            if name != "lat_bnds" and name != "lon_bnds":
                                del ds[name]
                        ds.to_zarr(store=store, mode="w")
                        ds.close()


                    ds = xr.open_mfdataset(chunk_paths, combine='by_coords').chunk(
                        {"time": chunk_size_time, "lat": chunk_size_lat, "lon": chunk_size_lon})

                    if "lat_bnds" in ds:
                        del ds["lat_bnds"]
                    if "lon_bnds" in ds:
                        del ds["lon_bnds"]

                    if creating:
                        ds.to_zarr(store=store,mode="a", encoding={field_name:{"compressor":compressor} for field_name in field_names})
                    else:
                        ds.to_zarr(store=store,mode="a", append_dim="time",
                                   encoding={})

                    ds.close()
                    logger.debug("Appended data to zarr store %s" % (path))
                    creating = False

                day += len(chunk_paths)

            partial_or_complete = "partial" if premature_end_of_data else "complete"
            logger.info("Completed %s publication for year %s to %s (processed %d days)"%(partial_or_complete,str(year),path,days_processed))

        if not climatology:
            try:
                ds = zarr.open(store=store,mode="a")
                ds.attrs["start_time"] = start_time
                ds.attrs["time_coverage_start"] = start_time
                stop_time = "%04d0101T000000Z"%(end_year+1)
                ds.attrs["stop_time"] = stop_time
                ds.attrs["time_coverage_end"] = stop_time
            except Exception as ex:
                print(ex)

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

    parser.add_argument('--climatology', action='store_true',
                        dest='climatology',
                        help='Run the climatology')

    parser.add_argument('--start-year', action='store',
                        dest='start_year',
                        type=int,
                        help='Specify the start year for reslicing SST data',
                        metavar="YEAR")

    parser.add_argument('--end-year', action='store',
                        dest='end_year',
                        type=int,
                        help='Specify the end year for reslicing SST data',
                        metavar="YEAR")

    parser.add_argument('--chunk-size', action='store',
                        dest='chunk_size',
                        type=str,
                        help='resolution of chunks, in format "days-size,lat-size,lon-size"',
                        default="5,360,720")

    parser.add_argument('--output-path', action='store',
                        dest='output_path',
                        type=str,
                        help='Specify the location of output zarr data, either local filesystem or s3://bucket/folder/filename.zarr',
                        default=DEFAULT_OUTPUT_PATH)

    parser.add_argument('--first-days-only', type=int,
                        help='for testing, include only this many days from the start of each year',
                        default=None)

    parser.add_argument('--verbose', action='store_true',
                        dest='verbose',
                        help='log debug output')

    args = parser.parse_args()
    if args.verbose:
        logger.setLevel(DEBUG)

    chunk_size = tuple(map(lambda n:int(n),args.chunk_size.split(",")))
    publisher = Publisher(args.input_cci, args.input_c3s, args.climatology_folder, chunk_size)

    if not args.climatology and (not args.start_year or not args.end_year):
        logger.warning("Please specify EITHER (--start-year <YEAR> AND --end-year <YEAR>) OR --climatology.")
        sys.exit(-1)

    start_year = None if args.climatology else args.start_year
    end_year = None if args.climatology else args.end_year

    publisher.publish(args.output_path,start_year,end_year,args.first_days_only)


