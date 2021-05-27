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
import math
from collections import defaultdict

from numcodecs import Zstd
compressor=Zstd(level=1)

import xarray as xr
import numpy as np

DEFAULT_SST_CCI_ANALYSIS_L4_PATH = "/Users/cv922550/Data/sst/data/CDR_v2/Analysis/L4/v2.1"  # "/neodc/esacci/sst/data/CDR_v2/Analysis/L4/v2.1"
DEFAULT_C3S_SST_ANALYSIS_L4_PATH = "/Users/cv922550/Data/sst/data/ICDR_v2/Analysis/L4/v2.0/"  # "/neodc/c3s_sst/data/ICDR_v2/Analysis/L4/v2.0/"
DEFAULT_SST_CCI_CLIMATOLOGY_PATH = "/Users/cv922550/Data/sst/data/CDR_v2/Climatology/L4/v2.1/"  # "/neodc/esacci/sst/data/CDR_v2/Climatology/L4/v2.1/"

DEFAULT_OUTPUT_FOLDER="/tmp" # "/group_workspaces/jasmin2/nceo_uor/niall/reslice"

INPUT_RESOLUTION=0.05

sst_field_names = ['analysed_sst', 'analysis_uncertainty', 'sea_ice_fraction']
climatology_field_names = ['sea_water_temperature']

logger = Logger("reslicer")
logger.setLevel(INFO)
sh = StreamHandler(sys.stdout)
formatter = Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
sh.setFormatter(formatter)
logger.addHandler(sh)

def make_json_compliant(v):
    if isinstance(v, np.int32) or isinstance(v, np.int16) or isinstance(v, np.int8):
        v = int(v)
    elif isinstance(v, np.float32):
        v = float(v)
    return v

class Publisher(object):

    @staticmethod
    def parseTime(s):
        return datetime.datetime.strptime(s, "%Y%m%dT%H%M%SZ")

    def __init__(self,input_cci_folder,input_c3s_folder,climatology_folder,chunk_size):
        self.input_cci_folder = input_cci_folder
        self.input_c3s_folder = input_c3s_folder
        self.climatology_folder = climatology_folder
        self.chunk_size = chunk_size

    def publish(self,output_parent_folder,year,first_days_only=None):

        if first_days_only is not None:
            days = first_days_only
        else:
            days = 365 if year == "climatology" else self.getDaysInYear(year)

        chunk_size_days = self.chunk_size[0]

        if year == "climatology":
            field_names = climatology_field_names
        else:
            field_names = sst_field_names

        path = os.path.join(output_parent_folder,str(year)+".zarr")

        rootgroup = zarr.open(path, mode='w')
        datasets = {}
        for field_name in field_names:
            datasets[field_name] = rootgroup.create_dataset(field_name,shape=(days,3600,7200),chunks=self.chunk_size, compressor=compressor, dtype='f4')

        logger.info("Starting reslicing for year %s"%(str(year)))

        day = 0
        self.resetProgress()

        premature_end_of_data = False
        days_processed = 0

        attrs = {}
        field_attrs = {}

        if year == 1981:
            day = 243 # only data from september 1st

        while day < days:
            chunk_day = 0
            chunk_data = defaultdict(lambda :[])
            while (chunk_day < chunk_size_days) and (not premature_end_of_data):
                if day + chunk_day >= days:
                    break
                if year == "climatology":
                    dt = datetime.datetime(2018, 1, 1) + datetime.timedelta(day+chunk_day)
                else:
                    dt = datetime.datetime(year, 1, 1) + datetime.timedelta(day+chunk_day)

                self.reportProgress(day,days,"reslice_days IN")

                if year == "climatology":
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
                ds = xr.open_dataset(in_path)
                print(ds)
                if len(attrs) == 0:
                    for key in ds.attrs:
                        attrs[key] = ds.attrs[key]
                else:
                    for key in ds.attrs:
                        if key in attrs:
                            if attrs[key] != ds.attrs[key]:
                                print(key,attrs[key])
                                del attrs[key]

                for field_name in field_names:
                    if field_name not in field_attrs:
                        field_attrs[field_name] = {}
                        for key in ds[field_name].attrs:
                            field_attrs[field_name][key] = ds[field_name].attrs[key]
                    else:
                        for key in ds[field_name].attrs:
                            if key in field_attrs[field_name]:
                                if field_attrs[field_name][key] != ds[field_name].attrs[key]:
                                    del field_attrs[field_name][key]
                    a = ds[field_name].data
                    a = np.ma.filled(a,math.nan)
                    a = a.reshape(1, 3600, 7200)
                    chunk_data[field_name].append(a)
                chunk_day += 1
                days_processed += 1
                ds.close()

            for field_name in field_names:
                chunk = np.concatenate(chunk_data[field_name], axis=0)
                datasets[field_name][day:day+chunk.shape[0], :, :] = chunk

            day += chunk_size_days

        for key in attrs:
            rootgroup.attrs[key] = make_json_compliant(attrs[key])

        for field_name in field_attrs:
            for key in field_attrs[field_name]:
                datasets[field_name].attrs[key] = make_json_compliant(field_attrs[field_name][key])

        partial_or_complete = "partial" if premature_end_of_data else "complete"
        logger.info("Completed %s reslicing for year %s to %s (processed %d days)"%(partial_or_complete,str(year),path,days_processed))


    def resetProgress(self):
        self.pct = None

    def reportProgress(self,n,total,label):
        pct = int((100 * n)/total)
        if self.pct is None or pct != self.pct:
            logger.info("%s progress %d percent" % (label,pct))
            self.pct = pct

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

    parser.add_argument('--year', action='store',
                        dest='year',
                        type=str,
                        help='Specify the year for reslicing SST data',
                        metavar="YEAR",
                        default="2018")

    parser.add_argument('--chunk_size', action='store',
                        dest='chunk_size',
                        type=str,
                        help='resolution of chunks, in format "days,degrees-lat,degrees-lon"',
                        default="7,360,720")

    parser.add_argument('--output-folder', action='store',
                        dest='output_folder',
                        type=str,
                        help='Specify the location of output files',
                        default=DEFAULT_OUTPUT_FOLDER)

    parser.add_argument('--first_days_only', type=int,
                        help='for testing, include only this many days from the start of the year',
                        default=None)

    parser.add_argument('--verbose', action='store_true',
                        dest='verbose',
                        help='log debug output')

    args = parser.parse_args()
    if args.verbose:
        logger.setLevel(DEBUG)

    chunk_size = tuple(map(lambda n:int(n),args.chunk_size.split(",")))
    publisher = Publisher(args.input_cci, args.input_c3s, args.climatology_folder,chunk_size)

    if not args.year:
        logger.warning("Please specify --year <YEAR> or --year climatology.")
        sys.exit(-1)

    year = args.year
    if year != "climatology":
        year = int(year)

    publisher.publish(args.output_folder,year,args.first_days_only)


