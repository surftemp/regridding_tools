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
"""

import os.path
import datetime
import sys
from logging import Logger, INFO, DEBUG, StreamHandler, Formatter
import zarr
import math

from numcodecs import Zstd
compressor=Zstd(level=1)

import cf
import numpy as np

DEFAULT_SST_ANALYSIS_L4_PATH = "/neodc/c3s_sst/data/ICDR_v2/Analysis/L4/v2.0/"
DEFAULT_SST_CCI_CLIMATOLOGY_PATH = "/neodc/esacci/sst/data/CDR_v2/Climatology/L4/v2.1/"

DEFAULT_OUTPUT_FOLDER="/gws/nopw/j04/nceo_uor/niall/reslice"

INPUT_RESOLUTION=0.05

sst_field_names = ['sea_water_temperature', 'sea_water_temperature standard_error', 'sea_ice_area_fraction']
climatology_field_names = ['sea_water_temperature']

logger = Logger("reslicer")
logger.setLevel(INFO)
sh = StreamHandler(sys.stdout)
formatter = Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
sh.setFormatter(formatter)
logger.addHandler(sh)

class Reslicer(object):

    @staticmethod
    def parseTime(s):
        return datetime.datetime.strptime(s, "%Y%m%dT%H%M%SZ")

    def __init__(self,input_sst_folder,climatology_folder,output_resolution):
        self.input_sst_folder = input_sst_folder
        self.climatology_folder = climatology_folder
        self.output_resolution = output_resolution

    def reslice(self,output_parent_folder,year):

        days = 365 if year == "climatology" else self.getDaysInYear(year)

        chunk_size_days = 7

        if year == "climatology":
            field_names = climatology_field_names
        else:
            field_names = sst_field_names

        path = os.path.join(output_parent_folder,str(year)+".zarr")
        spatial_chunk_size = round(self.output_resolution/INPUT_RESOLUTION)
        za = zarr.open(path, mode='w',shape=(days,3600,7200,len(field_names)),chunks=(chunk_size_days, spatial_chunk_size, spatial_chunk_size, len(field_names)), compressor=compressor, dtype='f4')

        logger.info("Starting reslicing for year %s"%(str(year)))

        day = 0
        self.resetProgress()

        premature_end_of_data = False
        days_processed = 0

        if year == 1981:
            day = 243 # only data from september 1st

        while day < days:
            chunk_day = 0
            chunk_data = []
            while (chunk_day < chunk_size_days) and (not premature_end_of_data):
                if day + chunk_day >= days:
                    break
                if year == "climatology":
                    dt = datetime.datetime(2018, 1, 1) + datetime.timedelta(day+chunk_day)
                else:
                    dt = datetime.datetime(year, 1, 1) + datetime.timedelta(day+chunk_day)

                self.reportProgress(day,days,"reslice_days IN")

                if year == "climatology":
                    in_path = os.path.join(self.climatology_folder,"D%03d-ESACCI-L4_GHRSST-SSTdepth-OSTIA-GLOB_CDR3.0-v02.0-fv01.0.nc"%(1+day+chunk_day))
                else:
                    in_path = os.path.join(self.input_sst_folder,str(year),"%02d"%(dt.month),"%02d"%(dt.day),
                        "%s-C3S-L4_GHRSST-SSTdepth-OSTIA-GLOB_ICDR3.0-v02.0-fv01.0.nc"%(dt.strftime("%Y%m%d120000")))

                if not os.path.exists(in_path):
                    msg = "could not locate input file %s"%(in_path)
                    logger.warning(msg)
                    # NOTE - not a problem if the data does not exist for the year
                    premature_end_of_data = True
                    break

                logger.debug("reading %s"%(in_path))
                in_fields = cf.read(in_path,extra='field_ancillary')

                arrs = []
                for field_idx in range(len(field_names)):
                    field_name = field_names[field_idx]
                    a = in_fields.select(field_name)[0].array
                    a = np.ma.filled(a,math.nan)
                    a = a.reshape(1, 3600, 7200, 1)
                    arrs.append(a)

                chunk_data.append(np.concatenate(arrs,axis=3))
                chunk_day += 1
                days_processed += 1
                in_fields.close()

            if chunk_data:
                chunk = np.concatenate(chunk_data, axis=0)
                za[day:day+chunk.shape[0], :, :, :] = chunk
                day += chunk_size_days

            if premature_end_of_data:
                break

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

    parser.add_argument('--input-sst', action='store',
                        dest='input_sst',
                        type=str,
                        help='Specify the location of input SST files',
                        default=DEFAULT_SST_ANALYSIS_L4_PATH)

    parser.add_argument('--input-climatology', action='store',
                        dest='climatology_folder',
                        type=str,
                        help='Specify the location of input climatology files',
                        default=DEFAULT_SST_CCI_CLIMATOLOGY_PATH)

    parser.add_argument('--year', action='store',
                        dest='year',
                        type=str,
                        help='Specify the year for reslicing SST data, or climatology',
                        metavar="YEAR",
                        default="2018")

    parser.add_argument('--output-resolution', action='store',
                        dest='output_resolution',
                        type=int,
                        help='Specify the spatial resolution of chunks, in degrees',
                        default=5)

    parser.add_argument('--output-folder', action='store',
                        dest='output_folder',
                        type=str,
                        help='Specify the location of output files',
                        default="/tmp")

    parser.add_argument('--verbose', action='store_true',
                        dest='verbose',
                        help='log debug output')

    args = parser.parse_args()
    if args.verbose:
        logger.setLevel(DEBUG)

    reslicer = Reslicer(args.input_sst, args.climatology_folder, args.output_resolution)

    if not args.year:
        logger.warning("Please specify --year <YEAR> or --year climatology.")
        sys.exit(-1)

    year = args.year
    if year != "climatology":
        year = int(year)


    reslicer.reslice(args.output_folder,year)


