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

def zdump(path):
    ds = zarr.open(path)
    print(ds.info)
    for key in ds.attrs:
        print("\t%s => %s"%(key,str(ds.attrs[key])))
    print(ds.tree())

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('--file',
                        type=str,
                        help='Specify the location of zarr file',default="/tmp/2018.zarr")

    args = parser.parse_args()

    zdump(args.file)


