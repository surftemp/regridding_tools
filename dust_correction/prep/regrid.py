#!/usr/bin/env python
# -*- coding: utf-8 -*-

# wrapper script for regridding using cfpython
#
# usage:
#
# python regrid.py <input_file.nc> <output_file.nc> <reference.nc>
#
# performs spherical regridding using the bilinear method
#

import cf
import sys

input_path = sys.argv[1]
output_path = sys.argv[2]
reference_file = sys.argv[3]

g = cf.read(reference_file)[0]

print("Converting %s -> %s"%(input_path,output_path))
f = cf.read(input_path)
regridded = []
for field in range(4):
    rg = f[field].regrids(g,method="bilinear",src_cyclic=True,dst_cyclic=True)
    print(".")
    regridded.append(rg)

cf.write(regridded,output_path)