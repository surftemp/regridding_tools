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

import xarray as xr
import zarr
import s3fs

import matplotlib.pyplot as plt

def zdump(path,plotvar,outpath):
    print(path)

    zarr_ds = xr.open_zarr(path)

    print(zarr_ds)

    if plotvar:
        plotvar = plotvar.split(":")
        name = plotvar[0]
        time = int(plotvar[1])
        lat_start = int(plotvar[2])
        lat_end = int(plotvar[3])
        lon_start = int(plotvar[4])
        lon_end = int(plotvar[5])
        zarr_ds[name].isel(time=time,lat=range(lat_start,lat_end),lon=range(lon_start,lon_end)).plot()
        plt.show()

    if outpath:
        from numcodecs import Zstd
        encoding = {}
        for name in zarr_ds.variables:
            encoding[name] = {"zlib": True, "complevel": 3}
        zarr_ds.to_netcdf(outpath,encoding=encoding)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('--path',
                        type=str,
                        help='Specify the location of zarr file',default="s3://surftemp.sst/zarr/2018.zarr")

    parser.add_argument('--plot',
                        type=str,
                        help='Specify the name of a variable to plot', default="analysed_sst:0:1000:1200:2000:2200")

    parser.add_argument('--output-path',
                        type=str,
                        help='Specify the location of a netcdf4 file to dump data to', default="")

    args = parser.parse_args()

    zdump(args.path,args.plot,args.output_path)


