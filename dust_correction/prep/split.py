# split yearly netcdf4 files into daily ones
#
# usage:
#
# python split.py <input-folder> <output-folder>
#
# for years 2003 to 2018, load the file regridded<year>_0p5.nc and output files yyyy/mm/dd/cams_ecmwf_yyyymmdd1200_dust.nc
#

import cf
import os
import sys
import datetime
import numpy as np

input_folder = sys.argv[1]
output_folder = sys.argv[2]

for year in range(2003,2020):

    dt = datetime.datetime(year,1,1)

    input_path = os.path.join(input_folder,"regridded%s_0p5.nc"% (year))
    f = cf.read(input_path)

    os.makedirs(output_folder,exist_ok=True)
    year_folder = os.path.join(output_folder,str(year))
    os.makedirs(year_folder,exist_ok=True)
    days = f[0].data.shape[0]
    print(days)

    for day in range(0,days):
        day_dt = dt+datetime.timedelta(days=day)
        month = "%02d"%(day_dt.month)
        month_day = "%02d"%(day_dt.day)
        month_folder = os.path.join(year_folder,month)
        os.makedirs(month_folder,exist_ok=True)
        day_folder = os.path.join(month_folder,month_day)
        os.makedirs(day_folder,exist_ok=True)

        filename = "cams_ecmwf_%d%s%s1200_dust.nc"%(year,month,month_day)
        outpath = os.path.join(day_folder,filename)
        dayfields = []
        for field in range(4):
            dayfields.append(f[field].subspace[day,...])
        for dayfield in dayfields:
            dayfield.data.nc_set_hdf5_chunksizes([1, 360, 720])
        cf.write(dayfields,outpath,datatype={np.dtype('float64'): np.dtype('float32')}, compress=1)
        for dayfield in dayfields:
            dayfield.close()
        print(outpath)

    f.close()