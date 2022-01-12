import xarray as xr
import numpy as np
import argparse
import os
import datetime
import sys
from logging import Logger, StreamHandler, Formatter, INFO
import netCDF4

_dustname = 'dustAdjustmentParametersAllYears_{:.1f}.nc'
_spikename = 'L4_spike_adjustment_{:.1f}.nc'

logger = Logger("dust_correction")
logger.setLevel(INFO)
sh = StreamHandler(sys.stdout)
formatter = Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
sh.setFormatter(formatter)
logger.addHandler(sh)

class DustAdjuster:
    """
    taken from https://github.com/surftemp/plots-dev/blob/master/projects/cci3/l4_adjust/adjust_md.py
    """

    def __init__(self, data_root, dust_power=0.5):
        self.data_dir = data_root
        self.dust_power = dust_power
        self.read_dust()
        self.read_spike()

    def read_dust(self):
        """Read the dust adjustment file and expand the array so we can use
        simple linear interpolation"""
        filename = os.path.join(self.data_dir, _dustname.format(self.dust_power))
        logger.info("Reading {}".format(filename))
        ds = xr.open_dataset(filename)
        # Indices for expansion. Lon wraps, lat uses nearest
        ilat = [0] + list(range(len(ds.lat))) + [-1]
        ilon = [-1] + list(range(len(ds.lon))) + [0]
        # Expand dust array
        ds = ds.isel(lat=ilat, lon=ilon)
        # Coordinates are Pandas IndexVariable which are immutable.
        # Must copy to a numpy array before we can modify them.
        new_lat = ds.lat.values
        new_lat[0] = -92.5
        new_lat[-1] = 92.5
        ds['lat'] = new_lat
        new_lon = ds.lon.values
        new_lon[0] -= 360
        new_lon[-1] += 360
        ds['lon'] = new_lon
        self.dust = ds

    def read_spike(self):
        """Read the spike adjustment file"""
        filename = os.path.join(self.data_dir, _spikename.format(self.dust_power))
        logger.info("Reading {}".format(filename))
        ds = xr.open_dataset(filename)
        # there are some NaNs at the end of the array (Dec 17-31 2018) - fill the values with zero
        self.spike = ds.fillna(0.0)

    def calculate(self, lon, lat, time):
        """Calcualte the L4 SST adjustment for given lat, lon, time"""
        iopts = {'fill_value': 0.0}  # Options passed to interpolator
        dsDi = self.dust.interp(lon=lon, lat=lat, time=time, kwargs=iopts)
        dsSi = self.spike.interp(time=time, kwargs=iopts)
        adjust = dsDi.scale * dsDi.dust**self.dust_power  # units of K
        uncert = dsDi.frUnc * adjust                      # units of K
        adjust += dsSi.global_l4_adj.values
        return adjust, uncert


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data-dir', type=str, help="directory containing dust data files", default="/Users/cv922550/Data/dust")
    parser.add_argument('--input-dir', type=str, help="input folder with SST files", default="/Users/cv922550/Data/sst/data/ICDR_v2/Analysis/L4/v2.0/")
    parser.add_argument('--input-year', type=int, help="year to process (defaults to all years)",default=-1)
    parser.add_argument('--input-month', type=int, help="input month to process (defaults to all months)", default=-1)
    parser.add_argument('--input-day', type=int, help="input day to process (defaults to all days)", default=-1)
    parser.add_argument('--output-dir', type=str, help="output folder to store adjusted SST files", default="/tmp")
    parser.add_argument('--cdr-version', type=str, help="set the CDR version (2.0, 2.1 etc)", default="2.1")

    args = parser.parse_args()

    # create an adjuster to perform the dust correction
    adjuster = DustAdjuster(args.data_dir,1.0)
    cdr_version = args.cdr_version

    # we expect the input directory to hold data files according to the standard directory structure YYYY/MM/DD

    # find all input files to correct by scanning the input directory, filtered according to any specified input year/month/day
    for year in sorted(os.listdir(args.input_dir)):
        year_dir = os.path.join(args.input_dir, year)
        if not os.path.isdir(year_dir):
            continue
        if args.input_year != -1 and int(year) != args.input_year:
            continue
        for month in sorted(os.listdir(year_dir)):
            month_dir = os.path.join(year_dir, month)
            if not os.path.isdir(month_dir):
                continue
            if args.input_month != -1 and int(month) != args.input_month:
                continue
            for day in sorted(os.listdir(month_dir)):
                day_dir = os.path.join(month_dir, day)
                if not os.path.isdir(day_dir):
                    continue
                if args.input_day != -1 and int(day) != args.input_day:
                    continue

                # build the name of the input file we expect, either C3S or CCI
                dt = datetime.datetime(year=int(year),month=int(month),day=int(day))
                if int(year) < 2017:
                    filename = "%s-ESACCI-L4_GHRSST-SSTdepth-OSTIA-GLOB_CDR%s-v02.0-fv01.0.nc"%(dt.strftime("%Y%m%d120000"),cdr_version)
                    sst_fieldname = 'analysed_sst'
                    sst_uncertainty_fieldname = 'analysed_sst_uncertainty'
                else:
                    filename = "%s-C3S-L4_GHRSST-SSTdepth-OSTIA-GLOB_ICDR%s-v02.0-fv01.0.nc"%(dt.strftime("%Y%m%d120000"),cdr_version)
                    sst_fieldname = 'analysed_sst'
                    if cdr_version == "2.0":
                        sst_uncertainty_fieldname = 'analysis_uncertainty'
                    else:
                        sst_uncertainty_fieldname = 'analysed_sst_uncertainty'
                in_path = os.path.join(day_dir,filename)

                # open the dataset
                logger.info("Reading %s"%(in_path))
                ds = xr.open_dataset(in_path)

                # apply the adjustment to the sst and uncertainty fields
                adj, unc = adjuster.calculate(ds.lon, ds.lat, ds.time)
                ds[sst_fieldname] += adj
                ds[sst_uncertainty_fieldname].values = np.sqrt(np.array(ds[sst_uncertainty_fieldname].values ** 2 + unc.values ** 2))

                # create a folder under the output directory of the same structure as the input, ie YYYY/MM/DD
                out_folder = os.path.join(args.output_dir,year,month,day)
                os.makedirs(out_folder,exist_ok=True)

                # write the adjusted dataset into the output folder using the same filename
                out_path = os.path.join(out_folder,filename)
                logger.info("Writing %s" % (out_path))

                if int(year) >= 2017:
                    # avoid strange issue with CF incompatible calendar (https://github.com/surftemp/sst-services/issues/10)
                    encoding = {"time":{"calendar":"gregorian","units":"seconds since 1981-01-01 00:00:00"}}
                    if "time_bnds" in ds:
                        encoding["time_bnds"] = {"calendar":"gregorian","units":"seconds since 1981-01-01 00:00:00"}
                    ds.to_netcdf(out_path,encoding=encoding)
                else:
                    ds.to_netcdf(out_path)


