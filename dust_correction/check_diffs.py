import xarray

import xarray as xr
import numpy as np
import argparse
import os
import datetime
import sys
from logging import Logger, StreamHandler, Formatter, INFO

_dustname = 'dustAdjustmentParametersAllYears_{:.1f}.nc'
_spikename = 'L4_spike_adjustment_{:.1f}.nc'

logger = Logger("dust_correction")
logger.setLevel(INFO)
sh = StreamHandler(sys.stdout)
formatter = Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
sh.setFormatter(formatter)
logger.addHandler(sh)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-dir', type=str, help="input folder with SST files", default="/Users/cv922550/Data/sst/data/ICDR_v2/Analysis/L4/v2.0/")
    parser.add_argument('--input-year', type=int, help="year to process (defaults to all years)",default=2018)
    parser.add_argument('--input-month', type=int, help="input month to process (defaults to all months)", default=12)
    parser.add_argument('--input-day', type=int, help="input day to process (defaults to all days)", default=31)
    parser.add_argument('--output-dir', type=str, help="output folder to store adjusted SST files", default="/tmp")

    max_diff = None
    max_diff_date = None
    min_diff = None
    min_diff_date = None

    max_uncertainty_diff = None
    max_uncertainty_diff_date = None
    min_uncertainty_diff = None
    min_uncertainty_diff_date = None

    args = parser.parse_args()
    # we expect the input directory to hold data files according to the standard directory structure YYYY/MM/DD

    # find all input files to correct by scanning the input directory, filtered according to any specified input year/month/day
    for year in sorted(os.listdir(args.input_dir)):

        max_year_diff = None
        max_year_diff_date = None
        min_year_diff = None
        min_year_diff_date = None

        max_year_uncertainty_diff = None
        max_year_uncertainty_diff_date = None
        min_year_uncertainty_diff = None
        min_year_uncertainty_diff_date = None

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
                    filename = "%s-ESACCI-L4_GHRSST-SSTdepth-OSTIA-GLOB_CDR2.1-v02.0-fv01.0.nc"%(dt.strftime("%Y%m%d120000"))
                    sst_fieldname = 'analysed_sst'
                    sst_uncertainty_fieldname = 'analysed_sst_uncertainty'
                else:
                    filename = "%s-C3S-L4_GHRSST-SSTdepth-OSTIA-GLOB_ICDR2.0-v02.0-fv01.0.nc"%(dt.strftime("%Y%m%d120000"))
                    sst_fieldname = 'analysed_sst'
                    sst_uncertainty_fieldname = 'analysis_uncertainty'
                in_path = os.path.join(day_dir,filename)

                out_folder = os.path.join(args.output_dir, year, month, day)
                out_path = os.path.join(out_folder, filename)

                # open the datasets
                logger.info("Reading %s"%(in_path))
                ds_orig = xr.open_dataset(in_path).fillna(0)
                ds_corrected  = xr.open_dataset(out_path).fillna(0)

                sst_diffs = ds_corrected[sst_fieldname].data - ds_orig[sst_fieldname].data
                uncertainty_diffs = ds_corrected[sst_uncertainty_fieldname].data - ds_orig[sst_uncertainty_fieldname].data

                max_sst_diff = np.max(sst_diffs)
                min_sst_diff = np.min(sst_diffs)

                # track global min/max of sst values
                if max_diff is None or max_sst_diff > max_diff:
                    max_diff = max_sst_diff
                    max_diff_date = year+"-"+month+"-"+day

                if min_diff is None or min_sst_diff < min_diff:
                    min_diff = min_sst_diff
                    min_diff_date = year + "-" + month + "-" + day

                # track yearly min/max of sst values
                if max_year_diff is None or max_sst_diff > max_year_diff:
                    max_year_diff = max_sst_diff
                    max_year_diff_date = year+"-"+month+"-"+day

                if min_year_diff is None or min_sst_diff < min_year_diff:
                    min_year_diff = min_sst_diff
                    min_year_diff_date = year + "-" + month + "-" + day

                max_sst_uncertainty_diff = np.max(uncertainty_diffs)
                min_sst_uncertainty_diff = np.min(uncertainty_diffs)

                # track global min/max of uncertainty values
                if max_uncertainty_diff is None or max_uncertainty_diff > max_sst_uncertainty_diff:
                    max_uncertainty_diff = max_sst_uncertainty_diff
                    max_uncertainty_diff_date = year + "-" + month + "-" + day

                if min_uncertainty_diff is None or min_uncertainty_diff > min_sst_uncertainty_diff:
                    min_uncertainty_diff = min_sst_uncertainty_diff
                    min_uncertainty_diff_date = year + "-" + month + "-" + day

                # track yearly min/max of uncertainty values
                if max_year_uncertainty_diff is None or max_sst_uncertainty_diff > max_year_uncertainty_diff:
                    max_year_uncertainty_diff = max_sst_uncertainty_diff
                    max_year_uncertainty_diff_date = year + "-" + month + "-" + day

                if min_year_uncertainty_diff is None or min_sst_uncertainty_diff < min_year_uncertainty_diff:
                    min_year_uncertainty_diff = min_sst_uncertainty_diff
                    min_year_uncertainty_diff_date = year + "-" + month + "-" + day

                ds_orig.close()
                ds_corrected.close()
        logger.info("YEAR %s: Min SST increase: %f %s, Max SST increase: %f %s, Min Uncertainty Increase: %f %s, Max Uncertainty Increase: %f %s"%(year,min_year_diff,min_year_diff_date,max_year_diff,max_year_diff_date,min_year_uncertainty_diff,min_year_uncertainty_diff_date,max_year_uncertainty_diff,max_year_uncertainty_diff_date))

    logger.info(
        "Global: Min SST increase: %f %s, Max SST increase: %f %s, Min Uncertainty Increase: %f %s, Max Uncertainty Increase: %f %s" % (
        min_diff, min_diff_date, max_diff, max_diff_date, min_uncertainty_diff,
        min_uncertainty_diff_date, max_uncertainty_diff, max_uncertainty_diff_date))




