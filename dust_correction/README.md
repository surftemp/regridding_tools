# Applying dust correction

The aim is to create new sets of C3S/CCI/Climatology files.  For more details of using these with the regridding service,
see the relevant section of the SETUP document: https://github.com/surftemp/sst-services/blob/master/docs/SETUP.md

The correction program [correct.py](correct.py) uses the `DustAdjuster` class from from https://github.com/surftemp/plots-dev/blob/master/projects/cci3/l4_adjust/adjust_md.py

The program to rebuild the climatology [climatology.py](climatology.py) is based on https://github.com/surftemp/common/blob/master/scripts/cci/create_climatology.py

## Setup

Obtain files `dustAdjustmentParametersAllYears_1.0.nc` and `L4_spike_adjustment_1.0.nc` (I've used the files from https://livereadingac-my.sharepoint.com/:f:/r/personal/c_j_merchant_reading_ac_uk/Documents/TempDataSharing?csf=1&web=1&e=rum2hJ)

### Setup a conda environment

```
conda env create -f environment.yml
conda activate xarr
```

## Applying the correction - manual usage

To correct a single year use:

```
python correct.py --data-dir `pwd` --input-year 2019 --input-dir /data/input/ICDR_v2/Analysis/L4/v2.0 --output-dir /data/tmp/dust
```

Note there are also `--input-month` and `--input-day` options for quickly converting a small number of files

## Applying the correction - batch usage on the regridding VM

Run the script `batch.sh` to apply the correction to all data, creating new datasets with the same folder structure and naming as the C3S and CCI input data.  If running in a different environment, you'll need to customise the script

## Rebuilding the climatology from years 1982 to 2010

Once a set of corrected SST files has been created, a climatology should be created.

Run the script `climatology.sh` to do this (set up with suitable paths for the regridding VM).  This uses the `climatology.py` program which has been tested to check that the settings are compatible with that used to create the existing climatology, namely:

* period from 1982 to 2010
* window size of 2
* in the first 2 days of a year, window includes last days from previous year
* in the last 2 days of a yearm window includes first days from next year
* ignore day 366 in leap years

## Checking differences

The script `check_diffs.py` can be used to compare the original and corrected SST files and find the largest and smallest increases in uncertainty and SST
by year and over all the data
