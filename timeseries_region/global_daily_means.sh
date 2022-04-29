#!/bin/bash

# compute the daily global mean SSTs

# choose input data - dust corrected or not

# original
# INPUT_SSTS=/data/input/reslice/

# dust corrected
INPUT_SSTS=/data2/input/dust/reslice

# configure output path
OUTPUT_PATH=global_sst_timeseries.csv

# activate environment for running
conda activate cfpy37

python ./maketimeseriesSSTs.py -180.0 180.0 -90.0 90.0 daily --start_year 1981 --end_year 2021 --start_month 9 --start_day 1 --end_month 12 --end_day 31 --f_max 1.0 --tau 3 --spatial_lambda 1.0 --in_path $INPUT_DATA --out_path $OUTPUT_PATH