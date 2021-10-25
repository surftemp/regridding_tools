#!/bin/bash 

#SBATCH --partition=long-serial 
#SBATCH --job-name=zarr_test
#SBATCH -o test2.out 
#SBATCH -e test2.err 
#SBATCH -t 96:00:00
#SBATCH --mem=10000

conda activate aws_env

START_YEAR=1990
END_YEAR=1999

for ((year=$START_YEAR;year<=$END_YEAR;year++)); do
    if [ "$year" = "1981" ]; then
	start_date="1981-09-01"
    else
        start_date="$year-01-01"
    fi
    end_date="$year-12-31"
    python ztest.py --start-date "$start_date" --end-date "$end_date" --verbose
done
