#!/bin/bash 

#SBATCH --partition=long-serial 
#SBATCH --job-name=zarr_test
#SBATCH -o test.out 
#SBATCH -e test.err 
#SBATCH -t 24:00:00
#SBATCH --mem=16384

conda activate aws_env

START_YEAR=1981
END_YEAR=1989

for ((year=$START_YEAR;year<=$END_YEAR;year++)); do
    if [ "$year" = "1981" ]; then
	start_date="1981-09-01"
    else
        start_date="$year-01-01"
    fi
    end_date="$year-12-31"
    python ztest.py --start-date "$start_date" --end-date "$end_date" --verbose
done
