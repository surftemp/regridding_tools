#!/bin/bash 

#SBATCH --partition=long-serial 
#SBATCH --job-name=zarr_upload
#SBATCH -o upload.out 
#SBATCH -e upload.err 
#SBATCH -t 128:00:00
#SBATCH --mem=100000

conda activate aws_env

# suggest uploading a decade at a time which should complete within a couple of days
# uploads must be consecutive in time
START_YEAR=1981
END_YEAR=1990

for ((year=$START_YEAR;year<=$END_YEAR;year++)); do
    python publish.py --start-year $START_YEAR --end-year $END_YEAR --verbose
done
