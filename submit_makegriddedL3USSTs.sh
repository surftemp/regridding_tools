#!/bin/bash
#SBATCH --partition=short-serial
#SBATCH -o makegriddedL3USSTs_%j.out
#SBATCH -e makegriddedL3USSTs_%j.err
#SBATCH --time=24:00:00

conda activate regridding_tools
python ./makegriddedL3USSTs.py 5 5 monthly AVHRR07_G --year 1982 --sst_l3u_path /neodc/esacci/sst/data/CDR_v2/AVHRR/L3U/v2.1 --out_path /gws/nopw/j04/esacci_sst/output/l3_regrid/CDR2.1
