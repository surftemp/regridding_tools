#!/bin/bash 

#SBATCH --partition=long-serial 
#SBATCH --job-name=1980s
#SBATCH -o 1980s.out 
#SBATCH -e 1980s.err 
#SBATCH -t 96:00:00

conda activate aws_env

python publish.py --start-year 1981 --end-year 1990 --verbose
