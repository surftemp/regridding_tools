#!/bin/bash 

#SBATCH --partition=long-serial 
#SBATCH --job-name=1980s
#SBATCH -o 1980s.out 
#SBATCH -e 1980s.err 
#SBATCH -t 96:00:00
#SBATCH --mem=8192

conda activate aws_env

python publish.py --start-year 1990 --end-year 1990 --verbose
python publish.py --start-year 1991 --end-year 1991 --verbose
python publish.py --start-year 1992 --end-year 1992 --verbose
python publish.py --start-year 1993 --end-year 1993 --verbose
python publish.py --start-year 1994 --end-year 1994 --verbose
python publish.py --start-year 1995 --end-year 1995 --verbose
python publish.py --start-year 1996 --end-year 1996 --verbose
python publish.py --start-year 1997 --end-year 1997 --verbose
python publish.py --start-year 1998 --end-year 1998 --verbose
python publish.py --start-year 1999 --end-year 1999 --verbose
