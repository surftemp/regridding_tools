#!/bin/bash 

#SBATCH --partition=long-serial 
#SBATCH --job-name=1980s
#SBATCH -o 1980s.out 
#SBATCH -e 1980s.err 
#SBATCH -t 96:00:00
#SBATCH --mem=8192

conda activate aws_env

python publish.py --start-year 1981 --end-year 1982 --verbose
python publish.py --start-year 1983 --end-year 1983 --verbose
python publish.py --start-year 1984 --end-year 1984 --verbose
python publish.py --start-year 1985 --end-year 1985 --verbose
python publish.py --start-year 1986 --end-year 1986 --verbose
python publish.py --start-year 1987 --end-year 1987 --verbose
python publish.py --start-year 1988 --end-year 1988 --verbose
python publish.py --start-year 1989 --end-year 1989 --verbose
