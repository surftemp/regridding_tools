import argparse
import glob
import os
import subprocess

_submission_script = """#!/bin/bash
#SBATCH --partition=short-serial
#SBATCH -o makegriddedL3USSTs_{0}_{1}.out
#SBATCH -e makegriddedL3USSTs_{0}_{1}.err
#SBATCH --time=12:00:00
#SBATCH --mem=2000

conda activate regridding_tools
python ./makegriddedL3USSTs.py 5 5 monthly {0} --year {1} --sst_l3u_path {2} --out_path {3}
"""
_submission_script_file = 'tmp_submit_makegriddedL3USSTs.sh'

def submit_makegriddedL3USSTs(sensor, sst_l3u_path, out_path):
    """
    Submit the L3U regridder for all the years for a particular sensor.
    :param sensor:
    :param sst_l3u_path:
    :param out_path:
    :return:
    """
    out_path = out_path + sensor
    os.makedirs(out_path, exist_ok=True)
    years = glob.glob(os.path.join(sst_l3u_path, sensor, '[0-9]' * 4))
    for year in years:
        with open(_submission_script_file, 'w') as f:
            f.write(_submission_script.format(sensor, year, sst_l3u_path, out_path))
        subprocess.run(['sbatch', _submission_script_file])

if __name__ == '__main__':
    parser = argparse.ArgumentParser('Submit the L3U regridder for all the years for a particular sensor.')
    parser.add_argument('sensor', help='The name of the sensor. The input files should be in a directory with this name.')
    parser.add_argument('sst_l3u_path', help='The top level directory containing the L3U input files.')
    parser.add_argument('out_path', help='The path to the directory in which to save the output files.')
    args = parser.parse_args()

    submit_makegriddedL3USSTs(args.sensor, args.sst_l3u_path, args.out_path)
