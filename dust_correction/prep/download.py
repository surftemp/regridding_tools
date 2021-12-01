#!/usr/bin/env python

import cdsapi

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--start-year", type=int, default=2020)
parser.add_argument("--end-year", type=int, default=2021)

args = parser.parse_args()

c = cdsapi.Client(progress=False)

for year in range(args.start_year, args.end_year):
    try:
        c.retrieve(
            'cams-global-reanalysis-eac4',
            {
                'date': '%04d-01-01/%04d-12-31' % (year, year),
                'format': 'netcdf',
                'variable': [
                    'dust_aerosol_optical_depth_550nm',
                    'vertically_integrated_mass_of_dust_aerosol_0.03-0.55um',
                    'vertically_integrated_mass_of_dust_aerosol_0.55-9um',
                    'vertically_integrated_mass_of_dust_aerosol_9-20um',
                ],
                'time': '12:00',
            },
            'dust_%04d.nc' % year)
    except Exception as ex:
        print(str(ex))

