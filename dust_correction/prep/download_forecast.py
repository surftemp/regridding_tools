#!/usr/bin/env python

import cdsapi

c = cdsapi.Client(progress=False)

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--start-date", type=str, default="2021-01-01", help="start date for download in YYYY-MM-DD format")
parser.add_argument("--end-date", type=str, default="2021-07-01",  help="end date for download in YYYY-MM-DD format")

args = parser.parse_args()

c.retrieve(
    'cams-global-atmospheric-composition-forecasts',
    {
        'date': '%s/%s'%(args.start_date,args.end_date),
        'type': 'forecast',
        'format': 'netcdf',
        'variable': [
            'dust_aerosol_optical_depth_550nm',
            'vertically_integrated_mass_of_dust_aerosol_0.03-0.55um',
            'vertically_integrated_mass_of_dust_aerosol_0.55-9um',
            'vertically_integrated_mass_of_dust_aerosol_9-20um'
        ],
        'time': '12:00',
        'leadtime_hour': '0',
    },
    'dust_forecast.nc')


