# dust contamination - task to prepare data

This task is to obtain netcdf4 files containing 4 dust aerosol optical depth (AOD) variables from the Atmospheric Data Store (ADS) and create a climatology.

The data will then be used to study how dust contamination affects sea surface temperature estimates.

The required variables are:

* duaod550   Dust Aerosol Optical Depth at 550nm
* aermssdul  Vertically integrated mass of dust aerosol (9 - 20 um)
* aermssdum  Vertically integrated mass of dust aerosol (0.55 - 9 um)
* aermssdus  Vertically integrated mass of dust aerosol (0.03 - 0.55 um)

(1) download data from the ADS using download.py

To register for access and set up a python environment and key, see https://ads.atmosphere.copernicus.eu/api-how-to

Then run:

```
python download.py --start-year YYYY --end-year YYYY
```

This should download files from the ADS containing the 4 dust measurements, regridded to 0.5 deg, at 12.00 every day.  One file is generated for each year.

```
dust_2003.nc
...
dust_2020.nc
```

(2) Download forecast data for the current year using download_forecast.py

```
python download_forecast.py --start-date YYYY-MM-DD --end-year YYYY-MM-DD
```

This will download the data to file `dust_forecast.nc`

(3) regrid/interpolate the files using cf-python to 0.5 degree 

You will need cfpython installed to run the `regrid.py` script.  I used the following to create a conda environment:

```
conda create -n cfpy python=3.8
conda activate cfpy
conda install udunits2=2.2.25
pip install cf-python
conda install -c conda-forge esmpy
```

You will also need a reference netCDF4 file with the desired lat/lon grid.  I used the reference file `20121011120000-ESACCI-L3C_GHRSST-SSTsubskin-AMSR2-CDR2.0_night-v02.0-fv01.0.nc`

Then run:
```
python regrid.py dust_2003.nc regridded2003_0p5.nc <reference file>
python regrid.py dust_2004.nc regridded2003_0p5.nc <reference file>
...
python regrid.py output_2018_0p5 regridded2018_0p5.nc <reference file>
```

(4) split regridded data into daily files

This converts the regridded*.nc files into daily files named daily/YEAR/day1.nc, daily/YEAR/day2.nc, etc

You will need cfpython installed to run the `split.py` script (see step 2 for install method)

Then run:

```
python split.py
```

(5) create climatology

The script `create_dust_climatology.py` works through groups of daily files to compute a climatology from the years 2003 to 2017 inclusive and a 5 day window.

You will need cfpython installed to run the `create_dust_climatology.py` script (see step 2 for install method)

This script is based on the original script https://github.com/surftemp/common/blob/master/scripts/cci/create_climatology.py where most changes relate to extra statistics added and a conversion to python3

The script will output climatology files, one per day of the year in the `climatology` output folder when run with the following command

```
mkdir climatology
python create_dust_climatology.py --start 2003 --stop 2017 --window 2 daily climatology
```
