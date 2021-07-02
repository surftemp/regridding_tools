# dust contamination - task to prepare data

This task is to obtain netcdf4 files containing 4 dust aerosol optical depth (AOD) variables from ECMWF and create a climatology.

The data will then be used to study how dust contamination affects sea surface temperature estimates.

The required variables are:

* duaod550   Dust Aerosol Optical Depth at 550nm
* aermssdul  Vertically integrated mass of dust aerosol (9 - 20 um)
* aermssdum  Vertically integrated mass of dust aerosol (0.55 - 9 um)
* aermssdus  Vertically integrated mass of dust aerosol (0.03 - 0.55 um)

(1) download data from ECMWF using download.py

See the following URL for helping to create the MARS request.

https://apps.ecmwf.int/data-catalogues/cams-reanalysis/?axis_param=209.210&class=mc&expver=eac4&stream=oper&type=an&year=2003&month=jan&levtype=sfc

A python environment with the ECMWF API installed is required.  https://www.ecmwf.int/en/forecasts/access-forecasts/ecmwf-web-api has more details on obtaining and installing a key enabling use of the API.

I used conda to create a python 3.5 environment for this using:

```
conda create -n dust2 python=3.5
conda activate dust2
pip install ecmwf-api-client
```

Then run:

```
python download.py
```

This should download files from ECMWF containing the 4 dust measurements, regridded to 0.5 deg, at 12.00 every day:

```
output2003_0p5.nc
...
output2018_0p5.nc
```

Note: the ECMWF API is no longer availale and the download will have to be gotten from the ADS.  Issue https://github.com/surftemp/general-tasks/issues/17 is created to track this.

(2) regrid the files using cfpython (spherical regridding / blinear)

The ECMWF downloaded files have different boundaries (lat from -90.0 to +90, lon from 0 to 359.5).  
We need instead the boundaries (lat from -89.75 to 89.75, lon from -179.75 to 179.75).

You will need cfpython installed to run the `regrid.py` script.  I used the following to create a conda environment:

```
conda create -n cfpy python=3.8
conda activate cfpy
conda install -c conda-forge cartopy
conda install udunits2=2.2.25
pip install cf-python cf-plot
conda install -c conda-forge esmpy
```

You will also need a reference netCDF4 file with the desired lat/lon grid.  I used the reference file `20121011120000-ESACCI-L3C_GHRSST-SSTsubskin-AMSR2-CDR2.0_night-v02.0-fv01.0.nc`

Then run:
```
python regrid.py output_2003_0p5 regridded2003_0p5.nc <reference file>
python regrid.py output_2003_0p5 regridded2003_0p5.nc <reference file>
...
python regrid.py output_2018_0p5 regridded2018_0p5.nc <reference file>
```

(3) split regridded data into daily files

This converts the regridded*.nc files into daily files named daily/YEAR/day1.nc, daily/YEAR/day2.nc, etc

You will need cfpython installed to run the `split.py` script (see step 2 for install method)

Then run:

```
python split.py
```

(4) create climatology

The script `create_dust_climatology.py` works through groups of daily files to compute a climatology from the years 2003 to 2017 inclusive and a 5 day window.

You will need cfpython installed to run the `create_dust_climatology.py` script (see step 2 for install method)

This script is based on the original script https://github.com/surftemp/common/blob/master/scripts/cci/create_climatology.py where most changes relate to extra statistics added and a conversion to python3

The script will output climatology files, one per day of the year in the `climatology` output folder when run with the following command

```
mkdir climatology
python create_dust_climatology.py --start 2003 --stop 2017 --window 2 daily climatology
```
