# global_regridding

Tools for regridding high resolution SST files to lower spatial and temporal resolutions over a global extent with a
file for each time slice.

## makegriddedSSTs.py

A tool to regrid daily 0.05 degree resolution L4 SST files to a lower spatial and/or temporal resolution. It will
process the data for a particular year. A start month, start day, end month and end day can be specified. These must be
consistent with the temporal resolution, which can be annual, monthly, nominal 10-day, nominal 5-day or N day from the
beginning of the year. The spatial resolution can be any multiple of 0.05 degrees that is also a factor of 180 for the
latitude or 360 for the longitude up to 5 degrees. This is because the regridding is done by averaging each box of cells
at the higher resolution that correspond to a single cell at the lower resolution. The outputs can be anomalies or
absolute temperatures. In order to calculate the former the L4 climatology is regridded to the same resolution as the
data and subtracted from the absolute temperature. The fraction of sea at the higher resolution in each cell is output.
By default the fraction of sea ice in each particular cell is output as well, but this can be optionally be omitted. A
fraction of sea ice above which to exclude a cell from the calculation can be specified. However, in calculating the
climatology this was assumed to be 1.0 and this cannot be changed. The regridded uncertainties are also calculated and
parameters specifying the spatial (lambda) and temporal (tau) scales on which the uncertainties are assumed to be fully
calculated can be specified. By default these are 3 degrees and 7 days respectively. Convenience fields for the year,
mmonth, day and day of year are also output. The paths to the top level directories for the input data and climatology
can be specified. These are assumed to be ordered by levels of subdirectories for year, month and day. The directory in
which to save the output files can be specified and the output files can optionally be zipped into a single file with a
user specifed name.

```
usage: makegriddedSSTs.py [-h] [--year YEAR] [--start_month START_MONTH] [--start_day START_DAY] [--end_month END_MONTH]
                          [--end_day END_DAY] [--anomalies] [--no_sea_ice_fraction] [--f_max F_MAX] [--tau TAU]
                          [--spatial_lambda SPATIAL_LAMBDA] [--sst_cci_analysis_l4_path SST_CCI_ANALYSIS_L4_PATH]
                          [--c3s_sst_analysis_l4_path C3S_SST_ANALYSIS_L4_PATH]
                          [--sst_cci_climatology_path SST_CCI_CLIMATOLOGY_PATH] [--out_path OUT_PATH]
                          [--zip_name ZIP_NAME]
                          lon_resolution lat_resolution time_resolution

Regrid L4 SST data for a particular year.

positional arguments:
  lon_resolution        The target longitude resolution. This must be a multiple of 0.05 degrees and a factor of 360
                        degrees. It must also not be greater than 5 degrees.
  lat_resolution        The target latitude resolution. This must be a multiple of 0.05 degrees and a factor of 180
                        degrees. It must also not be greater than 5 degrees.
  time_resolution       The target time resolution. This can be 'annual', 'monthly', '10-day' for dekads, '5-day' for
                        pentads or an integer for regular N day regridding aligned with the start of the year.

optional arguments:
  -h, --help            show this help message and exit
  --year YEAR           The year of the data to be regridded. Default is 1982.
  --start_month START_MONTH
                        The first month of the data to be regridded. Default is 1.
  --start_day START_DAY
                        The first day of the month of the data to be regridded. Default is 1.
  --end_month END_MONTH
                        The final month of the data to be regridded. Default is 12.
  --end_day END_DAY     The final day of the month of the data to be regridded. Default is 31.
  --anomalies           Output anomalies instead of absolute SSTs.
  --no_sea_ice_fraction
                        Do not output the sea ice fraction.
  --f_max F_MAX         The fraction of sea ice above which a cell is ignored in calculating the regridded SST SST.
                        Default is 1.0. When calculating anomalies, the climatology is always calculated with an f_max
                        of 1.0 regardless of this value.
  --tau TAU             Timescale within which errors are assumed to be fully correlated in days. Default is 7 days.
  --spatial_lambda SPATIAL_LAMBDA
                        Spatial scale within which errors are assumed to be fully correlated in degrees. Default is 3.0
                        degrees.
  --sst_cci_analysis_l4_path SST_CCI_ANALYSIS_L4_PATH
                        Path to the SST CCI Analysis Level 4 input data.
  --c3s_sst_analysis_l4_path C3S_SST_ANALYSIS_L4_PATH
                        Path to the C3S SST Analysis Level 4 input data.
  --sst_cci_climatology_path SST_CCI_CLIMATOLOGY_PATH
                        Path to the SST CCI Climatology input data.
  --out_path OUT_PATH   The path in which to write the output file of regridded data.
  --zip_name ZIP_NAME   Combine all output files into a single zip file with this name.
```

### Example

Regrid SST CCI files for 1982 to 0.25 x 0.25 degree spatial resolution and monthly temporal resolution with custom
paths.

```
python makegriddedSSTs.py 0.25 0.25 monthly --year 1982 --start_month 1 --start_day 1 --end_month 12 --end_day 31
    --sst_cci_analysis_l4_path /path/to/sst_cci/files --sst_cci_climatology_path /path/to/climatology/files
    --out_path /path/to/output/files
```

## makegriddedL3USSTs.py

A tool to regrid 0.05 degree resolution L3U SST files to lower spatial and temporal resolution. It will process the data
for a particular year. A start month and end month can be specified. These must be consistent with the temporal
resolution, which can be annual, monthly or nominal 10-day. The spatial resolution can be any multiple of 0.05 degrees
that is also a factor of 180 for the latitude or 360 for the longitude between 0.25 and 5 degrees. This is because
regridding is done by averaging each box of cells at the higher resolution that correspond to a single cell at the lower
resolution. The output can be anomalies or absolute temperatures. It is the anomalies that are regridded in this tool
and the regridded climatology is added to the regridded anomalies to obtain the regridded absolute temperatures. The
number of observations that contributed to each cell in the regridded SST is also output. A minimum quality level both
for the cell and the file can be specified below which the data is ignored. Convenience fields for the year, mmonth, day
and day of year are also output. The paths to the top level directories for the input data and climatology can be
specified. These are assumed to be ordered by levels of subdirectories for sensor, year, month and day for the L3U files
and year, month and day for the climatology. The directory in which to save the output files can be specified and the
output files can optionally be zipped into a single file with a user specified name.

```
usage: makegriddedL3USSTs.py [-h] [--year YEAR] [--start_month START_MONTH] [--end_month END_MONTH]
                             [--anomalies] [--min_qlevel MIN_QLEVEL] [--min_flevel MIN_FLEVEL]
                             [--sst_l3u_path SST_L3U_PATH]
                             [--sst_cci_climatology_path SST_CCI_CLIMATOLOGY_PATH]
                             [--out_path OUT_PATH] [--zip_name ZIP_NAME]
                             lon_resolution lat_resolution time_resolution sensor

Regrid L3U SST data for a particular year.

positional arguments:
  lon_resolution        The target longitude resolution. This must be a multiple of 0.05 degrees
                        and a factor of 360 degrees. It must also not be greater than 5 degrees.
  lat_resolution        The target latitude resolution. This must be a multiple of 0.05 degrees and
                        a factor of 180 degrees. It must also not be greater than 5 degrees.
  time_resolution       The target time resolution. This can be 'annual', 'monthly', '10-day' for
                        dekads.
  sensor                The sensor for which to perform regridding.

optional arguments:
  -h, --help            show this help message and exit
  --year YEAR           The year of the data to be regridded. Default is 1981.
  --start_month START_MONTH
                        The first month of the data to be regridded. Default is 1.
  --end_month END_MONTH
                        The final month of the data to be regridded. Default is 12.
  --anomalies           Output anomalies instead of absolute SSTs.
  --min_qlevel MIN_QLEVEL
                        The minimum quality level of the L3U cell SST for inclusion. Default is 4.
                        If 4 then quality levels 4 and 5 are included. Minimum is 3 such that
                        quality levels 3, 4 and 5 anre included. Maximum is 5 such that only
                        quality level 5 is included.
  --min_flevel MIN_FLEVEL
                        The minimum file quality level. Default is 3. An integer from 0 to 3.
  --sst_l3u_path SST_L3U_PATH
                        Path to the SST Level 3 Uncollated input data.
  --sst_cci_climatology_path SST_CCI_CLIMATOLOGY_PATH
                        Path to the SST CCI Climatology input data.
  --out_path OUT_PATH   The path in which to write the output file of regridded data.
  --zip_name ZIP_NAME   Combine all output files into a single zip file with this name.
```

### Example

Regrid files for the sensor AVHRRMTA_G to 0.5 x 0.5 degree spatial resolution and nominal 10-day temporal resolution
for March 2019 with custom paths.

```
python makegriddedL3USSTs.py 0.5 0.5 10-day AVHRRMTA_G --year 2019 --start_month 3 --end_month 3
    --sst_l3u_path /path/to/l3u/files --sst_cci_climatology_path /path/to/climatology/files
    --out_path /path/to/output/files
```
### submit_makegriddedL3USSTs.py

Submit batch jobs of the L3U regridder on JASMIN for all years for a particular sensor.

```
usage: Submit the L3U regridder for all the years for a particular sensor. [-h] sensor sst_l3u_path out_path

positional arguments:
  sensor        The name of the sensor. The input files should be in a directory with this name.
  sst_l3u_path  The top level directory containing the L3U input files.
  out_path      The path to the directory in which to save the output files.

optional arguments:
  -h, --help    show this help message and exit
```

### Example

Submit batch jobs for all years for the AVHRRMTA_G sensor. Files are expected to be in:
`/path/to/input/files/AVHRRMTA_G`

```
python submit_makegriddedL3USSTs.py AVHRRMTA_G /path/to/input/files /path/to/output/files
```