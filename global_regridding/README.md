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