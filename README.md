# regridding_tools
Tools for regridding high-resolution surface temperature data to more manageable resolutions.

## makegriddedSSTs.py
A tool to regrid daily 0.05 degree resolution L4 SST files to a lower spatial and/or temporal resolution. It will
process the data for a particular year. A start month, start day, end month and end day can be specified. These must be
consistent with the temporal resolution, which can be annual, monthly, nominal 10-day, nominal 5-day or N day from the
beginning of the year. The spatial resolution can be any multiple of 0.05 that is also a factor of 3600 for the
latitude or 7200 for the longitude up to 5 degrees. This is because the regridding is done by averaging each box of
cells at the higher resolution that correspond to a single cell at the lower resolution. The outputs can be anomalies or
absolute temperatures. In order to calculate the L4 climatology is regridded to the same resolution as the data and
subtracted from the absolute temperature. By default the fraction of sea ice in each particular cell is output as well,
but this can be optionally be omitted. A fraction of sea ice above which to exclude a cell from the calculation can
be specified. However, in calculating the climatology this was assumed to be 1.0 and this cannot be changed. The
regridded uncertainties are also calculated and parameters specifying the spatial (lambda) and temporal (tau) scales on
which the uncertainties are assumed to be fully calculated can be specified. By default these are 3 degrees and 7 days
respectively. The paths to the top level directories for the input data and climatology can be specified. These are
assumed to be ordered by levels of subdirectories for year, month and day. The directory in which to save the output
files can be specified and the output files can optionally be zipped into a single file with a user specifed name.