# Time Series Extraction for SST data

This task is described in [Local_Mean_SST_Timeseries_Service.docx](https://livereadingac-my.sharepoint.com/:w:/r/personal/c_j_merchant_reading_ac_uk/_layouts/15/Doc.aspx?sourcedoc=%7BC8F55706-C4D7-4E40-A886-4C56510499C1%7D&file=Local_Mean_SST_Timeseries_Service.docx)

This folder and sub-folders contain scripts for extracting time series data from SST for a given bounding box and temporal cycle

The code has been moved from the `general-tasks` repo.

Script `maketimeseriesSSTs.py` can be invoked on the command line (or invoked via its `makeTimeSeriesSSTs` function) for example:

```
python extract_timeseries.py 6 7 -33 -32 daily --start_year 1982 --end_year 1988  --out_path test.nc
makeTimeSeriesSSTs Processing 60%   ############################
makeTimeSeriesSSTs Processing 90%   ########################################
makeTimeSeriesSSTs Processing 99%   ############################################
```

The first four arguments represent the spatial bounding box, given as MIN-LON, MAX-LON, MIN-LAT, MAX-LAT
The fifth argument represents the temporal period one of 'daily','<N-days>',monthly','5-day','10-day'.  

Run the script with `--help` to get a listing of all supported options.

To obtain output in .csv format, supply a file with suffix `.csv` using the --out_path option.

The calculations yielding sea temperature, sea ice fraction, sea temperature anomaly and standard error are identical to that performed by the makegridedSST.py but use a different implementation.

## Reslicing to create the input data

This script uses as input data, SST data which has been transformed into a "resliced" dataset stored .zarr format, with one .zarr folder per year and a climatology.zarr to hold the climatology

The chunking system for the .zarr files is designed to facilitate the efficient retrieval of subsets of data needed to create time series.

The [reslicing folder](reslicing/README.md) contains a script and associated files for building this resliced data.  See its README for more information.

## Verification

To verify the results for time series extraction, we compare the values for the time series with that obtained for equivalent cells within a 
fully regridded dataset at matching spatial and temporal resolution.  The calculations specified in the task descriptions are the same but makegriddedSSTs.py and maketimeseriesSSTs.py 
have different implementations (they do both rely on the cfpython and numpy libraries).

Tests are implemented as part of the sst-services repo - see [https://github.com/surftemp/sst-services/tree/master/test](https://github.com/surftemp/sst-services/tree/master/test)
 