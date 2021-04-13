# Plotting

Tools for visualising regridded data.

## compare_L3_had.py

Create plots comparing regridded L3U data to HadSST and regridded L4 data. Panel plots of maps and histograms for each
month are created as well as whole series maps, histograms and time series.

```
usage: compare_L3_Had.py [-h] [--hadfile HADFILE] [--hadfileuncert HADFILEUNCERT] [--titlestr TITLESTR] [--umax UMAX]
                         [--xbins XBINS] [--xmin XMIN] [--xmax XMAX] [--vmax VMAX] [--nol4]
                         l3path l4path outPicPath

Create plots comparing regridded L3U files to regridded L4 and HadSST

positional arguments:
  l3path                Path to regridded L3U files.
  l4path                Path to regridded L4 files.
  outPicPath            Path to create output plots in.

optional arguments:
  -h, --help            show this help message and exit
  --hadfile HADFILE     Path to HadSST file with actuals and median.
  --hadfileuncert HADFILEUNCERT
                        Path to HadSST file with total uncertainty.
  --titlestr TITLESTR   A title string to prefix to the titles of the plots.
  --umax UMAX           Maximum uncertainty to include in best HadSST data. Default is 0.35
  --xbins XBINS         Number of x bins to use for histograms. Default is 80.
  --xmin XMIN           Minimum of x range for histograms. Default is -1.75.
  --xmax XMAX           Maximum of x range for histograms. Default is 1.75.
  --vmax VMAX           Maximum of colorscale for maps. Default is 0.8 degrees C.
  --nol4                Do not compare to the L4 dataset only to the HadSST.
```

### Example

Create plots comparing NOAA 6 data to HadSST specifying not to compare to the L4 files. The path to the L4 files is
still expected.

```
python compare_L3_Had.py /path/to/l3u/files/AVHRR06_G /path/to/l4/files /path/to/output/files
    --hadfile /path/to/HadSST4/HadSST.4.0.0.0_actuals_median.nc
    --hadfileuncert /path/to/HadSST4/HadSST.4.0.0.0_total_uncertainty.nc --titlestr "N06 v3.3.2" --nol4
```