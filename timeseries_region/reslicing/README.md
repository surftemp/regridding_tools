# Spatial reslicing task

This task is described in [Space-Sliced-SST-L4.docx](https://livereadingac-my.sharepoint.com/:w:/r/personal/c_j_merchant_reading_ac_uk/Documents/EO%20Scientific%20Programmers/Regridding/Space-Sliced-SST-L4.docx?d=w9c4652bfa2d04309bd96101696a42a7b&csf=1&web=1&e=jIymUt)

In summary, the objective is to reprocess the 0.05 degree daily CCI/C3S SST and climatology data to obtain equivalent data which can be 
more efficiently accessed to retrieve data for a particular area defined by a cell

This tool can reslice the data to provide output SST files in zarr format where each file holds a year's worth of daily 0.05 degree resolution data 
while storing the data in 5 degree * 5 degree * 7 day chunks.  

```
<output-folder>/2018.zarr
<output-folder>/2017.zarr
...
```

The tool also outputs climatology files according to a similar naming system.

```
<output-folder>/climatology.zarr
```

## Running on the regridding VM using supplied scripts

Run reslice.py configured for the file layout of the regridding VM using `./reslice.sh <START-YEAR> <END-YEAR>` (original SST files) and `./reslice_dust.sh <START-YEAR> <END-YEAR>` (dust corrected SST files).

## Options

Use the `--input-cci`, `--input-c3s` and `--input-climatology` to specify the locations of the various input file datasets

The defaults for these options should match the locations on CEDA/JASMIN under `/neodc`

Use the `--output-folder` to specify the folder under which the output files will be written.

Use the `--output-resolution` option to control the size of the spatial area (in degrees) described by each chunk of the output data file.  This defaults to 5 (degrees)

You must provide the `--year` option to specify the year or the climatology to process (see below)

### SST data reslicing options

Use `--year` to reslice a given week (1-52) within a year and output a resliced .zarr folder for that year.

Output folders will be named `<output-folder>/<year>.zarr`

### Climatology data reslicing options

Same as above for SST data but supply value `climatology` instead of the year to reslice the climatology.  

The output folder will be named `<output-folder>/climatology.zarr`

