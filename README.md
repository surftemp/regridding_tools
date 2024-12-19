# regridding_tools

Tools for regridding high-resolution surface temperature data to more manageable resolutions.

These tools were originally developed to work with the ESA CCI SST dataset (https://climate.esa.int/en/projects/sea-surface-temperature/) produced by the University of Reading, the Met Office and other partners, and are used by the University of Reading's surftemp data portal (https://surftemp.net).

Main software provided in this repo:

* [Global regridding tools](global_regridding/README.md)
* [Tools for visualising regridded data](plotting/README.md)
* [Publish L4 SST data to AWS opendata](publish_to_zarr/README.md)
* [Relice into zarr format for efficient timeseries computation](timeseries_region/README.md)
