# Publish L4 SST data to AWS opendata

For example see: https://registry.opendata.aws/mur/  

## In this folder

* publication program (publish.py)

Appends data for a given year to a zarr data store.

* test program (ztest.py)

Tests that the published zarr data store is consistent with the source netcdf4 data

* consumer sample code (zdump.py)

Checks that a consumer can connect to the zarr data store and extract data

