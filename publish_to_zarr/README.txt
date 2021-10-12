Sea Surface Temperature Daily Analysis: European Space Agency Climate Change Initiative product version 2.1
===========================================================================================================

Summary
-------

Global daily-mean sea surface temperatures, presented on a 0.05° latitude-longitude grid, with gaps between available daily observations filled by statistical means, spanning late 1981 to recent time.

Synthesised from multiple Earth orbiting satellites carrying infrared imagers from thousands of billions of individual measurements. Underlying observation resolution ranges from 1 to 20 km, and after gap filling the typical feature resolution is ~20 km. Suitable for large-scale oceanographic meteorological and climatological applications, such as evaluating or constraining environmental models or case-studies of marine heat wave events.

Adhering to community data standards and names. Includes temperature uncertainty information and auxiliary information about land-sea fraction and sea-ice coverage. To understand the data for your application, read the paper [1] using <a href="www.nature.com/articles/s41597-019-0236-x">www.nature.com/articles/s41597-019-0236-x</a> to cite in any published usage.

The v2.1 record is known to have biases associated with desert dust aerosol and erratic calibration of early-record sensors [1]. Adjustments to reduce these biases and include additional uncertainty in these effects have been developed, as described in [2] and are applied to this data. These adjustments operate on monthly and >5 degree time-space scales.

Data structure
--------------

The data is organised into the following variables

analysed_sst
    Analysed sea surface temperature (Kelvin)
analysis_uncertainty
    Estimated error standard deviation of analysed_sst (kelvin)
sea_ice_fraction
    The estimated fraction of the area covered by sea ice</td><td>-</td>
mask
    Bit mask (bit0:sea,bit1:land:bit2:lake,bit3:ice)

Accessing the data
------------------

An example of how to access the data using python and the xarray library is:

        import xarray as xr
        import matplotlib.pyplot as plt
        import s3fs
        s3 = s3fs.S3FileSystem(anon=True)
        store = s3fs.S3Map(root="s3://surftemp-sst/data/sst.zarr", s3=s3, create=False)
        sst_ds = xr.open_zarr(store)
        # select an area covering the great barrier reef at a particular date
        da = sst_ds["analysed_sst"].sel(time="2016-02-18",
                lat=slice(-15.6,-10.69867),lon=slice(141.4,153.6))
        da.plot(cmap="coolwarm")
        plt.show() # need this if not running in a jupyter notebook

References
----------

[1] Merchant, C.J., Embury, O., Bulgin, C.E., Block, T., Corlett, G.K., Fiedler, E., Good, S.A., Mittaz, J., Rayner, N.A., Berry, D., Eastwood, S., Taylor, M., Tsushima, Y., Waterfall, A., Wilson, R. and Donlon, C. (2019), Satellite-based time-series of sea-surface temperature since 1981 for climate applications. Scientific Data 6, 223, doi:10.1038/s41597-019-0236-x

[2] Merchant, C.J. and Embury, O. (2020) Adjusting for desert-dust-related biases in a climate data record of sea surface temperature. Remote Sensing, 12 (16). 2554. ISSN 2072-4292 doi:10.3390/rs12162554
