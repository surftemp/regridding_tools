#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 11:14:46 2021

@author: nn904972
"""
import os

import numpy as np
import scipy.optimize
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Define paths to data
l3path = "/Users/charles/Data/l3u_regrid_test/CDR2.1"
hadfiles = ["/Users/charles/Data/HadSST4/HadSST.4.0.0.0_actuals_median.nc",
            "//Users/charles/Data/HadSST4/HadSST.4.0.0.0_total_uncertainty.nc"]
l4path = "/Users/charles/Data/monthly_regridded_ssts_5_degrees"
outPicPath = "/Users/charles/Data/HadSSTComps/Plots"

titlestr = 'N07 '


def compare_L3_Had(hadfiles, l3path, l4path, outPicPath, titlestr="", umax=0.35, xbins=80, xrange=(-4, 4)):
    # Read the L3 data regridded to 5 degree monthly
    dsl3 = xr.open_mfdataset(os.path.join(l3path, '*.nc'), combine='by_coords')

    # Read the HadSST4 actuals dataset
    dsH = xr.open_dataset(hadfiles[0])
    dsHu = xr.open_dataset(hadfiles[1])

    # Read in the v2.1 L4 regridded to 5 degree for comparison
    dsl4 = xr.open_mfdataset(os.path.join(l4path, '*.nc'), combine='by_coords')

    # Rename some HadSST variables and make times compatible with L3 labelling
    dsH = dsH.rename_dims({'latitude': 'lat', 'longitude': 'lon'})
    dsH = dsH.rename({'latitude': 'lat', 'longitude': 'lon'})
    dsH = dsH.rename({'latitude_bnds': 'lat_bnds', 'longitude_bnds': 'lon_bnds', 'tos': 'sst'})
    dsH = dsH.where(dsH.sst < 9e36)
    dsH['sst'] += 273.15

    dsHu = dsHu.rename_dims({'latitude': 'lat', 'longitude': 'lon'})
    dsHu = dsHu.rename({'latitude': 'lat', 'longitude': 'lon'})
    dsHu = dsHu.rename({'latitude_bnds': 'lat_bnds', 'longitude_bnds': 'lon_bnds', 'tos_unc': 'usst'})
    dsHu = dsHu.where(dsHu.usst < 9e36)

    # Truncate the dsH record to the same months as the L3, then equalise the times
    dsH = dsH.where(dsH['time'] >= dsl3['time_bnds'][0, 0], drop=True)
    dsH = dsH.where(dsH['time'] <= dsl3['time_bnds'][-1, 1], drop=True)
    dsH['time'] = dsl3['time']

    dsHu = dsHu.where(dsHu['time'] >= dsl3['time_bnds'][0, 0], drop=True)
    dsHu = dsHu.where(dsHu['time'] <= dsl3['time_bnds'][-1, 1], drop=True)
    dsHu['time'] = dsl3['time']

    # Truncate the L4 record to the same months as the L3
    dsl4 = dsl4.where(dsl4['time'] >= dsl3['time_bnds'][0, 0], drop=True)
    dsl4 = dsl4.where(dsl4['time'] <= dsl3['time_bnds'][-1, 1], drop=True)

    # Mask all the datasets to effectively completely ocean cells
    # To avoid small-area low-N effects
    dsl4 = dsl4.where(dsl4.sea_fraction > 0.99).where(dsl4.sea_ice_fraction < 0.01)
    dsH = dsH.where(dsl4.sea_fraction > 0.99).where(dsl4.sea_ice_fraction < 0.01)
    dsHu = dsHu.where(dsl4.sea_fraction > 0.99).where(dsl4.sea_ice_fraction < 0.01)
    dsl3 = dsl3.where(dsl4.sea_fraction > 0.99).where(dsl4.sea_ice_fraction < 0.01)

    # Compare the  new regridded l3 and the previous l4
    # 0 --> a plot per month gather in annual plots
    # Histograms similarly, but on a common binsize/scale adn vertical axis for a given year
    (dsl3.sst - dsl4.sst)[0, :, :].plot(vmax=0.5)
    plt.show()
    (dsl3.sst - dsl4.sst)[0, :, :].plot.hist(bins=xbins, range=xrange)
    add_gaussian((dsl3.sst - dsl4.sst)[0, :, :], plt.gca(), xbins, xrange)

    plt.show()

    # Compare the l3 and HadSST4
    (dsl3.sst - dsH.sst)[0, :, :].plot()
    plt.show()
    (dsl3.sst - dsH.sst)[0, :, :].plot.hist(bins=xbins, range=xrange)
    add_gaussian((dsl3.sst - dsH.sst)[0, :, :], plt.gca(), xbins, xrange)
    plt.show()
    # ((dsl3.sst - dsH.sst)/dsHu.usst)[0,:,:].plot(vmax = 3)
    # ((dsl3.sst - dsH.sst)/dsHu.usst)[0,:,:].plot.hist()

    # Whole series plots, all data
    # Compare the  new regridded l3 and the previous l4
    (dsl3.sst - dsl4.sst).mean(dim='time').plot(vmax=1.0)
    plt.show()
    (dsl3.sst - dsl4.sst).plot.hist(bins=xbins, range=xrange)
    add_gaussian(dsl3.sst - dsl4.sst, plt.gca(), xbins, xrange)
    plt.show()
    # Compare the  new regridded l3 and the previous Had4
    (dsl3.sst - dsH.sst).mean(dim='time').plot(vmax=1.0)
    plt.show()
    (dsl3.sst - dsH.sst).plot.hist(bins=xbins, range=xrange)
    add_gaussian(dsl3.sst - dsH.sst, plt.gca(), xbins, xrange)
    plt.show()

    # Timeseries of statistics
    # First, using all common data
    dt = (dsl3.sst - dsl4.sst).mean(dim=('lat', 'lon'))
    sd = (dsl3.sst - dsl4.sst).std(dim=('lat', 'lon'))
    dtu = dt + sd
    dtl = dt - sd
    dt.plot(color='black')
    dtu.plot(color='red')
    dtl.plot(color='red')
    plt.title(titlestr + "L3 - L4, global mean +/- 1 sig")
    plt.show()

    dt = (dsl3.sst - dsH.sst).mean(dim=('lat', 'lon'))
    sd = (dsl3.sst - dsH.sst).std(dim=('lat', 'lon'))
    dtu = dt + sd
    dtl = dt - sd
    dt.plot(color='black')
    dtu.plot(color='red')
    dtl.plot(color='red')
    plt.title(titlestr + "L3 - Had4, all-data mean +/- 1 sig")
    plt.show()

    # Then the same things using only the low-uncertainty data in HadSST4
    # Mask the high uncertainty HadSST4
    dsHg = dsH.where(dsHu.usst < umax)
    dsHug = dsHu.where(dsHu.usst < umax)

    # Compare the l3 and HadSST4 best data in HadSST4
    (dsl3.sst - dsHg.sst)[0, :, :].plot(vmax=1)
    plt.show()
    (dsl3.sst - dsHg.sst)[0, :, :].plot.hist(bins=xbins, range=xrange)
    add_gaussian((dsl3.sst - dsHg.sst)[0, :, :], plt.gca(), xbins, xrange)
    plt.show()
    # ((dsl3.sst - dsHg.sst)/dsHug.usst)[0,:,:].plot(vmax = 3)
    # ((dsl3.sst - dsHg.sst)/dsHug.usst)[0,:,:].plot.hist()

    # Whole series plots, best Had data
    # Compare the  new regridded l3 and the previous Had4
    (dsl3.sst - dsHg.sst).mean(dim='time').plot(vmax=1.0)
    plt.show()
    (dsl3.sst - dsHg.sst).plot.hist(bins=xbins, range=xrange)
    add_gaussian(dsl3.sst - dsHg.sst, plt.gca(), xbins, xrange)
    plt.show()

    # using best data
    # add 1% and 99% quantiles, and median
    dt = (dsl3.sst - dsHg.sst).mean(dim=('lat', 'lon'))
    sd = (dsl3.sst - dsHg.sst).std(dim=('lat', 'lon'))
    dtu = dt + sd
    dtl = dt - sd
    q1 = (dsl3.sst - dsHg.sst).quantile(0.01, dim=('lat', 'lon'))
    q99 = (dsl3.sst - dsHg.sst).quantile(0.99, dim=('lat', 'lon'))

    dt.plot(color='black')
    dtu.plot(color='red')
    dtl.plot(color='red')
    q1.plot(color='orange')
    q99.plot(color='orange')
    plt.title(titlestr + "L3 - Had4, best in-situ mean +/- 1 sig")
    plt.legend(['Q99', '+1 sig', 'mean', '-1 sig', 'Q01'])
    plt.show()


def add_gaussian(data, axs, xbins, xrange):
    """
    Add a simple Gaussian to a plot.

    :param data:
        The data to fit the Gaussian to.
    :param axs:
        The axes of the plot to add the Gaussian to.
    :return:
    """
    def _gauss(x, a, mu, sigma):
        """Gaussian function"""
        return a * np.exp(-(x - mu) ** 2 / (2. * sigma ** 2))

    # Change as required - these are just the defaults from the PVIR plots
    hist, edge = np.histogram(data, xbins, xrange)
    xval = (edge[:-1] + edge[1:]) / 2
    # Find the best fit Gaussian to the histogram data
    # Alternatively we could use robust statistics. e.g. median and robust standard deviation
    try:
        p0 = [hist.max(), 0.0, 1.0]
        coef, vmat = scipy.optimize.curve_fit(_gauss, xval.astype(float), hist.astype(float), p0=p0)
        axs.plot(xval, _gauss(xval, *coef), color='black', linestyle='--')
    except RuntimeError:
        coef = None

    return coef


if __name__ == '__main__':
    compare_L3_Had(hadfiles, l3path, l4path, outPicPath, titlestr=titlestr)