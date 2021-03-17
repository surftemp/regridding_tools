#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 11:14:46 2021

@author: nn904972
"""
from calendar import month_name
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

titlestr = 'N07'


def compare_L3_Had(hadfiles, l3path, l4path, outPicPath, titlestr="", umax=0.35, xbins=80, xrange=(-4, 4)):
    """
    Compare regridded L3U data to HadSST.

    :param hadfiles:
    :param l3path:
    :param l4path:
    :param outPicPath:
    :param titlestr:
    :param umax:
    :param xbins:
    :param xrange:
    :return:
    """
    # Strip whitespace from either side of the title string.
    titlestr = titlestr.strip()

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
    diff = dsl3.sst - dsl4.sst
    plot_maps_and_hists(diff, titlestr + ' L3 - L4', 'l3-l4', outPicPath,
                        vmax=0.5, xbins=xbins, xrange=xrange)

    # Compare the l3 and HadSST4
    diff = dsl3.sst - dsH.sst
    plot_maps_and_hists(diff, titlestr + ' L3 - HadSST', 'l3-had', outPicPath,
                        xbins=xbins, xrange=xrange)

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
    diff = dsl3.sst - dsHg.sst
    plot_maps_and_hists(diff, titlestr + ' L3 - HadSST Best Data', 'l3-had_best', outPicPath,
                        vmax=1.0, xbins=xbins, xrange=xrange)

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
    med = (dsl3.sst - dsHg.sst).median(dim=('lat', 'lon'))
    q1 = (dsl3.sst - dsHg.sst).quantile(0.01, dim=('lat', 'lon'))
    q99 = (dsl3.sst - dsHg.sst).quantile(0.99, dim=('lat', 'lon'))

    q1.plot(color='orange', label='Q01')
    dtu.plot(color='red', label='+1 sig')
    dt.plot(color='black', label='mean')
    med.plot(color='darkgrey', label='median')
    dtl.plot(color='red', label='-1 sig')
    q99.plot(color='orange', label='Q99')
    plt.title(titlestr + "L3 - Had4, best in-situ mean +/- 1 sig")
    plt.legend()
    plt.show()


def plot_maps_and_hists(diff, titlestr, filestr, outPicPath, vmin=None, vmax=None, xbins=80, xrange=(-4, 4)):
    """
    Plot maps and histograms.

    :param diff:
    :param titlestr:
    :param filestr:
    :param outPicPath:
    :param vmin:
    :param vmax:
    :param xbins:
    :param xrange:
    :return:
    """
    yearly_diffs = diff.groupby('time.year')
    for year, yearly_diff in yearly_diffs:
        time = yearly_diff.time
        time_size = time.size

        if time_size == 1:
            p = yearly_diff.plot(vmin=vmin, vmax=vmax,
                                 transform=ccrs.PlateCarree(), subplot_kws=dict(projection=ccrs.PlateCarree()))
            p.axes.coastlines()
        else:
            p = yearly_diff.plot(x='lon', y='lat', col='time', col_wrap=None if time_size <= 4 else 4,
                                 vmin=vmin, vmax=vmax,
                                 aspect=yearly_diff['lon'].size / yearly_diff['lat'].size,
                                 transform=ccrs.PlateCarree(), subplot_kws=dict(projection=ccrs.PlateCarree()))
            for ax in p.axes.flat:
                ax.coastlines()

        plt.suptitle(titlestr + ' Maps for ' + str(year))
        plt.savefig(os.path.join(outPicPath, str(year) + '_' + filestr + '_maps.pdf'))
        plt.close()

        nrows = (time_size - 1) // 4 + 1
        ncols = time_size if time_size <= 4 else 4
        fig, axs = plt.subplots(nrows, ncols, figsize=(6.4 * ncols, 4.8 * nrows),
                                sharex=(time_size > ncols), sharey=(time_size > 1))
        if time_size == 1:
            axs = [axs]
        else:
            axs = axs.flat
        for i, ax in enumerate(axs):
            if i < time_size:
                yearly_diff[i, :, :].plot.hist(bins=xbins, range=xrange, ax=ax)
                ax.title.set_text('time = ' + time[i].dt.strftime('%Y-%m-%dT%X').values)
                coeff = add_gaussian(yearly_diff[i, :, :], ax, xbins, xrange)
                plt.text(0.75, 0.8, '$a = {0:.3f}$\n$\mu = {1:.3f}$\n$\sigma = {2:.3f}$'.format(*coeff),
                         transform=ax.transAxes)
            else:
                ax.set_axis_off()
        plt.suptitle(titlestr + ' Histograms for ' + str(year))
        plt.savefig(os.path.join(outPicPath, str(year) + '_' + filestr + '_hists.pdf'))
        plt.close()


def add_gaussian(data, ax, xbins, xrange):
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

    hist, edge = np.histogram(data, xbins, xrange)
    xval = (edge[:-1] + edge[1:]) / 2
    # Find the best fit Gaussian to the histogram data
    # Alternatively we could use robust statistics. e.g. median and robust standard deviation
    try:
        p0 = [hist.max(), 0.0, 1.0]
        coef, vmat = scipy.optimize.curve_fit(_gauss, xval.astype(float), hist.astype(float), p0=p0)
        ax.plot(xval, _gauss(xval, *coef), color='black', linestyle='--')
    except RuntimeError:
        coef = None

    return coef


if __name__ == '__main__':
    compare_L3_Had(hadfiles, l3path, l4path, outPicPath, titlestr=titlestr)