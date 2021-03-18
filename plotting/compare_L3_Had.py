#!/usr/bin/env python3
import argparse
import os

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
import xarray as xr

def compare_L3_Had(hadfiles, l3path, l4path, outPicPath, titlestr='', umax=0.35, xbins=80, xrange=(-2.5, 2.5)):
    """
    Compare regridded L3U data to L4 and HadSST.

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
    titlestr = format_titlestr(titlestr)

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
    diff = dsl3.sst - dsl4.sst
    create_plots(diff, titlestr + 'L3 - L4', 'l3-l4', outPicPath, xbins, xrange, vmax=0.5)

    # Compare the l3 and HadSST4
    diff = dsl3.sst - dsH.sst
    create_plots(diff, titlestr + 'L3 - HadSST', 'l3-had ', outPicPath, xbins, xrange)

    # Then the same things using only the low-uncertainty data in HadSST4
    # Mask the high uncertainty HadSST4
    dsHg = dsH.where(dsHu.usst < umax)
    dsHug = dsHu.where(dsHu.usst < umax)

    # Compare the l3 and HadSST4 best data in HadSST4
    diff = dsl3.sst - dsHg.sst
    create_plots(diff, titlestr + 'L3 - HadSST Best Data', 'l3-had_best', outPicPath, xbins, xrange, vmax=1.0)


def create_plots(diff, titlestr, filestr, outPicPath, xbins, xrange, vmin=None, vmax=None):
    """
    Create maps, histograms and time series of statistics for a particular difference between two datasets.

    :param diff:
    :param titlestr:
    :param filestr:
    :param outPicPath:
    :param xbins:
    :param xrange:
    :param vmin:
    :param vmax:
    :return:
    """
    titlestr = format_titlestr(titlestr)

    # Plot maps and histograms
    yearly_diffs = diff.groupby('time.year')
    for year, yearly_diff in yearly_diffs:
        time = yearly_diff.time
        time_size = time.size

        # Panel plots of maps for each month for each year
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

        plt.suptitle(titlestr + 'Maps for ' + str(year))
        plt.savefig(os.path.join(outPicPath, filestr + '_maps_' + str(year) + '.pdf'))
        plt.close()

        # The same thing for histograms
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
        plt.suptitle(titlestr + 'Histograms for ' + str(year))
        plt.savefig(os.path.join(outPicPath, filestr + '_hists_' + str(year) + '.pdf'))
        plt.close()

    # Whole series plots, all data
    p = diff.mean(dim='time').plot(vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree(),
                                          subplot_kws=dict(projection=ccrs.PlateCarree()))
    p.axes.coastlines()
    plt.suptitle(titlestr + 'Mean Map')
    plt.savefig(os.path.join(outPicPath, filestr + '_mean_map.pdf'))
    plt.close()

    diff.plot.hist(bins=xbins, range=xrange)
    coeff = add_gaussian(diff, plt.gca(), xbins, xrange)
    plt.text(0.75, 0.8, '$a = {0:.3f}$\n$\mu = {1:.3f}$\n$\sigma = {2:.3f}$'.format(*coeff), transform=ax.transAxes)
    plt.suptitle(titlestr + 'Whole Series Histogram')
    plt.savefig(os.path.join(outPicPath, filestr + '_whole_series_hist.pdf'))
    plt.close()

    # Time series of statistics
    dt = diff.mean(dim=('lat', 'lon'))
    sd = diff.std(dim=('lat', 'lon'))
    dtu = dt + sd
    dtl = dt - sd
    med = diff.median(dim=('lat', 'lon'))
    q1 = diff.quantile(0.01, dim=('lat', 'lon'))
    q99 = diff.quantile(0.99, dim=('lat', 'lon'))

    q1.plot(color='orange', label='Q01')
    dtu.plot(color='red', label='+1 sig')
    dt.plot(color='black', label='mean')
    med.plot(color='darkgrey', label='median')
    dtl.plot(color='red', label='-1 sig')
    q99.plot(color='orange', label='Q99')
    plt.title(titlestr + 'Mean and Median with +/- 1 Sig, 1st and 99th Percentile')
    plt.legend()
    plt.savefig(os.path.join(outPicPath, filestr + '_time_series.pdf'))
    plt.close()


def format_titlestr(titlestr):
    """Strip whitespace from either side of the title string and add a single space."""
    titlestr = titlestr.strip() + ' '
    return titlestr


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
    parser = argparse.ArgumentParser(description='Create plots comparing regridded L3U files to regridded L4 and HadSST')
    parser.add_argument('l3path', help='Path to regridded L3U files.')
    parser.add_argument('l4path', help='Path to regridded L4 files.')
    parser.add_argument('outPicPath', help='Path to create output plots in.')
    parser.add_argument('--hadfile', help='Path to HadSST file with actuals and median.',
                        default='/gws/nopw/j04/esacci_sst/validation/HadSST4/HadSST.4.0.0.0_actuals_median.nc')
    parser.add_argument('--hadfileuncert', help='Path to HadSST file with total uncertainty.',
                        default='/gws/nopw/j04/esacci_sst/validation/HadSST4/HadSST.4.0.0.0_total_uncertainty.nc')
    parser.add_argument('--titlestr', help='A title string to prefix to the titles of the plots.', default='')
    parser.add_argument('--umax', help='Maximum uncertainty to include in best HadSST data. Default is 0.35',
                        type=float, default=0.35)
    args = parser.parse_args()

    compare_L3_Had([args.hadfile, args.hadfileuncert], args.l3path, args.l4path, args.outPicPath,
                   titlestr=args.titlestr, umax=args.umax)