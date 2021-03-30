#!/usr/bin/env python3
import argparse
from calendar import month_name
import os

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
import xarray as xr

def compare_L3_Had(hadfiles, l3path, l4path, outPicPath, titlestr='', umax=0.35, xbins=80, xrange=(-1.75, 1.75),
                   vmax=0.8):
    """
    Compare regridded L3U data to L4 and HadSST.

    :param hadfiles:
    :param l3path:
    :param l4path:
    :param outPicPath:
    :param titlestr:
    :param umax:
    :param xbins:
    :param xrange_had:
    :param xrange_l4:
    :param vmax:
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

    # Truncate the dsH record to the same time period as the L3
    dsH = dsH.where(dsH['time'] >= dsl3['time_bnds'][0, 0], drop=True)
    dsH = dsH.where(dsH['time'] <= dsl3['time_bnds'][-1, 1], drop=True)

    dsHu = dsHu.where(dsHu['time'] >= dsl3['time_bnds'][0, 0], drop=True)
    dsHu = dsHu.where(dsHu['time'] <= dsl3['time_bnds'][-1, 1], drop=True)

    # Truncate the L4 record to the same time period as the L3
    dsl4 = dsl4.where(dsl4['time'] >= dsl3['time_bnds'][0, 0], drop=True)
    dsl4 = dsl4.where(dsl4['time'] <= dsl3['time_bnds'][-1, 1], drop=True)

    # Subset the dsH SST record to be on the same months as the L3 and then equalise the times
    dsl3_subset = [(year, month) in zip(dsl3.time.dt.year, dsl3.time.dt.month)
                   for year, month in list(zip(dsH.time.dt.year, dsH.time.dt.month))]
    dsH_sst = dsH.sst[dsl3_subset, :, :]
    dsH_subset = [(year, month) in zip(dsH_sst.time.dt.year, dsH_sst.time.dt.month)
                  for year, month in list(zip(dsl3.time.dt.year, dsl3.time.dt.month))]
    dsl3_sst = dsl3.sst[dsH_subset, :, :]
    dsH_sst['time'] = dsl3_sst['time']
    dsHu_usst = dsHu.usst[dsl3_subset, :, :]
    dsHu_usst['time'] = dsl3_sst['time']

    # Subset the L4 SST, sea fraction and sea ice fraction record to be on the same months as the L3
    subset = [(year, month) in zip(dsl3.time.dt.year, dsl3.time.dt.month)
              and (year, month) in zip(dsH.time.dt.year, dsH.time.dt.month)
              for year, month in list(zip(dsl4.time.dt.year, dsl4.time.dt.month))]
    dsl4_sst = dsl4.sst[subset, :, :]
    dsl4_sf = dsl4.sea_fraction[subset, :, :]
    dsl4_sif = dsl4.sea_ice_fraction[subset, :, :]
    assert (dsl4_sst.time == dsl3_sst.time).all()

    # Mask all the SST datasets to effectively completely ocean cells
    # To avoid small-area low-N effects
    dsl4_sst = dsl4_sst.where(dsl4_sf > 0.99).where(dsl4_sif < 0.01)
    dsH_sst = dsH_sst.where(dsl4_sf > 0.99).where(dsl4_sif < 0.01)
    dsHu_usst = dsHu_usst.where(dsl4_sf > 0.99).where(dsl4_sif < 0.01)
    dsl3_sst = dsl3_sst.where(dsl4_sf > 0.99).where(dsl4_sif < 0.01)

    # Compare the  new regridded l3 and the previous l4
    diff_l4 = dsl3_sst - dsl4_sst
    create_whole_series_plots(diff_l4, titlestr + 'L3 - L4', 'l3-l4', outPicPath, xbins, xrange, vmax=vmax)

    # Compare the l3 and HadSST4
    diff_had = dsl3_sst - dsH_sst
    create_whole_series_plots(diff_had, titlestr + 'L3 - HadSST', 'l3-had', outPicPath, xbins, xrange, vmax=vmax)

    # Then the same things using only the low-uncertainty data in HadSST4
    # Mask the high uncertainty HadSST4
    dsHg_sst = dsH_sst.where(dsHu_usst < umax)
    dsHug_usst = dsHu_usst.where(dsHu_usst < umax)

    # Compare the l3 and HadSST4 best data in HadSST4
    diff_had_best = dsl3_sst - dsHg_sst
    create_whole_series_plots(diff_had_best, titlestr + 'L3 - HadSST Best Data', 'l3-had_best', outPicPath, xbins,
                              xrange, vmax=vmax)

    # Create monthly plots
    for t1, t2, t3 in zip(diff_l4.groupby('time.year'),
                          diff_had.groupby('time.year'),
                          diff_had_best.groupby('time.year')):
        year, yearly_diff_l4 = t1
        year, yearly_diff_had = t2
        year, yearly_diff_had_best = t3
        print('Year:', year)
        n_plots = 4
        assert (12 // n_plots) * n_plots == 12, 'n_plots must be a factor of 12.'
        monthly_diffs_l4 = get_monthly_diffs(yearly_diff_l4)
        monthly_diffs_had = get_monthly_diffs(yearly_diff_had)
        monthly_diffs_had_best = get_monthly_diffs(yearly_diff_had_best)
        for month in range(1, 12, n_plots):
            create_monthly_maps(monthly_diffs_l4, monthly_diffs_had, monthly_diffs_had_best, month, month + n_plots,
                                year, titlestr, vmax, outPicPath)
            create_monthly_hists(monthly_diffs_l4, monthly_diffs_had, monthly_diffs_had_best, month, month + n_plots,
                                 year, titlestr, xbins, xrange, outPicPath)


def create_monthly_maps(monthly_diffs_l4, monthly_diffs_had, monthly_diffs_had_best, start_month, end_month, year,
                        titlestr, vmax, outPicPath):
    """
    Create panel plot of maps for the given year for the months from start month to end month - 1.

    :param monthly_diffs_l4:
    :param monthly_diffs_had:
    :param monthly_diffs_had_best:
    :param start_month:
    :param end_month:
    :param year:
    :param titlestr:
    :param vmax:
    :param outPicPath:
    :return:
    """
    print('Plotting maps...')
    ncols = end_month - start_month
    nrows = 3
    fig, axs = plt.subplots(nrows, ncols, figsize=(6.4 * ncols, 4.8 * nrows), sharex=True, sharey=True,
                            subplot_kw=dict(projection=ccrs.PlateCarree()))
    plot = None
    for month in range(start_month, end_month):
        print('Month:', month)
        col = month - start_month

        diff_l4 = monthly_diffs_l4.get(month, None)
        ax = axs[0, col]
        p = plot_monthly_map(ax, diff_l4, 'L3 - L4 for ' + month_name[month], vmax)
        if p is not None:
            plot = p

        diff_had = monthly_diffs_had.get(month, None)
        ax = axs[1, col]
        p = plot_monthly_map(ax, diff_had, 'L3 - HadSST for ' + month_name[month], vmax)
        if p is not None:
            plot = p

        diff_had_best = monthly_diffs_had_best.get(month, None)
        ax = axs[2, col]
        p = plot_monthly_map(ax, diff_had_best, 'L3 - HadSST Best Data for ' + month_name[month], vmax)
        if p is not None:
            plot = p

    if plot is not None:
        fig.subplots_adjust(right=0.825)
        cbar_ax = fig.add_axes([0.875, 0.15, 0.025, 0.7])
        fig.colorbar(plot, cax=cbar_ax)
        fig.suptitle(titlestr + 'Maps for ' + str(year) + ' from ' + month_name[start_month] + ' to '
                     + month_name[end_month - 1], fontsize='xx-large')
        plt.savefig(os.path.join(outPicPath, 'maps_' + str(year) + '_' + str(start_month) + '_' + str(end_month - 1)
                                 + '.pdf'))
        plt.close()


def create_monthly_hists(monthly_diffs_l4, monthly_diffs_had, monthly_diffs_had_best, start_month, end_month, year,
                         titlestr, xbins, xrange, outPicPath):
    """
    Create panel plot of histograms for the given year for the months from start month to end month - 1.

    :param monthly_diffs_l4:
    :param monthly_diffs_had:
    :param monthly_diffs_had_best:
    :param start_month:
    :param end_month:
    :param year:
    :param titlestr:
    :param xbins:
    :param xrange:
    :param outPicPath:
    :return:
    """
    print('Plotting histograms...')
    ncols = end_month - start_month
    nrows = 3
    fig, axs = plt.subplots(nrows, ncols, figsize=(6.4 * ncols, 4.8 * nrows), sharex=True, sharey=True)
    plot = None
    for month in range(start_month, end_month):
        print('Month:', month)
        col = month - start_month

        diff_l4 = monthly_diffs_l4.get(month, None)
        ax = axs[0, col]
        p = plot_monthly_hist(ax, diff_l4, 'L3 - L4 for ' + month_name[month], xbins, xrange)
        if p is not None:
            plot = p

        diff_had = monthly_diffs_had.get(month, None)
        ax = axs[1, col]
        p = plot_monthly_hist(ax, diff_had, 'L3 - HadSST for ' + month_name[month], xbins, xrange)
        if p is not None:
            plot = p

        diff_had_best = monthly_diffs_had_best.get(month, None)
        ax = axs[2, col]
        p = plot_monthly_hist(ax, diff_had_best, 'L3 - HadSST Best Data for ' + month_name[month], xbins, xrange)
        if p is not None:
            plot = p

    if plot is not None:
        fig.suptitle(titlestr + 'Histograms for ' + str(year) + ' from ' + month_name[start_month] + ' to '
                     + month_name[end_month - 1], fontsize='xx-large')
        plt.savefig(os.path.join(outPicPath, 'hists_' + str(year) + '_' + str(start_month) + '_' + str(end_month - 1)
                                 + '.pdf'))
        plt.close()


def plot_monthly_map(ax, diff, titlestr, vmax):
    """
    Plot an individual map as part of a panel plot of monthly data on the given axis.

    :param ax:
    :param diff:
    :param titlestr:
    :param vmax:
    :return:
    """
    if diff is not None:
        p = diff.plot(vmax=vmax, ax=ax, add_colorbar=False, transform=ccrs.PlateCarree())
        ax.title.set_text(titlestr)
        ax.coastlines()
        return p
    else:
        ax.set_axis_off()
        return None

def plot_monthly_hist(ax, diff, titlestr, xbins, xrange):
    """
    Plot an individual map as part of a panel plot of monthly data on the given axis.

    :param ax:
    :param diff:
    :param titlestr:
    :param xbins:
    :param xrange:
    :return:
    """
    if diff is not None:
        p = diff.plot.hist(bins=xbins, range=xrange, ax=ax)
        coef = add_gaussian(diff, ax, xbins, xrange)
        plt.text(0.75, 0.8, '$a = {0:.3f}$\n$\mu = {1:.3f}$\n$\sigma = {2:.3f}$'.format(*coef), transform=ax.transAxes)
        ax.title.set_text(titlestr)
        return p
    else:
        ax.set_axis_off()
        return None


def get_monthly_diffs(yearly_diff):
    """
    Create a dictionary of the monthly differences from a dataset of the difference for a particular year.
    :param yearly_diff:
    :return:
    """
    monthly_diffs = {int(month): diff for month, diff in yearly_diff.groupby('time.month')}
    return monthly_diffs


def create_whole_series_plots(diff, titlestr, filestr, outPicPath, xbins, xrange, vmin=None, vmax=None):
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

    # Whole series plots, all data
    p = diff.mean(dim='time').plot(vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree(),
                                          subplot_kws=dict(projection=ccrs.PlateCarree()))
    p.axes.coastlines()
    plt.suptitle(titlestr + 'Mean Map')
    plt.savefig(os.path.join(outPicPath, filestr + '_mean_map.pdf'))
    plt.close()

    diff.plot.hist(bins=xbins, range=xrange)
    ax = plt.gca()
    coef = add_gaussian(diff, ax, xbins, xrange)
    plt.text(0.75, 0.8, '$a = {0:.3f}$\n$\mu = {1:.3f}$\n$\sigma = {2:.3f}$'.format(*coef), transform=ax.transAxes)
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
    plt.title(titlestr + 'Mean and Median\nwith +/- 1 Sig, 1st and 99th Percentile')
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
    parser.add_argument('--xbins', help='Number of x bins to use for histograms. Default is 80.',
                        type=int, default=80)
    parser.add_argument('--xmin', help='Minimum of x range for histograms. Default is -1.75.',
                        type=float, default=-1.75)
    parser.add_argument('--xmax', help='Maximum of x range for histograms. Default is 1.75.',
                        type=float, default=1.75)
    parser.add_argument('--vmax', help='Maximum of colorscale for maps. Default is 0.8 degrees C.',
                        type=float, default=0.8)
    args = parser.parse_args()

    compare_L3_Had((args.hadfile, args.hadfileuncert), args.l3path, args.l4path, args.outPicPath,
                   titlestr=args.titlestr, umax=args.umax, xbins=args.xbins, xrange=(args.xmin, args.xmax),
                   vmax=args.vmax)
