#!/usr/bin/env python
# -*- coding: utf-8 -*-

# script for creating a climatology from dust corrected SST data
#
# Note in each input leap year, only 365 days are provided, December 31st is omitted
#
# input variables are:
#    analysed_sst
#
# output climatology variables provide the mean
#
# example usage:
#
# python create_climatology.py --start 2003 --stop 2017 --window 2 daily climatology
#
# based on:
#
# https://github.com/surftemp/common/blob/master/scripts/cci/create_climatology.py

__author__ = "Owen Embury <o.embury@reading.ac.uk>, Claire Macintosh <c.r.macintosh@reading.ac.uk>"
__license__ = "GNU GPL v3"
__version__ = "2.0"

import glob
import os.path
import datetime
import optparse
import uuid
import numpy as np
import netCDF4 as ncdf
import numpy.testing

ccitfmt = '%Y%m%dT%H%M%SZ'

input_pattern = '%(Y)d/%(m)02d/%(d)02d/%(Y)d%(m)02d%(d)02d120000-ESACCI-L4_GHRSST-SSTdepth-OSTIA-GLOB_CDR2.1-v02.0-fv01.0.nc'
output_pattern = 'D%03d-ESACCI-L4_GHRSST-SSTdepth-OSTIA-GLOB_CDR2.1-v02.0-fv01.0.nc'

fields = ["analysed_sst"]

def calcstats(files, name):
    """
    Calculate variable statistics (min,max,mean,sdev) across a list of netCDF objects.

    Parameters
    ----------
    ncs : list of netCDF4.Dataset
        List of netCDF.Datasets to read
    name : string
        netCDF variable name

    Returns
    -------
    Tuple of (minimum, maximum, mean, standard deviation)
    Each component has the same shape and type as the netCDF variable processed.
    """

    print("Reading %s"%(files[0]))
    nc0 = ncdf.Dataset(files[0])
    var = nc0.variables[name][:]
    vsum = var.copy()


    n = 1
    for file in files[1:]:
        print("Reading %s"%(file))
        nc = ncdf.Dataset(file)
        var = nc.variables[name][:]
        void = np.add(var, vsum, vsum)
        n = n + 1
        nc.close()
    scale = 1.0 / n
    void = np.multiply(scale, vsum, vsum)
    nc0.close()
    return vsum


def calcstats2(ncs, name):
    base_matrices = []
    for nc in ncs:
        base_matrices.append(nc.variables[name])
    m = np.concatenate(base_matrices)
    log_m = np.log(m)
    mean = np.mean(m, axis=0)
    median = np.median(m, axis=0)
    log_mean = np.mean(log_m, axis=0)
    log_sd = np.std(log_m, axis=0)
    return np.expand_dims(mean, axis=0), np.expand_dims(median, axis=0), np.expand_dims(log_mean,
                                                                                        axis=0), np.expand_dims(log_sd,
                                                                                                                axis=0)


def copyattr(src, dst, newattr={}, long_name_prefix=""):
    """
    Copy netCDF attributes from src variable/dataset to dst variable/dataset

    Parameters
    ----------
    src : netCDF4.Dataset or netCDF4.Variable
        Source variable or dataset
    dst : netCDF4.Dataset or netCDF4.Variable
        Desination variable or dataset
    newattr : dict, optional
        New attributes to add to destination variable/dataset. Will override
        attributes in source.

    Examples
    --------
    >>> nc  = ncdf.Dataset('original.nc')
    >>> out = ncdf.Dataset('copy.nc', 'w')
    >>> attr= {'history':'Copied file using Python'}
    >>> copyattr(nc, out, newattr=attr)  # Copy global attributes to new file.
    """
    new = newattr.copy()
    for attr in src.ncattrs():
        if attr == "_FillValue":
            continue  # can not set _FillValue as an attribute it seems?
        try:
            dst.setncattr(attr, new[attr])
            del new[attr]
        except KeyError:
            if attr == "long_name":
                new_long_name = long_name_prefix + " " + src.getncattr(attr)
                dst.setncattr(attr, new_long_name)
            else:
                dst.setncattr(attr, src.getncattr(attr))

    for attr in new:
        dst.setncattr(attr, new[attr])


def copydims(src, dst):
    """
    Copy netCDF dimensions from src dataset to dst dataset

    Parameters
    ----------
    src : netCDF4.Dataset
        Source dataset
    dst : netCDF4.Dataset
        Desination dataset
    """
    for dim in src.dimensions:
        if src.dimensions[dim].isunlimited():
            dst.createDimension(dim, None)
        else:
            dst.createDimension(dim, len(src.dimensions[dim]))


def copyvar(src, dst, newname=None, newopts={}, newattr={}, long_name_prefix=""):
    """
    Copy netCDF variable definition from src to dst

    Parameters
    ----------
    src : netCDF4.Variable
        Source variable
    dst : netCDF4.Dataset
        Desination dataset
    newname : string, optional
        Specify a different name for the destination variable
    newopts: dict, optional
        Specify different variable creation options (compression etc.)
        Will be passed to createVariable as keyword arguments
    newattr : dict, optional
        New attributes to add to destination variable. Will override
        attributes in source.

    Returns
    -------
    New netCDF variable in dst

    Examples
    --------
    >>> nc  = ncdf.Dataset('original.nc')
    >>> out = ncdf.Dataset('copy.nc', 'w')
    >>> opts = {'zlib':True, 'chunksizes':[1,900,1800]} # Enable compression
    >>> newvar = copyncvar(nc.variables['sst'], out, newopts=opts)
    >>> newvar[:] = nc.variables['sst'][:]              # Copy data
    """
    opts = src.filters()
    if opts is None:
        opts = {}
    if src.chunking() and src.chunking() != 'contiguous':
        opts['chunksizes'] = src.chunking()
    if hasattr(src, '_FillValue'):
        opts['fill_value'] = src._FillValue
    opts.update(newopts)
    if newname:
        name = newname
    else:
        name = src._name
    var = dst.createVariable(name, src.dtype, src.dimensions, **opts)
    copyattr(src, var, newattr, long_name_prefix)
    return var


def read(ncvar):
    """
    Read a netcdf variable to a numpy array using NaNs rather than as a masked
    array.
    """
    if not ncvar.maskandscale:
        return ncvar[:]
    ncvar.set_auto_maskandscale(False)
    data = ncvar[:]
    ncvar.set_auto_maskandscale(True)
    if hasattr(ncvar, '_FillValue'):
        mask = data == ncvar._FillValue
    elif hasattr(ncvar, 'missing_value'):
        mask = data == ncvar.missing_value
    else:
        mask = None
    if hasattr(ncvar, 'add_offset') and hasattr(ncvar, 'scale_factor'):
        data = data * ncvar.scale_factor + ncvar.add_offset
    elif hasattr(ncvar, 'scale_factor'):
        data = data * ncvar.scale_factor
    elif hasattr(ncvar, 'add_offset'):
        data += ncvar.add_offset
    if data.dtype.kind == 'f' and mask is not None:
        data[mask] = np.nan
    return data


def doytodate(year, doy):
    """
    Return a date for the given year and day-of-year ignoring leap years.
    Day-of-year will wrap - e.g. doytodate(2000,-1) = date(2000, 12, 30)
    """
    day_of_year = 1 + (doy - 1) % 365
    return datetime.date.fromordinal(day_of_year).replace(year=year)


def makefilelist(pattern, day_of_year, year_start, year_stop, window=0, path=''):
    files = []
    for year in range(year_start, year_stop + 1):
        base_dt = doytodate(year,day_of_year)
        for woff in range(-window, window + 1):
            dt = base_dt + datetime.timedelta(days=woff)
            fname = pattern%{"Y": dt.year, "m": dt.month, "d": dt.day}
            f = os.path.join(path, fname)
            files.extend(glob.glob(f))
    return files

def makefilelistX(pattern, day_of_year, year_start, year_stop, window=0, path=''):
    files = []
    for year in range(year_start, year_stop + 1):
        for doy in range(day_of_year - window, day_of_year + window + 1):
            adj_doy = 1 + ((doy - 1) % 365)
            dt = doytodate(year,doy)
            fname = pattern%{"Y": year, "m": dt.month, "d": dt.day}
            f = os.path.join(path, fname)
            files.extend(glob.glob(f))
    return files

def makeclimfile(day_of_year, year_start, year_stop, window=0, path='', outpath=''):
    if day_of_year < 1 or day_of_year > 365:
        raise Exception("day_of_year must be >= 1 and <= 365")
    files = makefilelist(input_pattern, day_of_year, year_start, year_stop, window, path)

    expected_file_count = (1 + (year_stop - year_start)) * (1 + 2 * window)
    if len(files) != expected_file_count:
        raise Exception("Incorrect number of input files located, found %d expected %d"%(len(files),expected_file_count))

    # print("Computing climatology for day %d with %d input days (%s)" % (day_of_year, len(files), str(files)))

    # ncs = [ncdf.Dataset(f) for f in files]

    refdate = doytodate(year_stop, day_of_year)
    file0 = os.path.join(path, input_pattern % {"Y": year_stop, "m": refdate.month, "d":refdate.day})
    try:
        nc0 = ncdf.Dataset(glob.glob(file0)[0])
    except IndexError:
        raise IOError(2, 'No such file or directory', file0)

    outfile = output_pattern % (day_of_year)

    newuuid = str(uuid.uuid4())
    newattr = {
        'title': u'Dust Corrected Climatology',
        'history': u'created with climatology.py',
        'netcdf_version_id': ncdf.__netcdf4libversion__,
        'uuid': newuuid,
        'tracking_id': newuuid,
        'date_created': datetime.datetime.now().strftime(ccitfmt),
        'start_time': refdate.replace(year=year_start).strftime(ccitfmt),
        'time_coverage_start': refdate.replace(year=year_start).strftime(ccitfmt),
        'stop_time': refdate.replace(year=year_stop).strftime(ccitfmt),
        'time_coverage_end': refdate.replace(year=year_stop).strftime(ccitfmt),
    }

    newopts = {'complevel': 9, 'shuffle': True}

    out = ncdf.Dataset(os.path.join(outpath, outfile), 'w', format='NETCDF4_CLASSIC')

    # Copy over existing attributes
    copyattr(nc0, out, newattr)
    copydims(nc0, out)

    skip = fields + ["analysed_sst_uncertainty", "sea_ice_fraction"]
    # Copy over existing variables
    for vname in nc0.variables:
        if vname in skip: continue
        var = copyvar(nc0.variables[vname], out)
        var[:] = nc0.variables[vname][:]

    for field in fields:
        vsum = calcstats(files, field)

        # mean, median, logmean, logsd = calcstats2(ncs, field)  # calculated in numpy in 1-step

        # print(str(field))

        # numpy.testing.assert_allclose(mean, vsum)  # sanity check of mean computed in two ways

        svar = nc0.variables[field]

        copyvar(svar, out, '' + field, newopts, long_name_prefix="")[:] = vsum
        # copyvar(svar, out, 'sdev_' + field, newopts, long_name_prefix="Standard Deviation of")[:] = vssq
        # copyvar(svar, out, 'mean_log_' + field, newopts, long_name_prefix="Mean of Log of")[:] = logmean
        # copyvar(svar, out, 'sdev_log_' + field, newopts, long_name_prefix="Standard Deviation of Log of")[:] = logsd
        # copyvar(svar, out, 'median_' + field, newopts, long_name_prefix="Median of")[:] = median

    out.close()


if __name__ == "__main__":
    usage = "usage: %prog [options]"
    parser = optparse.OptionParser(version=__version__, usage=usage)
    parser.add_option("--start", dest="start", type='int', help="Start year (default=1982)")
    parser.add_option("--stop", dest="stop", type='int', help="Stop year (default=2010)")
    parser.add_option("--window", dest="win", type='int', help="Smoothing window (default=2)")
    parser.add_option("--doy", dest="doy", type='int', help="Process single day of year (default=0 all)")

    parser.set_defaults(start=1982)
    parser.set_defaults(stop=2010)
    parser.set_defaults(win=2)
    parser.set_defaults(doy=0)

    (opts, args) = parser.parse_args()

    inpath = "/Users/cv922550/Data/sst/data/CDR_v2/Analysis/L4/v2.1" if len(args) < 1 else args[0]
    outpath = "/tmp/climatology" if len(args) < 2 else args[1]

    if opts.doy:
        makeclimfile(opts.doy, opts.start, opts.stop, opts.win, path=inpath, outpath=outpath)
    else:
        for doy in range(1, 366):
            makeclimfile(doy, opts.start, opts.stop, opts.win, path=inpath, outpath=outpath)
