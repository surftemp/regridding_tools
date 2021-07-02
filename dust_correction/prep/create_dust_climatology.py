#!/usr/bin/env python
# -*- coding: utf-8 -*-

# script for creating a climatology from dust aerosol data 
#
# 2003/day1.nc
# 2003/day2.nc
# ...
# 2018/day364.nc
# 2018/day365.nc
#
# Note in each input leap year, only 365 days are provided, December 31st is omitted
#
# input variables are:
#    duaod550   Dust Aerosol Optical Depth at 550nm
#    aermssdul  Vertically integrated mass of dust aerosol (9 - 20 um)
#    aermssdum  Vertically integrated mass of dust aerosol (0.55 - 9 um)
#    aermssdus  Vertically integrated mass of dust aerosol (0.03 - 0.55 um)
#
# output climatology variables provide the mean, standard deviation, mean(log), standard deviation(log) and median over the specified window and years
#
# example usage:
#
# python create_dust_climatology.py --start 2003 --stop 2017 --window 2 daily climatology
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

input_pattern = '%(YEAR)s/day%(DOY)s.nc'
output_pattern = '%03d-dust-climatology.nc'

fields = ["aermssdul","aermssdum","aermssdus","duaod550"]

def calcstats(ncs, name):
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
    var = ncs[0].variables[name][:]
    vmin = var.copy()
    vmax = var.copy()
    vsum = var.copy()

    vssq = np.square(var, var)

    n = 1
    for nc in ncs[1:]:
        var = nc.variables[name][:]
        void = np.fmin(var, vmin, vmin)
        void = np.fmax(var, vmax, vmax)
        void = np.add (var, vsum, vsum)
        void = np.add (np.square(var, var), vssq, vssq)
        n = n + 1
    scale = 1.0 / n
    void = np.multiply (scale, vsum, vsum)
    void = np.multiply (scale, vssq, vssq)
    void = np.subtract(vssq, np.square(vsum,var), vssq)
    void = np.sqrt(vssq, vssq)

    return vmin, vmax, vsum, vssq

def calcstats2(ncs, name):
    base_matrices = []
    for nc in ncs:
        base_matrices.append(nc.variables[name])
    m = np.concatenate(base_matrices)
    log_m = np.log(m)
    mean = np.mean(m,axis=0)
    median = np.median(m,axis=0)
    log_mean = np.mean(log_m,axis=0)
    log_sd = np.std(log_m,axis=0)
    return np.expand_dims(mean,axis=0), np.expand_dims(median,axis=0), np.expand_dims(log_mean,axis=0), np.expand_dims(log_sd,axis=0)

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
            continue # can not set _FillValue as an attribute it seems?
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
    >>> opts = {'zlib':True, 'chunksizes'=[1,900,1800]} # Enable compression
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
        mask = data==ncvar._FillValue
    elif hasattr(ncvar, 'missing_value'):
        mask = data==ncvar.missing_value
    else:
        mask = None
    if hasattr(ncvar, 'add_offset') and hasattr(ncvar, 'scale_factor'):
        data = data * ncvar.scale_factor + ncvar.add_offset
    elif hasattr(ncvar, 'scale_factor'):
        data = data * ncvar.scale_factor
    elif hasattr(ncvar, 'add_offset'):
        data += ncvar.add_offset
    if data.dtype.kind =='f' and mask is not None:
        data[mask] = np.nan
    return data


def doytodate(year, doy):
    """
    Return a date for the given year and day-of-year ignoring leap years.
    Day-of-year will wrap - e.g. doytodate(2000,-1) = date(2000, 12, 30)
    """
    day_of_year = 1 + (doy-1) % 365
    return datetime.date.fromordinal(day_of_year).replace(year=year)


def makefilelist(pattern, day_of_year, year_start, year_stop, window=0, path=''):
    files = []
    for year in range(year_start, year_stop+1):
        for doy in range(day_of_year-window, day_of_year+window+1):
            adj_doy = 1 + ((doy-1) % 365)
            fname = pattern % {"YEAR": year, "DOY": adj_doy}
            f = os.path.join(path,fname)
            files.extend(glob.glob(f))
    return files


def makeclimfile(day_of_year, year_start, year_stop, window=0, path='', outpath=''):
    if day_of_year < 1 or day_of_year > 365:
        raise Exception("day_of_year must be >= 1 and <= 365")
    files = makefilelist(input_pattern, day_of_year, year_start, year_stop, window, path)

    expected_file_count = (1+(year_stop - year_start))*(1+2*window)
    if len(files) != expected_file_count:
        raise Exception("Incorrect number of input files located")

    print("Computing climatology for day %d with %d input days (%s)"%(day_of_year,len(files),str(files)))


    ncs   = [ncdf.Dataset(f) for f in files]

    refdate = doytodate(year_stop, day_of_year)
    file0 = os.path.join(path,input_pattern%{"YEAR":year_stop,"DOY":day_of_year})
    try:
        nc0   = ncdf.Dataset(glob.glob(file0)[0])
    except IndexError:
        raise IOError(2, 'No such file or directory', file0)

    outfile = output_pattern %(day_of_year)
    
    newuuid = str(uuid.uuid4())
    newattr = {
            'title':u'Dust Climatology',
            'history':u'created with create_dust_climatology.py',
            'netcdf_version_id':ncdf.__netcdf4libversion__,
            'uuid':newuuid,
            'tracking_id':newuuid,
            'date_created':datetime.datetime.now().strftime(ccitfmt),
            'start_time':refdate.replace(year=year_start).strftime(ccitfmt),
            'time_coverage_start':refdate.replace(year=year_start).strftime(ccitfmt),
            'stop_time':refdate.replace(year=year_stop).strftime(ccitfmt),
            'time_coverage_end':refdate.replace(year=year_stop).strftime(ccitfmt),
            }

    newopts = { 'complevel':9, 'shuffle':True }

    out = ncdf.Dataset(os.path.join(outpath, outfile), 'w', format='NETCDF4_CLASSIC')

    #Copy over existing attributes
    copyattr(nc0, out, newattr)
    copydims(nc0, out)

    skip = fields
    #Copy over existing variables
    for vname in nc0.variables:
        if vname in skip: continue
        var = copyvar(nc0.variables[vname], out)
        var[:] = nc0.variables[vname][:]

    for field in fields:
        vmin, vmax, vsum, vssq = calcstats(ncs, field)

        mean, median, logmean, logsd = calcstats2(ncs,field) # calculated in numpy in 1-step

        print(str(field))

        numpy.testing.assert_allclose(mean,vsum) # sanity check of mean computed in two ways

        svar = nc0.variables[field]

        copyvar(svar, out, 'mean_'+field, newopts, long_name_prefix="Mean of")[:] = vsum
        copyvar(svar, out, 'sdev_'+field,newopts, long_name_prefix="Standard Deviation of")[:] = vssq
        copyvar(svar, out, 'mean_log_' + field, newopts, long_name_prefix="Mean of Log of")[:] = logmean
        copyvar(svar, out, 'sdev_log_' + field, newopts, long_name_prefix="Standard Deviation of Log of")[:] = logsd
        copyvar(svar, out, 'median_' + field, newopts, long_name_prefix="Median of")[:] = median

    out.close()
    nc0.close()

    for nc in ncs:
        nc.close()



if __name__ == "__main__":
    usage = "usage: %prog [options] path_in path_out"
    parser = optparse.OptionParser(version=__version__, usage=usage)
    parser.add_option("--start",  dest="start", type='int', help="Start year (default=1992)")
    parser.add_option("--stop",   dest="stop",  type='int', help="Stop year (default=2010)")
    parser.add_option("--window", dest="win",   type='int', help="Smoothing window (default=2)")
    parser.add_option("--doy",    dest="doy",   type='int', help="Process single day of year (default=0 all)")

    parser.set_defaults(start=1992)
    parser.set_defaults(stop =2010)
    parser.set_defaults(win  =0)
    parser.set_defaults(doy  =0)

    (opts, args) = parser.parse_args()
    
    if len(args) != 2:
        parser.error("incorrect number of arguments")

    if opts.doy:
        makeclimfile(opts.doy, opts.start, opts.stop, opts.win, path=args[0], outpath=args[1])
    else:
        for doy in range(1,366):
            makeclimfile(doy, opts.start, opts.stop, opts.win, path=args[0], outpath=args[1])
