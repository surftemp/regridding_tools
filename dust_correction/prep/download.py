#!/usr/bin/env python

from ecmwfapi import ECMWFDataServer

server = ECMWFDataServer()

for year in range(2003,2019):
    try:
        server.retrieve({
            "class": "mc",
            "dataset": "cams_reanalysis",
            "date": "%s-01-01/to/%s-12-31"%(year,year),
            "expver": "eac4",
            "levtype": "sfc",
            "param": "209.210/43.215/44.215/45.215",
            "stream": "oper",
            "time": "12:00:00",
            "type": "an",
            "grid": "0.5/0.5",
            "format": "netcdf",
            "target": "output%s_0p5.nc" % (year)
        })
    except Exception as ex:
        print(str(ex))

