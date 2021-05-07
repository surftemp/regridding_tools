import cf
import os.path
import datetime
import json
import numpy as np
import copy
from uuid import uuid4
from utils import TimeSeriesUtils

FIELD_PROPERTIES = json.loads(open(os.path.join(os.path.split(__file__)[0],"field_properties_template.json"),"r").read())
GLOBAL_PROPERTIES = json.loads(open(os.path.join(os.path.split(__file__)[0],"global_properties_template.json"),"r").read())


class NetCDF4Formatter(object):
    """
    Create a formatter for writing either timeseries or region data to a netcdf4 output file
    """

    def __init__(self,path="",time_resolution="daily",f_max=1.0,spatial_lambda=1.0,tau=3,output_anomaly=False,output_sea_ice=False,lon_min=0,lat_min=0,lon_max=5,lat_max=5,comment="",header=""):
        """
        Construct the netcdf4 formatter using options

        :param path: the output path
        :param time_resolution: the time resolution
        :param f_max: sea ice threshold
        :param spatial_lambda: parameter for uncertainty calculation
        :param tau: parameter for uncertainty calculation
        :param output_anomaly: output anomaly rather than absolute SST
        :param output_sea_ice: include sea ice fraction in output
        :param lon_min:  minimum longitude of box, must be aligned on 0.05 degree boundary
        :param lat_min:  minimum latitude of box, must be aligned on 0.05 degree boundary
        :param lon_max:  maximum longitude of box, must be aligned on 0.05 degree boundary
        :param lat_max:  maximum latitude of box, must be aligned on 0.05 degree boundary
        :param comment: an extra string to add to the output file metadata comment
        :param header: main text to use in the output file metadata comment
        """
        self.output_path = path
        self.time_resolution = time_resolution
        self.output_anomaly = output_anomaly
        self.output_sea_ice = output_sea_ice
        self.lon_min = lon_min
        self.lon_max = lon_max
        self.lat_min = lat_min
        self.lat_max = lat_max
        self.comment = comment
        self.f_max = f_max
        self.spatial_lambda = spatial_lambda
        self.tau = tau
        self.datapoints = []
        self.min_dt = None
        self.max_dt = None
        self.uuid = str(uuid4())
        self.header = header

    def write(self,start_dt,mid_dt,end_dt,sst_or_anomaly,uncertainty,sea_ice_fraction):
        """
        Write an entry to the output file covering a time period
        :param start_dt: start date of the period
        :param mid_dt: mid date of the period
        :param end_dt: end date of the period
        :param sst_or_anomaly: the absolute SST or anomaly SST value
        :param uncertainty: the estimated uncertainty standard error
        :param sea_ice_fraction: the sea ice fraction
        """
        self.datapoints.append((start_dt,mid_dt,end_dt,sst_or_anomaly,uncertainty,sea_ice_fraction))

    def createField(self,values,name,long_name,units,scalar_int=False):
        """
        Create a cf.Field object with a time axis and metadata
        :param values: the array of values to form the field's data
        :param name: the field short name
        :param long_name: the field long name
        :param units: the field's measurement units
        :param scalar_int: True iff the field contains integers
        :return: a cf.Field object
        """

        properties = copy.deepcopy(FIELD_PROPERTIES)
        field = cf.Field(properties=properties)

        for key in GLOBAL_PROPERTIES:
            value = GLOBAL_PROPERTIES[key]
            field.nc_set_global_attribute(key,value)

        field.nc_set_global_attribute('f_max', self.f_max)
        field.nc_set_global_attribute('lambda', self.spatial_lambda)
        field.nc_set_global_attribute('tau', np.int32(self.tau))

        field.nc_set_variable(name)
        field.set_property("standard_name",long_name)
        field.long_name = long_name

        # Note - domain axes are constructed based on running .creation_commands() on a regridded cf.Field

        # domain axis - time
        c = cf.DomainAxis(size=len(values))
        c.nc_set_dimension('time')
        field.set_construct(c, key='domainaxis0')

        array_width = 0
        array_height = 0
        if not scalar_int:
            array_height = values[0].shape[0]
            array_width = values[0].shape[1]

            # domain_axis - latitude
            c = cf.DomainAxis(size=array_height)
            c.nc_set_dimension('lat')
            field.set_construct(c, key='domainaxis1')

            # domain_axis - longitude
            c = cf.DomainAxis(size=array_width)
            c.nc_set_dimension('lon')
            field.set_construct(c, key='domainaxis2')

        # Set the field data from the passed values
        if not scalar_int:
            values = list(map(lambda arr: arr.reshape(-1,array_height,array_width),values))
            arr = np.concatenate(values,axis=0) # values is a list of np array of the correct shape
        else:
            arr = np.array(values).reshape((-1)) # values is a list of scalar values, convert to an array
        m_arr = np.ma.masked_invalid(arr)
        data = cf.Data(m_arr, units=units, dtype='i4' if scalar_int else "f4")

        if not scalar_int:
            field.set_data(data, axes=('domainaxis0', 'domainaxis1', 'domainaxis2'))
        else:
            field.set_data(data, axes=('domainaxis0'))

        # Note - dimension coordinates below are constructed based on running .creation_commands() on a regridded cf.Field

        # dimension_coordinate - Time
        c = cf.DimensionCoordinate()
        c.set_properties({'long_name': 'reference time of sst field', 'standard_name': 'time', 'axis': 'T',
                          'units': 'seconds since 1981-01-01 00:00:00', 'calendar': 'gregorian', 'comment': ''})
        c.nc_set_variable('time')
        data = cf.Data([mdt for (_,mdt,_) in self.times], units='seconds since 1981-01-01 00:00:00', calendar='gregorian', dtype='i4')
        c.set_data(data)
        b = cf.Bounds()
        b.set_properties({'long_name': 'Time cell boundaries',
                          'comment': 'Contains the start and end times for the time period the data represent.',
                          'units': 'seconds since 1981-01-01 00:00:00', 'calendar': 'gregorian'})
        b.nc_set_variable('time_bnds')
        data = cf.Data([[sdt, edt] for (sdt,_,edt) in self.times], units='seconds since 1981-01-01 00:00:00', calendar='gregorian',
                       dtype='i4')
        b.set_data(data)
        c.set_bounds(b)
        field.set_construct(c, axes=('domainaxis0',), key='dimensioncoordinate0', copy=False)

        if not scalar_int:
            # dimension_coordinate - Latitude
            c = cf.DimensionCoordinate()
            c.set_properties(
                {'standard_name': 'latitude', 'long_name': 'latitude', 'units': 'degrees_north', 'valid_min': -90.0,
                 'valid_max': 90.0, 'axis': 'Y', 'comment': ' Latitude geographical coordinates,WGS84 projection'})
            c.nc_set_variable('lat')
            cell_height = (self.lat_max-self.lat_min)/array_height
            data = cf.Data([np.float32(self.lat_min + i*cell_height + cell_height/2.0) for i in range(0,array_height)], units='degrees_north', dtype='f4')
            c.set_data(data)
            b = cf.Bounds()
            b.set_properties({'units': 'degrees_north'})
            b.nc_set_variable('lat_bnds')
            data = cf.Data([[np.float32(self.lat_min+i*cell_height), np.float32(self.lat_min+(i+1)*cell_height)] for i in range(0,array_height)], units='degrees_north', dtype='f4')
            b.set_data(data)
            c.set_bounds(b)
            field.set_construct(c, axes=('domainaxis1',), key='dimensioncoordinate1', copy=False)

            # dimension_coordinate - Longitude
            c = cf.DimensionCoordinate()
            c.set_properties(
                {'standard_name': 'longitude', 'long_name': 'longitude', 'units': 'degrees_east', 'valid_min': -180.0,
                 'valid_max': 180.0, 'axis': 'X', 'comment': ' Longitude geographical coordinates,WGS84 projection'})
            c.nc_set_variable('lon')
            cell_width = (self.lon_max - self.lon_min) / array_width
            data = cf.Data([np.float32(self.lon_min + i*cell_width + cell_width/2.0) for i in range(0,array_width)], units='degrees_east', dtype='f4')
            c.set_data(data)
            b = cf.Bounds()
            b.set_properties({'units': 'degrees_east'})
            b.nc_set_variable('lon_bnds')
            data = cf.Data([[np.float32(self.lon_min+i*cell_width), np.float32(self.lon_min+(i+1)*cell_width)] for i in range(0,array_width)], units='degrees_east', dtype='f4')
            b.set_data(data)
            c.set_bounds(b)
            field.set_construct(c, axes=('domainaxis2',), key='dimensioncoordinate2', copy=False)

        comment = self.header

        if self.comment:
            comment += self.comment

        if comment:
            field.nc_set_global_attribute("comment",comment)

        return field

    def close(self):
        """Close the formatter and flush all changes to disk"""
        years = []
        months = []
        days = []
        doys = []
        uncertainties = []
        sea_ice_fractions = []
        self.times = []
        ssts_or_anomalies = []
        self.min_dt = self.datapoints[0][0]
        self.max_dt = self.datapoints[-1][0]

        # go through each of the collected data points and append the fields to arrays, ready to package into cf.Field
        for datapoint in self.datapoints:
            (sdt,dt,edt,sst_or_anomaly,uncertainty,sea_ice_fraction) = datapoint
            self.times.append((TimeSeriesUtils.seconds_since_1981(sdt),TimeSeriesUtils.seconds_since_1981(dt),TimeSeriesUtils.seconds_since_1981(edt)))
            years.append(dt.year)
            months.append(dt.month)
            days.append(dt.day)
            doys.append(dt.timetuple().tm_yday)
            ssts_or_anomalies.append(sst_or_anomaly)
            uncertainties.append(uncertainty)
            sea_ice_fractions.append(sea_ice_fraction)

        # create the fields that will be written to file
        year_field = self.createField(years,"calendar_year","calendar year","1",scalar_int=True)
        month_field = self.createField(months, "calendar_month", "calendar month","1",scalar_int=True)
        day_field = self.createField(days, "day_of_month", "day of month","1",scalar_int=True)
        doy_field = self.createField(doys, "day_of_year", "day of year","1",scalar_int=True)
        sst_or_anomaly_field = self.createField(ssts_or_anomalies,
                    "sst_anomaly" if self.output_anomaly else "sst",
                    "sea_water_temperature_anomaly" if self.output_anomaly else "sea_water_temperature","K")
        uncertainty_field = self.createField(uncertainties, "sst_uncertainty", "sea_water_temperature uncertainty","K")
        fields = [year_field,month_field,day_field,doy_field,sst_or_anomaly_field,uncertainty_field]
        if self.output_sea_ice:
            sea_ice_field = self.createField(sea_ice_fractions, "sea_ice_area_fraction", "sea_ice_area_fraction","K")
            fields.append(sea_ice_field)

        # write out to file
        fl = cf.FieldList(fields)
        cf.write(fl,self.output_path,datatype={np.dtype('float64'): np.dtype('float32')}, compress=1)
        fl.close()
