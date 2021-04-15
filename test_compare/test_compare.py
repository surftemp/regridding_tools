# -*- coding: utf-8 -*-

#    sst-services
#    Copyright (C) 2020  National Centre for Earth Observation (NCEO)
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""this module automates the running of time series and regridding jobs detailed in test specifications

It will compares the results from different time series and regrid jobs to make sure that equivalent cell values in the
results match up.
Because the timeseries and regrid calculations are performed by independently coded solutions to a task specification,
agreement in the results provides some confidence that the calculations are being performed correctly"""

import argparse
import os.path
import csv
import glob

import cf
import sys
import numpy as np
import copy
from logging import Logger, INFO, StreamHandler, Formatter
import math

def divides(num,fac):
    rem = abs(num - fac*round(num/fac))
    return rem < 0.00001

from test_specs.test1_region_Dec2018_5day_005deg import spec as spec1
from test_specs.test2_timeseries_Dec2018_5day_5deg import spec as spec2

from test_spec import TestSpec

alltests = [spec1,spec2]

logger = Logger("verification")
logger.setLevel(INFO)
sh = StreamHandler(sys.stdout)
formatter = Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
sh.setFormatter(formatter)
logger.addHandler(sh)

# regrid and timeseries output datasets have different names
# map from the field names used in the regridded data to the field names used in the timeseries data

# mapping for csv output timeseries format
absolute_field_name_mapping_csv = {
    'sea_water_temperature':'mean temperature kelvin',
    'sea_water_temperature standard_error':'mean temperature uncertainty',
    'sea_ice_area_fraction':'fraction of sea-ice-covered ocean',
    'sea_area_fraction':'', # no equivalent in ts output
}

anomaly_field_name_mapping_csv = {
    'sea_water_temperature_anomaly': 'mean temperature anomaly kelvin',
    'sea_water_temperature standard_error':'mean temperature anomaly uncertainty',
    'sea_ice_area_fraction': 'fraction of sea-ice-covered ocean',
    'sea_area_fraction':'' # no equivalent in ts output
}

# mapping for netcdf4 output timeseries format
absolute_field_name_mapping_netcdf4 = {
    'sea_water_temperature':'sea_water_temperature',
    'sea_water_temperature standard_error':'sea_water_temperature uncertainty',
    'sea_ice_area_fraction':'sea_ice_area_fraction',
    'sea_area_fraction':'', # no equivalent in ts output
}

anomaly_field_name_mapping_netcdf4 = {
    'sea_ice_area_fraction':'sea_ice_area_fraction',
    'sea_area_fraction':'', # no equivalent in ts output
    'sea_water_temperature_anomaly': 'sea_water_temperature_anomaly',
    'sea_water_temperature standard_error': 'sea_water_temperature uncertainty'
}

# mapping for netcdf4 output region format
region_absolute_field_name_mapping_netcdf4 = {
    'sea_water_temperature':'sea_water_temperature',
    'sea_water_temperature standard_error':'sea_water_temperature uncertainty',
    'sea_ice_area_fraction':'sea_ice_area_fraction',
    'sea_area_fraction':'', # no equivalent in ts output
}

region_anomaly_field_name_mapping_netcdf4 = {
    'sea_ice_area_fraction':'sea_ice_area_fraction',
    'sea_area_fraction':'', # no equivalent in ts output
    'sea_water_temperature_anomaly': 'sea_water_temperature_anomaly',
    'sea_water_temperature standard_error': 'sea_water_temperature uncertainty'
}

class Verifier(object):

    """compare regrid results with one or more timeseries results for regrid and timeseries jobs with identical overlapping parameters"""

    def __init__(self,field_names,spec,regrid_folder,timeseries_locations,region_locations,tolerance_decimal_places=3):
        self.field_names = field_names
        self.spec = spec
        self.tolerance_decimal_places = tolerance_decimal_places
        self.failures = []
        self.regrid_folder = regrid_folder
        self.timeseries_locations = timeseries_locations
        self.region_locations = region_locations
        self.regrid_fields = None

    def open(self):
        self.regrid_fields = cf.read(os.path.join(self.regrid_folder,"*.nc"),aggregate={"relaxed_identities": True})

    def verify(self,name,is_anomaly):
        logger.info("Verifying results for test: %s",name)

        for (min_lon,min_lat,ts_path) in self.timeseries_locations:
            max_lat = min_lat + float(self.spec["latitude_step"])
            max_lon = min_lon + float(self.spec["longitude_step"])

            test = " at lat: (%0.2f,%0.2f) lon: (%0.2f,%0.2f)" % (min_lat, max_lat, min_lon, max_lon)
            logger.info("%s: Verifying timeseries"%(test))

            values = {} # map from field name to a list of values in the time series
            for field_name in self.field_names:
                ts_field_name = anomaly_field_name_mapping_csv[field_name] if is_anomaly else absolute_field_name_mapping_csv[field_name]
                if ts_field_name == '':
                    continue
                values[field_name] = []

            if ts_path.endswith(".nc"):
                ts_fields = cf.read(ts_path)
                for field_name in self.field_names:
                    ts_field_name = anomaly_field_name_mapping_netcdf4[field_name] if is_anomaly else absolute_field_name_mapping_netcdf4[field_name]
                    if ts_field_name == '':
                        continue
                    extracted_ts = ts_fields.select(ts_field_name)[0].array
                    if not np.ma.is_masked(extracted_ts) or not np.alltrue(extracted_ts.mask):
                        ts_values = np.ma.filled(extracted_ts, fill_value=math.nan)
                        values[field_name] = ts_values.flatten().tolist()
                    else:
                        # all values appear to be masked, create a list of NaNs
                        values[field_name] = [math.nan]*extracted_ts.shape[0]
                ts_fields.close()
            else: # CSV
                reader = csv.reader(open(ts_path))
                index = 0
                headers = {}

                for row in reader:
                    index += 1
                    if index < 3:
                        # skip first two rows containing comments
                        continue
                    if index == 3:
                        # read the field names from third row and setup
                        for idx in range(len(row)):
                            headers[row[idx]] = idx
                    else:
                        # read data rows from subsequent rows
                        for field_name in self.field_names:
                            ts_field_name = anomaly_field_name_mapping_csv[field_name] if is_anomaly else absolute_field_name_mapping_csv[field_name]
                            if ts_field_name == '':
                                continue
                            idx = headers[ts_field_name]
                            try:
                                value = float(row[idx])
                            except:
                                value = math.nan
                            values[field_name].append(value)

            for field_name in self.field_names:
                try:
                    if field_name in values:
                        self.testLocation(field_name + test, np.array(values[field_name]), min_lon, min_lat, field_name)
                except Exception as ex:
                    self.failures.append((field_name, min_lon, min_lat, max_lon, max_lat, ex))

        for (region_boundaries,region_paths) in self.region_locations:
            min_lon = region_boundaries[0][0]
            min_lat = region_boundaries[0][1]
            max_lon = region_boundaries[1][0]
            max_lat = region_boundaries[1][1]

            test = " at lat: (%0.2f,%0.2f) lon: (%0.2f,%0.2f)" % (min_lat, max_lat, min_lon, max_lon)
            logger.info("%s: Verifying region"%(test))

            values = {} # map from field name to the array

            region_fields = cf.read(region_paths)

            for field_name in self.field_names:
                region_field_name = region_anomaly_field_name_mapping_netcdf4[field_name] if is_anomaly else region_absolute_field_name_mapping_netcdf4[field_name]
                if region_field_name == '':
                    continue
                extracted_region = region_fields.select(region_field_name)[0].array
                region_values = np.ma.filled(extracted_region, fill_value=math.nan)
                values[field_name] = region_values

            region_fields.close()

            for field_name in self.field_names:
                try:
                    if field_name in values:
                        self.testRegion(field_name + test, values[field_name], min_lon, min_lat, max_lon, max_lat, field_name)
                except Exception as ex:
                    self.failures.append((field_name, min_lon, min_lat, max_lon, max_lat, ex))


        if len(self.failures)>0:
            for (field_name,min_lon,min_lat,max_lon,max_lat,ex) in self.failures:
                logger.error("Test %s failed field: %s lat: (%0.2f,%0.2f) lon:(%0.2f,%0.2f) error=%s",name,field_name,min_lat,max_lat,min_lon,max_lon,str(ex))

    def testLocation(self,test, ts_values,min_lon,min_lat,field_name):

        # obtain the time series from the regridded dataset
        lat_index = round((90 + min_lat) / float(self.spec["latitude_step"]))
        lon_index = round((180 + min_lon) / float(self.spec["longitude_step"]))

        regrid_field = self.regrid_fields.select(field_name)[0]
        regrid_ts = regrid_field.subspace[:, lat_index, lon_index].array

        regrid_values = np.ma.filled(regrid_ts,fill_value=math.nan)
        if field_name == "sea_ice_area_fraction":
            sea_area_fraction = self.regrid_fields.select("sea_area_fraction")[0].subspace[:, lat_index, lon_index].array
            sea_area_fraction = np.ma.filled(sea_area_fraction, fill_value=0.0)
            regrid_values = np.where(sea_area_fraction > 0,regrid_values,np.nan)


        regrid_values = regrid_values.flatten()
        # print("regridded:",regrid_filled.flatten().tolist())
        abs_diffs = np.abs(ts_values-regrid_values)
        # print("absdiffs:", abs_diffs.flatten().tolist())
        valid_count = np.count_nonzero(~np.isnan(ts_values))
        if valid_count > 0:
            logger.info("%s: Valid=%d Extract Mean Value=%0.4f, Biggest Difference=%0.8f" % (test, valid_count, np.nanmean(ts_values), np.nanmax(abs_diffs)))
        else:
            logger.info("%s: Valid=%d" % (test, valid_count))
        if regrid_ts.shape[0] != ts_values.shape[0]:
            msg = "Different number of values in series %s"%(test)
            logger.error(msg)
            raise Exception(msg)

        # compare data within tolerance
        np.testing.assert_array_almost_equal(regrid_values,ts_values,self.tolerance_decimal_places)

    def testRegion(self,test, region_values,min_lon,min_lat,max_lon,max_lat,field_name):
        # obtain the region from the regridded dataset
        lat_index_min = round((90 + min_lat) / 0.05)
        lon_index_min = round((180 + min_lon) / 0.05)
        lat_index_max = round((90 + max_lat) / 0.05)
        lon_index_max = round((180 + max_lon) / 0.05)

        regrid_field = self.regrid_fields.select(field_name)[0]

        regrid_area = regrid_field.subspace[:, lat_index_min:lat_index_max, lon_index_min:lon_index_max].array
        regrid_values = np.ma.filled(regrid_area, fill_value=math.nan)
        if field_name == "sea_ice_area_fraction":
            sea_area_fraction = self.regrid_fields.select("sea_area_fraction")[0].subspace[:, lat_index_min:lat_index_max, lon_index_min:lon_index_max].array
            sea_area_fraction = np.ma.filled(sea_area_fraction, fill_value=0.0)
            regrid_values = np.where(sea_area_fraction > 0,regrid_values,np.nan)

        if regrid_values.shape != region_values.shape:
            msg = "Different shapes %s"%(test)
            logger.error(msg)
            raise Exception(msg)

        abs_diffs = np.abs(region_values-regrid_values)

        valid_count = np.count_nonzero(~np.isnan(regrid_values))
        if valid_count > 0:
            logger.info("%s: Valid=%d Extract Mean Value=%0.4f, Biggest Difference=%0.8f" % (test, valid_count, np.nanmean(regrid_values), np.nanmax(abs_diffs)))
        else:
            logger.info("%s: Valid=%d" % (test, valid_count))

        # compare data within tolerance
        np.testing.assert_array_almost_equal(regrid_values,region_values,self.tolerance_decimal_places)

    def close(self):
        self.regrid_fields.close()

class TestRunner(object):
    """Handle the running of a suite of test cases - submitting the regrid and time series jobs, and checking the results"""

    def __init__(self,output_folder,cci_path,c3s_path,climatology_path,reslice_path,python_exe_path):
        self.cci_path = cci_path
        self.c3s_path = c3s_path
        self.climatology_path = climatology_path
        self.reslice_path = reslice_path
        self.output_folder = output_folder
        self.python_exe_path = python_exe_path

    def formatLon(self,lon):
        return ("E" if lon >= 0 else "W") + ("%0.2f"%(abs(lon))).replace(".","p")

    def formatLat(self,lat):
        return ("N" if lat >= 0 else "S") + ("%0.2f"%(abs(lat))).replace(".","p")

    def runTest(self,test):
        name = test["name"]
        spec = test["spec"]
        region_locations = test.get("region_tests", [])
        lat_step = float(spec["latitude_step"])
        lon_step = float(spec["longitude_step"])
        timeseries_locations = test.get("timeseries_tests", [])

        subfolder = os.path.join(self.output_folder, name)
        regrid_folder = os.path.join(subfolder, "regrid")

        if not os.path.exists(regrid_folder):
            os.makedirs(regrid_folder)
            try:
                self.runRegridTest(name,spec,regrid_folder)
            except Exception as ex:
                os.rmdir(regrid_folder)
                print(ex)

        # if os.path.exists(regrid_folder):
        #    self.checkRegridResults(name,regrid_folder,test["expected"])

        region_tests = []
        for (min_lonlat, max_lonlat) in region_locations:
            (min_lon,min_lat) = min_lonlat
            (max_lon,max_lat) = max_lonlat
            results_folder = "region_%s_%s_%s_%s" % (self.formatLon(min_lon), self.formatLat(min_lat),self.formatLon(max_lon), self.formatLat(max_lat) )
            region_folder = os.path.join(subfolder, results_folder)
            if not os.path.exists(region_folder):
                os.makedirs(region_folder)
                try:
                    self.runRegionTest(name,min_lon,min_lat,max_lon,max_lat,spec,region_folder)
                except Exception as ex:
                    region_tests.append(str(ex))
                    raise ex

            if os.path.exists(region_folder):
                self.checkRegionResults(name,min_lon,min_lat,max_lon,max_lat,region_folder,spec,test["expected"])
                region_tests.append((((min_lon, min_lat), (max_lon, max_lat)), glob.glob(region_folder+"/*")))

        timeseries_tests = []
        for (lon, lat) in timeseries_locations:
            if not divides(lon,lon_step):
                print("Skipping time series test at (%f,%f), longitude is not divisible by longitude step"%(lon,lat))
                continue
            if not divides(lat,lat_step):
                print("Skipping time series test at (%f,%f), latitude is not divisible by latitude step"%(lon,lat))
                continue
            output_formats = ["netcdf4","csv"]

            for output_format in output_formats:
                results_folder = "timeseries_%s_%s_%s_%s_%s" % (
                    self.formatLon(lon), self.formatLat(lat), self.formatLon(lon+lon_step), self.formatLat(lat+lat_step), output_format)
                timeseries_folder = os.path.join(subfolder, results_folder)
                if not os.path.exists(timeseries_folder):
                    os.makedirs(timeseries_folder)
                    try:
                        self.runTimeSeriesTest(name, lon, lat, spec, output_format,timeseries_folder)
                    except Exception as ex:
                        raise ex
                if os.path.exists(timeseries_folder):
                    self.checkTimeSeriesResults(name, lon, lat, lon_step, lat_step, timeseries_folder, output_format, spec, test["expected"])
                    timeseries_tests.append((lon, lat, glob.glob(timeseries_folder+"/*")[0]))

        if spec["anomaly_or_absolute"] == "anomaly":
            field_names = ["sea_water_temperature_anomaly"]
        else:
            field_names = ["sea_water_temperature"]

        field_names += ['sea_water_temperature standard_error', 'sea_ice_area_fraction','sea_area_fraction']

        v = Verifier(field_names, spec, regrid_folder, timeseries_tests, region_tests)
        v.open()
        v.verify(name, spec["anomaly_or_absolute"] == "anomaly")
        v.close()

    @staticmethod
    def findField(fields,field_name):
        for field in fields:
            if field.identity() == field_name:
                return field
        return None

    def runRegridTest(self,name,spec,output_folder):
        print("Running regrid test: %s"%(name))
        regrid_spec = copy.deepcopy(spec)
        regrid_spec["job_type"] = "regrid"
        ts = TestSpec(sst_cci_analysis_l4_path=self.cci_path,c3s_sst_analysis_l4_path=self.c3s_path,sst_cci_climatology_path=self.climatology_path,reslice_path=self.reslice_path,python_exe_path=self.python_exe_path)
        return ts.run(regrid_spec,"regrid",output_folder)

    def runTimeSeriesTest(self,name,lon,lat,spec,output_format,output_folder):
        print("Running time series test: %s (for lat=%f,lon=%f,output_format=%s)"%(name,lat,lon,output_format))
        timeseries_spec = copy.deepcopy(spec)
        timeseries_spec["job_type"] = "timeseries"
        lat_step = spec["latitude_step"]
        lon_step = spec["longitude_step"]
        del timeseries_spec["latitude_step"]
        del timeseries_spec["longitude_step"]
        timeseries_spec["timeseries_output_format"] = output_format
        timeseries_spec["bounding_box_longitude_min"] = str(lon)
        timeseries_spec["bounding_box_longitude_width"] = str(lon_step)
        timeseries_spec["bounding_box_latitude_min"] = str(lat)
        timeseries_spec["bounding_box_latitude_height"] = str(lat_step)
        ts = TestSpec(sst_cci_analysis_l4_path=self.cci_path, c3s_sst_analysis_l4_path=self.c3s_path,
                      sst_cci_climatology_path=self.climatology_path, reslice_path=self.reslice_path, python_exe_path=self.python_exe_path)
        ts.run(timeseries_spec,"timeseries",output_folder)
        return glob.glob(output_folder+"/*")[0]

    def runRegionTest(self,name,min_lon,min_lat,max_lon,max_lat,spec,output_folder):
        print("Running region test: %s (for area lat=%f-%f,lon=%f-%f)"%(name,min_lat,max_lat,min_lon,max_lon))
        region_spec = copy.deepcopy(spec)
        region_spec["job_type"] = "region"
        del region_spec["latitude_step"]
        del region_spec["longitude_step"]

        region_spec["bounding_box_longitude_min"] = str(min_lon)
        region_spec["bounding_box_longitude_width"] = str(max_lon-min_lon)
        region_spec["bounding_box_latitude_min"] = str(min_lat)
        region_spec["bounding_box_latitude_height"] = str(max_lat-min_lat)
        ts = TestSpec(sst_cci_analysis_l4_path=self.cci_path, c3s_sst_analysis_l4_path=self.c3s_path,
                      sst_cci_climatology_path=self.climatology_path, reslice_path=self.reslice_path,python_exe_path=self.python_exe_path)
        ts.run(region_spec, "region", output_folder)
        return glob.glob(output_folder+"/*")[0]



    def checkRegridResults(self,name,folder,expected):
        print("Checking regrid test results for: %s"%(name))

        expected_shape = expected["shape"]

        fields = cf.read(os.path.join(folder, "*.nc"))
        for fname in expected["fields"].keys():
            constraints = expected["fields"][fname]
            (range_min,range_max) = tuple(constraints["range"])
            field = TestRunner.findField(fields,fname)
            if field is None:
                raise Exception("Field %s not found"%(fname))

            print("\tChecking regrid field %s"%(fname))

            if field.array.shape != expected_shape:
                raise Exception("Field %s has unexpected shape: %s"%(fname,str(field.array.shape)))
            if not np.isfinite(field).all():
                raise Exception("Field %s has Inf or NaN values?"%(fname))
            minv = np.min(field.array)
            maxv = np.max(field.array)
            if minv < range_min:
                raise Exception("Field %s has value %f below expected minimum %f?" % (fname,minv,range_min))
            if maxv > range_max:
                raise Exception("Field %s has value %f above expected minimum %f?" % (fname,maxv,range_max))
        fields.close()

    def checkRegionResults(self,name,min_lon, min_lat,max_lon,max_lat,folder,spec,expected):
        print("Checking region test results for: %s"%(name))
        expected_shape = (expected["shape"][0],round((max_lat-min_lat)/0.05),round((max_lon-min_lon)/0.05))
        fields = cf.read(os.path.join(folder, "*.nc"))
        for fname in expected["fields"].keys():
            if spec["anomaly_or_absolute"] == "anomaly":
                mapped_fname = region_anomaly_field_name_mapping_netcdf4[fname]
            else:
                mapped_fname = region_absolute_field_name_mapping_netcdf4[fname]
            if mapped_fname == "":
                continue
            constraints = expected["fields"][fname]
            (range_min,range_max) = tuple(constraints["range"])
            field = TestRunner.findField(fields,mapped_fname)
            if field is None:
                raise Exception("Field %s not found"%(mapped_fname))
            print("\tChecking region field %s"%(mapped_fname))
            if field.array.shape != expected_shape:
                raise Exception("Field %s has unexpected shape: %s"%(fname,str(field.array.shape)))
            if not np.isfinite(field).all():
                raise Exception("Field %s has Inf or NaN values?"%(fname))
            minv = np.min(field.array)
            maxv = np.max(field.array)
            if minv < range_min:
                raise Exception("Field %s has value %f below expected minimum %f?" % (fname,minv,range_min))
            if maxv > range_max:
                raise Exception("Field %s has value %f above expected minimum %f?" % (fname,maxv,range_max))
        fields.close()

    def checkTimeSeriesResults(self,name,lon,lat,lon_step,lat_step,folder,output_format,spec,expected):
        print("Checking time series test results for: %s (for lat=%f,lon=%f,output_format=%s)"%(name,lat,lon,output_format))





if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('--tests', action='store',
                        dest='tests',
                        type=str,
                        help='Launch particular tests',
                        default="")

    parser.add_argument('--list', action='store_true',
                        dest='list_tests',
                        help='List available tests')

    parser.add_argument('--folder', action='store',
                        type=str,
                        default="/tmp",
                        help='specify parent folder for storing results')

    SST_CCI_ANALYSIS_L4_PATH = "/Users/cv922550/Data/sst/data/CDR_v2/Analysis/L4/v2.1/"
    C3S_SST_ANALYSIS_L4_PATH = "/Users/cv922550/Data/sst/data/ICDR_v2/Analysis/L4/v2.0/"
    SST_CCI_CLIMATOLOGY_PATH = "/Users/cv922550/Data/sst/data/CDR_v2/Climatology/L4/v2.1/"
    RESLICE_PATH = "/Users/cv922550/Data/reslice/"

    parser.add_argument('--reslice_path', default="/Users/cv922550/Data/reslice/",
                        help='Path to the resliced SST C3S/CCI Analysis Level 4 input data.')

    parser.add_argument('--cci_path', default="/Users/cv922550/Data/sst/data/CDR_v2/Analysis/L4/v2.1/",
                        help='Path to the SST CCI Analysis Level 4 input data.')

    parser.add_argument('--c3s_path', default="/Users/cv922550/Data/sst/data/ICDR_v2/Analysis/L4/v2.0/",
                        help='Path to the SST C3S Analysis Level 4 input data.')

    parser.add_argument('--climatology_path', default="/Users/cv922550/Data/sst/data/CDR_v2/Climatology/L4/v2.1/",
                        help='Path to the SST Climatology Analysis Level 4 input data.')

    parser.add_argument('--python_exe', default="/Users/cv922550/miniconda3/envs/cfpy37/bin/python",
                        help='Path to the python executable.')

    args = parser.parse_args()

    if args.list_tests:
        print("Available Tests")
        for test in alltests:
            print("%s (%s)"%(test["name"],test["description"]))
    else:
        tr = TestRunner(args.folder,args.cci_path,args.c3s_path,args.climatology_path,args.reslice_path,args.python_exe)
        if args.tests:
            test_names = args.tests.split(",")
            tests = list(filter(lambda t: t["name"] in test_names, alltests))
        else:
            tests = alltests
        for test in tests:
            print(test)
            tr.runTest(test)



