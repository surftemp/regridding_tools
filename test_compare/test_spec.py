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

"""
This module encapsulates the server side business logic for the regridding service.

The service supports two kinds of similar but distinct job: "regrid" and "timeseries" and the business logic varies accordingly
"""

import datetime
import os.path
import subprocess


class TestSpec(object):
    """
    Define the Specification for a regrid|timeseries job received from the web form and business logic for creating and executing
    tasks to accomplish the job

    The content of a specification is encoded as a python dict object
    """

    def __init__(self,sst_cci_analysis_l4_path, c3s_sst_analysis_l4_path, sst_cci_climatology_path, reslice_path, python_exe_path):
        self.sst_cci_analysis_l4_path = sst_cci_analysis_l4_path
        self.c3s_sst_analysis_l4_path = c3s_sst_analysis_l4_path
        self.sst_cci_climatology_path = sst_cci_climatology_path
        self.reslice_path = reslice_path
        self.python_exe_path = python_exe_path

    def run(self,spec,job_type,output_folder):
        """
        Construct and return a command line for executing a regridding or timeseries task from this spec, expressed as a list of strings

        ref https://github.com/surftemp/common/blob/master/scripts/cci/makegriddedSSTs.py
        """

        time_step = spec[TestSpec.TIME_STEP]

        start_date = datetime.datetime.strptime(spec[TestSpec.START_DATE],TestSpec.DATE_FORMAT)
        end_date = datetime.datetime.strptime(spec[TestSpec.END_DATE],TestSpec.DATE_FORMAT)

        anomaly_or_absolute = spec[TestSpec.ANOMALY_OR_ABSOLUTE]
        exclude_sea_ice = spec[TestSpec.EXCLUDE_SEA_ICE]
        generate_sea_ice_fraction = spec[TestSpec.GENERATE_SEA_ICE_FRACTION]
        tau = spec[TestSpec.TAU]
        spatial_lambda = spec[TestSpec.SPATIAL_LAMBDA]

        # if the user has specified a time step of N days, set the time step to the value of the
        # control specifying how many days
        if time_step == "N-daily":
            time_step = spec[TestSpec.N_DAILY_STEP]

        # work out the fraction of sea ice above which cells should be ignored
        if exclude_sea_ice:
            sea_ice_threshold = int(spec[TestSpec.SEA_ICE_THRESHOLD]) / 100.0
        else:
            sea_ice_threshold = 1.0

        if job_type == TestSpec.JOB_TYPE_REGRID:
            lat_step = spec[TestSpec.LATITUDE_STEP]
            lon_step = spec[TestSpec.LONGITUDE_STEP]
            positionals = [lon_step, lat_step, time_step]
        else:
            lon_min = float(spec[TestSpec.BOUNDING_BOX_LONGITUDE_MIN])
            lon_max = lon_min + float(spec[TestSpec.BOUNDING_BOX_LONGITUDE_WIDTH])
            lat_min = float(spec[TestSpec.BOUNDING_BOX_LATITUDE_MIN])
            lat_max = lat_min + float(spec[TestSpec.BOUNDING_BOX_LATITUDE_HEIGHT])
            positionals = ["%0.2f"%(lon_min),"%0.2f"%(lon_max),"%0.2f"%(lat_min),"%0.2f"%(lat_max), time_step]

        if job_type == "regrid" and start_date.year != end_date.year:
            raise Exception("start and end date must be in the same year")

        kwargs = []
        if job_type == TestSpec.JOB_TYPE_REGRID:
            kwargs += ["--year",str(start_date.year)]
        else:
            kwargs += ["--start_year", str(start_date.year)]
            kwargs += ["--end_year", str(end_date.year)]
        kwargs += ["--start_month",str(start_date.month),
                   "--start_day",str(start_date.day)]
        kwargs += ["--end_month", str(end_date.month),
                   "--end_day", str(end_date.day)]
        kwargs += ["--f_max", str(sea_ice_threshold)]
        kwargs += ["--tau", str(int(float(tau)))]
        kwargs += ["--spatial_lambda", spatial_lambda]
        if not generate_sea_ice_fraction:
            kwargs += ["--no_sea_ice_fraction"]
        if anomaly_or_absolute == TestSpec.ANOMALY:
            kwargs += ["--anomalies"]

        if job_type == TestSpec.JOB_TYPE_REGRID:
            kwargs += ["--sst_cci_analysis_l4_path", self.sst_cci_analysis_l4_path]
            kwargs += ["--c3s_sst_analysis_l4_path", self.c3s_sst_analysis_l4_path]
            kwargs += ["--sst_cci_climatology_path", self.sst_cci_climatology_path]
            output_zip_name = TestSpec.createOutputZipName(anomaly_or_absolute,lat_step,start_date,end_date)
            kwargs += ["--out_path",output_folder]
            cmdline_out = (positionals+kwargs, output_zip_name)
        else: # timeseries | region
            kwargs += ["--in_path", self.reslice_path]
            if job_type == TestSpec.JOB_TYPE_TIMESERIES:
                output_format = spec[TestSpec.TIMESERIES_OUTPUT_FORMAT]
                output_file_name = TestSpec.createOutputTimeseriesName(output_format,anomaly_or_absolute,lon_min,lon_max,lat_min,lat_max,start_date,end_date)
            else:
                output_file_name = TestSpec.createOutputRegionName(anomaly_or_absolute, lon_min, lon_max,
                                                                   lat_min, lat_max, start_date, end_date)
            kwargs += ["--out_path", os.path.join(output_folder,output_file_name)]
            cmdline_out = (positionals+kwargs, output_file_name)

        if job_type == TestSpec.JOB_TYPE_REGRID:
            script_path = "../global_regridding/makegriddedSSTs.py"
        elif job_type == TestSpec.JOB_TYPE_TIMESERIES:
            script_path = "../timeseries_region/maketimeseriesSSTs.py"
        elif job_type == TestSpec.JOB_TYPE_REGION:
            script_path = "../timeseries_region/makeregionSSTs.py"
        else:
            raise Exception("Invalid job type: "+job_type)

        (script_folder,script_filename) = os.path.split(script_path)
        cmd = " ".join([self.python_exe_path,script_filename] + cmdline_out[0])
        print("Running task: "+cmd)
        subp = subprocess.Popen(cmd, cwd=script_folder,shell=True)
        subp.wait()
        rcode = subp.returncode
        if rcode != 0:
            print("Task failed with return code %d"%(rcode))
            raise Exception("Failed")
        return cmdline_out[1]

    @staticmethod
    def createOutputZipName(anomaly_or_absolute, lat_step, start_date, end_date):
        """Create the fixed length name of a zip file incorporating important regrid job attributes"""
        # based on suggestion from Laura SST_XXXX_ZZZ_YYYYMMDD_YYYYMMDD_regridded.zip
        anon_abs = anomaly_or_absolute[:4].upper()
        spatial = lat_step.replace(".", "")
        while len(spatial) < 3:
            spatial += "0"
        start_s = start_date.strftime(TestSpec.FILENAME_DATE_FORMAT)
        end_s = end_date.strftime(TestSpec.FILENAME_DATE_FORMAT)
        return "SST_%s_%s_%s_%s_regridded.zip" % (anon_abs, spatial, start_s, end_s)

    @staticmethod
    def createOutputTimeseriesName(output_format,anomaly_or_absolute, lon_min,lon_max,lat_min,lat_max, start_date, end_date):
        """Create the name of a csv|netcdf4 file incorporating important timeseries job attributes"""
        # format will be SST_XXXX_MINLON_MAXLON_MINLAT_MAXLAT_YYYYMMDD_YYYYMMDD_timeseries.[csv|nc]
        anon_abs = anomaly_or_absolute[:4].upper()
        spatial = TestSpec.formatLon(lon_min) + "_" + TestSpec.formatLon(lon_max) + "_" + TestSpec.formatLat(lat_min) + "_" + TestSpec.formatLat(lat_max)
        start_s = start_date.strftime(TestSpec.FILENAME_DATE_FORMAT)
        end_s = end_date.strftime(TestSpec.FILENAME_DATE_FORMAT)
        return "SST_%s_%s_%s_%s_timeseries.%s" % (anon_abs, spatial, start_s, end_s,"csv" if output_format == "csv" else "nc")

    @staticmethod
    def createOutputRegionName(anomaly_or_absolute, lon_min, lon_max, lat_min, lat_max, start_date,
                                   end_date):
        """Create the name of a netcdf4 file incorporating important region job attributes"""
        # format will be SST_XXXX_MINLON_MAXLON_MINLAT_MAXLAT_YYYYMMDD_YYYYMMDD_region.nc
        anon_abs = anomaly_or_absolute[:4].upper()
        spatial = TestSpec.formatLon(lon_min) + "_" + TestSpec.formatLon(lon_max) + "_" + TestSpec.formatLat(
            lat_min) + "_" + TestSpec.formatLat(lat_max)
        start_s = start_date.strftime(TestSpec.FILENAME_DATE_FORMAT)
        end_s = end_date.strftime(TestSpec.FILENAME_DATE_FORMAT)
        return "SST_%s_%s_%s_%s_region.%s" % (anon_abs, spatial, start_s, end_s, "nc")

    @staticmethod
    def formatLat(lat):
        """format a latitude value for inclusion within a filename, eg 20 becomes 20N, -30.25 becomes 30:25S"""
        s1 = "%d"%abs(lat) if lat == int(lat) else ("%.2f"%(abs(lat))).replace(".",":")
        return s1 + ("N" if lat > 0 else "S")

    @staticmethod
    def formatLon(lon):
        """format a longitude value for inclusion within a filename, eg 10.85 becomes 10:85N"""
        s1 = "%d" % abs(lon) if lon == int(lon) else ("%.2f" % (abs(lon))).replace(".", ":")
        return s1 + ("E" if lon > 0 else "W")

    @staticmethod
    def getDatesFromOutputName(output_name):
        """extract and return the start and end date strings encoded into a file name by Spec.createOutputZipName or Spec.createOutputRegionName"""
        # all files end with YYYYMMDD_YYYYMMDD_<TYPE>.nc where TYPE is one of "timeseries,regridded,region"
        last_underscore = output_name.rindex("_")
        return (output_name[last_underscore-17:last_underscore-9],output_name[last_underscore-8:last_underscore])

    # fields from the user submitted form (matching the names in templates/index.html)
    # generally speaking the same field names are used in the spec objects used to define tasks and jobs

    JOB_TYPE = "job_type"
    JOB_TYPE_REGRID = "regrid"
    JOB_TYPE_TIMESERIES = "timeseries"
    JOB_TYPE_REGION = "region"

    TIME_STEP = "time_step"
    N_DAILY_STEP = "n_daily_step"
    START_DATE = "start_date"
    END_DATE = "end_date"
    EXCLUDE_SEA_ICE = "exclude_sea_ice"
    SEA_ICE_THRESHOLD = "sea_ice_threshold"
    EMAIL_ADDRESS = "email_address"
    ANOMALY_OR_ABSOLUTE = "anomaly_or_absolute"
    ANOMALY = "anomaly"
    ABSOLUTE = "absolute"
    GENERATE_SEA_ICE_FRACTION = "generate_sea_ice_fraction"
    TAU = "tau"
    SPATIAL_LAMBDA = "spatial_lambda"
    INCLUDE_BIAS_ADJUSTMENTS = "include_bias_adjustments"

    # regrid jobs only
    LONGITUDE_STEP = "longitude_step"
    LATITUDE_STEP = "latitude_step"

    # timeseries|region jobs only
    BOUNDING_BOX_LONGITUDE_MIN = "bounding_box_longitude_min"
    BOUNDING_BOX_LATITUDE_MIN = "bounding_box_latitude_min"
    BOUNDING_BOX_LONGITUDE_WIDTH = "bounding_box_longitude_width"
    BOUNDING_BOX_LATITUDE_HEIGHT = "bounding_box_latitude_height"

    # timeseries only
    TIMESERIES_OUTPUT_FORMAT = "timeseries_output_format"
    TIMESERIES_OUTPUT_FORMAT_CSV = "csv"
    TIMESERIES_OUTPUT_FORMAT_NETCDF4 = "netcdf4"

    # format for dates in START_DATE and END_DATE in the specification
    DATE_FORMAT = "%Y-%m-%d"

    # format for dates output in filename, eh YYYYMMDD
    FILENAME_DATE_FORMAT = "%Y%m%d"


    VALID_LAT_STEPS = ["0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.4", "0.45", "0.5", "0.6", "0.75", "0.8", "0.9", "1.0", "1.2", "1.25", "1.5", "1.8", "2.0", "2.25", "2.4", "2.5", "3.0", "3.6", "3.75", "4.0", "4.5", "5.0"]
    VALID_LON_STEPS = ["0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.4", "0.45", "0.5", "0.6", "0.75", "0.8", "0.9", "1.0", "1.2", "1.25", "1.5", "1.6", "1.8", "2.0", "2.25", "2.4", "2.5", "3.0", "3.6", "3.75", "4.0", "4.5", "4.8", "5.0"]

    # list of the time steps that are valid for regridding jobs
    VALID_TIME_STEPS_REGRID = ["annual","monthly","10-day","5-day","N-daily"]
    # list of the time steps that are valid for timeseries|region jobs
    VALID_TIME_STEPS_TIMESERIES = ["monthly", "10-day", "5-day", "N-daily", "daily"]
