spec = {
  "name": "test6",
  "description": "test regridding to 1 degrees 5-day for 1 year",
  "spec": {
    "time_step": "monthly",
    "n_daily_step": "1",
    "longitude_step": "1.0",
    "latitude_step": "1.0",
    "email_address": "test@test.com",
    "start_date": "2017-10-01",
    "end_date": "2017-12-31",
    "exclude_sea_ice": True,
    "sea_ice_threshold": "100",
    "anomaly_or_absolute": "absolute",
    "spatial_lambda": "1.0",
    "tau": "3",
    "generate_sea_ice_fraction": True,
    "include_bias_adjustments": False
  },
  "region_tests": [],
  "timeseries_tests": [(-10,85),(-10,80),(-10,75),(-10,70),(-10,65),(-10,60),(30,-60),(30,-50),(30,-40),(30,-30)],
  "expected": {
    "shape": (3,180,360),
    "fields": {
      "sea_area_fraction": {
        "range": [
          0.0,
          1.0
        ]
      },
      "sea_ice_area_fraction": {
        "range": [
          0.0,
          1.0
        ]
      },
      "sea_water_temperature": {
        "range": [
          269,
          310
        ]
      },
      "sea_water_temperature standard_error": {
        "range": [
          0,
          10.0
        ]
      }
    }
  }
}