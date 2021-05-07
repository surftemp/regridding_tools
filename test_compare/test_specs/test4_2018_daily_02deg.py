spec = {
  "name": "test4",
  "description": "test regridding to 0.2 degrees daily for 1 month",
  "spec": {
    "time_step": "N-daily",
    "n_daily_step": "1",
    "longitude_step": "0.2",
    "latitude_step": "0.2",
    "email_address": "test@test.com",
    "start_date": "2018-01-01",
    "end_date": "2018-12-31",
    "exclude_sea_ice": True,
    "sea_ice_threshold": "100",
    "anomaly_or_absolute": "absolute",
    "spatial_lambda": "1.0",
    "tau": "3",
    "generate_sea_ice_fraction": True,
    "include_bias_adjustments": True
  },
  "region_tests": [],
  "timeseries_tests": [(-10,85),(-10,80),(-10,75),(-10,70),(-10,65),(-10,60),(30,-60),(30,-50),(30,-40),(30,-30)],
  "expected": {
    "shape": (365,180*5,360*5),
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