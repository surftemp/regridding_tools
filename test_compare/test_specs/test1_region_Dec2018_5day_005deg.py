spec = {
  "name": "test1",
  "description": "test regridding vs region for 0.05 degrees, pentad, dec 2018",
  "spec": {
    "time_step": "5-day",
    "n_daily_step": "5",
    "longitude_step": "0.05",
    "latitude_step": "0.05",
    "email_address": "test@test.com",
    "start_date": "2018-12-01",
    "end_date": "2018-12-31",
    "exclude_sea_ice": True,
    "sea_ice_threshold": "100",
    "anomaly_or_absolute": "absolute",
    "spatial_lambda": "1.0",
    "tau": "3",
    "generate_sea_ice_fraction": True,
    "include_bias_adjustments": False
  },
  "region_tests": [((-10,-50),(-5,-40)),((100,-40),(120,-20)),((-30,60),(-10,70)),((-80,-50),(-70,-40)),((-5,-5),(5,5))],
  "timeseries_tests": [],
  "expected": {
    "shape": (6,3600,7200),
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