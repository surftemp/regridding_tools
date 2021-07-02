
import datetime

# values extracted from cfpython for 2003 .nc file using field.coord('T').array
# seem to be hours since 1900
val1 = 902892 # start 2003
val2 = 904476 # in the middle somewhere
val3 = 911628 # end 2003


def convert(hours_since_1900):
    """
    Convert datetime stored as integer into components

    Parameters
    ----------
    hours_since_1900 integer

    Returns
    -------

    integer array [year,month_of_year,day_of_month,day_of_year]
    """
    base_date = datetime.datetime(1900,1,1,00,00,00)
    time_delta = datetime.timedelta(hours=hours_since_1900)
    target_date = base_date + time_delta
    return [target_date.year,target_date.month,target_date.day,target_date.timetuple().tm_yday]

print(convert(val1))
print(convert(val2))
print(convert(val3))
