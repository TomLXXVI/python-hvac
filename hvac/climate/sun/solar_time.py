from datetime import datetime as DateTime
from datetime import time as Time
from datetime import date as Date
import math


def time_to_decimal_hour(t: Time) -> float:
    """Express time in decimal hours."""
    return t.hour + t.minute / 60.0 + t.second / 3600.0


def time_from_decimal_hour(t: float) -> Time:
    """Convert time from decimal hours to conventional time format."""
    # if t < 0: t += 24.0
    # if t > 24: t -= 24.0
    seconds = t * 3600.0
    hour = int(seconds // 3600)
    seconds = seconds % 3600
    minutes = int(seconds // 60)
    seconds = int(seconds % 60)
    return Time(hour, minutes, seconds)


def day_number(date: Date) -> int:
    """Get the number of the day in the year."""
    dt = date - Date(date.year, 1, 1)
    return dt.days + 1


def equation_of_time(day_num: int) -> float:
    """
    Returns the time difference in hours between solar noon and noon at local civil time due to the nature of the
    earth's orbital motion around the sun.
    """
    B = 2 * math.pi * (day_num - 81) / 364.0
    E_t = 9.87 * math.sin(2 * B) - 7.53 * math.cos(B) - 1.5 * math.sin(B)
    return E_t / 60  # expressed in decimal hours


def solar_time(datetime: DateTime, lon: float) -> Time:
    """
    Returns the solar time corresponding with the given datetime at the given longitude.
    """
    utc_offset = datetime.utcoffset().total_seconds() / 3600.0
    lon_ref = 15.0 * utc_offset
    n = day_number(datetime.date())
    t_civ = time_to_decimal_hour(datetime.time()) + (lon - lon_ref) / 15.0
    if t_civ < 0.0:
        t_civ += 24.0
        n -= 1
    if t_civ > 24.0:
        t_civ -= 24.0
        n += 1
    EoT = equation_of_time(n)
    t_sol = t_civ + EoT
    if t_sol < 0.0:
        t_sol += 24.0
    if t_sol > 24.0:
        t_sol -= 24.0
    return time_from_decimal_hour(t_sol)


def hour_angle(datetime: DateTime, lon: float) -> float:
    """
    Get solar hour angle in degrees at given date and time at the given longitude.
    """
    t_sol = time_to_decimal_hour(solar_time(datetime, lon))
    return 15 * (t_sol - 12)
