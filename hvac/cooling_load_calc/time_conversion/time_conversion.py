from datetime import date as Date
from datetime import time as Time
from datetime import datetime as DateTime
import pytz
from hvac import Quantity
from hvac.sun.time import (
    decimal_hour_to_time,
    solar_to_clock_time,
    clock_to_solar_time,
    time_to_decimal_hour
)


def convert_to_clock_time(
    t_sol_sec: float,
    date: Date,
    L_loc: Quantity,
    tz_loc: str
) -> tuple[DateTime, DateTime]:
    """Converts local solar time in seconds from midnight to local solar time in
    date-time format and to local standard time in date-time format.

    Parameters
    ----------
    t_sol_sec:
        Solar time in seconds from midnight (midnight is 0 s).
    date:
        The date of the day under consideration.
    L_loc:
        The longitude of the location under consideration. East from Greenwich
        is a positive angle, and west is a negative angle.
    tz_loc:
        The time zone of the location according to the tz database notation.

    Returns
    -------
    A 2-tuple with:
    -   Python datetime object representing the date and local standard time that
        corresponds with the given solar time in seconds.
    -   Python datetime object representing the date and local solar time that
        corresponds with the given solar time in seconds.
    """
    t_sol_hr = t_sol_sec / 3600
    sol_time = decimal_hour_to_time(t_sol_hr)
    sol_datetime = DateTime.combine(date, sol_time)  # naive date-time
    tz_loc = pytz.timezone(tz_loc)
    sol_datetime = tz_loc.localize(sol_datetime)  # timezone aware date-time
    std_datetime = solar_to_clock_time(sol_datetime, L_loc.to('deg').m)
    return std_datetime, sol_datetime


def convert_to_solar_seconds(
    clock_time: Time,
    date: Date,
    L_loc: Quantity,
    tz_loc: str
) -> float:
    """Converts a clock time (local standard time) in a given time zone on a
    given date to solar time in seconds from midnight (0 s).

    Parameters
    ----------
    clock_time:
        The clock time being read.
    date:
        The date at which the clock time is being read.
    L_loc:
        The longitude where the clock time is being read. East of Greenwich is
        a positive angle, west of greenwich is a negative angle.
    tz_loc:
        The timezone in which the clock time is being read conform the tz
        database notation.

    Returns
    -------
    A float, being the solar time in seconds from midnight (0 s).
    """
    dt_clock = DateTime.combine(
        date=date,
        time=clock_time
    )
    tz_loc = pytz.timezone(tz_loc)
    dt_clock = tz_loc.localize(dt_clock)
    dt_sol = clock_to_solar_time(dt_clock, L_loc.to('deg').m)
    t_sol = dt_sol.time()
    t_sol_sec = time_to_decimal_hour(t_sol) * 3600
    return t_sol_sec
