import numpy as np
import numpy.typing as npt
from datetime import time as Time
from datetime import date as Date
from datetime import datetime as DateTime
from datetime import timedelta


def time_to_decimal_hour(time: Time) -> float:
    """Converts time into decimal hours."""
    milliseconds = time.hour * 3.6e6
    milliseconds += time.minute * 6.0e4
    milliseconds += time.second * 1.0e3
    milliseconds += time.microsecond / 1.0e3
    hour = milliseconds / 3.6e6
    return hour


def decimal_hour_to_time(t_hr: npt.ArrayLike) -> Time:
    """Converts decimal hour into time."""
    milliseconds = t_hr * 3.6e6
    hours = milliseconds // 3.6e6
    milliseconds %= 3.6e6
    minutes = milliseconds // 6e4
    milliseconds %= 6e4
    seconds = milliseconds // 1e3
    milliseconds %= 1e3
    return Time(int(hours), int(minutes), int(seconds), int(milliseconds * 1e3))


def day_number(date: Date) -> int:
    """Returns the number of the day of the year of a given date."""
    time_delta = date - Date(date.year, 1, 1)
    return time_delta.days + 1


def day_number_to_date(n: int, year: int = 2023) -> Date:
    """Convert day number `n` to a Python date object."""
    date = Date(year, 1, 1) + timedelta(days=n - 1)
    return date


def hour_angle_to_solar_time(omega: npt.ArrayLike) -> npt.ArrayLike:
    """Returns solar time in decimal hours, given hour angle `omega`
    in radians.
    """
    t_sol = np.degrees(omega) / 15 + 12
    return t_sol


def _time_correction(dt_clock: DateTime, L_loc: float, n: int) -> float:
    """The time correction factor in minutes accounts for the variation of
    the local solar time within a given time zone due to the longitude
    variations within the time zone and also incorporates the equation of time,
    an empirical equation that corrects for the eccentricity of the Earth's orbit
    and the Earth's axial tilt.
    """
    if L_loc > 0.0:  # east of Greenwich
        L_loc = 360.0 - L_loc
    else:
        # Greenwich or west of Greenwich
        L_loc = abs(L_loc)

    utc_offset = dt_clock.utcoffset().total_seconds() / 3600.0
    L_st = 15.0 * utc_offset
    if L_st > 0.0:  # east of Greenwich
        L_st = 360.0 - L_st
    else:
        # Greenwich or west of Greenwich
        L_st = abs(L_st)

    B = np.radians((n - 1) * 360 / 365)
    E = (229.2 * (
        7.5e-5 + 1.868e-3 * np.cos(B) - 0.032077 * np.sin(B)
        - 0.014615 * np.cos(2 * B) - 0.04089 * np.sin(2 * B)
    ))
    tc = 4.0 * (L_st - L_loc) + E
    return tc


def clock_to_solar_time(dt_clock: DateTime, L_loc: float) -> DateTime:
    """Convert a time-zone aware local standard date-time `dt_clock` at longitude
    `L_loc` to solar date-time. Longitude `L_loc` in decimal degrees must be
    positive east of Greenwich and negative west of Greenwich.

    Notes
    -----
    To create a time-zone aware date-time object, we can use the package `pytz`.
    ```
    tz = pytz.timezone('Europe/Brussels')
    dt_naive = datetime(2023, 7, 17, 12, 0, 0)
    dt_std = tz.localize(dt_naive)
    ```
    """
    n = day_number(dt_clock.date())
    tc = _time_correction(dt_clock, L_loc, n)
    t_std = time_to_decimal_hour(dt_clock.time())
    t_solar = t_std + tc / 60.0
    if t_solar < 0.0:
        t_solar += 24.0
        n -= 1
    if t_solar >= 24.0:
        t_solar -= 24.0
        n += 1
    dt_solar = DateTime.combine(
        day_number_to_date(n, year=dt_clock.year),
        decimal_hour_to_time(t_solar),
        tzinfo=dt_clock.tzinfo
    )
    return dt_solar


def solar_to_clock_time(dt_solar: DateTime, L_loc: float) -> DateTime:
    """Convert time-zone aware solar date-time `dt_solar` at longitude `L_loc`
    to local standard date-time. Longitude `L_loc` is in degrees and must be
    positive east of UTC and negative west of UTC.
    """
    n = day_number(dt_solar.date())
    tc = _time_correction(dt_solar, L_loc, n)
    t_solar = time_to_decimal_hour(dt_solar.time())
    t_std = t_solar - tc / 60.0
    if t_std < 0.0:
        t_std += 24.0
        n -= 1
    if t_std >= 24.0:
        t_std -= 24.0
        n += 1
    dt_local = DateTime.combine(
        day_number_to_date(n, year=dt_solar.year),
        decimal_hour_to_time(t_std),
        tzinfo=dt_solar.tzinfo
    )
    return dt_local
