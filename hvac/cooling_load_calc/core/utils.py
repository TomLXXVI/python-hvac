from datetime import date as Date
from datetime import time as Time
from datetime import datetime as DateTime
import pytz
import numpy as np
from scipy.optimize import minimize
from hvac import Quantity
from hvac.sun.time import (
    decimal_hour_to_time,
    solar_to_clock_time,
    clock_to_solar_time,
    time_to_decimal_hour
)


Q_ = Quantity


class AirLayerTemperatureSolver:
    """Solves for the temperatures on both sides of an air layer, when the
    following conditions are known:
    - the exterior temperature
    - the unit thermal resistance between the exterior and the outer air layer
      side
    - the unit thermal resistance between the inner air layer side and the
      interior zone air
    - the zone air temperature.

    The solving method searches for the temperatures on the outer and the inner
    side of the air layer, so that the heat flow from the exterior to the
    outer side of the air layer balances with the heat flow from the inner side
    of the air layer to the interior zone.
    """

    def __init__(
        self,
        T_ext: Quantity,
        T_int: Quantity,
        R_ea: Quantity,
        R_ai: Quantity
    ) -> None:
        """Creates an `AirLayerTemperatureSolver` instance.

        Parameters
        ----------
        T_ext:
            Exterior temperature.
        T_int:
            Interior temperature.
        R_ea:
            Unit thermal resistance between exterior and outer air layer side.
        R_ai:
            Unit thermal resistance between inner air layer side and interior.
        """
        self.T_ext = T_ext.to('K').m
        self.T_int = T_int.to('K').m
        self.R_ea = R_ea.to('K * m**2 / W').m
        self.R_ai = R_ai.to('K * m**2 / W').m

    def _eq1(self, T_ae: float) -> float:
        """Heat flow from the exterior to the outer side of the air layer."""
        return (self.T_ext - T_ae) / self.R_ea

    def _eq2(self, T_ai: float) -> float:
        """Heat flow from the inner side of the air layer to the interior zone."""
        return (T_ai - self.T_int) / self.R_ai

    def _sys_eqs(self, unknowns: np.ndarray) -> float:
        T_ae = unknowns[0]
        T_ai = unknowns[1]
        Q_dot_ea = self._eq1(T_ae)
        Q_dot_az = self._eq2(T_ai)
        dev_Q_dot = abs(Q_dot_ea - Q_dot_az)
        return dev_Q_dot

    def solve(
        self,
        T_ae_guess: Quantity,
        T_ai_guess: Quantity
    ) -> tuple[Quantity, ...]:
        """Solves for the temperatures on the outer and the inner side of the
        air layer.

        Parameters
        ----------
        T_ae_guess:
            Initial guess for the temperature on the outer side of the air layer.
        T_ai_guess:
            Initial guess for the temperature on the inner side of the air layer.

        Returns
        -------
        1. temperature on the outer side of the air layer
        2. temperature on the inner side of the air layer
        3. temperature difference between outer and inner side of the air layer
        4. average temperature in the air layer
        """
        ini_T_ae = T_ae_guess.to('K').m
        ini_T_ai = T_ai_guess.to('K').m
        sol = minimize(self._sys_eqs, np.array([ini_T_ae, ini_T_ai]))
        T_ae = Q_(sol.x[0], 'K')
        T_ai = Q_(sol.x[1], 'K')
        dT = abs(T_ae - T_ai)
        T_avg = (T_ae + T_ai) / 2
        return T_ae, T_ai, dT, T_avg


def convert_to_clock_time(
    time_index: int,
    dt_hr: float,
    date: Date,
    L_loc: Quantity,
    tz_loc: str
) -> tuple[DateTime, DateTime]:
    """Converts a time index with given time step in local solar time and also
    to local standard time.

    Parameters
    ----------
    time_index:
        Time index indicating a moment in time of the day.
    dt_hr:
        The time step in decimal hours used to perform the dynamic heat transfer
        calculations of exterior building elements.
    date:
        The date of the day under consideration.
    L_loc:
        The longitude of the location under consideration. East from Greenwich
        is a positive angle, and west is a negative angle.
    tz_loc:
        The time zone of the location conform the tz database notation.

    Returns
    -------
    A 2-tuple with:
    -   Python datetime object representing the date and local standard time that
        corresponds with the given time index and time step.
    -   Python datetime object representing the date and local solar time that
        corresponds with the given time index and time step.
    """
    t_sol_hr_dec = time_index * dt_hr
    sol_time = decimal_hour_to_time(t_sol_hr_dec)
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
    """Converts a clock time (a local standard time) in a given time zone on a
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
    A float, being the solar time in seconds from midnight (0 s)
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
