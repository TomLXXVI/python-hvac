"""
Derive from measured insolation data the components of solar irradiance on the
horizontal surface at a given instance of time.
"""
from typing import Tuple, List
import math
from datetime import datetime as DateTime
from datetime import time as Time
from datetime import date as Date
import pandas as pd
from hvac import Quantity
from .location import Location
from .clear_sky import get_I0_norm
from .solar_time import day_number


Q_ = Quantity


def solar_declination(day_num: int) -> float:
    a = 23.45 * math.pi / 180.0
    b = math.radians(360 * (day_num + 10) / 365.25)
    sin_d = -math.sin(a) * math.cos(b)
    d = math.asin(sin_d)
    return d


class MonthlyCorrelationModel:
    """
    Estimate the average hourly beam, diffuse sky and global irradiance on the
    horizontal surface from the monthly mean daily global irradiation `H_glo_hor`
    on the horizontal surface.
    """

    def __init__(self, loc: Location, datetime: DateTime, H_glo_hor: Quantity):
        """
        Parameters
        ----------
        loc : Location
            The geographic location under consideration.
        datetime : DateTime
            The date and time for which the average hourly irradiance components
            need to be calculated.
        H_glo_hor : Quantity
            The monthly mean daily global irradiation on the horizontal surface
            at the location during the given month.
        """
        self.loc: Location = loc
        self._date = datetime.date()
        self._datetime = datetime
        self.H_glo_hor = H_glo_hor.to('J / m ** 2').m
        self._w_ss = self.loc.sunset_hour_angle(self._date).to('rad').m
        self._w = self.loc.hour_angle(datetime).to('rad').m

    def _calc_H0_glo_hor(self) -> float:
        # total daily extraterrestrial irradiation on a horizontal surface
        lat = self.loc.lat.to('rad').m
        d = solar_declination(day_number(self._date))
        tau_day = 86400
        PI = math.pi
        cos = math.cos
        sin = math.sin
        I0 = get_I0_norm(day_number(self._date))
        H0 = (tau_day / PI) * I0 * cos(lat) * cos(d) * (sin(self._w_ss) - self._w_ss * cos(self._w_ss))
        return H0

    def _calc_KT(self) -> float:
        # daily clearness index: measure for the relative daily amount of global irradiation on the horizontal surface
        # with reference to the extraterrestrial daily amount of irradiation.
        H0_glo_hor = self._calc_H0_glo_hor()
        return self.H_glo_hor / H0_glo_hor

    def _estimate_H_dif(self) -> float:
        # estimate monthly average daily diffuse irradiation on the horizontal surface
        KT = self._calc_KT()
        pi = math.pi
        cos = math.cos
        R = 0.775 + 0.347 * (self._w_ss - pi / 2) - (0.505 + 0.261 * (self._w_ss - pi / 2)) * cos(2.0 * (KT - 0.9))
        H_dif = R * self.H_glo_hor
        return H_dif

    def _r_dif(self) -> float:
        # correlation function `I_dif = f(H_dif)` for the diffuse irradiance on the horizontal surface
        tau_day = 86400
        pi = math.pi
        cos = math.cos
        sin = math.sin
        r_dif = (pi / tau_day)
        r_dif *= (cos(self._w) - cos(self._w_ss))
        r_dif /= (sin(self._w_ss) - self._w_ss * cos(self._w_ss))
        return r_dif

    def _r_glo(self) -> float:
        # correlation function `I_glo_hor = f(H_glo_hor)` for the global irradiance on the horizontal surface
        a = 0.4090 + 0.5016 * math.sin(self._w_ss - math.pi / 3.0)
        b = 0.6609 - 0.4767 * math.sin(self._w_ss - math.pi / 3.0)
        r_glo = (a + b * math.cos(self._w)) * self._r_dif()
        return r_glo

    def _calc_I_glo_hor(self) -> float:
        I_glo_hor = self._r_glo() * self.H_glo_hor
        return I_glo_hor

    def _calc_I_dif(self) -> float:
        H_dif = self._estimate_H_dif()
        I_dif = self._r_dif() * H_dif
        return I_dif

    def _calc_I_beam(self, I_glo_hor: float, I_dif: float) -> float:
        theta_s = self.loc.sun_position(self._datetime).zenith.to('rad').m
        I_beam = (I_glo_hor - I_dif) / math.cos(theta_s)
        return I_beam

    def estimate(self) -> Tuple[Quantity, ...]:
        """
        Returns an estimate for the beam, diffuse and global irradiance on the horizontal surface at the location on
        the specified date and time.

        Returns
        -------
        A tuple with 3 Quantity objects:
        1. Beam irradiance
        2. Diffuse irradiance
        3. Global irradiance
        """
        sunrise = self.loc.sunrise(self._date)
        sunset = self.loc.sunset(self._date)
        if sunrise.time() <= self._datetime.time() <= sunset.time():
            I_glo_hor = self._calc_I_glo_hor()
            I_dif = self._calc_I_dif()
            I_beam = self._calc_I_beam(I_glo_hor, I_dif)
            return Q_(I_beam, 'W / m ** 2'), Q_(I_dif, 'W / m ** 2'), Q_(I_glo_hor, 'W / m ** 2')
        else:
            return Q_(0.0, 'W / m ** 2'), Q_(0.0, 'W / m ** 2'), Q_(0.0, 'W / m ** 2')

    @classmethod
    def daily_profile(cls, loc: Location, date: Date, H_glo_hor: Quantity) -> pd.DataFrame:
        """
        Returns a Pandas Dataframe with the daily profile of hourly average beam, diffuse and global irradiance on the
        horizontal surface at the given location on the given date, knowing the monthly mean daily irradiation at this
        location.
        """
        t_ax = [Time(h, 0, 0) for h in range(24)]
        t_ax.append(Time(23, 59, 59))
        t_ax = [DateTime.combine(date, t) for t in t_ax]
        I_beam_profile: List[Quantity] = []; I_dif_profile: List[Quantity] = []; I_glo_hor_profile: List[Quantity] = []
        for t in t_ax:
            mcm = cls(loc, t, H_glo_hor)
            I_beam, I_dif, I_glo_hor = mcm.estimate()
            I_beam_profile.append(I_beam)
            I_dif_profile.append(I_dif)
            I_glo_hor_profile.append(I_glo_hor)
        d = {
            'time': t_ax,
            'I_beam': [I_beam.to('W / m**2').m for I_beam in I_beam_profile],
            'I_dif': [I_dif.to('W / m**2').m for I_dif in I_dif_profile],
            'I_glo_hor': [I_glo_hor.to('W / m**2').m for I_glo_hor in I_glo_hor_profile]
        }
        df = pd.DataFrame(d)
        df.set_index('time', inplace=True)
        return df


class HourlyCorrelationModel:
    """
    Estimate average hourly beam and diffuse sky irradiance on the horizontal surface from the measured hourly
    average global irradiance `I_glo_hor` on that horizontal surface.
    """

    def __init__(self, loc: Location, datetime: DateTime, I_glo_hor: Quantity):
        """
        Parameters
        ----------
        loc: Location
            The geographic location under consideration.
        datetime : DateTime
            The date and time the average hourly irradiance components need to be calculated.
        I_glo_hor : Quantity
             The hourly average global irradiance on the horizontal surface.
        """
        self.loc = loc
        self._date = datetime.date()
        self._datetime = datetime
        self.I_glo_hor = I_glo_hor.to('W / m**2').m
        # self._w_ss = self.loc.sunset_hour_angle(self._date)('rad')
        # self._w = self.loc.hour_angle(self._datetime)('rad')

    def _calc_KT(self) -> float:
        # hourly clearness index: measure for the relative amount of global irradiance on the horizontal surface
        # with reference to the extraterrestrial amount of irradiance.
        I0 = get_I0_norm(day_number(self._date))
        theta_s = self.loc.sun_position(self._datetime).zenith.to('rad').m
        KT = self.I_glo_hor / (I0 * math.cos(theta_s))
        return KT

    def _calc_I_dif(self) -> float:
        KT = self._calc_KT()
        if 0 <= KT < 0.22:
            R = 1.0 - 0.09 * KT
        elif 0.22 <= KT < 0.80:
            R = 0.9511 - 0.1604 * KT + 4.388 * KT ** 2 - 16.638 * KT ** 3 + 12.336 * KT ** 4
        else:
            R = 0.165
        I_dif = R * self.I_glo_hor
        return I_dif

    def _calc_I_beam(self, I_dif: float) -> float:
        theta_s = self.loc.sun_position(self._datetime).zenith.to('rad').m
        I_beam = (self.I_glo_hor - I_dif) / math.cos(theta_s)
        return I_beam

    def estimate(self) -> Tuple[Quantity, ...]:
        """
        Returns an estimate for the beam and sky diffuse irradiance on the horizontal surface.
        """
        sunrise = self.loc.sunrise(self._date)
        sunset = self.loc.sunset(self._date)
        if sunrise.time() <= self._datetime.time() <= sunset.time():
            I_dif = self._calc_I_dif()
            I_beam = self._calc_I_beam(I_dif)
            return Q_(I_beam, 'W / m**2'), Q_(I_dif, 'W / m**2')
        else:
            return Q_(0.0, 'W / m**2'), Q_(0.0, 'W / m**2')
