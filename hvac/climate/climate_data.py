from typing import Optional, Callable, Dict, List, Any
import math
from datetime import date as Date
from datetime import datetime as DateTime
from datetime import time as Time
from scipy.interpolate import interp1d
from hvac import Quantity
from hvac.fluids import HumidAir
from .sun import Location, ClearSkyModel, SunPath
from .sun.solar_time import time_to_decimal_hour

Q_ = Quantity


class _DesignDayTemperatureProfile:
    """
    Only for internal use. Used in `ClimateData`.
    Calculates the outdoor dry-bulb and wet-bulb temperature for each hour of the
    given design day based on the correlation of Erbs et al. (1983).
    Usage: for determination of peak loads and for selecting HVAC equipment.
    """
    def __init__(
        self,
        location: Location,
        date: Date,
        Tdb_max: Quantity,
        Tdb_range: Quantity,
        Twb_mc: Quantity
    ):
        """
        Create `DesignDayTemperatureProfile` object.

        Parameters
        ----------
        location :
            Geographic location for which the temperature profiles are to be
            calculated.
        date :
            Date of the design day.
        Tdb_max :
            Design value of outdoor air dry-bulb temperature.
        Tdb_range :
            Mean difference between daily maximum and minimum temperature.
        Twb_mc :
            Mean coincident wet-bulb temperature (i.e. the average value of the
            wet-bulb temperature at the design dry-bulb temperature).
        """
        W_des = HumidAir(Tdb=Tdb_max, Twb=Twb_mc).W
        self._design_day_T_profiles = []  # (*)
        Tdb_range = Tdb_range.to('delta_degC').m
        k = Tdb_max.to('degC').m - Tdb_range / 2.0
        # calculate db- and wb-temperature at each hour of the design day
        for t_hr in range(25):
            # calculate solar time
            if t_hr == 24:
                date_time = DateTime.combine(date, Time(23, 59, 59))
            else:
                date_time = DateTime.combine(date, Time(t_hr, 0, 0))
            t_sol = time_to_decimal_hour(location.solar_time(date_time))
            # calculate db-temperature at t_sol
            tau = 2.0 * math.pi * (t_sol - 1) / 24.0
            s = 0.4632 * math.cos(tau - 3.8051)
            s += 0.0984 * math.cos(2 * tau - 0.360)
            s += 0.0168 * math.cos(3 * tau - 0.822)
            s += 0.0138 * math.cos(4 * tau - 3.513)
            Tdb = Q_(k + Tdb_range * s, 'degC')
            # calculate wb-temperature at t_sol
            W_sat = HumidAir(Tdb=Tdb, RH=Q_(100, 'pct')).W
            W = min(W_des, W_sat)
            Twb = HumidAir(Tdb=Tdb, W=W).Twb
            # add db- and wb-temperature to design day temperature profile
            self._design_day_T_profiles.append((date_time, Tdb, Twb))

        # (*) List of tuples with 3 elements: datetime, dry-bulb temperature
        # *T_db*, and wet-bulb temperature *T_wb*.

    @property
    def Tdb_profile(self) -> Dict[str, List[Any]]:
        """
        Get daily profile of dry-bulb temperature at design day.

        Returns
        -------
        Pandas DataFrame of time moments 't' (`datetime` objects) and
        corresponding list of dry-bulb temperatures 'T' (`Quantity` objects).
        """
        d = {
            't': [tup[0] for tup in self._design_day_T_profiles],
            'T': [tup[1] for tup in self._design_day_T_profiles]
        }
        return d

    @property
    def Twb_profile(self) -> Dict[str, List[Any]]:
        """
        Get daily profile of wet-bulb temperature at design day.

        Returns
        -------
        Dictionary with list of time moments 't' (`datetime` objects) and
        corresponding list of wet-bulb temperatures 'T' (`Quantity` objects).
        """
        d = {
            't': [tup[0] for tup in self._design_day_T_profiles],
            'T': [tup[2] for tup in self._design_day_T_profiles]
        }
        return d

    @staticmethod
    def _interpolate_T_profile(T_profile) -> Callable[[float], float]:
        t_ax = [time_to_decimal_hour(dt.time()) * 3600.0 for dt in T_profile['t']]
        T_ax = [T.to('degC').m for T in T_profile['T']]
        interpolant = interp1d(x=t_ax, y=T_ax, kind='cubic')
        return interpolant

    @property
    def Tdb(self) -> Callable[[float], float]:
        """
        Get function object that accepts time of day *t* in seconds from
        00:00:00 and returns the corresponding dry-bulb temperature at that
        moment in degC.
        """
        return self._interpolate_T_profile(self.Tdb_profile)

    @property
    def Twb(self) -> Callable[[float], float]:
        """
        Get function object that accepts time of day in seconds from 00:00:00
        and returns the corresponding wet-bulb temperature at that moment in degC.
        """
        return self._interpolate_T_profile(self.Twb_profile)


class ClimateData:
    """
    Class that holds design data about the outdoor climate at the given location
    and calculates from this design data:
    (1) the daily hourly average dry-bulb and wet-bulb temperature based on a
    correlation implemented in `DesignDayTemperatureProfile`, and
    (2) the daily hourly average solar irradiance on the horizontal surface
    based on the clear-sky model.
    """

    def __init__(self):
        self.design_day: Optional[Date] = None
        self.location: Optional[Location] = None
        self.Tdb_avg: Optional[Quantity] = None
        self.Tdb_range: Optional[Quantity] = None
        self.Twb_mc: Optional[Quantity] = None
        self.tau_beam: Optional[float] = None
        self.tau_dif: Optional[float] = None
        self.rho: Optional[Quantity] = None
        self.cp: Optional[Quantity] = None
        self._design_day_temp_profile: Optional[_DesignDayTemperatureProfile] = None
        self._clear_sky_model: Optional[ClearSkyModel] = None

    @classmethod
    def create(
        cls,
        design_day: Date,
        location: Location,
        Tdb_avg: Quantity,
        Tdb_range: Quantity,
        Twb_mc: Quantity,
        tau_beam: float,
        tau_dif: float
    ) -> 'ClimateData':
        """
        Create `ClimateData` object.

        Parameters
        ----------
        design_day:
            Date for which the heat transfer through the exterior building
            element is to be determined.
        location :
            The geographic location for which the climate data applies to.
        Tdb_avg:
            Outside air dry-bulb design temperature.
        Tdb_range:
            Outside air dry-bulb temperature range.
        Twb_mc:
            Outside air wet-bulb design temperature (mean coincident wet-bulb
            temperature).
        tau_beam :
            Clear-sky optical depth of the atmosphere at the given location for
            solar beam radiation on the given design day. Used to calculate
            clear-sky solar irradiance.
        tau_dif :
            Clear-sky optical depth of the atmosphere at the given location for
            solar diffuse radiation on the given design day. Used to calculate
            clear-sky solar irradiance.
        """
        obj = cls()
        obj.design_day = design_day
        obj.location = location
        obj.Tdb_avg = Tdb_avg
        obj.Tdb_range = Tdb_range
        obj.Twb_mc = Twb_mc
        obj.tau_beam = tau_beam
        obj.tau_dif = tau_dif
        outdoor_air = HumidAir(Tdb=obj.Tdb_avg, Twb=obj.Twb_mc)
        obj.rho = outdoor_air.rho
        obj.cp = outdoor_air.cp
        # to get hourly db- and wb temperature profile on design day
        obj._design_day_temp_profile = _DesignDayTemperatureProfile(
            location=obj.location,
            date=obj.design_day,
            Tdb_max=obj.Tdb_avg,
            Tdb_range=obj.Tdb_range,
            Twb_mc=obj.Twb_mc
        )
        # to get hourly clear sky irradiance profile on horizontal surface
        obj._clear_sky_model = ClearSkyModel(
            location=obj.location,
            tau_beam=obj.tau_beam,
            tau_dif=obj.tau_dif
        )
        return obj

    def get_data_members(self) -> Dict[str, Any]:
        return {
            'design_day': self.design_day,
            'location': self.location,
            'T_db': self.Tdb_avg,
            'db_range': self.Tdb_range,
            'T_wb_mc': self.Twb_mc,
            'tau_beam': self.tau_beam,
            'tau_dif': self.tau_dif
        }

    def Tdb_ext(self, t: float) -> float:
        """
        Get outdoor air dry-bulb temperature in degC at time *t* in seconds
        from 00:00:00.
        """
        return float(self._design_day_temp_profile.Tdb(t))

    def Twb_ext(self, t: float) -> float:
        """
        Get the outdoor air wet-bulb temperature in degC at time *t* in seconds
        from 00:00:00.
        """
        return float(self._design_day_temp_profile.Twb(t))

    @property
    def Tdb_profile(self) -> Dict[str, List[Any]]:
        """
        Get daily profile of dry-bulb temperature at design day.

        Returns
        -------
        Dictionary with list of time moments 't' (`datetime` objects) and
        corresponding list of dry-bulb temperatures 'T' (`Quantity` objects).
        """
        return self._design_day_temp_profile.Tdb_profile

    @property
    def Twb_profile(self) -> Dict[str, List[Any]]:
        """
        Get daily profile of wet-bulb temperature at design day.

        Returns
        -------
        Dictionary with list of time moments 't' (`datetime` objects) and
        corresponding list of wet-bulb temperatures 'T' (`Quantity` objects).
        """
        return self._design_day_temp_profile.Twb_profile

    # noinspection PyTypeChecker
    @property
    def irr_profile(self) -> Dict[str, List[Any]]:
        """
        Get daily profile of solar irradiance based on ASHRAE's Clear-Sky Model.

        Returns
        -------
        Dictionary of time moments 't' (`datetime` objects) and
        corresponding lists of irradiance components: 'beam', 'dir_hor',
        'dif', and 'glo_hor' (`Quantity` objects).
        """
        return self._clear_sky_model.get_daily_profile(self.design_day)

    @property
    def sun_position_profile(self) -> Dict[str, List[Any]]:
        """
        Get daily profile of sun positions (for each hour of the day).

        Returns
        -------
        Dictionary {
            't': list of `datetime` objects
            'azimuth': list of azimuth angles of the sun (`Quantity` objects)
            'elevation': list of elevation angles of the sun (`Quantity` objects)
        }
        """
        sun_path = SunPath(self.location, self.design_day)
        t_ax, azi_ax, elev_ax = sun_path.axes
        return {
            't': t_ax,
            'azimuth': azi_ax,
            'elevation': elev_ax
        }
