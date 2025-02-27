from __future__ import annotations
from typing import Callable
from pathlib import Path
import math
from datetime import date as Date
from scipy.interpolate import interp1d
from hvac import Quantity
from hvac.fluids import HumidAir
from hvac.sun import tmy, correlation, clear_sky, Location, Surface
from hvac.sun.time import decimal_hour_to_time
from hvac.charts.matplotlibwrapper import LineChart


Q_ = Quantity
PathLike = Path | str
SolarRadiation = (
    tmy.SolarRadiation |
    correlation.SolarRadiation |
    clear_sky.SolarRadiation
)


daily_temp_rng_fractions = [
    0.82, 0.88, 0.92, 0.95, 0.98, 1.00,
    0.98, 0.91, 0.74, 0.55, 0.38, 0.23,
    0.13, 0.05, 0.00, 0.00, 0.06, 0.14,
    0.24, 0.39, 0.50, 0.59, 0.68, 0.75
]
# from ASHRAE Handbook-Fundamentals 2017, Chapter 14, Table 6


def dry_bulb_temperature(
    t_sol_dec: float,
    T_db_des: Quantity,
    T_db_rng: Quantity
) -> Quantity:
    """Returns the dry-bulb temperature at `t_sol_dec`, calculated according to
    ASHRAE Handbook-Fundamentals 2017, Chapter 14, §3.1, Temperatures.

    Parameters
    ----------
    t_sol_dec:
        Solar time in decimal hours.
    T_db_des:
        Dry-bulb design temperature.
    T_db_rng:
        Mean coincident daily dry-bulb temperature range.
    """
    T_db_des = T_db_des.to('degC').magnitude
    T_db_rng = T_db_rng.to('delta_degC').magnitude
    t_sol_hr = int(round(t_sol_dec))
    frac = daily_temp_rng_fractions[t_sol_hr]
    T_db = T_db_des - frac * T_db_rng
    return Q_(T_db, 'degC')


def wet_bulb_temperature(
    t_sol_dec: float,
    T_wb_mc: Quantity,
    T_wb_rng: Quantity
) -> Quantity:
    """Returns the wet-bulb temperature at `t_sol_dec`, calculated according to
    ASHRAE Handbook-Fundamentals 2017, Chapter 14, §3.1, Temperatures.

    Parameters
    ----------
    t_sol_dec:
        Solar time in decimal hours.
    T_wb_mc:
        The mean coincident wet-bulb temperature that goes with `T_db_des`.
    T_wb_rng:
        Mean coincident daily wet-bulb temperature range.
    """
    T_wb_mc = T_wb_mc.to('degC').magnitude
    T_wb_rng = T_wb_rng.to('delta_degC').magnitude
    t_sol_hr = int(round(t_sol_dec))
    frac = daily_temp_rng_fractions[t_sol_hr]
    T_wb = T_wb_mc - frac * T_wb_rng
    return Q_(T_wb, 'degC')


def synthetic_daily_temperature_profiles(
    T_db_des: Quantity,
    T_wb_mc: Quantity,
    T_db_rng: Quantity,
    T_wb_rng: Quantity
) -> tuple[list[Quantity], list[Quantity]]:
    """Returns a list with the dry-bulb temperatures and a list with the
    wet-bulb temperatures for each hour (0..23 h) of the selected day,
    calculated according to ASHRAE Handbook-Fundamentals 2017, Chapter 14, §3.1,
    Temperatures.

    Parameters
    ----------
    T_db_des:
        Dry-bulb design temperature.
    T_wb_mc:
        The mean coincident wet-bulb temperature that goes with `T_db_des`.
    T_db_rng:
        Mean coincident daily dry-bulb temperature range.
    T_wb_rng:
        Mean coincident daily wet-bulb temperature range.
    """
    T_db_prof = [
        dry_bulb_temperature(t_sol_dec, T_db_des, T_db_rng)
        for t_sol_dec in range(0, 24)
    ]
    T_wb_prof = [
        wet_bulb_temperature(t_sol_dec, T_wb_mc, T_wb_rng)
        for t_sol_dec in range(0, 24)
    ]
    return T_db_prof, T_wb_prof


class WeatherData:
    """Contains routines to create hourly outdoor temperature profiles (dry-bulb
    and wet-bulb) and hourly solar radiation profiles on a specified day of the
    year (i.e. the "design day").
    
    The class is also used to determine the solar radiation on exterior surfaces
    and through windows.

    Depending on the constructor method and its parameters, there three possible
    ways to create a `WeatherData` object.
    
    (1) From TMY-data.
    If the constructor method `create_from_tmy_data` is used, the solar radiation 
    calculations are done by an instance of class `sun.tmy.SolarRadiation` 
    using hourly TMY-data from a specified csv-file.
    
    (2) From climatic design data tables.
    If the constructor method `create_from_climatic_design_data` is used, the 
    solar radiation calculations can be done either by an instance of class
    `sun.correlation.SolarRadiation` or by an instance of class 
    `sun.clear_sky.SolarRadiation`.
    
    (2.a) `sun.correlation.SolarRadiation` 
    Is used when either the monthly average of the daily irradiation on a 
    horizontal surface `H_avg` is specified, or the monthly average of the daily
    clearness-index `K_T_avg`. 
    Solar radiation at a given time moment is then determined using a 
    statistical model. However, for a given monthly average, multiple 
    distributions of the solar radiation on a given day are possible. Therefore,
    different results may appear for the same day of the month.
    
    (2.b) `sun.clear_sky.SolarRadiation`
    If neither `H_avg`, nor `K_T_avg` is specified, the clear-sky model will be 
    used, which is implemented in class `sun.clear_sky.SolarRadiation`.
    """

    def __init__(self) -> None:
        self.location: Location | None = None
        self.sol_rad: SolarRadiation | None = None
        self.T_db_prof: list[Quantity] | None = None
        self.T_wb_prof: list[Quantity] | None = None
        self._T_db_interp: Callable[[float], float] | None = None
        self._T_wb_interp: Callable[[float], float] | None = None
        self.T_db_max: Quantity | None = None
        self.T_db_avg: Quantity | None = None

    @classmethod
    def create_from_tmy_data(
        cls,
        location: Location,
        date: Date,
        tmy_file: PathLike,
        tmy_params: tuple[str, str] = ('%Y%m%d:%H%M', 'UTC')
    ) -> WeatherData:
        """Creates a `WeatherData` object from TMY-data saved in a csv-file.

        Parameters
        ----------
        location:
            Geographic location for which the weather data is valid.
        date:
            Date for which the weather data is valid.
        tmy_file:
            Path to the csv-file with hourly TMY data.
        tmy_params:
            2-tuple with the date-time format used in the csv-file and the
            timezone in which the date-times are expressed in the csv-file,
            using the tz-database notation.
            (https://en.wikipedia.org/wiki/List_of_tz_database_time_zones).
        """
        obj = cls()
        obj.location = location
        obj.location.date = date
        obj.sol_rad = tmy.SolarRadiation(location, tmy_file, tmy_params)
        df_day = obj.sol_rad.tmy.get_day(date.month, date.day)
        obj.T_db_prof = [
            Q_(T_db, 'degC')
            for T_db in df_day['T2m'].values
        ]
        RH_prof = [
            Q_(RH, 'pct')
            for RH in df_day['RH'].values
        ]
        obj.T_wb_prof = [
            HumidAir(Tdb=T_db, RH=RH).Twb.to('degC')
            for T_db, RH in zip(obj.T_db_prof, RH_prof)
        ]
        interpol_fun = obj._interpolate_profiles()
        obj._T_db_interp = interpol_fun[0]
        obj._T_wb_interp = interpol_fun[1]
        obj.T_db_max = max(obj.T_db_prof)
        obj.T_db_avg = sum(T_db.to('K') for T_db in obj.T_db_prof) / len(obj.T_db_prof)
        return obj

    @classmethod
    def create_from_climatic_design_data(
        cls,
        location: Location,
        date: Date,
        T_db_des: Quantity,
        T_db_rng: Quantity,
        T_wb_mc: Quantity,
        T_wb_rng: Quantity,
        H_avg: Quantity | None = None,
        K_T_avg: float | None = None
    ) -> WeatherData:
        """Creates a `WeatherData` object from information that can be looked up 
        in e.g. ASHRAE's climatic design data tables.

        Parameters
        ----------
        location:
            Geographic location for which the weather data is valid.
        date:
            Date for which the weather data is valid, i.e. the "design day".
        T_db_des:
            Monthly design dry-bulb temperature, i.e., the maximum temperature
            on the selected day.
        T_db_rng:
            Mean daily temperature range, i.e., the difference between the
            maximum and the minimum temperature on the selected day.
        T_wb_mc:
            Mean coincident wet-bulb temperature.
        T_wb_rng:
            Mean coincident daily wet-bulb temperature range.
        H_avg: optional
            Monthly average daily global radiation on a horizontal surface.
        K_T_avg: optional
            Monthly average daily clearness index.

        Notes
        -----
        1. If neither `H_avg`, nor `K_T_avg` is specified, the clear-sky model
        will be used to determine solar radiation incident on surfaces.

        2. Hourly values of dry-bulb and wet-bulb temperature are determined
        with the correlation of Erbs et al. (1983).

        3. In case `H_avg` or `K_T_avg` is specified, correlations are used to
        determine a daily distribution of the solar radiation on the specified
        `date`. For a given monthly average of the daily irradiation `H_avg` or
        a monthly average of the daily clearness index `K_T_avg` multiple
        distributions of daily solar radiation are possible. This means that
        for a given `date` different results may appear each time the code is
        called.

        References
        ----------
        ASHRAE Handbook-Fundamentals 2017, Chapter 14 Climatic Design
        Information.
        """
        obj = cls()
        obj.location = location
        obj.location.date = date
        if H_avg is not None or K_T_avg is not None:
            if H_avg is not None: K_T_avg = float(H_avg / location.sun.H_o_avg)
            obj.sol_rad = correlation.SolarRadiation(location, date.month, K_T_avg)
        else:
            obj.sol_rad = clear_sky.SolarRadiation(location)
        obj.T_db_prof, obj.T_wb_prof = synthetic_daily_temperature_profiles(
            T_db_des=T_db_des,
            T_wb_mc=T_wb_mc,
            T_db_rng=T_db_rng,
            T_wb_rng=T_wb_rng
        )
        interpol_fun = obj._interpolate_profiles()
        obj._T_db_interp = interpol_fun[0]
        obj._T_wb_interp = interpol_fun[1]
        obj.T_db_max = max(obj.T_db_prof)
        obj.T_db_avg = sum(T_db.to('K') for T_db in obj.T_db_prof) / len(obj.T_db_prof)
        return obj

    @property
    def date(self) -> Date:
        return self.location.date

    def plot_daily_temperature_profile(self) -> None:
        """Shows a diagram with the hourly values (in solar time) of the dry-bulb
        and the wet-bulb temperature on the date (the design day) indicated when
        instantiating the `WeatherData` object.
        """
        chart = LineChart()
        chart.add_xy_data(
            label='dry-bulb temperature',
            x1_values=[hr for hr in range(24)],
            y1_values=[T_db.to('degC').m for T_db in self.T_db_prof],
            style_props={'marker': 'o', 'linestyle': '--'}
        )
        chart.add_xy_data(
            label='wet-bulb temperature',
            x1_values=[hr for hr in range(24)],
            y1_values=[T_wb.to('degC').m for T_wb in self.T_wb_prof],
            style_props={'marker': 'o', 'linestyle': '--'}
        )
        chart.x1.add_title('solar hour')
        chart.y1.add_title('temperature, °C')
        chart.add_legend()
        chart.show()

    def plot_daily_solar_radiation_profile(self) -> None:
        """Shows a diagram with the hourly values (in solar time) of the total
        irradiance or irradiation on a horizontal surface on the date
        (the design day) indicated when instantiating the `WeatherData` object.

        In case the `WeatherData` object was instantiated by specifying `H_avg`
        or `K_T_avg`, an instance of `sun.correlation.SolarRadiation` is used
        to perform the solar radiation calculations. This instance will
        return a bar chart of the hourly irradiation components (beam, diffuse
        and total) on the horizontal surface. `sun.tmy.SolarRadiation` and
        `sun.clear_sky.SolarRadiation` will return a line chart of the
        irradiance on a horizontal surface.
        """
        if isinstance(self.sol_rad, (tmy.SolarRadiation, clear_sky.SolarRadiation)):
            self.sol_rad.plot_irradiance(self.date.month, self.date.day)
        elif isinstance(self.sol_rad, correlation.SolarRadiation):
            self.sol_rad.plot_hourly_irradiation(self.date.day)

    def _interpolate_profiles(self):
        """Creates interpolation functions from the dry-bulb temperature profile 
        and the wet-bulb temperature profile on the design day.
        
        This allows to determine by linear interpolation the dry-bulb 
        temperature and the wet-bulb temperature at any time moment during the 
        design day.
        
        The interpolation functions expect to get solar time expressed in 
        decimal hours (float) and return temperature values (float) in units of
        degC.
        """
        x = [hr for hr in range(25)]
        y1 = [T_db.to('degC').m for T_db in self.T_db_prof]
        y1.append(y1[0])  # 24 h coincides with 0 h
        y2 = [T_wb.to('degC').m for T_wb in self.T_wb_prof]
        y2.append(y2[0])  # 24 h coincides with 0 h
        T_db_interp = interp1d(x, y1, kind='slinear')
        T_wb_interp = interp1d(x, y2, kind='slinear')
        return T_db_interp, T_wb_interp

    def T_db(self, t_sol_sec: float) -> Quantity:
        """Returns the outdoor air dry-bulb temperature at the specified solar 
        time of the design day.

        Parameters
        ----------
        t_sol_sec: float
            Solar time in seconds from midnight (= 0 s).
        """
        return Q_(self._T_db_interp(t_sol_sec / 3600), 'degC')

    def T_wb(self, t_sol_sec: float) -> Quantity:
        """Returns the outdoor air wet-bulb temperature at the specified solar 
        time of the design day.

        Parameters
        ----------
        t_sol_sec: float
            Solar time in seconds from midnight (= 0 s).
        """
        return Q_(self._T_wb_interp(t_sol_sec / 3600), 'degC')

    def I_T(
        self,
        hr1: int,
        surf: Surface,
        rho_g: Quantity = Q_(0.2, 'frac'),
        sky_model: str = ''
    ) -> tuple[Quantity, Quantity, Quantity] | None:
        """Returns the hourly solar irradiation components on a tilted surface 
        at the given solar time hour `hr1`, which specifies the start time of a 
        one hour-period.

        Parameters
        ----------
        hr1: int
            start hour in solar time of the 1-hour period
        surf: Surface
            tilted surface.
        rho_g: Quantity, default 0.2 frac
            ground reflectance
        sky_model: str, {'anisotropic.hdkr, 'isotropic'}, optional
            sky model to be used for calculating the irradiation on the tilted
            surface. By default, the anisotropic sky model according to Perez
            is used.

        Returns
        -------
        3-tuple:
        - hourly total irradiation
        - hourly beam irradiation
        - hourly diffuse sky and ground reflected irradiation
        """
        if isinstance(self.sol_rad, (tmy.SolarRadiation, clear_sky.SolarRadiation)):
            hic = self.sol_rad.tilted_hourly_irradiation(
                month=self.date.month,
                day=self.date.day,
                hr1=hr1,
                surf=surf,
                rho_g=rho_g,
                sky_model=sky_model
            )
            I_d = hic.I_d + hic.I_g  # sky diffuse and ground reflected irradiation
            return hic.I, hic.I_b, I_d
        elif isinstance(self.sol_rad, correlation.SolarRadiation):
            hic = self.sol_rad.tilted_hourly_irradiation(
                day=self.date.day,  # note: the month was already set on instantiation of `correlation.SolarRadiation`
                hr1=hr1,
                surf=surf,
                rho_g=rho_g,
                sky_model=sky_model
            )
            I_d = hic.I_d + hic.I_g
            return hic.I, hic.I_b, I_d


def sol_air_temperature(
    G_t: Quantity,
    T_db: Quantity,
    tilt_angle: Quantity,
    surface_color: str = 'dark',
    R_surf: Quantity | None = None,
    a_surf: Quantity | None = None
) -> Quantity:
    """Returns the sol-air temperature of an exterior surface.

    Parameters
    ----------
    G_t:
        Total irradiance on the surface.
    T_db:
        Outdoor air dry-bulb temperature.
    tilt_angle:
        The tilt or slope angle of the surface.
    surface_color: optional
        Either a 'dark' (default) or 'light' colored surface.
    R_surf: optional
        Thermal resistance of the exterior surface film.
    a_surf: optional
        Absorption factor of the exterior surface.

    Notes
    -----
    If `R_surf` and `a_surf` are specified, parameter `surface_color` is ignored.
    """
    G_t = G_t.to('W / m**2').m
    T_db = T_db.to('degC').m
    tilt_angle = tilt_angle.to('rad').m
    R_surf = R_surf.to('K * m**2 / W').m if R_surf is not None else None
    a_surf = a_surf.to('frac').m if a_surf is not None else None
    dT_sky = 3.9 * math.cos(tilt_angle)
    if R_surf is not None and a_surf is not None:
        T_sa = T_db + a_surf * R_surf * G_t - dT_sky
    elif surface_color == 'light':
        T_sa = T_db + 0.026 * G_t - dT_sky
    else:
        T_sa = T_db + 0.052 * G_t - dT_sky
    return Q_(T_sa, 'degC')


class ExteriorSurface(Surface):
    """Represents the exterior surface of opaque exterior building elements 
    and windows. This class is used internally inside the class 
    `exterior_building_element.ExteriorBuildingElement` and inside the class 
    `fenestration.Window`. It contains methods needed to calculate the exterior
    temperature and solar radiation incident on the surface at a given time 
    moment of the design day.
    """
    def __init__(
        self,
        weather_data: WeatherData,
        gamma: Quantity,
        beta: Quantity,
        surface_color: str = 'dark',
        R_surf: Quantity | None = None,
        a_surf: Quantity | None = None,
        rho_g: Quantity = Q_(0.2, 'frac'),
        sky_model: str = ''
    ) -> None:
        """Creates an `ExteriorSurface` object.

        Parameters
        ----------
        weather_data:
            Instance of class `WeatherData` responsible for delivering all
            the climatic design data (temperature and solar radiation) on the
            specified location and date.
        gamma:
            Azimuth angle of the surface. South = 0°, West = +90°, East = -90°.
        beta:
            Slope angle of the surface, i.e. the angle between the horizontal
            plane and the surface.
        surface_color: optional
            Either a 'dark' (default) or 'light' colored surface.
        R_surf: optional
            Thermal resistance of the exterior surface film.
        a_surf: optional
            Absorption factor of the exterior surface.
        rho_g: Quantity, default 0.2 frac
            Ground reflectance.
        sky_model: str, {'anisotropic.hdkr, 'isotropic'}, optional
            Sky model to be used for calculating the irradiation on the tilted
            surface. By default, the anisotropic sky model according to Perez
            is used.

        Notes
        -----
        If `R_surf` and `a_surf` are specified, parameter `surface_color` is
        ignored.
        """
        super().__init__(weather_data.location, gamma, beta)
        self.weather_data = weather_data
        self.surface_color = surface_color
        self.R_surf = R_surf
        self.a_surf = a_surf
        self.rho_g = rho_g
        self.sky_model = sky_model
        interpol_fun = self._interpolate_profiles()
        self._T_sa_interp = interpol_fun[0]
        self._I_T_interp = interpol_fun[1]
        self._I_Tb_interp = interpol_fun[2]
        self._I_Td_interp = interpol_fun[3]

    def _interpolate_profiles(self):
        """Creates interpolation functions of the sol-air temperature, the total
        solar irradiation, the solar beam irradiation, and the diffuse solar
        irradiation on the exterior surface on the design day.

        The interpolation functions expect to retrieve solar time expressed in 
        decimal hours (float) and return temperature values (float) in units of 
        degC and solar irradiation values (float) in units of J/m².
        """
        x = [hr for hr in range(25)]
        y3, y4, y5, y6 = [], [], [], []
        for hr in x[:-1]:
            I_T, I_Tb, I_Td = self.weather_data.I_T(hr, self, self.rho_g, self.sky_model)
            G_T_avg = I_T.to('J / m**2') / Q_(3600, 's')
            T_sa = sol_air_temperature(
                G_t=G_T_avg,
                T_db=self.weather_data.T_db(hr * 3600.0),
                tilt_angle=self.beta,
                surface_color=self.surface_color,
                R_surf=self.R_surf,
                a_surf=self.a_surf
            )
            y3.append(T_sa.to('degC').m)
            y4.append(I_T.to('J / m**2').m)
            y5.append(I_Tb.to('J / m**2').m)
            y6.append(I_Td.to('J / m**2').m)
        y3.append(y3[0])  # 24 h coincides with 0 h
        y4.append(y4[0])
        y5.append(y5[0])
        y6.append(y6[0])
        T_sa_interp = interp1d(x, y3, kind='slinear')
        I_T_interp = interp1d(x, y4, kind='slinear')
        I_Tb_interp = interp1d(x, y5, kind='slinear')
        I_Td_interp = interp1d(x, y6, kind='slinear')
        return (
            T_sa_interp,
            I_T_interp, I_Tb_interp, I_Td_interp
        )

    def T_sa(self, t_sol_sec: float) -> Quantity:
        """Returns the sol-air temperature of the exterior surface at the
        specified solar time of the design day.

        Parameters
        ----------
        t_sol_sec: float
            Solar time in seconds from midnight (= 0 s).
        """
        return Q_(self._T_sa_interp(t_sol_sec / 3600), 'degC')

    def T_db(self, t_sol_sec: float) -> Quantity:
        """Returns the outdoor air dry-bulb temperature at the specified solar 
        time of the design day.

        Parameters
        ----------
        t_sol_sec: float
            Solar time in seconds from midnight (= 0 s).
        """
        return self.weather_data.T_db(t_sol_sec)

    def theta_i(self, t_sol_sec: float) -> Quantity:
        """Returns the solar incidence angle on the exterior surface at the
        specified solar time of the design day.

        Parameters
        ----------
        t_sol_sec: float
            Solar time in seconds from midnight (= 0 s).
        """
        t_sol_hr = t_sol_sec / 3600
        sol_time = decimal_hour_to_time(t_sol_hr)
        self.location.solar_time = sol_time
        return self.theta

    def alpha_s(self, t_sol_sec: float) -> Quantity:
        """Returns the solar altitude angle at the specified solar time of the
        design day.

        Parameters
        ----------
        t_sol_sec: float
            Solar time in seconds from midnight (= 0 s).
        """
        t_sol_hr = t_sol_sec / 3600
        sol_time = decimal_hour_to_time(t_sol_hr)
        self.location.solar_time = sol_time
        return self.location.sun.alpha

    def gamma_ss(self, t_sol_sec: float) -> Quantity:
        """Returns the sun-surface azimuth angle at the specified solar time of
        the design day.

        Parameters
        ----------
        t_sol_sec: float
            Solar time in seconds from midnight (= 0 s).
        """
        t_sol_hr = t_sol_sec / 3600
        sol_time = decimal_hour_to_time(t_sol_hr)
        self.location.solar_time = sol_time
        gamma_s = self.location.sun.gamma
        return gamma_s - self.gamma

    def G_T(self, t_sol_sec: float) -> tuple[Quantity, Quantity, Quantity]:
        """Returns the hourly average irradiance components on the exterior
        surface at the specified solar time of the design day.

        Parameters
        ----------
        t_sol_sec: float
            Solar time in seconds from midnight (= 0 s).

        Returns
        -------
        - hourly average of total irradiance
        - hourly average of beam irradiance
        - hourly average of diffuse sky and ground reflected irradiance
        """
        I_T = self._I_T_interp(t_sol_sec / 3600)
        I_Tb = self._I_Tb_interp(t_sol_sec / 3600)
        I_Td = self._I_Td_interp(t_sol_sec / 3600)
        G_T = Q_(I_T / 3600, 'W / m**2')
        G_Tb = Q_(I_Tb / 3600, 'W / m**2')
        G_Td = Q_(I_Td / 3600, 'W / m**2')
        return G_T, G_Tb, G_Td
