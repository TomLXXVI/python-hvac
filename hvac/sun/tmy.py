from datetime import date as Date
from datetime import time as Time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.integrate import quad
from hvac import Quantity
from hvac.charts.matplotlibwrapper import LineChart
from hvac.sun.time import (
    clock_to_solar_time,
    time_to_decimal_hour,
)
from hvac.sun.surface import Location, hour_angle, Surface
from hvac.sun import transposition
from hvac.sun.utilizability import Utilizability
from hvac.sun.clear_sky import IrradianceComponents, HourlyIrradiationComponents

Q_ = Quantity


class TMY:
    """Class for reading TMY datafiles."""

    def __init__(
        self,
        file_path: str,
        dt_fmt: str = '%Y%m%d:%H%M',
        tz: str = 'UTC',
    ) -> None:
        """
        Creates a `TMY` instance by reading the csv-file pointed to by
        `file_path` into a Pandas DataFrame object.

        The implementation of this class is based on the csv-file with TMY-data
        that is downloaded from PVGIS.
        (https://re.jrc.ec.europa.eu/pvg_tools/en/#api_5.2)
        This means that date and time are grouped in a single column, which is
        the first column of the csv-file. Also, the cvs-file must only contain
        the actual TMY-data (other rows or columns need to be deleted from the
        csv-file).

        Parameters
        ----------
        file_path:
            The path to the csv-file containing the TMY-data.
        dt_fmt: default '%Y%m%d:%H%M'
            The format in which date-times are written in the csv-file.
        tz: default 'UTC'
            The timezone of the date-times in the csv-file, expressed in
            tz-database notation.
            (https://en.wikipedia.org/wiki/List_of_tz_database_time_zones).
        """
        self.file_path = file_path
        self.dataframe = self._read_csv_file(dt_fmt, tz)

    def _read_csv_file(self, dt_fmt: str, tz: str) -> pd.DataFrame:
        df = pd.read_csv(
            self.file_path,
            parse_dates=[0],
            dayfirst=True,
            date_format=dt_fmt
        )
        df = df.set_index(df.columns[0])
        df.index = df.index.tz_localize(tz, nonexistent='shift_forward')
        return df

    @property
    def header(self) -> pd.Index:
        """Returns the column names of the csv-file."""
        return self.dataframe.columns

    def get_month(self, month: int) -> pd.DataFrame:
        """Returns only the TMY-data for the given month, indicated by its
        integer index.
        """
        months = self.dataframe.groupby(self.dataframe.index.month)
        month = months.get_group(month)
        return month

    def get_day(self, month: int, day: int) -> pd.DataFrame:
        """Returns only the TMY-data for the given month and day, both indicated
        by their integer index.
        """
        month = self.get_month(month)
        days = month.groupby(month.index.day)
        day = days.get_group(day)
        return day

    def convert_to_timezone(self, tz: str) -> None:
        """Convert the date-time index to time zone `tz`."""
        self.dataframe.index = self.dataframe.index.tz_convert(tz)
        self.dataframe.index.name = f'Time ({tz})'
        self.dataframe.sort_index(inplace=True)

    def convert_to_solar_time(self, tz: str, L_loc: float) -> None:
        """Convert the date-time index to local solar time. Parameter `tz` is
        the local timezone with respect to GMT or UTC (e.g. 'Etc/GMT-1') and is
        used to determine the local standard time meridian of the time zone.
        Do not specify a timezone with DST, should it be in use at the location
        under consideration, as this may lead to errors at times when DST becomes
        active or inactive.
        Longitude `L_loc` in decimal degrees must be positive east of UTC and
        negative west of UTC.
        """
        if self.dataframe.index.tz != tz:
            self.convert_to_timezone(tz)
        self.dataframe.index = pd.DatetimeIndex([
            clock_to_solar_time(dt_local, L_loc)
            for dt_local in self.dataframe.index
        ])
        self.dataframe.index.name = 'Time (solar time)'


class SolarRadiation:
    """Estimation of solar radiation on horizontal and tilted surfaces at a
    given location using TMY data.
    """

    def __init__(
        self,
        location: Location,
        tmy_file: str | Path,
        tmy_params: tuple[str, str] = ('%Y%m%d:%H%M', 'UTC')
    ) -> None:
        """Creates a `SolarRadiation` instance based on TMY-data.

        Parameters
        ----------
        location: Location
            Location where the TMY-data applies.
        tmy_file: str | Path
            File path to csv-file with TMY data.
        tmy_params: optional
            2-tuple with the format in which date-times are written in the
            csv-file (default value: '%Y%m%d:%H%M') and the timezone of the
            date-times in the csv-file, expressed in its tz-database notation
            (default value: 'UTC').
        """
        self.location = location
        self.tmy = TMY(tmy_file, dt_fmt=tmy_params[0], tz=tmy_params[1])
        self.tmy.convert_to_solar_time(
            tz=location.timezone,
            L_loc=location.L_loc.to('deg').m
        )

    def plot_irradiance(self, month: int, day: int) -> None:
        """Plot a curve for the given month and day of the total irradiance on
        the horizontal surface."""
        tmy_data = self.tmy.get_day(month, day)
        G_data = tmy_data['G(h)']  # W / m ** 2
        omega_data = [
            hour_angle(time_to_decimal_hour(dt.time()))
            for dt in G_data.index
        ]
        G_fun = interp1d(omega_data, G_data.values, kind='slinear')
        omega_interp = np.linspace(omega_data[0], omega_data[-1], num=100)
        G_interp = G_fun(omega_interp)
        chart = LineChart()
        chart.add_xy_data(
            label='data',
            x1_values=omega_data,
            y1_values=G_data,
            style_props={'marker': 'o', 'linestyle': 'none'}
        )
        chart.add_xy_data(
            label='interpolated',
            x1_values=omega_interp,
            y1_values=G_interp,
            style_props={'linestyle': '--'}
        )
        chart.x1.add_title('hour angle, rad')
        chart.y1.add_title('irradiance, W/m²')
        chart.show()

    def tilted_hourly_irradiance_array(
        self,
        month: int,
        day: int,
        surf: Surface,
        rho_g: Quantity = Q_(0.2, 'frac'),
        sky_model: str = ''
    ) -> IrradianceComponents:
        """Returns the irradiance components on the tilted surface for each hour
        of the given day.

        Parameters
        ----------
        month: int
            numerical index of the month of the year
        day: int
            numerical index of the day of the month
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
        IrradianceComponents, a named tuple with fields:
        - G: `Quantity` array of the global irradiance on the tilted surface
             for each hour of the day (0..23).
        - G_b: Same, but with beam irradiance values instead.
        - G_d: Same, but with diffuse irradiance values instead.
        - G_g: Same, but with ground reflected irradiance values instead.
        """
        self.location.date = Date(2001, month, day)
        tmy_data = self.tmy.get_day(month, day)
        G_arr = np.array(tmy_data['G(h)'].values)
        G_d_arr = np.array(tmy_data['Gd(h)'].values)
        G_b_arr = G_arr - G_d_arr
        G_d_arr = Q_(G_d_arr, 'W / m ** 2')
        G_b_arr = Q_(G_b_arr, 'W / m ** 2')
        G_T_arr, G_Tb_arr, G_Td_arr, G_Tg_arr = [], [], [], []
        if sky_model.lower() == 'anisotropic.hdkr':
            SkyModel = transposition.AnisotropicSkyModel.HDKR
        elif sky_model.lower() == 'isotropic':
            SkyModel = transposition.IsotropicSkyModel
        else:
            SkyModel = transposition.AnisotropicSkyModel.Perez
        for G_b, G_d in zip(G_b_arr, G_d_arr):
            if G_b > 0 and G_d > 0:
                sky_model = SkyModel(
                    surface=surf,
                    rho_g=rho_g,
                    G_b=G_b,
                    G_d=G_d
                )
                G_T_arr.append(sky_model.G_T.to('W / m ** 2').m)
                G_Tb_arr.append(sky_model.G_Tb.to('W / m ** 2').m)
                G_Td_arr.append(sky_model.G_Td.to('W / m ** 2').m)
                G_Tg_arr.append(sky_model.G_Tg.to('W / m ** 2').m)
            else:
                G_T_arr.append(0.0)
                G_Tb_arr.append(0.0)
                G_Td_arr.append(0.0)
                G_Tg_arr.append(0.0)
        return IrradianceComponents(
            G=Q_(G_T_arr, 'W / m ** 2'),
            G_b=Q_(G_Tb_arr, 'W / m ** 2'),
            G_d=Q_(G_Td_arr, 'W / m ** 2'),
            G_g=Q_(G_Tg_arr, 'W / m ** 2')
        )

    def hourly_irradiation(
        self,
        month: int,
        day: int,
        hr1: int
    ) -> HourlyIrradiationComponents:
        """Returns the hourly irradiation components on a horizontal surface for
        the given hour on the given day in the given month. `hr1` specifies the
        start time of the 1 hour-period, assuming solar time.

        Returns
        -------
        HourlyIrradiationComponents, a named tuple with fields:
        - I: total irradiation on a horizontal surface at hour `hr1`
        - I_b: beam irradiation
        - I_d: sky-diffuse irradiation
        """
        self.location.date = Date(2001, month, day)
        tmy_data = self.tmy.get_day(month, day)
        G_data = np.array(tmy_data['G(h)'].values)  # W / m ** 2
        G_d_data = np.array(tmy_data['Gd(h)'].values)
        G_b_data = G_data - G_d_data
        omega_sr = self.location.omega_sr.to('rad').m
        omega_ss = self.location.omega_ss.to('rad').m
        omega_data = [
            hour_angle(time_to_decimal_hour(dt.time()))
            for dt in tmy_data.index
        ]
        G_interp = interp1d(omega_data, G_data, kind='slinear')
        G_b_interp = interp1d(omega_data, G_b_data, kind='slinear')
        G_d_interp = interp1d(omega_data, G_d_data, kind='slinear')
        omega_1 = hour_angle(hr1)
        omega_2 = hour_angle(hr1 + 1)
        if (omega_2 <= omega_sr) or (omega_1 >= omega_ss):
            return HourlyIrradiationComponents(
                I=Q_(0.0, 'J / m ** 2'),
                I_b=Q_(0.0, 'J / m ** 2'),
                I_d=Q_(0.0, 'J / m ** 2'),
                I_g=None
            )
        if omega_1 < omega_sr < omega_2:
            omega_1 = omega_sr
        if omega_1 < omega_ss < omega_2:
            omega_2 = omega_ss
        I = quad(lambda omega: G_interp(omega), omega_1, omega_2)[0]  # W * rad / m ** 2
        I_b = quad(lambda omega: G_b_interp(omega), omega_1, omega_2)[0]
        I_d = quad(lambda omega: G_d_interp(omega), omega_1, omega_2)[0]
        k = 1.0 / np.radians(15.0)  # hr / rad
        return HourlyIrradiationComponents(
            I=Q_(k * I, 'W * hr / m ** 2').to('J / m ** 2'),
            I_b=Q_(k * I_b, 'W * hr / m ** 2').to('J / m ** 2'),
            I_d=Q_(k * I_d, 'W * hr / m ** 2').to('J / m ** 2'),
            I_g=None
        )

    def avg_hourly_irradiation(self, month: int, hr1: int) -> HourlyIrradiationComponents:
        """Returns the monthly averages of the hourly irradiation components on
        a horizontal surface at the given hour `hr1` and the arrays of the hourly
        irradiation components on a horizontal surface at the same hour for
        every day of the given month.

        Returns
        -------
        HourlyIrradiationComponents, a named tuple with fields:
        - I: 2-tuple with the monthly average of total irradiation at `hr1` and
             a`Quantity`-array of the daily values of total irradiation at `hr1`.
        - I_b: same for beam irradiation
        - I_d: same for sky-diffuse irradiation
        """
        tmy_data = self.tmy.get_month(month)
        num_days = tmy_data.index.day.nunique()
        hic_lst = [
            self.hourly_irradiation(month, d, hr1) for d in range(1, num_days + 1)
        ]
        I_arr = np.array([hic.I.to('J / m**2').m for hic in hic_lst])
        I_avg = np.mean(I_arr)
        I_b_arr = np.array([hic.I_b.to('J / m**2').m for hic in hic_lst])
        I_b_avg = np.mean(I_b_arr)
        I_d_arr = np.array([hic.I_d.to('J / m**2').m for hic in hic_lst])
        I_d_avg = np.mean(I_d_arr)
        return HourlyIrradiationComponents(
            I=(Q_(I_avg, 'J / m ** 2'), Q_(I_arr, 'J / m ** 2')),
            I_b=(Q_(I_b_avg, 'J / m ** 2'), Q_(I_b_arr, 'J / m ** 2')),
            I_d=(Q_(I_d_avg, 'J / m ** 2'), Q_(I_d_arr, 'J / m ** 2')),
            I_g=None
        )

    def tilted_hourly_irradiation(
        self,
        month: int,
        day: int,
        hr1: int,
        surf: Surface,
        rho_g: Quantity = Q_(0.2, 'frac'),
        sky_model: str = ''
    ) -> Quantity:
        """Returns the hourly irradiation components on the tilted surface for
        the given hour on the given day in the given month. `hr1` specifies the
        start time of the 1 hour-period, assuming solar time.

        Parameters
        ----------
        month: int
            numerical index of the month of the year
        day: int
            numerical index of the day of the month
        hr1: int
            start hour (solar time) of the 1-hour period
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
        HourlyIrradiationComponents, a named tuple with fields:
        - I: total irradiation on the tilted surface at `hr1`
        - I_b: beam irradiation
        - I_d: sky-diffuse irradiation
        - I_g: ground reflected irradiation
        """
        I, I_b, I_d, _ = self.hourly_irradiation(month, day, hr1)
        if I > 0:
            hr2 = hr1 + 1
            t1, t2 = Time(hr1), Time(hr2)
            I_b, I_d = surf.location.sun.split_radiation(I=I, time_interval=(t1, t2))
            if sky_model.lower() == 'anisotropic.hdkr':
                SkyModel = transposition.AnisotropicSkyModel.HDKR
            elif sky_model.lower() == 'isotropic':
                SkyModel = transposition.IsotropicSkyModel
            else:
                SkyModel = transposition.AnisotropicSkyModel.Perez
            sky_model = SkyModel(
                surface=surf,
                rho_g=rho_g,
                G_b=I_b,
                G_d=I_d,
                time_interval=(t1, t2)
            )
            I_T = sky_model.G_T.to('J / m ** 2')
            I_Tb = sky_model.G_Tb.to('J / m**2')
            I_Td = sky_model.G_Td.to('J / m**2')
            I_Tg = sky_model.G_Tg.to('J / m**2')
            return HourlyIrradiationComponents(
                I=I_T,
                I_b=I_Tb,
                I_d=I_Td,
                I_g=I_Tg
            )
        else:
            return HourlyIrradiationComponents(
                *(Q_(0.0, 'J / m ** 2') for _ in range(4))
            )

    def tilted_avg_hourly_irradiation(
        self,
        month: int,
        hr1: int,
        surf: Surface,
        rho_g: Quantity = Q_(0.2, 'frac'),
        sky_model: str = ''
    ) -> HourlyIrradiationComponents:
        """Returns the monthly averages of the hourly irradiation components on
        the tilted surface at the given hour and the arrays of the hourly
        irradiation components on the tilted surface at the same hour for every
        day of the given month.

        Parameters
        ----------
        month: int
            numerical index of the month of the year
        hr1: int
            start hour (solar time) of the 1-hour period
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
        HourlyIrradiationComponents, a named tuple with fields:
        - I: 2-tuple with the monthly average of total irradiation at `hr1` and
             a`Quantity`-array of the daily values of total irradiation at `hr1`.
        - I_b: same for beam irradiation
        - I_d: same for sky-diffuse irradiation
        - I_g: same for ground reflected irradiation
        """
        tmy_data = self.tmy.get_month(month)
        num_days = tmy_data.index.day.nunique()
        hic_lst = [
            self.tilted_hourly_irradiation(d, hr1, surf, rho_g, sky_model)
            for d in range(1, num_days + 1)
        ]
        I_T_arr = np.array([hic.I.to('J / m**2').m for hic in hic_lst])
        I_T_avg = np.mean(I_T_arr)
        I_Tb_arr = np.array([hic.I_b.to('J / m**2').m for hic in hic_lst])
        I_Tb_avg = np.mean(I_Tb_arr)
        I_Td_arr = np.array([hic.I_d.to('J / m**2').m for hic in hic_lst])
        I_Td_avg = np.mean(I_Td_arr)
        I_Tg_arr = np.array([hic.I_g.to('J / m**2').m for hic in hic_lst])
        I_Tg_avg = np.mean(I_Tg_arr)
        return HourlyIrradiationComponents(
            I=(Q_(I_T_avg, 'J / m ** 2'), Q_(I_T_arr, 'J / m ** 2')),
            I_b=(Q_(I_Tb_avg, 'J / m ** 2'), Q_(I_Tb_arr, 'J / m ** 2')),
            I_d=(Q_(I_Td_avg, 'J / m ** 2'), Q_(I_Td_arr, 'J / m ** 2')),
            I_g=(Q_(I_Tg_avg, 'J / m ** 2'), Q_(I_Tg_arr, 'J / m ** 2'))
        )

    def daily_irradiation(self, month: int, day: int) -> Quantity:
        """Returns the daily irradiation `H` on the horizontal surface for the
        given day in the given month.
        """
        self.location.date = Date(2001, month, day)
        tmy_data = self.tmy.get_day(month, day)
        G_data = tmy_data['G(h)']  # W / m ** 2
        omega_sr = self.location.omega_sr.to('rad').m
        omega_ss = self.location.omega_ss.to('rad').m
        omega = [
            hour_angle(time_to_decimal_hour(dt.time()))
            for dt in G_data.index
        ]
        G_interp = interp1d(omega, G_data.values, kind='slinear')
        H = quad(lambda omega: G_interp(omega), omega_sr, omega_ss, limit=100)[0]  # W * rad / m ** 2
        k = 1.0 / np.radians(15.0)  # hr / rad
        return Q_(k * H, 'W * hr / m ** 2').to('J / m ** 2')

    def avg_daily_irradiation(self, month: int) -> tuple[Quantity, Quantity]:
        """Returns the average daily irradiation `H_avg` on the horizontal
        surface for the given month and the array of daily irradiation `H` for
        every day in the given month.
        """
        tmy_data = self.tmy.get_month(month)
        num_days = tmy_data.index.day.nunique()
        H_arr = np.array([
            self.daily_irradiation(month, day).m
            for day in range(1, num_days + 1)
        ])
        H_avg = np.mean(H_arr)
        return Q_(H_avg, 'J / m ** 2'), Q_(H_arr, 'J / m ** 2')

    def avg_daily_clearness_index(self, month: int) -> tuple[float, np.ndarray]:
        """Returns the average daily clearness index `K_T_avg` on the horizontal
        surface for the given month and the array of daily clearness indexes
        `K_T` for every day in the given month.
        """
        tmy_data = self.tmy.get_month(month)
        num_days = tmy_data.index.day.nunique()
        H_o_arr = []
        for day in range(1, num_days + 1):
            self.location.date = Date(2001, month, day)
            H_o_arr.append(self.location.sun.H_o.to('J / m ** 2').m)
        H_o_avg = np.mean(H_o_arr)
        H_avg, H_arr = self.avg_daily_irradiation(month)
        K_T_avg = H_avg.to('J / m ** 2').m / H_o_avg
        K_T_arr = H_arr.to('J / m ** 2').m / H_o_arr
        return K_T_avg, K_T_arr

    def tilted_daily_irradiation(
        self,
        month: int,
        day: int,
        surf: Surface,
        rho_g: Quantity = Q_(0.2, 'frac'),
        sky_model: str = ''
    ) -> Quantity:
        """Returns the daily irradiation `H_T` on the tilted surface for the
        given day in the given month.
        """
        H = self.daily_irradiation(month, day)
        H_b, H_d = surf.location.sun.split_radiation(H=H)
        if sky_model.lower() == 'anisotropic.hdkr':
            SkyModel = transposition.AnisotropicSkyModel.HDKR
        elif sky_model.lower() == 'isotropic':
            SkyModel = transposition.IsotropicSkyModel
        else:
            SkyModel = transposition.AnisotropicSkyModel.Perez
        sky_model = SkyModel(
            surface=surf,
            rho_g=rho_g,
            G_b=H_b,
            G_d=H_d,
            time_interval=(surf.sunrise, surf.sunset)
        )
        H_T = sky_model.G_T.to('J / m ** 2')
        return H_T

    def tilted_avg_daily_irradiation(
        self,
        month: int,
        surf: Surface,
        rho_g: Quantity = Q_(0.2, 'frac'),
        sky_model: str = ''
    ) -> tuple[Quantity, Quantity]:
        """Returns the average of the daily irradiation `H_T_avg` on the tilted
        surface for the given month and the array of daily irradiation `H_T` on
        the tilted surface for every day of the given month.

        Parameters
        ----------
        month: int
            numerical index of the month of the year
        surf: Surface
            tilted surface.
        rho_g: Quantity, default 0.2 frac
            ground reflectance
        sky_model: str, {'anisotropic.hdkr, 'isotropic'}, optional
            sky model to be used for calculating the irradiation on the tilted
            surface. By default, the anisotropic sky model according to Perez
            is used.
        """
        tmy_data = self.tmy.get_month(month)
        num_days = tmy_data.index.day.nunique()
        H_T_arr = np.array([
            self.tilted_daily_irradiation(month, day, surf, rho_g, sky_model).m
            for day in range(1, num_days + 1)
        ])
        H_T_avg = np.mean(H_T_arr)
        return Q_(H_T_avg, 'J / m ** 2'), Q_(H_T_arr, 'J / m ** 2')

    def monthly_utilizable_energy(
        self,
        month: int,
        surf: Surface,
        I_Tc: list[Quantity] | Quantity,
        rho_g: Quantity = Q_(0.2, 'frac'),
        sky_model: str = ''
    ) -> Quantity:
        """
        Returns the average utilizable energy for the given month.

        Parameters
        ----------
        month: int
            numerical index of the month of the year
        surf: Surface
            tilted surface.
        I_Tc:
            The critical radiation level of the tilted surface at the given
            solar time. Either a single quantity, if the critical level is
            constant in time, or a list of 24 quantities that corresponds with
            each hour of solar time.
        rho_g: Quantity, default 0.2 frac
            ground reflectance
        sky_model: str, {'anisotropic.hdkr, 'isotropic'}, optional
            sky model to be used for calculating the irradiation on the tilted
            surface. By default, the anisotropic sky model according to Perez
            is used.
        """
        hr_range = range(0, 24)
        lI_avg, lI_T_avg = [], []
        for hr in hr_range:
            I_data, *_ = self.avg_hourly_irradiation(month, hr)
            I_T_data, *_ = self.tilted_avg_hourly_irradiation(month, hr, surf, rho_g, sky_model)
            I_avg = I_data[0]
            I_T_avg = I_T_data[0]
            lI_avg.append(I_avg)
            lI_T_avg.append(I_T_avg)
        H_avg, _ = self.avg_daily_irradiation(month)
        util = Utilizability(surf, month, H_avg)
        ue = util.monthly_utilizable_energy(
            t=[Time(h) for h in hr_range],
            I_T_avg=lI_T_avg,
            I_Tc=I_Tc,
            I_avg=lI_avg
        )
        return ue


def frequency_table(
    dataset: np.ndarray,
    bins: np.ndarray | str | int = 'auto'
) -> pd.DataFrame:
    """Returns a Pandas DataFrame with the frequency table of the given dataset.

    Parameters
    ----------
    dataset:
        Numpy array of values.
    bins:
        See Numpy documentation about `numpy.histogram`:
        If bins is an int, it defines the number of equal-width bins in the
        given range. If bins is a sequence, it defines a monotonically
        increasing array of bin edges, including the rightmost
        edge, allowing for non-uniform bin widths. If bins is a string, it
        defines the method used to calculate the optimal bin width, as defined
        by `histogram_bin_edges`. By default, bins is set to 'auto', which
        means that Numpy will determine the optimal bin width and consequently
        the number of bins.
    """
    freq, bin_edges = np.histogram(dataset, bins=bins)
    bin_edges = bin_edges.round(2)
    bins = [
        f"[{bin_edges[i]}, {bin_edges[i+1]})"
        for i in range(len(bin_edges) - 2)
    ]
    bins += [f"[{bin_edges[-2]}, {bin_edges[-1]}]"]
    cum_freq = np.cumsum(freq)
    df = pd.DataFrame({
        'bin': bins,
        'mid': [(bin_edges[i] + bin_edges[i+1]) / 2 for i in range(len(bin_edges) - 1)],
        'freq': freq,
        'cum_freq': cum_freq
    })
    return df
