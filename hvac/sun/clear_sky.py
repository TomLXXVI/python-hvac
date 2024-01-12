from collections import namedtuple
from datetime import date as Date
from datetime import time as Time
import calendar
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad
from hvac import Quantity
from hvac.charts.matplotlibwrapper import LineChart
from hvac.sun.surface import Location, Surface, hour_angle
from hvac.sun import transposition
from hvac.sun.utilizability import Utilizability
from hvac.sun.correlation import HourlyIrradiationComponents

Q_ = Quantity


IrradianceComponents = namedtuple(
    'IrradianceComponents',
    ['G', 'G_b', 'G_d', 'G_g']
)


class SolarRadiation:

    def __init__(self, location: Location):
        self.location = location

    def _get_hourly_irradiance_arrays(
        self,
        month: int,
        day: int
    ) -> IrradianceComponents:
        """Returns for the given day and month the hourly values of the
        irradiance components on a horizontal surface.

        Parameters
        ----------
        month:
            Number of the month in a year.
        day:
            Number of the day in the selected month.

        Returns
        -------
        IrradianceComponents, a named tuple with fields:
        - G: `Quantity` array of the global irradiance on a horizontal surface
             for each hour of the day (0..23).
        - G_b: Same, but with beam irradiance values instead.
        - G_d: Same, but with diffuse irradiance values instead.
        """
        self.location.date = Date(2001, month, day)
        G_arr, G_b_arr, G_d_arr = [], [], []
        for hr in range(24):
            self.location.solar_time = Time(hr)
            if self.location.omega_sr < self.location.omega < self.location.omega_ss:
                G_b = self.location.sun.G_cb.to('W / m**2').magnitude
                G_d = self.location.sun.G_cd.to('W / m**2').magnitude
                G = G_b + G_d
                G_arr.append(G)
                G_b_arr.append(G_b)
                G_d_arr.append(G_d)
            else:
                G_arr.append(0.0)
                G_b_arr.append(0.0)
                G_d_arr.append(0.0)
        return IrradianceComponents(
            G=Q_(G_arr, 'W / m**2'),
            G_b=Q_(G_b_arr, 'W / m**2'),
            G_d=Q_(G_d_arr, 'W / m**2'),
            G_g=None
        )

    def plot_irradiance(self, month: int, day: int) -> None:
        """Plot a curve for the given month and day of the total irradiance on
        the horizontal surface.
        """
        G_data, *_ = self._get_hourly_irradiance_arrays(month, day)
        G_data = [G.to('W / m**2').m for G in G_data]
        omega_data = [hour_angle(hr) for hr in range(24)]
        G_fun = interp1d(omega_data, G_data, kind='slinear')
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
        G_arr, G_b_arr, G_d_arr, _ = self._get_hourly_irradiance_arrays(month, day)
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
        G_data, G_b_data, G_d_data, _ = self._get_hourly_irradiance_arrays(month, day)
        G_data = [G.to('W / m**2').m for G in G_data]
        G_b_data = [G_b.to('W / m**2').m for G_b in G_b_data]
        G_d_data = [G_d.to('W / m**2').m for G_d in G_d_data]
        omega_data = [hour_angle(hr) for hr in range(24)]
        omega_sr = self.location.omega_sr.to('rad').m
        omega_ss = self.location.omega_ss.to('rad').m
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
        num_days = calendar.monthrange(2001, month)[1]
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
    ) -> HourlyIrradiationComponents:
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
        num_days = calendar.monthrange(2001, month)[1]
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
        G_data, *_ = self._get_hourly_irradiance_arrays(month, day)
        G_data = [G.to('W / m**2').m for G in G_data]
        omega = [hour_angle(hr) for hr in range(24)]
        omega_sr = self.location.omega_sr.to('rad').m
        omega_ss = self.location.omega_ss.to('rad').m
        G_interp = interp1d(omega, G_data, kind='slinear')
        H = quad(lambda omega: G_interp(omega), omega_sr, omega_ss, limit=100)[0]  # W * rad / m ** 2
        k = 1.0 / np.radians(15.0)  # hr / rad
        return Q_(k * H, 'W * hr / m ** 2').to('J / m ** 2')

    def avg_daily_irradiation(self, month: int) -> tuple[Quantity, Quantity]:
        """Returns the average daily irradiation `H_avg` on the horizontal
        surface for the given month and the array of daily irradiation `H` for
        every day in the given month.
        """
        num_days = calendar.monthrange(2001, month)[1]
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
        num_days = calendar.monthrange(2001, month)[1]
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
        num_days = calendar.monthrange(2001, month)[1]
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
            I_avg = I_data[0]
            I_T_data, *_ = self.tilted_avg_hourly_irradiation(month, hr, surf, rho_g, sky_model)
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
