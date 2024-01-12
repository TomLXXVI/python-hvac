from collections import namedtuple
from datetime import date as Date
from datetime import time as Time
import numpy as np
import pandas as pd
from hvac import Quantity
from hvac.charts.matplotlibwrapper import BarChart
from hvac.sun.surface import Location, Surface
from hvac.sun.geometry import hour_angle
from hvac.sun.radiation import (
    ReferenceDates,
    cumulative_K_T_histogram,
    estimate_H_d_fraction,
    estimate_r_t,
    # estimate_r_d,
    estimate_I_d_fraction
)
from hvac.sun import transposition
from hvac.sun.utilizability import Utilizability

Q_ = Quantity


DailyIrradiationDistribution = namedtuple(
    'DailyIrradiationDistribution',
    ['H', 'H_b', 'H_d']
)
HourlyIrradiationDistribution = namedtuple(
    'HourlyIrradiationDistribution',
    ['I', 'I_b', 'I_d']
)
HourlyIrradiationComponents = namedtuple(
    'HourlyIrradiationComponents',
    ['I', 'I_b', 'I_d', 'I_g']
)


class SolarRadiation:
    """Estimation of solar radiation on horizontal and tilted surfaces at a
    given location during a given month if only the monthly average daily
    clearness index is known, using correlations from "SOLAR ENGINEERING OF
    THERMAL PROCESSES, PHOTOVOLTAICS AND WIND" (5th. Ed.), Chapters 1 and 2.
    """
    def __init__(self, location: Location, month: int, K_T_avg: float):
        """
        Creates a `SolarRadiation` instance.

        Parameters
        ----------
        location:
            Geographical location under consideration.
        month:
            Month in which solar radiation is to be estimated.
        K_T_avg:
            Monthly average daily clearness index for the given month.
        """
        self.location = location
        self.month = month
        self.ref_date = ReferenceDates.get_date_for(self.month)
        self.num_days = pd.Period(self.ref_date, freq='D').days_in_month
        self.K_T_avg = K_T_avg
        self.K_T_arr = self._get_K_T_distribution()
        self.H_arr, self.H_d_arr, self.H_b_arr = self._get_H_distribution(self.K_T_arr)
        self.H_avg = np.mean(self.H_arr)

    def _get_K_T_distribution(self) -> np.ndarray:
        """Get a distribution of the daily clearness index for each day of the
        month. Returns a Numpy array of which the number of elements equals the
        number of days in the month and the values of the elements are the
        daily clearness indexes.
        See "SOLAR ENGINEERING OF THERMAL PROCESSES, PHOTOVOLTAICS AND
        WIND" (5th. Ed.), p. 73, §2.9.
        """
        def _get_dist(bin_edges, hist) -> np.ndarray:
            """Create a distribution from a given histogram."""
            bin_data = [
                np.random.uniform(
                    low=bin_edges[i],
                    high=bin_edges[i+1],
                    size=int(hist[i])
                ) for i in range(len(hist))
            ]
            dist = np.concatenate(bin_data)
            np.random.shuffle(dist)
            return dist

        bin_width = 0.01
        # create a range of K_T bins
        K_T_edges = np.arange(0, 1 + bin_width, bin_width)
        # get a cumulative frequency distribution of days expressed as a
        # fraction of total number of days in the month
        cum_frac_arr = cumulative_K_T_histogram(K_T_edges, self.K_T_avg)
        # remove impossible values (negative values or values greater than 1)
        K_T_edges = K_T_edges[(cum_frac_arr > 0) & (cum_frac_arr <= 1.0)]
        cum_frac_arr = cum_frac_arr[(cum_frac_arr > 0) & (cum_frac_arr <= 1.0)]
        # convert to a cumulative frequency distribution of days in the month
        cum_days_arr = np.round(self.num_days * cum_frac_arr)
        # demount the cumulative frequency distribution into a frequency
        # distribution of days in the month
        days_arr = [cum_days_arr[0]]
        days_arr += [
            cum_days_arr[i] - cum_days_arr[i-1]
            for i in range(1, len(cum_days_arr))
        ]
        days_arr = np.array(days_arr)
        # add lower limit to K_T range
        K_T_edges = np.insert(K_T_edges, 0, K_T_edges[0] - bin_width)
        # generate a random distribution from the given histogram
        K_T_arr = _get_dist(K_T_edges, days_arr)
        return K_T_arr

    def _get_H_distribution(
        self,
        K_T_arr: np.ndarray
    ) -> DailyIrradiationDistribution:
        """Get the daily irradiation distribution for the whole month.
        Returns a tuple with the distribution of daily total irradiation,
        of daily diffuse irradiation, and of daily beam irradiation.
        A 'distribution' is a `Quantity` Numpy array of which the number of
        elements equals the number of days in the month and the values of the
        elements are the daily irradiation values.
        """
        num_days = K_T_arr.size
        dates = [Date(2001, self.month, d) for d in range(1, num_days + 1)]
        H_arr, H_d_arr = [], []
        for date, K_T in zip(dates, K_T_arr):
            self.location.date = date
            H_o = self.location.sun.H_o.to('J / m ** 2').m
            H = K_T * H_o
            H_arr.append(H)
            omega_ss = self.location.omega_ss.to('rad').m
            H_d_fr = estimate_H_d_fraction(K_T, omega_ss)
            H_d = H_d_fr * H
            H_d_arr.append(H_d)
        H_arr = Q_(np.array(H_arr), 'J / m ** 2')
        H_d_arr = Q_(np.array(H_d_arr), 'J / m ** 2')
        H_b_arr = H_arr - H_d_arr
        return DailyIrradiationDistribution(
            H=H_arr,
            H_b=H_b_arr,
            H_d=H_d_arr
        )

    def _get_I_distribution(self, day: int) -> HourlyIrradiationDistribution:
        """Returns the hourly irradiation distribution of a single day."""

        def _is_daylight() -> list[bool]:
            omega_sr = self.location.omega_sr
            omega_ss = self.location.omega_ss
            result = []
            for hr1 in range(0, 24):
                omega1 = Q_(hour_angle(hr1), 'rad')
                omega2 = Q_(hour_angle(hr1 + 1), 'rad')
                if omega1 <= omega_sr:
                    result.append(False)
                elif omega2 >= omega_ss:
                    result.append(False)
                else:
                    result.append(True)
            return result

        H = self.H_arr[day - 1]
        # H_d = self.H_d_arr[day - 1]
        date = Date(2001, self.month, day)
        self.location.date = date
        omega_ss = self.location.omega_ss.to('rad').m
        hr_rng = np.arange(0.5, 24.5, 1.0)
        omega_rng = hour_angle(hr_rng)
        r_t_rng = estimate_r_t(omega_rng, omega_ss)
        I_arr = r_t_rng * H
        t_rng = [
            (Time(hr), Time(hr + 1))
            for hr in range(0, 23)
        ]
        t_rng += [(Time(23), Time(0))]
        I_o_arr = [
            self.location.sun.I_o(t1, t2).to('J / m ** 2').m
            for t1, t2 in t_rng
        ]
        I_o_arr = Q_(np.array(I_o_arr), 'J / m ** 2')
        k_T_arr = [
            (I / I_o).m if I_o != 0 else 0
            for I, I_o in zip(I_arr, I_o_arr)
        ]
        I_d_fr_arr = np.array([estimate_I_d_fraction(k_T) for k_T in k_T_arr])
        I_d_arr = I_d_fr_arr * I_arr
        # r_d_rng = estimate_r_d(omega_rng, omega_ss)
        # I_d_arr = r_d_rng * H_d
        I_arr = np.where(_is_daylight(), I_arr, 0.0)
        I_d_arr = np.where(_is_daylight(), I_d_arr, 0.0)
        I_b_arr = I_arr - I_d_arr
        return HourlyIrradiationDistribution(
            I=I_arr,
            I_b=I_b_arr,
            I_d=I_d_arr
        )

    def get_I_distribution(
        self,
        day: int | None = None
    ) -> list[HourlyIrradiationDistribution] | HourlyIrradiationDistribution:
        """Returns the hourly irradiation distribution of a single day or,
        if `day` is `None`, returns a list with the hourly irradiation
        distribution of all days in the month.

        The hourly irradiation distribution of a single day is packed in a
        tuple. This tuple has 3 elements:
        0.  the hourly total irradiation distribution
        1.  the hourly diffuse irradiation distribution
        2.  the hourly beam irradiation distribution

        A 'distribution' is a `Quantity` Numpy array of which the number of
        elements equals the number of hours in the day and the values of the
        elements are the hourly irradiation values.
        """
        if day is None:
            ref_date = ReferenceDates.get_date_for(self.month)
            num_days = pd.Period(ref_date, freq='D').days_in_month
            irr_data_arr = [
                self._get_I_distribution(day)
                for day in range(1, num_days + 1)
            ]
            return irr_data_arr
        else:
            return self._get_I_distribution(day)

    def plot_hourly_irradiation(self, day: int) -> None:
        I_arr, I_b_arr, I_d_arr = self._get_I_distribution(day)
        hr_arr = [hr for hr in range(24)]
        bar_width = (hr_arr[1] - hr_arr[0]) / 2 * 0.8
        chart = BarChart()
        chart.add_xy_data(
            label='total irradiation',
            x1_values=hr_arr,
            y1_values=I_arr.to('kJ / m ** 2').m.tolist(),
            style_props={
                'width': bar_width,
                'align': 'edge'
            }
        )
        chart.add_xy_data(
            label='beam irradiation',
            x1_values=[hr + bar_width for hr in hr_arr],
            y1_values=I_b_arr.to('kJ / m ** 2').m.tolist(),
            style_props={
                'width': bar_width,
                'align': 'edge'
            }
        )
        chart.add_xy_data(
            label='diffuse irradiation',
            x1_values=[hr + bar_width for hr in hr_arr],
            y1_values=I_d_arr.to('kJ / m ** 2').m.tolist(),
            style_props={
                'width': bar_width,
                'align': 'edge',
                'bottom': I_b_arr.to('kJ / m ** 2').m.tolist()
            }
        )
        chart.x1.add_title('solar hour')
        chart.x1.scale(0, 24, 1)
        chart.y1.add_title('irradiation, kJ/m²')
        chart.add_legend(columns=3)

        date = Date(2001, self.month, day).strftime('%b %d')
        text = f"{date}\n" + "$K_{T}$ = " + f"{round(self.K_T_arr[day - 1], 2)}"
        chart.y1.axes.text(
            0.05, 0.95, text,
            transform=chart.y1.axes.transAxes,
            fontsize=12,
            verticalalignment='top',
            bbox={'boxstyle': 'round', 'facecolor': 'wheat', 'alpha': 0.5}
        )
        chart.show()

    def hourly_irradiation(self, day: int, hr1: int) -> HourlyIrradiationComponents:
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
        I_arr, I_b_arr, I_d_arr = self._get_I_distribution(day)
        return HourlyIrradiationComponents(
            I=I_arr[hr1],
            I_b=I_b_arr[hr1],
            I_d=I_d_arr[hr1],
            I_g=None
        )

    def avg_hourly_irradiation(self, hr1: int) -> HourlyIrradiationComponents:
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
        hic_lst = [
            self.hourly_irradiation(d, hr1) for d in range(1, self.num_days + 1)
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
        I, I_b, I_d, _ = self.hourly_irradiation(day, hr1)
        if I > 0:
            hr2 = hr1 + 1
            t1, t2 = Time(hr1), Time(hr2)
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
            I_T = sky_model.G_T.to('J / m**2')
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
        hr1: int,
        surf: Surface,
        rho_g: Quantity = Q_(0.2, 'frac'),
        sky_model: str = ''
    ) -> tuple[Quantity, Quantity]:
        """Returns the monthly averages of the hourly irradiation components on
        the tilted surface at the given hour and the arrays of the hourly
        irradiation components on the tilted surface at the same hour for every
        day of the given month.

        Parameters
        ----------
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
        hic_lst = [
            self.tilted_hourly_irradiation(d, hr1, surf, rho_g, sky_model)
            for d in range(1, self.num_days + 1)
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

    def daily_irradiation(self, day: int) -> Quantity:
        """Returns the daily irradiation `H` on the horizontal surface for the
        given day in the given month.
        """
        return self.H_arr[day - 1]

    def avg_daily_irradiation(self) -> tuple[Quantity, Quantity]:
        """Returns the average daily irradiation `H_avg` on the horizontal
        surface for the given month and the array of daily irradiation `H` for
        every day in the given month.
        """
        return self.H_avg, self.H_arr

    def avg_daily_clearness_index(self) -> tuple[float, np.ndarray]:
        """Returns the average daily clearness index `K_T_avg` on the horizontal
        surface for the given month and the array of daily clearness indexes
        `K_T` for every day in the given month.
        """
        return self.K_T_avg, self.K_T_arr

    def tilted_daily_irradiation(
        self,
        day: int,
        surf: Surface,
        rho_g: Quantity = Q_(0.2, 'frac'),
        sky_model: str = ''
    ) -> Quantity:
        """Returns the daily irradiation `H_T` on the tilted surface for the
        given day in the given month.
        """
        H_b = self.H_b_arr[day - 1]
        H_d = self.H_d_arr[day - 1]
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
        surf: Surface,
        rho_g: Quantity = Q_(0.2, 'frac'),
        sky_model: str = ''
    ) -> tuple[Quantity, Quantity]:
        """Returns the average of the daily irradiation `H_T_avg` on the tilted
        surface for the given month and the array of daily irradiation `H_T` on
        the tilted surface for every day of the given month.

        Parameters
        ----------
        surf: Surface
            tilted surface.
        rho_g: Quantity, default 0.2 frac
            ground reflectance
        sky_model: str, {'anisotropic.hdkr, 'isotropic'}, optional
            sky model to be used for calculating the irradiation on the tilted
            surface. By default, the anisotropic sky model according to Perez
            is used.
        """
        H_T_arr = np.array([
            self.tilted_daily_irradiation(d, surf, rho_g, sky_model).m
            for d in range(1, self.num_days + 1)
        ])
        H_T_avg = np.mean(H_T_arr)
        return Q_(H_T_avg, 'J / m ** 2'), Q_(H_T_arr, 'J / m ** 2')

    def monthly_utilizable_energy(
        self,
        surf: Surface,
        I_Tc: list[Quantity] | Quantity,
        rho_g: Quantity = Q_(0.2, 'frac'),
        sky_model: str = ''
    ) -> Quantity:
        """
        Returns the average utilizable energy for the given month.

        Parameters
        ----------
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
            I_data, *_ = self.avg_hourly_irradiation(hr)
            I_avg = I_data[0]
            I_T_data, *_ = self.tilted_avg_hourly_irradiation(hr, surf, rho_g, sky_model)
            I_T_avg = I_T_data[0]
            lI_avg.append(I_avg)
            lI_T_avg.append(I_T_avg)
        util = Utilizability(surf, self.month, self.H_avg)
        ue = util.monthly_utilizable_energy(
            t=[Time(h) for h in hr_range],
            I_T_avg=lI_T_avg,
            I_Tc=I_Tc,
            I_avg=lI_avg
        )
        return ue
