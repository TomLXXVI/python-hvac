from datetime import time as Time
import numpy as np
import pandas as pd
from hvac import Quantity
from hvac.sun.surface import Surface
from hvac.sun.radiation import ReferenceDates, estimate_r_t
from hvac.sun.geometry import hour_angle
from hvac.sun.time import time_to_decimal_hour


class Utilizability:

    def __init__(
        self,
        surface: Surface,
        month: int,
        H_avg: Quantity
    ) -> None:
        """
        Creates a `Utilizability` instance to estimate the amount of utilizable
        energy of a tilted surface (solar collector).

        Parameters
        ----------
        surface:
            Tilted surface for which the utilizable energy is to be estimated.
        month:
            The month for which the utilizability of the tilted surface is to be
            determined.
        H_avg:
            The monthly average daily irradiation for the given month.
        """
        self.surface = surface
        self.ref_date = ReferenceDates.get_date_for(month)
        self.surface.location.date = self.ref_date
        self.num_days = pd.Period(str(self.ref_date), freq='D').days_in_month
        self._H_avg = H_avg.to('J / m ** 2').m
        H_o_avg = self.surface.location.sun.H_o.to('J / m ** 2').m
        self.K_T_avg = self._H_avg / H_o_avg

    def _get_k_T_avg(self, omega: float) -> float:
        omega_ss = self.surface.location.omega_ss.to('rad').m
        a = 0.409 + 0.5016 * np.sin(omega_ss - np.radians(60))
        b = 0.6609 - 0.4767 * np.sin(omega_ss - np.radians(60))
        k_T_avg = (a + b * np.cos(omega)) * self.K_T_avg
        return k_T_avg

    def utilizability(
        self,
        t: Time,
        I_T_avg: Quantity,
        I_Tc: Quantity,
        I_avg: Quantity | None = None
    ) -> float:
        """
        Get utilizability `fi` at solar time `t` according to the algorithm
        of Clark et al. (1983).

        Utilizability is the ratio of the amount of total received radiation
        that is above the critical level (i.e. the utilizable energy) to the
        total radiation incident on the tilted surface. The utilizability for a
        particular solar time in a certain month, is determined based on the
        monthly average of the hourly radiation on the tilted surface at that
        solar time.

        Source: "SOLAR ENGINEERING OF THERMAL PROCESSES, PHOTOVOLTAICS AND WIND"
        (5th. Ed.) by John A. Duffie and William A. Beckman, p. 124

        Parameters
        ----------
        t:
            Local solar time for which the utilizability is to be determined.
        I_T_avg:
            Monthly average of the hourly radiation on tilted surface at the
            given solar time.
        I_Tc:
            The critical radiation level of the tilted surface at the given
            solar time.
        I_avg: optional
            Monthly average of the hourly radiation on the horizontal surface
            at the given solar time. If None, `I_avg` will be estimated from
            the monthly average daily irradiation.
        """
        I_T_avg = I_T_avg.to('J / m ** 2').m
        I_Tc = I_Tc.to('J / m ** 2').m
        if I_T_avg <= 0.0:
            return 0.0
        else:
            beta = self.surface.beta.to('rad').m
            delta = self.surface.location.delta.to('rad').m
            omega_ss = self.surface.location.omega_ss.to('rad').m
            omega = hour_angle(time_to_decimal_hour(t))
            if I_avg is None:
                r_t = estimate_r_t(omega, omega_ss)
                I_avg = max(0.0, r_t * self._H_avg)
            else:
                I_avg = I_avg.to('J / m ** 2').m
            if I_avg <= 0.0:
                return 0.0
            X_c = I_Tc / I_T_avg
            R_h_avg = I_T_avg / I_avg
            k_T_avg = self._get_k_T_avg(omega)
            a = R_h_avg / k_T_avg ** 2
            b = np.cos(beta) / k_T_avg ** 2
            c = k_T_avg / np.cos(delta) ** 2
            X_m = 1.85 + 0.169 * a - 0.0696 * b - 0.981 * c
            g = (X_m - 1) / (2 - X_m)
            if X_c >= X_m:
                fi = 0.0
            elif X_m == 2.0:
                fi = (1 - X_c / X_m) ** 2
            else:
                a = abs(g)
                b = (g ** 2 + (1 + 2 * g) * (1 - X_c / X_m) ** 2) ** 0.5
                fi = abs(a - b)
        return fi

    def hourly_utilizable_energy(
        self,
        t: Time,
        I_T_avg: Quantity,
        I_Tc: Quantity,
        I_avg: Quantity | None = None
    ) -> Quantity:
        """
        Returns the average utilizable energy for the given month at the given
        solar time `t`.

        Parameters
        ----------
        t:
            Local solar time for which the utilizability is to be determined.
        I_T_avg:
            Monthly average of the hourly radiation on tilted surface at the
            given solar time.
        I_Tc:
            The critical radiation level of the tilted surface at the given
            solar time.
        I_avg: optional
            Monthly average of the hourly radiation on the horizontal surface
            at the given solar time. If None, `I_avg` will be estimated from
            the monthly average daily irradiation.
        """
        fi = self.utilizability(t, I_T_avg, I_Tc, I_avg)
        ue = self.num_days * fi * I_T_avg
        return ue.to('MJ / m ** 2')

    def monthly_utilizable_energy(
        self,
        t: list[Time],
        I_T_avg: list[Quantity],
        I_Tc: list[Quantity] | Quantity,
        I_avg: list[Quantity] | None = None
    ) -> Quantity:
        """
        Returns the average utilizable energy for the given month.

        Parameters
        ----------
        t:
            List of local solar times for which the utilizability is to be
            determined.
        I_T_avg:
            List of monthly average of the hourly radiation on tilted surface at
            the given range of solar times.
        I_Tc:
            The critical radiation level of the tilted surface at the given
            solar time. Either a single quantity, if the critical level is
            constant in time, or a list of quantities that corresponds with the
            list of local solar times.
        I_avg: optional
            List of monthly average of the hourly radiation on the horizontal
            surface at the given range of solar times. If None, `I_avg` will be
            estimated from the monthly average daily irradiation.
        """
        try:
            iter(I_Tc)
        except TypeError:
            I_Tc = [I_Tc] * len(I_T_avg)
        if I_avg is None:
            _zip = zip(t, I_T_avg, I_Tc)
        else:
            _zip = zip(t, I_T_avg, I_Tc, I_avg)
        ue = sum(self.hourly_utilizable_energy(*tup) for tup in _zip)
        return ue
