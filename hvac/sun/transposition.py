from datetime import time as Time
import numpy as np
from hvac import Quantity
from hvac.sun.surface import Surface
from hvac.sun.radiation import air_mass, estimate_avg_H_d_fraction, ReferenceDates


Q_ = Quantity


class IsotropicSkyModel:
    """
    Implementation of the isotropic sky model according to Liu and Jordan.

    Use this class to estimate the radiation (irradiance) on a tilted
    surface when the radiation on the horizontal surface is known and
    isotropic sky conditions are assumed.

    The total radiation on a tilted surface is the sum of beam radiation,
    isotropic diffuse radiation from the sky area, circumsolar diffuse
    radiation with the same direction as beam radiation, diffuse radiation
    from the horizon, and reflected radiation from buildings, fields, etc.
    ```
    G_T = (
        (G_b * R_b) + (G_d_iso * F_cs) + (G_d_cs * R_b)
        + (G_d_hz * F_chz) + (G * rho_g * F_cg)
    )
    ```
    In the isotropic diffuse model the radiation on a tilted surface is
    considered to include only 3 components: beam, isotropic diffuse, and
    solar radiation diffusely reflected from the ground. The third and
    fourth term (circumsolar diffuse and horizon brightening) are taken as
    zero.

    The isotropic model is the simplest, giving the most conservative
    estimates of radiation on the tilted surface, and has been widely used.

    This model can be applied for instantaneous radiation (irradiance `G`)
    or hourly irradiation (`I`).
    """

    def __init__(
        self,
        surface: 'Surface',
        rho_g: Quantity = Q_(0.2, 'frac'),
        G_b: Quantity | None = None,
        G_d: Quantity | None = None,
        time_interval: tuple[Time, Time] | None = None
    ) -> None:
        """
        Creates `IsotropicSkyModel` instance.

        Parameters
        ----------
        surface:
            the tilted surface for which the total radiation is to be
            determined
        rho_g: default 0.2 frac
            the reflectance (the albedo) of the ground
        G_b: optional
            beam irradiance or hourly irradiation on the horizontal surface.
        G_d: optional
            diffuse irradiance or hourly irradiation on the horizontal surface.
        time_interval: optional
            solar time interval (begin and end time) for which the total
            radiation is to be determined. Needed when hourly radiation
            `I_T` on the tilted surface is to be determined.
        """
        self.surface = surface
        self.rho_g = rho_g

        if G_b.check('[mass] / [time] ** 3'):  # irradiance
            self.G_b = G_b.to_base_units()
            self.G_d = G_d.to_base_units()
            self.R_b = self.surface.R_b

        if G_b.check('[mass] / [time] ** 2'):  # irradiation
            self.G_b = G_b.to_base_units()
            self.G_d = G_d.to_base_units()
            t1 = time_interval[0]
            t2 = time_interval[1]
            self.R_b = self.surface.R_b_ave(t1, t2)

        # view factor of the tilted surface to the sky:
        self.F_cs = (1 + np.cos(self.surface.beta.to('rad'))) / 2

        # view factor of the tilted surface to the ground:
        self.F_cg = (1 - np.cos(self.surface.beta.to('rad'))) / 2

    @property
    def G_T(self) -> Quantity:
        """Total irradiance or hourly irradiation on tilted surface."""
        return self.G_Tb + self.G_Td + self.G_Tg

    @property
    def G_Tb(self) -> Quantity:
        """Beam irradiance or hourly irradiation on tilted surface."""
        return self.G_b * self.R_b

    @property
    def G_Td(self) -> Quantity:
        """Diffuse sky irradiance or hourly irradiation on tilted surface."""
        return self.G_d * self.F_cs

    @property
    def G_Tg(self) -> Quantity:
        """Ground reflected irradiance or hourly irradiation on tilted surface."""
        G = self.G_b + self.G_d
        return G * self.rho_g * self.F_cg

    @property
    def R(self) -> Quantity:
        """Ratio of total irradiance or hourly irradiation on tilted surface to
        total radiation or hourly irradiation on horizontal surface.
        """
        G = self.G_b + self.G_d
        R = self.G_T / G
        return R.to('frac')


class AnisotropicSkyModel:
    """For estimating the irradiance `G` or hourly radiation `I` on a tilted
    surface when the irradiance/hourly irradiation on the horizontal surface is
    known and anisotropic sky conditions are assumed.

    The anisotropic sky model extends the isotropic sky model by also taking
    circumsolar diffuse and horizon brightening into account.
    """

    class HDKR:
        """Implementation of the anisotropic sky model according to Hay,
        Davies, Klucher and Reindl.

        Notes
        -----
        The HDKR model produces results that are closer to measured values than
        the isotropic sky model.
        The HDKR model is suggested for surfaces sloped toward the equator.
        """

        def __init__(
            self,
            surface: Surface,
            rho_g: Quantity = Q_(0.2, 'frac'),
            G_b: Quantity | None = None,
            G_d: Quantity | None = None,
            time_interval: tuple[Time, Time] | None = None
        ) -> None:
            """Creates `AnisotropicSkyModel.HDKR` instance.

            Parameters
            ----------
            surface:
                the tilted surface for which the total radiation is to be
                determined
            rho_g: default 0.2 frac
                the reflectance (the albedo) of the ground
            G_b: optional
                beam irradiance or hourly irradiation on the horizontal surface.
            G_d: optional
                diffuse irradiance or hourly irradiation on the horizontal surface.
            time_interval: optional
                solar time interval (begin and end time) for which the total
                radiation is to be determined. Use this when irradiation is to
                be determined instead of irradiance.
            """
            self.surface = surface
            self.rho_g = rho_g

            if G_b.check('[mass] / [time] ** 3'):  # irradiance
                self.G_b = G_b.to_base_units()
                self.G_d = G_d.to_base_units()
                self.R_b = self.surface.R_b
                self.A_i = self.G_b / surface.location.sun.G_o  # anisotropy index

            if G_b.check('[mass] / [time] ** 2'):  # irradiation
                self.G_b = G_b.to_base_units()
                self.G_d = G_d.to_base_units()
                t1 = time_interval[0]
                t2 = time_interval[1]
                self.R_b = self.surface.R_b_ave(t1, t2)
                self.A_i = self.G_b / surface.location.sun.I_o(t1, t2)

            # modulating factor to account for cloudiness
            self.f = np.sqrt(self.G_b / (self.G_b + self.G_d))

            # view factor of the tilted surface to the sky:
            self.F_cs = (1 + np.cos(self.surface.beta.to('rad'))) / 2

            # view factor of the tilted surface to the ground:
            self.F_cg = (1 - np.cos(self.surface.beta.to('rad'))) / 2

        @property
        def G_T(self) -> Quantity:
            """Total irradiance (hourly irradiation) on tilted surface."""
            return self.G_Tb + self.G_Td + self.G_Tg

        @property
        def G_Tb(self) -> Quantity:
            """Beam irradiance (hourly irradiation) on tilted surface."""
            return (self.G_b + self.G_d * self.A_i) * self.R_b

        @property
        def G_Td(self) -> Quantity:
            """Diffuse sky irradiance (hourly irradiation) on tilted
            surface.
            """
            # correction factor for horizon brightening:
            f_hz = 1.0 + self.f * np.sin(self.surface.beta.to('rad') / 2) ** 3
            return self.G_d * (1.0 - self.A_i) * self.F_cs * f_hz

        @property
        def G_Tg(self) -> Quantity:
            """Ground reflected irradiance (hourly irradiation) on tilted
            surface."""
            return (self.G_b + self.G_d) * self.rho_g * self.F_cg

        @property
        def R(self) -> Quantity:
            """Ratio of total irradiance (hourly irradiation) on tilted
            surface to total irradiance (hourly irradiation) on horizontal
            surface.
            """
            G = self.G_b + self.G_d
            R = self.G_T / G
            return R.to('frac')

    class Perez:
        """Implementation of Perez' model for anisotropic sky conditions.

        Notes
        -----
        The Perez model is the least conservative model. It agrees the best
        by a small margin with measurements. For surfaces with azimuth angle
        far from 0° in the northern hemisphere or 180° in the southern
        hemisphere, the Perez model is suggested.
        """
        _f_dict = {
            (1.000, 1.065): {
                'f11': -0.008, 'f12': 0.588, 'f13': -0.062,
                'f21': -0.060, 'f22': 0.072, 'f23': -0.022
            },
            (1.065, 1.230): {
                'f11': 0.130, 'f12': 0.683, 'f13': -0.151,
                'f21': -0.019, 'f22': 0.066, 'f23': -0.029
            },
            (1.230, 1.500): {
                'f11': 0.330, 'f12': 0.487, 'f13': -0.221,
                'f21': 0.055, 'f22': -0.064, 'f23': -0.026
            },
            (1.500, 1.950): {
                'f11': 0.568, 'f12': 0.187, 'f13': -0.295,
                'f21': 0.109, 'f22': -0.152, 'f23': -0.014
            },
            (1.950, 2.800): {
                'f11': 0.873, 'f12': -0.392, 'f13': -0.362,
                'f21': 0.226, 'f22': -0.462, 'f23': 0.001
            },
            (2.800, 4.500): {
                'f11': 1.132, 'f12': -1.237, 'f13': -0.412,
                'f21': 0.288, 'f22': -0.823, 'f23': 0.056
            },
            (4.500, 6.200): {
                'f11': 1.060, 'f12': -1.600, 'f13': -0.359,
                'f21': 0.264, 'f22': -1.127, 'f23': 0.131
            },
            (6.200, float('inf')): {
                'f11': 0.678, 'f12': -0.327, 'f13': -0.250,
                'f21': 0.156, 'f22': -1.377, 'f23': 0.251
            }
        }
        # statistically derived coefficients for ranges of values of
        # clearness parameter `eta` that determine the brightness
        # coefficients `F1` and `F2`

        def __init__(
            self,
            surface: 'Surface',
            rho_g: Quantity = Q_(0.2, 'frac'),
            G_b: Quantity | None = None,
            G_d: Quantity | None = None,
            time_interval: tuple[Time, Time] | None = None
        ) -> None:
            """Creates `AnisotropicSkyModel.Perez` instance.

            Parameters
            ----------
            surface:
                the tilted surface for which the total radiation is to be
                determined
            rho_g: default 0.2 frac
                the reflectance (the albedo) of the ground
            G_b: optional
                beam radiation (irradiance) on the horizontal surface.
                If None, clear sky beam irradiance will be assumed.
            G_d: optional
                diffuse radiation (irradiance) on the horizontal surface.
                If None, clear sky diffuse irradiance will be assumed.
            time_interval: optional
                solar time interval (begin and end time) for which the total
                radiation is to be determined. Use this when irradiation is to
                be determined instead of irradiance.
            """
            self.surface = surface
            self.rho_g = rho_g

            if G_b.check('[mass] / [time] ** 3'):  # irradiance
                self.G_b = G_b.to_base_units()
                self.G_d = G_d.to_base_units()
                self.R_b = self.surface.R_b

            if G_b.check('[mass] / [time] ** 2'):  # irradiation
                self.G_b = G_b.to_base_units()
                self.G_d = G_d.to_base_units()
                t1 = time_interval[0]
                t2 = time_interval[1]
                self.R_b = self.surface.R_b_ave(t1, t2)

            # view factor of the tilted surface to the sky:
            self.F_cs = (1 + np.cos(self.surface.beta.to('rad'))) / 2

            # view factor of the tilted surface to the ground:
            self.F_cg = (1 - np.cos(self.surface.beta.to('rad'))) / 2

        @property
        def _eta(self) -> float:
            """Clearness parameter."""
            theta_z = self.surface.location.sun.theta.to('deg')
            G_bn = self.G_b / np.cos(theta_z.to('rad'))
            eta_n = (self.G_d + G_bn) / self.G_d + 5.535e-6 * theta_z ** 3
            eta_d = 1 + 5.535e-6 * theta_z ** 3
            eta = eta_n / eta_d
            return eta.to('frac').m

        @property
        def _Delta(self) -> float:
            """Brightness parameter."""
            m = air_mass(self.surface.location.sun.theta.to('rad'))
            if not self.G_d.check('[mass] / [time] ** 2'):
                # if G_d is a power quantity
                G_on = self.surface.location.sun.G_on
                Delta = m * self.G_d / G_on
            else:
                # if G_d is an energy quantity
                I_on = self.surface.location.sun.G_on * Q_(1, 'hr')
                Delta = m * self.G_d / I_on
            return Delta.to('frac').m

        @property
        def _F1(self) -> float:
            f11, f12, f13 = 0.0, 0.0, 0.0
            for rng in self._f_dict.keys():
                if rng[0] <= self._eta < rng[1]:
                    f11 = self._f_dict[rng]['f11']
                    f12 = self._f_dict[rng]['f12']
                    f13 = self._f_dict[rng]['f13']
                    break
            theta_z = self.surface.location.sun.theta.to('rad').m
            F1 = max(0, f11 + f12 * self._Delta + f13 * theta_z)
            return F1

        @property
        def _F2(self) -> float:
            f21, f22, f23 = 0.0, 0.0, 0.0
            for rng in self._f_dict.keys():
                if rng[0] <= self._eta < rng[1]:
                    f21 = self._f_dict[rng]['f21']
                    f22 = self._f_dict[rng]['f22']
                    f23 = self._f_dict[rng]['f23']
                    break
            theta_z = self.surface.location.sun.theta.to('rad').m
            F2 = f21 + f22 * self._Delta + f23 * theta_z
            return F2

        @property
        def G_T(self) -> Quantity:
            """Total irradiance (hourly irradiation) on tilted surface."""
            return self.G_Tb + self.G_Td + self.G_Tg

        @property
        def G_Tb(self) -> Quantity:
            """Beam irradiance (hourly irradiation) on tilted surface."""
            return self.G_b * self.R_b

        @property
        def G_Td(self) -> Quantity:
            """Diffuse sky irradiance (hourly irradiation) on tilted
            surface.
            """
            a = max(0, np.cos(self.surface.theta.to('rad').m))
            b = max(
                np.cos(np.radians(85.0)),
                np.cos(self.surface.location.sun.theta.to('rad').m)
            )
            G_Td = self.G_d * (
                    (1 - self._F1) * self.F_cs
                    + self._F1 * a / b
                    + self._F2 * np.sin(self.surface.beta.to('rad').m)
            )
            return G_Td

        @property
        def G_Tg(self) -> Quantity:
            """Ground reflected irradiance (hourly irradiation) on tilted
            surface.
            """
            return (self.G_b + self.G_d) * self.rho_g * self.F_cg

        @property
        def R(self) -> Quantity:
            """Ratio of total irradiance (hourly irradiation) on tilted
            surface to total irradiance (hourly irradiation) on horizontal
            surface.
            """
            G = self.G_b + self.G_d
            R = self.G_T / G
            return R.to('frac')


class KTMethod:
    """Implementation of the method of Klein and Theilacker for estimating the
    monthly average daily radiation `H_T_avg` on a tilted surface, when the
    monthly average daily radiation `H_avg` on the horizontal surface is known.
    """
    def __init__(
        self,
        month: int,
        surface: Surface,
        H_avg: Quantity | None,
        rho_g: Quantity = Q_(0.2, 'frac')
    ) -> None:
        """
        Creates `KTMethod` instance.

        Parameters
        ----------
        month:
            index of the month of the year.
        surface:
            the tilted surface for which the total radiation is to be
            determined
        H_avg:
            monthly average daily radiation on the horizontal surface
        rho_g: default 0.2 frac
            the reflectance (the albedo) of the ground

        Notes
        -----
        When calculating the monthly average daily radiation on the tilted
        surface, make sure that the date at the location of the surface
        has been set to the reference date of the month under consideration
        (see class `ReferenceDates` in module `radiation.py` to get the
        reference date of a given month).
        """
        self.surface = surface
        self.rho_g = rho_g.to('frac').m
        self.H_avg = H_avg.to_base_units()

        self.surface.location.date = ReferenceDates.get_date_for(month)
        H_o_avg = self.surface.location.sun.H_o.to_base_units()
        K_T_avg = self.H_avg / H_o_avg
        self._frac_H_d_avg = estimate_avg_H_d_fraction(
            K_T_avg.to('frac').m,
            omega_ss=self.surface.location.omega_ss.to('deg').m
        )

        # view factor of the tilted surface to the sky:
        self.F_cs = (1 + np.cos(self.surface.beta.to('rad').m)) / 2

        # view factor of the tilted surface to the ground:
        self.F_cg = (1 - np.cos(self.surface.beta.to('rad').m)) / 2

        beta = self.surface.beta.to('rad').m
        fi = self.surface.location.fi.to('rad').m
        gamma = self.surface.gamma.to('rad').m
        delta = self.surface.location.delta.to('rad').m
        omega_s = self.surface.location.omega_ss.to('rad').m
        self._a = (
            0.409
            + 0.5016 * np.sin(omega_s - np.radians(60))
        )
        self._b = (
            0.6609
            - 0.4767 * np.sin(omega_s - np.radians(60))
        )
        self._d = (
            np.sin(omega_s)
            - omega_s * np.cos(omega_s)
        )
        self._a_apo = self._a - self._frac_H_d_avg
        self._A = (
            np.cos(beta)
            + np.tan(fi) * np.cos(gamma) * np.sin(beta)
        )
        self._B = (
            np.cos(omega_s) * np.cos(beta)
            + np.tan(delta) * np.sin(beta) * np.cos(gamma)
        )
        self._C = (
            np.sin(beta) * np.sin(gamma)
            / np.cos(fi)
        )
        omega_sr = self.surface.omega_sr.to('rad').m
        omega_ss = self.surface.omega_ss.to('rad').m
        if omega_ss >= omega_sr:
            self._D = max(0.0, self._G(omega_ss, omega_sr))
        else:
            omega_s = self.surface.location.omega_ss.to('rad').m
            self._D = max(
                0.0,
                self._G(omega_ss, -omega_sr) + self._G(omega_s, -omega_sr)
            )

    def _omega_sr(self) -> float:
        """Sunrise hour angle of sloped surface."""
        omega_s = self.surface.location.omega_ss.to('rad').m
        a_n = (
            self._A * self._B
            + self._C
            * np.sqrt(self._A ** 2 - self._B ** 2 + self._C ** 2)
        )
        a_d = self._A ** 2 + self._C ** 2
        a = a_n / a_d
        omega_sr = min(omega_s, np.arccos(a))
        if (self._A > 0 and self._B > 0) or (self._A >= self._B):
            omega_sr = -omega_sr
        return omega_sr

    def _omega_ss(self) -> float:
        """Sunset hour angle of sloped surface."""
        omega_s = self.surface.location.omega_ss.to('rad').m
        a_n = (
            self._A * self._B
            - self._C
            * np.sqrt(self._A ** 2 - self._B ** 2 + self._C ** 2)
        )
        a_d = self._A ** 2 + self._C ** 2
        a = a_n / a_d
        omega_ss = -min(omega_s, np.arccos(a))
        if (self._A > 0 and self._B > 0) or (self._A >= self._B):
            omega_ss = -omega_ss
        return omega_ss

    def _G(self, omega1: float, omega2: float) -> float:
        G = (1 / (2 * self._d) * (
            (self._b * self._A / 2 - self._a_apo * self._B)
            * (omega1 - omega2)
            + (self._a_apo * self._A - self._b * self._B)
            * (np.sin(omega1) - np.sin(omega2))
            - self._a_apo * self._C
            * (np.cos(omega1) - np.cos(omega2))
            + (self._b * self._A / 2)
            * (
                np.sin(omega1) * np.cos(omega1)
                - np.sin(omega2) * np.cos(omega2)
            )
            + (self._b * self._C / 2)
            * (np.sin(omega1) ** 2 - np.sin(omega2) ** 2)
        ))
        return G

    @property
    def R_avg(self) -> Quantity:
        """Ratio of monthly average daily total radiation on the tilted surface
        to that on a horizontal surface.
        """
        R_avg = (
            self._D
            + self._frac_H_d_avg * self.F_cs
            + self.rho_g * self.F_cg
        )
        return Q_(R_avg, 'frac')

    @property
    def H_T_avg(self) -> Quantity:
        """Monthly average daily total radiation on tilted surface."""
        return self.R_avg * self.H_avg
