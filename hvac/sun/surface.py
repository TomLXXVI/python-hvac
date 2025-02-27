# Source:
# Duffie, J. A., Beckman, W. A., & Blair, N. (2020).
# SOLAR ENGINEERING OF THERMAL PROCESSES, PHOTOVOLTAICS AND WIND.
# John Wiley & Sons.

from datetime import date as Date
from datetime import time as Time
import numpy as np
import pandas as pd
from scipy.optimize import root_scalar
from scipy.integrate import quad
from .. import Quantity
from .time import (
    day_number,
    time_to_decimal_hour,
    hour_angle_to_solar_time,
    decimal_hour_to_time
)
from .geometry import (
    sunset_hour_angle,
    declination,
    hour_angle,
    zenith_angle,
    solar_azimuth_angle,
    solar_altitude_angle,
    incidence_angle,
    profile_angle,
    daylight_time
)
from .radiation import (
    ClimateType,
    ext_irradiance_normal,
    ext_irradiance_hor_surf,
    clear_sky_beam_transmittance,
    clear_sky_diffuse_transmittance,
    estimate_I_d_fraction,
    estimate_H_d_fraction,
    estimate_avg_H_d_fraction,
    estimate_r_t,
    estimate_r_d
)

Q_ = Quantity


class Location:

    def __init__(
        self,
        fi: Quantity,
        L_loc: Quantity | None = None,
        date: Date | None = None,
        solar_time: Time | None = None,
        altitude: Quantity = Q_(0, 'm'),
        climate_type: str = ClimateType.MID_LATITUDE_SUMMER,
        timezone: str = 'UTC'
    ) -> None:
        """
        Creates a `Location` object with which the solar path and the solar
        radiation can be determined at a given date and solar time.

        Parameters
        ----------
        fi:
            Latitude of the location; north positive; -pi/2 <= fi <= pi/2
        L_loc:
            Longitude of the location; positive east of Greenwich, negative
            west of Greenwich.
        date: Date, optional
            Current date at the location.
        solar_time: Time, optional
            Current solar time at the location.
        altitude: Quantity, optional
            Altitude of the location.
        climate_type:
            Type of climate as defined in class `ClimateType`
        timezone: optional, default UTC
            Tz-database identifier of the time zone of the location, indicating
            the offset from UTC (e.g. Etc/GMT-1).
            See: https://en.wikipedia.org/wiki/List_of_tz_database_time_zones
            The time zone is only used to determine the local standard time
            meridian. Do not specify a timezone with daylight saving time (DST),
            should it be in use at the location under consideration, as this may
            give errors at times when DST becomes active or inactive (e.g. use
            'Etc/GMT-1' instead of 'Europe/Brussels').
        """
        self._fi = fi.to('rad')  # `.to(...)` returns new object
        self._L_loc = L_loc
        self._date = date
        self._solar_time = Time(12, 0, 0) if solar_time is None else solar_time
        self._altitude = altitude
        self._climate_type = climate_type
        self._timezone = timezone
        self.sun = Sun(self)

    @property
    def date(self) -> Date:
        """Current date at the location."""
        return self._date

    @date.setter
    def date(self, date: Date):
        """Set the current date at the location."""
        self._date = date

    @property
    def solar_time(self) -> Time:
        """Current solar time at the location."""
        return self._solar_time

    @solar_time.setter
    def solar_time(self, time: Time):
        """Set the current solar time at the location."""
        self._solar_time = time

    @property
    def fi(self) -> Quantity:
        """Latitude of the location."""
        return self._fi

    @fi.setter
    def fi(self, v: Quantity) -> None:
        """Set latitude of the location."""
        self._fi = v.to('rad')
        # ensure that `fi` is always in radians

    @property
    def L_loc(self) -> Quantity:
        """Longitude of the location."""
        return self._L_loc

    @L_loc.setter
    def L_loc(self, v: Quantity) -> None:
        """Set longitude of the location."""
        self._L_loc = v

    @property
    def altitude(self) -> Quantity:
        """Altitude of location."""
        return self._altitude

    @altitude.setter
    def altitude(self, v: Quantity) -> None:
        """Set altitude of location."""
        self._altitude = v

    @property
    def climate_type(self) -> str:
        """Type of climate at the location."""
        return self._climate_type

    @climate_type.setter
    def climate_type(self, v: str) -> None:
        """Set type of climate at the location."""
        self._climate_type = v

    @property
    def timezone(self) -> str:
        """Timezone of the location."""
        return self._timezone

    @timezone.setter
    def timezone(self, v: str) -> None:
        """Set timezone of the location."""
        self._timezone = v

    @property
    def delta(self) -> Quantity:
        """Declination."""
        delta = declination(day_number(self._date))
        return Q_(delta, 'rad')

    @property
    def omega(self) -> Quantity:
        """Hour angle."""
        omega = hour_angle(time_to_decimal_hour(self.solar_time))
        return Q_(omega, 'rad')

    @property
    def omega_ss(self) -> Quantity:
        """Hour angle at sunset."""
        omega_ss = sunset_hour_angle(self._fi.m, self.delta.m)
        return Q_(omega_ss, 'rad')

    @property
    def sunset(self) -> Time:
        """Solar time at sunset."""
        t_sol = hour_angle_to_solar_time(self.omega_ss.m)
        return decimal_hour_to_time(t_sol)

    @property
    def omega_sr(self) -> Quantity:
        """Hour angle at sunrise."""
        return -self.omega_ss

    @property
    def sunrise(self) -> Time:
        """Solar time at sunrise."""
        t_sol = hour_angle_to_solar_time(self.omega_sr.m)
        return decimal_hour_to_time(t_sol)

    @property
    def daylight_time(self) -> Quantity:
        """Daylight time in decimal hours."""
        N = daylight_time(self.fi.m, self.delta.m)
        return Q_(N, 'hr')


class Sun:
    G_sc = Q_(1367.0, 'W / m ** 2')  # solar constant

    def __init__(self, location: Location):
        """
        Creates a `Sun` instance linked to a given geographic location.

        Parameters
        ----------
        location:
            Location where the solar geometry is to be determined.
        """
        self.location = location

    @property
    def theta(self) -> Quantity:
        """Zenith angle."""
        theta_z = zenith_angle(
            self.location.delta.m,
            self.location.fi.m,
            self.location.omega.m
        )
        return Q_(theta_z, 'rad')

    @property
    def gamma(self) -> Quantity:
        """
        Solar azimuth angle, with zero due south, east negative, and west 
        positive; -pi <= gamma <= pi.
        """
        gamma_s = solar_azimuth_angle(
            self.location.omega.m,
            self.theta.m,
            self.location.fi.m,
            self.location.delta.m
        )
        return Q_(gamma_s, 'rad')

    @property
    def alpha(self) -> Quantity:
        """Solar altitude angle."""
        alpha_s = solar_altitude_angle(
            self.location.delta.m,
            self.location.fi.m,
            self.location.omega.m
        )
        return Q_(alpha_s, 'rad')

    @property
    def G_on(self) -> Quantity:
        """
        Extraterrestrial irradiance incident on the plane normal to the
        radiation (beam radiation).
        """
        n = day_number(self.location.date)
        G_on = ext_irradiance_normal(n)
        return Q_(G_on, 'W / m ** 2')

    @property
    def G_o(self) -> Quantity:
        """
        Extraterrestrial irradiance incident on a horizontal plane (outside
        the atmosphere).
        """
        n = day_number(self.location.date)
        G_o = ext_irradiance_hor_surf(n, self.theta.to('rad').m)
        return Q_(G_o, 'W / m ** 2')

    @property
    def H_o(self) -> Quantity:
        """
        Integrated daily extraterrestrial irradiation (insolation)
        on a horizontal surface.
        """
        omega_sr = self.location.omega_sr.to('rad').m
        omega_ss = self.location.omega_ss.to('rad').m
        delta = self.location.delta.to('rad').m
        fi = self.location.fi.to('rad').m
        G_on = self.G_on.to('W / m ** 2').m

        # extraterrestrial irradiance on horizontal surface [W/m²] in function
        # of hour angle [rad]
        def _G_o(omega: float):
            G_o = G_on * np.cos(zenith_angle(delta, fi, omega))
            return G_o

        # integral of _G_o(omega) between sunrise and sunset hour angle has
        # units of W * rad / m²
        H_o = quad(_G_o, omega_sr, omega_ss)[0]

        # we need to scale the radians to an equivalent time duration in hours
        # to get H in units of W * hr / m²: rad * (scaling_factor in hr/rad) = hr
        t_sr = hour_angle_to_solar_time(omega_sr)
        t_ss = hour_angle_to_solar_time(omega_ss)
        k = (t_ss - t_sr) / (omega_ss - omega_sr)  # scaling factor in hr/rad
        return Q_(k * H_o, 'W * hr / m ** 2')

    @property
    def H_o_avg(self) -> Quantity:
        """
        Monthly mean daily extraterrestrial irradiation (insolation) on a
        horizontal surface for the current month in `self.location.date`.
        """
        num_days = pd.Period(self.location.date, freq='D').days_in_month
        location = Location(fi=self.location.fi)
        H_o_arr = []
        for day_num in range(num_days):
            location.date = Date(
                self.location.date.year,
                self.location.date.month,
                day_num + 1
            )
            H_o_arr.append(location.sun.H_o.to('J / m**2').magnitude)
        H_o_avg = np.mean(H_o_arr)
        return Q_(H_o_avg, 'J / m**2')

    def I_o(self, t1: Time, t2: Time) -> Quantity:
        """
        Extraterrestrial irradiation (insolation) on a horizontal surface
        between solar times `t1` and `t2`.
        """
        t1_hr = time_to_decimal_hour(t1)
        t2_hr = time_to_decimal_hour(t2)
        omega_1 = hour_angle(t1_hr)
        omega_2 = hour_angle(t2_hr)
        omega_sr = self.location.omega_sr.to('rad').m
        omega_ss = self.location.omega_ss.to('rad').m
        delta = self.location.delta.to('rad').m
        fi = self.location.fi.to('rad').m
        G_on = self.G_on.to('W / m ** 2').m

        def _G_o(omega: float):
            if omega_sr <= omega <= omega_ss:
                G_o = G_on * np.cos(zenith_angle(delta, fi, omega))
                return G_o
            else:
                return 0.0

        I_o = quad(_G_o, omega_1, omega_2)[0]
        k = (t2_hr - t1_hr) / (omega_2 - omega_1)
        return Q_(k * I_o, 'W * hr / m ** 2')

    @property
    def tau_b(self) -> Quantity:
        """
        Atmospheric transmittance for clear-sky beam radiation.

        This is the ratio of beam radiation to the extraterrestrial beam
        radiation on a horizontal or tilted surface.
        """
        A = self.location.altitude.to('km').m
        tau_b = clear_sky_beam_transmittance(
            A,
            self.location.climate_type,
            self.theta
        )
        return Q_(tau_b, 'frac')

    @property
    def G_cbn(self) -> Quantity:
        """Clear-sky beam normal irradiance."""
        return self.tau_b * self.G_on

    @property
    def G_cb(self) -> Quantity:
        """Clear-sky beam irradiance on horizontal surface."""
        return self.tau_b * self.G_o

    def I_cb(self, t1: Time, t2: Time) -> Quantity:
        """
        Clear-sky beam irradiation on horizontal surface between solar times
        `t1` and `t2`.
        """
        t1_hr = time_to_decimal_hour(t1)
        t2_hr = time_to_decimal_hour(t2)
        omega_1 = hour_angle(t1_hr)
        omega_2 = hour_angle(t2_hr)
        omega_sr = self.location.omega_sr.to('rad').m
        omega_ss = self.location.omega_ss.to('rad').m
        delta = self.location.delta.to('rad').m
        fi = self.location.fi.to('rad').m
        G_on = self.G_on.to('W / m ** 2').m

        def _G_cb(omega: float):
            if omega_sr <= omega <= omega_ss:
                theta_z = zenith_angle(delta, fi, omega)
                G_o = G_on * np.cos(theta_z)
                tau_b = clear_sky_beam_transmittance(
                    self.location.altitude.to('km').m,
                    self.location.climate_type,
                    theta_z
                )
                return tau_b * G_o
            else:
                return 0.0

        I_cb = quad(_G_cb, omega_1, omega_2)[0]
        k = (t2_hr - t1_hr) / (omega_2 - omega_1)
        return Q_(k * I_cb, 'W * hr / m ** 2')

    @property
    def tau_d(self) -> Quantity:
        """
        Atmospheric transmittance for clear-sky diffuse radiation.

        This is the ratio of diffuse radiation `G_cd` to the extraterrestrial
        beam radiation on the horizontal plane `G_o`.
        """
        tau_d = clear_sky_diffuse_transmittance(self.tau_b.to('frac').m)
        return Q_(tau_d, 'frac')

    @property
    def G_cd(self) -> Quantity:
        """Clear-sky diffuse irradiance on horizontal surface."""
        return self.tau_d * self.G_o

    def I_cd(self, t1: Time, t2: Time) -> Quantity:
        """
        Clear-sky diffuse irradiation on horizontal surface between solar
        times `t1` and `t2`.
        """
        t1_hr = time_to_decimal_hour(t1)
        t2_hr = time_to_decimal_hour(t2)
        omega_1 = hour_angle(t1_hr)
        omega_2 = hour_angle(t2_hr)
        omega_sr = self.location.omega_sr.to('rad').m
        omega_ss = self.location.omega_ss.to('rad').m
        delta = self.location.delta.to('rad').m
        fi = self.location.fi.to('rad').m
        G_on = self.G_on.to('W / m ** 2').m

        def _G_cd(omega: float):
            if omega_sr <= omega <= omega_ss:
                theta_z = zenith_angle(delta, fi, omega)
                G_o = G_on * np.cos(theta_z)
                tau_b = clear_sky_beam_transmittance(
                    self.location.altitude.to('km').m,
                    self.location.climate_type,
                    theta_z
                )
                tau_d = clear_sky_diffuse_transmittance(tau_b)
                return tau_d * G_o
            else:
                return 0.0

        I_cd = quad(_G_cd, omega_1, omega_2)[0]
        k = (t2_hr - t1_hr) / (omega_2 - omega_1)
        return Q_(k * I_cd, 'W * hr / m ** 2')

    def split_radiation(
        self,
        I: Quantity | None = None,
        H: Quantity | None = None,
        H_avg: Quantity | None = None,
        time_interval: tuple[Time, Time] | None = None
    ) -> tuple[Quantity, Quantity]:
        """
        Splits total radiation on horizontal surface into beam and diffuse
        radiation.

        Parameters
        ----------
        I: optional
            hourly total radiation on horizontal surface
        H: optional
            daily total radiation on horizontal surface
        H_avg: optional
            monthly average daily total radiation on horizontal surface
        time_interval: optional
            the time interval for which the hourly total radiation `I` is
            given; only used together with parameter `I`
        """
        if I is not None:
            I = I.to_base_units()
            I_o = self.I_o(t1=time_interval[0], t2=time_interval[1]).to_base_units()
            k_T = I / I_o
            frac_I_d = estimate_I_d_fraction(k_T)
            I_d = frac_I_d * I
            I_b = I - I_d
            return I_b, I_d
        if H is not None:
            H = H.to_base_units()
            H_o = self.H_o.to_base_units()
            K_T = H / H_o
            frac_H_d = estimate_H_d_fraction(
                K_T,
                self.location.omega_ss.to('rad').m
            )
            H_d = frac_H_d * H
            H_b = H - H_d
            return H_b, H_d
        if H_avg is not None:
            H_avg = H_avg.to_base_units()
            H_o = self.H_o.to_base_units()
            K_T_avg = H_avg / H_o
            frac_H_d_avg = estimate_avg_H_d_fraction(
                K_T_avg,
                self.location.omega_ss.to('rad').m
            )
            H_d_avg = frac_H_d_avg * H_avg
            H_b_avg = H_avg - H_d_avg
            return H_b_avg, H_d_avg


class Surface:

    def __init__(
        self,
        location: Location,
        gamma: Quantity | None,
        beta: Quantity | None
    ) -> None:
        """
        Creates a `Surface` instance at a given geographical location.

        Parameters
        ----------
        location:
            The location where the surface is situated.
        gamma:
            The azimuth angle of the surface.
        beta:
            The slope angle of the surface.
        """
        self.location = location
        if gamma is not None:
            self._gamma = gamma.to('rad')  # `.to(...)` returns new object
        else:
            self._gamma = None
        if beta is not None:
            self._beta = beta.to('rad')
        else:
            self._beta = None

    @property
    def gamma(self) -> Quantity:
        """
        Surface azimuth angle, with zero due south, east negative, and
        west positive; -pi <= gamma <= pi.
        """
        return self._gamma

    @gamma.setter
    def gamma(self, v: Quantity):
        """
        Set surface azimuth angle, with zero due south, east negative, and
        west positive; -pi <= gamma <= pi.
        """
        self._gamma = v.to('rad')
        # ensure that internally `gamma` is always in radians

    @property
    def beta(self) -> Quantity:
        """Slope angle of surface."""
        return self._beta

    @beta.setter
    def beta(self, v: Quantity) -> Quantity:
        """Set slope angle of surface."""
        self._beta = v.to('rad')
        # ensure that internally `beta` is always in radians

    @property
    def theta(self) -> Quantity:
        """Incidence angle"""
        theta = incidence_angle(
            self.beta.m, self.gamma.m, self.location.delta.m,
            self.location.fi.m, self.location.omega.m,
            self.location.sun.theta.m, self.location.sun.gamma.m
        )
        return Q_(theta, 'rad')

    @property
    def alpha_p(self) -> Quantity:
        """
        Profile angle of beam radiation.

        It is the projected solar altitude angle onto a vertical plane that
        is perpendicular to the surface in question.
        """
        alpha_p = profile_angle(
            self.location.sun.alpha.m,
            self.location.sun.gamma.m,
            self.gamma.m
        )
        return Q_(alpha_p, 'rad')

    def _find_omega_sr_ss(self, which: str) -> Quantity:
        """
        Hour angle at which sunrise or sunset occurs on the surface.

        Notes
        -----
        When sunrise occurs on the surface, the incidence angle is 90°, i.e.
        the sun rays are parallel with the plane of the surface.
        """

        def eq(omega: float) -> float:
            theta = incidence_angle(
                self.beta.m,
                self.gamma.m,
                self.location.delta.m,
                self.location.fi.m,
                omega,
                self.location.sun.theta.m,
                self.location.sun.gamma.m
            )
            return np.cos(theta)

        if which == 'sr':  # sunrise
            omega_0 = hour_angle(0)
            omega_1 = hour_angle(12)  # solar noon
        else:  # sunset
            omega_0 = hour_angle(12)
            omega_1 = hour_angle(24)
        sol = root_scalar(eq, bracket=[omega_0, omega_1])
        return Q_(sol.root, 'rad')

    @property
    def omega_sr(self) -> Quantity:
        """Hour angle at which sunrise occurs on the surface."""
        omega_sr_surf = self._find_omega_sr_ss('sr')
        omega_sr_loc = self.location.omega_sr
        return max(omega_sr_surf, omega_sr_loc)
        # the sunrise on the surface cannot happen earlier than the sunrise on
        # the location

    @property
    def omega_ss(self) -> Quantity:
        """Hour angle at which sunset occurs on the surface."""
        omega_ss_surf = self._find_omega_sr_ss('ss')
        omega_ss_loc = self.location.omega_ss
        return min(omega_ss_surf, omega_ss_loc)
        # if the sunset at the location happens earlier, this will also be the
        # actual sunset on the surface.

    @property
    def sunrise(self) -> Time:
        """Solar time at which sunrise occurs on the surface."""
        t_sol = hour_angle_to_solar_time(self.omega_sr.m)
        return decimal_hour_to_time(t_sol)

    @property
    def sunset(self) -> Time:
        """Solar time at which sunset occurs on the surface."""
        t_sol = hour_angle_to_solar_time(self.omega_ss.m)
        return decimal_hour_to_time(t_sol)

    @property
    def R_b(self) -> Quantity:
        """
        The ratio of beam radiation on the tilted surface to that on a
        horizontal plane, i.e. the gain in incident beam radiation by tilting
        a receiving surface, at the currently set date and solar time.
        """
        return np.cos(self.theta) / np.cos(self.location.sun.theta)

    def R_b_ave(self, t1: Time, t2: Time) -> Quantity:
        """
        The average ratio of beam radiation on the tilted surface to that on
        a horizontal plane between the given solar times `t1` and `t2` on the
        currently set date.
        """
        delta = self.location.delta.to('rad').m
        fi = self.location.fi.to('rad').m
        beta = self.beta.to('rad').m
        gamma = self.gamma.to('rad').m
        omega1 = hour_angle(time_to_decimal_hour(t1))
        omega2 = hour_angle(time_to_decimal_hour(t2))
        a = (
            np.sin(delta) * np.sin(fi) * np.cos(beta)
            - np.sin(delta) * np.cos(fi) * np.sin(beta) * np.cos(gamma)
        ) * (omega2 - omega1)
        a += (
            np.cos(delta) * np.cos(fi) * np.cos(beta)
            + np.cos(delta) * np.sin(fi) * np.sin(beta) * np.cos(gamma)
        ) * (np.sin(omega2) - np.sin(omega1))
        a -= (
            np.cos(delta) * np.sin(beta) * np.sin(gamma)
        ) * (np.cos(omega2) - np.cos(omega1))
        b = (
            np.cos(fi) * np.cos(delta) * (np.sin(omega2) - np.sin(omega1))
            + (np.sin(fi) * np.sin(delta) * (omega2 - omega1))
        )
        R_b_ave = max(0.0, a / b)
        return Q_(R_b_ave, 'frac')

    @property
    def G_Tcb(self) -> Quantity:
        """Clear-sky beam radiation (irradiance) on tilted surface."""
        return self.location.sun.G_cbn * np.cos(self.theta.to('rad').m)

    def estimate_I_Tb(self, H: Quantity, t1: Time, t2: Time) -> Quantity:
        """
        Estimate hourly beam radiation on tilted surface at solar times `t1`
        and `t2` when daily total radiation on horizontal surface `H` is given.
        """
        H = H.to('MJ / m ** 2').m
        H_o = self.location.sun.H_o.to('MJ / m ** 2').m
        K_T = H / H_o
        omega_ss = self.location.omega_ss.to('rad').m
        omega = self.location.omega.to('rad').m
        frac_H_d = estimate_H_d_fraction(
            K_T=K_T,
            omega_ss=omega_ss
        )
        H_d = frac_H_d * H
        r_t = estimate_r_t(
            omega=omega,
            omega_ss=omega_ss
        )
        r_d = estimate_r_d(
            omega=omega,
            omega_ss=omega_ss
        )
        I = Q_(r_t * H, 'MJ / m ** 2')
        I_d = Q_(r_d * H_d, 'MJ / m ** 2')
        I_b = I - I_d
        R_b = self.R_b_ave(t1, t2)
        I_bT = R_b * I_b
        return I_bT

    @property
    def H_To(self) -> Quantity:
        """
        Integrated daily extraterrestrial irradiation (insolation) on tilted 
        surface.
        """
        omega_sr = self.omega_sr.to('rad').m
        omega_ss = self.omega_ss.to('rad').m
        beta = self.beta.to('rad').m
        gamma = self.gamma.to('rad').m
        delta = self.location.delta.to('rad').m
        fi = self.location.fi.to('rad').m
        G_on = self.location.sun.G_on.to('W / m ** 2').m

        # extraterrestrial irradiance on tilted surface [W/m²] as function
        # of hour angle [rad]
        def _G_To(omega: float):
            theta_z = zenith_angle(delta, fi, omega)
            gamma_s = solar_azimuth_angle(omega, theta_z, fi, delta)
            G_To = G_on * np.cos(incidence_angle(
                beta, gamma, delta, fi, omega, theta_z, gamma_s
            ))
            return G_To

        # integral of _G_o(omega) between sunrise and sunset hour angle has
        # units of W * rad / m²
        H_To = quad(_G_To, omega_sr, omega_ss)[0]

        # we need to scale the radians to an equivalent time duration in hours
        # to get H in units of W * hr / m²: rad * (scaling_factor in hr/rad) = hr
        t_sr = hour_angle_to_solar_time(omega_sr)
        t_ss = hour_angle_to_solar_time(omega_ss)
        k = (t_ss - t_sr) / (omega_ss - omega_sr)  # scaling factor in hr/rad
        return Q_(k * H_To, 'W * hr / m ** 2')


class TrackingSurfaceEWS(Surface):
    """
    Tracking surface rotated about a horizontal east-west axis with a single
    daily adjustment so that beam radiation is normal to the surface at noon
    each day.
    """

    def __init__(self, location: Location):
        k = location.fi - location.delta
        beta = np.abs(k)
        gamma = Q_(np.pi, 'rad') if k <= 0 else Q_(0.0, 'rad')
        super().__init__(location, gamma, beta)


class TrackingSurfaceEWC(Surface):
    """
    Tracking surface rotated about a horizontal east-west axis with continuous 
    adjustment to minimize the angle of incidence.
    """

    def __init__(self, location: Location):
        super().__init__(location, None, None)

    @property
    def gamma(self) -> Quantity:
        """Surface azimuth angle."""
        gamma_s = np.abs(self.location.sun.gamma)
        if gamma_s < np.pi / 2:
            return Q_(0.0, 'rad')
        else:
            return Q_(np.pi, 'rad')

    @property
    def beta(self) -> Quantity:
        """Slope angle of surface."""
        beta = np.arctan(
            np.tan(self.location.sun.theta)
            * np.abs(np.cos(self.location.sun.gamma))
        )
        return beta


class TrackingSurfaceNSC(Surface):
    """
    Tracking surface rotated about a horizontal north-south axis with continuous
    adjustment to minimize the angle of incidence.
    """

    def __init__(self, location: Location):
        super().__init__(location, None, None)

    @property
    def gamma(self) -> Quantity:
        """Surface azimuth angle."""
        gamma_s = self.location.sun.gamma
        if gamma_s > 0:
            return Q_(np.pi / 2, 'rad')
        else:
            return Q_(-np.pi / 2, 'rad')

    @property
    def beta(self) -> Quantity:
        """Slope angle of surface."""
        beta = np.arctan(
            np.tan(self.location.sun.theta)
            * np.abs(np.cos(self.gamma - self.location.sun.gamma))
        )
        return beta


class TrackingSurfaceFSV(Surface):
    """
    Tracking surface with a fixed slope rotated about a vertical axis.
    The angle of incidence is minimized when the surface azimuth angle tracks
    the solar azimuth angle.
    """

    def __init__(self, location: Location, beta: Quantity):
        super().__init__(location, None, beta)

    @property
    def gamma(self) -> Quantity:
        """Surface azimuth angle."""
        return self.location.sun.gamma


class TrackingSurfaceNSE(Surface):
    """
    Tracking surface rotated about a north-south axis parallel to the
    earth's axis with continuous adjustment to minimize the incidence angle.
    """

    def __init__(self, location: Location):
        super().__init__(location, None, None)

    @property
    def gamma(self) -> Quantity:
        """Surface azimuth angle."""
        theta_q = (
                np.cos(self.location.sun.theta) * np.cos(self.location.fi)
                + np.sin(self.location.sun.theta) * np.sin(self.location.fi)
                * np.cos(self.location.sun.gamma)
        )
        k = np.arctan(
            np.sin(self.location.sun.theta)
            * np.sin(self.location.sun.gamma)
            / (np.cos(theta_q) * np.sin(self.location.fi))
        ) * self.location.sun.gamma
        if k >= 0.0:
            C1 = 0
        else:
            C1 = 1
        if self.location.sun.gamma >= 0.0:
            C2 = 1
        else:
            C2 = -1
        gamma = np.arctan(
            np.sin(self.location.sun.theta)
            * np.sin(self.location.sun.gamma)
            / (np.cos(theta_q) * np.sin(self.location.fi))
        ) + np.pi * C1 * C2
        return gamma

    @property
    def beta(self) -> Quantity:
        """Slope angle of surface."""
        beta = np.arctan(np.tan(self.location.fi) / np.cos(self.gamma))
        return beta


class TrackingSurface2Axes(Surface):
    """
    Tracking surface that is continuously tracking about two axes to minimize
    the angle of incidence.
    """
    def __init__(self, location: Location):
        super().__init__(location, None, None)

    @property
    def gamma(self) -> Quantity:
        """Surface azimuth angle."""
        return self.location.sun.gamma

    @property
    def beta(self) -> Quantity:
        """Slope angle of surface."""
        return self.location.sun.theta
