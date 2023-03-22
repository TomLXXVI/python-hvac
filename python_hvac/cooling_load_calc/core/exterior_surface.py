from typing import Any, Callable
import math
from scipy.interpolate import interp1d
from hvac import Quantity
from hvac.climate import ClimateData
from hvac.climate.sun import Surface, AnisotropicSkyModel
from hvac.climate.sun.solar_time import time_to_decimal_hour


Q_ = Quantity


class ExteriorSurface(Surface):
    """
    Intended for internal use only. It is an attribute of `ExteriorBuildingElement`
    and `Window`. Used for:
    - calculating the daily hourly values of global solar irradiance incident on
    the surface using the `AnisotropicSkyModel`
    - calculating the daily hourly average sol-air temperatures on the exterior
    surface of a building element
    """

    def __init__(
        self,
        azimuth: Quantity,
        tilt: Quantity,
        width: Quantity,
        height: Quantity,
        climate_data: ClimateData,
        surface_resistance: Quantity | None = None,
        surface_absorptance: Quantity | None = None,
        surface_color: str = 'dark-colored'
    ) -> None:
        """
        Create `ExteriorSurface` object.

        Parameters
        ----------
        azimuth: Quantity
            The azimuth angle of the surface measured clockwise from due North.
            (N = 0°, E = 90°, S = 180°, and W = 270°)
        tilt: Quantity
            The inclination angle of the surface with respect to the horizontal
            surface (0° <= tilt <= 90°).
        width: Quantity
            The width of the surface.
        height: Quantity
            The height of the surface.
        climate_data: ClimateData
            `ClimateData`-object encapsulating the daily temperature and solar
            irradiance profiles at the given design day on the location under
            investigation.
        surface_resistance: Quantity, optional
            The thermal resistance of the exterior surface film used to calculate
            the sol-air temperature of the surface.
        surface_absorptance: Quantity, optional
            The absorptance of the exterior surface film used to calculate the
            sol-air temperature of the surface.
        surface_color: ['dark-colored' (default), 'light-colored']
            Indicate if the surface is either dark-colored or light-colored. You
            can use this instead of specifying `surface_resistance` and
            `surface_absorptance`.
        """
        super().__init__(azimuth, tilt, width, height)
        self.climate_data = climate_data

        # Get the daily profile of the position of the sun (hourly values) as a
        # dict {
        #   't': List[datetime],
        #   'azimuth': List[Quantity],
        #   'elevation': List[Quantity]
        # }
        self.sun_pos_profile = self.climate_data.sun_position_profile

        # Create a function that takes the hour of the day in decimal hours
        # (float) and returns the elevation angle of the sun in degrees (float)
        # at that hour of the day.
        self._sun_alt_fun = self._interpolate_sun_altitude()

        # Create a function that takes the hour of the day in decimal hours
        # (float) and returns the solar-surface azimuth angle in degrees (float)
        # at that hour of the day.
        self._sun_surf_azi_fun = self._interpolate_solar_surface_azimuth()

        # Get the daily profile of the solar incidence angle and solar
        # irradiance on the surface (hourly values) as a dict {
        #   't': List[datetime],
        #   'theta_i': List[Quantity],
        #   'glo_sur': List[Quantity],
        #   'dir_sur': List[Quantity],
        #   'dif_sur': List[Quantity]
        # }
        self.irr_profile = AnisotropicSkyModel.daily_profile(
            location=climate_data.location,
            surface=self,
            irradiance_profile_hor=climate_data.irr_profile
        )

        # Create a function that takes the hour of the day in decimal hours
        # (float) and returns the incidence angle of the sun rays on the surface
        # in degrees (float) at that hour of the day.
        self._theta_i_fun = self._interpolate_incidence_angle()

        # Create a function that takes the hour of the day in decimal hours
        # (float) and returns the direct irradiance on the surface in W/m²
        # (float) at that hour of the day.
        self._I_dir_fun = self._interpolate_direct_irradiance()

        # Create a function that takes the hour of the day in decimal hours
        # (float) and returns the diffuse irradiance on the surface in W/m²
        # (float) at that hour of the day.
        self._I_dif_fun = self._interpolate_diffuse_irradiance()
        self._I_glo_fun = self._interpolate_global_irradiance()

        if None not in [surface_resistance, surface_absorptance]:
            # Calculate the daily profile of sol-air temperature of the surface
            # (hourly values) as a dict {
            #   't': List[datetime],
            #   'T': List[Quantity]
            # }
            self.T_sol_profile = self._get_sol_air_temperature_profile(
                R_surf=surface_resistance,
                a_surf=surface_absorptance
            )
        else:
            self.T_sol_profile = self._get_sol_air_temperature_profile(
                surface_color=surface_color
            )
        # Create a function that takes the hour of the day in decimal hours
        # (float) and returns the sol-air temperature in degC (float) at that
        # hour of the day.
        self._T_sol_fun = self._interpolate_sol_air_temperature()

    def _get_sol_air_temperature_profile(
        self,
        R_surf: Quantity | None = None,
        a_surf: Quantity | None = None,
        surface_color: str = 'dark-colored'
    ) -> dict[str, list[Any]]:
        tilt = self.tilt.to('rad').m
        dT_sky = 3.9 * math.cos(tilt)
        T_arr = [
            T.to('degC').m
            for T in self.climate_data.Tdb_profile['T']
        ]
        I_glo_arr = [
            I.to('W / m ** 2').m
            for I in self.irr_profile['glo_sur']
        ]
        if None not in [R_surf, a_surf]:
            R_surf = R_surf.to('m ** 2 * K / W').m
            a_surf = a_surf.to('frac').m
            Tsol_arr = [
                T + a_surf * R_surf * I_glo - dT_sky
                for T, I_glo in zip(T_arr, I_glo_arr)
            ]
        elif surface_color == 'light-colored':
            Tsol_arr = [
                T + 0.026 * I_glo - dT_sky
                for T, I_glo in zip(T_arr, I_glo_arr)
            ]
        else:
            Tsol_arr = [
                T + 0.052 * I_glo - dT_sky
                for T, I_glo in zip(T_arr, I_glo_arr)
            ]
        d = {
            't': self.irr_profile['t'],
            'T': [Q_(Tsol, 'degC') for Tsol in Tsol_arr]
        }
        return d

    def _interpolate_sol_air_temperature(self) -> Callable[[float], float]:
        t_ax = [
            time_to_decimal_hour(dt.time()) * 3600.0
            for dt in self.T_sol_profile['t']
        ]
        T_ax = [
            T.to('degC').m
            for T in self.T_sol_profile['T']
        ]
        interpolant = interp1d(x=t_ax, y=T_ax, kind='cubic')
        return interpolant

    def _interpolate_incidence_angle(self) -> Callable[[float], float]:
        t_ax = [
            time_to_decimal_hour(dt.time()) * 3600.0
            for dt in self.irr_profile['t']
        ]
        theta_i_ax = [
            theta_i.to('rad').m
            for theta_i in self.irr_profile['theta_i']
        ]
        interpolant = interp1d(x=t_ax, y=theta_i_ax)
        return interpolant

    def _interpolate_direct_irradiance(self) -> Callable[[float], float]:
        t_ax = [
            time_to_decimal_hour(dt.time()) * 3600.0
            for dt in self.irr_profile['t']
        ]
        I_dir_ax = [
            I_dir.to('W / m ** 2').m
            for I_dir in self.irr_profile['dir_sur']
        ]
        interpolant = interp1d(x=t_ax, y=I_dir_ax)
        return interpolant

    def _interpolate_diffuse_irradiance(self) -> Callable[[float], float]:
        t_ax = [
            time_to_decimal_hour(dt.time()) * 3600.0
            for dt in self.irr_profile['t']
        ]
        I_dif_ax = [
            I_dif.to('W / m ** 2').m
            for I_dif in self.irr_profile['dif_sur']
        ]
        interpolant = interp1d(x=t_ax, y=I_dif_ax)
        return interpolant

    def _interpolate_global_irradiance(self) -> Callable[[float], float]:
        t_ax = [
            time_to_decimal_hour(dt.time()) * 3600.0
            for dt in self.irr_profile['t']
        ]
        I_glo_ax = [
            I_glo.to('W / m ** 2').m
            for I_glo in self.irr_profile['glo_sur']
        ]
        interpolant = interp1d(x=t_ax, y=I_glo_ax)
        return interpolant

    def _interpolate_sun_altitude(self) -> Callable[[float], float]:
        t_ax = [
            time_to_decimal_hour(dt.time()) * 3600.0
            for dt in self.sun_pos_profile['t']
        ]
        sun_alt_ax = [
            alt.to('rad').m
            for alt in self.sun_pos_profile['elevation']
        ]
        interpolant = interp1d(x=t_ax, y=sun_alt_ax)
        return interpolant

    def _interpolate_solar_surface_azimuth(self) -> Callable[[float], float]:
        t_ax = [
            time_to_decimal_hour(dt.time()) * 3600.0
            for dt in self.sun_pos_profile['t']
        ]
        sun_surf_azi_ax = [
            (sun_azi - self.azimuth).to('rad').m
            for sun_azi in self.sun_pos_profile['azimuth']
        ]
        interpolant = interp1d(x=t_ax, y=sun_surf_azi_ax)
        return interpolant

    def T_sol(self, t: float) -> float:
        """
        Get sol-air temperature in degC at time t in seconds from 00:00:00.
        """
        return self._T_sol_fun(t)

    def theta_i(self, t: float) -> float:
        """
        Get incidence angle of sun rays on surface in radians at time t in
        seconds from 00:00:00.
        """
        return self._theta_i_fun(t)

    def I_dir(self, t: float) -> float:
        """
        Get solar direct irradiance on surface in W/m² at time t in seconds from
        00:00:00.
        """
        return float(self._I_dir_fun(t))

    def I_dif(self, t: float) -> float:
        """
        Get solar diffuse irradiance on surface in W/m² at time t in seconds from
        00:00:00.
        """
        return float(self._I_dif_fun(t))

    def I_glo(self, t: float) -> float:
        """
        Get solar global irradiance on surface in W/m² at time t in seconds from
        00:00:00.
        """
        return float(self._I_glo_fun(t))

    def sun_altitude(self, t: float) -> float:
        """
        Get sun altitude angle in radians at time t in seconds from 00:00:00.
        """
        return self._sun_alt_fun(t)

    def sun_surface_azimuth(self, t: float) -> float:
        """
        Get the solar-surface azimuth angle in radians at time t in seconds
        from 00:00:00.
        """
        return self._sun_surf_azi_fun(t)
