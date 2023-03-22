"""
Estimate solar radiation incident on surfaces with arbitrary tilt and orientation.

Solar radiation on tilted surfaces comes from 3 components:
- direct radiation
- diffuse sky radiation
- diffuse radiation reflected from the ground
"""

from typing import Optional, Dict, List, Any
import math
from datetime import datetime as DateTime
from hvac import Quantity
from .path import HorizonProfile
from .location import Location

Q_ = Quantity


class Surface:
    """
    Define a surface.

    Attributes
    ----------
    azimuth : Quantity
        Azimuth angle of the surface measured clockwise from North (0°).
    tilt : Quantity
        Tilt angle of the surface.
    width : Quantity, optional
        Width of the surface. Default value is 1 m.
    height : Quantity, optional
        Height of the surface. Default value is 1 m
    area : Quantity
        Area of the surface.
    horizon_profile : HorizonProfile, optional
        Horizon profile as seen from the position of the surface.
    """
    def __init__(
        self,
        azimuth: Quantity,
        tilt: Quantity,
        width: Quantity = Q_(1.0, 'm'),
        height: Quantity = Q_(1.0, 'm'),
        horizon_profile: Optional[HorizonProfile] = None,
        area: Optional[Quantity] = None
    ):
        self.azimuth = azimuth
        self.tilt = tilt
        self.width = width
        self.height = height
        self.horizon_profile = horizon_profile
        self.area = area or (width * height)


class TranspositionModel:
    """
    Calculate irradiance on a titled surface when solar irradiance on the
    horizontal surface is known.
    """
    F_sky: float = float('nan')  # correction factor for diffuse sky radiation
    F_grd: float = float('nan')  # correction factor for diffuse radiation reflected from the ground

    def __init__(
        self,
        loc: Location,
        datetime: DateTime,
        surface: Surface,
        I_beam: Quantity,
        I_dif: Quantity,
        I_glo_hor: Optional[Quantity] = None,
        rho_grd: float = 0.2,
        pv: bool = False
    ):
        """
        Parameters
        ----------
        loc: Location
            The geographic location under consideration.
        datetime : DateTime
            Local standard time at location where solar irradiance is investigated.
        surface : Surface
            Surface under investigation.
        I_beam : Quantity
            Solar irradiance along the direction of the sun rays.
        I_dif : Quantity
            Diffuse irradiance from the sky on the horizontal surface.
        I_glo_hor: Quantity, optional
            Global irradiance on the horizontal surface. If `None` it will be
            derived from `I_beam` and `I_dif`.
        rho_grd : float, optional
            Reflectivity of the ground (ground albedo). Default value is 0.2.
        pv : bool, optional
            To indicate that the surface is a PV panel. Default is `False`.
            If `True`, the global irradiance effectively absorbed by the PV
            panel is determined by using the incidence angle factors.
        """
        self.loc = loc
        self._date = datetime.date()
        self._datetime = datetime
        self.surface = surface
        self.I_beam = I_beam
        self.I_dif = I_dif
        self._I_glo_hor = I_glo_hor
        self.rho_grd = rho_grd
        self.pv = pv

    @property
    def theta_i(self) -> Quantity:
        """Incidence angle of the sun on tilted surface."""
        sp = self.loc.sun_position(self._datetime)
        fi_sun = sp.azimuth.to('rad').m
        theta_sun = sp.zenith.to('rad').m
        fi_sur = self.surface.azimuth.to('rad').m
        theta_sur = self.surface.tilt.to('rad').m
        a = math.sin(theta_sun) * math.sin(theta_sur) * math.cos(fi_sun - fi_sur)
        b = math.cos(theta_sun) * math.cos(theta_sur)
        cos_i = a + b
        if 0.0 <= cos_i <= 1.0:
            theta_i = math.acos(cos_i)
        else:
            theta_i = math.pi / 2.0
        return Q_(theta_i, 'rad')

    @property
    def I_glo_hor(self) -> Quantity:
        """Global irradiance on horizontal surface."""
        if self._I_glo_hor is None:
            theta_sun = self.loc.sun_position(self._datetime).zenith.to('rad').m
            I_glo_hor = self.I_beam * math.cos(theta_sun) + self.I_dif
            return I_glo_hor
        else:
            return self._I_glo_hor

    @property
    def I_dir_sur(self) -> Quantity:
        """Direct irradiance on tilted surface (i.e. normal to the surface)."""
        theta_i = self.theta_i.to('rad').m
        I_dir = self.I_beam * math.cos(theta_i)
        if self.pv is True:
            I_dir *= self.K_ta_dir
        return I_dir

    @property
    def I_dif_sky(self) -> Quantity:
        """Sky-reflected diffuse radiation on tilted surface."""
        I_dif_sky = self.F_sky * self.I_dif
        if self.pv is True:
            I_dif_sky *= self.K_ta_dif
        return I_dif_sky

    @property
    def I_dif_grd(self) -> Quantity:
        """Ground-reflected diffuse radiation on tilted surface."""
        I_dif_grd = self.F_grd * self.rho_grd * self.I_glo_hor
        if self.pv is True:
            I_dif_grd *= self.K_ta_grd
        return I_dif_grd

    @property
    def I_glo_sur(self) -> Quantity:
        """
        Global radiation on the tilted surface.

        This is the sum of the direct, sky-reflected diffuse and
        ground-reflected diffuse components.
        """
        I_glo_surf = self.I_dir_sur + self.I_dif_sky + self.I_dif_grd
        return I_glo_surf

    @classmethod
    def daily_profile(
        cls,
        location: Location,
        surface: Surface,
        irradiance_profile_hor: Dict[str, List[Any]],
        rho_grd: float = 0.2,
        pv: bool = False
    ) -> Dict[str, List[Any]]:
        """
        Returns daily profile of:
        - incidence angle of the sun on the tilted surface,
        - the direct, diffuse and global irradiance on the tilted surface
        when the daily profile of beam irradiance and diffuse
        irradiance on the horizontal surface is given.
        """
        theta_i_profile = []
        I_glo_sur_profile = []
        I_dir_sur_profile = []
        I_dif_sur_profile = []
        t_ax = irradiance_profile_hor['t']
        I_beam_ax = irradiance_profile_hor['beam']
        I_dif_ax = irradiance_profile_hor['dif']
        for t, I_beam, I_dif in zip(t_ax, I_beam_ax, I_dif_ax):
            tpm = cls(location, t, surface, I_beam, I_dif, rho_grd=rho_grd, pv=pv)
            theta_i_profile.append(tpm.theta_i)
            I_glo_sur_profile.append(tpm.I_glo_sur)
            I_dir_sur_profile.append(tpm.I_dir_sur)
            I_dif_sur_profile.append(tpm.I_dif_sky + tpm.I_dif_grd)
        d = {
            't': t_ax,
            'theta_i': theta_i_profile,
            'glo_sur': I_glo_sur_profile,
            'dir_sur': I_dir_sur_profile,
            'dif_sur': I_dif_sur_profile,
        }
        return d

    @staticmethod
    def _calculate_transmittance(theta_i: float):
        # theta_i = incidence angle of sun rays
        n_air = 1.0  # index of refraction of air
        n_pan = 1.1  # index of refraction of glass panel
        sin_theta_r = (n_air / n_pan) * math.sin(theta_i)
        theta_r = math.asin(sin_theta_r)  # angle of refraction
        A = math.sin(theta_r - theta_i) ** 2
        B = math.sin(theta_r + theta_i) ** 2
        C = math.tan(theta_r - theta_i) ** 2
        D = math.tan(theta_r + theta_i) ** 2
        try:
            tau_r = 1.0 - 0.5 * ((A / B) + (C / D))  # transmittance of glass cover
        except ZeroDivisionError:
            tau_r = 1.0
        K = 4.0  # 1/m
        t_pan = 0.002  # thickness of glass panel cover, m
        tau_a = math.exp(-K * t_pan / math.cos(theta_r))  # absorptance factor of glass
        tau = tau_a * tau_r  # fraction of solar radiation that reaches the PV cells
        return tau

    def _calculate_incidence_angle_factor(self, theta_i: Quantity) -> Quantity:
        tau_0 = self._calculate_transmittance(0.0)
        tau = self._calculate_transmittance(theta_i.to('rad').m)
        K = tau / tau_0
        return Q_(K, 'frac')

    @property
    def K_ta_dir(self) -> Quantity:
        """Incidence angle factor for direct radiation."""
        return self._calculate_incidence_angle_factor(self.theta_i)

    @property
    def K_ta_dif(self) -> Quantity:
        """Incidence angle factor for diffuse sky radiation."""
        return self._calculate_incidence_angle_factor(Q_(58.0, 'deg'))

    @property
    def K_ta_grd(self) -> Quantity:
        """Incidence angle factor for ground reflected radiation."""
        return self._calculate_incidence_angle_factor(Q_(58.0, 'deg'))


class IsotropicSkyModel(TranspositionModel):
    """
    Radiation from the sky is isotropic, i.e. sky-diffuse and ground-reflected solar components are uniformly
    distributed over the sky dome.
    """
    def __init__(
            self,
            loc: Location,
            datetime: DateTime,
            surface: Surface,
            I_beam: Quantity,
            I_dif: Quantity,
            I_glo_hor: Optional[Quantity] = None,
            rho_grd: float = 0.2,
            pv: bool = False
    ):
        super().__init__(loc, datetime, surface, I_beam, I_dif, I_glo_hor, rho_grd, pv)
        theta_p = self.surface.tilt.to('rad').m
        self.F_sky = (1 + math.cos(theta_p)) / 2.0
        self.F_grd = (1 - math.cos(theta_p)) / 2.0


class AnisotropicSkyModel(TranspositionModel):
    """
    Implements the ASHRAE anisotropic sky model, recommended for building load calculations. It applies only to
    clear-sky conditions and should not be used for cloudy conditions.
    """
    def __init__(
            self,
            loc: Location,
            datetime: DateTime,
            surface: Surface,
            I_beam: Quantity,
            I_dif: Quantity,
            I_glo_hor: Optional[Quantity] = None,
            rho_grd: float = 0.2,
            pv: bool = False
    ):
        super().__init__(loc, datetime, surface, I_beam, I_dif, I_glo_hor, rho_grd, pv)
        theta_p = self.surface.tilt.to('rad').m
        theta_i = self.theta_i.to('rad').m
        Y = max(0.45, 0.55 + 0.437 * math.cos(theta_i) + 0.313 * math.cos(theta_i) ** 2)
        if theta_p <= math.pi / 2:
            self.F_sky = Y * math.sin(theta_p) + math.cos(theta_p)
        else:
            self.F_sky = Y * math.sin(theta_p)
        self.F_grd = (1 - math.cos(theta_p)) / 2.0
