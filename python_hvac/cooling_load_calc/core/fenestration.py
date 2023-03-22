from typing import cast
from dataclasses import dataclass
import shelve
import numpy as np
import pandas as pd
from hvac import Quantity
from hvac.climate import ClimateData
from hvac.climate.sun.solar_time import time_from_decimal_hour
from hvac.cooling_load_calc.core.exterior_surface import ExteriorSurface


Q_ = Quantity


class WindowThermalProperties:
    db_path: str

    class SolarHeatGainCoefficient:
        """
        The solar heat gain coefficient (SHGC) takes into account the solar
        radiation that is transmitted directly through the window and the
        solar radiation that is first absorbed by the window and then released
        to the interior.
        """

        def __init__(self, cog_dir: dict[float, float], cog_dif: float, wnd: float):
            self._cog_dir = (list(cog_dir.keys()), list(cog_dir.values()))
            self.cog_dif = cog_dif
            self.wnd = wnd

        def cog_dir(self, theta_i: float) -> float:
            SHGC_dir = np.interp(theta_i, self._cog_dir[0], self._cog_dir[1])
            return SHGC_dir

    def __init__(
        self,
        ID: str,
        U: Quantity,
        SHGC_cog_dir: dict[float, float],
        SHGC_cog_dif: float,
        SHGC_wnd: float
    ) -> None:
        """
        Initialize a `WindowThermalProperties` instance.

        Parameters
        ----------
        ID:
            Identifier for the type of window
        U:
            Overall U value of the entire window, including edge effects and frame
            (see ASHRAE Fundamentals 2017, Ch. 15, Table 4).
        SHGC_cog_dir:
            Dictionary of which the keys are solar incidence angles and the
            values the corresponding SHGC value (see ASHRAE Fundamentals 2017,
            Ch. 15, Table 10).
        SHGC_cog_dif:
            The SHGC for diffuse solar radiation (see ASHRAE Fundamentals 2017,
            Ch. 15, Table 10).
        SHGC_wnd:
            The SHGC of the entire window at normal incidence (see ASHRAE
            Fundamentals 2017, Ch. 15, Table 10).
        """
        self.ID = ID
        self.U = U
        self.SHGC = self.SolarHeatGainCoefficient(SHGC_cog_dir, SHGC_cog_dif, SHGC_wnd)

    def save(self):
        with shelve.open(self.db_path) as shelf:
            shelf[self.ID] = self

    @classmethod
    def load(cls, ID: str) -> 'WindowThermalProperties':
        with shelve.open(cls.db_path) as shelf:
            return cast('WindowThermalProperties', shelf[ID])

    @classmethod
    def overview(cls) -> list[str]:
        """Returns a list with all the ID's of objects stored in the shelf."""
        with shelve.open(cls.db_path) as shelf:
            return list(shelf.keys())

    def __str__(self):
        str_ = f"Window-type: {self.ID}\n"
        str_ += '-' * len(str_[:-1]) + '\n'
        str_ += f"\tU entire window: {self.U.to('W / (m ** 2 * K)'):~P.1f}\n"
        str_ += f"\tSHGC(0°) direct radiation: {self.SHGC.cog_dir(0):.1f}\n"
        str_ += f"\tSHGC diffuse radiation: {self.SHGC.cog_dif:.1f}\n"
        str_ += f"\tSHGC(0°) entire window: {self.SHGC.wnd:.1f}\n"
        return str_


@dataclass
class ExteriorShadingDevice:
    """
    Overhang over window or a recessed window.

    Parameters
    ----------
    vertical_projection: default 0 m
        Projection distance of vertical part of exterior shading device.
    horizontal_projection: default 0 m
        Projection distance of horizontal part of exterior shading device.
    width_offset : default 0 m
        Distance between vertical window edge and vertical side of exterior
        shading device.
    height_offset : default 0 m
        Distance between window's upper edge and the horizontal underside of
        exterior shading device.
    """
    vertical_projection: Quantity = Q_(0, 'm')
    horizontal_projection: Quantity = Q_(0, 'm')
    width_offset: Quantity = Q_(0, 'm')
    height_offset: Quantity = Q_(0, 'm')


@dataclass
class InteriorShadingDevice:
    """
    Louvered shades, roller shades, draperies or insect screens
    (see ASHRAE Fundamentals 2017, Ch. 15, §5.2 and tables 14A to 14G).

    Parameters
    ----------
    IAC_dif:
        Interior attenuation coefficient for diffuse solar radiation.
    F_rad:
        Radiative fraction of solar heat gain through window in combination with
        interior shading device.
    IAC_0 : default None
        Interior attenuation coefficient for direct solar radiation at normal
        incidence angle (0°).
    IAC_60: default None
        Interior attenuation coefficient for direct solar radiation at 60°
        incidence angle.
    louver_orient: default 'horizontal'
        Sets the orientation of the louvers: either 'horizontal', or 'vertical'.

    Notes
    -----
    Only in case of louvered shades 3 IACs are considered (IAC_0, IAC_60, and
    IAC_dif). In all other cases only 1 IAC is provided that applies both to
    direct and diffuse solar radiation. In these cases it suffices to provide
    only a value to parameter `IAC_dif`.
    """
    IAC_dif: float
    F_rad: float
    IAC_0: float | None = None
    IAC_60: float | None = None
    louver_orient: str = 'horizontal'


class Window:

    def __init__(self):
        self.ID: str = ''
        self._exterior_surface: ExteriorSurface | None = None
        self.therm_props: WindowThermalProperties | None = None
        self.F_rad: float | None = None
        self._ext_shading_dev: ExteriorShadingDevice | None = None
        self._int_shading_dev: InteriorShadingDevice | None = None

    @classmethod
    def create(
        cls,
        ID: str,
        azimuth: Quantity,
        tilt: Quantity,
        width: Quantity,
        height: Quantity,
        climate_data: ClimateData,
        therm_props: WindowThermalProperties,
        F_rad: float = 0.46,
        ext_shading_dev: ExteriorShadingDevice | None = None,
        int_shading_dev: InteriorShadingDevice | None = None
    ) -> 'Window':
        """
        Create `Window` instance.

        Parameters
        ----------
        ID:
            Useful name to identify the `Window` instance.
        azimuth:
            Azimuth angle of the window.
        tilt:
            Tilt angle of the window.
        width:
            Width of the window opening.
        height:
            Height of the window opening.
        climate_data:
            Instance of `ClimateData` class, containing the climatic design data.
        therm_props:
            See class `WindowThermalProperties`. The thermal and solar properties
            of the window.
        F_rad: default 0.46
            The fraction of conductive and diffuse solar heat gain that is
            transferred by radiation to the interior thermal mass of the space.
            (see ASHRAE Fundamentals 2017, chapter 18, table 14).
        ext_shading_dev: default None
            See class `ExternalShadingDevice`. E.g. overhang or recessed window.
            In case a window is equipped with an external shading device, or in
            case of a recessed window, part of the window may be shaded depending
            on the position of the sun during the course of day.
        int_shading_dev: default None
            See class `InteriorShadingDevice`. E.g. louvered shades, roller
            shades, draperies, insect screens.
        """
        window = cls()
        window.ID = ID
        window._exterior_surface = ExteriorSurface(
            azimuth=azimuth,
            tilt=tilt,
            width=width,
            height=height,
            climate_data=climate_data,
            surface_resistance=None,
            surface_absorptance=None
        )
        window.therm_props = therm_props
        window.F_rad = F_rad
        window._ext_shading_dev = ext_shading_dev
        window._int_shading_dev = int_shading_dev
        return window

    @property
    def area(self) -> Quantity:
        return self._exterior_surface.area

    @property
    def UA(self) -> float:
        U = self.therm_props.U.to('W / (m ** 2 * K)').m
        A = self.area.to('m ** 2').m
        return U * A

    def T_ext(self, t: float) -> float:
        return self._exterior_surface.climate_data.Tdb_ext(t)

    def get_conductive_heat_gain(
        self,
        t: float,
        T_int: float
    ) -> dict[str, float]:
        """
        Get conductive heat gain through window (positive from outdoors to
        indoors) at time t seconds from 00:00:00.

        Returns a dict. The value of key `rad` is the radiative fraction of the
        conductive heat gain in Watts. The value of key `conv` is the convective
        fraction of the conductive heat gain in Watts.
        """
        U = self.therm_props.U.to('W / (m ** 2 * K)').m
        A = self.area.to('m ** 2').m
        T_ext = self._exterior_surface.climate_data.Tdb_ext(t)
        Q = U * A * (T_ext - T_int)
        F_rad = self.F_rad
        Q_rad = F_rad * Q
        Q_conv = Q - Q_rad
        return {'rad': Q_rad, 'conv': Q_conv, 'T_ext': T_ext}

    def _sunlit_area(self, t: float) -> float:
        wnd_width = self._exterior_surface.width.to('m').m
        wnd_height = self._exterior_surface.height.to('m').m
        sunlit_area = wnd_width * wnd_height
        if self._ext_shading_dev is not None:
            beta = self._exterior_surface.sun_altitude(t)
            gamma = self._exterior_surface.sun_surface_azimuth(t)
            omega = np.arctan(np.tan(beta) / np.cos(gamma))
            vp = self._ext_shading_dev.vertical_projection.to('m').m
            hp = self._ext_shading_dev.horizontal_projection.to('m').m
            wo = self._ext_shading_dev.width_offset.to('m').m
            ho = self._ext_shading_dev.height_offset.to('m').m
            shadow_width = vp * abs(np.tan(gamma))
            shadow_height = hp * np.tan(omega)
            sunlit_width = wnd_width + wo - shadow_width
            sunlit_height = wnd_height + ho - shadow_height
            if sunlit_width < 0.0:
                sunlit_width = 0.0
            elif sunlit_width > wnd_width:
                sunlit_width = wnd_width
            if sunlit_height < 0.0:
                sunlit_height = 0.0
            elif sunlit_height > wnd_height:
                sunlit_height = wnd_height
            sunlit_area = sunlit_width * sunlit_height
        return sunlit_area

    def _IAC_dir(self, t: float) -> float:
        if self._int_shading_dev.IAC_60 is not None:
            # louvred shade (slat-type sunshade)
            IAC_x = self._int_shading_dev.IAC_60 - self._int_shading_dev.IAC_0
            beta = self._exterior_surface.sun_altitude(t)
            gamma = self._exterior_surface.sun_surface_azimuth(t)
            omega = np.arctan(np.tan(beta) / np.cos(gamma))
            if self._int_shading_dev.louver_orient == 'horizontal':
                # horizontal louvred shade
                IAC_dir = self._int_shading_dev.IAC_0 + IAC_x * min(1.0, 0.02 * omega)
            else:
                # vertical louvred shade
                IAC_dir = self._int_shading_dev.IAC_0 + IAC_x * min(1.0, 0.02 * gamma)
            return IAC_dir
        else:
            return self._int_shading_dev.IAC_dif

    def get_solar_heat_gain(self, t: float) -> dict[str, float]:
        """
        Get solar heat gain through window (positive from outdoors to indoors)
        at time t seconds from 00:00:00.

        Returns a dict. The value of key `rad` is the radiative fraction of the
        conductive heat gain in Watts. The value of key `conv` is the convective
        fraction of the conductive heat gain in Watts.
        """
        theta_i = self._exterior_surface.theta_i(t) * 180.0 / np.pi
        SHGC_dir = self.therm_props.SHGC.cog_dir(theta_i)
        SHGC_dif = self.therm_props.SHGC.cog_dif
        SHGC_wnd = self.therm_props.SHGC.wnd
        f = SHGC_wnd / self.therm_props.SHGC.cog_dir(0.0)
        SHGC_dir = f * SHGC_dir
        SHGC_dif = f * SHGC_dif
        I_dir = self._exterior_surface.I_dir(t)
        I_dif = self._exterior_surface.I_dif(t)
        I_glo = self._exterior_surface.I_glo(t)
        A = self._exterior_surface.area.to('m ** 2').m
        A_sunlit = self._sunlit_area(t)

        if self._int_shading_dev is not None:
            IAC_dir = self._IAC_dir(t)
            IAC_dif = self._int_shading_dev.IAC_dif
            F_rad = self._int_shading_dev.F_rad
        else:
            IAC_dir = 1.0
            IAC_dif = 1.0
            F_rad = self.F_rad

        Q_dir = max(SHGC_dir * A_sunlit * I_dir * IAC_dir, 0.0)
        Q_dif = max(SHGC_dif * A * I_dif * IAC_dif, 0.0)
        Q_rad = F_rad * (Q_dir + Q_dif)
        Q_conv = (1.0 - F_rad) * (Q_dir + Q_dif)
        return {'rad': Q_rad, 'conv': Q_conv, 'I_glo': I_glo}

    def get_heat_transfer(
        self,
        T_int: Quantity,
        units: dict[str, str] | None = None
    ) -> pd.DataFrame:
        """
        Get the conductive and solar heat transfers through the window at each
        hour of the design day.

        Parameters
        ----------
        T_int
            Space air temperature.
        units:
            Dictionary to specify the desired units of the returned quantities.
            Use key 'T' to set the desired units for temperature (default
            unit is degC). Key 'I' is for irradiance (default unit is W/m²). And
            key 'Q' is for heat flow (default unit is W).
        """
        _units = {
            'T': 'degC',
            'I': 'W / m ** 2',
            'Q': 'W'
        }
        if units is not None:
            _units.update(units)
        T_int = T_int.to('degC').m
        time_range = range(24)
        Q_wnd_gain = []
        for t in time_range:
            Q_cond = self.get_conductive_heat_gain(t * 3600.0, T_int)
            Q_sol = self.get_solar_heat_gain(t * 3600.0)
            Q_wnd_gain.append([
                Q_(Q_cond['T_ext'], 'degC').to(_units['T']).m,
                Q_(Q_cond['conv'], 'W').to(_units['Q']).m,
                Q_(Q_cond['rad'], 'W').to(_units['Q']).m,
                Q_(Q_cond['conv'] + Q_cond['rad'], 'W').to(_units['Q']).m,
                Q_(Q_sol['I_glo'], 'W / m ** 2').to(_units['I']).m,
                Q_(Q_sol['conv'], 'W').to(_units['Q']).m,
                Q_(Q_sol['rad'], 'W').to(_units['Q']).m,
                Q_(Q_sol['conv'] + Q_sol['rad'], 'W').to(_units['Q']).m
            ])
        time_range = list(time_range)
        time_range = [time_from_decimal_hour(t) for t in time_range]
        Q_wnd_gain = pd.DataFrame(
            data=Q_wnd_gain,
            columns=pd.MultiIndex.from_arrays([
                ['Q_cond'] * 4 + ['Q_sol'] * 4,
                ['T_ext', 'conv', 'rad', 'tot', 'I_glo', 'conv', 'rad', 'tot']
            ]),
            index=time_range
        )
        return Q_wnd_gain
