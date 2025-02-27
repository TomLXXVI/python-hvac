import numpy as np
import typing
from dataclasses import dataclass
import shelve
from hvac import Quantity
from .weather_data import ExteriorSurface, WeatherData

if typing.TYPE_CHECKING:
    from .building_elements import ExteriorBuildingElement

Q_ = Quantity


class WindowThermalProperties:
    db_path: str

    class SolarHeatGainCoefficient:
        """The solar heat gain coefficient (SHGC) takes the solar radiation into
        account that is transmitted directly through the window and also the
        solar radiation that is first absorbed by the window and then released
        by convection to the interior.
        """
        def __init__(
            self,
            cog_dir: dict[float, float],
            cog_dif: float,
            wnd: float
        ) -> None:
            """Creates a `SolarHeatGainCoefficient` object.

            Parameters
            ----------
            cog_dir:
                Dictionary of which the keys are the incidence angles of direct
                sunlight in degrees and the values are the corresponding SHGC
                values of the window glazing.
            cog_dif:
                The SHGC that applies to incident diffuse solar radiation.
            wnd:
                The SHGC of the entire window at normal incidence.
            
            References
            ----------
            ASHRAE 2017, Chapter 15, Table 10.
            """
            self._cog_dir = (list(cog_dir.keys()), list(cog_dir.values()))
            self.cog_dif = cog_dif
            self.wnd = wnd

        def cog_dir(self, theta_i: float) -> float:
            """Returns the SHGC value for any incidence angle `theta_i` in
            degrees determined by linear interpolation.
            """
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
        """Creates a `WindowThermalProperties` object.

        Parameters
        ----------
        ID:
            Identifier for the type of window.
        U:
            Overall U value of the entire window, including edge effects and
            window frame (see ASHRAE Fundamentals 2017, Ch. 15, Table 4).
        SHGC_cog_dir:
            Dictionary of which the keys are solar incidence angles in degrees
            and the values the corresponding SHGC values (see ASHRAE
            Fundamentals 2017, Ch. 15, Table 10).
        SHGC_cog_dif:
            The SHGC for diffuse solar radiation (see ASHRAE Fundamentals 2017,
            Ch. 15, Table 10).
        SHGC_wnd:
            The SHGC of the entire window at normal incidence (see ASHRAE
            Fundamentals 2017, Ch. 15, Table 10).
        """
        self.ID = ID
        self.U = U
        self.SHGC = self.SolarHeatGainCoefficient(
            SHGC_cog_dir,
            SHGC_cog_dif,
            SHGC_wnd
        )

    def save(self):
        """Saves the `WindowThermalProperties` object to a shelf on disk.
        The file path to the shelf must be set first through class attribute
        `db_path` (str).
        """
        with shelve.open(self.db_path) as shelf:
            shelf[self.ID] = self

    @classmethod
    def load(cls, ID: str) -> 'WindowThermalProperties':
        """Loads the `WindowThermalProperties` object with the given ID
        from the shelf, whose file path was set through class attribute
        `db_path` (str).
        """
        with shelve.open(cls.db_path) as shelf:
            return typing.cast('WindowThermalProperties', shelf[ID])

    @classmethod
    def overview(cls) -> list[str]:
        """Returns a list with all the ID's of objects stored in the shelf,
        whose file path was set through class attribute `db_path` (str).
        """
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
    """Represents an overhang above a window or a recessed window.

    Parameters
    ----------
    vert_proj_dist: default 0 m
        Protrusion distance of the vertical surfaces of the exterior shading 
        device.
    hor_proj_dist: default 0 m
        Protrusion distance of the exterior shading device's horizontal surface 
        above the window.
    hor_offset : default 0 m
        Distance between the upright window edges and the vertical sides of the 
        exterior shading device.
    vert_offset : default 0 m
        Distance between the window's horizontal upper edge and the underside of
        the exterior shading device's horizontal surface.
    """
    vert_proj_dist: Quantity = Q_(0, 'm')
    hor_proj_dist: Quantity = Q_(0, 'm')
    hor_offset: Quantity = Q_(0, 'm')
    vert_offset: Quantity = Q_(0, 'm')


@dataclass
class InteriorShadingDevice:
    """Represents louvered shades, roller shades, draperies or insect screens
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
        Interior attenuation coefficient for direct solar radiation at an
        incidence angle of 60°.
    louver_orient: default 'horizontal'
        Sets the orientation of the louvers: either 'horizontal' or 'vertical'.

    Notes
    -----
    Only for louvered shades, three IAC values are being considered (IAC_0,
    IAC_60, and IAC_dif). In other cases there is only one IAC value, which
    applies both to direct and diffuse solar radiation. In these cases it
    suffices to provide only a value to parameter `IAC_dif`.
    """
    IAC_dif: float
    F_rad: float
    IAC_0: float | None = None
    IAC_60: float | None = None
    louver_orient: str = 'horizontal'


class Window:
    """Represents a window in an exterior building element."""

    def __init__(self):
        self.name: str = ''
        self.props: WindowThermalProperties | None = None
        self.width: Quantity | None = None
        self.height: Quantity | None = None
        self.area: Quantity | None = None
        self.F_rad: float | None = None
        self.ext_surf: ExteriorSurface | None = None
        self.ext_shading: ExteriorShadingDevice | None = None
        self.int_shading: InteriorShadingDevice | None = None
        self.parent: ExteriorBuildingElement | None = None
        
    @classmethod
    def create(
        cls,
        name: str,
        azimuth_angle: Quantity,
        slope_angle: Quantity,
        width: Quantity,
        height: Quantity,
        weather_data: WeatherData,
        props: WindowThermalProperties,
        F_rad: float = 0.46,
        ext_shading: ExteriorShadingDevice | None = None,
        int_shading: InteriorShadingDevice | None = None
    ) -> 'Window':
        """Creates a `Window` object.

        Parameters
        ----------
        name:
            Name to identify the `Window` object.
        azimuth_angle:
            The azimuth angle of the window. South = 0°, West = +90° and 
            East = -90°.
        slope_angle:
            Slope angle of the window. E.g., a vertical window has a slope angle
            of 90°.
        width:
            Width of the window opening.
        height:
            Height of the window opening.
        weather_data:
            Instance of `WeatherData` class, containing the climatic design
            information.
        props:
            See class `WindowThermalProperties`. The thermal and solar properties
            of the window.
        F_rad: default 0.46
            Radiative fraction of conduction heat gain through window to the 
            interior thermal mass of the space. (see ASHRAE Fundamentals 2017, 
            chapter 18, table 14).
        ext_shading: default None
            E.g. an overhang or recessed window. See: class `ExternalShadingDevice`.
            If a window has an external shading device, or if the window is
            recessed, part of the window may be shaded depending on the position
            of the sun during the day.
        int_shading: default None
            E.g. a louvered shade, roller shade, drapery, or insect screen.
            See: class `InteriorShadingDevice`.
        """
        window = cls()
        window.name = name
        window.ext_surf = ExteriorSurface(
            weather_data=weather_data,
            gamma=azimuth_angle,
            beta=slope_angle,
        )
        window.props = props
        window.width = width
        window.height = height
        window.area = width * height
        window.F_rad = F_rad
        window.ext_shading = ext_shading
        window.int_shading = int_shading
        return window
    
    def solar_heat_gain(self, t: Quantity) -> tuple[Quantity, ...]:
        """Returns the solar heat gain at solar time `t`.
        
        Returns
        -------
        Q_sol:
            The global solar heat gain through the window at solar time `t`.
        Q_sol_conv:
            The convective fraction of the solar heat gain.
        Q_sol_rad:
            The radiative fraction of the solar heat gain.
        """
        t = t.to('s').m
        theta_i = self.ext_surf.theta_i(t)
        SHGC_dir = self.props.SHGC.cog_dir(theta_i.to('deg').m)
        SHGC_dif = self.props.SHGC.cog_dif
        SHGC_wnd = self.props.SHGC.wnd
        f = SHGC_wnd / self.props.SHGC.cog_dir(0.0)
        SHGC_dir = f * SHGC_dir
        SHGC_dif = f * SHGC_dif
        G_T, G_Tb, G_Td = self.ext_surf.G_T(t)
        A = self.area.to('m**2')
        A_sunlit = self._sunlit_area(t)
        if self.int_shading is not None:
            IAC_dir = self._IAC_dir(t)
            IAC_dif = self.int_shading.IAC_dif
            F_rad = self.int_shading.F_rad
        else:
            IAC_dir = 1.0
            IAC_dif = 1.0
            F_rad = 1.0
        Q_dot_dir = max(SHGC_dir * A_sunlit * G_Tb * IAC_dir, Q_(0.0, 'W'))
        Q_dot_dif = max(SHGC_dif * A * G_Td * IAC_dif, Q_(0.0, 'W'))
        Q_dot_sol = Q_dot_dir + Q_dot_dif
        Q_dot_rad = F_rad * Q_dot_sol
        Q_dot_conv = (1.0 - F_rad) * Q_dot_sol
        return Q_dot_sol, Q_dot_conv, Q_dot_rad
    
    def _sunlit_area(self, t: float) -> Quantity:
        A_sunlit = self.area
        if self.ext_shading is not None:
            w_wnd = self.width.to('m').m
            h_wnd = self.height.to('m').m
            beta = self.ext_surf.alpha_s(t).to('rad').m
            gamma = self.ext_surf.gamma_ss(t).to('rad').m
            omega = np.arctan(np.tan(beta) / np.cos(gamma))
            vp = self.ext_shading.vert_proj_dist.to('m').m
            hp = self.ext_shading.hor_proj_dist.to('m').m
            wo = self.ext_shading.hor_offset.to('m').m
            ho = self.ext_shading.vert_offset.to('m').m
            w_shadow = vp * abs(np.tan(gamma))
            h_shadow = hp * np.tan(omega)
            w_sunlit = w_wnd + wo - w_shadow
            h_sunlit = h_wnd + ho - h_shadow
            if w_sunlit < 0.0:
                w_sunlit = 0.0
            elif w_sunlit > w_wnd:
                w_sunlit = w_wnd
            if h_sunlit < 0.0:
                h_sunlit = 0.0
            elif h_sunlit > h_wnd:
                h_sunlit = h_wnd
            A_sunlit = w_sunlit * h_sunlit
        return Q_(A_sunlit, 'm**2')

    def _IAC_dir(self, t: float) -> float:
        if self.int_shading.IAC_60 is not None:
            # louvred shade (slat-type sunshade)
            IAC_x = self.int_shading.IAC_60 - self.int_shading.IAC_0
            beta = self.ext_surf.alpha_s(t).to('rad').m  # solar altitude angle
            gamma = self.ext_surf.gamma_ss(t).to('rad').m  # surface-solar azimuth angle
            if self.int_shading.louver_orient == 'horizontal':
                # horizontal louvred shade
                omega = np.arctan(np.tan(beta) / np.cos(gamma))  # vertical profile angle
                IAC_dir = self.int_shading.IAC_0 + IAC_x * min(1.0, 0.02 * omega)
            else:
                # vertical louvred shade
                IAC_dir = self.int_shading.IAC_0 + IAC_x * min(1.0, 0.02 * gamma)
            return IAC_dir
        else:
            return self.int_shading.IAC_dif

    @property
    def UA(self) -> Quantity:
        """Returns the total transmittance of the window."""
        U = self.props.U.to('W / (m**2 * K)')
        A = self.area.to('m**2')
        return U * A
    
    def conduction_heat_gain(self, T_za: Quantity) -> tuple[Quantity, ...]:
        """Returns the hourly values of conductive heat gain through the window
        on the selected design day, together with its convective and radiative 
        components.

        Parameters
        ----------
        T_za:
            Zone-air temperature.

        Returns
        -------
        Q_dot_cond:
            `Quantity` array containing the hourly average values of conduction 
            heat gain on the selected design day and with the specified zone-air
            temperature.
        Q_dot_cond_conv:
            Part of the conduction heat gain which is transferred by convection 
            to the zone air.
        Q_dot_cond_rad:
            Part of the conduction heat gain which is transferred by radiation 
            to the interior thermal mass of the zone.
        """
        T_db = Quantity.from_list([
            self.ext_surf.T_db(hr * 3600) for hr in range(0, 24)
        ])
        T_za = Quantity.from_list([T_za] * len(T_db))
        Q_dot_cnd = self.UA * (T_db - T_za)
        Q_dot_cnd_rad = self.F_rad * Q_dot_cnd
        Q_dot_cnd_conv = Q_dot_cnd - Q_dot_cnd_rad
        return Q_dot_cnd, Q_dot_cnd_conv, Q_dot_cnd_rad
        