from typing import cast, Callable
from dataclasses import dataclass
import shelve
import numpy as np
from hvac import Quantity
from hvac.cooling_load_calc_old.core.weather_data import WeatherData, ExteriorSurface


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
        The file path to the shelf must be set through class attribute
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
            return cast('WindowThermalProperties', shelf[ID])

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
        Projection distance of the vertical part of the exterior shading device.
    hor_proj_dist: default 0 m
        Projection distance of the horizontal part of exterior shading device.
    hor_offset : default 0 m
        Distance between the vertical window edge and the vertical side of the
        exterior shading device.
    vert_offset : default 0 m
        Distance between the window's upper edge and the horizontal underside of
        exterior shading device.
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
        self.ID: str = ''
        self.T_zone: Callable[[float], Quantity] | None = None
        self.props: WindowThermalProperties | None = None
        self.width: Quantity | None = None
        self.height: Quantity | None = None
        self.area: Quantity | None = None
        self.F_rad: float | None = None
        self._ext_surf: ExteriorSurface | None = None
        self._ext_shading: ExteriorShadingDevice | None = None
        self._int_shading: InteriorShadingDevice | None = None

    @classmethod
    def create(
        cls,
        ID: str,
        T_zone: Callable[[float], Quantity],
        gamma: Quantity,
        beta: Quantity,
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
        ID:
            Useful name to identify the `Window` object.
        T_zone:
            The zone air temperature, being a function with signature
            `f(t_sol_sec: float) -> Quantity` which takes the solar time in
            seconds and returns the temperature in the zone as a `Quantity`
            object. This may allow a time-variable temperature in the zone.
        gamma:
            Azimuth angle of the window.
        beta:
            Tilt angle of the window.
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
            Radiative fraction of solar heat gain through window to the interior
            thermal mass of the space. (see ASHRAE Fundamentals 2017, chapter 18,
            table 14).
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
        window.ID = ID
        window.T_zone = T_zone
        window._ext_surf = ExteriorSurface(
            weather_data=weather_data,
            gamma=gamma,
            beta=beta,
        )
        window.props = props
        window.width = width
        window.height = height
        window.area = width * height
        window.F_rad = F_rad
        window._ext_shading = ext_shading
        window._int_shading = int_shading
        return window

    def conductive_heat_gain(
        self,
        dt_hr: float = 1.0,
    ) -> tuple[Quantity, Quantity, Quantity]:
        """Returns the conductive heat output at the interior side of the
        window, and also its convective and radiative component, for each time
        index k of the design day.

        Parameters
        ----------
        dt_hr:
            Time step expressed as a fraction of 1 hour, e.g., `dt_hr` = 1/4
            means the time step for the calculations is one quarter of an hour.
            The default value is 1 hour.

        Returns
        -------
        3-tuple with:
        - a `Quantity`-array with the conductive heat gains at each time index
        - a `Quantity`-array with the convective components
        - a `Quantity`-array with the radiative components
        """
        num_steps = int(round(24 / dt_hr))
        dt_sec = dt_hr * 3600
        U = self.props.U.to('W / (m**2 * K)')
        A = self.area.to('m**2')
        dT = Quantity.from_list([
            self._ext_surf.T_db(k * dt_sec).to('K') - self.T_zone(k * dt_sec).to('K')
            for k in range(num_steps)
        ])
        Q_dot_cond = U * A * dT
        Q_dot_rad = self.F_rad * Q_dot_cond
        Q_dot_conv = (1.0 - self.F_rad) * Q_dot_cond
        return Q_dot_cond, Q_dot_conv, Q_dot_rad

    def solar_heat_gain(
        self,
        dt_hr: float = 1.0
    ) -> tuple[Quantity, Quantity, Quantity]:
        """Returns the solar heat gain through the window, and also its
        radiative and convective components, for each time index k of the design
        day.

        Parameters
        ----------
        dt_hr:
            Time step expressed as a fraction of 1 hour, e.g., `dt_hr` = 1/4
            means the time step between the calculations is one quarter of an
            hour. The default value is 1 hour.

        Returns
        -------
        Returns
        -------
        3-tuple with:
        - a `Quantity`-array with the solar heat gains at each time index
        - a `Quantity`-array with the convective components
        - a `Quantity`-array with the radiative components
        """
        num_steps = int(round(24 / dt_hr))
        dt_sec = dt_hr * 3600
        Q_dot_sol, Q_dot_rad, Q_dot_conv = [], [], []
        for k in range(num_steps):
            theta_i = self._ext_surf.theta_i(k * dt_sec)
            SHGC_dir = self.props.SHGC.cog_dir(theta_i.to('deg').m)
            SHGC_dif = self.props.SHGC.cog_dif
            SHGC_wnd = self.props.SHGC.wnd
            f = SHGC_wnd / self.props.SHGC.cog_dir(0.0)
            SHGC_dir = f * SHGC_dir
            SHGC_dif = f * SHGC_dif
            G_T, G_Tb, G_Td = self._ext_surf.G_T(k * dt_sec)
            A = self.area.to('m**2')
            A_sunlit = self._sunlit_area(k * dt_sec)
            if self._int_shading is not None:
                IAC_dir = self._IAC_dir(k * dt_sec)
                IAC_dif = self._int_shading.IAC_dif
                F_rad = self._int_shading.F_rad
            else:
                IAC_dir = 1.0
                IAC_dif = 1.0
                F_rad = self.F_rad
            Q_dot_dir = max(SHGC_dir * A_sunlit * G_Tb * IAC_dir, Q_(0.0, 'W'))
            Q_dot_dif = max(SHGC_dif * A * G_Td * IAC_dif, Q_(0.0, 'W'))
            Q_dot_sol_ = Q_dot_dir + Q_dot_dif
            Q_dot_rad_ = F_rad * Q_dot_sol_
            Q_dot_conv_ = (1.0 - F_rad) * Q_dot_sol_
            Q_dot_sol.append(Q_dot_sol_.to('W'))
            Q_dot_rad.append(Q_dot_rad_.to('W'))
            Q_dot_conv.append(Q_dot_conv_.to('W'))
        Q_dot_sol = Quantity.from_list(Q_dot_sol)
        Q_dot_rad = Quantity.from_list(Q_dot_rad)
        Q_dot_conv = Quantity.from_list(Q_dot_conv)
        return Q_dot_sol, Q_dot_conv, Q_dot_rad

    def _sunlit_area(self, t_sol_sec: float) -> Quantity:
        A_sunlit = self.area
        if self._ext_shading is not None:
            w_wnd = self.width.to('m').m
            h_wnd = self.height.to('m').m
            beta = self._ext_surf.alpha_s(t_sol_sec).to('rad').m
            gamma = self._ext_surf.gamma_ss(t_sol_sec).to('rad').m
            omega = np.arctan(np.tan(beta) / np.cos(gamma))
            vp = self._ext_shading.vert_proj_dist.to('m').m
            hp = self._ext_shading.hor_proj_dist.to('m').m
            wo = self._ext_shading.hor_offset.to('m').m
            ho = self._ext_shading.vert_offset.to('m').m
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

    def _IAC_dir(self, t_sol_sec: float) -> float:
        if self._int_shading.IAC_60 is not None:
            # louvred shade (slat-type sunshade)
            IAC_x = self._int_shading.IAC_60 - self._int_shading.IAC_0
            beta = self._ext_surf.alpha_s(t_sol_sec).to('rad').m  # solar altitude angle
            gamma = self._ext_surf.gamma_ss(t_sol_sec).to('rad').m  # surface-solar azimuth angle
            if self._int_shading.louver_orient == 'horizontal':
                # horizontal louvred shade
                omega = np.arctan(np.tan(beta) / np.cos(gamma))  # vertical profile angle
                IAC_dir = self._int_shading.IAC_0 + IAC_x * min(1.0, 0.02 * omega)
            else:
                # vertical louvred shade
                IAC_dir = self._int_shading.IAC_0 + IAC_x * min(1.0, 0.02 * gamma)
            return IAC_dir
        else:
            return self._int_shading.IAC_dif

    @property
    def UA(self) -> Quantity:
        """Returns the total transmittance of the window."""
        U = self.props.U.to('W / (m**2 * K)')
        A = self.area.to('m**2')
        return U * A
