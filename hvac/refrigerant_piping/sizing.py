"""
SIZING REFRIGERATION LINES IN AIR CONDITIONING SYSTEMS.

This module implements the sizing procedures for suction lines, discharge lines,
and liquid lines made from copper piping as outlined in TRANE, Air Conditioning
Clinic: Refrigerant Piping (June 2011, TRG-TRC006-EN).
"""
from collections.abc import Iterable
import math
from abc import ABC, abstractmethod
from dataclasses import dataclass
from scipy.interpolate import interp1d
from hvac import Quantity
from hvac.fluids import Fluid
from hvac.heat_transfer.forced_convection.internal_flow import CircularTube
from hvac.refrigerant_piping.copper_tubing import CopperTube


Q_ = Quantity
g = Q_(9.81, 'm / s**2')


def _get_suction_u_min_limit_interpolation():
    """Returns an interpolation function that takes the nominal tube diameter
    in inch and returns the minimum allowable refrigerant flow velocity in
    feet / min in the suction line.
    """
    u_min_limit = [
        (3/8, 370),
        (1/2, 460),
        (5/8, 520),
        (3/4, 560),
        (7/8, 600),
        (1 + 1/8, 700),
        (1 + 3/8, 780),
        (1 + 5/8, 840),
        (2 + 1/8, 980),
        (2 + 5/8, 1080),
        (3 + 1/8, 1180),
        (3 + 5/8, 1270),
        (4 + 1/8, 1360)
    ]
    x, y = zip(*u_min_limit)
    interp = interp1d(x, y, bounds_error=False, fill_value='extrapolate')
    return interp


def _get_discharge_u_min_limit_interpolation():
    """Returns an interpolation function that takes the nominal tube diameter
    in inch and returns the minimum allowable refrigerant flow velocity in
    feet / min in the discharge line.
    """
    u_min_limit = [
        (5/16, 220),
        (3/8, 250),
        (1/2, 285),
        (5/8, 315),
        (3/4, 345),
        (7/8, 375),
        (1 + 1/8, 430),
        (1 + 3/8, 480),
        (1 + 5/8, 520),
        (2 + 1/8, 600),
        (2 + 5/8, 665),
        (3 + 1/8, 730)
    ]
    x, y = zip(*u_min_limit)
    interp = interp1d(x, y, bounds_error=False, fill_value='extrapolate')
    return interp


interp_suction_u_min_limit = _get_suction_u_min_limit_interpolation()
interp_discharge_u_min_limit = _get_discharge_u_min_limit_interpolation()


@dataclass
class RefrigerantCycleInfo:
    """Dataclass that holds information about the refrigeration cycle.

    Parameters
    ----------
    refrigerant:
        Type of refrigerant fluid.
    T_evp:
        Evaporation temperature.
    T_cnd:
        Condensing temperature.
    evp_superheat:
        Amount of refrigerant superheat at the evaporator outlet with respect
        to the evaporation temperature.
    cnd_subcooling:
        Amount of refrigerant subcooling at the condenser outlet with respect
        to the condensing temperature.
    cmp_superheat:
        Amount of refrigerant superheat at the condenser inlet with respect
        to the condensing temperature.
    """
    refrigerant: Fluid
    T_evp: Quantity
    T_cnd: Quantity
    evp_superheat: Quantity
    cnd_subcooling: Quantity
    cmp_superheat: Quantity

    def __post_init__(self):
        self.T_evp_out = self.T_evp.to('K') + self.evp_superheat.to('K')
        self.T_cnd_in = self.T_cnd.to('K') + self.cmp_superheat.to('K')
        self.T_cnd_out = self.T_cnd.to('K') - self.cnd_subcooling.to('K')
        self.P_evp = self.refrigerant(T=self.T_evp, x=Q_(0, 'frac')).P
        self.P_cnd = self.refrigerant(T=self.T_cnd, x=Q_(0, 'frac')).P
        self.rfg_cnd_sat_liq = self.refrigerant(T=self.T_cnd, x=Q_(0, 'frac'))
        self.rfg_cnd_sat_vap = self.refrigerant(T=self.T_cnd, x=Q_(1, 'frac'))
        self.rfg_cnd_in = self.refrigerant(P=self.P_cnd, T=self.T_cnd_in)
        self.rfg_cnd_out = self.refrigerant(P=self.P_cnd, T=self.T_cnd_out)
        self.rfg_evp_sat_liq = self.refrigerant(T=self.T_evp, x=Q_(0.0, 'frac'))
        self.rfg_evp_sat_vap = self.refrigerant(T=self.T_evp, x=Q_(1.0, 'frac'))
        self.rfg_evp_in = self.refrigerant(P=self.P_evp, h=self.rfg_cnd_out.h)
        self.rfg_evp_out = self.refrigerant(P=self.P_evp, T=self.T_evp_out)
        self.q_evp = self.rfg_evp_out.h - self.rfg_evp_in.h


class RefrigerantLine(ABC):
    u_max_limit: Quantity = None

    def __init__(
        self,
        cycle_info: RefrigerantCycleInfo,
        copper_tube: CopperTube,
        L: Quantity,
        Q_dot_evp_max: Quantity,
        elevation: Quantity = Q_(0.0, 'm'),
        Q_dot_evp_min: Quantity | None = None
    ) -> None:
        """
        Parameters
        ----------
        cycle_info:
            Instance of dataclass `RefrigerantCycleInfo` containing the
            operating design parameters of the refrigerant cycle.
        copper_tube:
            Instance of dataclass `CopperTube` containing the dimensions of the
            circular cross-section (nominal diameter, outside diameter, inside
            diameter).
        L:
            Straight length of the refrigerant line without fittings and other
            accessories.
        Q_dot_evp_max:
            Maximum evaporator load. This will determine the maximum
            flow velocity of the refrigerant in the tube.
        elevation:
            The difference in height between the outlet and the inlet of the
            refrigerant line. If the outlet is below the inlet, the elevation
            has a negative value.
        Q_dot_evp_min: optional, default None
            Minimum evaporator load. This will determine the minimum flow
            velocity of the refrigerant. Leave this `None` in case of an
            on/off-controlled compressor.
        """
        self.copper_tube = copper_tube
        self.L = L
        self.elevation = elevation
        self.cycle_info = cycle_info
        self.Q_dot_evp_max = Q_dot_evp_max
        self.Q_dot_evp_min = Q_dot_evp_min or Q_dot_evp_max

        self.m_dot_max = self._get_m_dot(self.Q_dot_evp_max)
        self.V_dot_max = self._get_V_dot(self.m_dot_max)
        self.m_dot_min = self._get_m_dot(self.Q_dot_evp_min)
        self.V_dot_min = self._get_V_dot(self.m_dot_min)

        self._tube: CircularTube | None = None
        self.fittings: Iterable[tuple[str, Quantity, int]] | None = None
        self.u_max: Quantity = Q_(float('nan'), 'm / s')
        self.u_min: Quantity = Q_(float('nan'), 'm / s')
        self.L_eq: Quantity = Q_(float('nan'), 'm')
        self.dP: Quantity = Q_(float('nan'), 'Pa')
        self.dT_sat: Quantity = Q_(float('nan'), 'K')
        self.is_double_riser: bool = False

        self.units = {'u': 'm / s', 'dP': 'kPa', 'dT_sat': 'K'}
        self.warning: str | None = None

    def _get_m_dot(self, Q_dot: Quantity) -> Quantity:
        """Returns the steady-state refrigerant mass flow rate based on the
        system's cooling capacity `Q_dot`.
        """
        return Q_dot / self.cycle_info.q_evp

    @abstractmethod
    def _get_V_dot(self, m_dot: Quantity) -> Quantity:
        """Returns the refrigerant volume flow rate that corresponds with the
        given mass flow rate `m_dot`.
        """
        ...

    @abstractmethod
    def flow_velocity(self) -> tuple:
        """Get refrigerant flow velocity."""
        ...

    def _check_u_max(self, u_max: Quantity) -> str:
        if u_max.to('feet / min') < self.u_max_limit:
            r = 'OK'
        else:
            r = 'TOO HIGH'
        return r

    def _check_u_min(self, u_min: Quantity) -> str:
        ...

    def _calculate_equivalent_length(
        self,
        fittings: Iterable[tuple[str, Quantity, int]] | None
    ) -> Quantity:
        """Calculates the equivalent length of the refrigeration line together
        with the fittings present in the line.
        """
        fittings_L_eq = Q_(0.0, 'm')
        if isinstance(fittings, Iterable):
            self.fittings = fittings
            fittings_L_eq = sum(
                n * fitting_L_eq
                for _, fitting_L_eq, n in fittings
            )
        self.L_eq = self.L + fittings_L_eq
        return self.L_eq

    @abstractmethod
    def pressure_drop(
        self,
        fittings: Iterable[tuple[str, Quantity, int]] | None = None,
        dP_add: Quantity | Iterable[Quantity] | None = None
    ) -> tuple[Quantity, Quantity]:
        """Returns the pressure loss across the refrigeration line at the
        maximum cooling load due to flow friction and elevation of the outlet
        with respect to the inlet of the line. However, in the case of a
        negative elevation, the associated pressure gain is ignored. Also,
        returns the associated change in saturation temperature that corresponds
        with the pressure drop.
        """
        ...

    def get_equivalent_length(self) -> Quantity | None:
        """Returns the total equivalent length of the line including all its
        fittings and other accessories. Note that before calling this method,
        the pressure drop along the line must have been calculated by calling
        method `pressure_drop(...)`.
        """
        if not math.isnan(self.dP.magnitude):
            f = self._tube.friction_factor()
            u = self._tube.mean_velocity()
            D = self._tube.Dh
            rho = self._tube.fluid.rho
            L_eq = 2 * D * self.dP / (f * rho * u ** 2)
            return L_eq.to('m')
        return None

    def __str__(self):
        u_u = self.units['u']
        u_dP = self.units['dP']
        u_dT_sat = self.units['dT_sat']
        fields = [
            f"{self.copper_tube.DN}",
            f"maximum flow velocity = {self.u_max.to(u_u):~P.1f}",
            f"minimum flow velocity = {self.u_min.to(u_u):~P.1f}",
            f"maximum pressure drop = {self.dP.to(u_dP):~P.3f} "
            f"({self.dT_sat.to(u_dT_sat):~P.3f})"
        ]
        if self.is_double_riser:
            fields.insert(0, "double riser")
        s = ' | '.join(fields)
        if self.warning is not None:
            s += f"\nWarning: {self.warning}"
        return s


class VaporLine(RefrigerantLine):

    @abstractmethod
    def _get_V_dot(self, m_dot: Quantity) -> Quantity:
        ...

    def flow_velocity(self) -> tuple[Quantity, str, Quantity, str]:
        """Returns the flow velocities of the refrigerant in the line that
        corresponds with the maximum and minimum cooling load of the system.
        This method also checks if maximum and minimum velocity are between the
        allowable upper and lower limit to avoid noise and to ensure proper oil
        return.

        Returns
        -------
        Tuple with 4 elements:
        1.  the flow velocity at maximum cooling load
        2.  string that indicates if the maximum flow velocity is 'OK' or
            'TOO HIGH'
        3.  the flow velocity at minimum cooling load
        4.  string that indicates if minimum flow velocity is 'OK', 'TOO LOW',
            or 'TOO HIGH'
        """
        A = math.pi * (self.copper_tube.D_int ** 2) / 4
        self.u_max = self.V_dot_max / A
        r_max = self._check_u_max(self.u_max)
        self.u_min = self.V_dot_min / A
        r_min = self._check_u_min(self.u_min)
        return self.u_max.to('m / s'), r_max, self.u_min.to('m / s'), r_min

    def _check_u_max(self, u_max: Quantity) -> str:
        if u_max.to('feet / min') < self.u_max_limit:
            r = 'OK'
        else:
            r = 'TOO HIGH'
        return r

    @abstractmethod
    def _get_u_min_limit(self, DN: str) -> Quantity:
        pass

    def _check_u_min(self, u_min: Quantity) -> str:
        u_min = u_min.to('feet / min')
        # minimum allowable flow velocity to ensure oil return:
        u_min_limit = self._get_u_min_limit(self.copper_tube.DN)
        if not self.elevation.magnitude > 0.0:
            # if not a vertical suction riser:
            u_min_limit *= 0.75
        if u_min_limit <= u_min < self.u_max_limit:
            r = 'OK'
        elif u_min < u_min_limit:
            r = 'TOO LOW'
        else:
            r = 'TOO HIGH'
        return r

    @abstractmethod
    def pressure_drop(
        self,
        fittings: Iterable[tuple[str, Quantity, int]] | None = None,
        dP_add: Quantity | Iterable[Quantity] | None = None
    ) -> tuple[Quantity, Quantity]:
        ...


class SuctionLine(VaporLine):
    u_max_limit = Q_(4000, 'ft / min')

    def _get_V_dot(self, m_dot: Quantity) -> Quantity:
        # volume flow rate of refrigerant at evaporator outlet
        return m_dot / self.cycle_info.rfg_evp_out.rho

    def _get_u_min_limit(self, DN: str) -> Quantity:
        DN = '+'.join(DN.split(' '))
        DN = eval(DN)
        u_min_limit = interp_suction_u_min_limit(DN)
        return Q_(u_min_limit, 'feet / min')

    def pressure_drop(
        self,
        fittings: Iterable[tuple[str, Quantity, int]] | None = None,
        dP_add: Quantity | Iterable[Quantity] | None = None
    ) -> tuple[Quantity, Quantity]:
        rfg = self.cycle_info.rfg_evp_out
        L_eq = self._calculate_equivalent_length(fittings)
        self._tube = CircularTube(self.copper_tube.D_int, L_eq, rfg, e=Q_(0.0015, 'mm'))
        self._tube.m_dot = self.m_dot_max
        dP_friction = self._tube.pressure_drop()
        if self.elevation.magnitude > 0.0:
            dP_elevation = g * rfg.rho * self.elevation
        else:
            dP_elevation = Q_(0.0, 'Pa')
        if dP_add is not None:
            try:
                len(dP_add)
            except TypeError:
                # `dP_add` is a single `Quantity`
                pass
            else:
                dP_add = sum(dP_add)
        else:
            dP_add = Q_(0.0, 'Pa')
        self.dP = dP_friction + dP_add
        dP = self.dP + dP_elevation
        P = self.cycle_info.P_evp - dP
        T_sat = self.cycle_info.refrigerant(P=P, x=Q_(1.0, 'frac')).T
        self.dT_sat = self.cycle_info.T_evp - T_sat
        return self.dP.to('Pa'), self.dT_sat.to('K')


class DischargeLine(VaporLine):
    u_max_limit = Q_(3500, 'ft / min')

    def _get_V_dot(self, m_dot: Quantity) -> Quantity:
        # volume flow rate of refrigerant at compressor outlet
        return m_dot / self.cycle_info.rfg_cnd_in.rho

    def _get_u_min_limit(self, DN: str) -> Quantity:
        DN = '+'.join(DN.split(' '))
        DN = eval(DN)
        u_min_limit = interp_suction_u_min_limit(DN)
        return Q_(u_min_limit, 'feet / min')

    def pressure_drop(
        self,
        fittings: Iterable[tuple[str, Quantity, int]] | None = None,
        dP_add: Quantity | Iterable[Quantity] | None = None
    ) -> tuple[Quantity, Quantity]:
        rfg = self.cycle_info.rfg_cnd_in
        L_eq = self._calculate_equivalent_length(fittings)
        self._tube = CircularTube(self.copper_tube.D_int, L_eq, rfg, e=Q_(0.0015, 'mm'))
        self._tube.m_dot = self.m_dot_max
        dP_friction = self._tube.pressure_drop()
        if self.elevation.magnitude > 0.0:
            dP_elevation = g * rfg.rho * self.elevation
        else:
            dP_elevation = Q_(0.0, 'Pa')
        if dP_add is not None:
            try:
                len(dP_add)
            except TypeError:
                # `dP_add` is a single `Quantity`
                pass
            else:
                dP_add = sum(dP_add)
        else:
            dP_add = Q_(0.0, 'Pa')
        self.dP = dP_friction + dP_add
        dP = self.dP + dP_elevation
        P = self.cycle_info.P_cnd + dP
        T_sat = self.cycle_info.refrigerant(P=P, x=Q_(1.0, 'frac')).T
        self.dT_sat = T_sat - self.cycle_info.T_cnd
        return self.dP.to('Pa'), self.dT_sat.to('K')


class LiquidLine(RefrigerantLine):
    u_max_limit = Q_(600, 'ft / min')
    # Minimum limit of remaining subcooling at expansion device:
    dT_sc_min = Q_(5, 'delta_degF')

    def _get_V_dot(self, m_dot: Quantity) -> Quantity:
        # volume flow rate of refrigerant at condenser outlet
        return m_dot / self.cycle_info.rfg_cnd_out.rho

    def flow_velocity(self) -> tuple[Quantity, str]:
        A = math.pi * (self.copper_tube.D_int ** 2) / 4
        self.u_max = self.V_dot_max / A
        r_max = self._check_u_max(self.u_max)
        return self.u_max.to('m / s'), r_max

    def pressure_drop(
        self,
        fittings: Iterable[tuple[str, Quantity, int]] = None,
        dP_add: Quantity | Iterable[Quantity] | None = None
    ) -> tuple[Quantity, Quantity]:
        rfg = self.cycle_info.rfg_cnd_out
        L_eq = self._calculate_equivalent_length(fittings)
        self._tube = CircularTube(self.copper_tube.D_int, L_eq, rfg, e=Q_(0.0015, 'mm'))
        self._tube.m_dot = self.m_dot_max
        dP_friction = self._tube.pressure_drop()
        if self.elevation.magnitude > 0.0:
            dP_elevation = g * rfg.rho * self.elevation
        else:
            dP_elevation = Q_(0.0, 'Pa')
        if dP_add is not None:
            try:
                len(dP_add)
            except TypeError:
                # `dP_add` is a single `Quantity`
                pass
            else:
                dP_add = sum(dP_add)
        else:
            dP_add = Q_(0.0, 'Pa')
        self.dP = dP_friction + dP_add
        dP = self.dP + dP_elevation
        P = self.cycle_info.P_cnd - dP
        T_sat = self.cycle_info.refrigerant(P=P, x=Q_(0.0, 'frac')).T
        self.dT_sat = self.cycle_info.T_cnd - T_sat
        if not self._check_remaining_subcooling():
            self.warning = "Flashing of refrigerant in the liquid line is likely to occur."
        return self.dP.to('Pa'), self.dT_sat.to('K')

    def _check_remaining_subcooling(self) -> bool:
        """Checks if enough subcooling of the refrigerant is still available at
        the entrance of the expansion device.
        """
        # Calculate the difference between the saturation temperature of the
        # subcooled refrigerant at the outlet of the condenser and the
        # saturation temperature of liquid refrigerant having the same enthalpy
        # as the refrigerant at the condenser outlet (i.e. where in the
        # log(p)/h-diagram the vertical line through the refrigerant state at
        # the condenser outlet intersects the liquid saturation line ). This is
        # the maximum saturation temperature drop available before flashing of
        # refrigerant will occur.
        h = self.cycle_info.rfg_cnd_out.h
        rfg_sat_liq = self.cycle_info.refrigerant(
            x=Q_(0.0, 'frac'),
            h=h,
            P=self.cycle_info.P_evp
        )
        dT_sat_max = self.cycle_info.T_cnd.to('K') - rfg_sat_liq.T.to('K')
        # The difference between the maximum available saturation temperature
        # drop and the actual saturation temperature drop must be greater than
        # the safety margin:
        diff = dT_sat_max - self.dT_sat
        if diff.to('K') < self.dT_sc_min.to('K'):
            return False
        return True


class RefrigerantLineSizer(ABC):
    """
    Abstract base class for sizing suction lines, discharge lines, and liquid
    lines.
    """
    _RFG_LINE = None

    def __init__(
        self,
        copper_tubes: Iterable[CopperTube | None, ...],
        cycle_info: RefrigerantCycleInfo,
        Q_dot_evp_max: Quantity,
        Q_dot_evp_min: Quantity | None = None,
    ) -> None:
        """
        Parameters
        ----------
        copper_tubes:
            List of copper tubes to select from.
        cycle_info:
            The operating parameters of the refrigeration cycle.
        Q_dot_evp_max:
            Maximum system cooling capacity.
        Q_dot_evp_min:
            Minimum system cooling capacity in case the circuit contains a
            compressor capable of unloading. If there is only one compressor
            that cycles on and off, `Q_dot_evp_min` should be left to `None`.
        """
        self.cycle_info = cycle_info
        self.Q_dot_evp_max = Q_dot_evp_max
        self.Q_dot_evp_min = Q_dot_evp_min
        # Sort copper tubes from small to large internal diameter.
        # noinspection PyTypeChecker
        self.copper_tubes = sorted(
            copper_tubes,
            key=lambda copper_tube: copper_tube.D_int.to('mm').m
        )

    def size(
        self,
        L: Quantity,
        elevation: Quantity = Q_(0.0, 'm')
    ) -> RefrigerantLine | tuple[RefrigerantLine, RefrigerantLine] | None:
        """Sizes the refrigerant line so that the refrigerant flow velocity is
        between the minimum and maximum allowable limit (to ensure proper oil
        return and to prevent noise and pipe wall erosion) while minimizing the
        pressure drop.
        If no suitable single pipe is found and the refrigerant line has a
        positive elevation, a double riser may be returned. In that case there
        should be no fittings in the riser. This means that in case the
        refrigerant line has a vertical riser, the horizontal/vertical drop
        sections and the vertical riser should be sized separately.

        Parameters
        ----------
        L:
            Straight length of the refrigerant line (or segment) without any
            fittings or accessories.
        elevation:
            Height of the line outlet with respect to the line inlet. Only
            relevant if the suction line has a vertical riser.

        Returns
        -------
        -   A single refrigerant line, or
        -   the small and large section of a double riser, or
        -   None in case no copper tube in the list of copper tubes satisfies
            the flow velocity criterion.
        """
        refrigerant_lines = []
        for copper_tube in self.copper_tubes:
            refrigerant_line = self._RFG_LINE(
                cycle_info=self.cycle_info,
                copper_tube=copper_tube,
                L=L,
                Q_dot_evp_max=self.Q_dot_evp_max,
                elevation=elevation,
                Q_dot_evp_min=self.Q_dot_evp_min
            )
            u_max, flag_max, u_min, flag_min = refrigerant_line.flow_velocity()
            if flag_max == 'OK' and flag_min == 'OK':
                dP, _ = refrigerant_line.pressure_drop()
                refrigerant_lines.append(refrigerant_line)
        if refrigerant_lines:
            dP_lst = Quantity.from_list([rl.dP for rl in refrigerant_lines])
            idx_dP_min = dP_lst.magnitude.argmin()
            return refrigerant_lines[idx_dP_min]
        if not refrigerant_lines and elevation.m > 0.0:
            # Try a double riser.
            small_section, large_section = self._create_double_riser(
                L=L,
                elevation=elevation,
            )
            return small_section, large_section
        return None

    def _create_double_riser(
        self,
        L: Quantity,
        elevation: Quantity,
    ) -> tuple[RefrigerantLine, RefrigerantLine]:
        small_section, large_section = None, None
        for copper_tube in self.copper_tubes:
            # The small riser piping is sized based on the minimum capacity of
            # the system.
            small_section_ = self._RFG_LINE(
                cycle_info=self.cycle_info,
                copper_tube=copper_tube,
                L=L,
                Q_dot_evp_max=self.Q_dot_evp_min,
                elevation=elevation,
            )
            u_max, flag_max, u_min, flag_min = small_section_.flow_velocity()
            if flag_min == 'OK':
                small_section = small_section_
            # The larger riser piping is sized based on a maximum capacity equal
            # to the difference between the maximum and minimum system capacity,
            # and a minimum capacity equal to the system's minimum capacity.
            large_section_ = self._RFG_LINE(
                cycle_info=self.cycle_info,
                copper_tube=copper_tube,
                L=L,
                Q_dot_evp_max=self.Q_dot_evp_max - self.Q_dot_evp_min,
                elevation=elevation,
                Q_dot_evp_min=self.Q_dot_evp_min
            )
            u_max, flag_max, u_min, flag_min = large_section_.flow_velocity()
            if flag_max == 'OK' and flag_min == 'OK':
                large_section = large_section_
        # Returns the small riser section with the largest diameter that
        # satisfies the minimum velocity criterion and the large riser section
        # with the largest diameter that satisfies both the minimum and maximum
        # velocity criterion.
        small_section.is_double_riser = True
        large_section.is_double_riser = True
        return small_section, large_section


class SuctionLineSizer(RefrigerantLineSizer):
    """Selects a copper tube with suitable diameter for a suction line."""
    _RFG_LINE = SuctionLine


class DischargeLineSizer(RefrigerantLineSizer):
    """Selects a copper tube with suitable diameter for a discharge line."""
    _RFG_LINE = DischargeLine


class LiquidLineSizer(RefrigerantLineSizer):
    """Selects a copper tube with suitable diameter for a liquid line."""
    _RFG_LINE = LiquidLine

    def size(
        self,
        L: Quantity,
        elevation: Quantity = Q_(0.0, 'm')
    ) -> RefrigerantLine | None:
        """Sizes the liquid line so that the refrigerant flow velocity is below
        the maximum allowable limit (to prevent noise and pipe wall erosion).
        (As oil and liquid refrigerant mix readily, oil return within the liquid
        line is not a concern.)
        Copper tubes are sorted from small to large internal diameter on
        instantiation of the liquid line sizer. The smallest copper tube that
        satisfies the maximum flow velocity criterion is returned as to minimize
        refrigerant charge in the liquid line.

        Parameters
        ----------
        L:
            Straight length of the refrigerant line (or segment) without any
            fittings or accessories.
        elevation:
            Height of the line outlet with respect to the line inlet. Only
            relevant if the suction line has a vertical riser.

        Returns
        -------
        A single refrigerant line, or `None` in case no copper tube in the list
        of available copper tubes satisfies the maximum flow velocity criterion.
        """
        for copper_tube in self.copper_tubes:
            refrigerant_line = self._RFG_LINE(
                cycle_info=self.cycle_info,
                copper_tube=copper_tube,
                L=L,
                Q_dot_evp_max=self.Q_dot_evp_max,
                elevation=elevation,
            )
            # Calculate the flow velocity in the liquid line that correspond
            # with the maximum cooling capacity:
            u_max, flag_max = refrigerant_line.flow_velocity()
            if flag_max == 'OK':
                return refrigerant_line
        return None
    