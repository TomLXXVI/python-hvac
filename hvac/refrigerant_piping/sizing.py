import math
from abc import ABC, abstractmethod
from dataclasses import dataclass
from hvac import Quantity
from hvac.fluids import Fluid
from hvac.refrigerant_piping.copper_tubing import CopperTube
from hvac.heat_transfer.forced_convection.internal_flow import CircularTube


Q_ = Quantity
g = Q_(9.81, 'm / s**2')


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
    u_min_limit: dict[str, Quantity] | None = None

    def __init__(
        self,
        copper_tube: CopperTube,
        elevation: Quantity,
        cycle_info: RefrigerantCycleInfo,
        Q_dot_evp_max: Quantity,
        Q_dot_evp_min: Quantity | None = None
    ) -> None:
        """
        Parameters
        ----------
        copper_tube:
            Instance of dataclass `CopperTube` containing the dimensions of the
            circular cross-section (nominal diameter, outside diameter, inside
            diameter).
        elevation:
            The difference in height between the outlet and inlet of the
            refrigerant line. If the outlet is below the inlet, the elevation
            has a negative value.
        cycle_info:
            Instance of dataclass `RefrigerantCycleInfo` containing operating
            parameters of the refrigerant cycle.
        Q_dot_evp_max:
            Maximum evaporator load. This will determine the maximum
            flow velocity of the refrigerant.
        Q_dot_evp_min: optional, default None
            Minimum evaporator load. This will determine the minimum
            flow velocity of the refrigerant. Leave this `None` in case of
            an on/off-controlled compressor.
        """
        self.copper_tube = copper_tube
        self.elevation = elevation
        self.cycle_info = cycle_info
        self.Q_dot_evp_max = Q_dot_evp_max
        self.Q_dot_evp_min = Q_dot_evp_min or Q_dot_evp_max
        self.m_dot_max = self._get_m_dot(self.Q_dot_evp_max)
        self.V_dot_max = self._get_V_dot(self.m_dot_max)
        self.m_dot_min = self._get_m_dot(self.Q_dot_evp_min)
        self.V_dot_min = self._get_V_dot(self.m_dot_min)

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

    @abstractmethod
    def pressure_drop(self, L_eq: Quantity) -> tuple[Quantity, Quantity]:
        """Returns the pressure loss across the refrigeration line at the
        maximum cooling load due to flow friction and elevation of the outlet
        with respect to the inlet of the line. However, in the case of a
        negative elevation, the associated pressure gain is ignored. Also,
        returns the associated change in saturation temperature that corresponds
        with the pressure drop.

        Parameters
        ----------
        L_eq:
            Total equivalent length of the refrigerant line, including the
            equivalent lengths of any fittings or accessories in the line.
        """
        ...


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
        u_max = self.V_dot_max / A
        r_max = self._check_u_max(u_max)
        u_min = self.V_dot_min / A
        r_min = self._check_u_min(u_min)
        return u_max.to('m / s'), r_max, u_min.to('m / s'), r_min

    def _check_u_max(self, u_max: Quantity) -> str:
        if u_max.to('feet / min') < self.u_max_limit:
            r = 'OK'
        else:
            r = 'TOO HIGH'
        return r

    def _check_u_min(self, u_min: Quantity) -> str:
        u_min = u_min.to('feet / min')
        # minimum allowable flow velocity to ensure oil return:
        u_min_limit = self.u_min_limit[self.copper_tube.DN]
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
    def pressure_drop(self, L_eq: Quantity) -> tuple[Quantity, Quantity]:
        ...


class SuctionLine(VaporLine):
    u_max_limit = Q_(4000, 'ft / min')
    # minimum allowable flow velocity for proper oil return in riser at
    # saturated suction temperature of 20 °F:
    u_min_limit = {
        k: Q_(v, 'feet / min') for k, v in
        [
            ('3/8', 370),
            ('1/2', 460),
            ('5/8', 520),
            ('3/4', 560),
            ('7/8', 600),
            ('1 1/8', 700),
            ('1 3/8', 780),
            ('1 5/8', 840),
            ('2 1/8', 980),
            ('2 5/8', 1080),
            ('3 1/8', 1180),
            ('3 5/8', 1270),
            ('4 1/8', 1360)
        ]
    }

    def _get_V_dot(self, m_dot: Quantity) -> Quantity:
        # volume flow rate of refrigerant at evaporator outlet
        return m_dot / self.cycle_info.rfg_evp_out.rho

    def pressure_drop(self, L_eq: Quantity) -> tuple[Quantity, Quantity]:
        rfg = self.cycle_info.rfg_evp_out
        tube = CircularTube(self.copper_tube.D_int, L_eq, rfg, e=Q_(0.0015, 'mm'))
        tube.m_dot = self.m_dot_max
        dP_friction = tube.pressure_drop()
        if self.elevation.magnitude > 0.0:
            dP_elevation = g * rfg.rho * self.elevation
        else:
            dP_elevation = Q_(0.0, 'Pa')
        dP = dP_friction + dP_elevation
        P = self.cycle_info.P_evp - dP
        T_sat = self.cycle_info.refrigerant(P=P, x=Q_(1.0, 'frac')).T
        dT_sat = self.cycle_info.T_evp - T_sat
        return dP.to('Pa'), dT_sat.to('K')


class DischargeLine(VaporLine):
    u_max_limit = Q_(3500, 'ft / min')
    # minimum allowable flow velocity for proper oil return in riser @
    # saturated condensing temperature of 80 °F:
    u_min_limit = {
        k: Q_(v, 'feet / min') for k, v in
        [
            ('5/16', 220),
            ('3/8', 250),
            ('1/2', 285),
            ('5/8', 315),
            ('3/4', 345),
            ('7/8', 375),
            ('1 1/8', 430),
            ('1 3/8', 480),
            ('1 5/8', 520),
            ('2 1/8', 600),
            ('2 5/8', 665),
            ('3 1/8', 730)
        ]
    }

    def _get_V_dot(self, m_dot: Quantity) -> Quantity:
        # volume flow rate of refrigerant at compressor outlet
        return m_dot / self.cycle_info.rfg_cnd_in.rho

    def pressure_drop(self, L_eq: Quantity) -> tuple[Quantity, Quantity]:
        rfg = self.cycle_info.rfg_cnd_in
        tube = CircularTube(self.copper_tube.D_int, L_eq, rfg, e=Q_(0.0015, 'mm'))
        tube.m_dot = self.m_dot_max
        dP_friction = tube.pressure_drop()
        if self.elevation.magnitude > 0.0:
            dP_elevation = g * rfg.rho * self.elevation
        else:
            dP_elevation = Q_(0.0, 'Pa')
        dP = dP_friction + dP_elevation
        P = self.cycle_info.P_cnd + dP
        T_sat = self.cycle_info.refrigerant(P=P, x=Q_(1.0, 'frac')).T
        dT_sat = T_sat - self.cycle_info.T_cnd
        return dP.to('Pa'), dT_sat.to('K')


class LiquidLine(RefrigerantLine):
    u_max_limit = Q_(600, 'ft / min')

    def _get_V_dot(self, m_dot: Quantity) -> Quantity:
        # volume flow rate of refrigerant at condenser outlet
        return m_dot / self.cycle_info.rfg_cnd_out.rho

    def flow_velocity(self) -> tuple[Quantity, str]:
        A = math.pi * (self.copper_tube.D_int ** 2) / 4
        u_max = self.V_dot_max / A
        r_max = self._check_u_max(u_max)
        return u_max.to('m / s'), r_max

    def pressure_drop(self, L_eq: Quantity) -> tuple[Quantity, Quantity]:
        rfg = self.cycle_info.rfg_cnd_out
        tube = CircularTube(self.copper_tube.D_int, L_eq, rfg, e=Q_(0.0015, 'mm'))
        tube.m_dot = self.m_dot_max
        dP_friction = tube.pressure_drop()
        if self.elevation.magnitude > 0.0:
            dP_elevation = g * rfg.rho * self.elevation
        else:
            dP_elevation = Q_(0.0, 'Pa')
        dP = dP_friction + dP_elevation
        P = self.cycle_info.P_cnd - dP
        T_sat = self.cycle_info.refrigerant(P=P, x=Q_(0.0, 'frac')).T
        dT_sat = self.cycle_info.T_cnd_out - T_sat
        return dP.to('Pa'), dT_sat.to('K')
