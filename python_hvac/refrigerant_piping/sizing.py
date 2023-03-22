import math
from abc import ABC, abstractmethod
from dataclasses import dataclass
from hvac import Quantity
from hvac.fluids import Fluid
from hvac.refrigerant_piping.copper_tubing import CopperTube
from hvac.fluid_flow import Pipe, Circular


Q_ = Quantity


@dataclass
class RefrigerantCycleData:
    rfg: Fluid
    T_eva: Quantity  # saturated suction temperature = evaporating temperature
    T_con: Quantity  # saturated condensing temperature
    ESH: Quantity    # evaporator superheat
    CSC: Quantity    # condensor subcooling
    CSH: Quantity    # compressor superheat

    def __post_init__(self):
        self.T_eva_o = self.T_eva + self.ESH
        self.T_con_i = self.T_con + self.CSH
        self.T_con_o = self.T_con - self.CSC
        self.P_eva = self.rfg(T=self.T_eva, x=Q_(0, 'frac')).P
        self.P_con = self.rfg(T=self.T_con, x=Q_(0, 'frac')).P
        self.rfg_con_i = self.rfg(P=self.P_con, T=self.T_con_i)
        self.rfg_con_o = self.rfg(P=self.P_con, T=self.T_con_o)
        self.rfg_eva_i = self.rfg(P=self.P_eva, h=self.rfg_con_o.h)
        self.rfg_eva_o = self.rfg(P=self.P_eva, T=self.T_eva_o)
        self.q_eva = self.rfg_eva_o.h - self.rfg_eva_i.h


class RefrigerantLine(ABC):
    vr_allow_max: Quantity = None
    vr_allow_min: dict[str, Quantity] | None = None

    def __init__(
        self,
        rfg_cycle_data: RefrigerantCycleData,
        Q_eva_max: Quantity,
        Q_eva_min: Quantity | None = None,
    ):
        """
        Parameters
        ----------
        rfg_cycle_data:
            Instance of dataclass RefrigerantCycleData containing the
            specifications of the refrigerant cycle.
        Q_eva_max:
            Maximum evaporator capacity. This will determine the maximum
            flow velocity of the refrigerant.
        Q_eva_min: optional, default None
            Minimum evaporator capacity. This will determine the minimum
            flow velocity of the refrigerant. Leave to default `None` in case of
            an ON/OFF-controlled compressor.
        """
        self.rcd = rfg_cycle_data
        self.Q_eva_max = Q_eva_max
        self.Q_eva_min = Q_eva_min or Q_eva_max
        # if Q_eva_min is None, set self.Q_eva_min equal to Q_eva_max
        self.mr_max = self._get_mr(self.Q_eva_max)
        self.Vr_max = self._get_Vr(self.mr_max)
        self.mr_min = self._get_mr(self.Q_eva_min)
        self.Vr_min = self._get_Vr(self.mr_min)

    def _get_mr(self, Q) -> Quantity:
        """Get mass flow rate of refrigerant in the system."""
        return Q / self.rcd.q_eva

    @abstractmethod
    def _get_Vr(self, mr) -> Quantity:
        """Get volume flow rate of refrigerant."""
        ...

    @abstractmethod
    def get_vr(self, *args, **kwargs) -> tuple:
        """Get refrigerant velocity."""
        ...

    def _check_vr_max(self, vr_max):
        if vr_max.to('feet / min') < self.vr_allow_max:
            r = 'OK'
        else:
            r = 'TOO HIGH'
        return r

    def _check_vr_min(self, vr_min, copper_tube, riser):
        ...

    @abstractmethod
    def get_dp(self, copper_tube: CopperTube, Leq: Quantity) -> Quantity:
        """Get pressure drop across refrigeration line."""
        ...


class VaporLine(RefrigerantLine):

    @abstractmethod
    def _get_Vr(self, mr) -> Quantity:
        ...

    def get_vr(self, copper_tube: CopperTube, riser: bool = True) -> tuple:
        """
        Get maximum and minimum flow velocity of refrigerant.
        This method also checks if maximum and minimum velocity are within the
        allowable upper and lower limits to avoid noise on one hand and to
        ensure proper oil return on the other hand.

        Parameters
        ----------
        copper_tube:
            Instance of dataclass CopperTube containing copper tube specs.
        riser: default True
            Indicate if the pipe is a vertical riser (True) or not (False).

        Returns
        -------
        Tuple with 4 elements:
        1. the flow velocity at maximum capacity
        2. string that indicates if maximum flow velocity is 'OK' or 'TOO HIGH'
        3. the flow velocity at minimum capacity
        4. string that indicates if minimum flow velocity is 'OK', 'TOO LOW', or
        'TOO HIGH'
        """
        A = math.pi * (copper_tube.di ** 2) / 4
        vr_max = self.Vr_max / A
        r_max = self._check_vr_max(vr_max)
        vr_min = self.Vr_min / A
        r_min = self._check_vr_min(vr_min, copper_tube, riser)
        return vr_max, r_max, vr_min, r_min

    def _check_vr_max(self, vr_max):
        if vr_max.to('feet / min') < self.vr_allow_max:
            r = 'OK'
        else:
            r = 'TOO HIGH'
        return r

    def _check_vr_min(self, vr_min, copper_tube, riser):
        vr_min = vr_min.to('feet / min')
        vr_allow_min = self.vr_allow_min[copper_tube.dn]
        # minimum allowable flow velocity to ensure oil return
        if not riser:
            vr_allow_min *= 0.75
        if vr_allow_min <= vr_min < self.vr_allow_max:
            r = 'OK'
        elif vr_min < vr_allow_min:
            r = 'TOO LOW'
        else:  # vr_min >= self.vr_allow_max:
            r = 'TOO HIGH'
        return r

    @abstractmethod
    def get_dp(self, copper_tube: CopperTube, Leq: Quantity) -> Quantity:
        ...


class SuctionLine(VaporLine):
    vr_allow_max = Q_(4000, 'ft / min')
    vr_allow_min = {
        # minimum allowable flow velocity for proper oil return in riser @
        # saturated suction temperature of 20 °F
        k: Q_(v, 'feet / min') for k, v in
        [
            ('3/8', 370), ('1/2', 460),
            ('5/8', 520), ('3/4', 560),
            ('7/8', 600), ('1 1/8', 700),
            ('1 3/8', 780), ('1 5/8', 840),
            ('2 1/8', 980), ('2 5/8', 1080),
            ('3 1/8', 1180), ('3 5/8', 1270),
            ('4 1/8', 1360)
        ]
    }

    def _get_Vr(self, mr):
        # volume flow rate of refrigerant at evaporator outlet
        return mr / self.rcd.rfg_eva_o.rho

    def get_dp(self, copper_tube: CopperTube, Leq: Quantity) -> Quantity:
        """
        Get pressure drop across suction line.

        Parameters
        ----------
        copper_tube:
            Instance of dataclass CopperTube containing copper tube specs.
        Leq:
            Equivalent length of suction line (including equivalent length
            of fittings).
        """
        pipe = Pipe.create(
            length=Leq,
            wall_roughness=Q_(0.0015, 'mm'),
            fluid=self.rcd.rfg(P=self.rcd.P_eva, x=Q_(1, 'frac')),
            cross_section=Circular.create(copper_tube.di),
            volume_flow_rate=self.Vr_max
        )
        return pipe.pressure_drop


class DischargeLine(VaporLine):
    vr_allow_max = Q_(3500, 'ft / min')
    vr_allow_min = {
        # minimum allowable flow velocity for proper oil return in riser @
        # saturated condensing temperature of 80 °F
        k: Q_(v, 'feet / min') for k, v in
        [
            ('5/16', 220), ('3/8', 250),
            ('1/2', 285), ('5/8', 315),
            ('3/4', 345), ('7/8', 375),
            ('1 1/8', 430), ('1 3/8', 480),
            ('1 5/8', 520), ('2 1/8', 600),
            ('2 5/8', 665), ('3 1/8', 730)
        ]
    }

    def _get_Vr(self, mr):
        # volume flow rate of refrigerant at condenser inlet
        return mr / self.rcd.rfg_con_i.rho

    def get_dp(self, copper_tube: CopperTube, Leq: Quantity) -> Quantity:
        """
        Get pressure drop across discharge line.

        Parameters
        ----------
        copper_tube:
            Instance of dataclass CopperTube containing copper tube specs.
        Leq:
            Equivalent length of discharge line (including equivalent length
            of fittings).
        """
        pipe = Pipe.create(
            length=Leq,
            wall_roughness=Q_(0.0015, 'mm'),
            fluid=self.rcd.rfg(P=self.rcd.P_con, x=Q_(1, 'frac')),
            cross_section=Circular.create(copper_tube.di),
            volume_flow_rate=self.Vr_max
        )
        return pipe.pressure_drop


class LiquidLine(RefrigerantLine):
    vr_allow_max = Q_(600, 'ft / min')

    def _get_Vr(self, mr) -> Quantity:
        return mr / self.rcd.rfg_con_o.rho

    def get_vr(self, copper_tube: CopperTube) -> tuple:
        A = math.pi * (copper_tube.di ** 2) / 4
        vr_max = self.Vr_max / A
        r_max = self._check_vr_max(vr_max)
        return vr_max, r_max

    def get_dp(self, copper_tube: CopperTube, Leq: Quantity) -> Quantity:
        """
        Get pressure drop across liquid line.

        Parameters
        ----------
        copper_tube:
            Instance of dataclass CopperTube containing copper tube specs.
        Leq:
            Equivalent length of discharge line (including equivalent length
            of fittings).
        """
        pipe = Pipe.create(
            length=Leq,
            wall_roughness=Q_(0.0015, 'mm'),
            fluid=self.rcd.rfg(P=self.rcd.P_con, x=Q_(0, 'frac')),
            cross_section=Circular.create(copper_tube.di),
            volume_flow_rate=self.Vr_max
        )
        return pipe.pressure_drop

    def get_dT(self, dP: Quantity, H: Quantity) -> Quantity:
        """
        Get equivalent change in saturated temperature.

        Parameters
        ----------
        dP:
            Total pressure drop across liquid line including fittings
            and accessoires.
        H:
            The elevation height between the outlet and inlet of the liquid
            line (if the outlet is below the inlet, `H` is negative).
        """
        rho = self.rcd.rfg_con_o.rho
        g = Q_(9.81, 'm / s ** 2')
        dP_elev = rho * g * H
        dP_tot = dP + dP_elev
        P_out = self.rcd.P_con - dP_tot  # pressure at liquid line outlet = TXV inlet
        T_sat = self.rcd.rfg(P=P_out, x=Q_(0, 'frac')).T  # saturation temperature @ P_out
        dT = self.rcd.T_con.to('K') - T_sat.to('K')
        return dT
