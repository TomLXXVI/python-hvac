from typing import Tuple
import math
from hvac import Quantity
from hvac.fluid_flow.conduit import Pipe
from ..abstract_fitting import AbstractFitting


Q_ = Quantity


class Tee(AbstractFitting):
    """
    Calculate resistance coefficient of a Tee or Wye according to Crane TP-410M.
    Straight leg and combined leg areas are equal.
    """
    def __init__(self, flow_pattern: str, combined_pipe: Pipe, branch_pipe: Pipe, theta: Quantity = Q_(90.0, 'deg'),
                 ID: str = ''):
        """
        Parameters
        ----------
        flow_pattern: str
            Possible values: 'diverging' or 'converging'.
        combined_pipe: Pipe
            Leg with the greatest volume flow rate.
        branch_pipe: Pipe
            Branch pipe.
        theta: Qty
            Branch leg angle.
        """
        super().__init__(ID)
        self._flow_pattern: str = flow_pattern
        self._D_br: float = branch_pipe.cross_section.hydraulic_diameter.to('mm').m
        self._D_comb: float = combined_pipe.cross_section.hydraulic_diameter.to('mm').m
        self._V_br: float = branch_pipe.volume_flow_rate.to('m ** 3 / s').m
        self._V_comb: float = combined_pipe.volume_flow_rate.to('m ** 3 /s').m
        self._pv_br: float = branch_pipe.velocity_pressure.to('Pa').m
        self._pv_comb: float = combined_pipe.velocity_pressure.to('Pa').m
        self._theta: float = theta.to('deg').m
        self._zeta_combined: float = 0.0
        self._zeta_branch: float = 0.0
        if self._flow_pattern == 'converging':
            self._zeta_combined, self._zeta_branch = self._tee_converging()
        if self._flow_pattern == 'diverging':
            self._zeta_combined, self._zeta_branch = self._tee_diverging()

    def _tee_converging(self) -> Tuple[float, float]:
        """Calculate resistance coefficient of straight leg and branch leg of a converging tee or wye."""
        beta = (self._D_br / self._D_comb) ** 2
        V_rat = self._V_br / self._V_comb
        cb = self._calc_c_branch(beta, V_rat)
        db, eb, fb, cr, dr, er, fr = self._get_coefficients(self._theta)
        zeta_branch = cb * (1.0 + db * (V_rat * 1.0 / beta) ** 2
                            - eb * (1.0 - V_rat) ** 2 - fb * 1.0 / beta * V_rat ** 2)
        if self._theta == 90.0:
            zeta_combined = 1.55 * V_rat - V_rat ** 2
        else:
            zeta_combined = cr * (1.0 + dr * (V_rat * 1.0 / beta) ** 2 - er * (1.0 - V_rat) ** 2
                                  - fr * 1.0 / beta * V_rat ** 2)
        # zeta_branch and zeta_combined are both referenced to the pipe section velocity in the combined leg, but
        # we want zeta_branch to be referenced to the pipe section velocity in the branch leg, therefore:
        zeta_branch = zeta_branch * beta * (1.0 / V_rat) ** 2
        return zeta_combined, zeta_branch

    @staticmethod
    def _calc_c_branch(beta: float, V_rat: float) -> float:
        cb = math.nan
        if beta <= 0.35:
            cb = 1.0
        elif beta > 0.35:
            if V_rat <= 0.4:
                cb = 0.9 * (1.0 - V_rat)
            elif V_rat > 0.4:
                cb = 0.55
        return cb

    @staticmethod
    def _get_coefficients(theta: float) -> Tuple[float, ...]:
        data = {
            30.0: [1.0, 2.0, 1.74, 1.0, 0.0, 1.0, 1.74],
            45.0: [1.0, 2.0, 1.41, 1.0, 0.0, 1.0, 1.41],
            60.0: [1.0, 2.0, 1.0, 1.0, 0.0, 1.0, 1.0],
            90.0: [1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        }
        return tuple(data[theta])

    def _tee_diverging(self) -> Tuple[float, float]:
        """Calculate resistance coefficient of straight leg and branch leg of a diverging tee or wye."""
        beta = (self._D_br / self._D_comb) ** 2
        V_rat = self._V_br / self._V_comb
        g = h = j = 0.0
        if 0.0 < self._theta <= 60.0:
            if beta <= 0.35:
                if V_rat <= 0.4:
                    g = 1.1 - 0.7 * V_rat
                elif V_rat > 0.4:
                    g = 0.85
            elif beta > 0.35:
                if V_rat <= 0.6:
                    g = 1.0 - 0.6 * V_rat
                elif V_rat > 0.6:
                    g = 0.6
            h = 1.0
            j = 2.0
        elif self._theta == 90.0:
            if math.sqrt(beta) <= 2.0 / 3.0:
                g = 1.0
                h = 1.0
                j = 2.0
            elif beta == 1.0 or V_rat * 1.0 / beta <= 2.0:
                g = 1.0 + 0.3 * V_rat ** 2
                h = 0.3
                j = 0.0
        zeta_branch = g * (1.0 + h * (V_rat * 1.0 / beta) ** 2
                           - j * (V_rat * 1.0 / beta) * math.cos(math.radians(self._theta)))
        m = 0.0
        if beta <= 0.4:
            m = 0.4
        elif beta > 0.4:
            if V_rat <= 0.5:
                m = 2.0 * (2.0 * V_rat - 1.0)
            elif V_rat > 0.5:
                m = 0.3 * (2.0 * V_rat - 1.0)
        zeta_combined = m * V_rat ** 2
        # zeta_branch and zeta_combined are both referenced to the pipe section velocity in the combined leg
        # we want zeta_branch to be referenced to the pipe section velocity in the branch leg, therefore:
        zeta_branch = zeta_branch * beta * (1.0 / V_rat) ** 2
        return zeta_combined, zeta_branch

    @property
    def zeta_b(self) -> float:
        """Get resistance coefficient (*float*) of branch leg."""
        return self._zeta_branch

    @property
    def zeta_c(self) -> float:
        """Get resistance coefficient (*float*) of combined leg."""
        return self._zeta_combined

    @property
    def pressure_drop_branch(self) -> Quantity:
        return Quantity(self._zeta_branch * self._pv_br, 'Pa')

    @property
    def pressure_drop_main(self) -> Quantity:
        return Quantity(self._zeta_combined * self._pv_comb, 'Pa')


class Reducer(AbstractFitting):
    """Calculate resistance coefficient of a reducer according to Crane TP-410M."""

    def __init__(self, large_pipe: Pipe, small_pipe: Pipe, L: Quantity, ID: str = ''):
        """
        Parameters
        ----------
        large_pipe: Pipe
            Pipe with large diameter.
        small_pipe: Pipe
            Pipe with small diameter.
        L: Quantity
            Length of reducer

        """
        super().__init__(ID)
        self._D_large: float = large_pipe.cross_section.hydraulic_diameter.to('mm').m  # mm
        self._D_small: float = small_pipe.cross_section.hydraulic_diameter.to('mm').m  # mm
        self._pv_large: float = large_pipe.velocity_pressure.to('Pa').m
        self._pv_small: float = small_pipe.velocity_pressure.to('Pa').m
        self._L: float = L.to('mm').m  # mm
        self._zeta_small, self._zeta_large = self._reducer()

    def _reducer(self) -> Tuple[float, float]:
        beta = self._D_small / self._D_large
        theta = 2.0 * math.atan((self._D_large - self._D_small) / (2.0 * self._L))
        zeta_small = math.nan
        zeta_large = math.nan
        if theta <= math.pi / 4.0:
            zeta_small = 0.8 * math.sin(theta / 2.0) * (1 - beta ** 2)
            zeta_large = zeta_small / beta ** 4
        elif math.pi / 4.0 < theta <= math.pi:
            zeta_small = 0.5 * math.sqrt(math.sin(theta / 2.0)) * (1 - beta ** 2)
            zeta_large = zeta_small / beta ** 4
        return zeta_small, zeta_large

    @property
    def zeta_small(self) -> float:
        """Get resistance coefficient (*float*) of small side of reducer."""
        return self._zeta_small

    @property
    def zeta_large(self):
        """Get resistance coefficient (*float*) of large side of reducer."""
        return self._zeta_large

    @property
    def pressure_drop(self) -> Quantity:
        return Quantity(self._zeta_small * self._pv_small, 'Pa')


class Enlarger(AbstractFitting):
    """Resistance coefficient of an Enlarger according to Crane TP-410M."""

    def __init__(self, large_pipe: Pipe, small_pipe: Pipe, L: Quantity, ID: str = ''):
        """
        Parameters
        ----------
        large_pipe: Pipe
            Large pipe.
        small_pipe: Pipe
            Small pipe.
        L: Quantity
            Length of enlarger.

        """
        super().__init__(ID)
        self._D_large: float = large_pipe.cross_section.hydraulic_diameter.to('mm').m  # mm
        self._D_small: float = small_pipe.cross_section.hydraulic_diameter.to('mm').m  # mm
        self._pv_large: float = large_pipe.velocity_pressure.to('Pa').m
        self._pv_small: float = small_pipe.velocity_pressure.to('Pa').m
        self._L: float = L.to('mm').m    # mm
        self._zeta_small, self._zeta_large = self._enlarger()

    def _enlarger(self) -> Tuple[float, float]:
        beta = self._D_small / self._D_large
        theta = 2 * math.atan((self._D_large - self._D_small) / (2 * self._L))
        zeta_small = math.nan
        zeta_large = math.nan
        if theta <= math.pi / 4.0:
            zeta_small = 2.6 * math.sin(theta / 2) * (1 - beta ** 2) ** 2
            zeta_large = zeta_small / beta ** 4
        elif math.pi / 4.0 < theta <= math.pi:
            zeta_small = (1 - beta ** 2) ** 2
            zeta_large = zeta_small / beta ** 4
        return zeta_small, zeta_large

    @property
    def zeta_small(self) -> float:
        """Get resistance coefficient (*float*) referred to small side of enlarger."""
        return self._zeta_small

    @property
    def zeta_large(self) -> float:
        """Get resistance coefficient (*float*) referred to large side of enlarger."""
        return self._zeta_large

    @property
    def pressure_drop(self) -> Quantity:
        return Quantity(self._zeta_large * self._pv_large, 'Pa')


class ValveReducedPortType1(AbstractFitting):
    """
    Resistance coefficient of a Reduced Port Valve Type 1 according to Crane TP-410M:

    - gate valve
    - ball valve

    """

    def __init__(self, pipe: Pipe, port_diameter: Quantity, theta: Quantity, zeta_unreduced: float, ID: str = ''):
        """
        Parameters
        ----------
        pipe: Pipe
            Pipe in which the valve is placed.
        port_diameter: Quantity
            Diameter of port, i.e. the smallest diameter.
        theta: Qty
            Angle of port reduction.
        zeta_unreduced: float
            Resistance coefficient of valve without reduced port.

        """
        super().__init__(ID)
        self._D_small: float = port_diameter.to('mm').m
        self._D_large: float = pipe.cross_section.hydraulic_diameter.to('mm').m
        self._pv: float = pipe.velocity_pressure.to('Pa').magnitude
        self._theta: float = theta.to('deg').m
        self._zeta_unreduced: float = zeta_unreduced
        self._zeta = self._valve_reduced_port_type1()

    def _valve_reduced_port_type1(self) -> float:
        beta = self._D_small / self._D_large
        zeta_reduced = math.nan
        theta = math.radians(self._theta)
        if beta < 1.0 and theta <= math.pi / 4.0:
            zeta_reduced = ((self._zeta_unreduced + math.sin(theta / 2.0)
                             * (0.8 * (1.0 - beta ** 2) + 2.6 * (1.0 - beta ** 2) ** 2)) / beta ** 4)
        elif beta < 1.0 and math.pi / 4.0 < theta <= math.pi:
            zeta_reduced = ((self._zeta_unreduced + 0.5 * math.sqrt(math.sin(theta / 2.0)) * (1.0 - beta ** 2)
                             + (1.0 - beta ** 2) ** 2) / beta ** 4)
        return zeta_reduced

    @property
    def zeta(self) -> float:
        """Get the resistance coefficient (*float*) of the reduced port valve."""
        return self._zeta

    @property
    def pressure_drop(self) -> Quantity:
        return Quantity(self._zeta * self._pv, 'Pa')


class GateValveReducedPort(ValveReducedPortType1):
    """Derived class of ValveReducedPortType1."""
    pass


class BallValveReducedPort(ValveReducedPortType1):
    """Derived class of ValveReducedPortType1."""
    pass


class ValveReducedPortType2(AbstractFitting):
    """
    Resistance coefficient of a Reduced Port Valve Type 2 according to Crane TP-410M:

    - globe valve
    - angle valve
    - check valve of lift type
    - check valve of stop type
    - plug valve
    - cock

    """
    def __init__(self, pipe: Pipe, port_diameter: Quantity, zeta_unreduced: float, ID: str = ''):
        """
        Parameters
        ----------
        pipe: Pipe
            Pipe in which the valve is placed.
        port_diameter: Quantity
            Diameter of port, i.e. the smallest diameter.
        zeta_unreduced: float
            Resistance coefficient of valve without reduced port.

        """
        super().__init__(ID)
        self._D_small: float = port_diameter.to('mm').m
        self._D_large: float = pipe.cross_section.hydraulic_diameter.to('mm').m
        self._pv: float = pipe.velocity_pressure.to('Pa').magnitude
        self._zeta_unreduced = zeta_unreduced
        self._zeta = self._valve_reduced_port_type2()

    def _valve_reduced_port_type2(self) -> float:
        beta = self._D_small / self._D_large
        zeta_reduced = (self._zeta_unreduced + beta * (0.5 * (1.0 - beta ** 2) + (1.0 - beta ** 2) ** 2)) / beta ** 4
        return zeta_reduced

    @property
    def zeta(self):
        """Get resistance coefficient of reduced port valve."""
        return self._zeta

    @property
    def pressure_drop(self) -> Quantity:
        return Quantity(self._zeta * self._pv, 'Pa')


class GlobeValveReducedPort(ValveReducedPortType2):
    pass


class AngleValveReducedPort(ValveReducedPortType2):
    pass


class LiftCheckValveReducedPort(ValveReducedPortType2):
    pass


class StopCheckValveReducedPort(ValveReducedPortType2):
    pass


class PlugValveReducedPort(ValveReducedPortType2):
    pass


class CockReducedPort(ValveReducedPortType2):
    pass
