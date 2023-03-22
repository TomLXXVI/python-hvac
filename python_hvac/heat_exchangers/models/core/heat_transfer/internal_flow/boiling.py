from typing import Optional
import warnings
import numpy as np
from scipy.optimize import fsolve
from hvac import Quantity
from hvac.fluids import Fluid

warnings.filterwarnings("ignore", category=RuntimeWarning)

Q_ = Quantity


class BoilingHeatTransfer:
    """
    Calculation of heat transfer coefficient
    in case of in-tube boiling heat transfer.
    """

    def __init__(
        self,
        Coolant: Fluid,
        mc: Quantity,
        Tc: Quantity,
        Ts: Optional[Quantity],
        Di: Quantity,
        g: Optional[Quantity] = None,
        x: Optional[Quantity] = None,
        tube_orientation: Optional[str] = None,
        correlation: Optional[str] = None,
        q: Optional[Quantity] = None
    ):
        self.Coolant = Coolant
        self.mc = mc.to('kg / s').m
        self.Tc = Tc.to('K').m
        self.Ts = Ts.to('K').m if Ts is not None else None
        self.Di = Di.to('m').m
        self.g = Q_(9.81, 'm / s ** 2').m if g is None else g.to('m / s ** 2').m
        self.x = Q_(0.5, 'frac').m if x is None else x.to('frac').m
        self.tube_orientation = 'vertical' if tube_orientation is None else tube_orientation
        self.correlation = 'shah' if correlation is None else correlation
        self.q = self._solve_q() if q is None else q.to('W / m ** 2').m

    def _hi_Shah(self, qi: float) -> float:
        """
        Saturated boiling heat transfer coefficient according to Shah (1982).
        [Nellis, G. F., & Klein, S. A. (2020). Introduction to Engineering Heat Transfer. Cambridge University Press.]
        """
        vapor = self.Coolant(T=Q_(self.Tc, 'K'), x=Q_(100, 'pct'))
        liquid = self.Coolant(T=Q_(self.Tc, 'K'), x=Q_(0, 'pct'))
        rho_v = vapor.rho.to('kg / m ** 3').m
        rho_f = liquid.rho.to('kg / m ** 3').m
        muf = liquid.mu.to('Pa * s').m
        cpf = liquid.cp.to('J / (kg * K)').m
        kf = liquid.k.to('W / (m * K)').m
        hfg = (vapor.h - liquid.h).to('J / kg').m
        A = np.pi * (self.Di ** 2) / 4
        G = self.mc / A
        Co = ((1 / self.x - 1) ** 0.8) * ((rho_v / rho_f) ** 0.5)  # convection number
        Bo = qi / (G * hfg)  # boiling number
        Frl = (G ** 2) / ((rho_f ** 2) * self.g * self.Di)  # Froude number for liquid
        Prf = cpf * muf / kf  # Prandtl number for liquid

        # superficial heat transfer coefficient of liquid alone.
        hL = 0.023 * ((G * (1 - self.x) * self.Di / muf) ** 0.8) * (Prf ** 0.4) * kf / self.Di

        if self.tube_orientation == 'vertical':
            N = Co
        else:  # horizontal
            N = Co if Frl > 0.04 else 0.38 * (Frl ** -0.3) * Co

        F = 14.7 if Bo >= 11e-4 else 15.43

        # convective boiling
        psi_cb = 1.8 / (N ** 0.8)
        if N > 1.0:
            # nucleate boiling regime
            psi_nb = 230 * Bo ** 0.5 if Bo >= 0.3e-4 else 1 + 46 * Bo ** 0.5
            psi = max(psi_nb, psi_cb)  # heat transfer enhancement factor
        elif 0.1 < N <= 1.0:
            # bubble suppression regime
            psi_bs1 = F * (Bo ** 0.5) * np.exp(2.74 * (N ** -0.1))
            psi = max(psi_cb, psi_bs1)
        else:
            psi_bs2 = F * (Bo ** 0.5) * np.exp(2.47 * (N ** -0.15))
            psi = max(psi_cb, psi_bs2)
        # two-phase heat transfer coefficient
        hTP = psi * hL
        return hTP

    def _solve_q(self) -> float:

        def _eq(q: np.ndarray) -> np.ndarray:
            q = q[0]
            out = self._hi_Shah(q) * (self.Ts - self.Tc) - q
            return np.array([out])

        roots = fsolve(_eq, np.array([1000.0]))
        q = float(roots[0])
        return q

    @property
    def h(self) -> Quantity:
        if self.correlation == 'shah':
            hi = self._hi_Shah(self.q)
            return Q_(hi, 'W / (m ** 2 * K)')

    def setCorrelation(self, corr: str):
        self.correlation = corr
