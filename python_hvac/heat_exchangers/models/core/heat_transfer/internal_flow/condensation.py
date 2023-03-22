from typing import Optional, Type
import numpy as np
from hvac import Quantity
from hvac.fluids import Fluid

Q_ = Quantity


class CondensationHeatTransfer:
    """
    Calculation of heat transfer coefficient
    in case of in-tube condensing heat transfer.
    """

    def __init__(
        self,
        Coolant: Optional[Fluid] = None,
        Tc: Optional[Quantity] = None,
        mc: Optional[Quantity] = None,
        Di: Optional[Quantity] = None,
        x: Optional[Quantity] = None,
        P: Optional[Quantity] = None,
        correlation: str = 'shah'
    ):
        """
        Parameters
        ----------
        Coolant: Type[Fluid]
            The type of coolant flowing inside the tubes (to be derived from class `Fluid`).
        Tc: Quantity
            Average coolant temperature at which the coolant properties will be determined.
        mc: Quantity
            The mass flow rate of coolant.
        Di: Quantity
            Inner diameter of the tubes.
        x: Quantity
            Vapor quality of coolant (refrigerant).
        P: Quantity
            Absolute pressure of coolant (refrigerant)
        correlation: str, ['shah' (default), 'ackers_crosser']
            The correlation to be used for the calculation of `hi`.
        """
        self.Coolant = Coolant
        self.Tc = Tc.to('K') if Tc is not None else None
        self.mc = mc.to('kg / s') if mc is not None else None
        self.Di = Di.to('m') if Di is not None else None
        self.x = x.to('frac') if x is not None else None
        self.P = P.to('Pa') if P is not None else None
        self.correlation = correlation

    def _h_AckersCrosser(self) -> Quantity:
        """
        Condensation heat transfer coefficient according to Dean Ackers and Crosser's correlation.
        [Arora, R. C., & Chandra, A. R. (2012). Refrigeration and Air Conditioning. PHI Learning. p. 530 & 544]
        """
        vapor = self.Coolant(T=self.Tc, x=Q_(100, 'pct'))  # saturated vapor
        liquid = self.Coolant(T=self.Tc, x=Q_(0, 'pct'))   # saturated liquid
        Re_g = 4 * self.mc / (np.pi * self.Di * vapor.mu)
        Re_f = 4 * self.mc / (np.pi * self.Di * liquid.mu)
        Pr_f = liquid.cp * liquid.mu / liquid.k
        Re_m = Re_f * (1 + np.sqrt(liquid.rho / vapor.rho))
        if Re_g < 5e4:
            Nu = 5.03 * (Re_m ** (1 / 3)) * (Pr_f ** (1 / 3))
        else:
            Nu = 0.0265 * (Re_m ** 0.8) * (Pr_f ** (1 / 3))
        hi = Nu * liquid.k / self.Di
        return hi

    def _h_Shah(self) -> Quantity:
        """
        Local condensation heat transfer coefficient according to Shah's correlation.
        [Arora, R. C., & Chandra, A. R. (2012). Refrigeration and Air Conditioning. PHI Learning. p. 529 & 544]
        """
        liquid = self.Coolant(T=self.Tc, x=Q_(0, 'pct'))  # saturated liquid
        P_crit = self.Coolant.critical_point.P
        Re_f = 4 * self.mc / (np.pi * self.Di * liquid.mu)
        Re_f.ito('')
        Pr_f = liquid.cp * liquid.mu / liquid.k
        Pr_f.ito('')
        h_f = 0.023 * (Re_f ** 0.8) * (Pr_f ** 0.4) * liquid.k / self.Di
        h_f.ito('W / (m**2 * K)')
        P_red = self.P / P_crit  # reduced pressure
        a = (1 - self.x) ** 0.8
        b = 3.8 * (self.x ** 0.76) * (1 - self.x) ** 0.04 / (P_red ** 0.38)
        hi = h_f * (a + b)
        return hi.to('W / (m**2 * K)')

    @property
    def h(self) -> Quantity:
        if self.correlation == 'shah':
            return self._h_Shah()
        if self.correlation == 'ackers_crosser':
            return self._h_AckersCrosser()

    def setCorrelation(self, corr: str):
        self.correlation = corr
