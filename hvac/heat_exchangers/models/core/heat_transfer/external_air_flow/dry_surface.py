"""Air side heat transfer (convection) coefficient for finned tube heat exchanger when surface is dry."""

from typing import Optional
from hvac import Quantity
from hvac.fluids import HumidAir


Q_ = Quantity


class DrySurfaceHeatTransfer:
    """
    Air-side heat transfer coefficient
    in case of dry surface.
    """

    def __init__(
        self,
        air: Optional[HumidAir] = None,
        vfa: Optional[Quantity] = None,
        Ac_to_Afa: Optional[float] = None,
        Ao_to_Afa: Optional[float] = None,
        Ap_to_Afa: Optional[float] = None,
        Do: Optional[Quantity] = None,
        Dh: Optional[Quantity] = None,
        correlation: str = 'mcquiston'
    ):
        """
        Parameters
        ----------
        air: HumidAir
            Air properties.
        vfa: Quantity
            Air face velocity.
        Ac_to_Afa: float
            Ratio of minimum flow area to face area.
        Ao_to_Afa: float
            Ratio of external area to face area.
        Ap_to_Afa: float
            Ratio of prime area to face area.
        Do: Quantity
            Outer diameter of tubes.
        Dh: Quantity
            Hydraulic diameter (acc. to property `Dh` of class `HeatExchangerArea`).
        correlation: str, ['mcquiston' (default), 'kays_london']
            The correlation to be used for the calculation of `h`.
        """
        self.air = air
        self.vfa = vfa.to('m / s') if vfa is not None else None
        self.Afa_to_Ac = 1 / Ac_to_Afa if Ac_to_Afa is not None else None
        self.Ao_to_Ap = Ao_to_Afa / Ap_to_Afa if (Ao_to_Afa is not None and Ap_to_Afa is not None) else None
        self.Do = Do.to('m') if Do is not None else None
        self.Dh = Dh.to('m') if Dh is not None else None
        self.correlation = correlation

    def _h_McQuiston(self) -> Quantity:
        """
        Air side heat transfer (convection) coefficient according to McQuiston.
        [Wang, S. K. (2001). Handbook of Air Conditioning and Refrigeration. McGraw-Hill Education. p. 15.41]
        """
        vmax = self.vfa * self.Afa_to_Ac
        Re = self.air.rho * vmax * self.Do / self.air.mu
        JP = (Re ** -0.4) * (self.Ao_to_Ap ** -0.15)
        js = 0.0014 + 0.2618 * JP
        Pr = (self.air.mu * self.air.cp / self.air.k).magnitude
        h = js * self.air.rho * self.air.cp * vmax * (Pr ** (-2 / 3))
        return h

    def _h_KaysLondon(self) -> Quantity:
        """
        Air side heat transfer coefficient according to Kays and London.
        [Arora, R. C., & Chandra, A. R. (2012). Refrigeration and Air Conditioning. PHI Learning. p. 524]
        """
        vmax = self.vfa * self.Afa_to_Ac
        Re = self.air.rho * vmax * self.Dh / self.air.mu
        Pr = self.air.mu * self.air.cp / self.air.k
        Nu = 0.117 * (Re ** 0.65) * (Pr ** (1 / 3))
        h = Nu * self.air.k / self.Dh
        return h

    @property
    def h(self) -> Quantity:
        if self.correlation == 'mcquiston':
            return self._h_McQuiston()
        if self.correlation == 'kays_london':
            return self._h_KaysLondon()

    def setCorrelation(self, corr: str):
        self.correlation = corr
