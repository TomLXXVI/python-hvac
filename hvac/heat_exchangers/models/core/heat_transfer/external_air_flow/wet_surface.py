"""Air side heat transfer coefficient for finned tube heat exchanger when surface is wet."""

from typing import Optional
from hvac import Quantity
from hvac.fluids import HumidAir, CP_HUMID_AIR

Q_ = Quantity


class WetSurfaceHeatTransfer:
    """
    Air-side heat transfer coefficient in case of a wet surface
    (due to condensation of water vapor from air).
    """
    def __init__(
        self,
        air: Optional[HumidAir] = None,
        vfa: Optional[Quantity] = None,
        Ac_to_Afa: Optional[float] = None,
        Ao_to_Afa: Optional[float] = None,
        Ap_to_Afa: Optional[float] = None,
        Do: Optional[Quantity] = None,
        sf: Optional[Quantity] = None,
        tf: Optional[Quantity] = None,
        condensation_regime: Optional[str] = None,
        correlation: Optional[str] = None
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
        sf: Quantity
            Centre-to-centre spacing between the fins.
        tf: Quantity
            Fin thickness.
        condensation_regime: str, ['film' (default), 'droplet']
            Kind of condensation that takes place on the heat transfer surface.
        correlation: str, ['mcquiston' (default), 'pacheco-vega']
            The correlation to be used for the calculation of `hco`.
        """
        self.air = air
        self.vfa = vfa.to('m / s') if vfa is not None else None
        self.Afa_to_Ac = 1 / Ac_to_Afa if Ac_to_Afa is not None else None
        self.Ao_to_Ap = Ao_to_Afa / Ap_to_Afa if (Ao_to_Afa is not None and Ap_to_Afa is not None) else None
        self.Do = Do.to('m') if Do is not None else None
        self.sf = sf.to('m') if sf is not None else None
        self.tf = tf.to('m') if tf is not None else None
        self.condensation_regime = 'film' if condensation_regime is None else condensation_regime
        self.correlation = 'mcquiston' if correlation is None else correlation

    def _hco_McQuiston(self) -> Quantity:
        """source: [Arturo Pacheco-Vega 2016 J. Phys.: Conf. Ser. 745 032052]"""
        vmax = self.vfa * self.Afa_to_Ac
        Re_Do = self.air.rho * vmax * self.Do / self.air.mu
        Re_sf = self.air.rho * vmax * self.sf / self.air.mu
        if self.condensation_regime == 'film':
            # film condensation
            fs = 0.84 + 4e-5 * (Re_sf ** 1.25)
        else:
            # drop-wise condensation
            dr = self.sf / (self.sf - self.tf)
            fs = (0.9 + 4.3e-5 * (Re_sf ** 1.25)) * (1 / dr)
        JP = (Re_Do ** -0.4) * (self.Ao_to_Ap ** -0.15) * fs
        js = 0.0014 + 0.2618 * JP
        Pr = (self.air.mu * self.air.cp / self.air.k).magnitude
        hco = js * self.air.rho * self.air.cp * vmax * (Pr ** (-2 / 3))
        return hco

    def _hco_PachecoVega(self) -> Quantity:
        """source: [Arturo Pacheco-Vega 2016 J. Phys.: Conf. Ser. 745 032052]"""
        vmax = self.vfa * self.Afa_to_Ac
        Re_Do = self.air.rho * vmax * self.Do / self.air.mu
        Re_sf = self.air.rho * vmax * self.sf / self.air.mu
        dr = self.sf / (self.sf - self.tf)
        if self.condensation_regime == 'film':
            # film condensation
            fs = (0.674 - 8.1e-4 * (Re_sf ** 0.737)) * (dr ** 8.6)
        else:
            # drop-wise condensation
            fs = (4.1 - 1119.4e-4 * (Re_sf ** 0.143)) * (dr ** -2.5)
        JP = (Re_Do ** 0.552) * (self.Ao_to_Ap ** 0.16) * fs
        js = 0.0169 - 6.52e-5 * JP
        Pr = (self.air.mu * self.air.cp / self.air.k).magnitude
        hco = js * self.air.rho * self.air.cp * vmax * (Pr ** (-2 / 3))
        return hco

    @property
    def h_co(self) -> Quantity:
        """
        Convection heat transfer coefficient between air and wetted surface.
        """
        if self.correlation == 'mcquiston':
            return self._hco_McQuiston()
        if self.correlation == 'pacheco-vega':
            return self._hco_PachecoVega()

    @property
    def h_tot(self) -> Quantity:
        """
        Total heat transfer coefficient between air and wetted surface (combined sensible heat and mass transfer).
        """
        ht = self.h_co / CP_HUMID_AIR
        return ht

    def setCorrelation(self, corr: str):
        self.correlation = corr
