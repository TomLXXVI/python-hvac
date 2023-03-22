"""Air-side Fanning friction factor of finned tube heat exchanger when surface is dry."""

from typing import Optional
from hvac import Quantity
from hvac.fluids import HumidAir

Q_ = Quantity


class DrySurfaceFriction:

    def __init__(
        self,
        air: Optional[HumidAir] = None,
        vfa: Optional[Quantity] = None,
        Ac_to_Afa: Optional[float] = None,
        Ao_to_Afa: Optional[float] = None,
        Ap_to_Afa: Optional[float] = None,
        sv: Optional[Quantity] = None,
        Do: Optional[Quantity] = None,
        sf: Optional[Quantity] = None,
        tf: Optional[Quantity] = None,
        correlation: str = 'mcquiston'
    ) -> None:
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
        sv: Quantity
            Vertical (transverse) spacing between tubes in a row.
        Do: Quantity
            Outer diameter of tubes.
        sf: Quantity
            Centre-to-centre spacing between the fins.
        tf: Quantity
            Fin thickness.
        correlation: str, ['mcquiston' (default)]
            The correlation to be used for the calculation of the Fanning friction factor.
        """
        self.air = air
        self.vfa = vfa.to('m / s') if vfa is not None else None
        self.Afa_to_Ac = 1 / Ac_to_Afa if Ac_to_Afa is not None else None
        self.Ao_to_Ap = Ao_to_Afa / Ap_to_Afa if (Ao_to_Afa is not None and Ap_to_Afa is not None) else None
        self.sv = sv.to('m') if sv is not None else None
        self.Do = Do.to('m') if Do is not None else None
        self.sf = sf.to('m') if sf is not None else None
        self.tf = tf.to('m') if tf is not None else None
        self.correlation = correlation

    def _f_fan_mcquiston(self) -> float:
        """Fanning friction factor from McQuiston's correlation.
        [F.C. McQuiston, Correlation of heat, mass and momentum transport
        coefficients for plate-fin–tube heat transfer surfaces with staggered
        tubes, ASHRAE Transact. 84 (1) (1978) 294–309]"""
        vmax = self.vfa * self.Afa_to_Ac
        Re_Do = self.air.rho * vmax * self.Do / self.air.mu
        Do_star = self.Ao_to_Ap * self.Do / (1 + (self.sv - self.Do) / self.sf)
        f1 = (Re_Do ** -0.25) * ((self.Do / Do_star) ** 0.25)
        f2 = (self.sv - self.Do) / (4 * (self.sf - self.tf))
        if self.sv > Do_star:
            f3 = self.sv / Do_star - 1
        else:
            f3 = 1
        FP = f1 * (f2 ** -0.4) * (f3 ** -0.5)
        f_fan = 4.094e-3 + 1.382 * FP ** 2
        return float(f_fan)

    @property
    def f_fan(self) -> float:
        """Get Fanning friction factor."""
        if self.correlation == 'mcquiston':
            return self._f_fan_mcquiston()

    @property
    def f_dar(self) -> float:
        """Get Darcy (aka Moody) friction factor."""
        f_dar = 4 * self.f_fan
        return f_dar

    def setCorrelation(self, corr: str):
        self.correlation = corr
