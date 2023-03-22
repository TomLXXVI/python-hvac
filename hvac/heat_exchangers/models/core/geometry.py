import numpy as np
from hvac import Quantity

Q_ = Quantity


class CoilGeometry:
    """
    Defines the coil core geometry and determines characteristic heat exchanger parameters for a single row.
    [source: Arora, R. C. & Chandra, A. R. (2012). Refrigeration and Air Conditioning. PHI Learning. p. 520]
    """
    def __init__(
        self,
        sv: Quantity,
        sh: Quantity,
        tf: Quantity,
        sf: Quantity,
        Do: Quantity,
        Di: Quantity,
        kf: Quantity,
        fin_shape: str,
        Afa: Quantity = Q_(1.0, 'm ** 2'),
        Nr: int | float = 1
    ) -> None:
        """
        Parameters
        ----------
        sv: Quantity
            Vertical (transverse) spacing between tubes in a row.
        sh: Quantity
            Horizontal (longitudinal) spacing between the tubes in different rows.
        tf: Quantity
            Thickness of the fins.
        sf: Quantity
            Centre-to-centre spacing between the fins.
        Do: Quantity
            Outer diameter of the tubes.
        Di: Quantity
            Inner diameter of the tubes.
        kf: Quantity
            Thermal conductivity of fin material.
        fin_shape: ['rectangular', 'circular', 'hexagonal']
            Shape of fins.
        Afa: Quantity
            Face area of coil.
        Nr: int
            Number of rows of the coil.

        Notes
        -----
        If Afa and Nr are left to the default values of 1 m² and 1, returned surface areas are valid per unit face area
        for a single row.
        """
        self.sv = sv
        self.sh = sh
        self.tf = tf
        self.sf = sf
        self.Do = Do
        self.Di = Di
        self.kf = kf
        self.fin_shape = fin_shape
        self.Afa = Afa
        self.Nr = Nr

        self._sv = sv.to('mm').m
        self._sh = sh.to('mm').m
        self._tf = tf.to('mm').m
        self._sf = sf.to('mm').m
        self._Do = Do.to('mm').m
        self._Di = Di.to('mm').m

        self.Ap_to_Afa = (self._sf - self._tf) / (self._sf * self._sv) * np.pi * self._Do
        self.Af_to_Afa = 2 / self._sf * (self._sh - np.pi * self._Do ** 2 / (4 * self._sv))
        self.Ac_to_Afa = (self._sf - self._tf) / self._sf * (1 - self._Do / self._sv)
        self.Ao_to_Afa = self.Ap_to_Afa + self.Af_to_Afa
        self.Ai_to_Afa = np.pi * self._Di / self._sv
        self._P = self.Ao_to_Afa / (self._sh / 1000)
        self._Dh = 4 * self.Ac_to_Afa / self._P

    @property
    def Ap(self) -> Quantity:
        """Bare tube area (prime area)."""
        return self.Ap_to_Afa * self.Afa.to('m ** 2') * self.Nr

    @property
    def Af(self) -> Quantity:
        """Fin area."""
        return self.Af_to_Afa * self.Afa.to('m ** 2') * self.Nr

    @property
    def Ac(self) -> Quantity:
        """Minimum flow area."""
        return self.Ac_to_Afa * self.Afa.to('m ** 2')

    @property
    def Ao(self) -> Quantity:
        """Total external heat transfer area."""
        return self.Ao_to_Afa * self.Afa.to('m ** 2') * self.Nr

    @property
    def Dh(self) -> Quantity:
        """Air-side hydraulic diameter."""
        return Q_(self._Dh, 'm')

    @property
    def Ai(self) -> Quantity:
        """Inside heat transfer area."""
        return self.Ai_to_Afa * self.Afa.to('m ** 2') * self.Nr

    @property
    def Ao_to_Ai(self) -> float:
        """Ratio of external to internal heat transfer area."""
        return self.Ao_to_Afa / self.Ai_to_Afa
