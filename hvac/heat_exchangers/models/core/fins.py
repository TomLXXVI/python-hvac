from typing import Optional
from abc import ABC
import numpy as np
from hvac import Quantity
from .geometry import CoilGeometry


class Fin(ABC):
    def __init__(
        self,
        ri: Optional[Quantity],
        t: Optional[Quantity],
        k: Optional[Quantity],
        h: Quantity,
        eta: Optional[Quantity] = None
    ) -> None:
        """
        Parameters
        ----------
        ri: Quantity, optional
            Outside radius of tube = inside radius of fin.
        t: Quantity, optional
            Fin thickness.
        k: Quantity, optional
            Thermal conductivity of fin material.
        h: Quantity
            Air-side heat transfer coefficient.
        eta: Quantity, optional
            Fin efficiency. Defaults to `None`, which means that fin efficiency
            needs to be calculated from the other parameters. Meant to be used
            only with derived class `DummyFin`.
        """
        self.ri = ri
        self.t = t
        self.k = k
        self.h = h
        self._eta = eta

    @property
    def _r(self) -> Optional[Quantity]:
        return None

    @property
    def efficiency(self) -> Quantity:
        """
        Returns fin efficiency calculated according to the correlation of Schmidt.
        [2017 ASHRAE Handbook - Fundamentals (SI), Chapter 4: Heat Transfer]
        """
        if self._eta is None:
            m = np.sqrt(2 * self.h / (self.k * self.t))
            fi = (self._r - 1) * (1 + 0.35 * np.log(self._r))
            L = self.ri * fi
            eta = np.tanh(m * L) / (m * L)
            return eta.to('frac')
        return self._eta.to('frac')


class CircularFin(Fin):
    def __init__(
        self,
        ri: Quantity,
        ro: Quantity,
        t: Quantity,
        k: Quantity,
        h: Quantity
    ) -> None:
        """
        Parameters
        ----------
        ri: Quantity
            Outside radius of tube.
        ro: Quantity
            Outside radius of fin
        t: Quantity
            Fin thickness.
        k: Quantity
            Thermal conductivity of fin material.
        h: Quantity
            Air-side heat transfer (convection) coefficient.
        """
        self.ro = ro
        super().__init__(ri, t, k, h)

    @property
    def _r(self) -> Quantity:
        _r = self.ro / self.ri
        return _r.to('frac')


class RectangularFin(Fin):
    def __init__(
        self,
        ri: Quantity,
        sv: Quantity,
        sh: Quantity,
        t: Quantity,
        k: Quantity,
        h: Quantity
    ) -> None:
        """
        Parameters
        ----------
        ri: Quantity
            Outside radius of tube.
        sv: Quantity
            Vertical spacing between tubes in a row.
        sh: Quantity
            Horizontal spacing between tubes in different rows.
        t: Quantity
            Fin thickness.
        k: Quantity
            Thermal conductivity of fin material.
        h: Quantity
            Air-side heat transfer (convection) coefficient.
        """
        self.Sv = sv
        self.Sh = sh
        super().__init__(ri, t, k, h)

    @property
    def _r(self) -> Quantity:
        L = max(self.Sv, self.Sh) / 2
        M = min(self.Sv, self.Sh) / 2
        psi = M / self.ri
        beta = L / M
        _r = 1.28 * psi * np.sqrt(beta - 0.2)
        return _r.to('frac')


class HexagonalFin(RectangularFin):
    @property
    def _r(self) -> Quantity:
        L = 0.5 * np.sqrt((self.Sv / 2) ** 2 + self.Sh ** 2)
        M = min(self.Sv / 2, self.Sh)
        psi = M / self.ri
        beta = L / M
        _r = 1.27 * psi * np.sqrt(beta - 0.3)
        return _r.to('frac')


class DummyFin(Fin):
    def __init__(self, h: Quantity, eta: Quantity) -> None:
        super().__init__(None, None, None, h, eta)


class FinnedSurface:
    def __init__(
        self,
        geometry: CoilGeometry,
        h: Quantity,
        eta_fin: Optional[Quantity] = None
    ) -> None:
        """
        Parameters
        ----------
        geometry: CoilGeometry
            Coil geometry object with the geometrical properties of the coil.
        h: Quantity
            Air-side heat transfer coefficient.
        eta_fin: Quantity, optional
            Fin efficiency. If `fin_shape` in `geometry` is None, you can create
            a "dummy fin" with the given fin efficiency `eta_fin`.
        """
        self.geometry = geometry
        self.h = h
        self.eta_fin = eta_fin
        self.fin = self._create_fin()

    @property
    def efficiency(self) -> Quantity:
        Af_to_Ao = self.geometry.Af_to_Afa * (1 / self.geometry.Ao_to_Afa)
        eta = 1 - Af_to_Ao * (1 - self.fin.efficiency)
        return eta

    def _create_fin(self) -> Optional[Fin]:
        fin = None
        if self.geometry.fin_shape == 'circular':
            fin = CircularFin(
                ri=self.geometry.Di / 2,
                ro=self.geometry.Do / 2,
                t=self.geometry.tf,
                k=self.geometry.kf,
                h=self.h
            )
        if self.geometry.fin_shape == 'rectangular':
            fin = RectangularFin(
                ri=self.geometry.Di / 2,
                sv=self.geometry.sv,
                sh=self.geometry.sh,
                t=self.geometry.tf,
                k=self.geometry.kf,
                h=self.h
            )
        if self.geometry.fin_shape == 'hexagonal':
            fin = HexagonalFin(
                ri=self.geometry.Di / 2,
                sv=self.geometry.sv,
                sh=self.geometry.sh,
                t=self.geometry.tf,
                k=self.geometry.kf,
                h=self.h
            )
        if self.geometry.fin_shape is None and self.eta_fin is not None:
            fin = DummyFin(
                h=self.h,
                eta=self.eta_fin
            )
        return fin
