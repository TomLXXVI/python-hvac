"""
VERTICAL UPWARD OR DOWNWARD, ISOTHERMAL OR NON-ISOTHERMAL FREE AIR JETS

References
----------
Hagström K., Sirén K., Zhivov A. M. (1999). Calculation Methods for Air Supply
Design in Industrial Facilities. Helsinki University of Technology, Laboratory
of Heating, Ventilating and Air Conditioning.
"""
import numpy as np
from hvac import Quantity
from hvac.fluids import FluidState
from .archimedes_number import archimedes_number


class CompactJet:
    """Jets formed by cylindrical tubes, nozzles, square or rectangular openings
    with a small aspect ratio, unshaded or shaded by perforated plates, grilles.
    Compact air jets are three-dimensional and axi-symmetric.
    """
    Pr = 0.75  # overall turbulent Prandtl number

    def __init__(
        self,
        A_o: Quantity,
        U_o: Quantity,
        supply_air: FluidState,
        room_air: FluidState,
        is_downward: bool = True,
        K1: float = 7.0,
    ) -> None:
        """Creates a `CompactJet` object.

        Parameters
        ----------
        A_o:
            Effective area of the diffuser opening.
        U_o:
            Average air velocity at diffuser outlet.
        supply_air:
            State of supply air.
        room_air:
            Mean state of room air (return air).
        is_downward:
            Indicates if the jet is projected downward (True) or upward (False).
        K1:
            Dynamic characteristic of the diffuser jet.
        """
        self._A_o = A_o.to('m**2')
        self._U_o = U_o.to('m / s')
        self._supply_air = supply_air
        self._room_air = room_air
        self._is_downward = is_downward
        self._K1 = K1
        # Thermal characteristic of the diffuser jet:
        self._K2 = np.sqrt((1 + self.Pr) / 2) * self._K1
        d_o = np.sqrt(4 * self._A_o / np.pi).to('m')
        dT_o = (self._supply_air.T - self._room_air.T).to('K')
        self._Ar_o = archimedes_number(d_o, self._U_o, dT_o, self._room_air)

    def _sign(self) -> int:
        if self._supply_air.T > self._room_air.T:
            # hot jet
            if self._is_downward:
                return -1
            else:
                return 1
        elif self._supply_air.T < self._room_air.T:
            # cold jet
            if self._is_downward:
                return 1
            else:
                return -1
        else:
            # isothermal jet
            return 0

    def _Kn(self, x: Quantity) -> float:
        # non-isothermal characteristic of the diffuser jet
        Kn = (
            1 + 2.5 * (self._K2 / self._K1 ** 2)
            * self._Ar_o * (x / np.sqrt(self._A_o)) ** 2
        )
        Kn = Kn ** (1 / 3)
        return Kn

    def maximum_velocity(self, x: Quantity) -> Quantity:
        upsilon_x = self._K1 * np.sqrt(self._A_o) * self._Kn(x) / x
        U_x = upsilon_x * self._U_o
        return U_x

    def maximum_temperature(self, x: Quantity) -> Quantity:
        theta_x = self._K2 * np.sqrt(self._A_o) / (self._Kn(x) * x)
        dT_o = (self._supply_air.T - self._room_air.T).to('K')
        T_r = self._room_air.T.to('K')
        T_x = T_r + theta_x * dT_o
        return T_x

    def maximum_travel(self) -> Quantity | None:
        """The maximum travel of a downward projected hot jet or an upward
        projected cold jet.
        """
        if self._sign() < 0:
            k = 0.63 * self._K1 / np.sqrt(self._K2 * self._Ar_o)
            d_o = np.sqrt(4 * self._A_o / np.pi)
            z_max = k * d_o
            return z_max
        return None


class LinearJet:
    """Jets formed by slots or rectangular openings with a large aspect ratio.
    The jet flow is approximately two dimensional and symmetric.
    """
    def __init__(
        self,
        h_o: Quantity,
        U_o: Quantity,
        supply_air: FluidState,
        room_air: FluidState,
        is_downward: bool = True,
        K1: float = 2.43,
        K2: float = 2.00
    ) -> None:
        """Creates a `LinearJet` object.

        Parameters
        ----------
        h_o:
            Height of the slot.
        U_o:
            Average air velocity at diffuser outlet.
        supply_air:
            State of supply air.
        room_air:
            Mean state of room air (return air).
        is_downward:
            Indicates if the jet is projected downward (True) or upward (False).
        K1:
            Dynamic characteristic of the diffuser jet.
        K2:
            Thermal characteristic of the diffuser jet.
        """
        self._h_o = h_o.to('m')
        self._U_o = U_o.to('m / s')
        self._supply_air = supply_air
        self._room_air = room_air
        self._is_downward = is_downward
        self._K1 = K1
        self._K2 = K2
        dT_o = (self._supply_air.T - self._room_air.T).to('K')
        self._Ar_o = archimedes_number(self._h_o, self._U_o, dT_o, self._room_air)

    def _sign(self) -> int:
        if self._supply_air.T > self._room_air.T:
            # hot jet
            if self._is_downward:
                return -1
            else:
                return 1
        elif self._supply_air.T < self._room_air.T:
            # cold jet
            if self._is_downward:
                return 1
            else:
                return -1
        else:
            # isothermal jet
            return 0

    def _Kn(self, x: Quantity) -> float:
        Kn = (
            1 + self._sign() * 1.8 * (self._K2 / self._K1 ** 2)
            * self._Ar_o * (x / np.sqrt(self._h_o)) ** 1.5
        )
        Kn = Kn ** (1 / 3)
        return Kn

    def maximum_velocity(self, x: Quantity) -> Quantity:
        upsilon_x = self._K1 * np.sqrt(self._h_o / x) * self._Kn(x)
        U_x = upsilon_x * self._U_o
        return U_x

    def maximum_temperature(self, x: Quantity) -> Quantity:
        theta_x = self._K2 * np.sqrt(self._h_o / x) / self._Kn(x)
        dT_o = (self._supply_air.T - self._room_air.T).to('K')
        T_r = self._room_air.T.to('K')
        T_x = T_r + theta_x * dT_o
        return T_x

    def maximum_travel(self) -> Quantity | None:
        """The maximum travel of a downward projected hot jet or an upward
        projected cold jet.
        """
        if self._sign() < 0:
            k = 0.67 * self._K1 ** (4 / 3) / (self._K2 * self._Ar_o) ** (2 / 3)
            z_max = k * self._h_o
            return z_max
        return None
