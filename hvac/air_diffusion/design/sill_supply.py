"""
SILL SUPPLY (HEATING AND COOLING)

References
----------
[1] Awbi, H. B. (2003). Ventilation of Buildings. Taylor & Francis. sec. 6.4.2
"""
from dataclasses import dataclass
import numpy as np
from scipy.optimize import root_scalar
from hvac import Quantity
from .side_wall_supply import RoomInfo


Q_ = Quantity
g = Q_(9.81, 'm / s**2')


@dataclass
class Output:
    """Groups the results returned by the design procedure for a sill supply.

    Attributes
    ----------
    V_dot:
        Required supply air volume flow rate.
    T_o:
        Required supply air temperature.
    U_o:
        Supply air velocity at the opening.
    A_o:
        Effective area of the supply opening.
    h:
        Effective width of the ATD.
    b:
        Effective length of the ATD.
    y_m:
        Maximum drop of the jet from the ceiling.
    """
    V_dot: Quantity
    T_o: Quantity
    U_o: Quantity
    A_o: Quantity
    h: Quantity
    b: Quantity
    y_m: Quantity

    def __post_init__(self) -> None:
        self.u = {
            'V_dot': ('m**3 / hr', 1),
            'T_o': ('degC', 0),
            'U_o': ('m / s', 2),
            'A_o': ('m**2', 4),
            'size': ('mm', 0),
            'drop': ('m', 2)
        }

    def __str__(self) -> str:
        out = [
            "supply air volume flow rate: "
            f"{self.V_dot.to(self.u['V_dot'][0]):~P.{self.u['V_dot'][1]}f}",
            "supply air temperature: "
            f"{self.T_o.to(self.u['T_o'][0]):~P.{self.u['T_o'][1]}f}",
            "supply air velocity: "
            f"{self.U_o.to(self.u['U_o'][0]):~P.{self.u['U_o'][1]}f}",
            "effective area of ATD: "
            f"{self.A_o.to(self.u['A_o'][0]):~P.{self.u['A_o'][1]}f}",
            "effective width of ATD: "
            f"{self.h.to(self.u['size'][0]):~P.{self.u['size'][1]}f}",
            "effective length of ATD: "
            f"{self.b.to(self.u['size'][0]):~P.{self.u['size'][1]}f}",
            "maximum jet drop: "
            f"{self.y_m.to(self.u['drop'][0]):~P.{self.u['drop'][1]}f}"
        ]
        out = '\n'.join(out)
        return out


class SillSupply:

    def __init__(
        self,
        room_info: RoomInfo,
        v_r: Quantity = Q_(0.2, 'm / s')
    ) -> None:
        """Creates a `SillSupply` object.

        Parameters
        ----------
        room_info:
            Instance of `RoomInfo`, containing the design information about the
            room.
        v_r: optional
            An acceptable mean room speed. Default is 0.2 m/s.
        """
        self._room = room_info
        self._v_r = v_r.to('m / s')
        self.output: Output | None = None

    def _supply_volume_flow_rate(self, dT_o: Quantity) -> Quantity:
        rho = self._room.air.rho
        cp = self._room.air.cp
        Q_dot = abs(self._room.Q_dot)
        dT_o = abs(dT_o.to('K'))
        V_dot = Q_dot / (rho * cp * dT_o)
        return V_dot.to('m ** 3 / s')

    def _supply_air_velocity(self, V_dot: Quantity) -> Quantity:
        B = self._room.B.magnitude
        H = self._room.H.magnitude
        v_r = self._v_r.magnitude
        V_dot = V_dot.magnitude
        rho = self._room.air.rho.magnitude
        U_o = Q_(19.75 * B * H * v_r ** 3 / (rho * V_dot), 'm / s')
        return U_o

    @staticmethod
    def _effective_grille_area(V_dot: Quantity, U_o: Quantity) -> Quantity:
        A_o = V_dot / U_o
        return A_o.to('m ** 2')

    @staticmethod
    def _effective_width(A_o: Quantity, b: Quantity) -> Quantity:
        h = A_o / b.to('m')
        return h

    def _archimedes_number(
        self,
        U_o: Quantity,
        h: Quantity,
        dT_o: Quantity
    ) -> float:
        g_ = g.magnitude
        beta = self._room.air.beta.magnitude
        dT_o = dT_o.to('K').magnitude
        h = h.magnitude
        U_o = U_o.magnitude
        Ar = g_ * beta * dT_o * h / U_o ** 2
        return Ar

    @staticmethod
    def _maximum_corner_velocity(
        H_sc: Quantity,  # height between grille and ceiling
        h: Quantity,
        Ar: float,
        U_o: Quantity
    ) -> Quantity:
        H_sc = H_sc.to('m').magnitude
        h = h.magnitude
        k = np.sqrt(5.4 / (H_sc / h)) + 2.15 * Ar * (H_sc / h - 5.4)
        U_c = k * U_o  # maximum jet velocity at the wall/ceiling corner
        return U_c

    @staticmethod
    def _dT_c(
        H_sc: Quantity,
        h: Quantity,
        dT_o: Quantity
    ) -> Quantity:
        H_sc = H_sc.to('m').magnitude
        h = h.magnitude
        k = 1.38 / (H_sc / h) ** 0.25
        dT_c = k * dT_o.to('K')
        # temperature difference between wall/ceiling corner and return air
        return dT_c

    def _maximum_jet_velocity(
        self,
        U_c: Quantity,
        b: Quantity,
        A_o: Quantity
    ) -> Quantity:
        B = self._room.B.magnitude
        b = b.to('m').magnitude
        L = self._room.L.magnitude
        A_o = A_o.magnitude
        k = 1 + 0.042 * (B / b) * (0.75 * L / np.sqrt(A_o) - 5)
        U_m = (1 / k) * U_c  # maximum jet velocity at 0.75 * L
        return U_m

    def _drop(self, dT_c: Quantity, A_o: Quantity, U_c: Quantity) -> Quantity:
        dT_c = abs(dT_c.magnitude)
        L = self._room.L.magnitude
        A_o = A_o.magnitude
        U_c = U_c.magnitude
        n = 8.27e-4 * dT_c * L ** 3
        d = np.sqrt(A_o) * U_c ** 2
        y_m = Q_(n / d, 'm')  # drop below the ceiling at which U_m occurs
        return y_m

    def calculate(self, b: Quantity, H_sc: Quantity, dT_o: Quantity) -> Output:
        """Calculates the air supply rate, the temperature difference between
        supply and room air, and the dimensions of the supply opening.
        The routine starts with the given initial guess of the supply/return air
        temperature difference `dT_o`, which is to be considered as a
        maximum value. It then tries to find the `dT_o` for which the
        maximum jet velocity equals 0.5 m/s at 0.75 x room length.

        Parameters
        ----------
        b:
            Effective length of the grille that's available.
        H_sc:
            Height between the sill and the ceiling.
        dT_o:
            Initial guess of the supply/return air temperature difference.
            Negative for cooling, positive for heating.

        Returns
        -------
        An instance of dataclass `Output` (see its docstring).
        """

        def _fun(dT_o: float) -> float:
            dT_o = Q_(dT_o, 'K')
            V_dot = self._supply_volume_flow_rate(dT_o)
            U_o = self._supply_air_velocity(V_dot)
            A_o = self._effective_grille_area(V_dot, U_o)
            h = self._effective_width(A_o, b)
            Ar = self._archimedes_number(U_o, h, dT_o)
            U_c = self._maximum_corner_velocity(H_sc, h, Ar, U_o)
            U_m = self._maximum_jet_velocity(U_c, b, A_o)
            dev = U_m.magnitude - 0.5
            return dev

        if dT_o.magnitude < 0:  # cooling
            bracket = (dT_o.to('K').m, -0.1)
        else:  # heating
            bracket = (0.5, dT_o.to('K').m)
        try:
            sol = root_scalar(_fun, bracket=bracket)
        except ValueError:
            if dT_o.magnitude < 0:  # cooling
                msg = (
                    "Cannot find a solution for `dT_o` "
                    f"equal or higher than {dT_o.to('K')}"
                )
            else:  # heating
                msg = (
                    "Cannot find a solution for `dT_o` "
                    f"equal or lower than {dT_o.to('K')}"
                )
            raise ValueError(msg) from None
        else:
            dT_o = Q_(sol.root, 'K')
            V_dot = self._supply_volume_flow_rate(dT_o)
            U_o = self._supply_air_velocity(V_dot)
            A_o = self._effective_grille_area(V_dot, U_o)
            h = self._effective_width(A_o, b)
            Ar = self._archimedes_number(U_o, h, dT_o)
            U_c = self._maximum_corner_velocity(H_sc, h, Ar, U_o)
            dT_c = self._dT_c(H_sc, h, dT_o)
            y_m = self._drop(dT_c, A_o, U_c)
            T_o = self._room.T_r + dT_o
            self.output = Output(V_dot, T_o, U_o, A_o, h, b, y_m)
            return self.output


def design_sill_supply(
    room_info: RoomInfo,
    b: Quantity,
    H_sc: Quantity,
    dT_o: Quantity,
    v_r: Quantity = Q_(0.2, 'm / s'),
) -> SillSupply:
    """Calculates the air supply rate, the temperature difference between
    supply and room air, and the dimensions of the supply opening.
    The routine starts with the given initial guess of the supply/return air
    temperature difference `dT_o`, which is to be considered as a
    maximum value. It then tries to find the `dT_o` for which the
    maximum jet velocity equals 0.5 m/s at 0.75 x room length.

    Parameters
    ----------
    room_info:
        Instance of `RoomInfo`, containing the design information about the
        room.
    b:
        Effective length of the grille that's available.
    H_sc:
        Height between the sill and the ceiling.
    dT_o:
        Initial guess of the supply/return air temperature difference.
    v_r: optional
        An acceptable mean room speed. Default is 0.2 m/s.

    Returns
    -------
    An instance of `SillSupply`. Get the design results through its `output`
    attribute (an instance of class `Output`).
    """
    sill_supply = SillSupply(room_info, v_r)
    sill_supply.calculate(b, H_sc, dT_o)
    return sill_supply
