"""
CEILING SUPPLY (ONLY COOLING)

References
----------
[1] Awbi, H. B. (2003). Ventilation of Buildings. Taylor & Francis. sec. 6.4.3
"""
import warnings
from dataclasses import dataclass
from enum import StrEnum
import numpy as np
from scipy.optimize import root_scalar
from hvac import Quantity
from hvac.fluids import Fluid, CoolPropWarning
from .side_wall_supply import RoomInfo

Q_ = Quantity

warnings.filterwarnings('ignore', category=CoolPropWarning)

Air = Fluid('Air')
P_atm = Q_(101_325, 'Pa')


@dataclass
class Output:
    """Groups the results returned by the design procedure for a ceiling supply.

    Attributes
    ----------
    dT_o:
        Supply/room air temperature difference.
    V_dot:
        Supply air volume flow rate.
    v_r:
        Mean room speed.
    U_o:
        Supply air velocity at diffuser opening.
    A_o:
        Effective area of circular diffuser. None for a linear diffuser.
    d_o:
        Effective outlet diameter of circular diffuser. None for a linear
        diffuser.
    h:
        Effective slot width of diffuser. For a centrally mounted two-way linear
        diffuser the width of each slot is 0.5 * h.
    T:
        Throw length to a velocity of 0.5 m/s.
    dT_r:
        Room temperature gradient, i.e. the vertical temperature difference
        between the floor and the ceiling.
    """
    dT_o: Quantity
    V_dot: Quantity
    v_r: Quantity
    U_o: Quantity
    A_o: Quantity | None
    d_o: Quantity | None
    h: Quantity
    T: Quantity
    dT_r: Quantity

    def __post_init__(self):
        self.u = {
            'dT_o': ('K', 2),
            'V_dot': ('m**3 / hr', 1),
            'v_r': ('m / s', 3),
            'U_o': ('m / s', 3),
            'A_o': ('m**2', 4),
            'd_o': ('mm', 1),
            'h': ('mm', 1),
            'T': ('m', 3),
            'dT_r': ('K', 1)
        }

    def __str__(self) -> str:
        out = [
            "supply/return air temperature difference: "
            f"{self.dT_o.to(self.u['dT_o'][0]):~P.{self.u['dT_o'][1]}f}",  # 0
            "supply air volume flow rate: "
            f"{self.V_dot.to(self.u['V_dot'][0]):~P.{self.u['V_dot'][1]}f}",  # 1
            "mean room speed: "
            f"{self.v_r.to(self.u['v_r'][0]):~P.{self.u['v_r'][1]}f}",  # 2
            "supply air velocity: "
            f"{self.U_o.to(self.u['U_o'][0]):~P.{self.u['U_o'][1]}f}",  # 3
            "effective slot width: "
            f"{self.h.to(self.u['h'][0]):~P.{self.u['h'][1]}f}",  # 6
            "throw length to 0.5 m/s: "
            f"{self.T.to(self.u['T'][0]):~P.{self.u['T'][1]}f}",  # 7
            "room temperature gradient: "
            f"{self.dT_r.to(self.u['dT_r'][0]):~P.{self.u['dT_r'][1]}f}",  # 8
        ]
        if self.A_o is not None:
            out[4:4] = [
                "effective area of supply opening: "
                f"{self.A_o.to(self.u['A_o'][0]):~P.{self.u['A_o'][1]}f}",  # 4
                "effective outlet diameter: "
                f"{self.d_o.to(self.u['d_o'][0]):~P.{self.u['d_o'][1]}f}",  # 5
            ]
        out = '\n'.join(out)
        return out


class CircularCeilingSupply:
    """
    Size the supply air terminal device in case of circular ceiling diffusers.
    """
    def __init__(self, room_info: RoomInfo) -> None:
        """Creates a `CircularCeiling` object.

        Parameters
        ----------
        room_info:
            Instance of `RoomInfo`, containing the design information about the
            room.
        """
        self._room = room_info
        self.output: Output | None = None

    def _mean_room_air_speed(self, M_o: Quantity) -> Quantity:
        Q_dot = abs(self._room.Q_dot.to('kW').magnitude)
        H = self._room.H.magnitude
        L = min(self._room.L.m, self._room.B.m)
        a = 1.2
        M_o = M_o.magnitude
        k1 = 0.22 * np.sqrt(0.19 * (Q_dot * H / M_o) ** 2 + 1)
        k2 = (L ** 2 / 4 + H ** 2) / (a ** 2 * M_o)
        v_r = Q_(np.sqrt(k1 / k2), 'm / s')
        return v_r

    def _supply_volume_flow_rate(self, dT_o: Quantity) -> Quantity:
        dT_o = abs(dT_o.to('K'))
        Q_dot = abs(self._room.Q_dot)
        rho = self._room.air.rho
        cp = self._room.air.cp
        V_dot = Q_dot / (rho * cp * dT_o)
        return V_dot.to('m ** 3 / s')

    def _effective_area(
        self,
        V_dot: Quantity,
        v_r: Quantity,
    ) -> Quantity:
        rho = self._room.air.rho.m
        V_dot = V_dot.m
        v_r = v_r.m
        A_o = Q_(0.0484 * rho * (V_dot / v_r) ** 2, 'm ** 2')
        return A_o

    @staticmethod
    def _supply_air_velocity(V_dot: Quantity, A_o: Quantity) -> Quantity:
        U_o = V_dot / A_o
        return U_o.to('m / s')

    def _supply_air_momentum(self, A_o: Quantity, U_o: Quantity):
        rho = self._room.air.rho
        M_o = rho * A_o * U_o ** 2
        return M_o.to('N')

    def _find_mean_room_air_speed(self, V_dot: Quantity) -> Quantity:

        def _initial_mean_room_speed() -> Quantity:
            L = min(self._room.L.m, self._room.B.m)
            H = self._room.H.m
            v_r = Q_(0.143 * L / np.sqrt(L ** 2 / 4 + H ** 2), 'm / s')
            return v_r

        v_r = _initial_mean_room_speed()
        tol = 0.001
        i_max = 20
        i = 0
        while i < i_max:
            A_o = self._effective_area(V_dot, v_r)
            U_o = self._supply_air_velocity(V_dot, A_o)
            M_o = self._supply_air_momentum(A_o, U_o)
            v_r_new = self._mean_room_air_speed(M_o)
            dev = abs(v_r_new.m - v_r.m)
            if dev <= tol:
                return v_r_new
            v_r = v_r_new
            i += 1
        else:
            raise ValueError(
                f'no solution for mean room speed after {i_max} iterations'
            )

    @staticmethod
    def _effective_slot_width(A_o: Quantity, d_o: Quantity) -> Quantity:
        h = A_o / (np.pi * d_o)
        return h.to('m')

    @staticmethod
    def _throw_length(
        U_o: Quantity,
        h: Quantity,
        K_v: float
    ) -> Quantity:
        # calculates the throw length to a velocity of 0.5 m/s
        U_o = U_o.magnitude
        h = h.magnitude
        T = Q_(4 * (K_v ** 2) * (U_o ** 2) * h, 'm')
        return T
    
    def _room_temperature_gradient(self, M_o: Quantity) -> Quantity:
        # vertical temperature difference between floor and ceiling
        M_o = M_o.magnitude
        V = self._room.V.magnitude
        H = self._room.H.magnitude
        Q_dot = abs(self._room.Q_dot.to('kW').magnitude)
        dT_r = Q_(4.41 * (M_o / V) * (Q_dot * H / M_o) ** 2, 'K')
        return dT_r

    def calculate(
        self,
        dT_o: Quantity,
        d_o: Quantity,
        K_v: float = 1.05
    ) -> Output:
        """Calculates the required air supply rate and the dimensions of the
        supply opening.

        The routine starts with the given initial guess for the diffuser outlet
        diameter `d_o`, which is to be considered as a maximum value. It then
        tries to find the actual outlet diameter for which the throw length
        equals 0.375 * room length.

        Parameters
        ----------
        dT_o:
            The selected supply/return air temperature difference.
            Negative for cooling.
        d_o:
            Initial guess of the diffuser outlet diameter.
        K_v: optional
            Throw constant. Default value applies to a circular cone diffuser.

        Returns
        -------
        An instance of dataclass `Output` that bundles the calculation
        results (see its docstring).
        """
        V_dot = self._supply_volume_flow_rate(dT_o)
        v_r = self._find_mean_room_air_speed(V_dot)
        if v_r > Q_(0.25, 'm / s'):
            warnings.warn(
                f"the mean room air speed {v_r:~P.2f} exceeds 0.25 m/s",
                category=RuntimeWarning
            )
        A_o = self._effective_area(V_dot, v_r)
        U_o = self._supply_air_velocity(V_dot, A_o)
        M_o = self._supply_air_momentum(A_o, U_o)
        dT_r = self._room_temperature_gradient(M_o)

        def _fun(d_o: float) -> float:
            d_o = Q_(d_o, 'm')
            h = self._effective_slot_width(A_o, d_o)
            T = self._throw_length(U_o, h, K_v)
            dev = T.m - 0.375 * self._room.L.m
            return dev

        try:
            sol = root_scalar(_fun, bracket=(0.01, d_o.to('m').m))
        except ValueError:
            msg = f"no solution for outlet diameter {d_o:~P.3f} or smaller"
            raise ValueError(msg) from None
        else:
            d_o = Q_(sol.root, 'm')
            h = self._effective_slot_width(A_o, d_o)
            T = self._throw_length(U_o, h, K_v)
            self.output = Output(dT_o, V_dot, v_r, U_o, A_o, d_o, h, T, dT_r)
            return self.output


def design_circular_ceiling_supply(
    room_info: RoomInfo,
    dT_o: Quantity,
    d_o: Quantity,
    K_v: float = 1.05
) -> CircularCeilingSupply:
    """Calculates the required air supply rate and the dimensions of the
    supply opening.

    The routine starts with the given initial guess for the diffuser outlet
    diameter `d_o`, which is to be considered as a maximum value. It then
    tries to find the actual outlet diameter for which the throw length
    equals 0.375 * room length.

    Parameters
    ----------
    room_info:
        Instance of `RoomInfo`, containing the design information about the
        room.
    dT_o:
        The selected supply/return air temperature difference.
        Negative for cooling.
    d_o:
        Initial guess of the diffuser outlet diameter.
    K_v: optional
        Throw constant. Default value applies to a circular cone diffuser.

    Returns
    -------
    An instance of `CircularCeilingSupply`. Get the design results through its
    `output` attribute (an instance of class `Output`).
    """
    ccs_obj = CircularCeilingSupply(room_info)
    ccs_obj.calculate(dT_o, d_o, K_v)
    return ccs_obj


class LinearCeilingSupply:
    """
    Size the supply air terminal device in case of linear ceiling diffusers.
    """
    def __init__(
        self,
        room_info: RoomInfo,
        position: str = 'central'
    ) -> None:
        """Creates a `LinearCeilingSupply` object.

        Parameters
        ----------
        room_info:
            Instance of `RoomInfo`, containing the design information about the
            room.
        position: {'end', 'central' (default)}
            Indicates if the linear diffuser is mounted at one end of the
            ceiling ('end'), discharging air horizontally in one direction, or
            centrally ('central'), discharging equal volumes of air in two
            horizontal directions.
        """
        self._room = room_info
        self.output: Output | None = None
        if position == 'central':
            self._L = self._room.L / 2
        else:
            self._L = self._room.L

    def _supply_volume_flow_rate(self, dT_o: Quantity) -> Quantity:
        q_dot = self._room.q_dot
        L = self._L
        rho = self._room.air.rho
        cp = self._room.air.cp
        dT_o = dT_o.to('K')
        V_dot = q_dot * L / (rho * cp * dT_o)
        return V_dot.to('m ** 3 / s / m')

    def _effective_slot_width(self, V_dot: Quantity, K_v: float) -> Quantity:

        def _mean_room_speed() -> Quantity:
            # calculates the mean room speed that would correspond with
            # a recommended throw length of 0.75 * room length
            T = (0.75 * self._L).magnitude
            B = self._room.B.magnitude
            rho = self._room.air.rho.magnitude
            v_r = Q_(np.sqrt(T * 0.0484 * rho * B / (4 * K_v ** 2)), 'm / s')
            return v_r

        v_r = _mean_room_speed()
        rho = self._room.air.rho.magnitude
        B = self._room.B.magnitude
        V_dot = V_dot.magnitude
        v_r = v_r.magnitude
        h = Q_(0.0484 * rho * B * (V_dot / v_r) ** 2, 'm')
        return h

    @staticmethod
    def _supply_air_velocity(V_dot: Quantity, h: Quantity) -> Quantity:
        U_o = V_dot / h
        return U_o.to('m / s')

    def _supply_air_momentum(self, U_o: Quantity, h: Quantity) -> Quantity:
        rho = self._room.air.rho
        B = self._room.B
        M_o = rho * B * h * U_o ** 2
        return M_o.to('N')

    def _mean_room_air_speed(self, M_o: Quantity) -> Quantity:
        # calculates the mean room air speed that depends on the supply air
        # momentum, the room load, and the height of the room.
        Q_dot = abs(self._room.Q_dot.to('kW').magnitude)
        H = self._room.H.magnitude
        L = self._room.L.magnitude
        a = 1.0
        M_o = M_o.magnitude
        k1 = 0.22 * np.sqrt(0.19 * (Q_dot * H / M_o) ** 2 + 1)
        k2 = (L ** 2 / 4 + H ** 2) / (a ** 2 * M_o)
        v_r = Q_(np.sqrt(k1 / k2), 'm / s')
        return v_r

    @staticmethod
    def _throw_length(U_o: Quantity, h: Quantity, K_v: float) -> Quantity:
        # calculates the throw length to a velocity of 0.5 m/s
        U_o = U_o.magnitude
        h = h.magnitude
        T = Q_(4 * (K_v ** 2) * (U_o ** 2) * h, 'm')
        return T

    def _room_temperature_gradient(self, M_o: Quantity) -> Quantity:
        # vertical temperature difference between floor and ceiling
        M_o = M_o.magnitude
        V = self._room.V.magnitude
        H = self._room.H.magnitude
        Q_dot = abs(self._room.Q_dot.to('kW').magnitude)
        dT_r = Q_(4.41 * (M_o / V) * (Q_dot * H / M_o) ** 2, 'K')
        return dT_r

    def calculate(self, dT_o: Quantity, K_v: float = 2.35) -> Output:
        """Calculates for a selected supply/return air temperature difference
        `dT_o`, the required supply air volume flow rate to compensate the room
        load, the required effective slot width of the linear diffuser to
        establish a throw to a velocity of 0.5 m/s at 0.75 * room_length, the
        resulting mean room speed, and the resulting vertical temperature
        difference between floor and ceiling in the room.

        Parameters
        ----------
        dT_o:
            Selected supply/return air temperature difference.
        K_v:
            Throw constant of the linear diffuser. Default is 2.35.

        Returns
        -------
        An instance of dataclass `Output` that bundles the calculation
        results (see its docstring).
        """
        dT_o = abs(dT_o.to('K'))
        V_dot = self._supply_volume_flow_rate(dT_o)
        h = self._effective_slot_width(V_dot, K_v)
        U_o = self._supply_air_velocity(V_dot, h)
        T = self._throw_length(U_o, h, K_v)
        M_o = self._supply_air_momentum(U_o, h)
        v_r = self._mean_room_air_speed(M_o)
        dT_r = self._room_temperature_gradient(M_o)
        self.output = Output(
            dT_o, V_dot * self._room.B, v_r,
            U_o, None, None, h, T, dT_r
        )
        return self.output


def design_linear_ceiling_supply(
    room_info: RoomInfo,
    dT_o: Quantity,
    K_v: float = 2.35,
    position: str = 'central'
) -> LinearCeilingSupply:
    """Calculates for a selected supply/return air temperature difference
    `dT_o`, the required supply air volume flow rate to compensate the room
    load, the required effective slot width of the linear diffuser to
    establish a throw to a velocity of 0.5 m/s at 0.75 * room_length, the
    resulting mean room speed, and the resulting vertical temperature
    difference between floor and ceiling in the room.

    Parameters
    ----------
    room_info:
        Instance of `RoomInfo`, containing the design information about the
        room.
    dT_o:
        Selected supply/return air temperature difference.
    K_v:
        Throw constant of the linear diffuser. Default is 2.35.
    position: {'end', 'central' (default)}
        Indicates if the linear diffuser is mounted at one end of the
        ceiling ('end'), discharging air horizontally in one direction, or
        centrally ('central'), discharging equal volumes of air in two
        horizontal directions.

    Returns
    -------
    An instance of `LinearCeilingSupply`. Get the design results through its
    `output` attribute (an instance of class `Output`).
    """
    lcs_obj = LinearCeilingSupply(room_info, position)
    lcs_obj.calculate(dT_o, K_v)
    return lcs_obj


class DiffuserType(StrEnum):
    CIRCULAR = 'circular'
    LINEAR_CONT = 'continuous slot'
    LINEAR_INTER = 'intermittent slot'


class VAVSupply:
    """In a VAV system the supply air volume flow rate is modulated to meet
    the room load at part-load conditions, while the supply air temperature
    `T_o` and consequently the supply air/room air temperature difference `dT_o`
    remain fixed.
    With this class the minimum room load can be determined that can be handled
    by the ceiling supply for a given `dT_o`, based on the criterion that the
    ratio `M_o / (Q_dot * H)` must be equal or greater than 0.07
    where `M_o` is the supply air momentum in Newtons, `Q_dot` is the actual
    room load in Watts, and `H` is the room height in metres. If this ratio
    should become smaller than 0.07, the convective air currents produced by
    internal heat loads dominate the room air movement and the temperature
    gradient in the occupied zone increases as a result.
    This class can also be used to determine:
    (1) the maximum local air velocity in the room for a given room load,
    (2) the mean air velocity in the room for a given room load,
    (3) the room load for which a certain upper limit of the maximum local
    air velocity in the room is reached,
    (4) the room load for which a certain lower or upper limit of the mean air
    velocity in the room is reached.
    """
    def __init__(
        self,
        room: RoomInfo,
        A_o: Quantity,
        dT_o: Quantity,
        diffuser_type: DiffuserType = DiffuserType.CIRCULAR,
        Z: Quantity = Q_(1.8, 'm'),
        r: float | None = None,
    ) -> None:
        """Creates a `VAVSupply` object.

        Parameters
        ----------
        room:
            Instance of `RoomInfo`, containing the design information about the
            room.
        A_o:
            Effective area of the supply air terminal.
        dT_o:
            Supply/room air temperature difference, defined as the difference
            between the supply air temperature and the room air temperature.
            - cooling: dT_o < 0
            - heating: dT_o > 0
        diffuser_type: optional
            Type of diffuser, see StrEnum class `DiffuserType`.
        Z: optional
            Height of the occupied zone.
        r: optional
            Only for intermittent slot diffusers: ratio of active length of
            diffuser to room length.
        """
        self.room = room
        self.L = room.L.to('m')
        self.B = room.B.to('m')
        self.H = room.H.to('m')
        self.Q_dot = room.Q_dot.to('W')
        self.T_r = room.T_r.to('K')
        self.A_o = A_o.to('m**2')
        self.dT_o = dT_o.to('K')
        self.diffuser_type = diffuser_type
        self.Z = Z.to('m')
        if self.diffuser_type is DiffuserType.LINEAR_INTER and r is None:
            raise AttributeError(
                'Parameter `r` must be specified in case '
                'of an intermittent slot diffuser.'
            )
        else:
            self.r = r

        self.tau = self._calculate_tau()
        self.Q_dot_min = self.minimum_room_load()
        self.room_air = Air(T=self.T_r, P=P_atm)
        self.T_o = self.T_r + dT_o
        self.supply_air = Air(T=self.T_o, P=P_atm)
        self.V_dot = self.supply_air_volume_flow_rate()
        self.M_o = self.supply_air_momentum()
        self.v_m = self.maximum_local_velocity()
        self.v_r = self.mean_room_air_speed()

    def supply_air_volume_flow_rate(self) -> Quantity:
        """Returns the supply air volume flow rate that corresponds with the
        given room load and the supply/room air temperature difference.
        """
        Q_dot = abs(self.Q_dot)
        dT_o = abs(self.dT_o)
        rho_o = self.supply_air.rho
        cp_o = self.supply_air.cp
        V_dot = Q_dot / (rho_o * cp_o * dT_o)
        return V_dot.to('m**3 / s')

    def supply_air_momentum(self) -> Quantity:
        """Returns the supply air momentum that corresponds with the
        given room load and the supply/room air temperature difference.
        """
        rho_o = self.supply_air.rho
        M_o = rho_o * (self.V_dot ** 2) / self.A_o
        return M_o.to('N')

    def _calculate_tau(self) -> Quantity:
        tau = np.sqrt(self.A_o) * abs(self.dT_o)
        return tau.to('K * m')

    def minimum_room_load(self) -> Quantity:
        """Returns the minimum room load the ceiling supply can handle,
        which depends on the effective area of the supply air opening, the
        supply/room air temperature difference, and the height of the room.
        """
        c_min = Q_(0.07, 'N / (kW * m)')
        std_air = Air(T=Q_(20, 'degC'), P=P_atm)
        Q_dot_min = c_min * std_air.rho * (std_air.cp ** 2) * (self.tau ** 2) * self.H
        return Q_dot_min.to('W')

    def maximum_local_velocity(self, Q_dot: Quantity | None = None) -> Quantity:
        """Returns the maximum local room air velocity that corresponds with the
        room load `Q_dot`. If `Q_dot` is None, the value of `self.Q_dot` is
        taken.
        """
        if self.L != self.B:
            L = np.sqrt(self.L * self.B)
        else:
            L = self.L
        a = (0.5 * L + self.H - self.Z).magnitude
        if Q_dot is None:
            Q_dot = self.Q_dot.to('kW').magnitude
        else:
            Q_dot = Q_dot.to('kW').magnitude
        tau = self.tau.magnitude
        L = L.magnitude
        v_m = float('nan')
        match self.diffuser_type:
            case DiffuserType.CIRCULAR:
                v_m = 0.9 * Q_dot / (tau * a)
            case DiffuserType.LINEAR_CONT:
                v_m = 1.06 * Q_dot / (tau * np.sqrt(L * a))
            case DiffuserType.LINEAR_INTER:
                v_m = 0.88 * Q_dot / (tau * np.sqrt(self.r * L * a))
        return Q_(v_m, 'm / s')

    def mean_room_air_speed(self, Q_dot: Quantity | None = None) -> Quantity:
        """Returns the mean room air velocity that corresponds with the room
        load `Q_dot`. If `Q_dot` is None, the value of `self.Q_dot` is taken.
        """
        a = 1.0
        if self.diffuser_type is DiffuserType.CIRCULAR:
            a = 1.2
        if Q_dot is None:
            Q_dot = self.Q_dot.to('kW').magnitude
        else:
            Q_dot = Q_dot.to('kW').magnitude
        H = self.H.magnitude
        tau = self.tau.magnitude
        L = self.L.magnitude
        n = 0.42 * a * Q_dot * (0.3 * (H * tau ** 2 / Q_dot) ** 2) ** 0.25
        d = tau * np.sqrt(0.25 * L ** 2 + H ** 2)
        v_r = Q_(n / d, 'm / s')
        return v_r

    def _find_room_load_01(self, v_m: Quantity) -> Quantity:
        """Searches for the room load for which the maximum local room air
        velocity is equal to `v_m`. The search range is between 0.5 x the
        minimum room load that the ceiling supply can handle and 2 x the given
        design room load. If no room load can be found within the search range
        NaN is returned.
        """
        v_m_target = v_m.to('m / s')

        def _fun(Q_dot: float) -> float:
            Q_dot = Q_(Q_dot, 'W')
            v_m = self.maximum_local_velocity(Q_dot)
            dev = (v_m - v_m_target).magnitude
            return dev

        Q_dot_min = 0.5 * self.minimum_room_load().to('W')
        Q_dot_max = 2.0 * self.Q_dot
        try:
            sol = root_scalar(_fun, bracket=(Q_dot_min.m, Q_dot_max.m))
        except ValueError:
            return Q_(float('nan'), 'W')
        else:
            Q_dot = Q_(sol.root, 'W')
            return Q_dot

    def _find_room_load_02(self, v_r: Quantity) -> Quantity:
        """Searches for the room load for which the mean room air velocity is
        equal to `v_r`. The search range is between 0.5 x the minimum room
        load that the ceiling supply can handle and 2 x the given design room
        load. If no room load can be found within the search range NaN is
        returned.
        """
        v_r_target = v_r.to('m / s')

        def _fun(Q_dot: float) -> float:
            Q_dot = Q_(Q_dot, 'W')
            v_r = self.maximum_local_velocity(Q_dot)
            dev = (v_r - v_r_target).magnitude
            return dev

        Q_dot_min = self.minimum_room_load().to('W')
        Q_dot_max = 2.0 * self.Q_dot
        try:
            sol = root_scalar(_fun, bracket=(Q_dot_min.m, Q_dot_max.m))
        except ValueError:
            return Q_(float('nan'), 'W')
        else:
            Q_dot = Q_(sol.root, 'W')
            return Q_dot

    def get_room_load_limits(
        self,
        v_m: Quantity = Q_(0.35, 'm / s'),
        v_r_min: Quantity = Q_(0.1, 'm / s'),
        v_r_max: Quantity = Q_(0.25, 'm / s')
    ) -> tuple[Quantity, ...]:
        """Returns the room loads that correspond with (1) the maximum local
        room velocity, (2) the minimum and (3) the maximum mean room air
        velocity.
        If a room load cannot be found within the search range of 0.5 x the
        minimum room load that the ceiling supply can handle and 2 x the given
        design room load, a NaN is returned.
        """
        Q_dot_1 = self._find_room_load_01(v_m)
        Q_dot_2 = self._find_room_load_02(v_r_min)
        Q_dot_3 = self._find_room_load_02(v_r_max)
        return Q_dot_1, Q_dot_2, Q_dot_3
