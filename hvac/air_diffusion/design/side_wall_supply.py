"""SIDE-WALL SUPPLY (ONLY COOLING)

References
----------
[1] Awbi, H. B. (2003). Ventilation of Buildings. Taylor & Francis.
    sec. 6.4.1
[2] Hagström K., Sirén K., Zhivov A. M. (1999). Calculation Methods for Air
    Supply Design in Industrial Facilities. Helsinki University of Technology,
    Laboratory of Heating, Ventilating and Air Conditioning.
"""
import warnings
from dataclasses import dataclass
import numpy as np
from scipy.optimize import root_scalar
from hvac import Quantity
from hvac.fluids import Fluid, FluidState, CoolPropWarning
from hvac.air_diffusion.air_jets.air_jets import CompactJet
from hvac.air_diffusion.air_jets.archimedes_number import archimedes_number

warnings.filterwarnings('ignore', category=CoolPropWarning)

Q_ = Quantity
Air = Fluid('Air')
P_atm = Q_(101_325, 'Pa')
g = Q_(9.81, 'm / s**2')


@dataclass
class RoomInfo:
    """Groups the design information about a room.

    Attributes
    ----------
    L:
        Length of the room, i.e. the dimension of the room in the jet flow
        direction.
    B:
        Width of the room, i.e. the dimension of the room normal to the jet
        flow direction.
    H:
        Height of the room.
    Z:
        Height of the occupied zone (1.8 m for standing or 1.2 m for seated
        occupancy).
    T_r:
        Room setpoint temperature.
    Q_dot:
        Sensible load of the room (which can be the design or maximum load).
        Must be positive for cooling load (heat gain), negative for heating
        load (heat loss).
    """
    L: Quantity
    B: Quantity
    H: Quantity
    Z: Quantity
    T_r: Quantity
    Q_dot: Quantity

    def __post_init__(self) -> None:
        self.L.ito('m')
        self.B.ito('m')
        self.H.ito('m')
        self.Z.ito('m')
        self.T_r.ito('K')
        self.air = Air(T=self.T_r, P=P_atm)
        self.Q_dot.ito('W')

    @property
    def V(self) -> Quantity:
        """Returns the volume of the room."""
        V = self.L * self.B * self.H
        return V

    @property
    def A(self) -> Quantity:
        """Returns the floor area of the room."""
        A = self.L * self.B
        return A

    @property
    def q_dot(self) -> Quantity:
        """Returns the room load per unit floor area."""
        q_dot = self.Q_dot / self.A
        return q_dot


@dataclass
class Output:
    """Groups the results returned by the design procedure of a side-wall
    supply.

    Attributes
    ----------
    N_min:
        Minimum room air change rate based on the room load.
    V_dot:
        Minimum required supply air volume flow rate.
    T_o:
        Required supply air temperature that goes with the minimum required
        volume flow rate.
    dT_o:
        Supply air/room air temperature difference that corresponds with the
        minimum required supply air volume flow rate.
    U_o:
        Supply air velocity at the opening.
    M_o:
        Supply air momentum.
    v_r:
        Mean room air speed.
    L_th:
        Throw length of the ATD to a jet centerline velocity of 0.5 m/s.
    A_o:
        Effective area of the supply opening.
    AR:
        Aspect ratio of the ATD, i.e. the ratio of its length to its width.
    b_o:
        Effective width of the ATD.
    h_o:
        Effective height of the ATD.
    dy_max:
        Tolerable jet drop.
    """
    N_min: Quantity
    V_dot: Quantity
    T_o: Quantity
    dT_o: Quantity
    U_o: Quantity
    M_o: Quantity
    v_r: Quantity
    L_th: Quantity
    A_o: Quantity
    AR: float
    b_o: Quantity
    h_o: Quantity
    dy_max: Quantity

    def __post_init__(self) -> None:
        self.u = {
            'N_min': ('1 / hr', 2),
            'V_dot': ('m**3 / hr', 1),
            'T_o': ('degC', 0),
            'dT_o': ('K', 1),
            'U_o': ('m / s', 2),
            'M_o': ('N', 1),
            'v_r': ('m / s', 2),
            'L_th': ('m', 1),
            'A_o': ('m**2', 4),
            'size': ('mm', 0),
            'dy_max': ('m', 2)
        }

    def __str__(self) -> str:
        out = [
            "minimum room air change rate: "
            f"{self.N_min.to(self.u['N_min'][0]):~P.{self.u['N_min'][1]}f}",
            "minimum supply air volume flow rate: "
            f"{self.V_dot.to(self.u['V_dot'][0]):~P.{self.u['V_dot'][1]}f}",
            "supply air temperature: "
            f"{self.T_o.to(self.u['T_o'][0]):~P.{self.u['T_o'][1]}f}",
            "supply air/room air temperature difference: "
            f"{self.dT_o.to(self.u['dT_o'][0]):~P.{self.u['dT_o'][1]}f}",
            "supply air velocity: "
            f"{self.U_o.to(self.u['U_o'][0]):~P.{self.u['U_o'][1]}f}",
            "supply air momentum: "
            f"{self.M_o.to(self.u['M_o'][0]):~P.{self.u['M_o'][1]}f}",
            "mean room air speed: "
            f"{self.v_r.to(self.u['v_r'][0]):~P.{self.u['v_r'][1]}f}",
            "throw to 0.5 m/s: "
            f"{self.L_th.to(self.u['L_th'][0]):~P.{self.u['L_th'][1]}f}",
            "effective area of ATD: "
            f"{self.A_o.to(self.u['A_o'][0]):~P.{self.u['A_o'][1]}f}",
            f"aspect ratio of ATD: {self.AR:.3f}",
            "effective length of ATD: "
            f"{self.b_o.to(self.u['size'][0]):~P.{self.u['size'][1]}f}",
            "effective height of ATD: "
            f"{self.h_o.to(self.u['size'][0]):~P.{self.u['size'][1]}f}",
            "maximum allowable air jet drop: "
            f"{self.dy_max.to(self.u['dy_max'][0]):~P.{self.u['dy_max'][1]}f}",
        ]
        out = '\n'.join(out)
        return out


class SideWallSupply:
    """
    Class for designing a single side-wall supply.

    Methods
    -------
    design:
        Runs the design procedure starting with an initial value for the
        effective height of the supply opening. The results are returned in
        an `Output` instance (see class `Output` of this module).
    analyze:
        For a supply opening of known size, determines the required minimum
        supply air volume flow rate based on a given value for the critical room
        Archimedes number and calculates the resulting performance
        characteristics. The results are returned in an `Output` instance (see
        class `Output` of this module).
    trajectory:
        Calculates the trajectory of the jet's centerline. It returns for a
        given horizontal distance `x` (or range of horizontal distances),
        besides the vertical distance `y` from the supply opening, the velocity
        and the temperature along the jet's centerline.
    """
    psi = 0.47

    def __init__(
        self,
        room: RoomInfo,
        d: Quantity,
        K1: float,
        Ar_r_crit: float = 1.e4
    ) -> None:
        """Creates a `SideWallSupply` object.

        Parameters
        ----------
        room:
            Instance of `RoomInfo` containing all the design information needed
            about the room.
        d:
            Distance between the upper edge of the supply opening and the
            ceiling.
        K1:
            Dynamic characteristic (throw constant) of the diffuser jet.
        Ar_r_crit:
            Critical room Archimedes number; depends on the room length to room
            height ratio, see e.g. ref. [1], p. 246, Table 6.2. The default
            value applies to a ratio of 2.
        """
        self._room = room
        self._d = d.to('m')
        self._K1 = K1
        self._Ar_r_crit = Ar_r_crit
        self._K2 = CompactJet.calculate_K2(K1)
        self._throw_ratio = 0.75
        self.compact_jet: CompactJet | None = None
        # Class `CompactJet` will be used to calculate the throw, the centerline
        # velocity and temperature along the jet's trajectory.
        self.output: Output | None = None
        # Class `Output` will be used to gather all the calculation results.

    def _supply_volume_flow_rate(self) -> Quantity:
        # Calculates the minimum required volume flow rate of supply air based
        # on the critical room Archimedes number, which depends on the ratio
        # of the room length L (i.e. in the direction of the jet) to the room
        # height.
        Q_dot = self._room.Q_dot.to('J / s')
        T_r = self._room.T_r.to('K')
        rho = self._room.air.rho
        cp = self._room.air.cp
        B = self._room.B.to('m')
        H = self._room.H.to('m')
        n = 2 * g * Q_dot * (B * H) ** 3
        d = self._Ar_r_crit * T_r * rho * cp * (B + H)
        V_dot = (n / d) ** (1 / 3)
        return V_dot.to('m**3 / s')

    def _air_temperature_difference(self, V_dot: Quantity) -> Quantity:
        # Returns the supply/return air temperature difference.
        rho = self._room.air.rho
        cp = self._room.air.cp
        Q_dot = self._room.Q_dot.to('W')
        dT_o = -Q_dot / (rho * cp * V_dot)
        # Q_dot > 0 = cooling load (heat gain), but for cooling
        # dT_o = T_o - T_r < 0.
        return dT_o.to('K')

    def _supply_air_state(self, dT_o: Quantity) -> FluidState:
        supply_air = Air(T=self._room.air.T + dT_o, P=P_atm)
        return supply_air

    def _required_supply_air_momentum(self, v_r: Quantity, rho_o: Quantity) -> Quantity:
        # Calculates the supply air momentum that limits the mean room air speed
        # to the value of `v_r` (see Ref. [2], p. 128, table 6.1).
        B = self._room.B
        H = self._room.H
        rho_r = self._room.air.rho
        M_r = rho_r * v_r ** 2 * (B * H)
        M_o = (1 / 0.78) ** 2 * (rho_o / rho_r) * M_r
        return M_o.to('N')

    @staticmethod
    def _supply_air_momentum(U_o: Quantity, A_o: Quantity, rho_o: Quantity) -> Quantity:
        M_o = rho_o * A_o * U_o ** 2
        return M_o.to('N')

    def _mean_room_air_speed(self, M_o: Quantity, rho_o: Quantity) -> Quantity:
        v_r = 0.78 * np.sqrt(M_o / (rho_o * self._room.B * self._room.H))
        return v_r.to('m / s')

    @staticmethod
    def _effective_area(V_dot: Quantity, M_o_req: Quantity, rho_o: Quantity) -> Quantity:
        # Calculates the required effective area of the supply opening based
        # on the given supply air volume flow rate and the supply air momentum
        # that limits the mean room air speed to the desired value.
        A_o = rho_o * V_dot ** 2 / M_o_req
        return A_o.to('m**2')

    def _maximum_tolerable_drop(self, h_o: Quantity) -> Quantity:
        # Calculates the available height between the center of the supply
        # opening and the upper edge of the occupied zone.
        H = self._room.H
        Z = self._room.Z
        h_o = h_o.to('m')
        d = self._d.to('m')
        delta = 0.075 * self._room.L  # jet spread
        dy_max = H - (Z + d + 0.5 * h_o + delta)
        return dy_max

    def _aspect_ratio(
        self,
        A_o: Quantity,
        dy_max: Quantity
    ) -> Quantity:
        # Calculates the required aspect ratio of the supply opening based on
        # the critical room Archimedes number, the allowable drop of the jet,
        # and the room dimensions. The aspect ratio is calculated so that the
        # trajectory of the jet cuts the upper edge of the occupied zone at a
        # horizontal distance from the supply opening equal to the room length.
        A_o = A_o.to('m ** 2')
        y_max = dy_max.to('m')
        d = self._d.to('m')
        L = self._room.L.to('m')
        B = self._room.B.to('m')
        H = self._room.H.to('m')
        a = self.psi * self._K2 / self._K1 ** 2
        b = d * A_o / y_max
        c = (B + H) / (2 * (B * H) ** 3) * self._Ar_r_crit * L ** 3
        AR = (a * b * c) ** (2 / 3)
        return AR

    @staticmethod
    def _supply_air_velocity(V_dot: Quantity, A_o: Quantity) -> Quantity:
        U_o = V_dot / A_o
        return U_o.to('m / s')

    def _throw(self, U_x: Quantity) -> Quantity:
        # Returns the throw length to velocity `U_x`
        L_th = self.compact_jet.throw_zone3(U_x)
        return L_th.to('m')

    def design(
        self,
        h_o: Quantity,
        v_r_minmax: tuple[Quantity, Quantity] = (Q_(0.1, 'm / s'), Q_(0.3, 'm / s')),
        L_th_frac: float = 0.75
    ) -> Output:
        """Calculates the supply air volume flow rate, the supply/room air
        temperature difference, and the effective size of the supply opening
        based on the critical room Archimedes number, a throw length to 0.5
        m/s being equal to `L_th_frac` x the room length, and an allowable range
        for the mean room air speed.

        Parameters
        ----------
        h_o:
            Initial value for the effective height of the ATD.
        v_r_minmax:
            Minimum and maximum allowable mean air room speed.
        L_th_frac:
            Target throw length to a jet centerline velocity of 0.5 m/s,
            expressed as a fraction of the room length.

        Returns
        -------
        An instance of dataclass `Output` (see its docstring).
        """
        # Calculate the required supply air volume flow rate and corresponding
        # supply air/return air temperature difference according to the critical
        # room Archimedes number and the room load:
        V_dot = self._supply_volume_flow_rate()
        dT_o = self._air_temperature_difference(V_dot)
        supply_air = self._supply_air_state(dT_o)

        # Find, between its minimum and maximum limit, the mean room air speed
        # for which the throw length to 0.5 m/s is equal to `L_th_frac` x the
        # room length.
        U_x = Q_(0.5, 'm / s')

        def _fun(v_r: float, L_th_frac: float) -> float | None:
            v_r = Q_(v_r, 'm / s')
            # Calculate the supply air momentum that corresponds with the
            # mean room air speed:
            M_o = self._required_supply_air_momentum(v_r, supply_air.rho)
            # Determine the effective area of the supply opening:
            A_o = self._effective_area(V_dot, M_o, supply_air.rho)
            # Calculate the jet air velocity at the supply outlet:
            U_o = self._supply_air_velocity(V_dot, A_o)
            if U_o > U_x:
                # Calculate the throw length to 0.5 m/s:
                # Use an instance of class `CompactJet` to calculate jet
                # characteristics (throw, centerline velocity, centerline
                # temperature).
                self.compact_jet = CompactJet(A_o, U_o, supply_air, self._room.air, self._K1)
                L_th = self._throw(U_x)
                dev = (L_th - L_th_frac * self._room.L).to('m').magnitude
                return dev
            return None

        # Check if the targeted throw length can be achieved between the minimum
        # and maximum mean room air speed.
        v_r_min, v_r_max = v_r_minmax
        dev_min = _fun(v_r_min, L_th_frac)
        dev_max = _fun(v_r_max, L_th_frac)
        if all((dev_min, dev_max)) and np.sign(dev_min) != np.sign(dev_max):
            sol = root_scalar(
                _fun,
                args=(L_th_frac,),
                bracket=(
                    v_r_min.to('m / s').m,
                    v_r_max.to('m / s').m
                )
            )
            v_r = Q_(sol.root, 'm/s')
            M_o = self._required_supply_air_momentum(v_r, supply_air.rho)
            A_o = self._effective_area(V_dot, M_o, supply_air.rho)
            U_o = self._supply_air_velocity(V_dot, A_o)
            self.compact_jet = CompactJet(A_o, U_o, supply_air, self._room.air, self._K1)
            L_th = self._throw(U_x)
        else:
            if np.sign(dev_min) < -0.5:
                raise ValueError(
                    f"The throw length factor of {L_th_frac} "
                    f"cannot be attained for mean room air velocities "
                    f"between {v_r_min:~P} and {v_r_max:~P}. "
                    "Try a smaller throw length factor..."
                )
            else:
                raise ValueError(
                    f"The throw length factor of {L_th_frac} "
                    f"cannot be attained for mean room air velocities "
                    f"between {v_r_min:~P} and {v_r_max:~P}. "
                    "Try a larger throw length factor..."
                )

        # Determine the effective height of the supply opening based on the
        # required effective area and the required aspect ratio of the supply
        # opening:
        i_max = 20
        i = 0
        tol = Q_(0.0001, 'm')
        while i < i_max:
            # Calculate the allowable drop of the air jet:
            dy_max = self._maximum_tolerable_drop(h_o)
            # Determine the required aspect ratio of the supply opening:
            AR = self._aspect_ratio(A_o, dy_max)
            # Determine a new value for the effective height:
            h_o_new = np.sqrt(A_o / AR)
            dev = abs(h_o_new - h_o)
            if dev <= tol:
                break
            h_o = h_o_new
            i += 1
        # Determine the effective length of the air supply outlet:
        # noinspection PyUnboundLocalVariable
        b_o = h_o * AR
        if self._d.to('m') > 0.5 * b_o:
            warnings.warn(
                f"the cold jet may not attach to the ceiling with "
                f"minimum length of ATD = {b_o.to('mm'):~P.0f}",
                category=RuntimeWarning
            )
        # Air change rate:
        N_min = V_dot / self._room.V
        # Gather results in an `Output` object:
        # noinspection PyUnboundLocalVariable
        self.output = Output(
            N_min, V_dot, supply_air.T,
            dT_o, U_o, M_o, v_r,
            L_th, A_o, AR, b_o, h_o, dy_max
        )
        return self.output

    def analyze(
        self,
        b_o: Quantity,
        h_o: Quantity,
        L_th_frac: float = 0.75
    ) -> Output:
        """For a given effective size of the supply opening, calculates the
        needed supply air volume flow rate and the needed supply/room air
        temperature difference based on a given critical room Archimedes number.

        Parameters
        ----------
        b_o:
            Effective length of the supply opening.
        h_o:
            Effective height of the supply opening.
        L_th_frac:
            Target throw length to a jet centerline velocity of 0.5 m/s,
            expressed as a fraction of the room length.

        Returns
        -------
        An instance of dataclass `Output` (see its docstring).
        """
        b_o = b_o.to('m')
        h_o = h_o.to('m')
        A_o = b_o * h_o
        AR = b_o / h_o
        # Calculate the required supply air volume flow rate and corresponding
        # supply air/return air temperature difference according to the critical
        # room Archimedes number and the room load:
        V_dot = self._supply_volume_flow_rate()
        dT_o = self._air_temperature_difference(V_dot)
        supply_air = self._supply_air_state(dT_o)
        # Calculate the jet air velocity at the supply outlet, the supply air
        # momentum, and the mean room air speed:
        U_o = self._supply_air_velocity(V_dot, A_o)
        M_o = self._supply_air_momentum(U_o, A_o, supply_air.rho)
        v_r = self._mean_room_air_speed(M_o, supply_air.rho)
        # Use an instance of class `CompactJet` to calculate jet characteristics
        # (throw, centerline velocity, centerline temperature).
        self.compact_jet = CompactJet(A_o, U_o, supply_air, self._room.air, self._K1)
        # Calculate the throw length to 0.5 m/s:
        U_x = Q_(0.5, 'm / s')
        if U_o > U_x:
            L_th = self._throw(U_x)
            L_th_frac_min = L_th_frac - 0.05
            L_th_frac_max = L_th_frac + 0.05
            if L_th < L_th_frac_min * self._room.L:
                warnings.warn(
                    f"throw {L_th:~P.1f} to {U_x:~P.2f} is shorter "
                    f"than {L_th_frac_min * self._room.L:~P.1f}",
                    category=RuntimeWarning
                )
            elif L_th > L_th_frac_max * self._room.L:
                warnings.warn(
                    f"throw {L_th:~P.1f} to {U_x:~P.2f} is longer "
                    f"than {L_th_frac_max * self._room.L.to('m'):~P.1f}",
                    category=RuntimeWarning
                )
        else:
            L_th = Q_(float('nan'), 'm')
        # Calculate the tolerable drop:
        dy_max = self._maximum_tolerable_drop(h_o)
        # Air change rate:
        N_min = (V_dot / self._room.V).to('1 / hr')
        # Gather results in an `Output` object:
        self.output = Output(
            N_min, V_dot, supply_air.T, dT_o, U_o, M_o,
            v_r, L_th, A_o, AR, b_o, h_o, dy_max
        )
        return self.output

    def trajectory(self, x: Quantity, h_o: Quantity | None = None) -> tuple[Quantity, ...]:
        """Given the horizontal distance `x` from the supply opening, returns
        the vertical distance `y` of the jet centerline, the jet centerline
        velocity, and the jet centerline temperature at this horizontal
        distance.

        If the side-wall supply hasn't been calculated yet, parameter `h_o` must
        be assigned an initial value for the effective height of the supply
        opening. Otherwise, this parameter can be left to `None`.

        Parameters
        ----------
        x:
            Horizontal distance from the supply opening.
        h_o:
            An initial value of the effective height of the ATD.
            Ignored when the `design` or `analyze` method has been called on
            this object before.

        Notes
        -----
        The trajectory is calculated based on a modified room Archimedes number
        that also takes into account the dimensions of the supply opening: see
        Awbi, H. B. (2003). Ventilation of Buildings. Taylor & Francis. sec.
        6.2.3, p. 249, eq. 6.47.
        """
        if self.output is None: self.design(h_o)
        dT_o = self.output.dT_o.to('K')
        d = self._d.to('m')
        h_o = self.output.h_o.to('m')
        AR = self.output.AR
        A_o = self.output.A_o.to('m**2')
        U_o = self.output.U_o.to('m / s')
        # Trajectory:
        a = self.psi * self._K2 / self._K1 ** 2
        b = d / (h_o * AR) ** 3
        Ar_o = archimedes_number(np.sqrt(A_o), U_o, dT_o, self._room.air)
        k = a * b * Ar_o
        y = k * (x ** 3)
        # Centerline velocity and centerline temperature:
        U_x = self.compact_jet.centerline_velocity_zone3(x)
        T_x = self.compact_jet.centerline_temperature_zone3(x)
        return y, U_x, T_x


def design_side_wall_supply(
    room_info: RoomInfo,
    h_o: Quantity,
    d: Quantity,
    K1: float,
    v_r_minmax: tuple[Quantity, Quantity] = (Q_(0.1, 'm / s'), Q_(0.3, 'm / s')),
    L_th_frac: float = 0.75,
    Ar_r_crit: float = 1.e4
) -> SideWallSupply:
    """Calculates the supply air volume flow rate, the supply/room air
    temperature difference, and the effective size of the supply opening
    based on a given critical room Archimedes number and an allowable range
    for the mean room air speed.

    Parameters
    ----------
    room_info:
        Instance of `RoomInfo`, containing all the design information needed
        about the room.
    h_o:
        An initial value for the effective height of the supply opening.
    d:
        Distance between the upper edge of the supply opening and the ceiling.
    K1:
        Dynamic characteristic (throw constant) of the diffuser air jet.
    v_r_minmax: optional
        Minimum and maximum allowable mean air room speed.
    L_th_frac: optional
        Target throw length to a jet centerline velocity of 0.5 m/s, expressed
        as a fraction of the room length.
    Ar_r_crit: optional
        Critical room Archimedes number; depends on the room length to room
        height ratio, see e.g. ref. [1], p. 246, Table 6.2. The default value
        applies to a ratio of 2.

    Returns
    -------
    An instance of `SideWallSupply`. Get the design results through its
    `output` attribute (an instance of class `Output`).
    
    Notes
    -----
    The room Archimedes number is the ratio of thermal buoyancy force (caused by 
    the temperature difference between room air and supply air) to the momentum 
    of the supply air jet. If the room Archimedes number is less than the 
    critical value, momentum forces dominate and the supply air jet mixes 
    effectively with the room air. Otherwise, if the room Archimedes number is 
    greater than the critical value, buoyancy forces dominate and the air may
    prematurely rise or fall due to temperature differences, leading to poor
    mixing or thermal discomfort.
    """
    sws_obj = SideWallSupply(room_info, d, K1, Ar_r_crit)
    sws_obj.design(h_o, v_r_minmax, L_th_frac)
    return sws_obj


def analyze_side_wall_supply(
    room_info: RoomInfo,
    b_o: Quantity,
    h_o: Quantity,
    d: Quantity,
    K1: float,
    L_th_frac: float = 0.75,
    Ar_r_crit: float = 1.e4
) -> SideWallSupply:
    """For a given effective size of the supply opening, calculates the
    needed supply air volume flow rate and the needed supply/room air
    temperature difference based on the critical room Archimedes number.

    Parameters
    ----------
    room_info:
        Instance of `RoomInfo`, containing all the design information needed
        about the room.
    b_o:
        Effective length of the supply opening.
    h_o:
        Effective height of the supply opening.
    d:
        Distance between the top of the supply opening and the ceiling.
    K1:
        Dynamic characteristic (throw constant) of the diffuser jet.
    L_th_frac: optional
        Target throw length to a jet centerline velocity of 0.5 m/s,
        expressed as a fraction of the room length.
    Ar_r_crit: optional
        Critical room Archimedes number; depends on the room length to room
        height ratio, see e.g. ref. [1], p. 246, Table 6.2. The default value
        applies to a ratio of 2.

    Returns
    -------
    An instance of `SideWallSupply`. Get the analysis results through its
    `output` attribute (an instance of class `Output`).
    """
    sws_obj = SideWallSupply(room_info, d, K1, Ar_r_crit)
    sws_obj.analyze(b_o, h_o, L_th_frac)
    return sws_obj
