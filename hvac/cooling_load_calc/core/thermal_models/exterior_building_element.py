"""LINEAR THERMAL NETWORK OF AN OPAQUE EXTERIOR BUILDING ELEMENT.
-----------------------------------------------------------------
Implementation of a linear thermal network (lumped capacitance model) using
a 2nd-order non-central backward finite difference approximation [1] for the
derivative dT/dt in the node equations.

This implementation is further used in subpackage `cooling_load_calc` for
calculating the conduction heat gain through opaque exterior building elements
in cooling load calculations (see class `ExteriorBuildingElement` in
`hvac.cooling_load_calc.core.building_element` and class `ConditionedZone` in
`hvac.cooling_load_calc.building.conditioned_zone`).

References
----------
[1] Kiusalaas, J. (2013).
    NUMERICAL METHODS IN ENGINEERING WITH PYTHON3.
    Cambridge University Press.
"""
from collections.abc import Sequence, Callable
from abc import ABC, abstractmethod
import numpy as np
import pandas as pd
from hvac import Quantity
from hvac.cooling_load_calc.core.construction_assembly import (
    ConstructionAssembly,
    ConstructionLayer
)


Q_ = Quantity
dt_sec: int = 3600  # time step in seconds


def set_time_step(dt_hr: float) -> float:
    """Sets the time step for solving the system of node equations in a linear
    thermal network.

    Parameters
    ----------
    dt_hr:
        Time step expressed as a fraction of 1 hour, e.g., `dt_hr` = 1/4 means
        the time step for the calculations is one quarter of an hour.

    Returns
    -------
    The time step converted to seconds.
    """
    global dt_sec
    dt_sec = dt_hr * 3600
    return dt_sec


class ExteriorBuildingElementNode(ABC):
    """Represents a general temperature node in the linear thermal network of an
    exterior building element.

    A temperature node has a unit thermal capacity `C`, expressed per unit of
    area, i.e., in SI units J/(K.m²).
    The node can be connected on two sides (the left and the right side) to
    either another temperature node in the linear network, or to the adjacent
    environment through a thermal resistor. The thermal resistor on the left is
    designated `R1`, while the resistor on the right is designated `R2`. The
    thermal resistance is expressed on the basis of one unit area, i.e., in SI
    units K.m²/W.
    The surface area, `A`, associated with a temperature node is considered
    separately.
    """
    def __init__(self):
        self.ID: str = ''
        self._C: float = 0.0              # J / K
        self._A: float = 1.0              # m**2
        self._R1: float = float('inf')    # K / W
        self._R2: float = float('inf')    # K / W

    def __str__(self) -> str:
        return (
            f"ID: {self.ID}\n"
            f"- A = {Q_(self._A, 'm**2'):~P.3g}\n"
            f"- C = {Q_(self._C, 'J / K'):~P.3g}\n"
            f"- R1 = {Q_(self._R1, 'K / W'):~P.3g}\n"
            f"- R2 = {Q_(self._R2, 'K / W'):~P.3g}"
        )

    @classmethod
    def create(
        cls,
        ID: str,
        C: Quantity = Q_(0.0, 'J / (K * m**2)'),
        A: Quantity = Q_(1.0, 'm**2'),
        R1: Quantity = Q_(float('inf'), 'K * m**2 / W'),
        R2: Quantity = Q_(float('inf'), 'K * m**2 / W')
    ) -> 'ExteriorBuildingElementNode':
        """Creates an `ExteriorBuildingElementNode` object.

        Parameters
        ----------
        ID:
            Name to identify the node in a linear thermal network.
        C:
            The thermal capacity of the node per unit area.
        A:
            The surface area associated with the node.
        R1:
            The unit thermal resistance between the preceding node and this node.
        R2:
            The unit thermal resistance between this node and the next node.
        """
        node = cls()
        node.ID = ID
        node._A = A.to('m**2').magnitude
        node._C = C.to('J / (K * m**2)').magnitude * node._A
        node._R1 = R1.to('K * m**2 / W').magnitude / node._A
        node._R2 = R2.to('K * m**2 / W').magnitude / node._A
        return node

    @abstractmethod
    def get_a_coefficients(self) -> list[float]:
        """Returns the coefficients of the node equation to put into the
        coefficient matrix A of the linear thermal network to be solved.
        """
        ...

    @abstractmethod
    def get_b_value(self, *args, **kwargs) -> float:
        """Returns the input value of the node equation to put in the input
        matrix B of the linear thermal network to be solved.
        """
        ...

    @property
    def A(self) -> Quantity:
        """Returns the surface area associated with the temperature node."""
        return Q_(self._A, 'm**2')

    @property
    def R1(self) -> Quantity | list[Quantity]:
        """Returns the unit thermal resistance R1."""
        if isinstance(self._R1, Sequence):
            return [Q_(r1 * self._A, 'K * m**2 / W') for r1 in self._R1]
        else:
            return Q_(self._R1 * self._A, 'K * m**2 / W')

    @property
    def R2(self) -> Quantity:
        """Returns the unit thermal resistance R2."""
        return Q_(self._R2 * self._A, 'K * m**2 / W')

    @property
    def C(self) -> Quantity:
        """Returns the thermal capacity per unit area C."""
        return Q_(self._C / self._A, 'J / (K * m**2)')


class ExteriorSurfaceNode(ExteriorBuildingElementNode):
    """Represents the first temperature node at the exterior side of an
    exterior building element. It is connected to the exterior temperature
    (sol-air temperature).
    """
    @property
    def _a1(self) -> float:
        """Coefficient in the node equation of the exterior surface node
        temperature.
        """
        a1 = 3 + (2 * dt_sec / self._C) * (1 / self._R1 + 1 / self._R2)
        return a1

    @property
    def _a2(self) -> float:
        """Coefficient in the node equation of the adjacent building mass node
        temperature.
        """
        a2 = -2 * dt_sec / (self._C * self._R2)
        return a2

    def _b(self, T: Sequence[float], T_ext: float) -> float:
        """Calculates the input-side (rhs) of the node equation at the current
        time moment t or time index k (t = k * dt).

        Parameters
        ----------
        T:
            List with the previous node temperature at time index k-1 and
            time index k-2.
        T_ext:
            The exterior sol-air temperature at the current time index k.
        """
        k = 2 * dt_sec / (self._C * self._R1)
        b = 4 * T[0] - T[1] + k * T_ext
        return b

    def get_a_coefficients(self) -> list[float]:
        return [self._a1, self._a2]

    def get_b_value(self, T: Sequence[float], T_ext: float) -> float:
        return self._b(T, T_ext)


class BuildingMassNode(ExteriorBuildingElementNode):
    """Represents a temperature node inside the exterior building element.
    It is connected on the left side with an exterior surface node or another
    building mass node. On the right side it is connected with an interior
    surface node or another building mass node.
    """
    @property
    def _a1(self) -> float:
        """Coefficient in the node equation of the preceding building mass node
        temperature or exterior surface node temperature.
        """
        a1 = -2 * dt_sec / (self._C * self._R1)
        return a1

    @property
    def _a2(self) -> float:
        """Coefficient in the node equation of this building mass node
        temperature.
        """
        a2 = 3 + (2 * dt_sec / self._C) * (1 / self._R1 + 1 / self._R2)
        return a2

    @property
    def _a3(self) -> float:
        """Coefficient in the node equation of the next building mass node
        temperature or interior surface node temperature.
        """
        a3 = -2 * dt_sec / (self._C * self._R2)
        return a3

    @staticmethod
    def _b(T: Sequence[float]) -> float:
        """Calculates the input-side (rhs) of the node equation at the current
        time moment t or time index k (t = k * dt).

        Parameters
        ----------
        T:
            List with the previous node temperature at time index k-1 and
            time index k-2.
        """
        b = 4 * T[0] - T[1]
        return b

    def get_a_coefficients(self) -> list[float]:
        return [self._a1, self._a2, self._a3]

    def get_b_value(self, T: Sequence[float]) -> float:
        return self._b(T)


class InteriorSurfaceNode(ExteriorBuildingElementNode):
    """Represents a temperature node on the interior surface of an exterior
    building element. It has no thermal capacity. On the right side the
    node is connected to the zone air temperature.
    As the interior surface node has no thermal capacity, the heat flow into
    the interior surface node equals the heat flow out to the zone air.
    """
    @property
    def _a1(self) -> float:
        """Coefficient in the node equation of the preceding building mass node
        temperature.
        """
        a1 = 1 / self._R1
        return a1

    @property
    def _a2(self) -> float:
        """Coefficient in the node equation of the internal surface node
        temperature.
        """
        a2 = -(1 / self._R1 + 1 / self._R2)
        return a2

    def get_a_coefficients(self) -> list[float]:
        return [self._a1, self._a2]

    def get_b_value(self, T_zone: float) -> float:
        return -T_zone / self._R2


class ExteriorBuildingElementLTN(list[ExteriorBuildingElementNode]):
    """Represents the linear thermal network model of an exterior building
    element. A `ExteriorBuildingElementLTN` object is a list of
    `ExteriorBuildingElementNode` objects.
    The first element must always be an `ExteriorSurfaceNode` object.
    The last element must always be an `InteriorSurfaceNode` object.
    Between the first and the last element, the list contains at least one or
    more `BuildingMassNode` objects.
    """
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._A: np.ndarray | None = None
        self._B: np.ndarray | None = None
        self._T_node_table: np.ndarray | None = None
        self._Q_dot_table: np.ndarray | None = None
        self._T_ext_series: list[float] | None = None
        self._T_zone_series: list[float] | None = None

    def __str__(self) -> str:
        s = ""
        for node in self:
            s += str(node) + '\n'
        return s

    def _init(
        self,
        init_values: list[list[Quantity]] | None = None,
        cycle: int | None = None
    ) -> None:
        """Initializes the temperature node table at each new calculation
        cycle.

        Parameters
        ---------
        init_values: optional
            2D-list with the initial temperatures (`Quantity` objects) of the
            nodes at time index -2 (first row) and at time index -1 (second
            row). The number of columns of this 2D-list must be equal to the
            number of nodes in the linear thermal network.
            If `init_values` is None, an array of zeros will be created.
        cycle:
            The index of the current cycle.
        """
        self._T_ext_series, self._T_zone_series = [], []
        if cycle is None or cycle == 0:
            if init_values is None:
                n_cols = len(self)
                n_rows = 2
                self._T_node_table = np.zeros((n_rows, n_cols))
            else:
                self._T_node_table = np.array([
                    [T.to('K').m for T in row]
                    for row in init_values
                ])
        else:
            # `cycle` > 0.
            # Initial values for the next cycle are the last two rows in
            # `_T_node_table` from the previous cycle:
            self._T_node_table = self._T_node_table[-2:, :]

    def _build_A_matrix(self) -> None:
        """Builds the coefficient matrix of the linear thermal network."""
        n = len(self)
        i_max = n - 1
        self._A = np.zeros((n, n))
        for i, node in enumerate(self):
            a = node.get_a_coefficients()
            if i == 0:
                # exterior surface node
                self._A[i, i] = a[0]
                self._A[i, i + 1] = a[1]
            elif 0 < i < i_max:
                # building mass node
                self._A[i, i - 1] = a[0]
                self._A[i, i] = a[1]
                self._A[i, i + 1] = a[2]
            elif i == i_max:
                # interior surface node
                self._A[i, i - 1] = a[0]
                self._A[i, i] = a[1]

    def _build_B_matrix(self, T_ext: float, T_zone: float) -> None:
        """Builds the input matrix of the linear thermal network at time moment
        t or time index k.

        Parameters
        ----------
        T_ext:
            Exterior temperature in Kelvin at time moment t or time index k at
            the side of the exterior surface node.
        T_zone:
            Zone air temperature in Kelvin at time moment t or time index k at
            the side of the interior surface node.
        """
        n = len(self)
        i_max = n - 1
        self._B = np.zeros((n,))
        for i, node in enumerate(self):
            if i == 0:
                # exterior surface node
                self._B[i] = node.get_b_value([
                    self._T_node_table[-1, i],
                    self._T_node_table[-2, i]],
                    T_ext
                )
            elif 0 < i < i_max:
                # building mass node
                self._B[i] = node.get_b_value([
                    self._T_node_table[-1, i],
                    self._T_node_table[-2, i]
                ])
            elif i == i_max:
                # interior surface node
                self._B[i] = node.get_b_value(T_zone)

    def _calculate_heat_flows(
        self,
        T_ext: float,
        T_node_row: np.ndarray,
    ) -> np.ndarray:
        """Calculates the heat flow in and out of the linear thermal network
        and between the intermediate nodes in the linear thermal network at
        the current time moment t or time index k.

        Parameters
        ----------
        T_ext: Quantity
            Exterior temperature in Kelvins at time moment t or time index k at
            the side of the exterior surface node.
        T_node_row: array of floats
            Numpy 1D-array of the node temperatures in units of Kelvin at time
            moment t or time index k.
        """
        n = len(self)
        i_max = n - 1
        Q_dot_row = np.zeros((n,))
        for i, node in enumerate(self):
            if i == 0:
                # heat flow into exterior surface node
                Q_dot_row[i] = (T_ext - T_node_row[i]) / node._R1
            elif 0 < i < i_max:
                # heat flow into building mass node
                Q_dot_row[i] = (T_node_row[i - 1] - T_node_row[i]) / node._R1
            elif i == i_max:
                # heat flow into interior surface node = heat flow out of
                # interior surface node (as C = 0)
                Q_dot_row[i] = (T_node_row[i - 1] - T_node_row[i]) / node._R1
        return Q_dot_row

    def solve_one_step(
        self,
        T_ext: float,
        T_zone: float
    ) -> np.ndarray:
        """Solves the linear thermal network for the node temperatures and heat
        flows at the current time moment or time index.

        Parameters
        ----------
        T_ext:
            Exterior temperature in Kelvins at the current time moment or index.
        T_zone:
            Interior zone air temperature in Kelvins at the current time moment
            or index.

        Returns
        -------
        Numpy 1D-array with the heat flows in units of Watts into the nodes at 
        the current time moment or time index. The last element is the heat that
        flows from the interior surface node toward the zone air.
        """
        self._T_ext_series.append(T_ext)
        self._T_zone_series.append(T_zone)
        self._build_B_matrix(T_ext, T_zone)
        T_node_row = np.linalg.solve(self._A, self._B)
        Q_dot_row = self._calculate_heat_flows(T_ext, T_node_row)
        T_node_row = T_node_row.reshape(1, len(T_node_row))
        self._T_node_table = np.append(self._T_node_table, T_node_row, axis=0)
        Q_dot_row = Q_dot_row.reshape(1, len(Q_dot_row))
        return Q_dot_row

    def solve(
        self,
        num_steps: int,
        T_ext: Callable[[float], Quantity],
        T_zone: Callable[[float], Quantity],
        init_values: list[list[Quantity]] | None = None,
        dt_hr: float = 1.0,
        num_cycles: int = 1
    ) -> None:
        """Solves the linear thermal network for the node temperatures in the
        course of time when the exterior temperature `T_ext` as a function of
        time and the zone air temperature `T_zone` as a function of time are
        known.

        Parameters
        ----------
        num_steps:
            The number of time steps in one cycle at which the node temperatures
            are calculated.
        T_ext:
            Function that takes time moment t in seconds from time zero as an
            argument and returns the exterior temperature (`Quantity` object) at
            that time.
        T_zone:
            Function that takes time moment t in seconds from time zero as an
            argument and returns the zone air temperature (`Quantity` object)
            at that time.
        init_values: optional
            2D-list with the initial temperatures (`Quantity` objects) of the
            nodes at time index -2 (first row) and at time index -1 (second
            row) before the first cycle. The number of columns of this 2D-list
            must be equal to the number of nodes in the linear thermal network.
            If `init_values` is None, an array of zeros will be created.
        dt_hr: optional
            The time step width in hours between two successive time moments
            at which the linear thermal network is solved. The default value
            is 1 hr. The product of the number of time steps per cycle and the
            time step width determines the duration (period) of one cycle.
        num_cycles:
            The number of times the cycle of `num_steps` calculations is repeated.
            Each new cycle starts with the node temperatures from the last two
            time indexes k-2 and k-1 of the previous cycle as the initial values.
            Only the last cycle is kept.

        Returns
        -------
        None
        """
        set_time_step(dt_hr)
        self._build_A_matrix()
        for cycle in range(num_cycles):
            self._init(init_values, cycle)
            self._Q_dot_table = np.zeros((num_steps, len(self)))
            for k in range(num_steps):
                t = k * dt_sec
                T_ext_ = T_ext(t).to('K').magnitude
                T_zone_ = T_zone(t).to('K').magnitude
                Q_dot_row = self.solve_one_step(T_ext_, T_zone_)
                self._Q_dot_table[k] = Q_dot_row

    def get_node_temperatures(self, unit: str = 'degC') -> pd.DataFrame:
        """Returns a Pandas DataFrame object of which the column titles are the
        ID's of the nodes in the linear thermal network, ordered from the
        exterior surface node on the left to the interior surface node on the
        right.
        Each row contains the values of the node temperatures calculated at
        each time index from the last calculation cycle (see parameter
        `num_cycles` of method `solve()`).
        """
        d = {
            'EXT': [
                Q_(T_ext, 'K').to(unit).magnitude
                for T_ext in self._T_ext_series
            ]
        }
        for c in range(self._T_node_table.shape[1]):
            node = self[c]
            d[node.ID] = [
                Q_(self._T_node_table[r, c], 'K').to(unit).magnitude
                for r in range(2, self._T_node_table.shape[0])
            ]
        d['ZONE'] = [
            Q_(T_zone, 'K').to(unit).magnitude
            for T_zone in self._T_zone_series
        ]
        df = pd.DataFrame(d)
        return df

    def get_heat_flows(self, unit: str = 'W') -> pd.DataFrame:
        """Returns a Pandas DataFrame object of which the column titles are the
        ID's of the nodes in the linear thermal network from the exterior
        surface node on the left to the interior surface node on the right.
        Each row contains the values of heat flows between the nodes calculated
        at each time index from the last calculation cycle (see parameter
        `num_cycles` of method `solve()`).
        """
        d = {}
        for c in range(self._Q_dot_table.shape[1]):
            node = self[c]
            d[node.ID] = [
                Q_(self._Q_dot_table[r, c], 'W').to(unit).magnitude
                for r in range(self._Q_dot_table.shape[0])
            ]
        df = pd.DataFrame(d)
        return df

    def get_coefficient_matrix(self) -> np.ndarray:
        """Returns the coefficient matrix A of the linear thermal network."""
        if self._A is None:
            self._build_A_matrix()
        return self._A


class ThermalStorageNode:
    """Represents the interior thermal mass in a zone."""
    def __init__(self):
        self.ID: str = ''
        self._A: float = 1.0            # m**2
        self._C: float = 0.0            # J / K
        self._R2: float = float('inf')  # K / W
        self._T_node_list: list[float] = []
        self._Q_dot_list: list[float] = []
        self._c: list[float] = []

    @classmethod
    def create(
        cls,
        ID: str,
        C: Quantity = Q_(0.0, 'J / (K * m**2)'),
        A: Quantity = Q_(1.0, 'm**2'),
        R2: Quantity = Q_(float('inf'), 'K * m**2 / W')
    ) -> 'ThermalStorageNode':
        """Creates a `ThermalStorageNode` object.

        Parameters
        ----------
        ID:
            Name to identify the node in a linear thermal network.
        C:
            The thermal capacity of the node per unit area.
        A:
            The surface area associated with the node.
        R2:
            The unit thermal resistance between this node and the zone air
            node.
        """
        node = cls()
        node.ID = ID
        node._A = A.to('m**2').magnitude
        node._C = C.to('J / (K * m**2)').magnitude * node._A
        node._R2 = R2.to('K * m**2 / W').magnitude / node._A
        return node

    def _calculate_coefficients(self):
        # Determine the coefficients of the node equation only once and keep
        # them in `self._c`:
        a = 3 + 2 * dt_sec / (self._R2 * self._C)
        b_lst = [4, -1, 2 * dt_sec / self._C, 2 * dt_sec / (self._R2 * self._C)]
        self._c = [b / a for b in b_lst]

    def _init(
        self,
        init_values: list[Quantity] | None = None,
        cycle: int | None = None
    ) -> None:
        """Initializes the temperature node list.

        Parameters
        ---------
        init_values:
            List with the 2 initial temperatures (`Quantity` objects) of the
            thermal storage node at time index -2 (first row) and at time index
            -1 (second row). If `init_values` is None, a list of zeros will be
            created.
        cycle:
            The index of the current cycle.
        """
        if cycle is None or cycle == 0:
            # Initialize the node temperature list with the initial values:
            if init_values is None:
                self._T_node_list = [Q_(0, 'K').m] * 2
            else:
                self._T_node_list = [T.to('K').m for T in init_values]
        else:
            # `cycle` > 0.
            # Initial values for the next cycle are the last two elements in
            # `_T_node_list` from the previous cycle:
            self._T_node_list = self._T_node_list[-2:]

    def solve_one_step(
        self,
        T_zone: float,
        Q_dot_rad: float
    ) -> float:
        """Solves the node equation for the node temperature at the current
        time moment or time index.
        """
        T = self._c[0] * self._T_node_list[-1]
        T += self._c[1] * self._T_node_list[-2]
        T += self._c[2] * Q_dot_rad
        T += self._c[3] * T_zone
        self._T_node_list.append(T)
        return T

    def solve(
        self,
        num_steps: int,
        Q_dot_rad: Callable[[float], Quantity],
        T_zone: Callable[[float], Quantity],
        init_values: list[Quantity] | None = None,
        dt_hr: float = 1.0,
        num_cycles: int = 1
    ) -> None:
        """Solves the thermal storage node for the node temperature in the
        course of time when the radiative heat flow input as a function of
        time and the zone air temperature `T_zone` as a function of time are
        known.

        Parameters
        ----------
        num_steps:
            The number of time steps in one cycle at which the node temperatures
            are calculated.
        Q_dot_rad:
            Function that takes time moment t in seconds from time zero as an
            argument and returns the radiative heat flow input (`Quantity`
            object) at that time.
        T_zone:
            Function that takes time moment t in seconds from time zero as an
            argument and returns the zone air temperature (`Quantity` object)
            at that time.
        init_values: optional
            1D-list with the 2 initial temperatures (`Quantity` objects) of the
            node at time index -2 (first element) and at time index -1 (second
            element) before the first cycle starts.
            If `init_values` is None, an list of zeros will be created.
        dt_hr: optional
            The time step width in hours between two successive time moments
            at which the node equation is solved. The default value
            is 1 hr. The product of the number of time steps per cycle and the
            time step width determines the duration (period) of one cycle.
        num_cycles:
            The number of times the cycle of `num_steps` calculations is repeated.
            Each new cycle starts with the node temperatures from the last two
            time indexes k-2 and k-1 of the previous cycle as the initial values.
            Only the last cycle is kept.

        Returns
        -------
        None
        """
        set_time_step(dt_hr)
        self._calculate_coefficients()
        for cycle in range(num_cycles):
            self._init(init_values, cycle)
            self._Q_dot_list = []
            for k in range(num_steps):
                t = k * dt_sec
                Q_dot_rad_ = Q_dot_rad(t).to('W').magnitude
                T_zone_ = T_zone(t).to('K').magnitude
                T_tmn = self.solve_one_step(T_zone_, Q_dot_rad_)
                Q_dot = (T_tmn - T_zone_) / self._R2
                self._Q_dot_list.append(Q_dot)

    def get_node_temperatures(self) -> Quantity:
        """Returns a `Quantity` array containing the values of the thermal
        storage node temperatures calculated at each time index from the last
        calculation cycle (see parameter `num_cycles` of method `solve()`).
        """
        return Q_(self._T_node_list, 'K')

    def get_heat_flows(self) -> Quantity:
        """Returns a `Quantity` array containing the values of the heat
        flow between the thermal storage node and the zone air node calculated
        at each time index from the last calculation cycle (see parameter
        `num_cycles` of method `solve()`).
        """
        return Q_(self._Q_dot_list, 'W')


class ExteriorBuildingElementLTNBuilder:
    """Takes the construction assembly of an exterior building element and
    returns the linear thermal network model of the exterior building element.

    The layers of a construction assembly can be further subdivided into a
    number of small slices. Each slice can be represented as a temperature node
    with an incoming thermal resistance R_in, an outgoing thermal resistance
    R_out, and a thermal capacity C (which can be zero). As such, the list of
    slices is equivalent to a list of nodes and each node in this list can be
    represented as a list [R_in, C, R_out]. Consequently, the linear thermal
    network of an exterior building element is a list with a pattern like:
    [R_in1, C1, R_out1, R_in2, C2, R_out2, ...].
    This list can be further reduced by adding the adjacent thermal resistors
    between every two thermal capacitors, which aren't zero and by omitting the
    zero-valued capacitors.
    Finally, this reduced list can be transformed into an instance of the
    `ExteriorBuildingElementLTN` class. This linear thermal network model of
    the exterior building element can then be solved for the node temperatures
    and the conduction heat flows between the nodes.
    """
    @staticmethod
    def build(
        constr_assem: ConstructionAssembly,
        A_net: Quantity
    ) -> ExteriorBuildingElementLTN | None:
        """Creates an `ExteriorBuildingElementLTN` object from a
        `ConstructionAssembly` object.

        Parameters
        ----------
        constr_assem:
            The construction assembly the exterior building element is made off.
        A_net:
            The net surface area of the exterior building element. This is the
            total surface area of the exterior building element minus the
            surface area of large openings like windows and doors.

        Returns
        -------
        Either an `ExteriorBuildingElementLTN` object or None.
        Returns `None` if the construction assembly has no construction layers
        (`constr_assem.layers is None` returns `True`).
        """
        if constr_assem.layers is not None:
            ltn_lst = ExteriorBuildingElementLTNBuilder._compose(constr_assem)
            ltn_red = ExteriorBuildingElementLTNBuilder._reduce(ltn_lst)
            ltn = ExteriorBuildingElementLTNBuilder._transform(ltn_red, A_net)
            return ltn
        return None

    @staticmethod
    def _compose(constr_assem: ConstructionAssembly) -> list[Quantity]:
        """Creates a list of thermal resistors and capacitors with a common
        pattern like:
        [R_in(1), C(1), R_out(1), R_in(2), C(2), R_out(2), ..., R_in(k), C(k), R_out(k)]
        """
        layers: list[ConstructionLayer] = list(constr_assem.layers.values())
        ltn = []
        # Iterate over all layers except the last one, which is always the
        # interior surface film:
        for layer in layers[:-1]:
            n = layer.num_slices
            R_slice = layer.R / (2 * n)
            C_slice = layer.C / n
            slices = [R_slice, C_slice, R_slice] * n
            ltn.extend(slices)
        # The last node in the linear thermal network is the interior surface
        # node, situated on the interior face of the exterior building element.
        # Its thermal capacity is zero, and its outgoing resistance is the
        # interior surface film resistance:
        ltn.append(Q_(0.0, 'J / (K * m**2)'))
        ltn.append(layers[-1].R)
        return ltn

    @staticmethod
    def _reduce(ltn: list[Quantity]) -> list[Quantity]:
        """Reduces `ltn` by adding resistors in series and omitting
        capacitors with a zero value (except the last C, which is the
        interior surface node).
        """
        # [R_in(1), C(1), R_(1+2), C(2), R_(2+3), ..., R_(k-1,k), C(k), R_out(k)]
        R_dummy = Q_(0, 'K * m**2 / W')
        ltn_red = []
        R = ltn[0]
        for i in range(1, len(ltn)):
            if R_dummy.check(ltn[i].dimensionality):
                # element i is a resistor: add it to R
                R += ltn[i]
            else:
                # element i is a capacitor (only keep C if it is greater than
                # zero or if it is the last C in the list).
                if ltn[i].m > 0.0 or i == len(ltn) - 2:
                    ltn_red.append(R)            # add R to the reduced list
                    ltn_red.append(ltn[i])       # add C to the reduced list
                    R = Q_(0.0, 'K * m**2 / W')  # reset R
        if R.magnitude > 0.0:
            ltn_red.append(R)
        return ltn_red

    @staticmethod
    def _transform(ltn_red: list[Quantity], A: Quantity) -> ExteriorBuildingElementLTN:
        """Transforms `ltn_red` into an `ExteriorBuildingElementLTN` object,
        being the linear thermal network model of the exterior building element.
        """
        n = len(ltn_red)
        esn_lst = ltn_red[:3]
        isn_lst = ltn_red[(n - 3):n]
        esn = ExteriorSurfaceNode.create(
            ID='ESN',
            A=A,
            C=esn_lst[1],
            R1=esn_lst[0],
            R2=esn_lst[2]
        )
        isn = InteriorSurfaceNode.create(
            ID='ISN',
            A=A,
            R1=isn_lst[0],
            R2=isn_lst[2]
        )
        bmn_lst = [
            BuildingMassNode.create(
                ID=f'BMN{k+1}',
                A=A,
                C=ltn_red[i],
                R1=ltn_red[i - 1],
                R2=ltn_red[i + 1]
            ) for k, i in enumerate(range(3, n - 2, 2))
        ]
        ltn = ExteriorBuildingElementLTN()
        ltn.append(esn)
        ltn.extend(bmn_lst)
        ltn.append(isn)
        return ltn
