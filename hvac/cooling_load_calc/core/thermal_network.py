from typing import Callable, TYPE_CHECKING
from abc import ABC, abstractmethod
import numpy as np
from hvac import Quantity
from .fenestration import Window
from .schedule import OnOffSchedule

if TYPE_CHECKING:
    from .building_element import (
        ExteriorBuildingElement,
        InteriorBuildingElement
    )


Q_ = Quantity


class AbstractNode(ABC):

    def __init__(
        self,
        ID: str,
        R_unit_lst: list[Quantity],
        C: Quantity | None = None,
        A: Quantity = Q_(1.0, 'm ** 2'),
        T_input: Callable[[float], float] | None = None,
        Q_input: Callable[[float, float | None], float] | None = None
    ):
        self.ID = ID
        self._A = A.to('m ** 2').m
        self.R_unit = [R.to('m ** 2 * K / W').m for R in R_unit_lst]
        self.R = [R / self._A for R in self.R_unit]
        self.C_unit = self.C = None
        if C is not None:
            self.C_unit = C.to('J / (K * m ** 2)').m 
            self.C = self.C_unit * self._A
        self.T_input = T_input
        self.Q_input = Q_input
        self.dt: float = 0.0

    @abstractmethod
    def get_coefficients(
        self,
        cooling_active: bool = True
    ) -> dict[str, float]:
        """Get coefficients of the LHS of the node equation."""
        ...

    @abstractmethod
    def get_input(
        self,
        k: int,
        T_node_prev: tuple[float, float] | None = None,
        cooling_active: bool = True
    ) -> float:
        """
        Get RHS of the node equation.

        Parameters
        ----------
        k:
            Current time index.
        T_node_prev:
            Tuple of the two previous node temperatures. First element
            (index 0 or -2) at time index k-2 and second element (index 1 or -1)
            at k-1.
        cooling_active:
            Indicate if cooling is ON (True) or OFF (False)
        """
        ...

    @property
    def A(self) -> float:
        return self._A

    @A.setter
    def A(self, v: Quantity) -> None:
        self._A = v.to('m ** 2').m
        self.R = [R / self._A for R in self.R_unit]
        if self.C_unit is not None: self.C = self.C_unit * self._A

    def __str__(self):
        l1 = f"Node {self.ID}:\n"
        l2 = "\tR-links: " + str([round(R, 2) for R in self.R]) + "\n"
        try:
            l3 = f"\tC: {self.C:.2f}\n"
        except TypeError:
            l3 = "\tC: None\n"
        return l1 + l2 + l3


class BuildingMassNode(AbstractNode):
    """
    A building mass node is situated within the mass of the construction assembly and
    is connected with a preceding node and a next node.

    Parameter `R_list` is a list of two resistances in the given order:
    1. the conduction resistance between this node and the preceding node, and
    2. the conduction resistance between this node and the next node.
    """
    def get_coefficients(
        self,
        cooling_active: bool = True
    ) -> dict[str, float]:
        return {
            'i, i-1': -2 * self.dt / (self.R[0] * self.C),
            'i, i': 3 + (2 * self.dt / self.C) * sum(1/R for R in self.R),
            'i, i+1': -2 * self.dt / (self.R[1] * self.C)
        }

    def get_input(
        self,
        k: int,
        T_node_prev: tuple[float, float] | None = None,
        cooling_active: bool = True
    ) -> float:
        return 4 * T_node_prev[-1] - T_node_prev[-2]


class ExteriorSurfaceNode(AbstractNode):
    """
    An exterior surface node is the first node of the 1D thermal network at the
    exterior side of the construction assembly that is connected to the sol-air
    temperature.

    Parameter `T_input` at instantiation of the node must be assigned a function
    that returns the sol-air temperature at each time moment k.

    Parameter `R_list` is a list of two resistances in the given order:
    1. The first resistance links the node to the sol-air temperature.
    2. The second resistance links the node to the next, internal node
    (which is an instance of class `BuildingMassNode`).
    """
    def get_coefficients(
        self,
        cooling_active: bool = True
    ) -> dict[str, float]:
        return {
            'i, i': 3 + (2 * self.dt / self.C) * sum(1/R for R in self.R),
            'i, i+1': -2 * self.dt / (self.R[1] * self.C)
        }

    def get_input(
        self,
        k: int,
        T_node_prev: tuple[float, float] | None = None,
        cooling_active: bool = True
    ) -> float:
        return (
            4 * T_node_prev[-1] - T_node_prev[-2]
            + 2 * self.dt / (self.C * self.R[0]) * self.T_input(k * self.dt)
        )


class InteriorSurfaceNode(AbstractNode):
    """
    An interior surface node is the last node of the 1-D thermal network of a
    construction assembly. It is situated on the surface of the construction
    assembly and as such it has no thermal capacity.

    The interior surface node can be connected to the zone air node through a
    convective resistance and/or to a thermal storage mass node within the zone
    through a radiative resistance.

    Parameter `T_input` at instantiation of the node must be assigned a function
    that returns the zone air temperature at each time moment *t* in units of
    seconds.

    Parameter `R_list` is a list of 2 or 3 resistances in the given order:
    1. the conductive resistance between the interior surface node and the
    preceding building mass node,
    2. the convective resistance between the interior surface node and the zone
    air node,
    3. the optionally radiative resistance between the interior surface node and
    the thermal storage node (if present).
    """
    def _get_coefficients_Tz_known(self) -> dict[str, float]:
        if len(self.R) == 3:
            # True if there is also a ThermalStorageNode present
            return {
                'i, i-1': 1 / self.R[0],
                'i, i': -sum(1/R for R in self.R),
                'i, -1': 1 / self.R[-1]
            }
        else:
            return {
                'i, i-1': 1 / self.R[0],
                'i, i': -sum(1/R for R in self.R)
            }

    def _get_input_Tz_known(
        self,
        k: int,
    ) -> float:
        return (-1 / self.R[1]) * self.T_input(k * self.dt)

    def _get_coefficients_Tz_unknown(self) -> dict[str, float]:
        if len(self.R) == 3:
            # True if there is also a ThermalStorageNode present
            return {
                'i, i-1': 1 / self.R[0],
                'i, i': -sum(1/R for R in self.R),
                'i, -2': 1 / self.R[-2],
                'i, -1': 1 / self.R[-1]
            }
        else:
            return {
                'i, i-1': 1 / self.R[0],
                'i, i': -sum(1/R for R in self.R),
                'i, -2': 1 / self.R[-1]
            }

    # noinspection PyMethodMayBeStatic
    def _get_input_Tz_unknown(self) -> float:
        return 0.0

    def get_coefficients(
        self,
        cooling_active: bool = True
    ) -> dict[str, float]:
        if cooling_active is True:
            return self._get_coefficients_Tz_known()
        else:
            return self._get_coefficients_Tz_unknown()

    def get_input(
        self,
        k: int,
        T_node_prev: tuple[float, float] | None = None,
        cooling_active: bool = True
    ) -> float:
        if cooling_active is True:
            return self._get_input_Tz_known(k)
        else:
            return self._get_input_Tz_unknown()


class ThermalStorageNode(AbstractNode):
    def __init__(
        self,
        ID: str,
        R_unit_lst: list[Quantity],
        C: Quantity | None = None,
        A: Quantity = Q_(1.0, 'm ** 2'),
        T_input: Callable[[float], float] | None = None,
        Q_input: Callable[[float], float] | None = None,
        windows: list[Window] | None = None,
        int_building_elements: list['InteriorBuildingElement'] | None = None,
        ext_doors: list['ExteriorBuildingElement'] | None = None,
        int_doors: list['InteriorBuildingElement'] | None = None
    ):
        super().__init__(ID, R_unit_lst, C, A, T_input, Q_input)
        if windows:
            self.UA_wnd = sum(
                wnd.F_rad * wnd.UA
                for wnd in windows
            )
            self.T_ext = windows[0].T_ext
        else:
            self.UA_wnd = 0.0
            self.T_ext = lambda t: 0.0
        if int_building_elements:
            self.int_building_elements = int_building_elements
            self.UA_ibe = sum(
                ibe.F_rad * ibe.UA
                for ibe in int_building_elements
            )
        else:
            self.int_building_elements = []
            self.UA_ibe = 0.0
        if ext_doors:
            self.UA_edr = sum(
                edr.F_rad * edr.UA
                for edr in ext_doors
            )
            self.T_sol = ext_doors[0].T_sol
        else:
            self.UA_edr = 0.0
        if int_doors:
            self.int_doors = int_doors
            self.UA_idr = sum(
                idr.F_rad * idr.UA
                for idr in int_doors
            )
        else:
            self.int_doors = []
            self.UA_idr = 0.0

    def _get_coefficients_Tz_known(self) -> dict[str, float]:
        a = {f'-1, {i + 1}': -2 * self.dt / (self.C * self.R[i]) for i in range(len(self.R) - 1)}
        a['-1, -1'] = 3 + (2 * self.dt / self.C) * sum(1/R for R in self.R)
        return a

    def _get_input_Tz_known(
        self,
        k: int,
        T_node_prev: tuple[float, float] | None = None,
    ) -> float:
        b = (
            4 * T_node_prev[-1] - T_node_prev[-2]
            + (2 * self.dt / (self.C * self.R[-1])) * self.T_input(k * self.dt)
        )
        if self.Q_input is None:
            return b
        else:
            return b + 2 * self.dt * self.Q_input(k * self.dt, None) / self.C

    def _get_coefficients_Tz_unknown(self) -> dict[str, float]:
        a = {f'-1, {i + 1}': -2 * self.dt / (self.C * self.R[i]) for i in range(len(self.R) - 1)}
        a['-1, -1'] = 3 + (2 * self.dt / self.C) * sum(1/R for R in self.R)
        a['-1, -2'] = (
            -2 * self.dt / (self.C * self.R[-1]) +
            2 * self.dt / self.C * (
                self.UA_wnd +
                self.UA_ibe +
                self.UA_edr +
                self.UA_idr
            )
        )
        return a

    def _get_input_Tz_unknown(
        self,
        k: int,
        T_node_prev: tuple[float, float]
    ) -> float:
        b = (
            4 * T_node_prev[-1] - T_node_prev[-2]
            + 2 * self.dt / self.C * (
                self.UA_wnd * self.T_ext(k * self.dt) +
                sum(ibe.UA * ibe.T_adj(k * self.dt) for ibe in self.int_building_elements) +
                self.UA_edr * self.T_sol(k * self.dt) +
                sum(idr.UA * idr.T_adj(k * self.dt) for idr in self.int_doors)
            )
        )
        return b

    def get_coefficients(
        self,
        cooling_active: bool = True
    ) -> dict[str, float]:
        if cooling_active is True:
            return self._get_coefficients_Tz_known()
        else:
            return self._get_coefficients_Tz_unknown()

    def get_input(
        self,
        k: int,
        T_node_prev: tuple[float, float] | None = None,
        cooling_active: bool = True
    ) -> float:
        if cooling_active is True:
            return self._get_input_Tz_known(k, T_node_prev)
        else:
            return self._get_input_Tz_unknown(k, T_node_prev)


class ZoneAirNode(AbstractNode):

    def __init__(
        self,
        ID: str,
        R_unit_lst: list[Quantity],
        A: Quantity = Q_(1, 'm ** 2'),
        Q_input: Callable[[float], float] | None = None,
        windows: list[Window] | None = None,
        int_building_elements: list['InteriorBuildingElement'] | None = None,
        ext_doors: list['ExteriorBuildingElement'] | None = None,
        int_doors: list['InteriorBuildingElement'] | None = None
    ):
        super().__init__(ID, R_unit_lst, None, A, None, Q_input)
        if windows:
            self.UA_wnd = sum(
                (1.0 - wnd.F_rad) * wnd.UA
                for wnd in windows
            )
            self.T_ext = windows[0].T_ext
        else:
            self.UA_wnd = 0.0
            self.T_ext = lambda t: 0.0
        if int_building_elements:
            self.int_building_elements = int_building_elements
            self.UA_ibe = sum(
                (1.0 - ibe.F_rad) * ibe.UA
                for ibe in int_building_elements
            )
        else:
            self.int_building_elements = []
            self.UA_ibe = 0.0
        if ext_doors:
            self.UA_edr = sum(
                (1.0 - edr.F_rad) * edr.UA
                for edr in ext_doors
            )
            self.T_sol = ext_doors[0].T_sol
        else:
            self.UA_edr = 0.0
        if int_doors:
            self.int_doors = int_doors
            self.UA_idr = sum(
                (1.0 - idr.F_rad) * idr.UA
                for idr in int_doors
            )
        else:
            self.int_doors = []
            self.UA_idr = 0.0

    def get_coefficients(
        self,
        cooling_active: bool = True
    ) -> dict[str, float]:
        a = {f'-2, {i + 1}': 1 / self.R[i] for i in range(len(self.R) - 1)}
        a['-2, -2'] = -(
            sum(1 / R for R in self.R) +
            self.UA_wnd +
            self.UA_ibe +
            self.UA_edr +
            self.UA_idr
        )
        a['-2, -1'] = 1 / self.R[-1]
        return a

    def get_input(
        self,
        k: int,
        T_node_prev: tuple[float, float] | None = None,
        cooling_active: bool = True
    ) -> float:
        b = (
            -self.UA_wnd * self.T_ext(k * self.dt)
            - sum(ibe.UA * ibe.T_adj(k * self.dt) for ibe in self.int_building_elements)
            - self.UA_edr * self.T_sol(k * self.dt)
            - sum(idr.UA * idr.T_adj(k * self.dt) for idr in self.int_doors)
        )
        if self.Q_input is not None:
            return b - self.Q_input(k * self.dt, None)
        else:
            return b


class ThermalNetwork:
    """Represents the linear thermal network model of a construction assembly."""

    def __init__(
        self,
        nodes: list[AbstractNode],
        int_surf_node_indices: list[int] | None = None
    ) -> None:
        self.nodes = nodes
        self.int_surf_node_indices = int_surf_node_indices
        self._T_ext: Callable[[float], float] | None = None
        self._T_int: Callable[[float], float] | None = None
        self._T_node_table: list[list[float]] | None = None

    @property
    def T_ext(self) -> Callable[[float], float] | None:
        return self._T_ext

    @T_ext.setter
    def T_ext(self, fun: Callable[[float], float]) -> None:
        """Set the temperature input function at the exterior side."""
        self._T_ext = fun
        self.nodes[0].T_input = self._T_ext

    @property
    def T_int(self) -> Callable[[float], float] | None:
        return self._T_int

    @T_int.setter
    def T_int(self, fun: Callable[[float], float]) -> None:
        """Set the temperature input function at the interior side."""
        self._T_int = fun
        self.nodes[-1].T_input = self._T_int

    @property
    def T_node_table(self) -> list[list[Quantity]]:
        tbl = [[Q_(T, 'degC') for T in row] for row in self._T_node_table]
        return tbl

    @T_node_table.setter
    def T_node_table(self, tbl: list[list[float]]) -> None:
        self._T_node_table = tbl

    @property
    def R_ext(self) -> Quantity:
        """Get the thermal resistance between the outdoor environment and
        exterior surface node."""
        R_ext = self.nodes[0].R[0]
        return Q_(R_ext, 'K / W')

    @property
    def R_int(self) -> Quantity:
        """Get the thermal resistance between the interior surface node and
        the indoor environment."""
        R_int = self.nodes[-1].R[-1]
        return Q_(R_int, 'K / W')

    def __str__(self):
        _str = ''
        for node_str in (str(node) for node in self.nodes):
            _str += node_str
        return _str


class ThermalNetworkSolver:
    """
    Solve a linear thermal network model. The solve method can be applied
    to the thermal network model of construction assembly, but also to the
    thermal network of a space.
    """
    dt: float = 0.0
    k_max: int = 0
    Tn_table: list[list[float]] = None

    class _CoolingOn:
        nodes: list[AbstractNode] = None
        n: int = 0
        dt: float = 0.0
        A: np.ndarray = None
        B: np.ndarray = None

        @classmethod
        def init(cls, nodes: list[AbstractNode]):
            cls.nodes = nodes
            cls.n = len(nodes)
            cls.dt = ThermalNetworkSolver.dt
            cls.B = np.zeros((cls.n, 1))
            cls.A = cls._build_matrix_A()

        @classmethod
        def _build_matrix_A(cls) -> np.ndarray:
            """
            Build the coefficient matrix A from the thermal network model.
            """
            A = np.zeros((cls.n, cls.n))
            int_surf_node_indexes = []
            for i, node in enumerate(cls.nodes):
                node.dt = cls.dt
                a = node.get_coefficients(cooling_active=True)
                if isinstance(node, ExteriorSurfaceNode):
                    A[i, i] = a['i, i']
                    A[i, i + 1] = a['i, i+1']
                elif isinstance(node, BuildingMassNode):
                    A[i, i - 1] = a['i, i-1']
                    A[i, i] = a['i, i']
                    A[i, i + 1] = a['i, i+1']
                elif isinstance(node, InteriorSurfaceNode):
                    A[i, i - 1] = a['i, i-1']
                    A[i, i] = a['i, i']
                    try:
                        A[i, -1] = a['i, -1']
                    except KeyError:
                        # the internal surface node is not connected to a
                        # thermal storage node
                        pass
                    else:
                        int_surf_node_indexes.append(i)
                elif isinstance(node, ThermalStorageNode):
                    A[i, i] = a['-1, -1']
                    for p, j in enumerate(int_surf_node_indexes):
                        A[i, j] = a[f'-1, {p + 1}']
                else:
                    pass
            return A

        @classmethod
        def _update_matrix_B(
            cls,
            k: int,
            Tn_prev: tuple[list[float], list[float]],
            cooling_active: bool
        ) -> None:
            """
            Update the input matrix B at time index k.
            """
            for i, node in enumerate(cls.nodes):
                cls.B[i] = node.get_input(
                    k,
                    (Tn_prev[-1][i], Tn_prev[-2][i]),
                    cooling_active
                )

        @classmethod
        def solve(
            cls,
            k: int,
            Tn_prev: tuple[list[float], list[float]],
        ):
            cls._update_matrix_B(k, Tn_prev, cooling_active=True)
            Tn_array = np.linalg.solve(cls.A, cls.B)
            Tn_list = np.transpose(Tn_array)[0].tolist()
            return Tn_list

    class _CoolingOff(_CoolingOn):

        @classmethod
        def _build_matrix_A(cls) -> np.ndarray:
            A = np.zeros((cls.n, cls.n))
            int_surf_node_indexes = []
            for i, node in enumerate(cls.nodes):
                node.dt = cls.dt
                a = node.get_coefficients(cooling_active=False)
                if isinstance(node, ExteriorSurfaceNode):
                    A[i, i] = a['i, i']
                    A[i, i + 1] = a['i, i+1']
                elif isinstance(node, BuildingMassNode):
                    A[i, i - 1] = a['i, i-1']
                    A[i, i] = a['i, i']
                    A[i, i + 1] = a['i, i+1']
                elif isinstance(node, InteriorSurfaceNode):
                    A[i, i - 1] = a['i, i-1']
                    A[i, i] = a['i, i']
                    A[i, -2] = a['i, -2']  # zone air node column
                    A[i, -1] = a['i, -1']  # thermal storage node column
                    int_surf_node_indexes.append(i)
                elif isinstance(node, ZoneAirNode):
                    for p, j in enumerate(int_surf_node_indexes):
                        A[-2, j] = a[f'-2, {p + 1}']
                    A[-2, -2] = a['-2, -2']
                    A[-2, -1] = a['-2, -1']
                elif isinstance(node, ThermalStorageNode):
                    for p, j in enumerate(int_surf_node_indexes):
                        A[-1, j] = a[f'-1, {p + 1}']
                    A[-1, -2] = a['-1, -2']
                    A[-1, -1] = a['-1, -1']
                else:
                    pass
            return A

        @classmethod
        def solve(
            cls,
            k: int,
            Tn_prev: tuple[list[float], list[float]],
            cooling_active: bool = False
        ):
            cls._update_matrix_B(k, Tn_prev, cooling_active=False)
            Tn_array = np.linalg.solve(cls.A, cls.B)
            Tn_list = np.transpose(Tn_array)[0].tolist()
            return Tn_list

    @classmethod
    def _init(cls, n: int, dt_hr: float):
        cls.dt = dt_hr * 3600
        cls.k_max = int(24.0 / dt_hr)
        cls.Tn_table = [
            [0.0] * n,
            [0.0] * n
        ]

    @classmethod
    def solve(
        cls,
        thermal_networks: tuple[ThermalNetwork, ThermalNetwork | None],
        cooling_schedule: OnOffSchedule | None = None,
        dt_hr: float = 1.0,
        n_cycles: int = 6
    ) -> ThermalNetwork:
        if cooling_schedule is None:
            thermal_network = thermal_networks[0]
            cls._init(len(thermal_network.nodes), dt_hr)
            cls._CoolingOn.init(thermal_network.nodes)
            for i in range(n_cycles):
                for k in range(cls.k_max):
                    Tn_list = cls._CoolingOn.solve(
                        k,
                        (cls.Tn_table[-1], cls.Tn_table[-2])
                    )
                    cls.Tn_table.append(Tn_list)
                if i < n_cycles - 1:
                    cls.Tn_table = cls.Tn_table[-2:]
            cls.Tn_table = cls.Tn_table[2:]
            thermal_network.T_node_table = cls.Tn_table
            return thermal_network
        else:
            cls._init(len(thermal_networks[1].nodes), dt_hr)
            cls._CoolingOn.init(thermal_networks[0].nodes)   # cooling ON: T_int = setpoint temperature
            cls._CoolingOff.init(thermal_networks[1].nodes)  # cooling OFF: T_int is floating
            for i in range(n_cycles):
                for k in range(cls.k_max):
                    if cooling_schedule(k * cls.dt):  # cooling ON
                        Tn_list = cls._CoolingOn.solve(
                            k,
                            (cls.Tn_table[-1], cls.Tn_table[-2]),
                        )
                        tmn = thermal_networks[0].nodes[-1]
                        Tz = tmn.T_input(k * cls.dt)
                        Tn_list.insert(-1, Tz)  # insert zone air temperature in temperature node table
                    else:  # cooling OFF
                        Tn_list = cls._CoolingOff.solve(
                            k,
                            (cls.Tn_table[-1], cls.Tn_table[-2]),
                        )
                    cls.Tn_table.append(Tn_list)
                if i < n_cycles - 1:
                    cls.Tn_table = cls.Tn_table[-2:]
            cls.Tn_table = cls.Tn_table[2:]
            thermal_networks[1].T_node_table = cls.Tn_table
            return thermal_networks[1]
