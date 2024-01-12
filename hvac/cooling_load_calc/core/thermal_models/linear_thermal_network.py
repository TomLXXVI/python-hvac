"""LINEAR THERMAL NETWORK OF AN OPAQUE EXTERIOR BUILDING ELEMENT.

Implementation of a linear thermal network (lumped capacitance model) using
`scipy.solve_ivp` to solve the system of differential node equations.
"""
from collections.abc import Callable, Sequence
import numpy as np
from scipy.integrate import solve_ivp
import pandas as pd
from hvac import Quantity


Q_ = Quantity


class TemperatureNode:
    """Represents a temperature node in a linear thermal network.

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

    If the temperature node is connected on one or on both sides to an
    environment, the known environmental temperature as a function of time can
    be connected to the other end of the thermal resistor present on that side
    of the temperature node: on the left side through parameter `T_fun1`, and
    on the right side through parameter `T_fun2`. To return the known
    temperature at a given time moment `t`, these functions need to be called
    by passing this time moment `t` to them. Time moment `t` is the time
    measured from a certain zero point in time, and must be expressed in seconds.

    Instead of a thermal resistor in combination with a temperature function, it
    is also possible to connect a function to the node that directly returns the
    heat flow into or out of the node at a given time moment `t`; the parameter
    `Q_dot_fun1` on the left side of the node is considered to be heat input,
    while the parameter `Q_dot_fun2` is considered to be heat output (so heat
    flows from left to right through the linear thermal network).
    """
    def __init__(self):
        """Creates an 'empty' `TemperatureNode` object."""
        self.ID: str = ''
        self._R1: float = float('inf')
        self._R2: float = float('inf')
        self._C: float = 0.0
        self._A: float = 1.0
        self.T_fun1: Callable[[float], Quantity] | None = None
        self.T_fun2: Callable[[float], Quantity] | None = None
        self.Q_dot_fun1: Callable[[float], Quantity] | None = None
        self.Q_dot_fun2: Callable[[float], Quantity] | None = None

    @classmethod
    def create(
        cls,
        ID: str = '',
        R1: Quantity = Q_(float('inf'), 'K * m**2 / W'),
        R2: Quantity = Q_(float('inf'), 'K * m**2 / W'),
        C: Quantity = Q_(0.0, 'J / (K * m**2)'),
        A: Quantity = Q_(1.0, 'm**2'),
        T_fun1: Callable[[float], Quantity] | None = None,
        T_fun2: Callable[[float], Quantity] | None = None,
        Q_dot_fun1: Callable[[float], Quantity] | None = None,
        Q_dot_fun2: Callable[[float], Quantity] | None = None
    ) -> 'TemperatureNode':
        """Creates a configured `TemperatureNode` object.

        Parameters
        ----------
        ID: str
            Identifier for the node in the linear thermal network.
        R1: Quantity
            Unit thermal resistance between the preceding node and this node in
            the linear thermal network.
        R2: Quantity
            Unit thermal resistance between this node and the next node in the
            linear thermal network.
        C: Quantity
            Unit thermal capacitance of this node.
        A: Quantity
            Surface area associated with this node.
        T_fun1: optional, default None
            A function of time that returns a temperature value (`Quantity`
            object) at a given time `t` (float) present on the left-hand side of
            this node. Time variable `t` should be of type float, being the time
            measured in seconds since an arbitrary initial time moment that is
            set to zero.
        T_fun2: optional, default None
            Same as `T_fun1`, but present on the right-hand side of this node.
        Q_dot_fun1: optional, default None
            A function of time, similar to `T_fun1`, that returns a thermal
            power value (`Quantity` object) at a given time `t` in seconds
            flowing into this node on the left-hand side.
        Q_dot_fun2: optional, default None
            A function of time, similar to `Q_fun1` that returns a thermal power
            value at a given time `t` flowing out of this node on the right-hand
            side.

        Returns
        -------
        TemperatureNode

        Notes
        -----
        If a temperature node is connected on both sides to another temperature
        node, the functions `T_fun1`, `T_fun2`, `Q_dot_fun1`, and `Q_dot_fun2`
        should all be `None`.

        If `Q_dot_fun1` is given, then any value assigned to thermal resistor
        `R1` will be ignored. The same applies to `R2` when `Q_dot_fun2` is
        given.
        """
        node = cls()
        node.ID = ID
        node._R1 = R1.to('K * m**2 / W').magnitude
        node._R2 = R2.to('K * m**2 / W').magnitude
        node._C = C.to('J / (K * m**2)').magnitude
        node._A = A.to('m**2').magnitude
        node.T_fun1 = T_fun1
        node.T_fun2 = T_fun2
        node.Q_dot_fun1 = Q_dot_fun1
        node.Q_dot_fun2 = Q_dot_fun2
        return node

    @property
    def R1(self) -> Quantity:
        """Returns the unit thermal resistance on the left side of this node."""
        return Q_(self._R1, 'K * m**2 / W')

    @property
    def R2(self) -> Quantity:
        """Returns the unit thermal resistance on the right side of this node."""
        return Q_(self._R2, 'K * m**2 / W')

    @property
    def C(self) -> Quantity:
        """Returns the unit thermal capacitance of this node."""
        return Q_(self._C, 'J / (K * m**2)')

    @property
    def A(self) -> Quantity:
        """Returns the surface area associated with this node."""
        return Q_(self._A, 'm**2')

    def Q_dot_1(self, t: float, T: float, T1: float | None = None) -> float:
        """Returns the thermal power input in units of Watt on the left
        side of this node at the time moment `t` measured in seconds from time
        zero. Parameter `T` is the node temperature in units of Kelvin of this
        node at time `t`. Parameter `T1` is the temperature in units of Kelvin
        at the same time `t` of the preceding node to the left from this node in
        a linear thermal network, if it is present, otherwise `T1` should be set
        to `None`.
        """
        if self.Q_dot_fun1 is not None:
            return self.Q_dot_fun1(t).to('W').magnitude
        elif self.T_fun1 is not None:
            T1 = self.T_fun1(t).to('K').magnitude
            R1 = self._R1 / self._A
            Q_dot_1 = (T1 - T) / R1
            return Q_dot_1
        elif T1 is not None:
            R1 = self._R1 / self._A
            Q_dot_1 = (T1 - T) / R1
            return Q_dot_1
        else:
            return 0.0

    def Q_dot_2(self, t: float, T: float, T2: float | None = None) -> float:
        """Returns the thermal power output in units of Watt on the right
        side of this node at the time moment `t` measured in seconds from time
        zero. Parameter `T` is the node temperature in units of Kelvin of this
        node at time `t`. Parameter `T2` is the temperature in units of Kelvin
        at the same time `t` of the next node to the right of this node in a
        linear thermal network, if it is present, otherwise `T2` should be set
        to `None`.
        """
        if self.Q_dot_fun2 is not None:
            return self.Q_dot_fun2(t).to('W').magnitude
        elif self.T_fun2 is not None:
            T2 = self.T_fun2(t).to('K').magnitude
            R2 = self._R2 / self._A
            Q_dot_2 = (T - T2) / R2
            return Q_dot_2
        elif T2 is not None:
            R2 = self._R2 / self._A
            Q_dot_2 = (T - T2) / R2
            return Q_dot_2
        else:
            return 0.0

    def node_equation(
        self,
        t: float,
        T: float,
        T1: float | None = None,
        T2: float | None = None
    ) -> float:
        """Returns the solution of the right-hand side of the heat balance
        equation:
            dT/dt = (1/C) * (Q_dot_1 - Q_dot_2)
        with dT/dt the change per unit time of the node temperature at the
        time moment `t` measured in seconds from time zero.
        The right-hand side of the equation represents the difference between
        the heat flowing into the node and the heat flowing out of the node at
        the time moment `t`.

        Parameters
        ----------
        t:
            The time moment measured in seconds from time zero for which the
            rhs of the heat balance equation of the node is solved.
        T:
            The temperature in units of Kelvin of this node at time moment `t`.
        T1:
            The temperature in units of Kelvin of the node to the left of this
            node at time moment `t`. If no node is present to the left of this
            node, `T1` should be set to `None`.
        T2:
            The temperature in units of Kelvin of the node to the right of this
            node at time moment `t`. If no node is present to the right of this
            node, `T2` should be set to `None`.
        """
        C = self._C * self._A
        out = (1 / C) * (self.Q_dot_1(t, T, T1) - self.Q_dot_2(t, T, T2))
        return out

    def solve(
        self,
        time_span: Sequence[float],
        T_init: Quantity,
        time_eval: Sequence[float] | None = None,
        continuous: bool = False
    ) -> tuple[list[float], list[Quantity]]:
        """Solves the heat balance equation (ODE) for the node temperature as a
         function of time using Scipy's function `solve_ivp()`.

        Parameters
        ----------
        time_span: 2-member sequence
            The initial and final moment between which the node temperatures are
            solved measured in seconds from time zero.
        T_init: Quantity
            The node temperature at the initial moment `t` = 0 s, i.e., time
            zero.
        time_eval: Sequence
            Times at which the node temperatures are returned, sorted between
            the initial and the final time moment, expressed in seconds from
            time zero.
        continuous: bool, default False
            Whether to compute a continuous solution.

        Returns
        -------
        A 2-member tuple containing:
        - A list of the time moments at which the node temperature is
          calculated.
        - A list with the corresponding node temperatures.
        """
        T_init = T_init.to('K').magnitude
        sol = solve_ivp(
            self.node_equation, time_span, [T_init],
            t_eval=time_eval, dense_output=continuous
        )
        return sol.t.tolist(), [Q_(T, 'K') for T in sol.y[0]]


class LinearThermalNetwork(list[TemperatureNode]):
    """Represents a linear thermal network as a list of temperature nodes
    ordered from left to right.

    Once a temperature node has been created, it can be appended to a linear
    thermal network, or multiple temperature nodes, ordered from left to right
    that constitute a linear thermal network, can be passed in an iterable to
    the `__init__()` method of class `LinearThermalNetwork`.

    Attributes
    ----------
    time_axis: list[float]
        List of time moments at which the node temperatures are calculated,
        measured in seconds from time zero.
    T_node_table: list[list[Quantity]]
        List with for each node in the linear thermal network, going from left
        to right, a list with the node temperatures at each time moment for
        which they were calculated (see `time_axis`). So, the number of rows
        (lists) is equal to the number of nodes and the number of columns
        (number of elements in each list) is equal to the number of time moments
        in list `time_axis`.
    Q_dot_table: list[list[Quantity]]
        List with the lists of heat flows between the nodes in the linear
        thermal network at the time moments in list `time_axis`. The number of
        rows (lists) is equal to the number of nodes in the linear thermal
        network plus one. The number of columns (number of elements in each
        list) is equal to the number of time moments in list `time_axis`.
        The first row is the thermal power that flows from the adjacent
        environment on the left into the first, left-most node in the linear
        thermal network. The last row is the thermal power that flows out of the
        last, right-most node in the linear thermal network and into the
        adjacent environment on the right.
    """
    def __init__(self, *args, **kwargs):
        """Creates a `LinearThermalNetwork` object, being a list of
        `TemperatureNode` objects.
        """
        super().__init__(*args, **kwargs)
        self.time_axis: list[float] = []
        self.T_node_table: list[list[Quantity]] = []
        self.Q_dot_table: list[list[Quantity]] = []

    def _network_equation(
        self,
        t: float,
        T_node_seq: Sequence[float]
    ) -> list[float]:
        """Solves the rhs of the heat balance equation of each node in the
        linear thermal network at time moment `t`, measured in seconds from
        time zero.

        Parameters
        ----------
        t:
            Time moment measured in seconds from time zero at which the rhs of
            the heat balance equation of each node in the linear thermal network
            is solved.
        T_node_seq:
            List of the node temperatures in units of Kelvin at time moment `t`.

        Returns
        -------
        A list with the solutions at time moment `t` of the rhs of the heat
        balance equations of the nodes in the linear thermal network.
        """
        out = []
        i_max = len(T_node_seq)
        for i in range(i_max):
            if i == 0:
                T = T_node_seq[i]
                T_next = T_node_seq[i+1]
                out.append(self[i].node_equation(t, T, T2=T_next))
            elif 0 < i < i_max-1:
                T_prev = T_node_seq[i-1]
                T = T_node_seq[i]
                T_next = T_node_seq[i+1]
                out.append(self[i].node_equation(t, T, T1=T_prev, T2=T_next))
            else:  # i == i_max - 1
                T_prev = T_node_seq[i-1]
                T = T_node_seq[i]
                out.append(self[i].node_equation(t, T, T1=T_prev))
        return out

    def _calculate_heat_flows(
        self,
        time_axis: np.ndarray,
        T_node_table: np.ndarray
    ) -> np.ndarray:
        """Calculates the heat flow into and out of the linear thermal network
        and also between the intermediate nodes in the linear thermal network
        for each time moment in the 1D-array `time_axis` after the node
        temperatures have been calculated.
        """
        rows = T_node_table.shape[0] + 1
        cols = T_node_table.shape[1]
        Q_dot_table = np.zeros((rows, cols))
        for j, t in enumerate(time_axis):
            for i in range(len(self)):
                if i == 0:
                    T = T_node_table[i, j]
                    Q_dot_table[i, j] = self[i].Q_dot_1(t, T)
                elif 0 < i < len(self) - 1:
                    T_prev = T_node_table[i-1, j]
                    T = T_node_table[i, j]
                    Q_dot_table[i, j] = self[i].Q_dot_1(t, T, T_prev)
                else:
                    T_prev = T_node_table[i-1, j]
                    T = T_node_table[i, j]
                    Q_dot_table[i, j] = self[i].Q_dot_1(t, T, T_prev)
                    Q_dot_table[i+1, j] = self[i].Q_dot_2(t, T)
        return Q_dot_table

    def solve(
        self,
        time_span: Sequence,
        T_node_init_seq: Sequence[Quantity],
        time_eval: Sequence | None = None,
        continuous: bool = False
    ) -> tuple[list[float], list[list[Quantity]], list[list[Quantity]]]:
        """Solves the linear thermal network for the node temperatures in the
        course of time.

        Parameters
        ----------
        time_span: 2-member sequence
            The initial and final moment between which the node temperatures are
            solved, measured in seconds from time zero.
        T_node_init_seq: Sequence of Quantity
            The node temperatures at time zero, ordered from the
            left-most to the right-most node in the linear thermal network.
        time_eval: Sequence
            Times at which the node temperatures are calculated, ordered from
            the initial towards the final time moment, measured in seconds from
            time zero.
        continuous: bool, default False
            Indicate whether to compute a continuous solution or not.

        Returns
        -------
        A 3-member tuple:
        - A list with the times at which the node temperatures were calculated,
          measured in seconds from time zero.
        - A list of lists with for each node in the linear thermal network the
          calculated node temperatures (`Quantity` objects) at each time moment.
        - A list of lists with for each node in the linear thermal network the
          calculated flow of heat into the node (`Quantity` objects) at each
          time moment. For the last node in the linear thermal network (i.e.,
          the right-most node), also the flow of heat out of the node is added.
        """
        # Convert the initial node temperatures to Kelvin:
        T_node_init_seq = [
            T_node_init.to('K').magnitude
            for T_node_init in T_node_init_seq
        ]
        # Solve the linear thermal network for the node temperatures as a
        # function of time:
        sol = solve_ivp(
            self._network_equation, time_span, T_node_init_seq,
            method='LSODA',
            t_eval=time_eval, dense_output=continuous
        )
        # Solve the linear thermal network for the heat flows between nodes as
        # a function of time:
        Q_dot_table = self._calculate_heat_flows(sol.t, sol.y)
        # Set the list with the calculated time moments to attribute `time_axis`:
        self.time_axis = sol.t.tolist()
        # Set the solved node temperatures in a 2D-list as `Quantity` objects to
        # attribute `T_node_table`:
        self.T_node_table = [
            [Q_(T, 'K') for T in sol.y[i]]
            for i in range(sol.y.shape[0])
        ]
        # Set the solved heat flows between the nodes in the linear thermal
        # network in a 2D-list as `Quantity` objects to attribute `Q_dot_table`:
        self.Q_dot_table = [
            [Q_(Q_dot, 'W') for Q_dot in Q_dot_table[i]]
            for i in range(Q_dot_table.shape[0])
        ]
        return self.time_axis, self.T_node_table, self.Q_dot_table

    def get_T_node_table(self, unit: str = 'degC') -> pd.DataFrame:
        """Returns the solution of the node temperatures in the units specified
        by parameter `unit` (by default °C) in a Pandas DataFrame.
        The column titles of the dataframe are the IDs of the nodes in the
        linear thermal network.
        The row-index corresponds with the time-axis (the time moments at which
        the node temperatures were calculated) expressed in hours from time zero.
        """
        d = {'time (hrs)': [t / 3600 for t in self.time_axis]}
        d.update({
            node.ID: [T.to(unit).m for T in self.T_node_table[i]]
            for i, node in enumerate(self)
        })
        df = pd.DataFrame(d)
        df = df.set_index('time (hrs)', drop=True)
        return df

    def get_Q_dot_table(self, unit: str = 'W') -> pd.DataFrame:
        """Returns the solution of the heat flows between the nodes in the
        linear thermal network in the units specified by parameter `unit` (by
        default W) in a Pandas DataFrame.
        The column titles of the dataframe are the IDs of the nodes in the
        linear thermal network. The column values are the heat flows into
        the nodes (the heat flow into one node is the same as the heat flow out
        of its preceding node). The last column 'OUT' contains the values of the
        heat flow out of the last, right-most node toward the interior
        environment.
        The row-index corresponds with the time-axis (the time moments at which
        the node temperatures were calculated) expressed in hours from time zero.
        """
        d = {'time (hrs)': [t / 3600 for t in self.time_axis]}
        d.update({
            node.ID: [Q_dot.to(unit).m for Q_dot in self.Q_dot_table[i]]
            for i, node in enumerate(self)
        })
        d.update({'OUT': [Q_dot.to(unit).m for Q_dot in self.Q_dot_table[-1]]})
        df = pd.DataFrame(d)
        df = df.set_index('time (hrs)', drop=True)
        return df
