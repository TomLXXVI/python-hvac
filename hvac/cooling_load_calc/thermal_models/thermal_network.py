from __future__ import annotations

import warnings
from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import Enum
import numpy as np
import control as ct
import pandas as pd
from hvac import Quantity
from .units import Units, UNITS


class ThermalNetworkComponent(ABC):
    preferred_unit: UNITS.Unit

    def __init__(self, value: Quantity) -> None:
        self.value = value.to(self.preferred_unit)

    @classmethod
    @abstractmethod
    def update(cls) -> None:
        ...

    def __repr__(self):
        return f"{self.value:~P.3g}"


class Resistor(ThermalNetworkComponent):
    """Represents a thermal resistor between two temperature nodes in a linear 
    thermal network.
    """
    preferred_unit = Units.unit_R

    @classmethod
    def update(cls) -> None:
        cls.preferred_unit = Units.unit_R

    def __add__(self, other: Resistor) -> Resistor:
        """Adds the values of two resistors in series together and returns a 
        single new resistor.
        """
        return Resistor(self.value + other.value)
    
    def __radd__(self, other: Resistor | int) -> Resistor:
        if isinstance(other, int) and other == 0:
            return self
        elif isinstance(other, Resistor):
            return self.__add__(other)
        else:
            raise NotImplementedError(
                f"Cannot add {type(other)} and {type(self)}."
            ) from None
    
    def __floordiv__(self, other: Resistor) -> Resistor:
        new_value = 1 / (1 / self.value + 1 / other.value)
        return Resistor(new_value)
    

class Capacitor(ThermalNetworkComponent):
    """Represents the thermal capacitor of a temperature node in a linear 
    thermal network.
    """
    preferred_unit = Units.unit_C

    @classmethod
    def update(cls) -> None:
        cls.preferred_unit = Units.unit_C


class HeatFlow(ThermalNetworkComponent):
    """Represents heat flow to or from a temperature node in a linear thermal 
    network.
    """
    preferred_unit = Units.unit_Q

    @classmethod
    def update(cls) -> None:
        cls.preferred_unit = Units.unit_Q


Units.register(Resistor, Capacitor)


class ConnectionSide(Enum):
    """Enum class to specify the connection side of a temperature node, either 
    a connection on the heat flow input side, or a connection on the heat flow 
    output side.
    """
    IN = 'in'
    OUT = 'out'
    
    def reverse(self) -> ConnectionSide:
        """Reverse the connection side."""
        if self.name == self.IN:
            return self.OUT
        return self.IN


class TemperatureNode:
    """Represents a node in a linear thermal network. For each node in the 
    network a heat balance equation can be written. The set of heat balance 
    equations of all the nodes in a network can be written in a state-space 
    representation (system matrix A, input matrix B, output matrix C, and 
    feedforward matrix D).
    The heat balance equation of each node in the network corresponds with a row
    in the system or state matrix of the network, and a row in the input matrix 
    of the network. 
    The coefficients in the system matrix link node temperatures (state  
    variables) with each other, while the coefficients in the input matrix link 
    any external temperatures or external heat flows (input variables) to the 
    node temperatures.
    """
    def __init__(
        self,
        name: str,
        network_name: str,
        capacitor: Capacitor,
    ) -> None:
        """Creates a `TemperatureNode` object.
        
        Parameters
        ----------
        name:
            Unique name to identify the node in the linear thermal network.
        network_name:
            Name of the linear thermal network to which the node belongs.
        capacitor:
            Thermal capacitor connected to the node, an instance of class 
            `Capacitor`.
        """
        self.name = f"{name}@{network_name}"  # name of the node object
        self.network_name = network_name  # name of the network the node belongs to
        self.state_var_name = self.name  # name of the node's state variable (i.e. the node temperature)
        self.capacitor = capacitor
        self.resistors_in: list[tuple[str | TemperatureNode, Resistor]] | None = None
        self.resistors_out: list[tuple[str | TemperatureNode, Resistor]] | None = None
        self.heat_flows_in: list[str] | None = None
        self.heat_flows_out: list[str] | None = None
        self.num_state_vars: int = 0  # number of state variables (node temperatures) in the node heat balance equation
        self.state_var_names: list[str] = []  # names of the state variables present in the node equation
        self.num_input_vars: int = 0  # number of input variables in the node equation
        self.input_vars_names: list[str] = []  # names of the input variables present in the node equation
        self.state_dict: dict[str, float] = {}  # key = node temperature name, value = its coeff. in the node equation
        self.input_dict: dict[str, float] = {}  # key = input variable name, value = its coeff. in the node equation
    
    def connect_resistor(
        self, 
        other_terminal: str | TemperatureNode, 
        resistor: Resistor,
        this_side: ConnectionSide
    ) -> None:
        """Connects a thermal resistor to the temperature node. A thermal 
        resistor connects the temperature node with another temperature node in 
        the network, or with a terminal on the other end of the resistor which 
        is external to the network.
        
        Parameters
        ----------
        other_terminal:
            In case the other end of the resistor is connected to a temperature
            node in the network, that `TemperatureNode` object. Otherwise, a
            unique string to identify the external terminal.
        resistor:
            The connecting thermal resistor, an instance of class `Resistor`.
        this_side: {ConnectionSide.IN, ConnectionSide.OUT}
            Indicates the side of the temperature node where the resistor needs
            to be connected to. A temperature node has two sides: on one side,
            heat can flow into the node, while on the other side heat can flow 
            out of the node. The side of the node is indicated with the enum 
            class `ConnectionSide`. 
        """
        if isinstance(other_terminal, str) and '@' not in other_terminal:
            other_terminal = f"{other_terminal}@{self.network_name}"
        if this_side == ConnectionSide.OUT:
            if self.resistors_out is None:
                self.resistors_out = [(other_terminal, resistor)]
            else:
                self.resistors_out.append((other_terminal, resistor))
        if this_side == ConnectionSide.IN:
            if self.resistors_in is None:
                self.resistors_in = [(other_terminal, resistor)]
            else:
                self.resistors_in.append((other_terminal, resistor))
    
    def add_heat_flow(
        self,
        name: str,
        side: ConnectionSide
    ) -> None:
        """Adds a heat flow to the temperature node. A heat flow is always 
        external to the network.
        
        Parameters
        ----------
        name:
            Name to identify the heat flow.
        side: {ConnectionSide.IN, ConnectionSide.OUT}
            Indicates the side of the temperature node where the heat flow needs
            to be connected to. A temperature node has two sides: on one side,
            heat can flow into the node, while on the other side heat can flow 
            out of the node. The side is specified with the enum class 
            `ConnectionSide` 
        """
        if isinstance(name, str): 
            name = f"{name}@{self.network_name}"
        if side == ConnectionSide.IN:
            if self.heat_flows_in is None:
                self.heat_flows_in = [name]
            else:
                self.heat_flows_in.append(name)
        if side == ConnectionSide.OUT:
            if self.heat_flows_out is None:
                self.heat_flows_out = [name]
            else:
                self.heat_flows_out.append(name)
    
    def _get_states(self) -> tuple[int, list[str]]:
        """Returns the names of the state variables (node temperatures) in the
        heat balance equation of this node.
        """
        states = []
        # iterate over the heat input resistors connected to the node; if the
        # terminal on the other side of the resistor (tup[0]) is a temperature
        # node, its temperature is a state variable in the heat balance equation
        # of this node.
        if self.resistors_in is not None:
            for tup in self.resistors_in:
                if isinstance(tup[0], TemperatureNode):
                    states.append(f"{tup[0].state_var_name}")
        # the temperature of this node is always a state variable in its own
        # heat balance equation
        states.append(f"{self.state_var_name}")
        # iterate over the heat output resistors...
        if self.resistors_out is not None:
            for tup in self.resistors_out:
                if isinstance(tup[0], TemperatureNode):
                    states.append(f"{tup[0].state_var_name}")
        num_states = len(states)
        return num_states, states
    
    def _get_inputs(self) -> tuple[int, list[str]]:
        """Returns the names of the input variables (external temperatures or
        heat flows) in the heat balance equation of this node.
        """
        inputs = []
        # Iterate over the heat input resistors connected to the node; if the
        # terminal on the other side of the resistor (tup[0]) is not a 
        # temperature node (and therefore a string), its temperature is an input
        # variable in the heat balance equation of this node.
        if self.resistors_in is not None:
            for tup in self.resistors_in:
                if isinstance(tup[0], str):
                    inputs.append(f"{tup[0]}")
        # Iterate over the heat output resistors...
        if self.resistors_out is not None:
            for tup in self.resistors_out:
                if isinstance(tup[0], str):
                    inputs.append(f"{tup[0]}")
        # Iterate over the input heat flows...
        if self.heat_flows_in is not None:
            for label in self.heat_flows_in:
                inputs.append(f"{label}")
        # Iterate over the output heat flows...
        if self.heat_flows_out is not None:
            for label in self.heat_flows_out:
                inputs.append(f"{label}")
        num_inputs = len(inputs)
        return num_inputs, inputs
    
    def _get_coefficient(
        self, 
        resistor: Resistor | None = None,
        heat_flow: bool = False
    ) -> float:
        """Calculates a coefficient in the heat balance equation of this node
        given either the connecting resistor, or the heat flow.
        """
        C = self.capacitor.value.m
        if resistor is not None:
            R = resistor.value.m
            coeff = 1 / (R * C)
        elif heat_flow:
            coeff = 1 / C
        else:
            raise ValueError('Either `resistor` or `heat_flow` must be given.')
        return float(coeff)
    
    def _get_system_middle_coeff(self) -> list[float]:
        """Calculates the coefficient that goes with the temperature of this 
        node in the heat balance equation of this node.
        """
        def _calc_coeff(resistor_tuples) -> float:
            _, resistors = zip(*resistor_tuples)
            resistors = list(resistors)
            coeff_ = sum(self._get_coefficient(resistor=r) for r in resistors)
            return coeff_
        
        coeff = 0.0
        if self.resistors_in is not None:
            coeff += _calc_coeff(self.resistors_in)
        if self.resistors_out is not None:
            coeff += _calc_coeff(self.resistors_out)
        return [-1 * coeff]
    
    def _get_system_left_coeffs(self) -> list[float]:
        """Calculates the coefficients that go with the temperatures of the 
        nodes on the heat input side of this node in the heat balance equation
        of this node. These coefficients are placed to the left of the 
        coefficient of this node in the row of the system matrix that 
        corresponds with this node.
        """
        coeffs = []
        if self.resistors_in is not None:
            for tup in self.resistors_in:
                if isinstance(tup[0], TemperatureNode):
                    coeffs.append(self._get_coefficient(resistor=tup[1]))
        return coeffs
    
    def _get_system_right_coeffs(self) -> list[float]:
        """Calculates the coefficients that go with the temperatures of the 
        nodes on the heat output side of this node in the heat balance equation
        of this node. These coefficients are placed to the right of the 
        coefficient of this node in the row of the system matrix that 
        corresponds with this node.
        """
        coeffs = []
        if self.resistors_out is not None:
            for tup in self.resistors_out:
                if isinstance(tup[0], TemperatureNode):
                    coeffs.append(self._get_coefficient(resistor=tup[1]))
        return coeffs
    
    def _get_state_coeffs(self) -> list[float]:
        """Collects the coefficients which are placed in the row of the system
        matrix that corresponds with this node. Coefficients in the system 
        matrix belong to state variables, which are temperatures of connected
        nodes in the network.
        """
        coeffs = []
        left_coeffs = self._get_system_left_coeffs()
        middle_coeff = self._get_system_middle_coeff()
        right_coeffs = self._get_system_right_coeffs()
        if left_coeffs:
            coeffs.extend(left_coeffs)
        if middle_coeff:
            coeffs.extend(middle_coeff)
        if right_coeffs:
            coeffs.extend(right_coeffs)
        if not coeffs:
            coeffs = [0.0]
        return coeffs
    
    def _get_input_coeffs(self) -> list[float]:
        """Collects the coefficients which are placed in the row of the 
        state-space input matrix that corresponds with this node. Coefficients 
        in the input matrix belong to input variables, which are external 
        temperatures or heat flows to which no other nodes in the network are 
        attached.
        """
        coeffs = []
        if self.resistors_in is not None:
            for tup in self.resistors_in:
                if isinstance(tup[0], str):
                    coeffs.append(self._get_coefficient(resistor=tup[1]))
        if self.resistors_out is not None:
            for tup in self.resistors_out:
                if isinstance(tup[0], str):
                    coeffs.append(self._get_coefficient(resistor=tup[1]))
        if self.heat_flows_in is not None:
            for _ in self.heat_flows_in:
                coeffs.append(self._get_coefficient(heat_flow=True))
        if self.heat_flows_out is not None:
            for _ in self.heat_flows_out:
                coeffs.append(-1 * self._get_coefficient(heat_flow=True))
        if not coeffs:
            coeffs = [0.0]
        return coeffs
    
    @staticmethod
    def _clean_up(
        resistors: list[tuple[str | TemperatureNode, Resistor]]
    ) -> list[tuple[str | TemperatureNode, Resistor]] | None:
        """Checks if multiple thermal resistors are connected to the same 
        `other_terminal`. If true, combine these parallel resistors into one
        equivalent resistor.
        """
        if resistors is not None:
            resistors_dict = {}
            for tup in resistors:
                if tup[0] not in resistors_dict.keys():
                    resistors_dict[tup[0]] = tup[1]
                else:
                    R = resistors_dict[tup[0]]
                    R = R // tup[1]
                    resistors_dict[tup[0]] = R
            resistors = [(k, v) for k, v in resistors_dict.items()]
            return resistors
        return None
        
    def equation(self) -> tuple[dict[str, float], ...]:
        """Returns the state and input part of the heat balance equation of the 
        node. The state and input part of the heat balance equation are each 
        represented by a dict. The state-dict maps the names of the state 
        variables (node temperatures) to their associated coefficients in 
        the node equation. The input-dict maps the names of input variables to 
        their associated coefficients in the node equation.
        """
        self.resistors_in = self._clean_up(self.resistors_in)
        self.resistors_out = self._clean_up(self.resistors_out)
        self.num_state_vars, self.state_var_names = self._get_states()
        A = self._get_state_coeffs()
        self.state_dict = {lbl: a for lbl, a in zip(self.state_var_names, A)}
        self.num_input_vars, self.input_vars_names = self._get_inputs() 
        B = self._get_input_coeffs()
        self.input_dict = {lbl: b for lbl, b in zip(self.input_vars_names, B)}
        return self.state_dict, self.input_dict
    
    def __repr__(self) -> str:
        return self.name


@dataclass
class Connection:
    """Data class to define a connection between two nodes that belong to two 
    different linear thermal networks. The connection is made by placing a 
    thermal resistor between the nodes.
    
    Attributes
    ----------
    node: TemperatureNode
        Node in the source network where the connection leaves.
    other_network: LinearThermalNetwork
        The other network where the connection arrives at a node.
    other_node: TemperatureNode
        Node in the other network where the connection arrives.
    resistor: Resistor
        Thermal resistor of the connection.
    other_side: {ConnectionSide.IN, ConnectionSide.OUT}
        The side of the other node where the connection ends. If 
        `ConnectionSide.IN`, it is assumed that heat flow is directed to the 
        other node. If `ConnectionSide.OUT`, it is assumed that heat flow is 
        directed away from the other node. In both cases the connection endpoint
        on the source node must be on the opposite side. 
    """
    node: TemperatureNode
    other_network: LinearThermalNetwork
    other_node: TemperatureNode
    resistor: Resistor
    other_side: ConnectionSide


class LinearThermalNetwork:
    """Represents a linear thermal network of temperature nodes interconnected 
    by thermal resistors. The set of heat balance node equations in the linear 
    thermal network can be described mathematically by a system using 
    state-space representation.
    """
    def __init__(
        self, 
        name: str, 
        nodes: tuple[TemperatureNode, ...],
    ) -> None:
        """Creates a `LinearThermalNetwork` object.
        
        Parameters
        ----------
        name:
            Name to identify the linear thermal network.
        nodes:
            Tuple with all the nodes that constitute the linear thermal network. 
            The nodes must be ordered according to their spatial position in the
            linear network.
        """
        self.name = name
        self.nodes = nodes
        self._state_var_mapping: dict[str, int] | None = None
        self._input_var_mapping: dict[str, int] | None = None
        self._output_nodes: tuple[TemperatureNode, ...] | None = None
        self._A: np.ndarray[float] | None = None
        self._B: np.ndarray[float] | None = None
        self._C: np.ndarray[float] | None = None
        self._D: np.ndarray[float] | None = None
    
    def _get_state_var_mapping(self) -> dict[str, int]:
        """Loops over the nodes in the linear thermal network. Reads the names 
        of the state variables (node temperatures) in the node equation. Each 
        state variable in the global system of the network is assigned a unique 
        index which will determine its row and/or column position in the
        system matrix.
        Returns a dictionary which maps the names of the state variables to 
        their index.
        """
        state_indices = {}
        i = 0
        for node in self.nodes:
            states_dict, _ = node.equation()
            for key in states_dict.keys():
                if key not in state_indices.keys():
                    state_indices[key] = i
                    i += 1
        return state_indices
    
    def _get_input_var_mapping(self) -> dict[str, int]:
        """Loops over the nodes in the linear thermal network. Reads the
        names of the input variables (external temperatures or heat flows) in 
        the node equation. Each input variable in the global system of the 
        network is assigned a unique index which will determine its column 
        position in the input matrix.
        Returns a dictionary which maps the names of the input variables to
        their index.
        """
        input_indices = {}
        i = 0
        for node in self.nodes:
            _, inputs_dict = node.equation()
            for key in inputs_dict.keys():
                if key not in input_indices.keys():
                    input_indices[key] = i
                    i += 1
        return input_indices
    
    def _create_system_matrix(self) -> np.ndarray:
        """Creates the square system matrix of the linear thermal network. The
        number of rows and columns equal the number of nodes in the network.
        Starts with a zero system matrix. Loops over each node in the network. 
        Reads the state variables and their associated coefficients in the node 
        equation. The row position of each coefficient in the system matrix is 
        determined by the assigned index of the node in the linear network; its 
        column position is determined by the assigned index of its state 
        variable.
        """
        n = len(self.nodes)
        A = np.zeros((n, n))
        for node in self.nodes:
            i = self._state_var_mapping[node.state_var_name]
            states_dict, _ = node.equation()
            for key, coeff in states_dict.items():
                j = self._state_var_mapping[key]
                A[i, j] = coeff
        return A
    
    def _create_input_matrix(self) -> np.ndarray:
        """Creates the input matrix of the linear thermal network. The number
        of rows equals the number of nodes in the network. The number of columns
        equals the number of input variables (external temperatures or heat 
        flows) in the network.
        Starts with a zero input matrix. Loops over each node in the network.
        Reads the input variables and their associated coefficients in the node 
        equation. The row position of each coefficient in the input matrix is 
        determined by the assigned index of the node in the linear network; its 
        column position is determined by the index of its input variable.
        """
        n = len(self.nodes)
        m = len(self._input_var_mapping)
        B = np.zeros((n, m))
        for node in self.nodes:
            i = self._state_var_mapping[node.state_var_name]
            _, inputs_dict = node.equation()
            for key, coeff in inputs_dict.items():
                j = self._input_var_mapping[key]
                B[i, j] = coeff
        return B
    
    def _create_output_matrix(self, output_nodes: tuple[TemperatureNode, ...]) -> np.ndarray:
        """Creates the output matrix of the linear thermal network. The number
        of rows equals the number of outputs. The number of columns equals the 
        number of state variables (node temperatures) in the network. 
        """
        p = len(output_nodes)
        n = len(self.nodes)
        C = np.zeros((p, n))
        for i, output_node in enumerate(output_nodes):
            try:
                output_index = self._state_var_mapping[output_node.state_var_name]
            except KeyError:
                raise KeyError(
                    f"Node '{output_node.name}' is not in this network."
                ) from None
            C[i, output_index] = 1
        return C
    
    def _create_feedforward_matrix(self, num_output_nodes: int = 1) -> np.ndarray:
        """Creates the feedforward matrix of the linear thermal network. The
        number of rows equals the number of output variables. The number of 
        columns equals the number of input variables (external temperatures or 
        heat flows) in the network. The feedforward matrix is always a zero 
        matrix (input variables are not coupled to output variables in the 
        system).
        """
        p = num_output_nodes
        m = len(self._input_var_mapping)
        D = np.zeros((p, m))
        return D
    
    def create_system(
        self, 
        output_nodes: tuple[TemperatureNode, ...] | None = None,
        reduced_order: int | None = None
    ) -> ct.StateSpace:
        """Creates the system of the linear thermal network in state-space 
        representation.
        
        Parameters
        ----------
        output_nodes:
            Nodes in the linear thermal network whose temperature will be the 
            output of the system. If `None`, the last node in the `nodes` list
            is taken to be the output node.
        reduced_order:
            The order to which the system should be reduced. If `None`, the 
            order of the system won't be reduced.
            As the system is build from a big number of temperature nodes, the
            system order will also be high. However, without noticeable loss of 
            accuracy, the order can normally be reduced to an order of 8 or 6
            (but not below the number of system outputs; should `reduced_order`
            be less than the number of system outputs, it will be limited to the
            number of system outputs). 
            Note that without system reduction, it may also happen that some 
            internal mathematical operations raise an exception (e.g. 
            `FloatingPointError`).
        
        Returns
        -------
        `ct.StateSpace` object.
        """
        if output_nodes is None: 
            self._output_nodes = (self.nodes[-1],)
        else:
            self._output_nodes = output_nodes
        self._state_var_mapping = self._get_state_var_mapping()
        self._input_var_mapping = self._get_input_var_mapping()
        self._A = self._create_system_matrix()
        self._B = self._create_input_matrix()
        self._C = self._create_output_matrix(self._output_nodes)
        self._D = self._create_feedforward_matrix(len(self._output_nodes))
        system = ct.ss(self._A, self._B, self._C, self._D)
        system.name = self.name
        system.set_inputs(list(self._input_var_mapping.keys()))
        system.set_states(list(self._state_var_mapping.keys()))
        system.set_outputs([output_node.state_var_name for output_node in self._output_nodes])
        if reduced_order is not None and system.nstates > reduced_order:
            red_system = ct.balanced_reduction(system, reduced_order, method='matchdc')
            red_system.name = self.name
            red_system.set_inputs(list(self._input_var_mapping.keys()))
            red_system.set_outputs([output_node.state_var_name for output_node in self._output_nodes])
            return red_system
        return system
    
    @property
    def A(self) -> pd.DataFrame | None:
        """Returns the system matrix of the linear thermal network as a
        Pandas DataFrame object.
        """
        if self._A is not None:
            df = pd.DataFrame(
                data=self._A,
                index=list(self._state_var_mapping.keys()),
                columns=list(self._state_var_mapping.keys())
            )
            return df
        else:
            warnings.warn(
                f"`A` is None. Call `create_system(...)` first.",
                category=UserWarning
            )
    
    @property
    def B(self) -> pd.DataFrame | None:
        """Returns the input matrix of the linear thermal network as a Pandas
        DataFrame object.
        """
        if self._B is not None:
            df = pd.DataFrame(
                data=self._B,
                index=list(self._state_var_mapping.keys()),
                columns=list(self._input_var_mapping.keys())
            )
            return df
        else:
            warnings.warn(
                f"`B` is None. Call `create_system(...)` first.",
                category=UserWarning
            )
    
    @property
    def C(self) -> pd.DataFrame | None:
        """Returns the output matrix of the linear thermal network as a Pandas
        DataFrame object.
        """
        if self._C is not None:
            df = pd.DataFrame(
                data=self._C,
                columns=list(self._state_var_mapping.keys())
            )
            return df
        else:
            warnings.warn(
                f"`C` is None. Call `create_system(...)` first.",
                category=UserWarning
            )
    
    @property
    def D(self) -> pd.DataFrame | None:
        """Returns the feedforward matrix of the linear thermal network as a 
        Pandas DataFrame object.
        """
        if self._D is not None:
            df = pd.DataFrame(
                data=self._D,
                columns=list(self._input_var_mapping.keys())
            )
            return df
        else:
            warnings.warn(
                f"`D` is None. Call `create_system(...)` first.",
                category=UserWarning
            )
    
    @property
    def inputs(self) -> tuple[str, ...] | None:
        """Returns the names of the input variables of the system."""
        if self._input_var_mapping is not None:
            return tuple(self._input_var_mapping.keys())
        else:
            warnings.warn(
                f"`inputs` is None. Call `create_system(...)` first."
            )
    
    @property
    def outputs(self) -> tuple[str, ...] | None:
        """Returns the names of the output variables of the system."""
        if self._output_nodes is not None:
            return tuple(output_node.name for output_node in self._output_nodes)
        else:
            warnings.warn(
                f"`output` is None. Call `create_system(...)` first."
            )
    
    def get_input_index(self, input_name: str) -> int:
        try:
            input_index = self._input_var_mapping[input_name]
        except KeyError:
            raise KeyError(
                f"Input '{input_name}' is not an input of the system."
            ) from None
        return input_index
    
    def connect(
        self, 
        connections: list[Connection],
        name: str
    ) -> LinearThermalNetwork:
        """Connects this linear thermal network with one or more other linear 
        thermal networks and returns the new resulting linear thermal network.
        
        Parameters
        ----------
        connections:
            Tuple of `Connection` objects. See docstring of data class 
            `Connection` on how to define a connection.
        name:
            Name for the resulting linear thermal network.
            
        Returns
        -------
        LinearThermalNetwork
        """
        def _apply_connection(conn: Connection):
            # Add the connecting resistor on the specified connection side of
            # the node in the other network to the node in this network.
            conn.other_node.connect_resistor(conn.node, conn.resistor, conn.other_side)
            # Also add the connecting resistor on the opposite connection side
            # of the node in this network to the node in the other network.
            conn.node.connect_resistor(conn.other_node, conn.resistor, conn.other_side.reverse())
        
        # Apply all connections.
        for connection in connections: _apply_connection(connection)
        
        # Collect the other networks to which a connection was made. A dict is
        # used to make sure that another network is not collected more than once.
        other_networks_dict = {}
        for connection in connections:
            other_network = connection.other_network
            other_networks_dict.setdefault(other_network.name, other_network)
        
        # Take the nodes of all networks and return the new network.
        nodes = []
        for nw in other_networks_dict.values(): nodes.extend(list(nw.nodes))
        nodes.extend(self.nodes)
        return LinearThermalNetwork(name, tuple(nodes))
