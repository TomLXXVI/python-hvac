import warnings
from typing import Dict, List, Optional, Tuple, Union
from abc import ABC
# import threading
from copy import copy
import csv
import pandas as pd
import dill as pickle
from hvac import Quantity
from hvac.fluids import Fluid, FluidState
from ..conduit.conduit import TConduit, Node, NodeArrow, Loop, Conduit, PseudoConduit, FlowSign
from ..conduit.cross_section import Rectangular, FlatOval, Circular, TCrossSection
from ..fittings.pipe import BalancingValve, ControlValve
from ..fittings.duct import ObstructionA15A
from ..schedule import PipeSchedule, DuctSchedule, PipeScheduleFactory, DuctScheduleFactory
from ..utils import SystemCurve


Q_ = Quantity
TSchedule = Union[PipeSchedule, DuctSchedule]
TVolumeDamper = Union[ObstructionA15A]
Air = Fluid('Air')
STANDARD_AIR = Air(T=Q_(20, 'degC'), P=Q_(101_325, 'Pa'))


class AnalysisWarning(Warning):
    """Warning for when the maximum number of iterations has been reached during
    the analysis of a network with the Hardy Cross method.
    """
    pass


def pretty_unit(unit: str) -> str:
    """Returns the prettified form of the given unit."""
    q = Quantity(0, unit)
    q = f"{q:~P}".split(' ')
    return q[1]


class FlowPath(List[TConduit]):
    """Class derived from `list` representing a flow path between the start
    and end node of a network.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.ambient_air: FluidState = STANDARD_AIR
        self._conduits: List[Conduit] = []

    def __str__(self):
        return '|'.join(conduit.ID for conduit in self.conduits)

    @property
    def conduits(self) -> List[Conduit]:
        """Get the conduits that belong to the flow path."""
        if not self._conduits:
            self._conduits = [
                conduit for conduit in self
                if not isinstance(conduit, PseudoConduit)
            ]
        return self._conduits

    @property
    def dyn_pressure_diff(self) -> Quantity:
        """Get the total pressure loss between the start and end node of the
        flow path due to fluid flow friction and fitting losses.
         """
        return sum(conduit.pressure_drop for conduit in self.conduits)

    @property
    def elev_pressure_diff(self) -> Quantity:
        """Returns the elevation pressure difference (aka thermal gravity effect
        or chimney effect) due to any height difference between the start and
        end node of the flow path.
        """
        try:
            first_node = self.conduits[0].start_node
        except IndexError as err:
            print('Oops')
            raise err
        last_node = self.conduits[-1].end_node
        z1 = first_node.height.to('m')
        z2 = last_node.height.to('m')
        rho = self.conduits[0].fluid.rho.to('kg / m ** 3')
        g = Q_(9.81, 'm / s ** 2')
        rho_a = self.ambient_air.rho.to('kg / m ** 3')
        return g * (rho - rho_a) * (z2 - z1)

    @property
    def tot_pressure_diff(self) -> Quantity:
        """Returns the total pressure difference between the start and end node
        of the flow path. It is the sum of the elevation pressure difference
        and the dynamic pressure difference.
        """
        dp_elev = self.elev_pressure_diff
        dp_dyn = self.dyn_pressure_diff
        dp_tot = dp_elev + dp_dyn  # dp_tot = (tp1 - tp2) + dp_pump
        return dp_tot

    def get_system_curve(
        self,
        V_wp: Quantity,
        dp_wp: Quantity,
        p_g1: Optional[Quantity] = Q_(0.0, 'Pa')
    ) -> SystemCurve:
        """Returns a `SystemCurve` object, representing the system curve of the
        flow path.

        Parameters
        ----------
        V_wp:
            Volume flow rate through the flow path.
        dp_wp:

        p_g1:
            Static pressure at the start node of the flow path.

        Returns
        -------
        SystemCurve
        """
        dp_dyn = self.dyn_pressure_diff
        dp_elev = self.elev_pressure_diff
        dp_loss = dp_dyn + dp_wp
        R_hyd = dp_loss / (V_wp ** 2)  # hydraulic resistance of the flow path
        p_v1 = self.conduits[0].velocity_pressure
        p_t1 = p_g1 + p_v1  # total pressure at start node of the flow path
        p_t2 = p_t1 - dp_elev - dp_dyn  # total pressure at end node
        # p_v2 = self.conduits[-1].velocity_pressure
        # p_g2 = p_t2 - p_v2  # static pressure at end node
        # dp_g = p_g2 - p_g1  # static pressure difference between end and start node
        # dp_v = p_v2 - p_v1  # velocity pressure difference between end and start node
        # dp_t = dp_g + dp_v  # difference between total pressure at end and start node
        dp_t = p_t2 - p_t1
        sys_curve = SystemCurve.create(R_hyd=R_hyd, dP_tot=dp_t, dP_elev=dp_elev)
        return sys_curve


class Network(ABC):
    default_units = {
        'length': 'm',
        'diameter': 'mm',
        'dim_width': 'mm',
        'dim_height': 'mm',
        'volume_flow_rate': 'm ** 3 / s',
        'velocity': 'm/s',
        'pressure': 'Pa',
        'specific_pressure': 'Pa / m',
        'pressure_head': 'm',
        'height': 'm',
        'wall_roughness': 'mm',
        'machine_coefficients': ['Pa', 'Pa / (m ** 3 / s)', 'Pa / (m ** 3 / s) ** 2']
    }

    def __init__(self):
        self.ID: str = ''
        self.start_node: Node = Node()
        self.end_node: Node = Node()
        self.fluid: FluidState = STANDARD_AIR
        self.ambient_air: FluidState = STANDARD_AIR
        self.wall_roughness: Quantity = Q_(0.09, 'mm')
        self.schedule: Optional[TSchedule] = None
        self.units = self.default_units
        self.nodes: Dict[str, Node] = {}
        self.conduits: Dict[str, TConduit] = {}
        self._flow_paths: List['FlowPath'] = []
        self.loops: Dict[str, Loop] = {}

    @classmethod
    def create(
        cls,
        ID: str,
        fluid: FluidState = STANDARD_AIR,
        ambient_air: FluidState = STANDARD_AIR,
        wall_roughness: Quantity = Q_(0.09, 'mm'),
        schedule: Optional[TSchedule] = None,
        start_node_ID: str = '',
        start_node_height: Quantity = Q_(0.0, 'm'),
        end_node_ID: str = '',
        end_node_height: Quantity = Q_(0.0, 'm'),
        units: Optional[Dict[str, str]] = None
    ) -> 'Network':
        """Creates an instance of a derived class of abstract class `Network`
        (see `PipeNetwork` for pipe systems and `DuctNetwork` for duct systems).

        Parameters
        ----------
        ID:
            Name of the network.
        fluid: default STANDARD_AIR
            The (incompressible) fluid that flows in the network. The default
            is standard air, which is dry air with a temperature of 20 °C at
            standard atmospheric pressure (101,325 Pa).
        ambient_air: default STANDARD_AIR
            The fluid of the environment where the network is located. Used for
            open networks to calculate the elevation pressure difference
            (gravity effect or chimney effect) along a flow path.
        wall_roughness: default 0.09 mm
            The absolute wall roughness of the pipes or ducts. The default value
            applies to medium smooth galvanized air ducts.
        schedule:
            The schedule contains the commercially available sizes of a type
            of duct or pipe. See module schedule.duct_schedule.py and module
            schedule.pipe_schedule.py
        start_node_ID:
            ID for the start node of the network.
        start_node_height: default 0 m
            Height of the start node with respect to a reference plane to
            which the heights of all network nodes are referred to.
        end_node_ID:
            ID of the network's end node.
        end_node_height: default 0 m
            Height of the end node referred to the same reference plane as the
            network's start node.
        units:
            Dictionary of which the keys are the names of the quantities and the
            values are the units to be used for these quantities. See class
            attribute `default_units` for a list of the quantities.
            Also, when reading a network configuration from csv-file, the values
            in the table will be assumed to be expressed in the units set through
            this parameter.
        """
        obj = cls()
        obj.ID = ID
        obj.fluid = fluid
        obj.ambient_air = ambient_air
        obj.wall_roughness = wall_roughness
        obj.schedule = schedule

        if units is not None:
            obj.units = dict(cls.default_units, **units)
        else:
            obj.units = cls.default_units

        obj.start_node = Node(start_node_ID, start_node_height)
        obj.nodes[obj.start_node.ID] = obj.start_node
        if end_node_ID is not None:
            obj.end_node = Node(end_node_ID, end_node_height)
            obj.nodes[obj.end_node.ID] = obj.end_node
        return obj

    def add_conduit(
        self,
        conduit: TConduit,
        conduit_ID: str,
        start_node_ID: str,
        end_node_ID: str,
        start_node_height: Quantity = Q_(0.0, 'm'),
        end_node_height: Quantity = Q_(0.0, 'm'),
        loop_ID: Optional[Union[str, Tuple[str, str]]] = None,
        zeta: Optional[float] = None
    ) -> None:
        """Adds a new conduit to the network.

        Parameters
        ----------
        conduit:
            Object of type `TConduit`, which is a union of classes `Conduit`,
            `Pipe`, `Duct`, and `PseudoConduit`. See module conduit.py.
        conduit_ID:
            Identifier for the conduit in the network.
        start_node_ID:
            ID for the start node of the conduit (node where the conduit fluid
            leaves the node).
        end_node_ID:
            ID for the end node of the conduit (node where the conduit fluid
            arrives at the node)
        start_node_height: default 0 m
            Height of the conduit's start node above the common reference plane
            of the network.
        end_node_height: default 0 m
            Height of the conduit's end node above the common reference plane
            of the network.
        loop_ID: str, tuple[str, str]
            Identifier to identify the loop or loops to which a conduit belongs
            in a network. A conduit can belong to 2 loops maximum.
        zeta: optional
            The sum of the resistance coefficients of all fittings and valves
            present in the conduit.
        """
        conduit.ID = conduit_ID
        if isinstance(zeta, float) and isinstance(conduit, Conduit):
            conduit.add_fitting(zeta)
        sn = self.nodes.setdefault(
            start_node_ID,
            Node(start_node_ID, start_node_height)
        )
        en = self.nodes.setdefault(
            end_node_ID,
            Node(end_node_ID, end_node_height)
        )
        sn.connect(conduit, NodeArrow.OUTGOING)
        en.connect(conduit, NodeArrow.INCOMING)
        conduit.start_node = sn
        conduit.end_node = en
        conduit.loop_ID = loop_ID
        self.conduits[conduit_ID] = conduit

    def _create_loops(self) -> None:
        # (re)create the loops for analysis with Hardy Cross method
        self.loops = {}
        for conduit in self.conduits.values():
            loop_ID = conduit.loop_ID
            if isinstance(loop_ID, tuple) and len(loop_ID) == 2:
                # conduit belongs to 2 loops
                loop = self.loops.setdefault(
                    loop_ID[0],
                    Loop(loop_ID[0])
                )
                other_loop = self.loops.setdefault(
                    loop_ID[1],
                    Loop(loop_ID[1])
                )
                loop.append(conduit)
                conduit.loops = [
                    self.loops[loop_ID[0]],
                    self.loops[loop_ID[1]]
                ]
                if isinstance(conduit, Conduit):
                    conduit_duplicate = Conduit.duplicate(
                        conduit,
                        flow_sign=conduit.flow_sign.reverse()
                        # flow sign in the other loop is opposite to the first loop
                    )
                else:
                    conduit_duplicate = PseudoConduit.duplicate(
                        conduit,
                        flow_sign=conduit.flow_sign.reverse()
                    )
                other_loop.append(conduit_duplicate)
                conduit_duplicate.loops = [
                    self.loops[loop_ID[1]],
                    self.loops[loop_ID[0]]
                ]
            if isinstance(loop_ID, str) and len(loop_ID) > 0:
                # conduit belongs to only 1 loop
                loop = self.loops.setdefault(
                    loop_ID,
                    Loop(loop_ID)
                )
                loop.append(conduit)
                conduit.loops = [self.loops[loop_ID]]

    def _hardy_cross_algorithm(self) -> None:
        # calculate the loop correction terms of all loops in the network
        for loop in self.loops.values():
            loop.calculate_correction_term()

        # apply the correction terms to the flow rate of each section in each
        # loop
        for loop in self.loops.values():
            for conduit in loop:
                # correction terms cannot be applied to pseudo-sections as they
                # have no flow
                if not isinstance(conduit, PseudoConduit):
                    # check if the section is common to 2 loops
                    if len(conduit.loops) == 2:
                        # if a section belongs to 2 loops: get the correction
                        # term of the other loop
                        if conduit.loops[0].ID != loop.ID:
                            other_loop = conduit.loops[0]
                        else:
                            other_loop = conduit.loops[1]
                        correction_term = loop.correction_term - other_loop.correction_term
                    else:
                        correction_term = loop.correction_term
                    # calculate new flow rate of the section
                    if conduit.flow_sign == FlowSign.COUNTERCLOCKWISE:
                        conduit.volume_flow_rate += correction_term
                    else:
                        conduit.volume_flow_rate -= correction_term

    def _check_loops(self, tolerance: Quantity) -> bool:
        # check if the pressure drop along each loop is smaller than the allowable tolerance in Pa.
        if False in [abs(loop.pressure_drop) < tolerance for loop in self.loops.values()]:
            return False
        return True

    def analyze(
        self,
        tolerance: Quantity = Quantity(1.0, 'Pa'),
        i_max: int = 100
    ) -> int:
        """Analyzes a network using the Hardy-Cross method after it has been
        configured with initial guesses for the volume flow rates in the
        conduits of the network. The Hardy-Cross method will determine the
        actual volume flow rates in the network.

        Parameters
        ----------
        tolerance: default 1 Pa
            The actual volume flow rates will be attained when the sum of
            conduit pressure losses along each loop has become zero. As it
            would take many iterations for the Hardy-Cross method to reach
            zero pressure drop in every loop of the network, while volume
            flow rates hardly change anymore, parameter `tolerance` allows
            to stop iteration when the pressure drop of all loops is less
            than the pressure drop assigned to `tolerance`.
        i_max: default 100
            In case the tolerance condition is not met after `i_max` iterations
            the routine will be terminated and an `OverflowError` exception
            will be raised.

        Notes
        -----
        When assigning initial guesses to the volume flow rates, care should be
        taken that the sum of flow rates that enter a node is equal to the
        sum of flow rates that leave the node.
        """
        self._create_loops()
        i = 0
        while not self._check_loops(tolerance):
            self._hardy_cross_algorithm()
            i += 1
            if i > i_max:
                warnings.warn(
                    'Tolerance could not be reached for '
                    f'all loops after {i_max} iterations',
                    category=AnalysisWarning
                )
                break
        return i

    # def _find_flow_paths(self):
    #     path = FlowPath()
    #     path.ambient_air = self.ambient_air
    #     self._flow_paths.append(path)
    #     self._search_flow_paths(self.start_node, path)
    #
    # def _search_flow_paths(self, node: Node, path: FlowPath):
    #     loop_condition = lambda n: len(n.outgoing) >= 1 if self.end_node is None else n.ID != self.end_node.ID
    #     while loop_condition(node):
    #         if len(node.outgoing) > 1:
    #             for conduit in node.outgoing[1:]:
    #                 new_path = FlowPath(path)
    #                 new_path.ambient_air = self.ambient_air
    #                 self._flow_paths.append(new_path)
    #                 new_path.append(conduit)
    #                 thread = threading.Thread(target=self._search_flow_paths, args=(conduit.end_node, new_path))
    #                 thread.start()
    #                 thread.join()
    #         path.append(node.outgoing[0])
    #         node = path[-1].end_node

    def _search_flow_paths(self) -> None:
        node = self.start_node
        path = FlowPath()
        path.ambient_air = self.ambient_air
        self._flow_paths.append(path)
        self._recursive_path_search(node, path)

    def _recursive_path_search(self, node: Node, path: FlowPath) -> None:
        while node.ID != self.end_node.ID:
            if len(node.outgoing) > 1:
                for conduit in node.outgoing[1:]:
                    if conduit not in path and not isinstance(conduit, PseudoConduit):
                        new_flow_path = FlowPath(path)
                        new_flow_path.ambient_air = self.ambient_air
                        self._flow_paths.append(new_flow_path)
                        new_flow_path.append(conduit)
                        next_node = self.nodes[conduit.end_node.ID]
                        self._recursive_path_search(next_node, new_flow_path)
                    else:
                        return
            try:
                conduit = node.outgoing[0]
            except IndexError:
                return
            else:
                if conduit not in path and not isinstance(conduit, PseudoConduit):
                    path.append(conduit)
                    node = self.nodes[conduit.end_node.ID]
                else:
                    return

    @property
    def flow_paths(self) -> List[FlowPath]:
        """Returns a list of `FlowPath` objects, representing the flow paths
        between the start and end node of the network.
        """
        if not self._flow_paths:
            self._search_flow_paths()
        return self._flow_paths

    @property
    def critical_path(self) -> FlowPath:
        """Returns the `FlowPath` object that represents the flow path between
        the start node and end node of the network that has the greatest
        pressure loss.
        """
        dp_tots = [fp.tot_pressure_diff for fp in self.flow_paths]
        index_max = dp_tots.index(max(dp_tots))
        return self._flow_paths[index_max]

    @property
    def volume_flow_rate(self) -> Quantity:
        """Get the volume flow rate that enters the network through its start
        node (and leaves at the same rate through the end node). This is
        applicable to a network with one entry (start node) and one exit
        (end node).
        """
        return sum(
            conduit.volume_flow_rate
            for conduit in self.start_node.outgoing
        )

    @property
    def total_pressure_difference(self) -> Quantity:
        """Get the total pressure loss along the critical path of the network."""
        return self.critical_path.tot_pressure_diff

    @property
    def hydraulic_resistance(self) -> Quantity:
        """Get the hydraulic resistance of the network, calculated as
        R = dP_crit / V² with dP_crit the pressure loss along the critical
        path and V the volume flow rate that enters and leaves the network.
        """
        dp_loss = self.critical_path.dyn_pressure_diff.to('Pa')
        V = self.volume_flow_rate.to('m ** 3 / s')
        return dp_loss / (V ** 2)

    def get_flow_path_table(self) -> pd.DataFrame:
        """Returns a Pandas DataFrame object with an overview of the flow
        paths from start node to end node in the network.
        """
        headers = [
            "path",
            f"Δp-elev. [{self.units['pressure']}]",
            f"Δp-dyn. [{self.units['pressure']}]",
            f"Δp-tot. [{self.units['pressure']}]",
            f"Δp-deficit [{self.units['pressure']}]"
        ]
        table = {header: [] for header in headers}
        for fp in self.flow_paths:
            table[headers[0]].append(str(fp))
            table[headers[1]].append(fp.elev_pressure_diff.to(self.units['pressure']).magnitude)
            table[headers[2]].append(fp.dyn_pressure_diff.to(self.units['pressure']).magnitude)
            table[headers[3]].append(fp.tot_pressure_diff.to(self.units['pressure']).magnitude)
            dp_deficit = self.total_pressure_difference - fp.tot_pressure_diff
            table[headers[4]].append(dp_deficit.to(self.units['pressure']).magnitude)
        return pd.DataFrame(table)

    def get_fitting_table(self) -> pd.DataFrame:
        """Returns a Pandas DataFrame object with an overview of all the
        fittings in the network."""
        fitting_tables = [
            conduit.get_fitting_table(self.units['pressure'])
            for conduit in self.conduits.values()
            if not isinstance(conduit, PseudoConduit)
        ]
        fitting_table = pd.concat(fitting_tables, ignore_index=True)
        return fitting_table

    def _clear(self) -> None:
        # Resets the attributes of the `Network` object.
        self.nodes = {}
        self.start_node = Node(self.start_node.ID, self.start_node.height)
        self.nodes[self.start_node.ID] = self.start_node
        if self.end_node is not None:
            self.end_node = Node(self.end_node.ID, self.end_node.height)
            self.nodes[self.end_node.ID] = self.end_node
        self.conduits = {}
        self._flow_paths = []
        self.loops = {}

    def _build(self, data: Union[csv.DictReader, pd.DataFrame]) -> None:
        # Build the network from a `csv.DictReader` or Pandas `DataFrame`.
        if isinstance(data, pd.DataFrame):
            data = data.to_dict(orient='records')
        for row in data:
            # every row in `csv.DictReader` is a dict of which the keys are
            # the column titles in the table
            if not row.get('fixed_pressure_difference'):
                # create cross-section
                cross_section = self._create_cross_section(
                    diameter=row.get('diameter'),
                    nominal_diameter=row.get('nominal_diameter'),
                    width=row.get('width'),
                    height=row.get('height'),
                    shape=row.get('shape'),
                    schedule=self._get_schedule(row.get('schedule')) or self.schedule,
                    units=self.units
                )
                # create conduit with cross-section
                conduit = Conduit.create(
                    length=self._quantify(
                        row['length'],
                        self.units['length']
                    ),
                    wall_roughness=self._extract_wall_roughness(
                        row.get('wall_roughness'),
                        self.units['wall_roughness']
                    ),
                    fluid=self.fluid,
                    cross_section=cross_section,
                    volume_flow_rate=self._quantify(
                        row.get('volume_flow_rate'),
                        self.units['volume_flow_rate']
                    ),
                    pressure_drop=self._quantify(
                        row.get('pressure_drop'),
                        self.units['pressure']
                    ),
                    specific_pressure_drop=self._quantify(
                        row.get('specific_pressure_drop'),
                        self.units['specific_pressure']
                    ),
                    machine_coefficients=self._extract_machine_coefficients(
                        row.get('machine_coefficients')
                    )
                )
            else:
                # create a `PseudoConduit` object with a fixed pressure drop.
                conduit = PseudoConduit.create(
                    fixed_pressure_drop=self._quantify(
                        row['fixed_pressure_difference'],
                        self.units['pressure']
                    )
                )
            # add `Conduit` (`Pipe` or `Duct`) object or `PseudoConduit`
            # object to network
            self.add_conduit(
                conduit=conduit,
                conduit_ID=row['conduit_ID'],
                start_node_ID=row['start_node_ID'],
                start_node_height=self._quantify(
                    row.get('start_node_height', '0.0'),
                    self.units['height']
                ),
                end_node_ID=row['end_node_ID'],
                end_node_height=self._quantify(
                    row.get('end_node_height', '0.0'),
                    self.units['height']
                ),
                loop_ID=self._extract_loop_id(row.get('loop_ID')),
                zeta=self._floatify(row.get('zeta'))
            )

    def load_from_csv(self, file_path: str):
        """
        Loads a network configuration from a csv-file and creates the network.

        Available column titles:
        ** 'conduit_ID': text, mandatory
                The ID of the conduit in the network.
        ** 'start_node_ID': text, mandatory
                The ID of the start node of the conduit. The start node is
                the endpoint of the conduit where fluid enters the conduit.
        - 'start_node_height': number, optional
                The height of the start node above a fixed reference plane
                (the height of all nodes in the network must be referred to this
                reference plane). The measuring unit is to be specified via
                parameter `units` when calling the class method `create` of the
                `PipeNetwork` or `DuctNetwork` class.
                The height of the nodes is only meaningful when designing or
                analyzing so called "open networks".
        ** 'end_node_ID': text, mandatory
                The ID of the end node of the conduit. The end node is
                the endpoint of the conduit where fluid leaves the conduit.
        - 'end_node_height': number, optional
                See 'start_node_height'
        - 'loop_ID': text, optional
                The ID of the loop to which the conduit belongs. In case the
                conduit belongs to two loops, this must be specified like so:
                (L1, L2). The loop IDs must be placed between parentheses and
                are separated by a comma (just like a tuple in Python).
                The loop IDs must be specified if you want to analyze a pipe or
                duct network.
        - 'zeta': number, optional
                The global resistance coefficient of all fittings and valves
                in the conduit. The zeta-value of a fitting relates the pressure
                drop across the fitting to the volume flow rate through the 
                fitting in the equation `dP = zeta * rho * v**2 / 2` where `dP` 
                is the pressure drop expressed in 'Pa', `rho` is the mass 
                density of the fluid expressed in 'kg / m**3', and `v` is the 
                flow velocity expressed in 'm / s'.
        - 'shape': text, optional {'circular', 'rectangular', or 'flat-oval'}
                The shape of the cross-section. If `shape` is not given a
                circular cross-section is assumed.
        - 'diameter': number, optional
                The inner diameter of the circular conduit. The measuring unit 
                is to be specified via parameter `units` when calling the class 
                method `create` of the `PipeNetwork` or `DuctNetwork` class.
        - 'nominal_diameter': number, optional
                The nominal diameter of the circular conduit. The measuring unit 
                is to be specified via parameter `units` when calling the class 
                method `create` of the `PipeNetwork` or `DuctNetwork` class.
                When the nominal diameter is specified, the pipe or duct 
                schedule must also be specified, either individually for the
                conduit, or globally for all conduits in the network (see more
                on this below in the explanation about column 'schedule').
        - 'width': number, optional
                The width of the cross-section in case of a rectangular or 
                flat-oval conduit. The measuring unit is to be specified via 
                parameter `units` when calling the class method `create` of the 
                `PipeNetwork` or `DuctNetwork` class.
        - 'height': number, optional
                The height of the cross-section in case of a rectangular or 
                flat-oval conduit. The measuring unit is to be specified via 
                parameter `units` when calling the class method `create` of the 
                `PipeNetwork` or `DuctNetwork` class.
        ** 'length': number, mandatory
                The length of the conduit. The measuring unit is to be specified
                via parameter `units` when calling the class method `create` of 
                the `PipeNetwork` or `DuctNetwork` class.
        - 'wall_roughness': number, optional
                The absolute conduit wall roughness. The measuring unit is to be
                specified via parameter `units` when calling the class method 
                `create` of the `PipeNetwork` or `DuctNetwork` class. If the
                same wall roughness applies to all conduits in the network, this
                value can also be specified via the class method `create` of the 
                `PipeNetwork` or `DuctNetwork` class.
        - 'volume_flow_rate': number, optional
                The volume flow rate through the conduit. The measuring unit is 
                to be specified via parameter `units` when calling the class 
                method `create` of the `PipeNetwork` or `DuctNetwork` class.
                If the volume flow rate through the conduits is still unknown,
                this column must be omitted in the csv-file.
                When analyzing the network with the Hardy Cross method, the
                sign of the volume flow rates is coupled to the loop sense,
                which by convention is taken to be clockwise. If the flow sense
                in a conduit is opposite to the loop sense, the volume flow rate
                must be a negative number. In case of a conduit that is common
                to two loops, the sign of the volume flow rate through the
                conduit is referred to the loop first mentioned.
        - 'pressure_drop': number, optional
                The pressure loss across the conduit. The measuring unit is
                to be specified via parameter `units` when calling the class
                method `create` of the `PipeNetwork` or `DuctNetwork` class.
                If the pressure drop across the conduits is still unknown,
                this column must be omitted in the csv-file.
        - 'specific_pressure_drop': number, optional
                The pressure loss per unit length of the conduit. The measuring
                unit is to be specified via parameter `units` when calling the
                class method `create` of the `PipeNetwork` or `DuctNetwork`
                class.
                The specific pressure drop is used by the equal friction method
                to size ducts given the design volume flow rates through the
                ducts.
        - 'machine_coefficients': text, optional
                Pump or fan curves are modeled using a second order polynomial
                `dP = a0 + a1 * V_dot + a2 * V_dot**2`. In case a conduit has
                a pump or fan, the polynomial coefficients (which can be
                determined by curve fitting; see class `PumpCurve` in module
                `utils.py`) can be specified like so: (a0, a1, a2). The
                coefficients must be placed between parentheses and are
                separated by a comma (just like a tuple in Python). The
                coefficients must be ordered from left to right from zero order
                a0 to second order a2.
        - 'fixed_pressure_difference': number, optional
                This only applies to "pseudo conduits", which are used to close
                loops when analyzing open networks with the Hardy Cross method.
                The measuring unit is to be specified via parameter `units` when
                calling the class method `create` of the `PipeNetwork` or
                `DuctNetwork` class.
        - 'schedule': text, optional
                ID of the pipe or duct schedule used to size the pipe or
                duct. This schedule must already be defined before calling the
                method `load_from_csv`. The schedules are held in the
                `PipeScheduleFactory` (in case a `PipeNetwork` instance will
                be created) or in the `DuctScheduleFactory` (in case a
                `DuctNetwork` instance will be created).
                If all the conduits in the network have the same schedule,
                column `schedule` can be omitted and the schedule can be
                assigned when calling the class method `create` of the
                `PipeNetwork` or `DuctNetwork` class.

        Notes
        -----
        The decimal sign used in a csv-file must always be a decimal point,
        never a comma.
        """
        self._clear()
        with open(file_path) as f:
            reader = csv.DictReader(f)
            self._build(reader)

    def load_from_dataframe(self, df: pd.DataFrame):
        """Loads network configuration from a Pandas `DataFrame` object and
        creates the network.
        """
        self._clear()
        self._build(df)

    @staticmethod
    def _extract_loop_id(
        text_item: Optional[str]
    ) -> Optional[Union[str, Tuple[str, ...]]]:
        if text_item is not None:
            text_item = text_item.strip()
            if text_item.startswith('('):
                t = text_item.strip('()')
                t = t.split(',')
                loop_ids = tuple([loop_id.strip() for loop_id in t])
                return loop_ids
            elif len(text_item) > 0:
                loop_id = text_item
                return loop_id
        return None

    @staticmethod
    def _quantify(text_item: Optional[str], unit: str) -> Optional[Quantity]:
        if text_item:
            magnitude = float(text_item)
            return Q_(magnitude, unit)
        return None

    @staticmethod
    def _floatify(text_item: Optional[str]) -> Optional[float]:
        if text_item:
            return float(text_item)
        return None

    def _extract_wall_roughness(self, text_item: Optional[str], unit: str):
        if text_item:
            wall_roughness = Network._quantify(text_item, unit)
            return wall_roughness
        else:
            return self.wall_roughness

    @staticmethod
    def _create_cross_section(
        diameter: Optional[str],
        nominal_diameter: Optional[str],
        width: Optional[str],
        height: Optional[str],
        shape: Optional[str],
        schedule: Optional[TSchedule],
        units: Dict[str, str]
    ) -> TCrossSection:
        if shape == 'circular' or shape is None:
            if diameter is not None:
                return Circular.create(
                    internal_diameter=Q_(float(diameter), units['diameter'])
                )
            elif nominal_diameter is not None and schedule is not None:
                return Circular.create(
                    nominal_diameter=Q_(float(nominal_diameter), units['diameter']),
                    schedule=schedule
                )
            elif schedule is not None:
                return Circular.create(schedule=schedule)
            else:
                return Circular.create()
        elif shape == 'rectangular':
            if width is not None and height is not None:
                return Rectangular.create(
                    height=Q_(float(height), units['dim_height']),
                    width=Q_(float(width), units['dim_width']),
                )
            elif height is not None:
                return Rectangular.create(
                    height=Q_(float(height), units['dim_height']),
                    schedule=schedule
                )
        elif shape == 'flat-oval':
            if width is not None and height is not None:
                return FlatOval.create(
                    height=Q_(float(height), units['dim_height']),
                    width=Q_(float(width), units['dim_width']),
                )
            elif height is not None:
                return FlatOval.create(
                    height=Q_(float(height), units['dim_height']),
                    schedule=schedule
                )

    def _extract_machine_coefficients(
        self,
        text_item: Optional[str]
    ) -> Optional[List[Quantity]]:
        if text_item is not None and text_item.startswith('('):
            t = text_item.strip('()')
            machine_coefficients = t.split(',')
            machine_coefficients = [float(c) for c in machine_coefficients]
            machine_coefficients = [
                Q_(c, u)
                for c, u in zip(machine_coefficients, self.default_units['machine_coefficients'])
            ]
            return machine_coefficients
        return None

    @staticmethod
    def _get_schedule(ID: Optional[str] = None) -> Optional[TSchedule]:
        pass


class PipeNetwork(Network):
    """Class derived from abstract base class `Network` that represents a pipe
    network composed of `Pipe` objects (which are also of class `Conduit`).
    """

    def __init__(self):
        super().__init__()
        self.balancing_valves: Dict[str, BalancingValve] = {}
        self.control_valves: Dict[str, ControlValve] = {}

    def get_pipe_table(self) -> pd.DataFrame:
        """Returns a Pandas `DataFrame` object with an overview of the pipes in
        the network.
        """
        units = self.units
        headers = [
            "pipe ID",
            f"L [{units['length']}]",
            f"DN [{units['diameter']}]",
            f"Di [{units['diameter']}]",
            f"V [{pretty_unit(units['volume_flow_rate'])}]",
            f"v [{pretty_unit(units['velocity'])}]",
            f"zeta",
            "Re",
            f"Δp-dyn. [{units['pressure']}]"
        ]
        table = {header: [] for header in headers}
        for conduit in filter(lambda c: not isinstance(c, PseudoConduit), self.conduits.values()):
            table[headers[0]].append(conduit.ID)
            table[headers[1]].append(conduit.length.to(units['length']).magnitude)
            table[headers[2]].append(conduit.cross_section.nominal_diameter.to(units['diameter']).magnitude)
            table[headers[3]].append(conduit.cross_section.equivalent_diameter.to(units['diameter']).magnitude)
            table[headers[4]].append(conduit.volume_flow_rate.to(units['volume_flow_rate']).magnitude)
            table[headers[5]].append(conduit.velocity.to(units['velocity']).magnitude)
            table[headers[6]].append(conduit.fittings_coefficient)
            table[headers[7]].append(conduit.reynolds_number)
            table[headers[8]].append(conduit.pressure_drop.to(units['pressure']).magnitude)
        return pd.DataFrame(table)

    def add_balancing_valve(
        self,
        cross_over_ID: str,
        pressure_drop_full_open: Quantity = Q_(3, 'kPa')
    ) -> float:
        """Adds a balancing valve to a conduit (cross-over) in the network.

        Parameters
        ----------
        cross_over_ID:
            ID of the conduit in the network to which the balancing valve is
            added. (A cross-over is a conduit that runs between the supply
            and return pipe of a hydronic network to feed e.g. radiator panels.)
        pressure_drop_full_open: default 3 kPa
            Pressure drop across the balancing valve when it is fully open.

        Returns
        -------
        Returns the calculated fully-open valve coefficient (Kvs-value) of the
        balancing valve to attain the given pressure drop.
        """
        balancing_valve = BalancingValve(
            pipe=self.conduits[cross_over_ID],
            ID=f'BAL-VLV-{cross_over_ID}',
            pressure_drop_full_open=pressure_drop_full_open
        )
        self.balancing_valves[cross_over_ID] = balancing_valve
        return balancing_valve.calculate_preliminary_Kvs()

    def set_balancing_valve_Kvs(self, cross_over_ID: str, Kvs: float) -> None:
        """Set the actual fully-open valve coefficient of a commercially
        available balancing valve.
        """
        balancing_valve = self.balancing_valves[cross_over_ID]
        balancing_valve.set_Kvs(Kvs)
        self.conduits[cross_over_ID].add_fitting(balancing_valve.zeta, balancing_valve.ID)

    def add_control_valve(
        self,
        cross_over_ID: str,
        target_authority: float = 0.5
    ) -> float:
        """Adds a control valve to a conduit (cross-over) of the network.

        Parameters
        ----------
        cross_over_ID:
            ID of the conduit in the network to which the control valve is
            added. (A cross-over is a conduit that runs between the supply
            and return pipe of a hydronic network to feed e.g. radiator panels.)
        target_authority: default 0.5
            The control valve authority aimed at. Control valve authority is
            defined as the ratio of the pressure drop across the fully open valve
            to the pressure drop across the fully closed valve. Practically, the
            pressure drop across the fully closed valve is determined as the
            pressure drop across the cross-over when the design volume flow
            rate is flowing through the cross-over.

        Returns
        -------
        Returns the calculated fully-open valve coefficient (Kvs-value) of the
        control valve to attain the given valve authority.
        """
        control_valve = ControlValve(
            pipe=self.conduits[cross_over_ID],
            ID=f'CTRL-VLV-{cross_over_ID}',
            target_authority=target_authority,
            pressure_drop_crit_path=self.critical_path.tot_pressure_diff
        )
        self.control_valves[cross_over_ID] = control_valve
        return control_valve.calculate_preliminary_Kvs()

    def set_control_valve_Kvs(self, cross_over_ID: str, Kvs: float) -> None:
        """Set the actual fully-open valve coefficient of a commercially
        available control valve.
        """
        control_valve = self.control_valves[cross_over_ID]
        control_valve.set_Kvs(Kvs)
        self.conduits[cross_over_ID].add_fitting(control_valve.zeta, control_valve.ID)

    def set_control_valve_opening(
        self,
        cross_over_ID: str,
        percent_open: int
    ) -> None:
        """Sets the control valve's opening position.

        Parameters
        ----------
        cross_over_ID:
            ID of the conduit in the network to which the control valve is
            added. (A cross-over is a conduit that runs between the supply
            and return pipe of a hydronic network to feed e.g. radiator panels.)
        percent_open:
            The opening position of the control valve expressed as a percentage.

        Notes
        -----
        This method can be used to analyze a pipe network at different opening
        positions of a control valve in the network.
        """
        control_valve = self.control_valves[cross_over_ID]
        zeta = control_valve.set_valve_opening(percent_open)
        self.conduits[cross_over_ID].add_fitting(zeta, control_valve.ID)

    def balance_network_at_design(self) -> None:
        """Runs the routine to determine the required, calculated flow
        coefficient settings of the balancing valves in the network so that
        the pressure drop along each flow path of the network becomes equal
        to the pressure drop along the critical path.
        """
        dp_crit = self.critical_path.tot_pressure_diff
        flow_paths = copy(self.flow_paths)
        for balancing_valve in self.balancing_valves.values():
            for i, flow_path in enumerate(flow_paths):
                if balancing_valve.pipe in flow_path:
                    dp_flow_path = flow_path.tot_pressure_diff
                    dp_deficit = dp_crit - dp_flow_path
                    balancing_valve.calculate_Kv_setting(dp_deficit)
                    balancing_valve.pipe.add_fitting(balancing_valve.zeta, balancing_valve.ID)
                    flow_paths.pop(i)
                    break

    def get_balancing_valve_table(self) -> pd.DataFrame:
        """Returns a Pandas `DataFrame` object with an overview of the balancing
        valves in the network.
        """
        table = {
            'pipe ID': [],
            'valve ID': [],
            'Kvs': [],
            'Kvr': []
        }
        for balancing_valve in self.balancing_valves.values():
            table['pipe ID'].append(balancing_valve.pipe.ID)
            table['valve ID'].append(balancing_valve.ID)
            table['Kvs'].append(balancing_valve.Kvs)
            table['Kvr'].append(balancing_valve.Kvr)
        return pd.DataFrame(table)

    def get_control_valve_table(self) -> pd.DataFrame:
        """Returns a Pandas `DataFrame` object with an overview of the control
        valves in the network.
        """
        dp_crit = self.critical_path.tot_pressure_diff
        table = {
            'pipe ID': [],
            'valve ID': [],
            'Kvs': [],
            'authority': []
        }
        for control_valve in self.control_valves.values():
            table['pipe ID'].append(control_valve.pipe.ID)
            table['valve ID'].append(control_valve.ID)
            table['Kvs'].append(control_valve.Kvs)
            table['authority'].append(control_valve.get_valve_authority(dp_crit))
        return pd.DataFrame(table)

    @staticmethod
    def _get_schedule(name: Optional[str] = None) -> Optional[PipeSchedule]:
        """Get the pipe schedule identified by `name` from the
        `PipeScheduleFactory`.
        """
        if name is None:
            return None
        else:
            return PipeScheduleFactory.get(name)


class DuctNetwork(Network):
    """Class derived from abstract base class `Network` that represents a duct
    network composed of `Duct` objects (which are also of class `Conduit`).
    """

    def __init__(self):
        super().__init__()
        self.balancing_dampers: Dict[str, TVolumeDamper] = {}

    def get_duct_table(self) -> pd.DataFrame:
        """Returns a Pandas `DataFrame` object with an overview of the ducts in
        the network.
        """
        units = self.units
        headers = [
            "duct ID",
            f"L [{units['length']}]",
            f"Deq [{units['diameter']}]",
            f"width [{units['dim_width']}]",
            f"height [{units['dim_height']}]",
            f"V [{pretty_unit(units['volume_flow_rate'])}]",
            f"v [{pretty_unit(units['velocity'])}]",
            "Re",
            f"Δp-dyn. [{units['pressure']}]"
        ]
        table = {header: [] for header in headers}
        conduits = (
            conduit for conduit in self.conduits.values()
            if not isinstance(conduit, PseudoConduit)
        )
        for conduit in conduits:
            table[headers[0]].append(conduit.ID)
            table[headers[1]].append(conduit.length.to(units['length']).magnitude)
            table[headers[2]].append(conduit.cross_section.equivalent_diameter.to(units['diameter']).magnitude)
            if isinstance(conduit.cross_section, (Rectangular, FlatOval)):
                table[headers[3]].append(conduit.cross_section.width.to(units['dim_width']).magnitude)
                table[headers[4]].append(conduit.cross_section.height.to(units['dim_height']).magnitude)
            else:
                table[headers[3]].append(float('nan'))
                table[headers[4]].append(float('nan'))
            if conduit.flow_sign == FlowSign.CLOCKWISE:
                table[headers[5]].append(conduit.volume_flow_rate.to(units['volume_flow_rate']).magnitude)
            else:
                table[headers[5]].append(-conduit.volume_flow_rate.to(units['volume_flow_rate']).magnitude)
            table[headers[6]].append(conduit.velocity.to(units['velocity']).magnitude)
            table[headers[7]].append(conduit.reynolds_number)
            table[headers[8]].append(conduit.pressure_drop.to(units['pressure']).magnitude)
        return pd.DataFrame(table)

    def add_balancing_damper(
        self,
        duct_ID: str,
        damper_type: str = 'A15A'
    ) -> None:
        """Adds a fully open balancing damper to a duct (branch duct) in the
        network.

        Parameters
        ----------
        duct_ID:
            ID of the conduit in the network to which the balancing damper is
            added.
        damper_type: default 'A15A'
            The type of balancing damper.
            At this stage, only 1 type of balancing damper has been implemented
            in the program. 'A15A' refers to the table in appendix A of the
            SMACNA guide "HVAC SYSTEMS DUCT DESIGN".

        Returns
        -------
        Returns the calculated fully-open valve coefficient (Kvs-value) of the
        balancing valve to attain the given pressure drop.
        """
        duct = self.conduits[duct_ID]
        if damper_type == 'A15B':
            volume_damper = None  # TODO: extend with other types of volume damper
        else:
            volume_damper = ObstructionA15A(
                duct,
                theta=Q_(0.0, 'deg'),
                ID=f'damper_{duct.ID}'
            )
        self.balancing_dampers[volume_damper.ID] = volume_damper
        duct.add_fitting(zeta=volume_damper.zeta, ID=volume_damper.ID)

    def balance_network_at_design(self) -> None:
        """Runs the routine to determine the required, calculated angle settings
        of the balancing dampers in the network so that the pressure drop along
        each flow path of the network becomes equal to the pressure drop along
        the critical path.
        """
        dp_crit = self.critical_path.tot_pressure_diff
        flow_paths = copy(self.flow_paths)
        for balancing_damper in self.balancing_dampers.values():
            for i, flow_path in enumerate(flow_paths):
                if balancing_damper.duct in flow_path:
                    dp_flow_path = flow_path.tot_pressure_diff
                    dp_deficit = dp_crit - dp_flow_path
                    dp_tot = balancing_damper.pressure_drop + dp_deficit
                    dp_tot = dp_tot.to('Pa').magnitude
                    pv = balancing_damper.duct.velocity_pressure.to('Pa').magnitude
                    balancing_damper.zeta = dp_tot / pv
                    balancing_damper.duct.add_fitting(balancing_damper.zeta, balancing_damper.ID)
                    flow_paths.pop(i)
                    break

    def get_balancing_damper_table(self) -> pd.DataFrame:
        """Returns a Pandas `DataFrame` object with an overview of the balancing
        dampers in the network and their calculated settings.
        """
        table = {
            'volume damper ID': [],
            'set angle': [],
            'zeta': []
        }
        for damper in self.balancing_dampers.values():
            table['volume damper ID'].append(damper.ID)
            table['set angle'].append(damper.theta)
            table['zeta'].append(damper.zeta)
        return pd.DataFrame(table)

    @staticmethod
    def _get_schedule(name: Optional[str] = None) -> Optional[DuctSchedule]:
        """Get the duct schedule identified by `name` from the
        `DuctScheduleFactory`.
        """
        if name is None:
            return None
        else:
            return DuctScheduleFactory.get(name)


def save_network(
    network: Union[PipeNetwork, DuctNetwork],
    file_path: str
) -> None:
    """Saves a `Network` object to disk file. For example, after designing a
     network, you can save it to disk, to open it in another script for
     analyzing the network.
    """
    with open(file_path, 'wb') as fh:
        pickle.dump(network, fh)


def load_network(file_path: str) -> Union[PipeNetwork, DuctNetwork]:
    """Loads a `Network` object from a disk file."""
    with open(file_path, 'rb') as fh:
        network = pickle.load(fh)
    return network
