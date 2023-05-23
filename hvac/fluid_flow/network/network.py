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


def pretty_unit(unit: str) -> str:
    q = Quantity(0, unit)
    q = f"{q:~P}".split(' ')
    return q[1]


class FlowPath(List[TConduit]):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.ambient_air: FluidState = STANDARD_AIR
        self._conduits: List[Conduit] = []

    def __str__(self):
        return '|'.join(conduit.ID for conduit in self.conduits)

    @property
    def conduits(self) -> List[Conduit]:
        if not self._conduits:
            self._conduits = [conduit for conduit in self if not isinstance(conduit, PseudoConduit)]
        return self._conduits

    @property
    def dynamic_pressure_difference(self) -> Quantity:
        return sum(conduit.pressure_drop for conduit in self.conduits)

    @property
    def elevation_pressure_difference(self) -> Quantity:
        # aka thermal gravity effect or chimney effect
        first_node = self.conduits[0].start_node
        last_node = self.conduits[-1].end_node
        z1 = first_node.height.to('m')
        z2 = last_node.height.to('m')
        rho = self.conduits[0].fluid.rho.to('kg / m ** 3')
        g = Q_(9.81, 'm / s ** 2')
        rho_a = self.ambient_air.rho.to('kg / m ** 3')
        return g * (rho - rho_a) * (z2 - z1)

    @property
    def total_pressure_difference(self) -> Quantity:
        dp_elev = self.elevation_pressure_difference
        dp_dyn = self.dynamic_pressure_difference
        dp_tot = dp_elev + dp_dyn  # dp_tot = (tp1 - tp2) + dp_pump
        return dp_tot

    def get_system_curve(
            self,
            V_wp: Quantity,
            dp_wp: Quantity,
            p_g1: Optional[Quantity] = Q_(0.0, 'Pa'),
    ) -> SystemCurve:
        dp_dyn = self.dynamic_pressure_difference
        dp_elev = self.elevation_pressure_difference
        dp_loss = dp_dyn + dp_wp
        R_hyd = dp_loss / (V_wp ** 2)
        p_v1 = self.conduits[0].velocity_pressure
        p_t1 = p_g1 + p_v1
        p_v2 = self.conduits[-1].velocity_pressure
        p_t2 = p_t1 - dp_elev - dp_dyn
        p_g2 = p_t2 - p_v2
        dp_g = p_g2 - p_g1
        dp_v = p_v2 - p_v1
        dp_t = dp_g + dp_v
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
        self.end_node: Optional[Node] = None
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
            end_node_ID: Optional[str] = None,
            end_node_height: Quantity = Q_(0.0, 'm'),
            units: Optional[Dict[str, str]] = None
    ) -> 'Network':
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
    ):
        conduit.ID = conduit_ID
        if isinstance(zeta, float) and isinstance(conduit, Conduit):
            conduit.add_fitting(zeta)
        sn = self.nodes.setdefault(start_node_ID, Node(start_node_ID, start_node_height))
        en = self.nodes.setdefault(end_node_ID, Node(end_node_ID, end_node_height))
        sn.connect(conduit, NodeArrow.OUTGOING)
        en.connect(conduit, NodeArrow.INCOMING)
        conduit.start_node = sn
        conduit.end_node = en
        conduit.loop_ID = loop_ID
        self.conduits[conduit_ID] = conduit

    def _create_loops(self):
        # (re)create the loops for analysis with Hardy Cross method
        self.loops = {}
        for conduit in self.conduits.values():
            loop_ID = conduit.loop_ID
            if isinstance(loop_ID, tuple) and len(loop_ID) == 2:
                loop = self.loops.setdefault(loop_ID[0], Loop(loop_ID[0]))
                other_loop = self.loops.setdefault(loop_ID[1], Loop(loop_ID[1]))

                conduit.loops = []
                conduit_duplicate = Conduit.duplicate(conduit, flow_sign=conduit.flow_sign.reverse())

                loop.append(conduit)
                other_loop.append(conduit_duplicate)

                conduit.loops = [self.loops[loop_ID[0]], self.loops[loop_ID[1]]]
                conduit_duplicate.loops = [self.loops[loop_ID[1]], self.loops[loop_ID[0]]]

            if isinstance(loop_ID, str) and len(loop_ID) > 0:
                loop = self.loops.setdefault(loop_ID, Loop(loop_ID))
                loop.append(conduit)
                conduit.loops = [self.loops[loop_ID]]

    def _hardy_cross_algorithm(self):
        # calculate the loop correction terms of all loops in the network
        for loop in self.loops.values():
            loop.calculate_correction_term()

        # apply the correction terms to the flow rate of each section in each loop
        for loop in self.loops.values():
            for conduit in loop:
                # correction terms cannot be applied to pseudo-sections as they have no flow
                if not isinstance(conduit, PseudoConduit):
                    # check if the section is common to 2 loops
                    if len(conduit.loops) == 2:
                        # if a section belongs to 2 loops: get the correction term of the other loop
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

    def analyze(self, tolerance: Quantity = Quantity(1.0, 'Pa'), i_max: int = 100) -> int:
        self._create_loops()
        i = 0
        while not self._check_loops(tolerance):
            self._hardy_cross_algorithm()
            i += 1
            if i > i_max:
                raise OverflowError(f'no solution after {i_max} iterations')
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

    def _search_flow_paths(self):
        node = self.start_node
        path = FlowPath()
        path.ambient_air = self.ambient_air
        self._flow_paths.append(path)
        self._recursive_path_search(node, path)

    def _recursive_path_search(self, node: Node, path: FlowPath):
        while True:
            if len(node.outgoing) > 1:
                for conduit in node.outgoing[1:]:
                    if conduit not in path:
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
                if conduit not in path:
                    path.append(conduit)
                    node = self.nodes[conduit.end_node.ID]
                else:
                    return

    @property
    def flow_paths(self) -> List[FlowPath]:
        if not self._flow_paths:
            self._search_flow_paths()
        return self._flow_paths

    @property
    def critical_path(self) -> FlowPath:
        dp_tots = [fp.total_pressure_difference for fp in self.flow_paths]
        index_max = dp_tots.index(max(dp_tots))
        return self._flow_paths[index_max]

    @property
    def volume_flow_rate(self) -> Quantity:
        return sum(conduit.volume_flow_rate for conduit in self.start_node.outgoing)

    @property
    def total_pressure_difference(self) -> Quantity:
        return self.critical_path.total_pressure_difference

    @property
    def hydraulic_resistance(self) -> Quantity:
        dp_loss = self.critical_path.dynamic_pressure_difference.to('Pa')
        V = self.volume_flow_rate.to('m ** 3 / s')
        return dp_loss / (V ** 2)

    def get_flow_path_table(self) -> pd.DataFrame:
        headers = [
            "path",
            f"Δp-elev. [{self.units['pressure']}]",
            f"Δp-dyn. [{self.units['pressure']}]",
            f"Δp-tot. [{self.units['pressure']}]",
            f"Δp-deficit [{self.units['pressure']}]"
        ]
        table = {header: [] for header in headers}
        for flow_path in self.flow_paths:
            table[headers[0]].append(str(flow_path))
            table[headers[1]].append(flow_path.elevation_pressure_difference.to(self.units['pressure']).magnitude)
            table[headers[2]].append(flow_path.dynamic_pressure_difference.to(self.units['pressure']).magnitude)
            table[headers[3]].append(flow_path.total_pressure_difference.to(self.units['pressure']).magnitude)
            dp_deficit = self.total_pressure_difference - flow_path.total_pressure_difference
            table[headers[4]].append(dp_deficit.to(self.units['pressure']).magnitude)
        return pd.DataFrame(table)

    def get_fitting_table(self) -> pd.DataFrame:
        fitting_tables = [
            conduit.get_fitting_table(self.units['pressure'])
            for conduit in self.conduits.values()
            if not isinstance(conduit, PseudoConduit)
        ]
        fitting_table = pd.concat(fitting_tables, ignore_index=True)
        return fitting_table

    def _clear(self):
        self.nodes = {}
        self.start_node = Node(self.start_node.ID, self.start_node.height)
        self.nodes[self.start_node.ID] = self.start_node
        if self.end_node is not None:
            self.end_node = Node(self.end_node.ID, self.end_node.height)
            self.nodes[self.end_node.ID] = self.end_node
        self.conduits = {}
        self._flow_paths = []
        self.loops = {}

    def _build(self, data: Union[csv.DictReader, pd.DataFrame]):
        if isinstance(data, pd.DataFrame):
            data = data.to_dict(orient='records')
        for row in data:
            if not row.get('fixed_pressure_difference'):
                cross_section = self._create_cross_section(
                    diameter=row.get('diameter'),
                    nominal_diameter=row.get('nominal_diameter'),
                    width=row.get('width'),
                    height=row.get('height'),
                    shape=row.get('shape'),
                    schedule=self._get_schedule(row.get('schedule')) or self.schedule,
                    units=self.units
                )
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
                conduit = PseudoConduit.create(
                    fixed_pressure_drop=self._quantify(
                        row['fixed_pressure_difference'],
                        self.units['pressure']
                    )
                )
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
        Load a network configuration from a csv-file.

        Available column titles:
        - 'conduit_ID': mandatory
        - 'start_node_ID': mandatory
        - 'start_node_height': optional
        - 'end_node_ID': mandatory
        - 'end_node_height': optional
        - 'loop_ID': optional
        - 'zeta': optional
        - 'shape': optional (possible values: 'circular', 'rectangular', or 'flat-oval')
        - 'diameter': optional
        - 'nominal_diameter': optional
        - 'width': optional
        - 'height': optional
        - 'length': mandatory
        - 'wall_roughness': optional
        - 'volume_flow_rate': optional
        - 'pressure_drop': optional
        - 'specific_pressure_drop': optional
        - 'machine_coefficients': optional
        - 'fixed_pressure_difference': optional
        - 'schedule': optional
        """
        self._clear()
        with open(file_path) as f:
            reader = csv.DictReader(f)
            self._build(reader)

    def load_from_dataframe(self, df: pd.DataFrame):
        self._clear()
        self._build(df)

    @staticmethod
    def _extract_loop_id(text_item: Optional[str]) -> Optional[Union[str, Tuple[str, ...]]]:
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

    def _extract_machine_coefficients(self, text_item: Optional[str]) -> Optional[List[Quantity]]:
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

    def __init__(self):
        super().__init__()
        self.balancing_valves: Dict[str, BalancingValve] = {}
        self.control_valves: Dict[str, ControlValve] = {}

    def get_pipe_table(self) -> pd.DataFrame:
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

    def add_balancing_valve(self, cross_over_ID: str, pressure_drop_full_open: Quantity = Q_(3, 'kPa')) -> float:
        balancing_valve = BalancingValve(
            pipe=self.conduits[cross_over_ID],
            ID=f'BAL-VLV-{cross_over_ID}',
            pressure_drop_full_open=pressure_drop_full_open
        )
        self.balancing_valves[cross_over_ID] = balancing_valve
        return balancing_valve.calculate_preliminary_Kvs()

    def set_balancing_valve_Kvs(self, cross_over_ID: str, Kvs: float) -> None:
        balancing_valve = self.balancing_valves[cross_over_ID]
        balancing_valve.set_Kvs(Kvs)
        self.conduits[cross_over_ID].add_fitting(balancing_valve.zeta, balancing_valve.ID)

    def add_control_valve(self, cross_over_ID: str, target_authority: float = 0.5) -> float:
        control_valve = ControlValve(
            pipe=self.conduits[cross_over_ID],
            ID=f'CTRL-VLV-{cross_over_ID}',
            target_authority=target_authority,
            pressure_drop_crit_path=self.critical_path.total_pressure_difference
        )
        self.control_valves[cross_over_ID] = control_valve
        return control_valve.calculate_preliminary_Kvs()

    def set_control_valve_Kvs(self, cross_over_ID: str, Kvs: float) -> None:
        control_valve = self.control_valves[cross_over_ID]
        control_valve.set_Kvs(Kvs)
        self.conduits[cross_over_ID].add_fitting(control_valve.zeta, control_valve.ID)

    def set_control_valve_opening(self, cross_over_ID: str, percent_open: int) -> None:
        control_valve = self.control_valves[cross_over_ID]
        zeta = control_valve.set_valve_opening(percent_open)
        self.conduits[cross_over_ID].add_fitting(zeta, control_valve.ID)

    def balance_network_at_design(self) -> None:
        dp_crit = self.critical_path.total_pressure_difference
        flow_paths = copy(self.flow_paths)
        for balancing_valve in self.balancing_valves.values():
            for i, flow_path in enumerate(flow_paths):
                if balancing_valve.pipe in flow_path:
                    dp_flow_path = flow_path.total_pressure_difference
                    dp_deficit = dp_crit - dp_flow_path
                    balancing_valve.calculate_Kv_setting(dp_deficit)
                    balancing_valve.pipe.add_fitting(balancing_valve.zeta, balancing_valve.ID)
                    flow_paths.pop(i)
                    break

    def get_balancing_valve_table(self) -> pd.DataFrame:
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
        dp_crit = self.critical_path.total_pressure_difference
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
        if name is None:
            return None
        else:
            return PipeScheduleFactory.get(name)


class DuctNetwork(Network):

    def __init__(self):
        super().__init__()
        self.balancing_dampers: Dict[str, TVolumeDamper] = {}

    def get_duct_table(self) -> pd.DataFrame:
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
        conduits = (conduit for conduit in self.conduits.values() if not isinstance(conduit, PseudoConduit))
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

    def add_balancing_damper(self, duct_ID: str, damper_type: str = 'A15A') -> None:
        duct = self.conduits[duct_ID]
        if damper_type == 'A15B':
            volume_damper = None  # TODO: extend with other types of volume damper
        else:
            volume_damper = ObstructionA15A(duct, theta=Q_(0.0, 'deg'), ID=f'damper_{duct.ID}')
        self.balancing_dampers[volume_damper.ID] = volume_damper
        duct.add_fitting(zeta=volume_damper.zeta, ID=volume_damper.ID)

    def balance_network_at_design(self) -> None:
        dp_crit = self.critical_path.total_pressure_difference
        flow_paths = copy(self.flow_paths)
        for balancing_damper in self.balancing_dampers.values():
            for i, flow_path in enumerate(flow_paths):
                if balancing_damper.duct in flow_path:
                    dp_flow_path = flow_path.total_pressure_difference
                    dp_deficit = dp_crit - dp_flow_path
                    dp_tot = balancing_damper.pressure_drop + dp_deficit
                    dp_tot = dp_tot.to('Pa').magnitude
                    pv = balancing_damper.duct.velocity_pressure.to('Pa').magnitude
                    balancing_damper.zeta = dp_tot / pv
                    balancing_damper.duct.add_fitting(balancing_damper.zeta, balancing_damper.ID)
                    flow_paths.pop(i)
                    break

    def get_balancing_damper_table(self) -> pd.DataFrame:
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
        if name is None:
            return None
        else:
            return DuctScheduleFactory.get(name)


def save_network(network: Union[PipeNetwork, DuctNetwork], file_path: str) -> None:
    with open(file_path, 'wb') as fh:
        pickle.dump(network, fh)


def load_network(file_path: str) -> Union[PipeNetwork, DuctNetwork]:
    with open(file_path, 'rb') as fh:
        network = pickle.load(fh)
    return network
