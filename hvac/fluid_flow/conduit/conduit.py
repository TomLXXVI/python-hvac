from typing import Optional, List, Dict, Union, Tuple
from abc import ABC, abstractmethod
from copy import deepcopy
from dataclasses import dataclass
from enum import IntEnum
import math
import pandas as pd
from hvac import Quantity, UNITS
from hvac.fluids import Fluid, FluidState
from .cross_section import TCrossSection
from . import exceptions

u = UNITS
Water = Fluid('Water')


@dataclass
class _Fitting:
    ID: str = ''
    zeta: float = float('nan')


class NodeArrow(IntEnum):
    INCOMING = 1
    OUTGOING = -1


class FlowSign(IntEnum):
    CLOCKWISE = 1
    COUNTERCLOCKWISE = -1

    def reverse(self) -> 'FlowSign':
        if self.value == FlowSign.CLOCKWISE:
            return FlowSign.COUNTERCLOCKWISE
        else:
            return FlowSign.CLOCKWISE


class Node:

    def __init__(self, ID: str = '', height: Quantity = 0 * u.m):
        self.ID = ID
        self.height = height
        self.incoming: List['Conduit'] = []
        self.outgoing: List['Conduit'] = []
        self.conduit_IDs: List[str] = []

    def connect(self, conduit: 'Conduit', node_arrow: NodeArrow):
        if conduit.ID not in self.conduit_IDs:
            self.conduit_IDs.append(conduit.ID)
            if node_arrow == NodeArrow.INCOMING:
                self.incoming.append(conduit)
            if node_arrow == NodeArrow.OUTGOING:
                self.outgoing.append(conduit)

    @property
    def volume_flow_rate(self) -> Quantity:
        V_in = sum(conduit.volume_flow_rate for conduit in self.incoming if not isinstance(conduit, PseudoConduit))
        V_out = sum(conduit.volume_flow_rate for conduit in self.outgoing if not isinstance(conduit, PseudoConduit))
        return V_out - V_in


class AbstractConduit(ABC):

    def __init__(self):
        self.ID: str = ''
        self.loop_ID: Optional[Union[str, Tuple[str, str]]] = None
        self.loops: List['Loop'] = []
        self.start_node: Optional[Node] = None
        self.end_node: Optional[Node] = None
        self.flow_sign: FlowSign = FlowSign.CLOCKWISE

    @classmethod
    @abstractmethod
    def create(cls, *args, **kwargs) -> 'AbstractConduit':
        ...

    @staticmethod
    @abstractmethod
    def duplicate(conduit: 'AbstractConduit', **kwargs) -> 'AbstractConduit':
        ...

    @property
    @abstractmethod
    def numerator(self) -> Quantity:
        ...

    @property
    @abstractmethod
    def denominator(self) -> Quantity:
        ...


class PseudoConduit(AbstractConduit):
    # only used for analysis

    def __init__(self):
        super().__init__()
        self._fixed_pressure_drop: Quantity = float('nan') * u.Pa

    @classmethod
    def create(cls, fixed_pressure_drop: Quantity) -> 'PseudoConduit':
        obj = cls()
        obj._fixed_pressure_drop = fixed_pressure_drop
        return obj

    @staticmethod
    def duplicate(conduit: 'PseudoConduit', **kwargs) -> 'PseudoConduit':
        obj = deepcopy(conduit)
        obj._fixed_pressure_drop = kwargs.get('fixed_pressure_drop', obj._fixed_pressure_drop)
        obj.ID = kwargs.get('ID', obj.ID)
        obj.start_node = kwargs.get('start_node', obj.start_node)
        obj.end_node = kwargs.get('end_node', obj.end_node)
        obj.flow_sign = kwargs.get('flow_sign', obj.flow_sign)
        obj.loops = kwargs.get('loops', obj.loops)
        return obj

    @property
    def numerator(self) -> Quantity:
        return self._fixed_pressure_drop.to('Pa')

    @property
    def denominator(self) -> Quantity:
        return Quantity(0.0, 'Pa * s / m ** 3')

    @property
    def pressure_drop(self) -> Quantity:
        return self._fixed_pressure_drop


class Conduit(AbstractConduit):
    MAX_ITERATIONS = 30
    SOLUTION_TOLERANCE = 1e-5

    def __init__(self):
        super().__init__()
        self.length: Quantity = float('nan') * u.m
        self.wall_roughness: Quantity = float('nan') * u.mm
        self.fluid: FluidState = Water(T=Quantity(15.0, 'degC'), P=Quantity(101.325, 'kPa'))
        self._cross_section: Optional[TCrossSection] = None
        self._volume_flow_rate: Quantity = Quantity(float('nan'), 'm ** 3 / s')
        self._pressure_drop: Quantity = float('nan') * u.Pa
        self._fittings: Dict[str, _Fitting] = {}
        self.machine_coefficients: Optional[List[Quantity]] = None  # pump or fan curve coefficients

    @classmethod
    def create(
        cls,
        length: Quantity,
        wall_roughness: Quantity,
        fluid: FluidState,
        cross_section: TCrossSection,
        volume_flow_rate: Optional[Quantity] = None,
        pressure_drop: Optional[Quantity] = None,
        specific_pressure_drop: Optional[Quantity] = None,
        machine_coefficients: Optional[List[Quantity]] = None
    ) -> 'Conduit':
        obj = cls()
        obj.length = length
        obj.wall_roughness = wall_roughness
        obj.fluid = fluid
        obj._cross_section = cross_section
        if (volume_flow_rate is not None) and (volume_flow_rate < 0.0):
            obj._volume_flow_rate = abs(volume_flow_rate)
            obj.flow_sign = FlowSign.COUNTERCLOCKWISE
        else:
            obj._volume_flow_rate = volume_flow_rate
            obj.flow_sign = FlowSign.CLOCKWISE
        # obj._specific_pressure_drop = specific_pressure_drop
        if isinstance(specific_pressure_drop, Quantity):
            obj._pressure_drop = specific_pressure_drop * length
        else:
            obj._pressure_drop = pressure_drop
        obj.machine_coefficients = machine_coefficients
        flow_rate_given = isinstance(obj._volume_flow_rate, Quantity)
        pressure_drop_given = isinstance(obj._pressure_drop, Quantity)
        diameter_given = not math.isnan(obj._cross_section.hydraulic_diameter.magnitude)
        if flow_rate_given and pressure_drop_given:
            obj._calculate_equivalent_diameter()
        elif pressure_drop_given and diameter_given:
            obj._calculate_volume_flow_rate()
        elif flow_rate_given and diameter_given:
            obj._calculate_pressure_drop()
        else:
            raise exceptions.ConduitConfigurationError('arguments missing')
        return obj

    # @staticmethod
    # def duplicate(conduit: 'Conduit', **kwargs) -> 'Conduit':
    #     obj = deepcopy(conduit) --> throws TypeError as deepcopy fails on CoolProp.CoolProp.AbstractState
    #     obj.length = kwargs.get('length', obj.length)
    #     obj.wall_roughness = kwargs.get('wall_roughness', obj.wall_roughness)
    #     obj.fluid = kwargs.get('fluid', obj.fluid)
    #     obj._cross_section = kwargs.get('cross_section', obj._cross_section)
    #     obj._volume_flow_rate = kwargs.get('volume_flow_rate', obj._volume_flow_rate)
    #     obj._pressure_drop = kwargs.get('pressure_drop', obj._pressure_drop)
    #     if (specific_pressure_drop := kwargs.get('specific_pressure_drop')) is not None:
    #         obj._pressure_drop = specific_pressure_drop * obj.length
    #     obj._fittings = kwargs.get('fittings', obj._fittings)
    #     obj.ID = kwargs.get('ID', obj.ID)
    #     obj.start_node = kwargs.get('start_node', obj.start_node)
    #     obj.end_node = kwargs.get('end_node', obj.end_node)
    #     obj.flow_sign = kwargs.get('flow_sign', obj.flow_sign)
    #     obj.loops = kwargs.get('loops', obj.loops)
    #     obj.machine_coefficients = kwargs.get('machine_coefficients', obj.machine_coefficients)
    #     return obj

    @staticmethod
    def duplicate(conduit: 'Conduit', **kwargs) -> 'Conduit':
        obj = Conduit.create(
            length=deepcopy(conduit.length),
            wall_roughness=deepcopy(conduit.wall_roughness),
            fluid=conduit.fluid.fluid(T=conduit.fluid.T, P=conduit.fluid.P),
            cross_section=deepcopy(conduit.cross_section),
            volume_flow_rate=deepcopy(conduit.volume_flow_rate),
            pressure_drop=deepcopy(conduit._pressure_drop),
            specific_pressure_drop=None,
            machine_coefficients=deepcopy(conduit.machine_coefficients)
        )
        obj.length = kwargs.get('length', obj.length)
        obj.wall_roughness = kwargs.get('wall_roughness', obj.wall_roughness)
        obj.fluid = kwargs.get('fluid', obj.fluid)
        obj._cross_section = kwargs.get('cross_section', obj._cross_section)
        obj._volume_flow_rate = kwargs.get('volume_flow_rate', obj._volume_flow_rate)
        obj._pressure_drop = kwargs.get('pressure_drop', obj._pressure_drop)
        if (specific_pressure_drop := kwargs.get('specific_pressure_drop')) is not None:
            obj._pressure_drop = specific_pressure_drop * obj.length
        obj._fittings = kwargs.get('fittings', deepcopy(conduit._fittings))
        obj.ID = kwargs.get('ID', deepcopy(conduit.ID))
        obj.start_node = kwargs.get("start_node", conduit.start_node)
        obj.end_node = kwargs.get("end_node", conduit.end_node)
        obj.flow_sign = kwargs.get('flow_sign', deepcopy(conduit.flow_sign))
        obj.loops = kwargs.get('loops', deepcopy(conduit.loops))
        obj.machine_coefficients = kwargs.get('machine_coefficients', obj.machine_coefficients)
        return obj

    def _calculate_darcy_friction_factor(self, Re: float, e: float) -> float:
        f = self._serghide(Re, e)
        return f

    @staticmethod
    def _churchill(Re: float, e: float) -> float:
        d = (7 / Re) ** 0.9 + 0.27 * e
        A = (2.457 * math.log(1 / d)) ** 16
        B = (37.530 / Re) ** 16
        f = 8 * ((8 / Re) ** 12 + 1 / (A + B) ** 1.5) ** (1 / 12)
        return f

    @staticmethod
    def _serghide(Re: float, e: float) -> float:
        # equation of Serghide
        f1 = -2.0 * math.log10(e / 3.7 + 12.0 / Re)
        f2 = -2.0 * math.log10(e / 3.7 + 2.51 * f1 / Re)
        f3 = -2.0 * math.log10(e / 3.7 + 2.51 * f2 / Re)
        f = (f1 - (f2 - f1) ** 2.0 / (f3 - 2.0 * f2 + f1)) ** -2.0
        return f

    def _calculate_equivalent_diameter(self):
        l = self.length.to('m').magnitude
        dp = self._pressure_drop.to('Pa').magnitude  # only friction pressure drop
        V = self._volume_flow_rate.to('m ** 3 / s').magnitude
        rho = self.fluid.rho.to('kg / m ** 3').magnitude
        mu = self.fluid.mu.to('Pa * s').magnitude
        e_abs = self.wall_roughness.to('m').magnitude
        f = 0.03
        i = 0
        while i <= self.MAX_ITERATIONS:
            d_eq = (f * (l / dp) * rho * (8.0 / math.pi ** 2) * V ** 2) ** (1 / 5)
            A = math.pi * d_eq ** 2 / 4
            v = V / A
            Re = v * d_eq * rho / mu
            e = e_abs / d_eq
            f_new = self._calculate_darcy_friction_factor(Re, e)
            if abs(f_new - f) <= self.SOLUTION_TOLERANCE:
                break
            f = f_new
            i += 1
        else:
            raise exceptions.ConduitIterationOverflow('no solution for diameter')
        self._cross_section.equivalent_diameter = d_eq * u.m  # this assignment will change the cross-section's size
        if self._cross_section.schedule:  # a commercially available size will be determined
            self._calculate_pressure_drop()  # a change of the cross-section's size, will alter the pressure drop

    def _calculate_volume_flow_rate(self):
        l = self.length.to('m').magnitude
        dp = self._pressure_drop.to('Pa').magnitude
        rho = self.fluid.rho.to('kg / m ** 3').magnitude
        mu = self.fluid.mu.to('Pa * s').magnitude
        zeta = self.fittings_coefficient
        e_abs = self.wall_roughness.to('m').magnitude
        dh = self._cross_section.hydraulic_diameter.to('m').magnitude
        e = e_abs / dh
        f = 1.325 / math.log10(e / 3.7) ** 2.0
        i = 0
        while i <= self.MAX_ITERATIONS:
            v = math.sqrt(2 * dp / ((f * (l / dh) + zeta) * rho))
            Re = v * dh * rho / mu
            f_new = self._calculate_darcy_friction_factor(Re, e)
            if abs(f_new - f) <= self.SOLUTION_TOLERANCE:
                break
            f = f_new
            i += 1
        else:
            raise exceptions.ConduitIterationOverflow('no solution for volume flow rate')
        A = self._cross_section.area.to('m ** 2').magnitude
        self._volume_flow_rate = Quantity(v * A, 'm ** 3 / s')

    def _calculate_pressure_drop(self):
        l = self.length.to('m').magnitude
        rho = self.fluid.rho.to('kg / m ** 3').magnitude
        mu = self.fluid.mu.to('Pa * s').magnitude
        zeta = self.fittings_coefficient
        e_abs = self.wall_roughness.to('m').magnitude
        dh = self._cross_section.hydraulic_diameter.to('m').magnitude
        e = e_abs / dh
        V = self._volume_flow_rate.to('m ** 3 / s').magnitude
        A = self._cross_section.area.to('m ** 2').magnitude
        v = V / A
        Re = v * dh * rho / mu
        f = self._calculate_darcy_friction_factor(Re, e)
        dp = (f * (l / dh) + zeta) * rho * v ** 2 / 2
        self._pressure_drop = dp * u.Pa

    @property
    def volume_flow_rate(self) -> Quantity:
        return self._volume_flow_rate

    @volume_flow_rate.setter
    def volume_flow_rate(self, value: Quantity):
        # when the flow rate through a conduit changes, the pressure drop across
        # the conduit will also change
        if value < 0.0:
            self._volume_flow_rate = abs(value)
            self.flow_sign = self.flow_sign.reverse()
        else:
            self._volume_flow_rate = value
        self._calculate_pressure_drop()

    @property
    def pressure_drop(self) -> Quantity:
        dp = self._pressure_drop
        if self.machine_coefficients is not None:
            V = self.volume_flow_rate
            m = self.machine_coefficients
            dp_machine = sum(m[i] * V ** i for i in range(len(m)))
            dp -= dp_machine
        return dp

    @pressure_drop.setter
    def pressure_drop(self, value: Quantity):
        # when the pressure across a conduit changes, the flow rate through the
        # conduit will also change
        self._pressure_drop = value
        self._calculate_volume_flow_rate()

    @property
    def specific_pressure_drop(self) -> Quantity:
        return self._pressure_drop / self.length

    @specific_pressure_drop.setter
    def specific_pressure_drop(self, value: Quantity):
        self.pressure_drop = value * self.length

    @property
    def cross_section(self) -> TCrossSection:
        return self._cross_section

    @cross_section.setter
    def cross_section(self, value: TCrossSection):
        # when the cross-section of a conduit changes, the pressure drop across
        # the conduit will also change for a given network flow rate.
        self._cross_section = value
        self._calculate_pressure_drop()

    @property
    def velocity(self) -> Quantity:
        return self._volume_flow_rate / self._cross_section.area

    @property
    def reynolds_number(self) -> float:
        rho = self.fluid.rho.to('kg / m ** 3').magnitude
        mu = self.fluid.mu.to('Pa * s').magnitude
        dh = self._cross_section.hydraulic_diameter.to('m').magnitude
        V = self._volume_flow_rate.to('m ** 3 / s').magnitude
        A = self._cross_section.area.to('m ** 2').magnitude
        v = V / A
        return v * dh * rho / mu

    @property
    def velocity_pressure(self) -> Quantity:
        rho = self.fluid.rho.to('kg / m ** 3').magnitude
        V = self._volume_flow_rate.to('m ** 3 / s').magnitude
        A = self._cross_section.area.to('m ** 2').magnitude
        v = V / A
        pv = rho * v ** 2 / 2
        return pv * u.Pa

    @property
    def fittings_coefficient(self) -> float:
        return sum(fitting.zeta for fitting in self._fittings.values())

    def add_fitting(self, zeta: float, ID: str = ''):
        self._fittings[ID] = _Fitting(ID, zeta)
        self._calculate_pressure_drop()

    def remove_fitting(self, ID: str) -> Optional[_Fitting]:
        try:
            fitting = self._fittings.pop(ID)
        except KeyError:
            return None
        else:
            self._calculate_pressure_drop()
            return fitting

    def get_fitting_table(self, pressure_unit: str) -> Optional[pd.DataFrame]:
        if len(self._fittings):
            header = [
                'conduit ID',
                'fitting ID',
                'zeta',
                f'pressure drop [{pressure_unit}]'
            ]
            table = {
                header[0]: [],
                header[1]: [],
                header[2]: [],
                header[3]: []
            }
            for fitting in self._fittings.values():
                table[header[0]].append(self.ID)
                table[header[1]].append(fitting.ID)
                table[header[2]].append(fitting.zeta)
                table[header[3]].append(fitting.zeta * self.velocity_pressure.to(pressure_unit).magnitude)
            return pd.DataFrame(table)

    @property
    def numerator(self) -> Quantity:
        dp = self.pressure_drop.to('Pa')
        n = self.flow_sign * dp
        return n

    @property
    def denominator(self) -> Quantity:
        dp = self._pressure_drop
        V = self._volume_flow_rate
        d = 2.0 * dp / V
        if self.machine_coefficients is not None:
            m = self.machine_coefficients
            der_dp_machine = sum(m[i] * i * V ** (i - 1) for i in range(1, len(m)))
            d -= der_dp_machine
        return d


class Pipe(Conduit):
    pass


class Duct(Conduit):
    pass


TConduit = Union[PseudoConduit, Conduit, Pipe, Duct]


class Loop(List[AbstractConduit]):

    def __init__(self, ID: str = '', *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.ID = ID
        self.correction_term: Quantity = Quantity(float('nan'), 'm ** 3 / s')

    def calculate_correction_term(self):
        n = 0.0
        d = 0.0
        for conduit in self:
            n += conduit.numerator
            d += conduit.denominator
        self.correction_term = n / d

    @property
    def pressure_drop(self) -> Quantity:
        return sum(conduit.numerator for conduit in self)
