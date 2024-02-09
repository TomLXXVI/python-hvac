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
from .exceptions import ConduitConfigurationError, ConduitIterationOverflow

u = UNITS
Q_ = Quantity
Water = Fluid('Water')


@dataclass
class Fitting:
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
    """Class that represents a network node. Any conduit in the network needs
    to be assigned a start and an end node. Nodes are created internally when
    creating a network and adding conduits to a network.
    """

    def __init__(self, ID: str = '', height: Quantity = 0 * u.m) -> None:
        """Creates a `Node` instance.

        Parameters
        ----------
        ID:
            Identifier to identify the specific node in the network.
        height:
            The height of the node above a given reference plane to which the
            heights of all nodes in the network are referenced.
        """
        self.ID = ID
        self.height = height
        self.incoming: List['Conduit'] = []
        self.outgoing: List['Conduit'] = []
        self.conduit_IDs: List[str] = []
        self._flow_sign_in_init: Dict[str, Quantity] = {}
        self._flow_sign_out_init: Dict[str, Quantity] = {}

    def connect(self, conduit: 'Conduit', node_arrow: NodeArrow) -> None:
        """Connect a conduit to the node.

        Parameters
        ----------
        conduit: Conduit, {Pipe, Duct}
            The conduit to connect to the node.
        node_arrow: NodeArrow, {NodeArrow.INCOMING, NodeArrow.OUTGOING}
            Indicates if the conduit arrives at the node or leaves from the
            node.
        """
        if conduit.ID not in self.conduit_IDs:
            self.conduit_IDs.append(conduit.ID)
            if node_arrow == NodeArrow.INCOMING:
                self.incoming.append(conduit)
                self._flow_sign_in_init[conduit.ID] = conduit.flow_sign
            if node_arrow == NodeArrow.OUTGOING:
                self.outgoing.append(conduit)
                self._flow_sign_out_init[conduit.ID] = conduit.flow_sign

    @property
    def volume_flow_rate(self) -> Quantity:
        """Returns the net volume flow rate into the node. A negative sign means
        that the net volume flow rate is leaving the node.

        The net volume flow rate of a node should be (very close) to zero if
        this node has no external connections, else the net volume flow rate
        should be equal to the known external flow rate that enters
        (positive sign) or leaves (negative sign) the node (see physical law of
        continuity).
        """
        V_in = [Q_(0, 'm ** 3  / s')]
        V_out = [Q_(0, 'm ** 3 / s')]
        for conduit in self.incoming:
            if not isinstance(conduit, PseudoConduit):
                flow_sign_in_init = self._flow_sign_in_init[conduit.ID]
                if conduit.flow_sign != flow_sign_in_init:
                    V_out.append(conduit.volume_flow_rate)
                else:
                    V_in.append(conduit.volume_flow_rate)
        for conduit in self.outgoing:
            if not isinstance(conduit, PseudoConduit):
                flow_sign_out_init = self._flow_sign_out_init[conduit.ID]
                if conduit.flow_sign != flow_sign_out_init:
                    V_in.append(conduit.volume_flow_rate)
                else:
                    V_out.append(conduit.volume_flow_rate)
        return sum(V_out) - sum(V_in)


class AbstractConduit(ABC):
    """Parent class from which all types of conduit are derived. It contains
    the attributes that all types of conduit have in common and that are
    needed to design and analyze a network of conduits. This class also declares
    a number of abstract methods that must be implemented by its derived classes.

    Attributes
    ----------
    ID: str
        Identifier to identify a conduit in a network.
    loop_ID: str, tuple[str, str]
        Identifier to identify the loop or loops to which a conduit belongs in
        a network. A conduit can belong to 2 loops maximum.
    loops: list[Loop]
        The `Loop` object(s) to which a conduit belongs in a network.
    start_node: Node
        The start node is the node where the conduit leaves.
    end_node: Node
        The end node is the node where the conduit arrives.
    """
    def __init__(self) -> None:
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
    """A pseudo conduit is a fictitious conduit used for network analysis.
    It represents a fixed pressure difference between two nodes of a network.
    """

    def __init__(self):
        super().__init__()
        self._fixed_pressure_drop: Quantity = float('nan') * u.Pa

    @classmethod
    def create(cls, fixed_pressure_drop: Quantity) -> 'PseudoConduit':
        """Creates a `PseudoConduit` with a fixed pressure difference."""
        obj = cls()
        obj._fixed_pressure_drop = fixed_pressure_drop
        return obj

    @staticmethod
    def duplicate(conduit: 'PseudoConduit', **kwargs) -> 'PseudoConduit':
        """Returns a duplicate (deep copy) of the conduit. Through **kwargs
        attributes of the copied conduit can be modified.

        kwargs
        ------
        ID: str
        start_node: Node
        end_node: Node
        flow_sign: FlowSign
        loops: list[Loop]
        fixed_pressure_drop: Quantity
        """
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
        """Used internally when analyzing a network with the Hardy-Cross method."""
        return self._fixed_pressure_drop.to('Pa')

    @property
    def denominator(self) -> Quantity:
        """Used internally when analyzing a network with the Hardy-Cross method."""
        return Q_(0.0, 'Pa * s / m ** 3')

    @property
    def pressure_drop(self) -> Quantity:
        """Returns the fixed pressure difference between the start and the end
        node of the pseudo conduit.
        """
        return self._fixed_pressure_drop


class Conduit(AbstractConduit):
    MAX_ITERATIONS = 30
    SOLUTION_TOLERANCE = 1e-5

    def __init__(self):
        super().__init__()
        self.length: Quantity = float('nan') * u.m
        self.wall_roughness: Quantity = float('nan') * u.mm
        self.fluid: FluidState = Water(T=Q_(15.0, 'degC'), P=Q_(101.325, 'kPa'))
        self._cross_section: Optional[TCrossSection] = None
        self._volume_flow_rate: Quantity = Quantity(float('nan'), 'm ** 3 / s')
        self._pressure_drop: Quantity = float('nan') * u.Pa
        self._fittings: Dict[str, Fitting] = {}
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
        """Creates an instance of one of the classes derived from class
        `Conduit` (see class `Pipe` and class `Duct`).

        Parameters
        ----------
        length:
            Length of the conduit.
        wall_roughness:
            Absolute wall roughness of the conduit.
        fluid:
            (Incompressible) fluid that flows in the conduit, used internally to
            get fluid properties, like mass density and viscosity of the fluid.
        cross_section:
            Cross-section of the conduit (see module cross_section.py).
        volume_flow_rate: optional
            Volume flow rate of the fluid.
        pressure_drop: optional
            Pressure drop across the conduit.
        specific_pressure_drop: optional
            Pressure drop per unit length of the conduit (this is used, e.g., in
            the equal friction method to size air ducts).
        machine_coefficients: optional
            List of polynomial coefficients that describe the pump or fan curve
            if a pump or fan should be present in the conduit.

        Notes
        -----
        Depending on the given input parameters, the appropriate output is
        calculated on instantiation of the `Conduit` object and accessible after
        instantiation.
        1.  If flow rate and (specific) pressure drop are given, the output
            will be the conduit's cross-section size.
        2.  If (specific) pressure drop and the conduit's cross-section size are
            given, the output will be the volume flow rate.
        3.  If volume flow rate and the conduit's cross-section size are given,
            the output will be the pressure drop across the conduit.

        Returns
        -------
        Conduit
        """
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
        if isinstance(specific_pressure_drop, Quantity):
            obj._pressure_drop = specific_pressure_drop * length
        else:
            obj._pressure_drop = pressure_drop
        obj.machine_coefficients = machine_coefficients
        # depending on the available input parameters, calculate the appropriate
        # output
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
            raise ConduitConfigurationError('arguments missing')
        return obj

    @staticmethod
    def duplicate(conduit: 'Conduit', **kwargs) -> 'Conduit':
        """Returns a duplicate (deep copy) of the conduit. Through **kwargs
        attributes of the copied conduit can be modified.

        kwargs
        ------
        length: Quantity
        wall_roughness: Quantity
        fluid: FluidState
        cross-section: {Circular, Rectangular, FlatOval}
        volume_flow_rate: Quantity
        pressure_drop: Quantity
        specific_pressure_drop: Quantity
        fittings:
        ID: str
        start_node: Node
        end_node: Node
        flow_sign: FlowSign
        loops: list[Loop]
        machine_coefficients: list[Quantity]
        """
        obj = deepcopy(conduit)
        obj.length = kwargs.get('length', obj.length)
        obj.wall_roughness = kwargs.get('wall_roughness', obj.wall_roughness)
        obj.fluid = kwargs.get('fluid', obj.fluid)
        obj._cross_section = kwargs.get('cross_section', obj._cross_section)
        obj._volume_flow_rate = kwargs.get('volume_flow_rate', obj._volume_flow_rate)
        obj._pressure_drop = kwargs.get('pressure_drop', obj._pressure_drop)
        if (specific_pressure_drop := kwargs.get('specific_pressure_drop')) is not None:
            obj._pressure_drop = specific_pressure_drop * obj.length
        obj._fittings = kwargs.get('fittings', obj._fittings)
        obj.ID = kwargs.get('ID', obj.ID)
        obj.start_node = kwargs.get('start_node', obj.start_node)
        obj.end_node = kwargs.get('end_node', obj.end_node)
        obj.flow_sign = kwargs.get('flow_sign', obj.flow_sign)
        obj.loops = kwargs.get('loops', obj.loops)
        obj.machine_coefficients = kwargs.get('machine_coefficients', obj.machine_coefficients)
        return obj

    # @staticmethod
    # def duplicate(conduit: 'Conduit', **kwargs) -> 'Conduit':
    #     obj = Conduit.create(
    #         length=deepcopy(conduit.length),
    #         wall_roughness=deepcopy(conduit.wall_roughness),
    #         fluid=conduit.fluid.fluid(T=conduit.fluid.T, P=conduit.fluid.P),
    #         cross_section=deepcopy(conduit.cross_section),
    #         volume_flow_rate=deepcopy(conduit.volume_flow_rate),
    #         pressure_drop=deepcopy(conduit._pressure_drop),
    #         specific_pressure_drop=None,
    #         machine_coefficients=deepcopy(conduit.machine_coefficients)
    #     )
    #     obj.length = kwargs.get('length', obj.length)
    #     obj.wall_roughness = kwargs.get('wall_roughness', obj.wall_roughness)
    #     obj.fluid = kwargs.get('fluid', obj.fluid)
    #     obj._cross_section = kwargs.get('cross_section', obj._cross_section)
    #     obj._volume_flow_rate = kwargs.get('volume_flow_rate', obj._volume_flow_rate)
    #     obj._pressure_drop = kwargs.get('pressure_drop', obj._pressure_drop)
    #     if (specific_pressure_drop := kwargs.get('specific_pressure_drop')) is not None:
    #         obj._pressure_drop = specific_pressure_drop * obj.length
    #     obj._fittings = kwargs.get('fittings', deepcopy(conduit._fittings))
    #     obj.ID = kwargs.get('ID', deepcopy(conduit.ID))
    #     obj.start_node = kwargs.get('start_node', deepcopy(conduit.start_node))
    #     obj.end_node = kwargs.get('end_node', deepcopy(conduit.end_node))
    #     obj.flow_sign = kwargs.get('flow_sign', deepcopy(conduit.flow_sign))
    #     obj.loops = kwargs.get('loops', deepcopy(conduit.loops))
    #     obj.machine_coefficients = kwargs.get('machine_coefficients', obj.machine_coefficients)
    #     return obj

    def _calculate_darcy_friction_factor(self, Re: float, e: float) -> float:
        """Returns the Darcy friction factor that corresponds with the given
        Reynolds number `Re` and wall roughness `e`.
        """
        f = self._serghide(Re, e)
        return f

    @staticmethod
    def _churchill(Re: float, e: float) -> float:
        """Darcy friction factor acc. to correlation of Churchill."""
        d = (7 / Re) ** 0.9 + 0.27 * e
        A = (2.457 * math.log(1 / d)) ** 16
        B = (37.530 / Re) ** 16
        f = 8 * ((8 / Re) ** 12 + 1 / (A + B) ** 1.5) ** (1 / 12)
        return f

    @staticmethod
    def _serghide(Re: float, e: float) -> float:
        """Darcy friction factor acc. to correlation of Serghide."""
        if Re > 0.0:
            f1 = -2.0 * math.log10(e / 3.7 + 12.0 / Re)
            try:
                f2 = -2.0 * math.log10(e / 3.7 + 2.51 * f1 / Re)
            except ValueError as err:
                print('Oops!')
                raise err
            f3 = -2.0 * math.log10(e / 3.7 + 2.51 * f2 / Re)
            f = (f1 - (f2 - f1) ** 2.0 / (f3 - 2.0 * f2 + f1)) ** -2.0
            return f
        else:
            # no flow means no friction
            return 0.0

    def _calculate_equivalent_diameter(self):
        """Calculates the equivalent diameter of the cross-section based on the
        volume flow rate and pressure drop.
        """
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
            raise ConduitIterationOverflow('no solution for diameter')
        self._cross_section.equivalent_diameter = d_eq * u.m  # this assignment may also change the cross-section's size
        if self._cross_section.schedule:
            # if also schedule is set, a commercially available size will be
            # determined automatically, which makes it necessary to recalculate
            # the pressure drop based on the real equivalent diameter of the
            # commercially available size.
            self._calculate_pressure_drop()

    def _calculate_volume_flow_rate(self):
        """Calculates the volume flow rate through the conduit based on the
        pressure drop across the conduit and the hydraulic diameter of the
        cross-section.
        """
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
            raise ConduitIterationOverflow('no solution for volume flow rate')
        A = self._cross_section.area.to('m ** 2').magnitude
        self._volume_flow_rate = Quantity(v * A, 'm ** 3 / s')

    def _calculate_pressure_drop(self):
        """Calculates the pressure drop across the conduit due to friction and
        fittings in the conduit based on the volume flow rate through the
        conduit and the hydraulic diameter of the cross-section.
        The total pressure loss is calculated, without considering possible
        pressure gain due to the presence of a pump or fan in the conduit.
        """
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
        """Get the volume flow rate through the conduit."""
        self._calculate_volume_flow_rate()
        return self._volume_flow_rate

    @volume_flow_rate.setter
    def volume_flow_rate(self, value: Quantity):
        """Set the volume flow rate through the conduit."""
        if value < 0.0:
            self._volume_flow_rate = abs(value)
            self.flow_sign = self.flow_sign.reverse()
            # `flow_sign` is only used internally when analyzing a network with
            # the Hardy-Cross method.
        else:
            self._volume_flow_rate = value
        # when the flow rate through a conduit changes, the pressure drop across
        # the conduit will also change
        self._calculate_pressure_drop()

    @property
    def pressure_drop(self) -> Quantity:
        """Get the net pressure drop (pressure loss) across the conduit, also
        taking possible pressure gain due to a pump or fan into account.
        """
        dp = self._pressure_drop
        # Internally, `_pressure_drop` is the total pressure drop without
        # considering possible pressure gain due to the presence of a pump
        # or fan in the conduit. In that case, the real pressure drop across
        # the conduit is returned after subtracting the pressure gain.
        if self.machine_coefficients is not None:
            V = self.volume_flow_rate
            m = self.machine_coefficients
            dp_machine = sum(m[i] * V ** i for i in range(len(m)))
            dp -= dp_machine
        return dp

    @pressure_drop.setter
    def pressure_drop(self, value: Quantity):
        """Set the pressure loss across the conduit (not taking any pressure
        gain due to pump or fan into account).
        """
        self._pressure_drop = value
        # when the pressure across a conduit changes, the flow rate through the
        # conduit will also change
        self._calculate_volume_flow_rate()

    @property
    def specific_pressure_drop(self) -> Quantity:
        """Get the specific pressure drop along the conduit (not taking any
        pressure gain due to pump or fan into account).
        """
        return self._pressure_drop / self.length

    @specific_pressure_drop.setter
    def specific_pressure_drop(self, value: Quantity):
        """Set the specific pressure drop across the conduit (not taking any
        pressure gain due to pump or fan into account).
        """
        self.pressure_drop = value * self.length

    @property
    def cross_section(self) -> TCrossSection:
        """Get the cross-section of the conduit."""
        return self._cross_section

    @cross_section.setter
    def cross_section(self, value: TCrossSection):
        """Set the cross-section of the conduit."""
        self._cross_section = value
        # when the cross-section of a conduit changes, the pressure drop across
        # the conduit will also change for a given network flow rate.
        self._calculate_pressure_drop()

    @property
    def velocity(self) -> Quantity:
        """Get the mean flow velocity in the conduit."""
        return self._volume_flow_rate / self._cross_section.area

    @property
    def reynolds_number(self) -> float:
        """Get the Reynolds number."""
        rho = self.fluid.rho.to('kg / m ** 3').magnitude
        mu = self.fluid.mu.to('Pa * s').magnitude
        dh = self._cross_section.hydraulic_diameter.to('m').magnitude
        V = self._volume_flow_rate.to('m ** 3 / s').magnitude
        A = self._cross_section.area.to('m ** 2').magnitude
        v = V / A
        return v * dh * rho / mu

    @property
    def velocity_pressure(self) -> Quantity:
        """Get the velocity pressure in the conduit."""
        rho = self.fluid.rho.to('kg / m ** 3').magnitude
        V = self._volume_flow_rate.to('m ** 3 / s').magnitude
        A = self._cross_section.area.to('m ** 2').magnitude
        v = V / A
        pv = rho * v ** 2 / 2
        return pv * u.Pa

    @property
    def fittings_coefficient(self) -> float:
        """Get the total resistance coefficient of any fittings in the conduit."""
        return sum(fitting.zeta for fitting in self._fittings.values())

    def add_fitting(self, zeta: float, ID: str = ''):
        """Add a fitting to the conduit.

        Parameters
        ----------
        zeta:
            Resistance coefficient of the fitting.
        ID:
            Identifier of the fitting, making it possible to modify the
            resistance coefficient of this fitting or remove the fitting
            afterward.
        """
        self._fittings[ID] = Fitting(ID, zeta)
        self._calculate_pressure_drop()

    def remove_fitting(self, ID: str) -> Optional[Fitting]:
        """Removes the fitting with the given ID from the conduit."""
        try:
            fitting = self._fittings.pop(ID)
        except KeyError:
            return None
        else:
            self._calculate_pressure_drop()
            return fitting

    def get_fitting_table(self, pressure_unit: str) -> Optional[pd.DataFrame]:
        """Returns a Pandas DataFrame object with an overview of all the fittings
        in the conduit and having 4 columns: 'conduit ID', 'fitting ID', 'zeta'
        and 'pressure drop [<pressure_unit>]' with the pressure drop across the
        fittings expressed in the given pressure unit.
        """
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
                zeta = float(fitting.zeta)
                table[header[0]].append(self.ID)
                table[header[1]].append(fitting.ID)
                table[header[2]].append(zeta)
                table[header[3]].append(zeta * self.velocity_pressure.to(pressure_unit).magnitude)
            return pd.DataFrame(table)

    @property
    def numerator(self) -> Quantity:
        """Used internally when analyzing a network with the Hardy-Cross method."""
        dp = self.pressure_drop.to('Pa')
        n = self.flow_sign * dp
        return n

    @property
    def denominator(self) -> Quantity:
        """Used internally when analyzing a network with the Hardy-Cross method."""
        V = self._volume_flow_rate
        if V.magnitude != 0.0:
            dp = self._pressure_drop
            d = 2.0 * dp / V
            if self.machine_coefficients is not None:
                m = self.machine_coefficients
                der_dp_machine = sum(m[i] * i * V ** (i - 1) for i in range(1, len(m)))
                d -= der_dp_machine
            return d
        else:
            return Quantity(0.0, 'Pa / (m ** 3 / s)')


class Pipe(Conduit):
    """Conduit that represents a pipe in a pipe network."""
    pass


class Duct(Conduit):
    """Conduit that represents a duct in a duct network."""
    pass


# Define a general `TConduit` annotation type for type annotations in the source
# code.
TConduit = Union[PseudoConduit, Conduit, Pipe, Duct]


class Loop(List[AbstractConduit]):
    """Class derived from `list` that represents a `Loop` in a network of conduits.
    Used for analyzing a network with the Hardy-Cross method.
    """
    def __init__(self, ID: str = '', *args, **kwargs):
        """Creates a `Loop` instance.

        Parameters
        ----------
        ID:
            Identifier to identify the specific loop in the network.
        """
        super().__init__(*args, **kwargs)
        self.ID = ID
        self.correction_term: Quantity = Quantity(float('nan'), 'm ** 3 / s')

    def calculate_correction_term(self):
        """Used internally when analyzing a network with the Hardy-Cross method."""
        n = 0.0
        d = 0.0
        for conduit in self:
            n += conduit.numerator
            d += conduit.denominator
        try:
            self.correction_term = n / d
        except ZeroDivisionError as err:
            print('Oops!')
            raise err

    @property
    def pressure_drop(self) -> Quantity:
        """Returns the pressure drop along the loop. After analysis with the
        Hardy-Cross method, the loop pressure drop should be nearly equal to
        zero (as any point in the loop can have only one pressure, which is
        also a consequence of the physical law of conservation of energy).
        """
        return sum(conduit.numerator for conduit in self)
