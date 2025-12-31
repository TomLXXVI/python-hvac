import warnings
from collections.abc import Sequence
from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import StrEnum
from itertools import combinations
import numpy as np
from scipy.optimize import fsolve
import CoolProp
import CoolProp.CoolProp as CP
from .. import Quantity
from .exceptions import CoolPropError, CoolPropMixtureError, CoolPropWarning

Q_ = Quantity


class CoolPropBackend(StrEnum):
    HEOS = 'HEOS'
    TTSE = 'TTSE&HEOS'  # tabular interpolation
    BICUBIC = 'BICUBIC&HEOS'  # tabular interpolation
    SRK = 'SRK'  # cubic equation of state Soave-Redlich-Kwong
    PR = 'PR'  # cubic equation of state Peng-Robinson
    INCOMP = 'INCOMP'  # incompressible pure fluids and mixtures


class ReferenceState(StrEnum):
    IIR = 'IIR'
    ASHRAE = 'ASHRAE'
    NBP = 'NBP'
    DEF = 'DEF'


@dataclass
class _CoolPropParameter:
    ID: int
    unit: str


CPP_ = _CoolPropParameter


# State variables defined in CoolProp
CPStateParameters: dict[str, _CoolPropParameter] = {
    'P': CPP_(CoolProp.iP, 'Pa'),
    'T': CPP_(CoolProp.iT, 'K'),
    'rho': CPP_(CoolProp.iDmass, 'kg / m**3'),
    'x': CPP_(CoolProp.iQ, 'frac'),
    'u': CPP_(CoolProp.iUmass, 'J / kg'),
    'h': CPP_(CoolProp.iHmass, 'J / kg'),
    's': CPP_(CoolProp.iSmass, 'J / kg / K'),
    'cp': CPP_(CoolProp.iCpmass, 'J / kg / K'),
    'cv': CPP_(CoolProp.iCvmass, 'J / kg / K'),
    'c': CPP_(CoolProp.ispeed_sound, 'm / s'),
    'k': CPP_(CoolProp.iconductivity, 'W / m / K'),
    'mu': CPP_(CoolProp.iviscosity, 'Pa * s'),
    'phase': CPP_(CoolProp.iPhase, '')
}


CPFluidParameters: dict[str, _CoolPropParameter] = {
    'P_crit': CPP_(CoolProp.iP_critical, 'Pa'),
    'T_crit': CPP_(CoolProp.iT_critical, 'K'),
    'P_trip': CPP_(CoolProp.iP_triple, 'Pa'),
    'T_trip': CPP_(CoolProp.iT_triple, 'K'),
    'P_min': CPP_(CoolProp.iP_min, 'Pa'),
    'P_max': CPP_(CoolProp.iP_max, 'Pa'),
    'T_min': CPP_(CoolProp.iT_min, 'K'),
    'T_max': CPP_(CoolProp.iT_max, 'K'),
    'molar_mass': CPP_(CoolProp.imolar_mass, 'kg / mol'),
    # incompressible solutions only:
    'T_freeze': CPP_(CoolProp.iT_freeze, 'K'),
    'fraction_min': CPP_(CoolProp.ifraction_min, 'frac'),
    'fraction_max': CPP_(CoolProp.ifraction_max, 'frac')
}


def _get_phase_description(phase_index: Quantity) -> str:
    if phase_index is not None:
        phase_index = phase_index.m
    else:
        phase_index = CoolProp.iphase_unknown
    match phase_index:
        case CoolProp.iphase_liquid:
            return 'sub-critical liquid'
        case CoolProp.iphase_supercritical:
            return 'super-critical (P > P_cnd, T > T_cnd)'
        case CoolProp.iphase_supercritical_gas:
            return 'super-critical gas (P < P_cnd, T > T_cnd)'
        case CoolProp.iphase_supercritical_liquid:
            return 'super-critical liquid (P > P_cnd, T < T_cnd)'
        case CoolProp.iphase_critical_point:
            return 'at the critical point'
        case CoolProp.iphase_gas:
            return 'sub-critical gas'
        case CoolProp.iphase_twophase:
            return 'two phase'
        case CoolProp.iphase_unknown:
            return 'unknown'
        case CoolProp.iphase_not_imposed:
            return 'not imposed'
        case _:
            return 'unknown'


@dataclass
class FluidState(ABC):
    name: str
    backend: str
    reference_state: str
    _state_dict: dict[str, Quantity]
    
    @property
    @abstractmethod
    def fluid(self) -> 'Fluid':
        ...


@dataclass
class PureFluidState(FluidState):

    def __post_init__(self):
        for k, v in self._state_dict.items():
            if k == 'phase':
                v = _get_phase_description(v)
                self._state_dict[k] = v
            setattr(self, k, v)

    def __str__(self) -> str:
        s = f"{self.name}\n"
        for k, v in self._state_dict.items():
            if isinstance(v, Quantity):
                s += f"{k} = {v:~P}\n"
            else:
                s += f"{k} = {v}\n"
        return s

    @property
    def fluid(self) -> 'PureFluid':
        return PureFluid(
            self.name,
            self.backend,
            self.reference_state
        )


class FractionType(StrEnum):
    MOLE = 'mole'
    MASS = 'mass'
    VOLUME = 'volume'


@dataclass
class Constituent:
    name: str
    mole_fraction: Quantity | None = None
    mass_fraction: Quantity | None = None
    volume_fraction: Quantity | None = None

    def __post_init__(self):
        # Determine the type of fraction (either mole, mass, or volume):
        if isinstance(self.mole_fraction, Quantity):
            self.fraction_type = FractionType.MOLE
        elif isinstance(self.mass_fraction, Quantity):
            self.fraction_type = FractionType.MASS
        elif isinstance(self.volume_fraction, Quantity):
            self.fraction_type = FractionType.VOLUME
        else:
            self.fraction_type = None
        # This will refer to the `PureFluid` instance of the constituent,
        # created when instantiating the mixture:
        self.fluid: PureFluid | None = None


@dataclass
class InteractionParams(dict):
    betaT: float | None = None
    gammaT: float | None = None
    betaV: float | None = None
    gammaV: float | None = None

    def __post_init__(self):
        for k, v in self.__dict__.items():
            self[k] = v
            

@dataclass
class MixtureState(FluidState):
    constituents: Sequence[Constituent]
    interaction_params: tuple[InteractionParams] | str | None

    def __post_init__(self):
        frac_type = self.constituents[0].fraction_type
        if frac_type is FractionType.MOLE:
            labels = [f"{c.name}[{c.mole_fraction}]" for c in self.constituents]
        elif frac_type is FractionType.MASS:
            labels = [f"{c.name}[{c.mass_fraction}]" for c in self.constituents]
        else:
            labels = [f"{c.name}[{c.volume_fraction}]" for c in self.constituents]
        self.constituent_names = ' & '.join(labels) + " " + str(frac_type.value)
        for k, v in self._state_dict.items():
            if k == 'phase':
                v = _get_phase_description(v)
                self._state_dict[k] = v
            setattr(self, k, v)

    def __str__(self) -> str:
        s = f"{self.name} | {self.constituent_names}\n"
        for k, v in self._state_dict.items():
            if isinstance(v, Quantity):
                s += f"{k} = {v:~P}\n"
            else:
                s += f"{k} = {v}\n"
        return s
    
    @property
    def fluid(self) -> 'Mixture':
        return Mixture(
            self.name,
            self.constituents,
            self.backend,
            self.reference_state,
            self.interaction_params
        )


class Fluid(ABC):

    def __init__(
        self,
        name: str,
        backend: str = CoolPropBackend.HEOS.value,
        reference_state: str = ReferenceState.DEF.value
    ) -> None:
        self.name = name
        self.backend = backend
        self.ref_state = reference_state
        self.abstract_state: CP.AbstractState | None = None

    @abstractmethod
    def _create_abstract_state(self, *args, **kwargs) -> CoolProp.AbstractState:
        ...

    @staticmethod
    def _generate_update_pair(**input_params):
        input_keys = list(input_params.keys())
        input_qties = list(input_params.values())
        update_pair = []
        for i in range(2):
            key = input_keys[i]
            cp_param = CPStateParameters[key]
            input_qty = input_qties[i]
            update_pair.extend([
                cp_param.ID,
                input_qty.to(cp_param.unit).magnitude
            ])
        up = CP.generate_update_pair(*update_pair)
        return up

    def _create_quantity(self, key: str) -> Quantity:
        cp_param = CPStateParameters[key]
        try:
            val = self.abstract_state.keyed_output(cp_param.ID)
        except ValueError as err:
            warnings.warn(
                f"Property '{key}' could not be determined. "
                f"CoolProp says: {err}",
                category=CoolPropWarning
            )
            qty = Q_(float('nan'), cp_param.unit)
        else:
            qty = Q_(val, cp_param.unit)
        return qty

    # noinspection PyTypeChecker
    @staticmethod
    def _get_phase(phase: str | None = None) -> int | None:
        match phase:
            case 'liquid':
                return CoolProp.iphase_liquid
            case 'gas':
                return CoolProp.iphase_gas
            case 'two_phase':
                return CoolProp.iphase_twophase
            case 'supercritical_liquid':
                return CoolProp.iphase_supercritical_liquid
            case 'supercritical_gas':
                return CoolProp.iphase_supercritical_gas
            case 'supercritical':
                return CoolProp.iphase_supercritical
            case None:
                return CoolProp.iphase_not_imposed
        return None
    
    @abstractmethod
    def __call__(self, **input_params: Quantity) -> FluidState:
        ...
    
    @abstractmethod
    def __deepcopy__(self, memo):
        ...
    
    @property
    def critical_point(self) -> FluidState:
        T_crit = self.abstract_state.T_critical()
        P_crit = self.abstract_state.p_critical()
        cpp_T_crit = CPFluidParameters['T_crit']
        cpp_P_crit = CPFluidParameters['P_crit']
        state = self.__call__(
            T=Q_(T_crit, cpp_T_crit.unit),
            P=Q_(P_crit, cpp_P_crit.unit)
        )
        return state

    @property
    def triple_point(self) -> FluidState:
        cpp_T_trip = CPFluidParameters['T_trip']
        T_trip = self.abstract_state.keyed_output(cpp_T_trip.ID)
        state = self.__call__(
            T=Q_(T_trip, cpp_T_trip.unit),
            x=Q_(0, 'frac')
        )
        return state

    def get_pressure_limits(self) -> tuple[Quantity, Quantity]:
        cpp_P_min = CPFluidParameters['P_min']
        cpp_P_max = CPFluidParameters['P_max']
        P_min = self.abstract_state.keyed_output(cpp_P_min.ID)
        P_max = self.abstract_state.keyed_output(cpp_P_max.ID)
        return Q_(P_min, cpp_P_min.unit), Q_(P_max, cpp_P_max.unit)

    def get_temperature_limits(self) -> tuple[Quantity, Quantity]:
        cpp_T_min = CPFluidParameters['T_min']
        cpp_T_max = CPFluidParameters['T_max']
        T_min = self.abstract_state.keyed_output(cpp_T_min.ID)
        T_max = self.abstract_state.keyed_output(cpp_T_max.ID)
        return Q_(T_min, cpp_T_min.unit), Q_(T_max, cpp_T_max.unit)

    @property
    def molar_mass(self) -> Quantity:
        cpp_mol_mass = CPFluidParameters['molar_mass']
        mol_mass = self.abstract_state.keyed_output(cpp_mol_mass.ID)
        return Q_(mol_mass, cpp_mol_mass.unit)


class PureFluid(Fluid):

    def __init__(
        self,
        name: str,
        backend: str = CoolPropBackend.HEOS.value,
        reference_state: str = ReferenceState.DEF.value
    ) -> None:
        """Creates an instance of `PureFluid`.

        Parameters
        ----------
        name:
            Name of the fluid as known by CoolProp.
        backend: optional
            Backend used by CoolProp to compute the state properties.
        reference_state: optional
            Reference state of the fluid.
        """
        super().__init__(name, backend, reference_state)
        self.abstract_state = self._create_abstract_state()

    def _create_abstract_state(self) -> CoolProp.AbstractState:
        CP.set_reference_state(
            f"{self.backend}::{self.name}",
            self.ref_state
        )
        return CoolProp.AbstractState(self.backend, self.name)

    def __call__(self, **input_params: Quantity) -> PureFluidState:
        phase = self._get_phase(input_params.pop('phase', None))
        update_pair = self._generate_update_pair(**input_params)
        try:
            self.abstract_state.specify_phase(phase)
            self.abstract_state.update(*update_pair)
        except ValueError as err:
            raise CoolPropError(err) from None
        state = PureFluidState(
            name=self.name,
            backend=self.backend,
            reference_state=self.ref_state,
            _state_dict={
                k: self._create_quantity(k) 
                for k in CPStateParameters.keys()
            }
        )
        return state
    
    def __deepcopy__(self, memo):
        return PureFluid(
            name=self.name,
            backend=self.backend,
            reference_state=self.ref_state
        )


class Mixture(Fluid):

    def __init__(
        self,
        name: str,
        constituents: Sequence[Constituent],
        backend: str = CoolPropBackend.HEOS.value,
        reference_state: str = ReferenceState.DEF.value,
        interaction_params: tuple[InteractionParams] | str | None = None
    ) -> None:
        """Creates a `Mixture` instance.
        
        Parameters
        ----------
        name:
            Name given to the mixture.
        constituents:
            The pure constituents of the mixture (`Constituent` objects).
        backend:
            CoolProp backend for the equation of state to be used (see 
            `CoolPropBackend` enum).
        reference_state:
            Pre-defined reference state or "zero-state" of the mixture (see 
            `ReferenceState` enum).
        interaction_params:
            The interaction parameters between each pair of constituents (see
            CoolProp's online documentation about mixtures). Either a tuple of
            `InteractionParams` objects (one for each pair of constituents), or 
            either one of the strings 'linear' or 'Lorentz-Berthelot'.
        """
        super().__init__(name, backend, reference_state)
        self.constituents = self._create_pure_fluids(constituents)
        self.interaction_params = interaction_params
        self.abstract_state = self._create_abstract_state(interaction_params)
        self._set_fractions()

    def _create_pure_fluids(
        self,
        constituents: Sequence[Constituent]
    ) -> Sequence[Constituent]:
        for constituent in constituents:
            fluid = PureFluid(constituent.name, self.backend, self.ref_state)
            constituent.fluid = fluid
        return constituents

    def _create_abstract_state(
        self,
        interaction_params: tuple[InteractionParams] | str | None
    ) -> CoolProp.AbstractState:
        if interaction_params is not None:
            CAS_numbers = tuple(
                CP.get_fluid_param_string(c.name, 'CAS')
                for c in self.constituents
            )
            CAS_pairs = self._get_binary_combinations(CAS_numbers)
            if isinstance(interaction_params, str):
                for cp in CAS_pairs: 
                    CP.apply_simple_mixing_rule(cp[0], cp[1], interaction_params)
            # Dummy entry in the library of interaction parameters.
            if isinstance(interaction_params, tuple):
                for cp in CAS_pairs:
                    CP.apply_simple_mixing_rule(cp[0], cp[1], 'linear')
            # Set configuration flag to allow binary interaction parameters to 
            # be overwritten.
            CP.set_config_bool(CP.OVERWRITE_BINARY_INTERACTION, True)

        name = '&'.join([c.name for c in self.constituents])
        abstract_state = CoolProp.AbstractState(self.backend, name)

        if isinstance(interaction_params, tuple):
            constituent_indices = tuple(i for i in range(len(self.constituents)))
            index_pairs = self._get_binary_combinations(constituent_indices)
            for i, ip in enumerate(index_pairs):
                interaction_params_pair_i = interaction_params[i] 
                for k, v in interaction_params_pair_i.items():
                    if isinstance(v, float):
                        abstract_state.set_binary_interaction_double(ip[0], ip[1], k, v)

        return abstract_state
    
    @staticmethod
    def _get_binary_combinations(tuple_input):
        pairs = list(combinations(tuple_input, 2))
        return pairs
    
    def _set_fractions(self) -> None:
        fraction_types = [c.fraction_type for c in self.constituents]
        frac_ref_type = fraction_types[0]
        if all([frac_type == frac_ref_type for frac_type in fraction_types]):
            if frac_ref_type is FractionType.MOLE:
                fractions = [
                    c.mole_fraction.to('frac').m
                    for c in self.constituents
                ]
                self.abstract_state.set_mole_fractions(fractions)
            elif frac_ref_type is FractionType.MASS:
                fractions = [
                    c.mass_fraction.to('frac').m
                    for c in self.constituents
                ]
                self.abstract_state.set_mass_fractions(fractions)
            elif frac_ref_type is FractionType.VOLUME:
                fractions = [
                    c.volume_fraction.to('frac').m
                    for c in self.constituents
                ]
                self.abstract_state.set_volu_fractions(fractions)
        else:
            raise ValueError('all fractions must be of the same type')
   
    def _solve_iteratively(self, **input_params: Quantity) -> dict[str, Quantity]:
        keys = list(input_params.keys())
        val_key_1 = keys[0] if keys[0] in ('P', 'T', 'x') else keys[1]
        val_key_2 = keys[1] if keys[1] not in ('P', 'T', 'x') else keys[0]
        search_key = keys[2]
        search_unit = CPStateParameters[search_key].unit
        search_value_ini = input_params[search_key].to(search_unit).m
        qty_1 = input_params[val_key_1]
        val_2 = input_params[val_key_2].to_base_units().m
        
        def _objective(search_val_: np.ndarray) -> np.ndarray:
            search_qty_ = Q_(search_val_[0], search_unit)
            state_ = self(**{val_key_1: qty_1, search_key: search_qty_})
            # noinspection PyProtectedMember
            _temp_val = state_._state_dict[val_key_2].to_base_units().m
            deviation = val_2 - _temp_val
            return np.array([deviation])
        
        search_val = fsolve(_objective, np.array([search_value_ini]))[0]
        search_qty = Q_(search_val, search_unit)
        return {val_key_1: qty_1, search_key: search_qty}
   
    def __call__(self, **input_params: Quantity) -> MixtureState:
        """Determines the full mixture state when known state variables are
        given.
        
        Parameters
        ----------
        input_params:
            Keyword arguments `key=value, ...` where `key` is the symbol of the
            known state variable and `value` its corresponding value. Normally,
            two known state variables suffice to determine the full state. 
            However, in case of mixtures the allowable combinations of state 
            variables are limited to: 
            - pressure (`P`) and temperature (`T`), 
            - `P` and vapor quality (`x`),
            - `T` and `x`
            To bypass this limitation somehow, a third keyword argument can be
            passed to determine the mixture state iteratively. The first 
            (or second) keyword argument must be either `P`, `T`, or `x`. The 
            second (or first) keyword argument can be any known state variable 
            with its value. The third keyword argument must be again either 
            `P`, `T`, or `x` with an initially guessed value being assigned to 
            it for starting the iterative routine.
        
        Raises
        ------
        CoolPropMixtureError:
            If the mixture state could not be determined. 
        """
        phase = self._get_phase(input_params.pop('phase', None))
        update_pair = self._generate_update_pair(**input_params)
        try:
            self.abstract_state.specify_phase(phase)
            self.abstract_state.update(*update_pair)
        except ValueError:
            keys = list(input_params.keys())[:2]
            mask = [k in ('P', 'T', 'x') for k in keys]
            if all(mask) or not any(mask) or len(input_params.keys()) < 3:
                raise CoolPropMixtureError(
                    "The mixture state can only be solved directly for P & T, "
                    "P & x, or T & x."
                )
            else:
                input_params = self._solve_iteratively(**input_params)
                update_pair = self._generate_update_pair(**input_params)
                try:
                    self.abstract_state.update(*update_pair)
                except ValueError as err:
                    raise CoolPropMixtureError(err)
        state = MixtureState(
            name=self.name,
            backend=self.backend,
            reference_state=self.ref_state,
            constituents=self.constituents,
            interaction_params=self.interaction_params,
            _state_dict={k: self._create_quantity(k) for k in CPStateParameters.keys()}
        )
        return state

    def get_fraction_limits(self) -> tuple:
        cpp_frac_min = CPFluidParameters['fraction_min']
        cpp_frac_max = CPFluidParameters['fraction_max']
        frac_min = self.abstract_state.keyed_output(cpp_frac_min.ID)
        frac_max = self.abstract_state.keyed_output(cpp_frac_max.ID)
        return Q_(frac_min, cpp_frac_min.unit), Q_(frac_max, cpp_frac_max.unit)

    def get_freezing_temp(self) -> Quantity:
        cpp_T_freeze = CPFluidParameters['T_freeze']
        T_freeze = self.abstract_state.keyed_output(cpp_T_freeze.ID)
        return Q_(T_freeze, cpp_T_freeze.unit)
    
    def __deepcopy__(self, memo):
        return Mixture(
            name=self.name,
            constituents=self.constituents,
            backend=self.backend,
            reference_state=self.ref_state,
            interaction_params=self.interaction_params
        )
