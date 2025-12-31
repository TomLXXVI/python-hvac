import warnings
from typing import List, Optional, Dict, Tuple
from dataclasses import dataclass
import CoolProp
import CoolProp.CoolProp as CP
import numpy as np
from scipy.optimize import fsolve
from .exceptions import CoolPropWarning, CoolPropError, CoolPropMixtureError
from .. import Quantity


Q_ = Quantity


@dataclass
class FluidState:
    fluid_attrs: Dict[str, Quantity]
    state_dict: Dict[str, Quantity]

    def __post_init__(self):
        for k, v in self.state_dict.items():
            if k == 'phase':
                v = self._get_phase_description(v)
            setattr(self, k, v)

    @staticmethod
    def _get_phase_description(phase_index: Quantity) -> str | None:
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

        return None

    @property
    def fluid(self) -> 'Fluid':
        # workaround for the fact that a `Fluid` object cannot be pickled
        # (not even with `dill`) due to CoolProp's `AbstractState` object.
        return Fluid(**self.fluid_attrs)


class Fluid:
    _coolprop_qties: Dict[str, Tuple[int, str]] = {
        'T': (CoolProp.iT, 'K'),
        'P': (CoolProp.iP, 'Pa'),
        'rho': (CoolProp.iDmass, 'kg / m **3'),
        'h': (CoolProp.iHmass, 'J / kg'),
        's': (CoolProp.iSmass, 'J / kg / K'),
        'cp': (CoolProp.iCpmass, 'J / kg / K'),
        'cv': (CoolProp.iCvmass, 'J / kg / K'),
        'x': (CoolProp.iQ, 'frac'),
        'k': (CoolProp.iconductivity, 'W / m / K'),
        'mu': (CoolProp.iviscosity, 'Pa * s'),  # dynamic or absolute viscosity
        'beta': (CoolProp.iisobaric_expansion_coefficient, '1 / K'),
        'sigma': (CoolProp.isurface_tension, 'N / m'),
        'phase': (CoolProp.iPhase, ''),
        'speed_of_sound': (CoolProp.ispeed_sound, 'm / s'),
    }

    def __init__(
        self,
        name: str,
        backend: str = 'HEOS',
        mass_fractions: Optional[List[Quantity]] = None,
        vol_fractions: Optional[List[Quantity]] = None,
        reference: str = 'DEF'
    ):
        """
        Creates a `Fluid`-instance.

        Parameters
        ----------
        name: str
            Name of the fluid or mixture. CoolProp must know the fluid or the
            constituents of the mixture. In case of a mixture, the constituents
            are separated by an ampersand (&) without spaces.
        backend: str, default: 'HEOS'
            The backend CoolProp must use to perform state calculations. See
            CoolProp's documentation for which backends are possible.
        mass_fractions: List[Quantity], default `None`
            Only for mixtures or incompressible fluids which are mass-based
            binary mixtures (i.e., water-based mixtures). The mass fractions of
            the constituents in the mixture in the same order as the
            constituents given in `name`. In case of a binary, water-based
            mixture only the mass fraction of the constituent which is not water
            needs to be specified.
        vol_fractions: List[Quantity], default `None`
            Only for mixtures or incompressible fluids which are volume-based
            binary mixtures (i.e., water-based mixtures). The volume fractions
            of the constituents in the same order as the constituents given in
            `name`. In case of a binary, water-based mixture only the volume
            fraction of the constituent which is not water needs to be specified.
        reference: str, default 'DEF'
            Determines the reference state for enthalpy and entropy. See
            CoolProp's documentation for the available possibilities.
        """
        self.fluid_name = name
        self.backend = backend
        self.mass_fractions = [mf.to('frac').m for mf in mass_fractions] if mass_fractions else None
        self.vol_fractions = [vf.to('frac').m for vf in vol_fractions] if vol_fractions else None
        self.reference = reference
        self._constituents: List[str] = []
        self._create_state_object()

    def _create_state_object(self) -> None:
        """Create CoolProp's AbstractState-object (we will call it the state
        object) with the given backend for the given fluid.
        """
        # Set the reference state for enthalpy and entropy.
        if self.reference != 'DEF':
            if not self._is_mixture():
                CP.set_reference_state(self.fluid_name, self.reference)
            else:
                for constituent in self._constituents:
                    CP.set_reference_state(constituent, self.reference)

        # Create the state object.
        self._state = CoolProp.AbstractState(self.backend, self.fluid_name)

        # In case of mixture: set the mass fractions of the constituents
        if self.mass_fractions is not None:
            self._state.set_mass_fractions(self.mass_fractions)

        # In case of mixture: set the volume fractions of the constituents
        if self.vol_fractions is not None:
            self._state.set_volu_fractions(self.vol_fractions)

    def _is_mixture(self) -> bool:
        """Check whether the fluid is a mixture."""
        self._constituents = self.fluid_name.split('&')
        if len(self._constituents) > 1:
            return True
        return False

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

    def _update(self, phase: str | None = None, **input_qties: Quantity) -> None:
        """Updates the state object of the `Fluid`-instance based on the
        given state variables in `input_qties`. Parameter `phase` can be used
        to impose the phase ('liquid', 'gas', 'two_phase', 'supercritical_liquid',
        'supercritical_gas', or 'supercritical').
        """
        phase = self._get_phase(phase)
        if len(input_qties) == 2:  # normal case
            input_qty_names = list(input_qties.keys())
            input_qty_values = list(input_qties.values())

            # First input state variable:
            # Get the CoolProp-ID of this variable and the unit that CoolProp uses
            # for this variable.
            coolprop_qty1, coolprop_unit1 = self._coolprop_qties[input_qty_names[0]]
            # Convert the value of the first state variable to the CoolProp-unit.
            val1 = input_qty_values[0].to(coolprop_unit1).m

            # Do the same with the second input state variable.
            coolprop_qty2, coolprop_unit2 = self._coolprop_qties[input_qty_names[1]]
            val2 = input_qty_values[1].to(coolprop_unit2).m

            # Use CoolProp's functions to update the state object.
            inputs = CoolProp.CoolProp.generate_update_pair(
                coolprop_qty1, val1,
                coolprop_qty2, val2
            )
            try:
                self._state.specify_phase(phase)
                self._state.update(*inputs)
            except ValueError as err:
                if self._is_mixture():
                    raise CoolPropMixtureError(err) from None  # signal to caller that mixture state cannot be solved
                else:
                    raise CoolPropError(err) from None

        if len(input_qties) == 3:  # special case for mixture
            self._find_mixture_state(**input_qties)

    def _find_mixture_state(self, **qties) -> None:
        """Workaround for mixtures when method `update` of `self._state`
        (CoolProp `AbstractState`-object) throws a `ValueError` exception because
        it cannot solve for the given combination of input state variables (the
        only combinations it can solve for are P and T or P and x or T and x).

        The workaround is based on a root finding algorithm. This implies that
        a third input state variable must be given of which the given value is
        just an initial guess for the root finding algorithm. The first input
        state variable must be either P, T or x (this is a restriction, but there
        is no other way). The second input state variable can be any valid state
        variable (see the keys of `_coolprop_qties`). The third input state
        variable must be again either P, T or x (but not the same one used for
        the first input state variable).

        The goal is to find the correct value of the third state variable so
        that at the state defined by the value of this third state variable and
        the value of the first one (which will be a combination of P and T or
        P and x or T and x) the value of the second state variable is equal to
        the value that we passed in for this second state variable.
        """
        qties_names = list(qties.keys())
        qties_units = [qty.units for qty in qties.values()]
        qties = list(qties.values())

        def eq(unknowns: np.ndarray) -> np.ndarray:
            qty2_ = Q_(unknowns[0], qties_units[2])
            input_qties_ = {qties_names[0]: qties[0], qties_names[2]: qty2_}
            qty1_new = self.__call__(**input_qties_).state_dict[qties_names[1]]
            out = (qty1_new - qties[1]).to(qties_units[1]).m
            return np.array([out])

        roots = fsolve(eq, np.array([qties[2].to(qties_units[2]).m]))
        qty2 = Q_(roots[0], qties_units[2])
        input_qties = {qties_names[0]: qties[0], qties_names[2]: qty2}
        self._update(**input_qties)

    def _get_quantity(self, qty_name: str) -> Quantity:
        """Returns the state variable denoted by `qty_name`. `qty_name` can be
        any of the keys defined in the dict `_coolprop_qties`."""
        try:
            qty = Q_(
                self._state.keyed_output(self._coolprop_qties[qty_name][0]),  # get the quantity magnitude from Coolprop
                self._coolprop_qties[qty_name][1]  # get the corresponding Coolprop-unit
            )
            return qty
        except ValueError as err:
            warnings.warn(
                f"CoolProp could not solve for quantity '{qty_name}: {err}'",
                category=CoolPropWarning
            )
            return None

    def _get_state(self) -> FluidState:
        """Get the current state of the fluid wrapped in a new `FluidState`
        instance."""
        fluid_attrs = {
            'name': self.fluid_name,
            'backend': self.backend,
            'mass_fractions': self.mass_fractions,
            'vol_fractions': self.vol_fractions,
            'reference': self.reference
        }
        fluid_state = {
            k: self._get_quantity(k)
            for k in self._coolprop_qties.keys()
        }
        return FluidState(fluid_attrs, fluid_state)

    def __call__(self, phase: str | None = None, **input_qties: Quantity) -> FluidState:
        """Pass the input state variables that change the fluid's current state
        and get the new state wrapped in a `FluidState` instance.
        Normally, only two input state variables are needed to define the state.
        But in the case of mixtures, the possible combinations are restricted.
        Only for the following combinations CoolProp can determine the state of
        a mixture:
        - `P` and `T`, or
        - `P` and `x`, or
        - `T` and `x`
        However it is possible to find the state of a mixture for other
        combinations if one of the state variables is P, T or x and the second
        is any other valid state variable. Therefore, it is, however, needed
        that an initial guess is also given for a third state variable, which
        must be either P, T or x, but cannot be equal to the first state
        variable.

        Example
        -------
        ```
        evaporator_in = R454B(
            T=Q_(-10, 'degC'),
            h=condenser_out.h,
            x=Q_(0, 'frac')
        )
        ```
        """
        self._update(**input_qties)
        return self._get_state()

    def __deepcopy__(self, memo):
        # needed this to solve copy-trouble with CoolProp: create a new instance
        # with the same attributes as the instance to be copied.
        new_fluid = type(self)(
            self.fluid_name,
            self.backend,
            self.mass_fractions,
            self.vol_fractions,
            self.reference
        )
        return new_fluid

    @property
    def critical_point(self) -> FluidState:
        """Returns state of fluid at critical point."""
        T_crit = Q_(CP.PropsSI('Tcrit', self.fluid_name), 'K')
        P_crit = Q_(CP.PropsSI('Pcrit', self.fluid_name), 'Pa')
        self._update(T=T_crit, P=P_crit)
        return self._get_state()

    @property
    def triple_point(self) -> FluidState:
        """Returns state of fluid at triple point."""
        T_triple = Q_(CP.PropsSI('Ttriple', self.fluid_name), 'K')
        P_triple = Q_(CP.PropsSI('ptriple', self.fluid_name), 'Pa')
        self._update(T=T_triple, P=P_triple)
        return self._get_state()

    @property
    def coolprop_abstract_state(self):
        """For internal use only."""
        return self._state
