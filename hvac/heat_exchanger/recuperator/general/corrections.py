"""
Nusselt number and friction factor corrections.

References
----------
[1] Shah, R. K., & Sekulic, D. P. Fundamentals of Heat Exchanger Design
    (John Wiley & Sons, 2003).
"""
import math
from hvac import Quantity
from hvac.fluids import FluidState, HumidAir


Q_ = Quantity


# property ratio method exponents according to [1], table 7.12
laminar_flow = {
    'heating': {
        'gas': {
            'n': 0.0,
            'm': 1.00
        },
        'liquid': {
            'n': -0.14,
            'm': 0.58
        }
    },
    'cooling': {
        'gas': {
            'n': 0.0,
            'm': 0.81
        },
        'liquid': {
            'n': -0.14,
            'm': 0.54
        }
    }
}


# property ratio method exponents according to [1], table 7.13
turbulent_flow = {
    'heating': {
        'gas': {
            'n': lambda Tw_on_Tm: -math.log10(Tw_on_Tm) ** (1 / 4) + 0.3,
            'm': -0.1
        },
        'liquid': {
            'n': -0.11,
            'm': 0.25
        }
    },
    'cooling': {
        'gas': {
            'n': 0,
            'm': -0.1
        },
        'liquid': {
            'n': -0.25,
            'm': 0.24
        }
    }
}


def correct_nusselt_number(
    Nu_cp: float,
    T_w: Quantity,
    flow_regime: str,
    thermal_regime: str,
    fluid: FluidState | HumidAir,
) -> float:
    """Returns the corrected Nusselt number acc. to [1], eqs. 7.157 and 7.158.

    Parameters
    ----------
    Nu_cp:
        Nusselt number evaluated at bulk temperature.
    T_w:
        Wall temperature.
    flow_regime: {'laminar', 'turbulent'}
        Flow regime of fluid.
    thermal_regime: {'cooling', 'heating'}
        Thermal regime of fluid.
    fluid:
        Thermophysical state of fluid.
    """
    if isinstance(fluid, FluidState):
        if 'liquid' in fluid.phase:
            fluid_phase = 'liquid'
        elif 'gas' in fluid.phase:
            fluid_phase = 'gas'
        else:
            raise ValueError('the fluid should be in a liquid or gaseous phase.')
    else:
        # in case of humid air
        fluid_phase = 'gas'
    if flow_regime == 'laminar':
        n = laminar_flow[thermal_regime][fluid_phase]['n']
    else:  # flow regime is turbulent
        n = turbulent_flow[thermal_regime][fluid_phase]['n']
    if fluid_phase == 'gas':
        if isinstance(fluid, FluidState):
            T_m = fluid.T.to('K').m
        else:
            T_m = fluid.Tdb.to('K').m
        T_w = T_w.to('K').m
        if flow_regime == 'turbulent' and thermal_regime == 'heating':
            n = n(T_w / T_m)
        Nu = (T_w / T_m) ** n * Nu_cp
    else:  # fluid phase is liquid
        mu_m = fluid.mu.to('Pa * s').m
        fluid_w = fluid.fluid(T=T_w, P=Q_(1, 'atm'))
        mu_w = fluid_w.mu.to('Pa * s').m
        Nu = (mu_w / mu_m) ** n * Nu_cp
    return Nu


def correct_friction_factor(
    f_cp: float,
    T_w: Quantity,
    flow_regime: str,
    thermal_regime: str,
    fluid: FluidState | HumidAir
) -> float:
    """Returns the corrected friction factor acc. to [1], eqs. 7.157 and 7.158.

    Parameters
    ----------
    f_cp:
        Friction factor evaluated at bulk temperature.
    T_w:
        Wall temperature.
    flow_regime: {'laminar', 'turbulent'}
        Flow regime of fluid.
    thermal_regime: {'cooling', 'heating'}
        Thermal regime of fluid.
    fluid:
        Thermophysical state of fluid.
    """
    if isinstance(fluid, FluidState):
        if 'liquid' in fluid.phase:
            fluid_phase = 'liquid'
        elif 'gas' in fluid.phase:
            fluid_phase = 'gas'
        else:
            raise ValueError('the fluid should be in a liquid or gaseous phase.')
    else:
        # fluid is humid air = always a gas
        fluid_phase = 'gas'
    if flow_regime == 'laminar':
        m = laminar_flow[thermal_regime][fluid_phase]['m']
    else:  # flow regime is turbulent
        m = turbulent_flow[thermal_regime][fluid_phase]['m']
    if fluid_phase == 'gas':
        T_m = fluid.T.to('K').m if isinstance(fluid, FluidState) else fluid.Tdb.to('K').m
        T_w = T_w.to('K').m
        f = (T_w / T_m) ** m * f_cp
    else:  # fluid phase is liquid
        mu_m = fluid.mu.to('Pa * s').m
        if isinstance(fluid, FluidState):
            fluid_w = fluid.fluid(T=T_w, P=Q_(1, 'atm'))
        else:
            fluid_w = HumidAir(Tdb=T_w, P=Q_(1, 'atm'))
        mu_w = fluid_w.mu.to('Pa * s').m
        f = (mu_w / mu_m) ** m * f_cp
    return f
