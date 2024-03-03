from hvac import Quantity
from hvac.fluids import FluidState

Q_ = Quantity


g = Q_(9.81, 'm / s**2')


def grashof_number(
    L_char: Quantity,
    dT: Quantity,
    fluid: FluidState
) -> float:
    """Dimensionless number that represents the ratio of the buoyancy force to
    the viscous force acting on a unit volume of an air jet.
    """
    beta = fluid.beta
    mu = fluid.mu
    rho = fluid.rho
    nu = mu / rho
    Gr = g * beta * dT * L_char ** 3 / (nu ** 2)
    return Gr.to_base_units().magnitude


def reynolds_number(
    L_char: Quantity,
    U_o: Quantity,
    fluid: FluidState
) -> float:
    """Dimensionless number that represents the ratio of the inertia force to
    the viscous force acting on a unit volume of an air jet.
    """
    mu = fluid.mu
    rho = fluid.rho
    nu = mu / rho
    Re = U_o * L_char / nu
    return Re.to_base_units().magnitude


def archimedes_number(
    L_char: Quantity,
    U: Quantity,
    dT: Quantity,
    fluid: FluidState
) -> float:
    """Dimensionless number that represents the ratio of the buoyancy force to
    the inertia force acting on a unit volume of an air jet.
    """
    Gr = grashof_number(L_char, dT, fluid)
    Re = reynolds_number(L_char, U, fluid)
    Ar = Gr / (Re ** 2)
    return Ar


def room_archimedes_number(
    B: Quantity,
    H: Quantity,
    T_r: Quantity,
    V_dot: Quantity,
    dT_o: Quantity
) -> float:
    """Returns the room Archimedes number.

    Parameters
    ----------
    B:
        Width of the room, i.e. the dimension of the room normal to the
        direction of the air jet.
    H:
        Height of the room.
    T_r:
        Average room temperature in the occupied zone.
    V_dot:
        Volume flow rate of supply air through the outlet.
    dT_o:
        Temperature difference between the supply air and the room air in the
        occupied zone.
    """
    B = B.to('m')
    H = H.to('m')
    T_r = T_r.to('K')
    V_dot = V_dot.to('m**3 / s')
    dT_o = dT_o.to('K')
    n = 2.0 * g * dT_o * (B * H) ** 3
    d = T_r * V_dot ** 2 * (B + H)
    Ar = n / d
    return Ar
