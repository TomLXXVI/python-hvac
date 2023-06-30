import numpy as np
from hvac import Quantity

Q_ = Quantity

g = Q_(9.81, 'm / s ** 2')


def characteristic_velocity(
    L_char: Quantity,
    beta: Quantity,
    T_surf: Quantity,
    T_inf: Quantity,
    g_: Quantity | None = None
) -> Quantity:
    """Returns the characteristic buoyancy velocity for a natural convection
    problem

    Parameters
    ----------
    L_char: Quantity
        Characteristic length of the object.
    beta: Quantity
        Volumetric thermal expansion coefficient, evaluated at the average
        fluid temperature. For an ideal gas, it is the inverse of its absolute
        temperature.
    T_surf: Quantity
        Surface temperature of the object.
    T_inf: Quantity
        Fluid temperature far away from the object.
    g_: Quantity, default None
        Projection of the gravity vector; used with tilted plates.

    Returns
    -------
    u_char: Quantity
        Characteristic buoyancy velocity.
    """
    if g_ is None: g_ = g
    u_char = np.sqrt(g_ * L_char * beta * (T_surf.to('K') - T_inf.to('K')))
    return u_char.to('m / s')


def reynolds_number(
    rho: Quantity,
    mu: Quantity,
    L_char: Quantity,
    u_char: Quantity
) -> float:
    """Calculates the Reynolds number.

    Parameters
    ----------
    rho: Quantity
        Mass density of the fluid.
    mu: Quantity
        Absolute or dynamic viscosity of the fluid.
    u_char: Quantity
        Characteristic buoyancy velocity.
    L_char: Quantity
        Characteristic length of the object over which the flowing fluid is
        passing.

    Returns
    -------
    Re: float
        Reynolds number.
    """
    rho = rho.to('kg / m ** 3').m
    mu = mu.to('Pa * s').m
    L_char = L_char.to('m').m
    u_char = u_char.to('m / s').m
    Re = rho * L_char * u_char / mu
    return Re


def grashof_number(Re: float) -> float:
    """Calculates the Grashof number.

    Parameters
    ----------
    Re: float
        Reynolds number

    Returns
    -------
    Gr: float
        Grashof number
    """
    Gr = Re ** 2
    return Gr


def prandtl_number(
    rho: Quantity,
    mu: Quantity,
    k: Quantity,
    cp: Quantity
) -> float:
    """Calculates the Prandtl number of a fluid.

    Parameters
    ----------
    rho: Quantity
        Mass density of the fluid.
    mu: Quantity
        Absolute or dynamic viscosity of the fluid.
    k: Quantity
        Thermal conductivity.
    cp: Quantity
        Specific heat of the fluid (at constant pressure).

    Returns
    -------
    Pr: float
        Prandtl number.
    """
    rho = rho.to('kg / m ** 3').m
    mu = mu.to('Pa * s').m
    k = k.to('W / (m * K)').m
    cp = cp.to('J / (kg * K)').m
    nu = mu / rho
    alpha = k / (rho * cp)
    Pr = nu / alpha
    return Pr


def rayleigh_number(Gr: float, Pr: float) -> float:
    """Calculates the Rayleigh number.

    Parameters
    ----------
    Gr: float
        Grashof number
    Pr: float
        Prandtl number

    Returns
    -------
    Ra: float
        Rayleigh number
    """
    Ra = Gr * Pr
    return Ra
