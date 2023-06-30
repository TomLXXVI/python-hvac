from hvac import Quantity

Q_ = Quantity

Re_crit = 5e5


def reynolds_number(
    rho: Quantity,
    mu: Quantity,
    u_inf: Quantity,
    L_char: Quantity
) -> float:
    """Calculates the Reynolds number.

    Parameters
    ----------
    rho: Quantity
        Mass density of the fluid.
    mu: Quantity
        Absolute or dynamic viscosity of the fluid.
    u_inf: Quantity
        Free stream velocity of the flowing fluid.
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
    u_inf = u_inf.to('m / s').m
    L_char = L_char.to('m').m
    Re = rho * u_inf * L_char / mu
    return Re


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


def get_flow_condition(Re: float) -> str:
    """Returns the flow condition of a fluid depending on the Reynolds number.
    The flow condition is either laminar or turbulent.
    """
    if Re < Re_crit:
        return 'laminar'
    else:
        return 'turbulent'


def shear_stress(Cf: float, rho: Quantity, u_inf: Quantity) -> Quantity:
    """Calculates the shear stress exerted by the moving fluid along the surface
    of the object.

    Parameters
    ----------
    Cf: float
        Local or average friction coefficient.
    rho: Quantity
        Mass density of the fluid.
    u_inf: Quantity
        Free stream velocity of the flowing fluid.

    Returns
    -------
    tau_s: Quantity
        Local or average shear stress, depending on whether the local or average
        friction coefficient is given.
    """
    rho = rho.to('kg / m ** 3')
    u_inf = u_inf.to('m / s')
    tau_s = Cf * rho * (u_inf ** 2) / 2
    return tau_s.to('Pa')


def drag_force(Cd: float, rho: Quantity, u_inf: Quantity, A_p: Quantity) -> Quantity:
    """Calculates the drag force exerted by the moving fluid on the object.

    Parameters
    ----------
    Cd: float
        Drag coefficient.
    rho:
        Mass density of fluid.
    u_inf:
        Free stream velocity of the flowing fluid.
    A_p: Quantity
        Projected area of the object perpendicular to the direction of the flow.

    Returns
    -------
    Fd: Quantity
        Drag force.
    """
    rho = rho.to('kg / m ** 3')
    u_inf = u_inf.to('m / s')
    A_p = A_p.to('m ** 2')
    Fd = Cd * rho * (u_inf ** 2) * A_p / 2
    return Fd.to('N')
