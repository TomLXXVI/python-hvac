from typing import Callable
import numpy as np
from scipy.integrate import quad, solve_ivp
from hvac import Quantity

Q_ = Quantity


def get_fluid_temp_from_heat_flux(
    x: Quantity,
    q_fun: Callable[[Quantity], Quantity],
    cp: Quantity,
    m_dot: Quantity,
    P_h: Quantity,
    T_fl_in: Quantity
) -> Quantity:
    """Returns the mean fluid temperature at a given position x along a tube, if
    the heat flux as a function of position is given.
    
    Parameters
    ----------
    x:
        Position along the tube where the fluid temperature has to be 
        determined. At the entrance of the tube x = 0.
    q_fun:
        Function that takes the position x along the tube and returns the 
        heat flux at this position.
    cp:
        The (average) specific heat of the fluid.
    m_dot:
        The mass flow rate of the fluid through the tube.
    P_h:
        The perimeter associated with the heat flux (this is not necessarily the
        full perimeter of the tube).
    T_fl_in:
        Temperature of the fluid entering the tube.
    """
    x_f = x.to('m').m
    cp = cp.to('J / (kg * K)').m
    m_dot = m_dot.to('kg / s').m
    P_h = P_h.to('m').m
    T_fl_in = T_fl_in.to('K').m

    # "De-quantify" the heat flux function to use it with Scipy's integrating
    # function:
    def _q(x: float) -> float:
        x = Q_(x, 'm')
        q_val = q_fun(x)
        return q_val.to('W / m ** 2').m

    # Take the integral of the heat-flux function from the entrance x = 0 to
    # the final position x_f (SI-units W / m):
    I = quad(_q, 0, x_f)[0]
    # Multiplying with the perimeter leads to the total heat transfer along the
    # tube's surface (SI-units W):
    Q_dot = I * P_h
    # The fluid temperature at position x can now be determined from the energy
    # balance:
    T_fl = T_fl_in + Q_dot / (m_dot * cp)
    return Q_(T_fl, 'K')


def get_wall_temp_from_heat_flux(
    q: Quantity,
    h: Quantity,
    T_fl: Quantity
) -> Quantity:
    """Returns the (inner) surface wall temperature at a given position along
    the tube.

    Parameters
    ----------
    q:
        The heat flux at the given position.
    h:
        The heat transfer coefficient at the given position.
    T_fl:
        The fluid temperature.
    """
    h = h.to('W / (m ** 2 * K)')
    q = q.to('W / m ** 2')
    T_fl = T_fl.to('K')
    T_w = T_fl + q / h
    return T_w


def get_fluid_temp_from_wall_temp(
    x: Quantity,
    T_w_fun: Callable[[Quantity], Quantity],
    cp: Quantity,
    h: Quantity,
    m_dot: Quantity,
    P_h: Quantity,
    T_fl_in: Quantity
) -> Quantity:
    """Returns the (mean) fluid temperature at a given position x along a tube,
    if the wall temperature as a function of position is given.

    Parameters
    ----------
    x:
        Position along the tube where the fluid temperature has to be
        determined. At the entrance of the tube x = 0.
    T_w_fun:
        Function that takes the position x along the tube and returns the
        wall temperature at this position.
    cp:
        The (average) specific heat of the fluid.
    h:
        The average heat transfer coefficient of the fluid between the entrance
        of the tube (x = 0) and the position x.
    m_dot:
        The mass flow rate of the fluid through the tube.
    P_h:
        The perimeter associated with the heat flux (this is not necessarily the
        full perimeter of the tube).
    T_fl_in:
        Temperature of the fluid entering the tube.
    """
    x_f = x.to('m').m
    cp = cp.to('J / (kg * K)').m
    h = h.to('W / (m ** 2 * K)').m
    m_dot = m_dot.to('kg / s').m
    P_h = P_h.to('m').m
    T_fl_in = T_fl_in.to('K').m

    def _T_w(x: float) -> float:
        x = Q_(x, 'm')
        T_w_val = T_w_fun(x)
        return T_w_val.to('K').m

    def f(x: float, T_fl: list[float]):
        return P_h * h * (_T_w(x) - T_fl[0]) / (m_dot * cp)

    sol = solve_ivp(f, [0, x_f], [T_fl_in], t_eval=[x_f])
    T_fl = Q_(sol.y.flatten()[0], 'K')
    return T_fl


def get_fluid_temp_from_ext_const_temp(
    x: Quantity,
    T_ext: Quantity,
    cp: Quantity,
    m_dot: Quantity,
    R_tot: Quantity,
    T_fl_in: Quantity,
    L: Quantity
) -> Quantity:
    """Returns the (mean) fluid temperature at a given position x along a tube
    if a constant external temperature is given.

    Parameters
    ----------
    x:
        Position along the tube where the fluid temperature has to be
        determined. At the entrance of the tube x = 0.
    T_ext:
        The spatial constant external temperature.
    cp:
        The (average) specific heat of the fluid.
    m_dot:
        The mass flow rate of the fluid through the tube.
    R_tot:
        The (average) total thermal resistance between the external temperature
        and the mean fluid temperature inside the tube along the entrance of the
        tube (x = 0) to the given position x (i.e. along a distance L).
    T_fl_in:
        Temperature of the fluid entering the tube.
    L:
        The distance between the entrance of the tube (x = 0) to the position x.
    """
    x = x.to('m')
    T_ext = T_ext.to('K')
    cp = cp.to('J / (kg * K)')
    m_dot = m_dot.to('kg / s')
    R_tot = R_tot.to('K / W')
    T_fl_in = T_fl_in.to('K')
    L = L.to('m')
    T_fl = T_ext - (T_ext - T_fl_in) * np.exp(-x / (m_dot * cp * R_tot * L))
    return T_fl
