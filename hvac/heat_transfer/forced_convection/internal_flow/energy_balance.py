from typing import Callable
import numpy as np
from scipy.integrate import quad, solve_ivp
from hvac import Quantity

Q_ = Quantity


def get_fluid_temp_from_heat_flx(
    x: Quantity,
    q_fun: Callable[[Quantity], Quantity],
    cp: Quantity,
    m_dot: Quantity,
    Ph: Quantity,
    T_fl_in: Quantity
) -> Quantity:
    x_f = x.to('m').m
    cp = cp.to('J / (kg * K)').m
    m_dot = m_dot.to('kg / s').m
    Ph = Ph.to('m').m
    T_fl_in = T_fl_in.to('K').m

    def _q(x: float) -> float:
        x = Q_(x, 'm')
        q_val = q_fun(x)
        return q_val.to('W / m ** 2').m

    I = quad(_q, 0, x_f)[0]
    T_fl = T_fl_in + Ph * I / (m_dot * cp)
    return Q_(T_fl, 'K')


def get_wall_temp_from_heat_flx(
    q: Quantity,
    h: Quantity,
    T_fl: Quantity
) -> Quantity:
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
    Ph: Quantity,
    T_fl_in: Quantity
) -> Quantity:
    x_f = x.to('m').m
    cp = cp.to('J / (kg * K)').m
    h = h.to('W / (m ** 2 * K)').m
    m_dot = m_dot.to('kg / s').m
    Ph = Ph.to('m').m
    T_fl_in = T_fl_in.to('K').m

    def _T_w(x: float) -> float:
        x = Q_(x, 'm')
        T_w_val = T_w_fun(x)
        return T_w_val.to('K').m

    def f(x: float, T_fl: list[float]):
        return Ph * h * (_T_w(x) - T_fl[0]) / (m_dot * cp)

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
    x = x.to('m')
    T_ext = T_ext.to('K')
    cp = cp.to('J / (kg * K)')
    m_dot = m_dot.to('kg / s')
    R_tot = R_tot.to('K / W')
    T_fl_in = T_fl_in.to('K')
    L = L.to('m')
    T_fl = T_ext - (T_ext - T_fl_in) * np.exp(-x / (m_dot * cp * R_tot * L))
    return T_fl
