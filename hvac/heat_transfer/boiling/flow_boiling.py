"""
In this module the Shah correlation is implemented for estimating the local heat
transfer coefficient for flow boiling inside a tube.

The Shah correlation was developed for saturated flow boiling at subcritical
heat fluxes and can be used for a wide range of vapor qualities, from saturated
liquid (x = 0) to the liquid-deficient and dry-out regimes that occur at
qualities of 0.8 or higher.

The correlation was taken from:
Nellis G.F., & Klein S.A. (2021) INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.
"""
import numpy as np
from scipy.interpolate import interp1d
from hvac import Quantity
from hvac.fluids import Fluid
from ..forced_convection.internal_flow import turbulent_flow as single_phase


Q_ = Quantity


def convection_number(x: float, rho_f: float, rho_g: float) -> float:
    Co = ((1 / x - 1) ** 0.8) * ((rho_g / rho_f) ** 0.5)
    return Co


def boiling_number(q: float, G: float, h_fg: float) -> float:
    Bo = q / (G * h_fg)
    return Bo


def froude_number(
    G: float,
    rho_f: float,
    Dh: float,
    g: Quantity = Q_(9.81, 'm / s ** 2')
) -> float:
    g = g.to('m / s ** 2').m
    Fr = (G ** 2) / ((rho_f ** 2) * g * Dh)
    return Fr


def N(Co: float, Fr: float, tube_orient: str) -> float:
    if tube_orient == 'vertical':
        return Co
    if tube_orient == 'horizontal' and Fr > 0.04:
        return Co
    if tube_orient == 'horizontal' and Fr <= 0.04:
        return 0.38 * Co * (Fr ** -0.3)
    return None


def h_tilde_cb(N: float) -> float:
    return 1.8 * (N ** -0.8)


def h_tilde_nb(Bo: float) -> float:
    if Bo >= 0.3e-4:
        return 230 * (Bo ** 0.5)
    else:
        return 1 + 46 * (Bo ** 0.5)


def h_tilde_bs1(Bo: float, N: float) -> float:
    k = (Bo ** 0.5) * np.exp(2.74 * (N ** -0.1))
    if Bo >= 11.0e-4:
        return 14.70 * k
    else:
        return 15.43 * k


def h_tilde_bs2(Bo: float, N: float) -> float:
    k = (Bo ** 0.5) * np.exp(2.74 * (N ** -0.15))
    if Bo >= 11.0e-4:
        return 14.70 * k
    else:
        return 15.43 * k


def h_tilde(N: float, h_cb: float, h_bs1: float, h_bs2: float, h_nb: float) -> float:
    if N <= 0.1:
        return max(h_cb, h_bs2)
    elif 0.1 < N <= 1.0:
        return max(h_cb, h_bs1)
    else:  # if N > 1.0:
        return max(h_cb, h_nb)


def h_l(f_f: float, Re_f: float, Pr_f: float, k_f: float, Dh: float) -> float:
    n = (f_f / 8) * (Re_f - 1000.0) * Pr_f
    d = 1 + 12.7 * ((Pr_f ** (2 / 3)) - 1) * ((f_f / 8) ** 0.5)
    h_l = (n / d) * (k_f / Dh)
    return h_l


def h(h_tilde: float, h_l: float) -> float:
    return h_tilde * h_l


def reynolds_number(G: float, x: float, Dh: float, mu_f: float) -> float:
    Re_f = G * (1 - x) * Dh / mu_f
    return Re_f


def prandtl_number(rho_f: float, mu_f: float, k_f: float, cp_f: float) -> float:
    nu_f = mu_f / rho_f
    alpha_f = k_f / (rho_f * cp_f)
    Pr = nu_f / alpha_f
    return Pr


def friction_factor(Re_f: float) -> float:
    n = 1.0
    d = (0.790 * np.log(Re_f) - 1.64) ** 2.0
    return n / d


class Tube:

    def __init__(
        self,
        Dh: Quantity,
        A: Quantity,
        fluid: Fluid,
        tube_orient: str = 'horizontal'
    ) -> None:
        """
        Creates `Tube` instance.

        Parameters
        ----------
        Dh: PlainQuantity
            Hydraulic diameter of the tube.
        A: PlainQuantity
            The cross-sectional area of the tube.
        fluid: Fluid
            The fluid that runs through the tube.
        tube_orient: str, {'vertical', 'horizontal' (default)}
            The orientation of the tube.
        """
        self.Dh = Dh.to('m').m
        self.A = A.to('m ** 2').m
        self.fluid = fluid
        self.tube_orient = tube_orient
        self._m_dot: float | None = None
        self._G: float | None = None
        self.Re = None

    @property
    def m_dot(self) -> Quantity:
        """Returns the mass flow rate of fluid through the tube."""
        if self._m_dot is not None:
            return Q_(self._m_dot, 'kg / s')
        return None

    @m_dot.setter
    def m_dot(self, v: Quantity) -> None:
        """Set the mass flow rate of fluid through the tube."""
        self._m_dot = v.to('kg / s').m
        self._G = self._m_dot / self.A  # mass velocity = mass flux

    def _heat_trf_coeff_Shah(
        self,
        x: float,
        T_sat: Quantity,
        q_s: Quantity,
        g: Quantity = Q_(9.81, 'm / s ** 2')
    ) -> float:
        # implementation of Shah correlation
        q_s = q_s.to('W / m ** 2').m
        sat_liq = self.fluid(T=T_sat, x=Q_(0.0, 'frac'))
        sat_vap = self.fluid(T=T_sat, x=Q_(1.0, 'frac'))
        Co = convection_number(
            x,
            sat_liq.rho.to('kg / m ** 3').m,
            sat_vap.rho.to('kg / m ** 3').m
        )
        Bo = boiling_number(
            q_s,
            self._G,
            sat_vap.h.to('J / kg').m - sat_liq.h.to('J / kg').m
        )
        Fr = froude_number(
            self._G,
            sat_liq.rho.to('kg / m ** 3').m,
            self.Dh,
            g
        )
        N_ = N(Co, Fr, tube_orient=self.tube_orient)
        h_tilde_ = h_tilde(
            N_,
            h_tilde_cb(N_),
            h_tilde_bs1(Bo, N_),
            h_tilde_bs2(Bo, N_),
            h_tilde_nb(Bo)
        )
        Re = reynolds_number(
            self._G,
            x,
            self.Dh,
            sat_liq.mu.to('Pa * s').m
        )
        Pr = prandtl_number(
            sat_liq.rho.to('kg / m ** 3').m,
            sat_liq.mu.to('Pa * s').m,
            sat_liq.k.to('W / (m * K)').m,
            sat_liq.cp.to('J / (kg * K)').m
        )
        f_f = friction_factor(Re)
        h_l_ = h_l(
            f_f,
            Re,
            Pr,
            sat_liq.k.to('W / (m * K)').m,
            self.Dh
        )
        h_ = h(h_tilde_, h_l_)
        return max(0.0, h_)

    def _heat_trf_coeff_sat_vapor(self, T_sat: Quantity) -> float:
        # returns the heat transfer coefficient for saturated vapor single-phase
        # flow assuming fully developed turbulent flow and smooth tube.
        sat_vap = self.fluid(T=T_sat, x=Q_(1.0, 'frac'))
        rho_g = sat_vap.rho.to('kg / m ** 3').m
        mu_g = sat_vap.mu.to('Pa * s').m
        k_g = sat_vap.k.to('W / (m * K)').m
        cp_g = sat_vap.cp.to('J / (kg * K)').m
        nu_g = mu_g / rho_g
        alpha_g = k_g / (rho_g * cp_g)
        u_g = self._G / rho_g
        Re_g = u_g * self.Dh * rho_g / mu_g
        Pr_g = nu_g / alpha_g
        f_fd = single_phase.friction_factor.SmoothTube.fully_developed_local_friction_factor(Re_g)
        Nu = single_phase.nusselt_number.fully_developed_local_nusselt_number(f_fd, Re_g, Pr_g)
        h = Nu * k_g / self.Dh
        return h

    def heat_trf_coeff(
        self,
        x: Quantity,
        T_sat: Quantity,
        q_s: Quantity,
        g: Quantity = Q_(9.81, 'm / s ** 2')
    ) -> Quantity:
        """
        Get the local heat transfer coefficient for flow boiling acc. to the
        Shah correlation.

        Parameters
        ----------
        x: PlainQuantity
            The vapor quality of the fluid.
        T_sat: PlainQuantity
            The saturation temperature of the fluid.
        q_s: PlainQuantity
            The heat flux at the tube wall added to the fluid.
        g: PlainQuantity, default 9.81 m/s²
            Acceleration of gravity

        Returns
        -------
        h: PlainQuantity
            The local heat transfer coefficient.
        """
        x = x.to('frac').m
        x_min = 1.0e-12
        if x == 0.0: x = x_min
        sat_liq = self.fluid(T=T_sat, x=Q_(0.0, 'frac'))
        self.Re = reynolds_number(self._G, x, self.Dh, sat_liq.mu.to('Pa * s').m)
        if self.Re <= 2300:
            # calculate the vapor quality that corresponds with Re = 2300
            x_2300 = max(x_min, 1 - 2300 * sat_liq.mu.to('Pa * s').m / (self._G * self.Dh))
            h_2300 = self._heat_trf_coeff_Shah(x_2300, T_sat, q_s, g)
            # calculate h for single-phase flow of saturated vapor
            h_sat_vap = self._heat_trf_coeff_sat_vapor(T_sat)
            # find h by linear interpolation between the 2 previous values
            interp = interp1d([x_2300, 1.0], [h_2300, h_sat_vap])
            h = interp(x)
            return Q_(h, 'W / (m ** 2 * K)')
        else:
            h = self._heat_trf_coeff_Shah(x, T_sat, q_s, g)
            return Q_(h, 'W / (m ** 2 * K)')
