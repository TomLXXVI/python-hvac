"""
In this module the correlation of Dobson and Chato (1998) is implemented for
estimating the local heat transfer coefficient for flow condensation inside a
horizontal tube.

This correlation is known to over-predict the heat transfer coefficient for
refrigerants that condense at high pressure (e.g. R125, R32, and R410a).

The correlation was taken from:
Nellis G.F., & Klein S.A. (2021) INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.
"""

import numpy as np
from hvac import Quantity
from hvac.fluids import Fluid
from ..boiling.flow_boiling import reynolds_number, prandtl_number, froude_number

Q_ = Quantity


def _lockhart_martinelli_parameter(
    x: float,
    rho_f: float,
    mu_f: float,
    rho_g: float,
    mu_g: float
) -> float:
    a = (rho_g / rho_f) ** 0.5
    b = (mu_f / mu_g) ** 0.1
    c = ((1 - x) / x) ** 0.9
    X_tt = a * b * c
    return X_tt


def _modified_froude_number(Re_f: float, Ga: float, X_tt: float) -> float:
    if Re_f <= 1250:
        Fr_mod = 0.025 * (Re_f ** 1.59) / (Ga ** 0.5)
    else:
        Fr_mod = 1.26 * (Re_f ** 1.04) / (Ga ** 0.5)
    z = ((1 + 1.09 * (X_tt ** 0.039)) / X_tt) ** 1.5
    Fr_mod *= z
    return Fr_mod


def _galileo_number(
    D: float,
    rho_f: float,
    mu_f: float,
    rho_g: float,
    g: Quantity = Q_(9.81, 'm / s ** 2')
) -> float:
    g = g.to('m / s ** 2').m
    n = g * rho_f * (rho_f - rho_g) * (D ** 3)
    d = mu_f ** 2
    Ga = n / d
    return Ga


def void_fraction(x: float, rho_f: float, rho_g: float) -> float:
    v_f = 1 / (1 + (1 - x) / x * (rho_g / rho_f) ** (2 / 3))
    return v_f


def _nusselt_number_forced_convection(
    Re_f: float,
    Pr_f: float,
    Fr: float,
    X_tt: float
) -> float:
    if Fr <= 0.7:
        C1 = 4.172 + 5.48 * Fr - 1.564 * Fr ** 2
        C2 = 1.773 - 0.169 * Fr
    else:
        C1 = 7.242
        C2 = 1.655
    z = np.sqrt(1.376 + C1 / (X_tt ** C2))
    Nu_fc = 0.0195 * (Re_f ** 0.8) * (Pr_f ** 0.4) * z
    return Nu_fc


def _h_wavy_flow(
    D: float,
    G: float,
    T_sat: float,
    T_surf: float,
    k_f: float,
    cp_f: float,
    Pr_f: float,
    mu_g: float,
    h_fg: float,
    X_tt: float,
    Ga: float,
    vf: float,
    Nu_fc: float
) -> float:
    h_fg_star = h_fg + 0.68 * cp_f * (T_sat - T_surf)
    A = np.arccos(2 * vf - 1) / np.pi
    a = k_f / D
    b = 0.23 / (1 + 1.11 * (X_tt ** 0.58))
    c = (G * D / mu_g) ** 0.12
    d = (h_fg_star / (cp_f * (T_sat - T_surf))) ** 0.25
    e = (Ga ** 0.25) * (Pr_f ** 0.25)
    f = A * Nu_fc
    h = a * (b * c * d * e + f)
    return h


def _h_annular_flow(
    D: float,
    k_f: float,
    Re_f: float,
    Pr_f: float,
    X_tt: float
) -> float:
    a = 0.023 * (Re_f ** 0.8) * (Pr_f ** 0.4)
    b = 1 + 2.22 / (X_tt ** 0.89)
    Nu_D = a * b
    h = k_f * Nu_D / D
    return h


class HorizontalTube:

    def __init__(
        self,
        D: Quantity,
        fluid: Fluid,
    ) -> None:
        self._D = D.to('m').m
        self.fluid = fluid
        self._m_dot: float | None = None
        self._G: float | None = None

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
        self._G = 4 * self._m_dot / (np.pi * self._D ** 2)  # mass velocity = mass flux

    def heat_trf_coeff(
        self,
        x: Quantity,
        T_sat: Quantity,
        T_surf: Quantity,
        g: Quantity = Q_(9.81, 'm / s ** 2')
    ) -> Quantity:
        if self._G is not None:
            x = x.to('frac').m
            liq_sat = self.fluid(T=T_sat, x=Q_(0.0, 'frac'))
            vap_sat = self.fluid(T=T_sat, x=Q_(1.0, 'frac'))
            k_f = liq_sat.k.to('W / (m * K)').m
            mu_f = liq_sat.mu.to('Pa * s').m
            rho_f = liq_sat.rho.to('kg / m ** 3').m
            cp_f = liq_sat.cp.to('J / (kg * K)').m
            Re_f = reynolds_number(self._G, x, self._D, mu_f)
            Pr_f = prandtl_number(rho_f, mu_f, k_f, cp_f)
            h_f = liq_sat.h.to('J / kg').m
            rho_g = vap_sat.rho.to('kg / m ** 3').m
            mu_g = vap_sat.mu.to('Pa * s').m
            h_g = vap_sat.h.to('J / kg').m
            h_fg = h_g - h_f
            X_tt = _lockhart_martinelli_parameter(x, rho_f, mu_f, rho_g, mu_g)
            Ga = _galileo_number(self._D, rho_f, mu_f, rho_g, g)
            Fr_mod = _modified_froude_number(Re_f, Ga, X_tt)
            if self._G < 500 and Fr_mod < 20:
                # wavy flow
                _T_sat = T_sat.to('K').m
                _T_surf = T_surf.to('K').m
                vf = void_fraction(x, rho_f, rho_g)
                Fr = froude_number(self._G, rho_f, self._D, g)
                Nu_fc = _nusselt_number_forced_convection(Re_f, Pr_f, Fr, X_tt)
                h = _h_wavy_flow(
                    self._D, self._G, _T_sat, _T_surf, k_f,
                    cp_f, Pr_f, mu_g, h_fg, X_tt, Ga, vf, Nu_fc
                )
            else:
                # annular turbulent flow
                h = _h_annular_flow(self._D, k_f, Re_f, Pr_f, X_tt)
            return Q_(h, 'W / (m ** 2 * K)')
        return None
