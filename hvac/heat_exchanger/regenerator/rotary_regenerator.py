"""
HEAT TRANSFER ANALYSIS OF A COUNTERFLOW ROTARY REGENERATOR USING THE
EPSILON-NTUo METHOD (ONLY FOR SENSIBLE HEAT TRANSFER).

Class `CounterFlowRotaryRegenerator` can be used to calculate the effectiveness
of a counterflow rotary regenerator, given the needed specifications of the
matrix, and of the hot and cold fluid (gas) at the hot/cold entrance of the
regenerator. Besides the effectiveness, the class also returns the heat
transfer rate between the hot and cold fluid, and the outlet temperatures of
the hot and the cold fluid.

References
----------
Shah, R. K., & Sekulic, D. P. (2003). Fundamentals of Heat Exchanger Design.
John Wiley & Sons. Chapter 5: Thermal design theory for regenerators, § 5.2.
"""
import warnings
from dataclasses import dataclass
import numpy as np
from hvac import Quantity
from . import epsilon_ntu
from .definitions import Fluid, Matrix


class RegeneratorWarning(Warning):
    pass


@dataclass
class Fluid(Fluid):
    """Groups the input data about the hot/cold fluid.

    Parameters
    ----------
    m_dot:
        Mass flow rate of fluid.
    cp:
        Specific heat of fluid.
    h:
        Heat transfer coefficient between fluid and matrix.
    T_i:
        Fluid temperature at the entry of the matrix.

    Attributes
    ----------
    C_dot:
        Heat capacity rate.
    """
    m_dot: Quantity
    cp: Quantity
    h: Quantity
    T_i: Quantity

    def __post_init__(self):
        self.C_dot = Fluid.heat_capacity_rate(self.m_dot, self.cp)


@dataclass
class Matrix(Matrix):
    """Groups the input data about the matrix.

    Parameters
    ----------
    D: optional
        Outer wheel (rotor) diameter.
    d: optional
        Inner wheel diameter, i.e. the shaft diameter.
    L:
        Rotor length.
    rho: optional
        Density of matrix material.
    k: optional
        Heat conduction transfer coefficient of the matrix material. This
        is used to determine the effect of longitudinal heat conduction between
        the hot and cold side of the regenerator on the effectiveness of the
        regenerator.
    sigma: optional
        Matrix porosity.
    beta:
        Packing density (surface-to-volume ratio).
    seal_fraction:
        Fraction of rotor face area covered by radial seals.
    theta_h:
        Disk sector angle or fraction of the frontal area through which hot
        fluid flows .
    theta_c:
        Disk sector angle or fraction of the frontal area through which cold
        fluid flows.
    M: optional
        Matrix mass. If not specified, parameters `rho` and `sigma` must be
        specified.
    A_fr: optional.
        Frontal area of the matrix. If not specified, parameters `D` and `d`
        must be specified.
    num_disks: optional
        Number of parallel disks the rotary regenerator is made off. The
        frontal area is multiplied with this number.
    t_w: optional
        Wall thickness of the matrix channels. This is used to calculate the
        thermal resistance of the matrix wall in case the effect of transverse
        wall heat conduction needs to be taken into account.

    Attributes
    ----------
    C_r_dot:
        Matrix wall heat capacity rate during the cold/hot period.
    A:
        Total matrix heat transfer surface area.
    A_h:
        Heat transfer area of the hot sector of the regenerator matrix wall.
    A_c:
        Heat transfer area of the cold sector of the regenerator matrix wall.
    A_kt:
        Total surface area for longitudinal heat conduction from the hot to the
        cold side of the regenerator.
    """
    L: Quantity
    c: Quantity
    N: Quantity
    beta: Quantity
    seal_fraction: Quantity
    theta_h: Quantity
    theta_c: Quantity
    D: Quantity | None = None
    d: Quantity | None = None
    A_fr: Quantity | None = None
    M: Quantity | None = None
    rho: Quantity | None = None
    sigma: Quantity | None = None
    num_disks: int = 1
    k: Quantity | None = None
    t_w: Quantity | None = None

    def __post_init__(self):
        if self.A_fr is None:
            if all((self.D, self.d)) is False:
                raise ValueError(
                    'the frontal area is not specified and '
                    'cannot be determined'
                )
            self.A_fr = np.pi * (self.D ** 2 - self.d ** 2) / 4 * self.num_disks
        if self.M is None:
            if all((self.rho, self.sigma)) is False:
                raise ValueError(
                    'the matrix mass is not specified and '
                    'cannot be determined'
                )
            self.M = Matrix.mass(self.A_fr, self.L, self.rho, self.sigma)
        self.C_r_dot = Matrix.heat_capacity_rate(
            m_w=self.M,
            c_w=self.c,
            N=self.N
        )
        self.A = Matrix.total_heat_transfer_area(
            self.A_fr, self.L, self.beta,
            self.seal_fraction
        )
        self.A_h = Matrix.sector_heat_transfer_area(
            A=self.A,
            theta_j=self.theta_h
        )
        self.A_c = Matrix.sector_heat_transfer_area(
            A=self.A,
            theta_j=self.theta_c
        )
        if self.sigma is None:
            self.sigma = Matrix.porosity(self.A, self.A_fr)
        if self.k is not None:
            self.A_kt = epsilon_ntu.A_kt(self.A_fr, self.sigma)


class CounterFlowRotaryRegenerator:
    """Represents a sensible heat recovery wheel.

    This class can be used to determine the sensible heat transfer rate from the
    hot to the cold fluid between both sides of the wheel when inlet
    temperatures of the hot and cold fluid are given.
    """
    def __init__(
        self,
        hot_fluid: Fluid,
        cold_fluid: Fluid,
        matrix: Matrix
    ) -> None:
        self.hot_fluid = hot_fluid
        self.cold_fluid = cold_fluid
        self.matrix = matrix
        self.C_dot_min = min(self.hot_fluid.C_dot, self.cold_fluid.C_dot)

    def _determine_C_dot_star(self) -> Quantity:
        C_dot_max = max(self.hot_fluid.C_dot, self.cold_fluid.C_dot)
        C_dot_star = epsilon_ntu.C_dot_star(self.C_dot_min, C_dot_max)
        return C_dot_star

    def _determine_hA_star(self) -> Quantity:
        hA_star = epsilon_ntu.hA_star(
            self.hot_fluid.h,
            self.hot_fluid.C_dot,
            self.matrix.theta_h,
            self.cold_fluid.h,
            self.cold_fluid.C_dot,
            self.matrix.theta_c
        )
        return hA_star

    def _determine_C_r_dot_star(self) -> Quantity:
        C_r_dot_star = epsilon_ntu.C_r_dot_star(
            self.matrix.C_r_dot,
            self.C_dot_min
        )
        return C_r_dot_star

    def _determine_NTU_o(self) -> Quantity:
        UA_o = epsilon_ntu.UA_o(
            self.hot_fluid.h,
            self.cold_fluid.h,
            self.matrix.A_h,
            self.matrix.A_c,
            self.matrix.t_w,
            self.matrix.k
        )
        NTU_o = epsilon_ntu.NTU_o(UA_o, self.C_dot_min)
        return NTU_o

    @property
    def eps(self) -> Quantity:
        """Regenerator effectiveness."""
        hA_star = self._determine_hA_star()
        if not 0.25 <= hA_star.m < 4.0:
            warnings.warn(
                "The convection conductance ratio (hA)* is outside the "
                "range [0.25, 4.0[. Calculation of regenerator effectiveness "
                "may not be accurate.",
                category=RegeneratorWarning
            )
        NTU_o = self._determine_NTU_o()
        C_dot_star = self._determine_C_dot_star()
        C_r_dot_star = self._determine_C_r_dot_star()
        if self.matrix.k is not None:
            lamda = epsilon_ntu.lamda(
                self.matrix.k, self.matrix.A_kt,
                self.matrix.L, self.C_dot_min
            )
            Fi = epsilon_ntu.Fi(lamda, NTU_o)
            C_lamda = epsilon_ntu.C_lamda(NTU_o, lamda, Fi)
        else:
            C_lamda = None
        eps = epsilon_ntu.eps(NTU_o, C_dot_star, C_r_dot_star, C_lamda)
        return eps

    @property
    def Q_dot_max(self) -> Quantity:
        """Theoretical maximum heat transfer rate in the regenerator."""
        Q_dot_max = epsilon_ntu.Q_dot_max(
            self.C_dot_min,
            self.hot_fluid.T_i,
            self.cold_fluid.T_i
        )
        return Q_dot_max

    @property
    def Q_dot(self) -> Quantity:
        """Heat transfer rate between hot and cold fluid in the regenerator."""
        Q_dot = epsilon_ntu.Q_dot(self.eps, self.Q_dot_max)
        return Q_dot

    @property
    def T_hot_out(self) -> Quantity:
        """Hot fluid outlet temperature."""
        T_h_o = self.hot_fluid.T_i.to('K') - self.Q_dot / self.hot_fluid.C_dot
        return T_h_o.to('K')

    @property
    def T_cold_out(self) -> Quantity:
        """Cold fluid outlet temperature."""
        T_c_o = self.cold_fluid.T_i.to('K') + self.Q_dot / self.cold_fluid.C_dot
        return T_c_o.to('K')
