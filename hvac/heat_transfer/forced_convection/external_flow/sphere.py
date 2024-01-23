"""Correlations for the drag coefficient and the average Nusselt number in case
of external forced convective flow of a fluid past a sphere.

The correlations were taken from:
Nellis G. F., & Klein S. A. (2021)
INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.
"""
from math import pi as PI
import warnings
from hvac import Quantity
from hvac.fluids import FluidState
from .general import reynolds_number, prandtl_number, drag_force


def drag_coefficient(Re_D: float) -> float:
    """Calculates the drag coefficient acc. to the correlation of White (1991)

    Parameters
    ----------
    Re_D:
        Reynolds number

    Returns
    -------
    C_D: float
        Drag coefficient

    Notes
    -----
    Valid for Re_D_o < 2.e5
    """
    if Re_D >= 2.e5:
        warnings.warn(f"Reynolds number {Re_D} is greater than 2.e5")
    C_D = 24 / Re_D + 6 / (1 + Re_D ** 0.5) + 0.40
    return C_D


def average_nusselt_number(Re_D: float, Pr: float) -> float:
    """Calculates the average Nusselt number acc. to the correlation of Whitaker
    (1972).

    Parameters
    ----------
    Re_D:
        Reynolds number
    Pr:
        Prandtl number

    Returns
    -------
    Nu_D_avg: float
        Average Nusselt number

    Notes
    -----
    The correlation is valid for 3.5 < Re_D_o < 8.e4 and 0.7 < Pr < 380.
    """
    if not (3.5 < Re_D < 8.e4):
        warnings.warn(f"Reynolds number {Re_D} out of range [3.5, 8.e4]")
    if not (0.7 < Pr < 380):
        warnings.warn(f"Prandtl number {Pr} out of range [0.7, 380]")
    Nu_D_avg = 2 + (0.4 * Re_D ** 0.5 + 0.06 * Re_D ** (2 / 3)) * (Pr ** 0.4)
    return Nu_D_avg


class Sphere:
    """Class for calculating the drag force exerted on a sphere by a moving
    fluid and the average heat transfer coefficient between the sphere and
    the moving fluid.
    """
    def __init__(self, D: Quantity, fluid: FluidState | None) -> None:
        """Creates a `Cylinder` instance.

        Parameters
        ----------
        D: Quantity
            Diameter of the sphere.
        fluid: FluidState
            The fluid moving along the surface of the sphere. A `FluidState`
            object groups the thermophysical properties of a fluid in a given
            state.
        """
        self.D = D.to('m')
        self.fluid = fluid
        self._u_inf = None

    @property
    def u_inf(self) -> Quantity:
        """Returns the free stream velocity of the flowing fluid with respect
        to the sphere.
        """
        return self._u_inf

    @u_inf.setter
    def u_inf(self, v: Quantity) -> None:
        """Sets the free stream velocity of the flowing fluid with respect to the
        sphere.
        """
        self._u_inf = v.to('m / s')

    @property
    def Re(self) -> float:
        Re_D = reynolds_number(
            self.fluid.rho, self.fluid.mu,
            self.u_inf, self.D
        )
        return Re_D

    @property
    def Pr(self) -> float:
        Pr = prandtl_number(
            self.fluid.rho, self.fluid.mu,
            self.fluid.k, self.fluid.cp
        )
        return Pr

    @property
    def volume(self) -> Quantity:
        return 4 * PI / 3 * (self.D / 2) ** 3

    @property
    def area(self) -> Quantity:
        return 4 * PI * (self.D / 2) ** 2

    def drag_force(self) -> Quantity:
        """Returns the drag force exerted by the moving fluid on the sphere."""
        C_D = drag_coefficient(self.Re)
        F_D = drag_force(C_D, self.fluid.rho, self.u_inf, PI * self.D ** 2 / 4)
        return F_D

    def avg_heat_trf_coeff(self) -> Quantity:
        """Returns the average heat transfer coefficient between the sphere and
        the moving fluid.
        """
        Nu_D_avg = average_nusselt_number(self.Re, self.Pr)
        h_avg = Nu_D_avg * self.fluid.k / self.D
        return h_avg.to('W / (m ** 2 * K)')
