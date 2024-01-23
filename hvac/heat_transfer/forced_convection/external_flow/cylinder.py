"""
Correlations for the drag coefficient and the average Nusselt number of
forced convective cross flow of a fluid past a cylinder.

The correlations were taken from:
Nellis G. F., & Klein S. A. (2021)
INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.
"""

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
    """
    C_D = (
        1.18
        + 6.8 / (Re_D ** 0.89)
        + 1.96 / (Re_D ** 0.5)
        - 0.0004 * Re_D / (1 + 3.64e-7 * Re_D ** 2)
    )
    return C_D


def average_nusselt_number(Re_D: float, Pr: float) -> float:
    """Calculates the average Nusselt number acc. to the correlation of Churchill
    and Bernstein (1977).

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
    The correlation is valid for Re_D_o * Pr >= 0.2.
    """
    Nu_D_avg = 0.3
    n = 0.62 * (Re_D ** 0.5) * (Pr ** (1 / 3))
    d = (1 + (0.4 / Pr) ** (2 / 3)) ** 0.25
    f = (1 + (Re_D / 2.82e5) ** 0.625) ** 0.80
    Nu_D_avg += (n / d) * f
    return Nu_D_avg


class Cylinder:
    """Class for calculating the drag force exerted on a cylinder by a moving
    fluid and the average heat transfer coefficient between the cylinder and
    the moving fluid.
    """
    def __init__(self, D: Quantity, L: Quantity, fluid: FluidState | None) -> None:
        """Creates a `Cylinder` instance.

        Parameters
        ----------
        D: Quantity
            Outer diameter of the cylinder.
        L: Quantity
            Length of the cylinder.
        fluid: FluidState
            The fluid moving along the surface of the cylinder. A `FluidState`
            object groups the thermophysical properties of a fluid in a given
            state.
        """
        self.D = D.to('m')
        self.L = L.to('m')
        self.fluid = fluid
        self._u_inf = None

    @property
    def u_inf(self) -> Quantity:
        """Returns the free stream velocity of the flowing fluid with respect
        to the cylinder.
        """
        return self._u_inf

    @u_inf.setter
    def u_inf(self, v: Quantity) -> None:
        """Sets the free stream velocity of the flowing fluid with respect to the
        cylinder.
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

    def drag_force(self) -> Quantity:
        """Returns the drag force exerted by the moving fluid on the cylinder."""
        C_D = drag_coefficient(self.Re)
        F_D = drag_force(C_D, self.fluid.rho, self.u_inf, self.L * self.D)
        return F_D

    def avg_heat_trf_coeff(self) -> Quantity:
        """Returns the average heat transfer coefficient between the cylinder and
        the moving fluid.
        """
        Nu_D_avg = average_nusselt_number(self.Re, self.Pr)
        h_avg = Nu_D_avg * self.fluid.k / self.D
        return h_avg.to('W / (m ** 2 * K)')
