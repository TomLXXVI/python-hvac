"""
Correlations for the average Nusselt number in case of external free convective
fluid flow past a horizontal cylinder.

The correlations were taken from:
Nellis G. F. , & Klein S. A.  (2021).
INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.
"""
import warnings
from hvac import Quantity
from hvac.fluids import FluidState
from hvac.heat_transfer.free_convection.general import (
    characteristic_velocity,
    reynolds_number,
    rayleigh_number,
    prandtl_number,
    grashof_number
)


def average_nusselt_number(Ra_D: float, Pr: float) -> float:
    """Calculates the average Nusselt number acc. to the correlation of
    Churchill and Chu.

    Parameters
    ----------
    Ra_D: float
        Rayleigh number
    Pr: float
        Prandtl number

    Returns
    -------
    Nu_D_avg: float
        Average Nusselt number

    Notes
    -----
    The correlation is valid for:
    - 1.e-5 < Ra < 1.e12
    """
    if not (1.e-5 < Ra_D < 1.e12):
        warnings.warn(f"Rayleigh number {Ra_D} is out of range ]1.e-5, 1.e12[")
    Nu_D_avg = 0.60
    n = 0.387 * Ra_D ** (1 / 6)
    d = (1 + (0.559 / Pr) ** (9 / 16)) ** (8 / 27)
    Nu_D_avg += n / d
    Nu_D_avg = Nu_D_avg ** 2.0
    return Nu_D_avg


class Cylinder:
    """Horizontal cylinder."""

    def __init__(
        self,
        D: Quantity,
        L: Quantity,
        T_surf: Quantity | None,
        T_inf: Quantity,
        fluid: FluidState | None
    ) -> None:
        self.D = D
        self.L = L
        self.T_surf = T_surf.to('K') if T_surf else None
        self.T_inf = T_inf.to('K')
        self.fluid = fluid

    @property
    def u_char(self) -> Quantity:
        """Characteristic buoyancy velocity."""
        u_char = characteristic_velocity(
            self.D, self.fluid.beta,
            self.T_surf, self.T_inf
        )
        return u_char

    @property
    def Re(self) -> float:
        """Reynolds number."""
        Re = reynolds_number(
            self.fluid.rho, self.fluid.mu,
            self.D, self.u_char
        )
        return Re

    @property
    def Pr(self) -> float:
        """Prandtl number."""
        Pr = prandtl_number(
            self.fluid.rho, self.fluid.mu,
            self.fluid.k, self.fluid.cp
        )
        return Pr

    @property
    def Gr(self) -> float:
        """Grashof number."""
        Gr = grashof_number(self.Re)
        return Gr

    @property
    def Ra(self) -> float:
        """Rayleigh number."""
        Ra = rayleigh_number(self.Gr, self.Pr)
        return Ra

    def avg_heat_trf_coeff(self) -> Quantity:
        """Returns the average heat transfer coefficient between the cylinder and
        the fluid.
        """
        Nu_D_avg = average_nusselt_number(self.Ra, self.Pr)
        h_avg = Nu_D_avg * self.fluid.k / self.D
        return h_avg.to('W / (m ** 2 * K)')
 