"""
Correlations for the average Nusselt number in case of external free convective
fluid flow past a vertical cylinder.

The correlations were taken from:
Nellis G. F. , & Klein S. A.  (2021).
INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.
"""
import numpy as np
from hvac import Quantity
from hvac.fluids import FluidState
from .vertical_plate import average_nusselt_number as vp_average_nusselt_number
from ..general import (
    characteristic_velocity,
    reynolds_number,
    rayleigh_number,
    prandtl_number,
    grashof_number
)


def average_nusselt_number(Ra_L: float, Pr: float, L_to_D: float) -> float:
    """Calculates the average Nusselt number acc. to the correlation of
    Sparrow and Gregg.

    Parameters
    ----------
    Ra_L: float
        Rayleigh number
    Pr: float
        Prandtl number
    L_to_D: float
        Ratio of cylinder length to diameter.

    Returns
    -------
    Nu_L_avg: float
        Average Nusselt number
    """
    Nu_L_avg_vp = vp_average_nusselt_number(Ra_L, Pr)
    zeta = 1.8 / Nu_L_avg_vp * L_to_D
    k = zeta / np.log(1 + zeta)
    Nu_L_avg = k * Nu_L_avg_vp
    return Nu_L_avg


class Cylinder:
    """Vertical cylinder."""

    def __init__(
        self,
        D: Quantity,
        L: Quantity,
        T_surf: Quantity | None,
        T_inf: Quantity,
        fluid: FluidState | None
    ) -> None:
        self.D = D.to('m')
        self.L = L.to('m')
        self.T_surf = T_surf.to('K') if T_surf else None
        self.T_inf = T_inf.to('K')
        self.fluid = fluid

    @property
    def u_char(self) -> Quantity:
        """Characteristic buoyancy velocity."""
        u_L = characteristic_velocity(
            self.L, self.fluid.beta,
            self.T_surf, self.T_inf
        )
        return u_L

    @property
    def Re(self) -> float:
        """Reynolds number."""
        Re = reynolds_number(
            self.fluid.rho, self.fluid.mu,
            self.L, self.u_char
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
        """Returns the average heat transfer coefficient between the cylinder
        and the fluid.
        """
        L_to_D = self.L.to('m') / self.D.to('m')
        Nu_L_avg = average_nusselt_number(self.Ra, self.Pr, L_to_D)
        h_avg = Nu_L_avg * self.fluid.k / self.L
        return h_avg.to('W / (m ** 2 * K)')
