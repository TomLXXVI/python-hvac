"""
Correlations for the average Nusselt number in case of external free convective
fluid flow past a vertical plate.

The correlations were taken from:
Nellis G. F. , & Klein S. A.  (2021).
INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.
"""
from hvac import Quantity
from hvac.fluids import FluidState
from hvac.heat_transfer.free_convection.general import (
    characteristic_velocity,
    reynolds_number,
    grashof_number,
    prandtl_number,
    rayleigh_number
)


def average_nusselt_number(Ra_L: float, Pr: float) -> float:
    """Calculates the average Nusselt number associated with an isothermal,
    vertical heated/cooled plate acc. to the correlation of Churchill and Chu.

    Parameters
    ----------
    Ra_L: float
        Rayleigh number
    Pr: float
        Prandtl number

    Returns
    -------
    Nu_L_avg: float
        Average Nusselt number
    """
    Nu_L_avg = 0.825
    n = 0.387 * Ra_L ** (1 / 6)
    d = (1 + (0.492 / Pr) ** (9 / 16)) ** (8 / 27)
    Nu_L_avg += n / d
    Nu_L_avg = Nu_L_avg ** 2
    return Nu_L_avg


class Plate:
    """Isothermal vertical heated/cooled plate."""

    def __init__(
        self,
        L: Quantity,
        T_surf: Quantity | None,
        T_inf: Quantity,
        fluid: FluidState | None
    ) -> None:
        """Creates `Plate` instance.

        Parameters
        ----------
        L: Quantity
            Vertical dimension of the plate.
        T_surf: Quantity
            Surface temperature of the plate.
        T_inf: Quantity
            Fluid temperature far away from the plate.
        fluid: FluidState
            Object that groups the thermophysical properties of the fluid.
        """
        self.L = L.to('m')
        self.T_surf = T_surf.to('K') if T_surf else None
        self.T_inf = T_inf.to('K')
        self.fluid = fluid

    @property
    def u_char(self) -> Quantity:
        """Characteristic buoyancy velocity."""
        u_char = characteristic_velocity(
            self.L, self.fluid.beta,
            self.T_surf, self.T_inf
        )
        return u_char

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
        """Returns the average heat transfer coefficient between the plate and
        the fluid.
        """
        Nu_L_avg = average_nusselt_number(self.Ra, self.Pr)
        h_avg = Nu_L_avg * self.fluid.k / self.L
        return h_avg.to('W / (m ** 2 * K)')
