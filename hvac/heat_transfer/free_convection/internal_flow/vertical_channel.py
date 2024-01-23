"""
Correlation for the average Nusselt number in case of internal, free convective,
fully developed fluid flow through an isothermal, vertical channel formed by
two vertical plates.

The correlations were taken from:
Nellis G. F. , & Klein S. A.  (2021).
INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.
"""
import numpy as np
from hvac import Quantity
from hvac.fluids import FluidState
from hvac.heat_transfer.free_convection.general import (
    characteristic_velocity,
    reynolds_number,
    grashof_number,
    prandtl_number,
    rayleigh_number
)


def average_nusselt_number(Ra_S: float, L_to_S: float) -> float:
    """Calculates the average Nusselt number associated with fully developed
    flow induced in an isothermal, vertical channel formed by two vertical
    plates acc. to the correlation of Elenbaas.

    Parameters
    ----------
    Ra_S: float
        Rayleigh number based on plate spacing S
    L_to_S: float
        Ratio of plate length (height) to plate spacing

    Returns
    -------
    Nu_S_avg: float
        Average Nusselt number
    """
    Nu_S_avg = (
        (Ra_S / 24) * (1 / L_to_S)
        * (1 - np.exp(-(35 / Ra_S) * L_to_S)) ** 0.75
    )
    return Nu_S_avg


class Channel:
    """Vertical channel formed by two parallel plates.

    Notes
    -----
    If the plates are sufficiently short, such that the flow does not become
    fully developed, then the flow can be modeled using the correlations for
    a vertical plate.
    """

    def __init__(
        self,
        L: Quantity,
        S: Quantity,
        T_surf: Quantity | None,
        T_inf: Quantity,
        fluid: FluidState | None
    ) -> None:
        self.L = L.to('m')
        self.S = S.to('m')
        self.T_surf = T_surf.to('K') if T_surf else None
        self.T_inf = T_inf.to('K')
        self.fluid = fluid

    def avg_heat_trf_coeff(self) -> Quantity:
        """Returns the average heat transfer coefficient between the parallel
        plates and the fluid when the flow is fully developed.
        """
        u_S = characteristic_velocity(
            self.S, self.fluid.beta,
            self.T_surf, self.T_inf
        )
        Re_S = reynolds_number(
            self.fluid.rho, self.fluid.mu,
            self.S, u_S
        )
        Pr = prandtl_number(
            self.fluid.rho, self.fluid.mu,
            self.fluid.k, self.fluid.cp
        )
        Gr_S = grashof_number(Re_S)
        Ra_S = rayleigh_number(Gr_S, Pr)
        L_to_S = self.L.to('m') / self.S.to('m')
        Nu_S_avg = average_nusselt_number(Ra_S, L_to_S)
        h_avg = Nu_S_avg * self.fluid.k / self.S
        return h_avg.to('W / (m ** 2 * K)')
