"""
Correlations for the average Nusselt number in case of external free convective
fluid flow past a horizontal plate.

The correlations were taken from:
Nellis G. F. , & Klein S. A.  (2021).
INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.
"""
import warnings
import numpy as np
from hvac import Quantity
from hvac.fluids import FluidState
from hvac.heat_transfer.free_convection.general import (
    characteristic_velocity,
    reynolds_number,
    rayleigh_number,
    prandtl_number,
    grashof_number
)


class LaminarFlow:

    @staticmethod
    def average_nusselt_number(Ra: float, Pr: float) -> float:
        """Calculates the average Nusselt number for laminar flow.

        Parameters
        ----------
        Ra: float
            Raleigh number.
        Pr: float
            Prandtl number.

        Returns
        -------
        Nu_L_avg: float
            Average Nusselt number.
        """
        C_lam = 0.671 / ((1 + (0.492 / Pr) ** (9 / 16)) ** (4 / 9))
        n = 1.4
        d1 = 0.835 * C_lam * (Ra ** 0.25)
        d2 = np.log(1 + n / d1)
        Nu_L_avg = n / d2
        return Nu_L_avg


class TurbulentFlow:

    @staticmethod
    def average_nusselt_number(Ra: float, Pr: float) -> float:
        """Calculates the average Nusselt number for turbulent flow.

        Parameters
        ----------
        Ra: float
            Raleigh number.
        Pr: float
            Prandtl number.

        Returns
        -------
        Nu_L_avg: float
            Average Nusselt number.
        """
        n = 1 + 0.0107 * Pr
        d = 1 + 0.01 * Pr
        C_tur = 0.14 * (n / d)
        Nu_L_avg = C_tur * (Ra ** (1 / 3))
        return Nu_L_avg


def average_nusselt_number_upwards(Ra: float, Pr: float) -> float:
    """Calculates the average Nusselt number for a horizontal plate, facing
    upward when heating, or facing downward when cooling acc. to correlation
    of Churchill and Usagi.

    Parameters
    ----------
    Ra: float
        Raleigh number.
    Pr: float
        Prandtl number.

    Returns
    -------
    Nu_L_avg: float
        Average Nusselt number.

    Notes
    -----
    The correlation is valid for 1.0 < Ra < 1.e10.
    """
    if not (1.0 < Ra < 1.e10):
        warnings.warn(f"Rayleigh number {Ra} out of range [1, 1.e10]")
    Nu_L_avg_lam = LaminarFlow.average_nusselt_number(Ra, Pr)
    Nu_L_avg_tur = TurbulentFlow.average_nusselt_number(Ra, Pr)
    Nu_L_avg = ((Nu_L_avg_lam ** 10) + (Nu_L_avg_tur ** 10)) ** (1 / 10)
    return Nu_L_avg


def average_nusselt_number_downwards(Ra: float, Pr: float) -> float:
    """Calculates the average Nusselt number for a horizontal plate, facing
    downward when heating, or facing upward when cooling.

    Parameters
    ----------
    Ra: float
        Raleigh number.
    Pr: float
        Prandtl number.

    Returns
    -------
    Nu_L_avg: float
        Average Nusselt number.

    Notes
    -----
    The correlation is valid for 1.e3 < Ra < 1.e10.
    """
    if not (1.e3 < Ra < 1.e10):
        warnings.warn(f"Rayleigh number {Ra} out of range [1.e3, 1.e10]")
    d2 = (1 + (1.9 / Pr) ** 0.9) ** (2 / 9)
    d1 = 2.5 / (0.527 * (Ra ** 0.20))
    d = np.log(1 + d1 * d2)
    n = 2.5
    Nu_L_avg = n / d
    return Nu_L_avg


class Plate:
    """Horizontal plate."""

    def __init__(
        self,
        L: Quantity,
        W: Quantity | None,
        T_surf: Quantity | None,
        T_inf: Quantity,
        fluid: FluidState | None,
        configuration: str = 'upwards-heated'
    ) -> None:
        """Creates `Plate` instance.

        Parameters
        ----------
        L: Quantity
            Length or diameter of plate
        W: Quantity
            Width of plate. Set to `None` if the plate has a circular shape.
        T_surf: Quantity
            Surface temperature of the plate.
        T_inf: Quantity
            Fluid temperature far away from the plate.
        fluid: FluidState
            Object that groups the thermophysical properties of the fluid.
        configuration: str
            - 'upwards-heated' (default)
            - 'downwards-cooled'
            - 'downwards-heated'
            - 'upwards-cooled'
        """
        self.L = L.to('m')
        self.W = W.to('m')
        self.T_surf = T_surf.to('K') if T_surf else None
        self.T_inf = T_inf.to('K')
        self.fluid = fluid
        self.configuration = configuration

    @property
    def L_char(self) -> Quantity:
        """Characteristic length of the plate."""
        if self.W is None:
            A = np.pi * (self.L ** 2) / 4
            P = np.pi * self.L
        else:
            A = self.L * self.W
            P = 2 * (self.L + self.W)
        L_char = A / P
        return L_char

    @property
    def u_char(self) -> Quantity:
        """Characteristic buoyancy velocity."""
        u_char = characteristic_velocity(
            self.L_char, self.fluid.beta,
            self.T_surf, self.T_inf
        )
        return u_char

    @property
    def Re(self) -> float:
        """Reynolds number."""
        Re = reynolds_number(
            self.fluid.rho, self.fluid.mu,
            self.L_char, self.u_char
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
        if self.configuration in ['upwards-heated', 'downwards-cooled']:
            Nu_L_avg = average_nusselt_number_upwards(self.Ra, self.Pr)
        else:
            Nu_L_avg = average_nusselt_number_downwards(self.Ra, self.Pr)
        h_avg = Nu_L_avg * self.fluid.k / self.L_char
        return h_avg.to('W / (m ** 2 * K)')
