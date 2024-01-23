"""
Correlations for the average Nusselt number in case of internal, free convective,
fluid flow inside a rectangular enclosure formed by two vertical plates with
the separation distance S between the plates being much less than either of the
other two dimensions (length L and width W).

Examples:
    - the glass covering over solar collectors
    - multi-pane windows.

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
    grashof_number,
    prandtl_number,
    rayleigh_number
)


def _average_nusselt_number_vertical1(
    Ra_S: float,
    Pr: float,
    L_to_S: float
) -> float:
    """Average Nusselt number for tilt angle 90° (vertical up) at low Rayleigh
    number and with W >> L.

    Parameters
    ----------
    Ra_S: float
        Rayleigh number associated with the spacing S between the plates of
        the rectangular enclosure (i.e. the thickness of the enclosure).
    Pr: float
        Prandtl number of the fluid.
    L_to_S: float
        The ratio of the length of the enclosure to the spacing (thickness) of
        the enclosure.

    Returns
    -------
    Nu_S_avg: float
        Average Nusselt number.

    Notes
    -----
    The correlation is valid if:
    - 10 < L/S < 40
    - 1 < Pr < 2.e4
    - 1.e4 < Ra_S < 1.e7
    """
    if not (10 < L_to_S < 40):
        warnings.warn(f"Ratio L/S {L_to_S} out of range ]10, 40[")
    if not (1 < Pr < 2.e4):
        warnings.warn(f"Prandtl number {Pr} out of range ]1, 2.e4[")
    if not (1.e4 < Ra_S < 1.e7):
        warnings.warn(f"Rayleigh number {Ra_S} out of range ]1.e4, 1.e7[")
    Nu_S_avg = 0.42 * (Ra_S ** 0.25) * (Pr ** 0.012) * (L_to_S ** -0.3)
    return Nu_S_avg


def _average_nusselt_number_vertical2(
    Ra_S: float,
    Pr: float,
    L_to_S: float
) -> float:
    """Average Nusselt number for tilt angle 90° (vertical up) at higher Rayleigh
    number and with W >> L.

    Parameters
    ----------
    Ra_S: float
        Rayleigh number associated with the spacing S between the plates of
        the rectangular enclosure (i.e. the thickness of the enclosure).
    Pr: float
        Prandtl number of the fluid.
    L_to_S: float
        The ratio of the length of the enclosure to the spacing (thickness) of
        the enclosure.

    Returns
    -------
    Nu_S_avg: float
        Average Nusselt number.

    Notes
    -----
    The correlation is valid if:
    - 1 < L/S < 40
    - 1 < Pr < 20
    - 1.e6 < Ra_S < 1.e9
    """
    if not (1 < L_to_S < 40):
        warnings.warn(f"Ratio L/S {L_to_S} out of range ]1, 40[")
    if not (1 < Pr < 20):
        warnings.warn(f"Prandtl number {Pr} out of range ]1, 20[")
    if not (1.e6 < Ra_S < 1.e9):
        warnings.warn(f"Rayleigh number {Ra_S} out of range ]1.e6, 1.e9[")
    Nu_S_avg = 0.046 * Ra_S ** (1 / 3)
    return Nu_S_avg


def _average_nusselt_number_vertical(
    Ra_S: float,
    Pr: float,
    L_to_S: float
) -> float:
    """Average Nusselt number for tilt angle 90° (vertical up).

    Parameters
    ----------
    Ra_S: float
        Rayleigh number associated with the spacing S between the plates of
        the rectangular enclosure (i.e. the thickness of the enclosure).
    Pr: float
        Prandtl number of the fluid.
    L_to_S: float
        The ratio of the length of the enclosure to the spacing (thickness) of
        the enclosure.

    Returns
    -------
    Nu_S_avg: float
        Average Nusselt number.
    """
    if Ra_S <= 1.e6:
        Nu_S_avg = _average_nusselt_number_vertical1(Ra_S, Pr, L_to_S)
    else:
        Nu_S_avg = _average_nusselt_number_vertical2(Ra_S, Pr, L_to_S)
    return Nu_S_avg


def _average_nusselt_number_tilted1(
    Ra_S: float,
    theta: float,
    L_to_S: float
) -> float:
    """Average Nusselt number for rectangular enclosure acc. to correlation of
    Hollands et al. (1976) when 0 <= tilt angle <= 1.22 rad.

    Parameters
    ----------
    Ra_S: float
        Rayleigh number associated with the spacing S between the plates of
        the rectangular enclosure (i.e. the thickness of the enclosure).
    theta: float
        Tilt angle in radians.
    L_to_S: float
        The ratio of the length of the enclosure to the spacing (thickness) of
        the enclosure.

    Returns
    -------
    Nu_S_avg: float
        Average Nusselt number.

    Notes
    -----
    The correlation is valid if:
    - L/S and W/s > 12
    """
    if not(L_to_S > 12):
        warnings.warn(f"Ratio L/S {L_to_S} is smaller than 12")
    Nu_S_avg = (
        1
        + 1.44 * max(0, 1 - 1708 / (Ra_S * np.cos(theta)))
        * max(0, 1 - 1708 * (np.sin(1.8 * theta) ** 1.6) / (Ra_S * np.cos(theta)))
        + max(0, (Ra_S * np.cos(theta) / 5830) ** (1 / 3) - 1)
    )
    return Nu_S_avg


def _average_nusselt_number_tilted2(
    Ra_S: float,
    Pr: float,
    theta: float,
    L_to_S: float
) -> float:
    """Average Nusselt number for rectangular enclosure acc. to correlation of
    Ayyaswamy and Catton (1973) when 1.22 < tilt angle < pi/2 rad.

    Parameters
    ----------
    Ra_S: float
        Rayleigh number associated with the spacing S between the plates of
        the rectangular enclosure (i.e. the thickness of the enclosure).
    Pr: float
        Prandtl number of the fluid.
    theta: float
        Tilt angle in radians.
    L_to_S: float
        The ratio of the length of the enclosure to the spacing (thickness) of
        the enclosure.

    Returns
    -------
    Nu_S_avg: float
        Average Nusselt number.
    """
    Nu_S_avg = _average_nusselt_number_vertical(Ra_S, Pr, L_to_S)
    Nu_S_avg *= (np.sin(theta) ** 0.25)
    return Nu_S_avg


def _average_nusselt_number_tilted3(
    Ra_S: float,
    Pr: float,
    theta: float,
    L_to_S: float
) -> float:
    """Average Nusselt number for rectangular enclosure acc. to correlation of
    Arnold et al. (1975) when pi/2 < tilt angle < pi rad.

    Parameters
    ----------
    Ra_S: float
        Rayleigh number associated with the spacing S between the plates of
        the rectangular enclosure (i.e. the thickness of the enclosure).
    Pr: float
        Prandtl number of the fluid.
    theta: float
        Tilt angle in radians.
    L_to_S: float
        The ratio of the length of the enclosure to the spacing (thickness) of
        the enclosure.

    Returns
    -------
    Nu_S_avg: float
        Average Nusselt number.
    """
    Nu_S_avg_ver = _average_nusselt_number_vertical(Ra_S, Pr, L_to_S)
    Nu_S_avg = 1 + (Nu_S_avg_ver - 1) * np.sin(theta)
    return Nu_S_avg


def average_nusselt_number(
    Ra_S: float,
    Pr: float,
    theta: float,
    L_to_S: float
) -> float:
    """Average Nusselt number for rectangular enclosure when
    0 <= tilt angle <= pi rad.

    Parameters
    ----------
    Ra_S: float
        Rayleigh number associated with the spacing S between the plates of
        the rectangular enclosure (i.e. the thickness of the enclosure).
    Pr: float
        Prandtl number of the fluid.
    theta: float
        Tilt angle in radians.
    L_to_S: float
        The ratio of the length of the enclosure to the spacing (thickness) of
        the enclosure.

    Returns
    -------
    Nu_S_avg: float
        Average Nusselt number.

    Notes
    -----
    - If tilt angle `theta` = 0 rad, then the cooler surface lies above the hotter
    surface.
    - If tilt angle `theta` = pi/2 rad, then the two surfaces are parallel to
    the gravity vector.
    - If tilt angle `theta` = pi rad, the hotter surface lies above the cooler
    surface.
    """
    if 0 <= theta <= 1.22:
        Nu_S_avg = _average_nusselt_number_tilted1(Ra_S, theta, L_to_S)
    elif 1.22 < theta < np.pi / 2:
        Nu_S_avg = _average_nusselt_number_tilted2(Ra_S, Pr, theta, L_to_S)
    elif theta == np.pi / 2:
        Nu_S_avg = _average_nusselt_number_vertical(Ra_S, Pr, L_to_S)
    else:  # if pi / 2 < theta <= np.pi
        Nu_S_avg = _average_nusselt_number_tilted3(Ra_S, Pr, theta, L_to_S)
    return Nu_S_avg


class Enclosure:
    """Rectangular enclosure formed by two vertical plates with the separation
    distance S between the plates being much less than either of the other two
    dimensions (length L and width W).
    """
    def __init__(
        self,
        L: Quantity,
        S: Quantity,
        theta: Quantity,
        T_hot: Quantity,
        T_cold: Quantity,
        fluid: FluidState | None
    ) -> None:
        """Creates `Enclosure` instance.

        Parameters
        ----------
        L: Quantity
            Length of the rectangular enclosure (i.e. the dimension which is also
            a leg of the tilted angle).
        S: Quantity
            Spacing between the hot and cold plate.
        theta: Quantity
            Tilt angle of the enclosure.
        T_hot: Quantity
            Temperature of the hot side of the enclosure.
        T_cold: Quantity
            Temperature of the cold side of the enclosure.
        fluid: FluidState
            Object that groups the thermophysical properties of the fluid.

        Notes
        -----
        - If tilt angle `theta` = 0°, then the cooler surface lies above the
        hotter surface.
        - If tilt angle `theta` = 90°, then the two surfaces are parallel to
        the gravity vector.
        - If tilt angle `theta` = 180°, the hotter surface lies above the cooler
        surface.
        """
        self.L = L.to('m')
        self.S = S.to('m')
        self.theta = theta.to('rad')
        self.T_hot = T_hot.to('K') if T_hot else None
        self.T_cold = T_cold.to('K') if T_cold else None
        self.fluid = fluid

    def avg_heat_trf_coeff(self) -> Quantity:
        """Returns the average heat transfer coefficient of the convective heat
        transfer between the parallel hot and cold plates due to the fluid
        present in the enclosure.

        Notes
        -----
        If the length `L` based Rayleigh number is less than 1000, the buoyancy
        force is insufficient to overcome the viscous force and the fluid remains
        stagnant. In this limit, heat transfer will happen by conduction only
        and therefore the Nusselt number will be 1.
        """
        Ra_L = self._rayleigh_number_L_based()
        if Ra_L >= 1.e3:
            u_S = characteristic_velocity(
                self.S, self.fluid.beta,
                self.T_hot, self.T_cold
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
            Nu_S_avg = average_nusselt_number(Ra_S, Pr, self.theta.to('rad').m, L_to_S)
            h_avg = Nu_S_avg * self.fluid.k / self.S
            return h_avg.to('W / (m ** 2 * K)')
        else:
            h_avg = self.fluid.k / self.S
            return h_avg.to('W / (m ** 2 * K)')

    def _rayleigh_number_L_based(self) -> float:
        u_L = characteristic_velocity(
            self.L, self.fluid.beta,
            self.T_hot, self.T_cold
        )
        Re_L = reynolds_number(
            self.fluid.rho, self.fluid.mu,
            self.L, u_L
        )
        Pr = prandtl_number(
            self.fluid.rho, self.fluid.mu,
            self.fluid.k, self.fluid.cp
        )
        Gr_L = grashof_number(Re_L)
        Ra_L = rayleigh_number(Gr_L, Pr)
        return Ra_L
