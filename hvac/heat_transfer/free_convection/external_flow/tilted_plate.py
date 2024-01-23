"""
Correlations for the average Nusselt number in case of external free convective
fluid flow past a tilted plate.

The correlations were taken from:
Nellis G. F. , & Klein S. A.  (2021).
INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.
"""
import numpy as np
from hvac import Quantity
from hvac.fluids import FluidState
from hvac.heat_transfer.free_convection.general import (
    g,
    characteristic_velocity,
    reynolds_number,
    rayleigh_number,
    prandtl_number,
    grashof_number
)
from ..external_flow.vertical_plate import average_nusselt_number as ver_avg_nusselt_number
from ..external_flow.horizontal_plate import average_nusselt_number_upwards as hor_uw_avg_nusselt_number
from ..external_flow.horizontal_plate import average_nusselt_number_downwards as hor_dw_avg_nusselt_number


class Plate:
    """Tilted plate."""

    def __init__(
        self,
        L: Quantity,
        W: Quantity | None,
        theta: Quantity,
        T_surf: Quantity | None,
        T_inf: Quantity,
        fluid: FluidState | None,
    ) -> None:
        """Creates `Plate` instance.

        Parameters
        ----------
        L: Quantity
            Length or diameter of plate
        W: Quantity
            Width of plate. Set to `None` if the plate has a circular shape.
        theta: Quantity
            Tilt angle of the plate with respect to the horizontal plane.
        T_surf: Quantity
            Surface temperature of the plate.
        T_inf: Quantity
            Fluid temperature far away from the plate.
        fluid: FluidState
            Object that groups the thermophysical properties of the fluid.

        Notes
        -----
        - When the tilt angle = 0°, the plate is heated upward facing.
        - When the tilt angle = 90°, the plate is vertical.
        - When the tilt angle = 180°, the plate is heated downward facing.
        """
        self.L = L.to('m')
        self.W = W.to('m')
        self.theta = theta.to('rad')
        self.T_surf = T_surf.to('K') if T_surf else None
        self.T_inf = T_inf.to('K')
        self.fluid = fluid

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
        u_char_ver = characteristic_velocity(
            self.L_char, self.fluid.beta, self.T_surf, self.T_inf,
            g_=g * np.sin(self.theta)
        )
        if self.theta < np.pi / 2:
            u_char_hor = characteristic_velocity(
                self.L_char, self.fluid.beta, self.T_surf, self.T_inf,
                g_=g * np.cos(self.theta)
            )
        else:
            u_char_hor = characteristic_velocity(
                self.L_char, self.fluid.beta, self.T_surf, self.T_inf,
                g_=g * np.cos(np.pi - self.theta)
            )
        u_char = max(u_char_ver, u_char_hor)
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

    def _avg_heat_trf_coeff_vertical(self) -> Quantity:
        # vertical plate with g * np.sin(self.theta)
        u_char = characteristic_velocity(
            self.L_char, self.fluid.beta, self.T_surf, self.T_inf,
            g_=g * np.sin(self.theta)
        )
        Re = reynolds_number(self.fluid.rho, self.fluid.mu, self.L_char, u_char)
        Gr = grashof_number(Re)
        Ra = rayleigh_number(Gr, self.Pr)
        Nu_L_avg = ver_avg_nusselt_number(Ra, self.Pr)
        h_avg = Nu_L_avg * self.fluid.k / self.L_char
        return h_avg.to('W / (m ** 2 * K)')

    def _avg_heat_trf_coeff_horizontal_up(self) -> Quantity:
        # horizontal plate upwards with g * np.cos(self.theta)
        u_char = characteristic_velocity(
            self.L_char, self.fluid.beta, self.T_surf, self.T_inf,
            g_=g * np.cos(self.theta)
        )
        Re = reynolds_number(self.fluid.rho, self.fluid.mu, self.L_char, u_char)
        Gr = grashof_number(Re)
        Ra = rayleigh_number(Gr, self.Pr)
        Nu_L_avg = hor_uw_avg_nusselt_number(Ra, self.Pr)
        h_avg = Nu_L_avg * self.fluid.k / self.L_char
        return h_avg.to('W / (m ** 2 * K)')

    def _avg_heat_trf_coeff_horizontal_down(self) -> Quantity:
        # horizontal plate downwards with g * np.cos(pi - self.theta)
        u_char = characteristic_velocity(
            self.L_char, self.fluid.beta, self.T_surf, self.T_inf,
            g_=g * np.cos(np.pi - self.theta)
        )
        Re = reynolds_number(self.fluid.rho, self.fluid.mu, self.L_char, u_char)
        Gr = grashof_number(Re)
        Ra = rayleigh_number(Gr, self.Pr)
        Nu_L_avg = hor_dw_avg_nusselt_number(Ra, self.Pr)
        h_avg = Nu_L_avg * self.fluid.k / self.L_char
        return h_avg.to('W / (m ** 2 * K)')

    def avg_heat_trf_coeff(self) -> Quantity:
        """Returns the average heat transfer coefficient between the tilted
        plate and the fluid.
        """
        h_avg_ver = self._avg_heat_trf_coeff_vertical()
        if self.theta < np.pi / 2:
            h_avg_hor = self._avg_heat_trf_coeff_horizontal_up()
        else:
            h_avg_hor = self._avg_heat_trf_coeff_horizontal_down()
        h_avg = max(h_avg_ver, h_avg_hor)
        return h_avg
