"""Correlations for the pressure drop and the average heat transfer coefficient
for external flow across a bank of tubes (e.g. in a shell-and-tube heat
exchanger). A distinction needs to be made between an inline tube bank and a
staggered tube bank.

The correlations were taken from:
Nellis G. F. , & Klein S. A.  (2021).
INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.

Friction factor data and tube row correction factors were taken from third-party
Python library `ht` (see https://github.com/CalebBell/ht).
"""
from abc import ABC, abstractmethod
import numpy as np
from scipy.interpolate import RectBivariateSpline
from hvac import Quantity
from hvac.fluids import FluidState
from hvac.heat_transfer.forced_convection.external_flow.general import prandtl_number
from ._tube_bank_data import (
    InlineTubeBankFrictionData,
    StaggeredTubeBankFrictionData
)

Q_ = Quantity


class TubeBank(ABC):

    def __init__(
        self,
        N_l: int,
        D_o: Quantity,
        s_l: Quantity,
        s_t: Quantity,
        fluid: FluidState
    ) -> None:
        """Creates a `TubeBank` instance.

        Parameters
        ----------
        N_l: int
            Number of tube rows in the flow direction (longitudinal direction).
        D_o: Quantity
            Outer tube diameter.
        s_l: Quantity
            Longitudinal pitch, i.e. the tube-to-tube spacing in the flow
            direction.
        s_t: Quantity
            Transverse pitch, i.e. the tube-to-tube spacing perpendicular to
            the flow direction.
        fluid: FluidState
            Fluid that flows across the bank of tubes.
        """
        self.N_l = N_l
        self.D_o = D_o.to('m')
        self.s_l = s_l.to('m')
        self.s_t = s_t.to('m')
        self.fluid = fluid
        self._u_inf = None
        self._T_surf = None

    @property
    def approach_velocity(self) -> Quantity:
        """Uniform free stream velocity of fluid approaching the tube bundle."""
        return self._u_inf

    @approach_velocity.setter
    def approach_velocity(self, v: Quantity) -> None:
        self._u_inf = v.to('m / s')

    @property
    def surface_temperature(self) -> Quantity:
        """Surface temperature of the tube."""
        return self._T_surf

    @surface_temperature.setter
    def surface_temperature(self, v: Quantity) -> None:
        self._T_surf = v.to('K')

    @property
    @abstractmethod
    def _u_max(self) -> Quantity:
        """Maximum flow velocity of fluid between the tubes."""
        ...

    def _reynolds_number(self) -> float:
        """Returns Reynolds number."""
        rho = self.fluid.rho.to('kg / m ** 3').m
        mu = self.fluid.mu.to('Pa * s').m
        D_o = self.D_o.m
        u_max = self._u_max.m
        Re = rho * D_o * u_max / mu
        return Re

    @abstractmethod
    def _friction_factor(self) -> float:
        """Returns the average friction factor acc. to Zukauska (1987)."""
        ...

    @property
    def pressure_drop(self) -> Quantity:
        """Get the pressure drop across the tube bank."""
        f = self._friction_factor()
        rho = self.fluid.rho.to('kg / m ** 3').m
        u_max = self._u_max.m
        dP = 0.5 * self.N_l * f * rho * (u_max ** 2)
        return Q_(dP, 'Pa')

    def _nusselt_number(self, C: float, m: float) -> float:
        """Returns the average Nusselt number acc. to Zukauskas (1972) /
        Bergman et al. (2011).
        """
        Pr = prandtl_number(  # Prandtl number at the bulk temperature of the fluid
            self.fluid.rho,
            self.fluid.mu,
            self.fluid.k,
            self.fluid.cp
        )
        Re = self._reynolds_number()
        Nu = C * (Re ** m) * (Pr ** 0.36)
        if self._T_surf is not None:
            fluid_surf = self.fluid.fluid(T=self._T_surf, P=self.fluid.P)
            Pr_surf = prandtl_number(  # Prandtl number at the surface temperature of the tube.
                fluid_surf.rho,
                fluid_surf.mu,
                fluid_surf.k,
                fluid_surf.cp
            )
            Nu *= (Pr / Pr_surf) ** 0.25
        cf = self._tube_row_correction()
        Nu *= cf
        return Nu

    @abstractmethod
    def _tube_row_correction(self) -> float:
        ...

    @abstractmethod
    def avg_heat_trf_coeff(self) -> Quantity:
        """Average heat transfer coefficient."""
        ...


class InlineTubeBank(TubeBank):
    """Class that represents a tube bank with an inline tube arrangement
    (each tube is directly behind the tube in the adjacent row).
    """
    def __init__(
        self,
        N_l: int,
        D_o: Quantity,
        s_l: Quantity,
        s_t: Quantity,
        fluid: FluidState
    ) -> None:
        """Creates an `InlineTubeBank` instance.

        Parameters
        ----------
        N_l: int
            Number of tube rows in the flow direction (longitudinal direction).
        D_o: Quantity
            Outer tube parameter
        s_l: Quantity
            Longitudinal pitch, i.e. the tube-to-tube spacing in the flow
            direction.
        s_t: Quantity
            Transverse pitch, i.e. the tube-to-tube spacing perpendicular to
            the flow direction.
        fluid: FluidState
            Fluid that flows across the bank of tubes.
        """
        super().__init__(N_l, D_o, s_l, s_t, fluid)
        # interpolation function returning friction factor for inline square
        self._f_interp = RectBivariateSpline(
            x=np.array(list(InlineTubeBankFrictionData.f_ax_dict.keys())),
            y=InlineTubeBankFrictionData.Re_ax,
            z=np.array(list(InlineTubeBankFrictionData.f_ax_dict.values()))
        )
        # interpolation function returning correction factor for friction factor
        # if tubes are not in a square arrangement
        self._cf_interp = RectBivariateSpline(
            x=np.array(list(InlineTubeBankFrictionData.cf_ax_dict.keys())),
            y=InlineTubeBankFrictionData.x_ax,
            z=np.array(list(InlineTubeBankFrictionData.cf_ax_dict.values()))
        )
        # correction factors for Nusselt number if N_l < 20
        self._cf_arr = [
            0.6768, 0.8089, 0.8687, 0.9054, 0.9303, 0.9465,
            0.9569, 0.9647, 0.9712, 0.9766, 0.9811, 0.9847,
            0.9877, 0.99, 0.992, 0.9937, 0.9953, 0.9969,
            0.9986
        ]

    @property
    def _u_max(self) -> Quantity:
        if self._u_inf is not None:
            u_max = self._u_inf * self.s_t / (self.s_t - self.D_o)
            return u_max.to('m / s')
        else:
            raise ValueError('no approach velocity set yet.')

    def _friction_factor(self) -> float:
        s_l_norm = (self.s_l / self.D_o).m
        Re = self._reynolds_number()
        f = self._look_up_square_friction_factor(Re, s_l_norm)
        if self.s_t != self.s_l:
            chi = self._look_up_correction_friction_factor()
            f *= chi
        return f

    def _look_up_square_friction_factor(self, Re: float, s_l_norm: float) -> float:
        f = self._f_interp(s_l_norm, Re)
        return f

    def _look_up_correction_friction_factor(self) -> float:
        Re = self._reynolds_number()
        x = (self.s_t / self.D_o - 1) * (self.s_l / self.D_o - 1)
        cf = self._cf_interp(Re, x)
        return cf

    def _tube_row_correction(self) -> float:
        if self.N_l <= 19:
            N_l = int(self.N_l)
            if N_l < 1:
                N_l = 1
            cf = self._cf_arr[N_l - 1]
            return cf
        else:
            return 1.0

    def avg_heat_trf_coeff(self) -> Quantity:
        Re = self._reynolds_number()
        if Re < 100:
            C, m = 0.9, 0.4
        elif Re < 1000:
            C, m = 0.52, 0.05  # shouldn't m be 0.5 instead?
        elif Re < 2e5:
            C, m = 0.27, 0.63
        else:
            C, m = 0.033, 0.8
        Nu = self._nusselt_number(C, m)
        k = self.fluid.k.to('W / (m * K)').m
        D_o = self.D_o.to('m').m
        h_avg = Nu * k / D_o
        return Q_(h_avg, 'W / (m ** 2 * K)')


class StaggeredTubeBank(TubeBank):
    """Class that represents a tube bank with a staggered tube arrangement
    (each row of tubes is offset by half of the transverse pitch relative
    to its adjacent rows).
    """
    def __init__(
        self,
        N_l: int,
        D_o: Quantity,
        s_l: Quantity,
        s_t: Quantity,
        fluid: FluidState
    ) -> None:
        """Creates an `StaggeredTubeBank` instance.

        Parameters
        ----------
        N_l: int
            Number of tube rows in the flow direction (longitudinal direction).
        D_o: Quantity
            Outer tube parameter
        s_l: Quantity
            Longitudinal pitch, i.e. the tube-to-tube spacing in the flow
            direction.
        s_t: Quantity
            Transverse pitch, i.e. the tube-to-tube spacing perpendicular to
            the flow direction.
        fluid: FluidState
            Fluid that flows across the bank of tubes.
        """
        super().__init__(N_l, D_o, s_l, s_t, fluid)
        # center-to-center distance of tubes in adjacent rows
        self.s_d = (self.s_l ** 2 + (self.s_t / 2) ** 2) ** 0.5
        # interpolation function returning friction factor for inline square
        self._f_interp = RectBivariateSpline(
            x=np.array(list(StaggeredTubeBankFrictionData.f_ax_dict.keys())),
            y=StaggeredTubeBankFrictionData.Re_ax,
            z=np.array(list(StaggeredTubeBankFrictionData.f_ax_dict.values()))
        )
        # interpolation function returning correction factor for friction factor
        # if tubes are not in a square arrangement
        self._cf_interp = RectBivariateSpline(
            x=np.array(list(StaggeredTubeBankFrictionData.cf_ax_dict.keys())),
            y=StaggeredTubeBankFrictionData.x_ax,
            z=np.array(list(StaggeredTubeBankFrictionData.cf_ax_dict.values()))
        )
        # correction factors for Nusselt number if N_l < 20
        self._cf_arr_low_Re = [
            0.8295, 0.8792, 0.9151, 0.9402, 0.957, 0.9677,
            0.9745, 0.9785, 0.9808, 0.9823, 0.9838, 0.9855,
            0.9873, 0.9891, 0.991, 0.9929, 0.9948, 0.9967,
            0.9987
        ]
        self._cf_arr_high_Re = [
            0.6273, 0.7689, 0.8473, 0.8942, 0.9254, 0.945,
            0.957, 0.9652, 0.9716, 0.9765, 0.9803, 0.9834,
            0.9862, 0.989, 0.9918, 0.9943, 0.9965, 0.998,
            0.9986
        ]

    @property
    def _u_max(self) -> Quantity:
        if self._u_inf is not None:
            c = self.s_d - self.D_o < (self.s_t - self.D_o) / 2
            if c:
                u_max = self._u_inf * self.s_t / (2 * (self.s_d - self.D_o))
            else:
                u_max = self._u_inf * self.s_t / (self.s_t - self.D_o)
            return u_max.to('m / s')
        else:
            raise ValueError('no approach velocity set yet.')

    def _friction_factor(self) -> float:
        s_t_norm = (self.s_t / self.D_o).m
        Re = self._reynolds_number()
        f = self._look_up_square_friction_factor(Re, s_t_norm)
        if self.s_t != self.s_d:
            chi = self._look_up_correction_friction_factor()
            f *= chi
        return f

    def _look_up_square_friction_factor(self, Re: float, s_t_norm: float) -> float:
        f = self._f_interp(s_t_norm, Re)
        return f

    def _look_up_correction_friction_factor(self) -> float:
        Re = self._reynolds_number()
        x = self.s_t / self.s_l
        cf = self._cf_interp(Re, x)
        return cf

    def _tube_row_correction(self) -> float:
        if self.N_l <= 19:
            Re = self._reynolds_number()
            N_l = int(self.N_l)
            if N_l < 1:
                N_l = 1
            if Re < 1000:
                cf = self._cf_arr_low_Re[N_l - 1]
            else:
                cf = self._cf_arr_high_Re[N_l - 1]
            return cf
        else:
            return 1.0

    def avg_heat_trf_coeff(self) -> Quantity:
        Re = self._reynolds_number()
        if Re < 500:
            C, m = 1.04, 0.4
        elif Re < 1000:
            C, m = 0.71, 0.5
        elif Re < 2e5:
            C, m = 0.35, 0.60
            f = (self.s_t / self.s_l) ** 0.2
            C *= f
        else:
            C, m = 0.031, 0.8
            f = (self.s_t / self.s_l) ** 0.2
            C *= f
        Nu = self._nusselt_number(C, m)
        k = self.fluid.k.to('W / (m * K)').m
        D_o = self.D_o.to('m').m
        h_avg = Nu * k / D_o
        return Q_(h_avg, 'W / (m ** 2 * K)')
