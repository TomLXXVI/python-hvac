"""
In this module some correlations are implemented that predict the behavior of
nucleate boiling and the limit of this mode which is the critical heat flux.

The correlations were taken from:
Nellis G. F., & Klein S. A. (2021)
INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.
"""
from math import pi as PI
from abc import ABC, abstractmethod
import warnings
from hvac import Quantity
from hvac.fluids import Fluid
from ..free_convection.general import prandtl_number


Q_ = Quantity


class Geometry(ABC):

    def __init__(
        self,
        L_char: float,
        L_nb: float,
        A_s: float | None = None
    ) -> None:
        self.L_tilde = L_char / L_nb
        self.L_nb = L_nb
        self.A_s = A_s

    @property
    @abstractmethod
    def C_crit(self) -> float:
        ...


class FlatPlate(Geometry):

    @property
    def C_crit(self) -> float:
        if self.L_tilde < 27.0:
            if self.L_tilde < 9.0:
                warnings.warn(f"L_char smaller than 9.0")
            elif 20.0 < self.L_tilde < 27.0:
                warnings.warn(f"L_char between 20.0 and 27.0")
            return 9 * PI * self.L_nb ** 2 / (5 * self.A_s)
        else:
            return 0.15


class HorizontalCylinder(Geometry):

    @property
    def C_crit(self) -> float:
        if self.L_tilde < 1.2:
            if self.L_tilde < 0.15:
                warnings.warn(f"L_char smaller than 0.15")
            return 0.12 / (self.L_tilde ** 0.25)
        else:
            return 0.12


class Sphere(Geometry):

    @property
    def C_crit(self) -> float:
        if self.L_tilde < 4.26:
            if self.L_tilde < 0.15:
                warnings.warn(f"L_char smaller than 0.15")
            return 0.227 / (self.L_tilde ** 0.5)
        else:
            return 0.11


class NucleateBoiling:
    geometries: dict[str, type[Geometry]] = {
        'plate': FlatPlate,
        'cylinder': HorizontalCylinder,
        'sphere': Sphere
    }
    g: Quantity = Q_(9.81, 'm / s ** 2')

    def __init__(
        self,
        fluid: Fluid,
        P_sat: Quantity,
        geometry: str | None = None,
        L_char: Quantity | None = None,
        A_s: Quantity | None = None
    ) -> None:
        """
        Creates a `NucleateBoiling` instance.

        Parameters
        ----------
        fluid: Fluid
            The fluid being considered.
        P_sat: PlainQuantity
            The pressure at which the fluid is boiling.
        geometry: str, {'plate', 'cylinder', 'sphere'}, optional
            Geometry of the heated surface in contact with the boiling fluid.
            If not specified, a large, finite body is assumed.
        L_char: PlainQuantity, optional
            Characteristic length of the heated surface in contact with the
            boiling fluid. This parameter must be set together with parameter
            `geometry`. In case `geometry` is 'plate', the characteristic length
            is the width or diameter of the plate; for a 'cylinder' or 'sphere'
            it is the radius of the cylinder or sphere.
        A_s: PlainQuantity, optional
            Surface area. This parameter is only needed if `geometry` is 'plate'.
        """
        self.fluid = fluid
        self.P_sat = P_sat
        if geometry is not None:
            self.geometry = self._create_geometry(geometry, L_char, A_s)
        else:
            self.geometry = None

    def _create_geometry(
        self,
        geometry: str,
        L_char: Quantity,
        A_s: Quantity | None = None
    ) -> Geometry:
        L_char = L_char.to('m').m
        fluid_liq_sat = self.fluid(P=self.P_sat, x=Q_(0, 'frac'))
        fluid_vap_sat = self.fluid(P=self.P_sat, x=Q_(1, 'frac'))
        sigma = fluid_liq_sat.sigma.to('N / m').m
        rho_liq_sat = fluid_liq_sat.rho.to('kg / m ** 3').m
        rho_vap_sat = fluid_vap_sat.rho.to('kg / m ** 3').m
        g = self.g.to('m / s ** 2').m
        L_nb = (sigma / (g * (rho_liq_sat - rho_vap_sat))) ** 0.5
        _Geometry = self.geometries[geometry]
        A_s = A_s.to('m ** 2').m if A_s is not None else None
        geometry = _Geometry(L_char, L_nb, A_s)
        return geometry

    def heat_flux(
        self,
        dT_e: Quantity,
        n: float = 1.0,
        C_nb: float = 0.013
    ) -> Quantity:
        """
        Get the heat flux in the nucleate boiling region of the given fluid as
        function of the excess temperature difference acc. to the correlation
        of Rohsenow (1952).

        Parameters
        ----------
        dT_e: PlainQuantity
            Excess temperature difference (i.e. the difference between the
            surface temperature and the saturation temperature of the fluid).
        n: float, {1.0, 1.7}
            Exponent, 1.0 for water and 1.7 for other fluids.
        C_nb: float, default 0.013
            Empirical constant related to the surface-fluid combination.
            See Nellis, INTRODUCTION TO ENGINEERING HEAT TRANSFER, table 11.1.

        Returns
        -------
        q_nb: PlainQuantity
            Heat flux in the nucleate boiling region.
        """
        dT_e = dT_e.to('K').m
        fluid_liq_sat = self.fluid(P=self.P_sat, x=Q_(0, 'frac'))
        Pr_liq_sat = prandtl_number(
            fluid_liq_sat.rho, fluid_liq_sat.mu,
            fluid_liq_sat.k, fluid_liq_sat.cp
        )
        mu_liq_sat = fluid_liq_sat.mu.to('Pa * s').m
        rho_liq_sat = fluid_liq_sat.rho.to('kg / m ** 3').m
        cp_liq_sat = fluid_liq_sat.cp.to('J / (kg * K)').m
        h_liq_sat = fluid_liq_sat.h.to('J / kg').m
        sigma = fluid_liq_sat.sigma.to('N / m').m

        fluid_vap_sat = self.fluid(P=self.P_sat, x=Q_(1, 'frac'))
        rho_vap_sat = fluid_vap_sat.rho.to('kg / m ** 3').m
        h_vap_sat = fluid_vap_sat.h.to('J / kg').m

        h_fg = h_vap_sat - h_liq_sat
        g = self.g.to('m / s ** 2').m
        q_nb = (
            mu_liq_sat * h_fg * (g * (rho_liq_sat - rho_vap_sat) / sigma) ** 0.5
            * ((cp_liq_sat * dT_e) / (C_nb * h_fg * Pr_liq_sat ** n)) ** 3
        )
        return Q_(q_nb, 'W / m ** 2')

    def heat_trf_coeff(
        self,
        dT_e: Quantity,
        n: float = 1.0,
        C_nb: float = 0.013
    ):
        """
        Get the heat transfer coefficient for nucleate boiling at the given
        excess temperature difference.

        Parameters
        ----------
        dT_e: PlainQuantity
            Excess temperature difference (i.e. the difference between the
            surface temperature and the saturation temperature of the fluid).
        n: float, {1.0, 1.7}
            Exponent, 1.0 for water and 1.7 for other fluids.
        C_nb: float, default 0.013
            Empirical constant related to the surface-fluid combination.
            See Nellis, INTRODUCTION TO ENGINEERING HEAT TRANSFER table 11.1.

        Returns
        -------
        h_nb: PlainQuantity
            Heat transfer coefficient for nucleate boiling.
        """
        dT_e = dT_e.to('K')
        q_nb = self.heat_flux(dT_e, n, C_nb)
        h_nb = q_nb / dT_e
        return h_nb.to('W / (m ** 2 * K)')

    def critical_heat_flux(self) -> Quantity:
        """Returns the heat flux at the burnout point.

        The burnout point is reached when vapor, which has a much lower
        conductivity than liquid, is produced at a rate that is so high that it
        hinders liquid to re-wet the surface. The result is a drastic and sudden
        increase in the excess temperature difference that may cause damage or
        even melt the material (i.e. boiling crisis). In most devices, it is
        important that the heat flux be kept below the critical heat flux to
        avoid the boiling crisis and the device operates safely in the
        nucleate boiling region.
        """
        fluid_liq_sat = self.fluid(P=self.P_sat, x=Q_(0, 'frac'))
        rho_liq_sat = fluid_liq_sat.rho.to('kg / m ** 3').m
        h_liq_sat = fluid_liq_sat.h.to('J / kg').m
        sigma = fluid_liq_sat.sigma.to('N / m').m

        fluid_vap_sat = self.fluid(P=self.P_sat, x=Q_(1, 'frac'))
        rho_vap_sat = fluid_vap_sat.rho.to('kg / m ** 3').m
        h_vap_sat = fluid_vap_sat.h.to('J / kg').m

        C_crit = self.geometry.C_crit if self.geometry is not None else 0.12
        h_fg = h_vap_sat - h_liq_sat
        g = self.g.to('m / s ** 2').m
        q_crit = (
            C_crit * h_fg * rho_vap_sat
            * (
                (sigma * g * (rho_liq_sat - rho_vap_sat))
                / (rho_vap_sat ** 2)
            ) ** 0.25
        )
        return Q_(q_crit, 'W / m ** 2')
