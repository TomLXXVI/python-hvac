"""
Correlations for local and average friction factor and Nusselt number in case
of external forced convective fluid flow parallel along a flat plate.

The correlations were taken from:
Nellis G. F., & Klein S. A. (2021)
INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.
"""
import math
from scipy.optimize import root_scalar
from scipy.integrate import quad
from hvac import Quantity
from hvac.fluids import FluidState
from .general import (
    reynolds_number,
    prandtl_number,
    Re_crit,
    get_flow_condition,
    shear_stress
)

Q_ = Quantity


class LaminarFlow:
    """Correlations for laminar flow over a flat plate."""

    @staticmethod
    def local_friction_factor(Re_x: float) -> float:
        """Calculates the local friction factor at position x along the
        characteristic length of the plate.

        Parameters
        ----------
        Re_x: float
            Reynolds number at position x along the characteristic length of
            the plate.

        Returns
        -------
        Cf: float
            Local friction factor.
        """
        Cf = 0.664 / math.sqrt(Re_x)
        return Cf

    @staticmethod
    def average_friction_factor(Re_x: float) -> float:
        """Calculates the average friction factor between the leading edge
        of the plate and position x along the characteristic length of the
        plate.

        Parameters
        ----------
        Re_x: float
            Reynolds number at position x along the characteristic length of
            the plate.

        Returns
        -------
        Cf_avg: float
            The average friction factor.
        """

        def f(Re_x_: float) -> float:
            Cf_lam = LaminarFlow.local_friction_factor(Re_x_)
            return Cf_lam

        I = quad(f, 0, Re_x)[0]
        Cf_avg = 1 / Re_x * I
        return Cf_avg

    class UniformTemperature:
        """Applies to a plate with constant, uniform temperature."""
        @staticmethod
        def local_nusselt_number(Re_x: float, Pr: float) -> float:
            """Calculates the local Nusselt number at position x along the
            characteristic length of the plate.

            Parameters
            ----------
            Re_x: float
                Reynolds number at position x along the characteristic length of
                the plate.
            Pr: float
                Prandtl number of the fluid.

            Returns
            -------
            Nu_x: float
                Local Nusselt number.
            """
            n = 0.3387 * Re_x ** (1 / 2) * Pr ** (1 / 3)
            d = ((1 + 0.0468 / Pr) ** (2 / 3)) ** (1 / 4)
            Nu_x = n / d
            return Nu_x

        @staticmethod
        def average_nusselt_number(Re_x: float, Pr: float) -> float:
            """Calculates the average Nusselt number between the leading edge
            of the plate and position x along the characteristic length of the
            plate.

            Parameters
            ----------
            Re_x: float
                Reynolds number at position x along the characteristic length of
                the plate.
            Pr: float
                Prandtl number of the fluid.

            Returns
            -------
            Nu_x_avg: float
                Average Nusselt number.
            """
            def f(Re_x_: float) -> float:
                LUT = LaminarFlow.UniformTemperature
                Nu_x = LUT.local_nusselt_number(Re_x_, Pr)
                return Nu_x / Re_x_

            Nu_x_avg = quad(f, 0, Re_x)[0]
            return Nu_x_avg
            # n = 0.6774 * Pr ** (1 / 3) * Re_x ** 0.5
            # d = (1 + (0.0468 / Pr) ** (2 / 3)) ** (1 / 4)
            # Nu_x_avg = n / d
            # return Nu_x_avg

    class UnheatedStartingLength:
        """Applies to a plate with an adiabatic section at its leading edge
        where the flow develops hydrodynamically before it begins to develop
        thermally.
        """
        @staticmethod
        def local_nusselt_number(
            x: Quantity,
            Re_x: float,
            Pr: float,
            L_uh: Quantity
        ) -> float:
            """Calculates the local Nusselt number at position x along the
            characteristic length of the plate.

            Parameters
            ----------
            x: Quantity
                Position x along the characteristic length of the plate where
                the local Nusselt number is to be determined.
            Re_x: float
                Reynolds number at position x along the characteristic length of
                the plate.
            Pr: float
                Prandtl number of the fluid.
            L_uh: Quantity
                Length of the unheated section starting from the leading edge
                in the direction of the flow.

            Returns
            -------
            Nu_x: float
                Local Nusselt number.
            """
            LUT = LaminarFlow.UniformTemperature
            Nu_x_uh = LUT.local_nusselt_number(Re_x, Pr)
            x = x.to('m').m
            L_uh = L_uh.to('m').m
            Nu_x = Nu_x_uh / ((1 - (L_uh / x) ** 0.75) ** (1 / 3))
            return Nu_x

        @staticmethod
        def average_nusselt_number(
            x: Quantity,
            Re_x: float,
            Pr: float,
            L_uh: Quantity
        ) -> float:
            """Calculates the average Nusselt number between the leading edge
            of the plate and position x along the characteristic length of the
            plate.

            Parameters
            ----------
            x: Quantity
                Position x along the characteristic length of the plate where
                the average Nusselt number is to be determined starting from
                the leading edge of the plate.
            Re_x: float
                Reynolds number at position x along the characteristic length of
                the plate.
            Pr: float
                Prandtl number of the fluid.
            L_uh: Quantity
                Length of the unheated section starting from the leading edge
                in the direction of the flow.

            Returns
            -------
            Nu_x_avg: float
                Average Nusselt number.
            """
            LUT = LaminarFlow.UniformTemperature
            Nu_x_avg_uh = LUT.average_nusselt_number(Re_x, Pr)
            x = x.to('m').m
            L_uh = L_uh.to('m').m
            Nu_x_avg = Nu_x_avg_uh
            p = 2
            Nu_x_avg *= x / (x - L_uh)
            Nu_x_avg *= (1 - (L_uh / x) ** ((p + 1) / (p + 2))) ** (p / (p + 1))
            return Nu_x_avg

    class UniformHeatFlux:
        """Applies to a plate heated (or cooled) with constant, uniform
        heat flux.
        """
        @staticmethod
        def local_nusselt_number(Re_x: float, Pr: float) -> float:
            """Calculates the local Nusselt number at position x along the
            characteristic length of the plate.

            Parameters
            ----------
            Re_x: float
                Reynolds number at position x along the characteristic length of
                the plate.
            Pr: float
                Prandtl number of the fluid.

            Returns
            -------
            Nu_x: float
                Local Nusselt number.

            Notes
            -----
            Correlation from Kays et al. (2005) valid for Prandtl numbers
            greater than 0.6.
            """
            Nu_x = 0.453 * Re_x ** (1 / 2) * Pr * (1 / 3)
            return Nu_x

        @staticmethod
        def average_nusselt_number(Re_x: float, Pr: float) -> float:
            """Calculates the average Nusselt number between the leading edge
            of the plate and position x along the characteristic length of the
            plate.

            Parameters
            ----------
            Re_x: float
                Reynolds number at position x along the characteristic length of
                the plate.
            Pr: float
                Prandtl number of the fluid.

            Returns
            -------
            Nu_x_avg: float
                Average Nusselt number.
            """
            Nu_x_avg = 0.680 * Re_x ** (1 / 2) * Pr ** (1 / 3)
            return Nu_x_avg


class TurbulentFlow:
    """Correlations for turbulent flow over a flat plate."""

    @staticmethod
    def _local_friction_factor_smooth(Re_x: float) -> float:
        """Calculates the local friction factor at position x along the
        characteristic length of the plate in case of a smooth plate
        (relative roughness er_x = e / x < 1.e-6 with e the absolute surface
        roughness of the plate).

        Parameters
        ----------
        Re_x: float
            Reynolds number at position x along the characteristic length of
            the plate.

        Returns
        -------
        Cf: float
            Local friction factor.
        """
        # when
        Cf = 0.455 / (math.log(0.06 * Re_x) ** 2)
        return Cf

    @staticmethod
    def _local_friction_factor_rough(Re_x: float, er_x: float) -> float:
        """Calculates the local friction factor at position x along the
        characteristic length of the plate in case of a rough plate
        (relative roughness er_x = e / x >= 1.e-6 with e the absolute surface
        roughness of the plate).

        Parameters
        ----------
        Re_x: float
            Reynolds number at position x along the characteristic length of
            the plate.
        er_x: float
            Relative roughness of the plate at position x along the
            characteristic length of the plate.

        Returns
        -------
        Cf: float
            Local friction factor.
        """
        kappa = 0.41  # Von Kármán constant

        def eq(Cf: float) -> float:
            Z = kappa * math.sqrt(2 / Cf)
            e_plus = Re_x * er_x / math.sqrt(2 / Cf)
            Re_x_new = 1.73 * (1 + 0.3 * e_plus) * math.exp(Z)
            Re_x_new *= (
                Z ** 2 - 4 * Z + 6
                - 0.3 * e_plus * (Z - 1) / (1 + 0.3 * e_plus)
            )
            return Re_x_new - Re_x

        Cf_smooth = TurbulentFlow._local_friction_factor_smooth(Re_x)
        Cf_max = 0.012
        # see Introduction to Engineering Heat Transfer (Nellis, G.F.), fig. 8.2
        Cf = root_scalar(eq, bracket=[Cf_smooth, Cf_max]).root
        return Cf

    @staticmethod
    def local_friction_factor(Re_x: float, er_x: float = 0.0):
        """Calculates the local friction factor at position x along the
        characteristic length of the plate. The characteristic length of a flat
        plate is the dimension of the plate in the direction of the flow.

        Parameters
        ----------
        Re_x: float
            Reynolds number at position x along the characteristic length of
            the plate.
        er_x: float
            Relative roughness of the plate at position x along the
            characteristic length of the plate. It is the ratio of the absolute
            surface roughness e to the distance x measured from the leading edge
            of the plate.

        Returns
        -------
        Cf: float
            Local friction factor.
        """
        if er_x < 1.e-6:
            return TurbulentFlow._local_friction_factor_smooth(Re_x)
        else:
            return TurbulentFlow._local_friction_factor_rough(Re_x, er_x)

    @staticmethod
    def average_friction_factor(
        x: Quantity,
        e: Quantity,
        u_inf: Quantity,
        rho: Quantity,
        mu: Quantity
    ) -> float:
        """Calculates the average friction factor between the leading edge
        of the plate and position x along the characteristic length of the
        plate.

        Parameters
        ----------
        x: Quantity
            Position x along the characteristic length of the plate where the
            average friction factor is to be determined.
        e: Quantity
            Absolute surface roughness of the plate.
        u_inf: Quantity
            Free stream velocity of flowing fluid along the characteristic
            length of the plate.
        rho: Quantity
            Mass density of fluid.
        mu: Quantity
            Absolute or dynamic viscosity of the fluid.

        Returns
        -------
        Cf_avg: float
            The average friction factor.
        """
        Re_x = reynolds_number(rho, mu, u_inf, x)
        e = e.to('m').m
        rho = rho.to('kg / m ** 3').m
        mu = mu.to('Pa * s').m
        u_inf = u_inf.to('m / s').m

        def f1(Re_x_: float) -> float:
            Cf_lam = LaminarFlow.local_friction_factor(Re_x_)
            return Cf_lam

        def f2(Re_x_: float) -> float:
            x_ = Re_x_ * mu / (rho * u_inf)
            er_x_ = e / x_
            Cf_tur = TurbulentFlow.local_friction_factor(Re_x_, er_x_)
            return Cf_tur

        I1 = quad(f1, 0, Re_crit)[0]
        I2 = quad(f2, Re_crit, Re_x)[0]
        Cf_avg = (1 / Re_x) * (I1 + I2)
        return Cf_avg

    class UniformTemperature:
        """Plate at constant, uniform temperature."""
        @staticmethod
        def _local_nusselt_number_smooth(Re_x: float, Pr: float) -> float:
            """Calculates the local Nusselt number at position x along the
            characteristic length of the plate in case of a smooth plate.
            (relative roughness er_x = e / x < 1.e-6 with e the absolute surface
            roughness of the plate).

            Parameters
            ----------
            Re_x: float
                Reynolds number at position x along the characteristic length of
                the plate.
            Pr: float
                Prandtl number.

            Returns
            -------
            Nu_x: float
                Local Nusselt number

            Notes
            -----
            The correlation is based on the modified Reynolds analogy (Colburn
            analogy) and valid for Re_crit < Re_x < 1.e8 and 0.5 < Pr < 60.
            """
            Nu_x = 0.0135 * Re_x ** (6 / 7) * Pr ** (1 / 3)
            return Nu_x

        # noinspection PyUnusedLocal
        @staticmethod
        def _local_nusselt_number_rough(Re_x: float, Pr: float, er_x: float) -> float:
            # not implemented: acc. to Nellis G.F., Introduction to Engineering
            # Heat Transfer, § 8.2.2. there is no generally acceptable correlation
            # for the Nusselt number associated with flow over a rough plate.
            return float('nan')

        @staticmethod
        def local_nusselt_number(Re_x: float, Pr: float, er_x: float = 0.0) -> float:
            """Calculates the local Nusselt number at position x along the
            characteristic length of the plate.

            Parameters
            ----------
            Re_x: float
                Reynolds number at position x.
            Pr: float
                Prandtl number of the fluid.
            er_x: float
                Relative roughness at position x.

            Returns
            -------
            Nu_x: float
                Local Nusselt number.
            """
            TUT = TurbulentFlow.UniformTemperature
            if er_x < 1.e-6:
                return TUT._local_nusselt_number_smooth(Re_x, Pr)
            else:
                return TUT._local_nusselt_number_rough(Re_x, Pr, er_x)

        @staticmethod
        def average_nusselt_number(
            x: Quantity,
            e: Quantity,
            Pr: float,
            u_inf: Quantity,
            rho: Quantity,
            mu: Quantity
        ) -> float:
            """Calculates the average Nusselt number between the leading edge
            of the plate and position x along the characteristic length of the
            plate.

            Parameters
            ----------
            x: Quantity
                Position x along the characteristic length of the plate.
            e: Quantity
                Absolute surface roughness of the plate.
            Pr: float
                Prandtl number of the fluid.
            u_inf: Quantity
                Free stream velocity of the fluid.
            rho: Quantity
                Mass density of fluid.
            mu: Quantity
                Absolute or dynamic viscosity of the fluid.

            Returns
            -------
            Nu_x_avg: float
                Average Nusselt number.
            """
            Re_x = reynolds_number(rho, mu, u_inf, x)
            e = e.to('m').m
            rho = rho.to('kg / m ** 3').m
            mu = mu.to('Pa * s').m
            u_inf = u_inf.to('m / s').m
            LUT = LaminarFlow.UniformTemperature
            TUT = TurbulentFlow.UniformTemperature

            def f1(Re_x_: float) -> float:
                Nu_x_lam = LUT.local_nusselt_number(Re_x_, Pr)
                return Nu_x_lam / Re_x_

            def f2(Re_x_: float) -> float:
                x_ = Re_x_ * mu / (rho * u_inf)
                er_x_ = e / x_
                Nu_x_tur = TUT.local_nusselt_number(Re_x_, Pr, er_x_)
                return Nu_x_tur / Re_x_

            I1 = quad(f1, 0, Re_crit)[0]
            I2 = quad(f2, Re_crit, Re_x)[0]
            Nu_x_avg = I1 + I2
            return Nu_x_avg

    class UnheatedStartingLength:
        """Applies to a plate with an adiabatic section at its leading edge
        where the flow develops hydrodynamically before it begins to develop
        thermally.
        """
        @staticmethod
        def local_nusselt_number(
            x: Quantity,
            Re_x: float,
            Pr: float,
            L_uh: Quantity,
            e: Quantity = Q_(0.0, 'mm')
        ) -> float:
            """Calculates the local Nusselt number at position x along the
            characteristic length of the plate.

            Parameters
            ----------
            x: Quantity
                Position x along the characteristic length of the plate where
                the local Nusselt number is to be determined.
            Re_x: float
                Reynolds number at position x along the characteristic length of
                the plate.
            Pr: float
                Prandtl number of the fluid.
            L_uh: Quantity
                Length of the unheated section starting from the leading edge
                in the direction of the flow.
            e: Quantity
                Absolute surface roughness of the plate.

            Returns
            -------
            Nu_x: float
                Local Nusselt number.
            """
            TUT = TurbulentFlow.UniformTemperature
            x = x.to('m').m
            L_uh = L_uh.to('m').m
            e = e.to('m').m
            er_x = e / x
            Nu_x_uh = TUT.local_nusselt_number(Re_x, Pr, er_x)
            Nu_x = Nu_x_uh / ((1 - (L_uh / x) ** 0.90) ** (1 / 9))
            return Nu_x

        @staticmethod
        def average_nusselt_number(
            x: Quantity,
            e: Quantity,
            L_uh: Quantity,
            Pr: float,
            u_inf: Quantity,
            rho: Quantity,
            mu: Quantity
        ) -> float:
            """Calculates the average Nusselt number between the leading edge
            of the plate and position x along the characteristic length of the
            plate.

            Parameters
            ----------
            x: Quantity
                Position x along the characteristic length of the plate.
            e: Quantity
                Absolute surface roughness of the plate.
            L_uh: Quantity
                Length of the unheated section starting from the leading edge
                in the direction of the flow.
            Pr: float
                Prandtl number of the fluid.
            u_inf: Quantity
                Free stream velocity of the fluid.
            rho: Quantity
                Mass density of fluid.
            mu: Quantity
                Absolute or dynamic viscosity of the fluid.

            Returns
            -------
            Nu_x_avg: float
                Average Nusselt number.
            """
            TUT = TurbulentFlow.UniformTemperature
            Nu_x_avg_uh = TUT.average_nusselt_number(x, e, Pr, u_inf, rho, mu)
            x = x.to('m').m
            L_uh = L_uh.to('m').m
            Nu_x_avg = Nu_x_avg_uh
            p = 8
            Nu_x_avg *= x / (x - L_uh)
            Nu_x_avg *= (1 - (L_uh / x) ** ((p + 1) / (p + 2))) ** (p / (p + 1))
            return Nu_x_avg

    class UniformHeatFlux:
        """Applies to a plate heated (or cooled) with constant, uniform
        heat flux.
        """
        @staticmethod
        def local_nusselt_number(Re_x: float, Pr: float) -> float:
            """Calculates the local Nusselt number at position x along the
            characteristic length of the plate.

            Parameters
            ----------
            Re_x: float
                Reynolds number at position x along the characteristic length of
                the plate.
            Pr: float
                Prandtl number of the fluid.

            Returns
            -------
            Nu_x: float
                Local Nusselt number.

            Notes
            -----
            Correlation from Kays et al. (2005) valid for Prandtl numbers
            greater than 0.6 and smaller than 60.
            """
            Nu_x = 0.0308 * Re_x ** 0.8 * Pr * (1 / 3)
            return Nu_x

        @staticmethod
        def average_nusselt_number(
            x: Quantity,
            e: Quantity,
            Pr: float,
            u_inf: Quantity,
            rho: Quantity,
            mu: Quantity
        ) -> float:
            """Calculates the average Nusselt number between the leading edge
            of the plate and position x along the characteristic length of the
            plate.

            Parameters
            ----------
            x: Quantity
                Position x along the characteristic length of the plate.
            e: Quantity
                Absolute surface roughness of the plate.
            Pr: float
                Prandtl number of the fluid.
            u_inf: Quantity
                Free stream velocity of the fluid.
            rho: Quantity
                Mass density of fluid.
            mu: Quantity
                Absolute or dynamic viscosity of the fluid.

            Returns
            -------
            Nu_x_avg: float
                Average Nusselt number.

            Notes
            -----
            According to Nellis, G.F. Introduction to Engineering Heat Transfer
            the average Nusselt number for uniform heat flux agrees closely with
            the uniform temperature result.
            """
            TUT = TurbulentFlow.UniformTemperature
            Nu_x_avg = TUT.average_nusselt_number(x, e, Pr, u_inf, rho, mu)
            return Nu_x_avg


class Plate:
    """Class for calculating the average shear stress exerted on a flat plate
    by a moving fluid and the average heat transfer coefficient between the
    plate and the moving fluid.
    """
    def __init__(
        self,
        L_char: Quantity,
        b: Quantity,
        fluid: FluidState | None,
        e: Quantity = Q_(0.0, 'mm'),
        L_uh: Quantity = Q_(0.0, 'm')
    ) -> None:
        """Creates a `Plate` instance with the given fixed parameters.

        Parameters
        ----------
        L_char: Quantity
            Characteristic length of the flat plate, i.e. the dimension of
            the plate in the direction of the flow.
        b: Quantity
            The other dimension of the flat plate.
        fluid: FluidState
            The fluid moving along the surface of the plate. A `FluidState`
            object groups the thermophysical properties of a fluid in a given
            state.
        e: Quantity
            Absolute surface roughness of the plate.
        L_uh: Quantity
            Unheated starting length of the plate in the direction of flow.
        """
        self.L_char = L_char.to('m')
        self.b = b.to('m')
        self.fluid = fluid
        self.e = e.to('m')
        self.L_uh = L_uh.to('m')
        self._u_inf = None

    @property
    def u_inf(self) -> Quantity:
        """Returns the free stream velocity of the flowing fluid with respect
        to the plate.
        """
        return self._u_inf

    @u_inf.setter
    def u_inf(self, v: Quantity) -> None:
        """Sets the free stream velocity of the flowing fluid with respect to the
        plate.
        """
        self._u_inf = v.to('m / s')

    @property
    def Re(self) -> float:
        Re_L = reynolds_number(
            self.fluid.rho, self.fluid.mu,
            self.u_inf, self.L_char
        )
        return Re_L

    @property
    def Pr(self) -> float:
        Pr = prandtl_number(
            self.fluid.rho, self.fluid.mu,
            self.fluid.k, self.fluid.cp
        )
        return Pr

    def avg_shear_stress(self) -> Quantity:
        """Returns the average shear stress exerted on the surface of the plate
        due to the moving fluid.
        """
        Re_L = self.Re
        flow_condition = get_flow_condition(Re_L)
        if flow_condition == 'laminar':
            Cf_avg = LaminarFlow.average_friction_factor(Re_L)
        else:
            Cf_avg = TurbulentFlow.average_friction_factor(
                x=self.L_char,
                e=self.e,
                u_inf=self._u_inf,
                rho=self.fluid.rho,
                mu=self.fluid.mu
            )
        tau_s = shear_stress(Cf_avg, self.fluid.rho, self._u_inf)
        return tau_s

    def avg_heat_trf_coeff(
        self,
        therm_bound_cond: str = 'ust'
    ) -> Quantity:
        """Returns the average heat transfer coefficient between the plate and
        the moving fluid.

        Parameters
        ----------
        therm_bound_cond: str, {'ust', 'uhs', 'uhf'}, default 'ust'
            Thermal boundary condition.
            - 'ust': plate with uniform surface temperature
            - 'uhs': plate with unheated starting length
            - 'uhf': plate heated/cooled with uniform heat flux
        """
        Re_L = self.Re
        Pr = self.Pr
        flow_condition = get_flow_condition(Re_L)
        L_UST = LaminarFlow.UniformTemperature
        L_UHS = LaminarFlow.UnheatedStartingLength
        L_UHF = LaminarFlow.UniformHeatFlux
        T_UST = TurbulentFlow.UniformTemperature
        T_UHS = TurbulentFlow.UnheatedStartingLength
        T_UHF = TurbulentFlow.UniformHeatFlux
        condition = flow_condition + '&' + therm_bound_cond
        Nu_avg = None
        match condition:
            case 'laminar&ust':
                Nu_avg = L_UST.average_nusselt_number(Re_L, Pr)
            case 'laminar&uhs':
                Nu_avg = L_UHS.average_nusselt_number(
                    self.L_char, Re_L,
                    Pr, self.L_uh
                )
            case 'laminar&uhf':
                Nu_avg = L_UHF.average_nusselt_number(Re_L, Pr)
            case 'turbulent&ust':
                Nu_avg = T_UST.average_nusselt_number(
                    self.L_char, self.e, Pr, self._u_inf,
                    self.fluid.rho, self.fluid.mu
                )
            case 'turbulent&uhs':
                Nu_avg = T_UHS.average_nusselt_number(
                    self.L_char, self.e, self.L_uh, Pr, self._u_inf,
                    self.fluid.rho, self.fluid.mu
                )
            case 'turbulent&uhf':
                Nu_avg = T_UHF.average_nusselt_number(
                    self.L_char, self.e, Pr, self._u_inf,
                    self.fluid.rho, self.fluid.mu
                )
        h_avg = Nu_avg * self.fluid.k / self.L_char
        return h_avg.to('W / (m ** 2 * K)')
