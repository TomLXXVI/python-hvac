"""Correlations for the drag coefficient and the average Nusselt number in case
of forced convective cross flow of a gas past non-circular cylinders.

The correlations for average Nusselt number were taken from:
Incropera, F. P., DeWitt, D. P., Bergman, T. L., & Lavine, A. S. (2019).
FUNDAMENTALS OF HEAT AND MASS TRANSFER.
Wiley.

The drag coefficients were taken from:
White, F. M. (2016)
FLUID MECHANICS.
New York: McGraw-Hill Education.
"""
import warnings
from abc import ABC
from hvac import Quantity
from hvac.fluids import FluidState
from .general import reynolds_number, prandtl_number, drag_force


class NonCircularCylinder(ABC):
    C: float = None
    m: float = None
    C_D: float = None
    Re_D_rng: tuple[float, float] = None

    def __init__(self, D: Quantity, L: Quantity, fluid: FluidState) -> None:
        """Creates a `NonCircularCylinder` instance.

        Parameters
        ----------
        D: Quantity
            Characteristic dimension of the cylinder, i.e. the dimension
            perpendicular to the direction of the flow.
        L: Quantity
            Length of the cylinder.
        fluid: FluidState
            The fluid moving along the surface of the cylinder. A `FluidState`
            object groups the thermophysical properties of a fluid in a given
            state.
        """
        self.D = D.to('m')
        self.L = L.to('m')
        self.fluid = fluid
        self._u_inf = None

    @property
    def u_inf(self) -> Quantity:
        """Returns the free stream velocity of the flowing fluid with respect
        to the cylinder.
        """
        return self._u_inf

    @u_inf.setter
    def u_inf(self, v: Quantity) -> None:
        """Sets the free stream velocity of the flowing fluid with respect to the
        cylinder.
        """
        self._u_inf = v.to('m / s')

    @property
    def Re(self) -> float:
        Re_D = reynolds_number(
            self.fluid.rho, self.fluid.mu,
            self.u_inf, self.D
        )
        return Re_D

    @property
    def Pr(self) -> float:
        Pr = prandtl_number(
            self.fluid.rho, self.fluid.mu,
            self.fluid.k, self.fluid.cp
        )
        return Pr

    def drag_force(self) -> Quantity:
        """Returns the drag force exerted by the moving fluid on the cylinder."""
        C_D = self.drag_coefficient(self.Re)
        F_D = drag_force(C_D, self.fluid.rho, self.u_inf, self.L * self.D)
        return F_D

    def avg_heat_trf_coeff(self) -> Quantity:
        """Returns the average heat transfer coefficient between the cylinder and
        the moving fluid.
        """
        Nu_D_avg = self.average_nusselt_number(self.Re, self.Pr)
        h_avg = Nu_D_avg * self.fluid.k / self.D
        return h_avg.to('W / (m ** 2 * K)')

    def drag_coefficient(self, Re_D: float) -> float:
        """Returns the drag coefficient C_D valid for Re_D_o >= 1.e4."""
        self._check_drag_coeff_validity(Re_D)
        return self.C_D

    def average_nusselt_number(self, Re_D: float, Pr: float) -> float:
        """Calculates the average Nusselt number acc. to the correlation of
        Jakob (1949).

        Parameters
        ----------
        Re_D:
            Reynolds number
        Pr:
            Prandtl number

        Returns
        -------
        Nu_D_avg: float
            Average Nusselt number

        Notes
        -----
        The correlation is valid for Pr >= 0.7 and within a certain range of
        Re_D_o which depends on the cross-sectional shape of the cylinder:
        - rhombus : 6000 - 60,000
        - square : 5000 - 60,000
        - hexagon, flat side towards fluid stream : 5200 - 105,000
        - hexagon, edge towards fluid stream : 4500 - 90,700
        - thin plate : 10,000 - 50,000
        """
        self._check_nusselt_validity(Re_D, Pr)
        Nu_D_avg = self.C * (Re_D ** self.m) * (Pr ** (1 / 3))
        return Nu_D_avg

    def _check_nusselt_validity(self, Re_D: float, Pr: float) -> None:
        if not (self.Re_D_rng[0] <= Re_D <= self.Re_D_rng[1]):
            warnings.warn(
                f"Reynolds number {round(Re_D, 3)} out of range "
                f"[{self.Re_D_rng[0]}, {self.Re_D_rng[1]}]"
            )
        if not (Pr >= 0.7):
            warnings.warn(
                f"Prandtl number {round(Pr, 3)} is less than 0.7"
            )

    @staticmethod
    def _check_drag_coeff_validity(Re_D: float) -> None:
        if not Re_D >= 1.e4:
            warnings.warn(
                f"Reynolds number {round(Re_D)} is less than 10,000"
            )


class RhombusCylinder(NonCircularCylinder):
    C = 0.304
    m = 0.59
    C_D = 1.6
    Re_D_rng = (6000, 60_000)


class SquareCylinder(NonCircularCylinder):
    C = 0.158
    m = 0.66
    C_D = 2.1
    Re_D_rng = (5000, 60_000)


class HexagonCylinder:

    class FlatSideToStream(NonCircularCylinder):
        C = (0.164, 0.039)
        m = (0.638, 0.78)
        C_D = 0.7
        Re_D_rng = (5200, 105_000)

        def average_nusselt_number(self, Re_D: float, Pr: float) -> float:
            self._check_nusselt_validity(Re_D, Pr)
            if 5200 <= Re_D < 20_400:
                C, m = self.C[0], self.m[0]
            else:
                C, m = self.C[1], self.m[1]
            Nu_D_avg = C * (Re_D ** m) * (Pr ** (1 / 3))
            return Nu_D_avg

    class EdgeToStream(NonCircularCylinder):
        C = 0.150
        m = 0.638
        C_D = 1.0
        Re_D_rng = (4500, 90_700)


class ThinPlate(NonCircularCylinder):
    C = 0.667
    m = 0.5
    C_D = 2.0
    Re_D_rng = (10_000, 50_000)
