import numpy as np
from scipy.special import iv, kv
from abc import ABC, abstractmethod
from hvac import Quantity


Q_ = Quantity
bessel_i = iv   # modified Bessel function of the first kind of real order
bessel_k = kv   # modified Bessel function of the second kind of real order


class Fin(ABC):
    """General base class to represent a fin."""

    def __init__(self) -> None:
        """Creates a `Fin` instance."""
        self._h_avg: float = 0.0
        self._T_b: float = 0.0
        self._T_fluid: float = 0.0

    @property
    def h_avg(self) -> Quantity:
        """Gets the average heat transfer coefficient between the fin and the
        surrounding fluid.
        """
        return Q_(self._h_avg, 'W / (m ** 2 * K)')

    @h_avg.setter
    def h_avg(self, v: Quantity):
        """Sets the average heat transfer coefficient between the fin and the
        surrounding fluid.
        """
        self._h_avg = v.to('W / (m ** 2 * K)').m

    @property
    def T_b(self) -> Quantity:
        """Gets the base temperature of the fin."""
        return Q_(self._T_b, 'K')

    @T_b.setter
    def T_b(self, v: Quantity) -> None:
        """Sets the base temperature of the fin."""
        self._T_b = v.to('K').m

    @property
    def T_fluid(self) -> Quantity:
        """Gets the fluid temperature (at some distance from the fin)."""
        return Q_(self._T_fluid, 'K')

    @T_fluid.setter
    def T_fluid(self, v: Quantity) -> None:
        """Sets the fluid temperature (at some distance from the fin)."""
        self._T_fluid = v.to('K').m

    @property
    @abstractmethod
    def _m(self) -> float:
        ...

    @property
    @abstractmethod
    def surface_area(self) -> Quantity:
        """Gets the surface area of the fin involved in the heat transfer between
        the fin and the surrounding fluid (fin tip included).
        """
        ...

    @property
    @abstractmethod
    def efficiency(self) -> float:
        """Gets the fin efficiency, i.e. the ratio of the actual heat transfer
        rate between the fin and the surrounding fluid to the heat transfer rate
        that would occur if the fin is at a uniform temperature equal to the
        base temperature."""
        ...

    @property
    def resistance(self) -> Quantity:
        """Gets the fin thermal resistance."""
        eta_fin = self.efficiency
        A_f = self.surface_area.m
        R_f = 1 / (eta_fin * self._h_avg * A_f)
        return Q_(R_f, 'K / W')

    @property
    def heat_transfer(self) -> Quantity:
        """Gets the heat transfer between the fin and the surrounding
        fluid.
        """
        R = self.resistance.m
        Q = (self._T_b - self._T_fluid) / R
        return Q_(Q, 'W')

    @property
    @abstractmethod
    def base_area(self) -> Quantity:
        """Gets the cross-sectional area of the fin base."""
        ...

    @property
    def effectiveness(self) -> float:
        """Gets the fin effectiveness, i.e. the ratio of the actual heat transfer
        rate to the heat transfer rate that would occur without fin."""
        A_f = self.surface_area.m
        A_b = self.base_area.m
        Q_id_fin = A_f * self._h_avg * (self._T_b - self._T_fluid)
        Q_no_fin = A_b * self._h_avg * (self._T_b - self._T_fluid)
        Q_fin = self.efficiency * Q_id_fin
        eps_fin = Q_fin / Q_no_fin
        return eps_fin


class StraightFin(Fin):
    """General base class to represent straight fins."""

    def __init__(
        self,
        L: Quantity,
        W: Quantity,
        t: Quantity,
        k: Quantity
    ) -> None:
        """
        Creates a `StraightFin` instance.

        Parameters
        ----------
        L: PlainQuantity
            The length of the fin.
        W: PlainQuantity
            The width of the fin.
        t: PlainQuantity
            The thickness of the fin (at the base of the fin).
        k: PlainQuantity
            The conductivity of the fin material.
        """
        self._L = L.to('m').m
        self._W = W.to('m').m
        self._t = t.to('m').m
        self._k = k.to('W / (m * K)').m
        super().__init__()

    @property
    def _m(self) -> float:
        return np.sqrt(2 * self._h_avg / (self._k * self._t))

    @property
    @abstractmethod
    def surface_area(self) -> Quantity:
        ...

    @property
    @abstractmethod
    def efficiency(self) -> float:
        ...

    @property
    def base_area(self) -> Quantity:
        A_b = self._W * self._t
        return Q_(A_b, 'm ** 2')


class StraightRectangularFin(StraightFin):

    @property
    def surface_area(self) -> Quantity:
        L_cor = self._L + self._t / 2
        A_f = 2 * L_cor * self._W
        return Q_(A_f, 'm ** 2')

    @property
    def efficiency(self) -> float:
        L_cor = self._L + self._t / 2
        n = np.tanh(self._m * L_cor)
        d = self._m * L_cor
        eta_fin = n / d
        return eta_fin


class StraightTriangularFin(StraightFin):

    @property
    def surface_area(self) -> Quantity:
        A_f = 2 * self._W * np.sqrt(self._L ** 2 + (self._t / 2) ** 2)
        return Q_(A_f, 'm ** 2')

    @property
    def efficiency(self) -> float:
        n = bessel_i(1, 2 * self._m * self._L)
        d = self._m * self._L * bessel_i(0, 2 * self._m * self._L)
        eta_fin = n / d
        return eta_fin


class StraightParabolicFin(StraightFin):

    @property
    def surface_area(self) -> Quantity:
        C1 = np.sqrt(1 + (self._t / self._L) ** 2)
        a = C1 * self._L
        b = (self._L ** 2) / self._t
        c = np.log(self._t / self._L + C1)
        A_f = self._W * (a + b * c)
        return Q_(A_f, 'm ** 2')

    @property
    def efficiency(self) -> float:
        n = 2
        d = np.sqrt(4 * (self._m * self._L) ** 2 + 1) + 1
        eta_fin = n / d
        return eta_fin


class SpineFin(Fin):
    """General base class to represent spine fins or pin fins."""

    def __init__(
        self,
        L: Quantity,
        D: Quantity,
        k: Quantity
    ) -> None:
        """
        Creates a `SpineFine` instance.

        Parameters
        ----------
        L: PlainQuantity
            The length of the fin.
        D: PlainQuantity
            The diameter of the fin (at the base of the fin).
        k: PlainQuantity
            The conductivity of the fin material.
        """
        self._L = L.to('m').m
        self._D = D.to('m').m
        self._k = k.to('W / (m * K)').m
        super().__init__()

    @property
    def _m(self) -> float:
        return np.sqrt(4 * self._h_avg / (self._k * self._D))

    @property
    @abstractmethod
    def surface_area(self) -> Quantity:
        ...

    @property
    @abstractmethod
    def efficiency(self) -> float:
        ...

    @property
    def base_area(self) -> Quantity:
        A_b = np.pi * self._D ** 2 / 4
        return Q_(A_b, 'm ** 2')


class SpineRectangularFin(SpineFin):

    @property
    def surface_area(self) -> Quantity:
        L_cor = self._L + self._D / 4
        A_f = np.pi * self._D * L_cor
        return Q_(A_f, 'm ** 2')

    @property
    def efficiency(self) -> float:
        L_cor = self._L + self._D / 4
        n = np.tanh(self._m * L_cor)
        d = self._m * L_cor
        eta_fin = n / d
        return eta_fin


class SpineTriangularFin(SpineFin):

    @property
    def surface_area(self) -> Quantity:
        a = np.pi * self._D / 2
        b = np.sqrt(self._L ** 2 + (self._D / 2) ** 2)
        A_f = a * b
        return Q_(A_f, 'm ** 2')

    @property
    def efficiency(self) -> float:
        n = 2 * bessel_i(2, 2 * self._m * self._L)
        d = self._m * self._L * bessel_i(1, 2 * self._m * self._L)
        eta_fin = n / d
        return eta_fin


class SpineParabolicFin(SpineFin):

    @property
    def surface_area(self) -> Quantity:
        C1 = 1 + 2 * (self._D / self._L) ** 2
        C2 = np.sqrt(1 + (self._D / self._L) ** 2)
        a = np.pi * (self._L ** 3) / (8 * self._D)
        b = C1 * C2
        c = self._L / (2 * self._D)
        d = np.log(2 * self._D * C2 / self._L + C1)
        A_f = a * (b - c * d)
        return Q_(A_f, 'm ** 2')

    @property
    def efficiency(self) -> float:
        n = 2
        d = np.sqrt(4 / 9 * (self._m * self._L) ** 2 + 1) + 1
        eta_fin = n / d
        return eta_fin


class CircularRectangularFin(Fin):
    """Class that represents a circular fin with a rectangular cross-section
    around a circular tube.
    """
    def __init__(
        self,
        r_i: Quantity,
        r_o: Quantity,
        t: Quantity,
        k: Quantity
    ) -> None:
        """Creates an `CircularRectangularFin` instance.

        Parameters
        ----------
        r_i: PlainQuantity
            The inner radius of the fin.
        r_o: PlainQuantity
            The outer radius of the fin.
        t: PlainQuantity
            The thickness of the fin.
        k: PlainQuantity
            The conductivity of the fin material.
        """
        self._r_i = r_i.to('m').m
        self._r_o = r_o.to('m').m
        self._t = t.to('m').m
        self._k = k.to('W / (m * K)').m
        super().__init__()

    @property
    def _m(self) -> float:
        return np.sqrt(2 * self._h_avg / (self._k * self._t))

    @property
    def surface_area(self) -> Quantity:
        r_o_cor = self._r_o + self._t / 2
        A_f = 2 * np.pi * (r_o_cor ** 2 - self._r_i ** 2)
        return Q_(A_f, 'm ** 2')

    @property
    def efficiency(self) -> float:
        r_o_cor = self._r_o + self._t / 2
        n1 = 2 * self._m * self._r_i
        d1 = (self._m * r_o_cor) ** 2 - (self._m * self._r_i) ** 2
        n2 = (
            bessel_k(1, self._m * self._r_i) * bessel_i(1, self._m * r_o_cor)
            - bessel_i(1, self._m * self._r_i) * bessel_k(1, self._m * r_o_cor)
        )
        d2 = (
            bessel_i(0, self._m * self._r_i) * bessel_k(1, self._m * r_o_cor)
            + bessel_k(0, self._m * self._r_i) * bessel_i(1, self._m * r_o_cor)
        )
        eta_fin = (n1 / d1) * (n2 / d2)
        return eta_fin

    @property
    def base_area(self) -> Quantity:
        A_b = 2 * np.pi * self._r_i * self._t
        return Q_(A_b, 'm ** 2')


class ConstantCrossSectionFin(Fin):
    """Class that represents a fin with constant cross-section."""

    def __init__(
        self,
        L: Quantity,
        P: Quantity,
        A_c: Quantity,
        k: Quantity
    ) -> None:
        """Creates a `ConstantCrossSectionFin` instance.

        Parameters
        ----------
        L: PlainQuantity
            The length of the fin.
        P: PlainQuantity
            The perimeter of the fin.
        A_c: PlainQuantity
            The cross-section area of the fin.
        k: PlainQuantity
            The conductivity of the fin material.
        """
        self._L = L.to('m').m
        self._P = P.to('m').m
        self._A_c = A_c.to('m ** 2').m
        self._k = k.to('W / (m * K)').m
        super().__init__()

    @property
    def _m(self) -> float:
        return np.sqrt(self._P * self._h_avg / (self._k * self._A_c))

    @property
    def surface_area(self) -> Quantity:
        L_cor = self._L + self._A_c / self._P
        A_f = self._P * L_cor
        return Q_(A_f, 'm ** 2')

    @property
    def efficiency(self) -> float:
        L_cor = self._L + self._A_c / self._P
        n = np.tanh(self._m * L_cor)
        d = self._m * L_cor
        eta_fin = n / d
        return eta_fin

    @property
    def base_area(self) -> Quantity:
        return Q_(self._A_c, 'm ** 2')


class PlainContinuousFin(Fin):
    """Represents a metal sheet pierced by round tubes in either a square inline
    or equilateral triangular staggered array.

    Fin efficiency is calculated according to Schmidt. The metal sheet is
    divided in "fin areas" around each tube which can be of rectangular or
    hexagonal shape.

    References
    ----------
    [1]     Thulukkanam, K. (2013). Heat Exchanger Design Handbook (2nd ed.).
            CRC Press. (§ 4.7.2.2)
    """
    def __init__(
        self,
        r_i: Quantity,
        t: Quantity,
        k: Quantity,
        M: Quantity,
        L: Quantity,
        fin_area_shape: str = 'hexagonal'
    ) -> None:
        """Creates an `PlainContinuousFin` instance. The fin efficiency is
        calculated according to the Schmidt method.

        Parameters
        ----------
        r_i: PlainQuantity
            The inner radius of the fin (= outer radius of the tube).
        t: PlainQuantity
            The thickness of the fin.
        k: PlainQuantity
            The conductivity of the fin material.
        M: PlainQuantity
            Half-height of rectangular fin area in case of an inline square
            tube array, or half-height of hexagonal fin area in case of
            staggered equilateral triangular tube array.
        L: PlainQuantity
            Half-height of rectangular fin area in case of inline square tube
            array, or distance between center and one side of a hexagonal fin
            area in case of staggered equilateral triangular tube array.
        fin_area_shape: {'rectangular', 'hexagonal' (default)}
            Shape of the fin area: 'rectangular' if inline square tube array,
            or 'hexagonal' if staggered equilateral triangular tube array.

        Notes
        -----
        The accuracy is especially good for square or regular hexagonal fins,
        for which L = M.
        """
        self._r_i = r_i.to('m').m
        self._t = t.to('m').m
        self._k = k.to('W / (m * K)').m
        self._M = M.to('m').m
        self._L = L.to('m').m
        self._lambda = self._M / self._r_i
        self._beta = self._L / self._M
        self._fin_area_shape = fin_area_shape
        super().__init__()

    @property
    def _m(self) -> float:
        return np.sqrt(2 * self._h_avg / (self._k * self._t))

    @property
    def _rho_eq(self) -> float:
        if self._fin_area_shape == 'rectangular':
            rho_eq = 1.28 * self._lambda * np.sqrt(self._beta - 0.2)
        else:  # self._fin_area_shape == 'hexagonal'
            rho_eq = 1.27 * self._lambda * np.sqrt(self._beta - 0.3)
        return rho_eq

    @property
    def _r_o_eq(self) -> float:
        """Equivalent fin radius."""
        r_o_eq = self._rho_eq * self._r_i
        return r_o_eq

    @property
    def efficiency(self) -> float:
        L_cor = (
            (self._r_o_eq - self._r_i)
            * (1 + self._t / (2 * (self._r_o_eq - self._r_i)))
            * (1 + 0.35 * np.log(self._rho_eq))
        )
        n = np.tanh(self._m * L_cor)
        d = self._m * L_cor
        eta_fin = n / d
        return eta_fin

    @property
    def surface_area(self) -> Quantity:
        r_o_cor = self._r_o_eq + self._t / 2
        A_f = 2 * np.pi * (r_o_cor ** 2 - self._r_i ** 2)
        return Q_(A_f, 'm ** 2')

    @property
    def base_area(self) -> Quantity:
        """Gets the cross-sectional area of the fin base."""
        A_b = 2 * np.pi * self._r_i * self._t
        return Q_(A_b, 'm ** 2')
