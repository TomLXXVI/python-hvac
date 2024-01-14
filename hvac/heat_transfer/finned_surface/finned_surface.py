from hvac import Quantity
from .fins import Fin

Q_ = Quantity


class FinnedSurface:
    """Class that represents a finned surface."""

    def __init__(self, A_b: Quantity, fin: Fin, N_fins: int) -> None:
        """Creates a `FinnedSurface` instance.

        Parameters
        ----------
        A_b: PlainQuantity
            The base area of the finned surface.
        fin: Fin
            The type of fin the finned surface is made of.
        N_fins: int
            The number of fins on the surface.
        """
        self.A_b = A_b.to('m ** 2')
        self.A_p = self.A_b - N_fins * fin.base_area  # prime area
        self.fin = fin
        self.A_f = N_fins * fin.surface_area
        self.A_t = self.A_p + self.A_f
        self._h_avg: float = 0.0
        self._T_b: float = 0.0
        self._T_fluid: float = 0.0

    @property
    def h_avg(self) -> Quantity:
        """Gets the average heat transfer coefficient between the finned surface
        and the surrounding fluid.
        """
        return Q_(self._h_avg, 'W / (m ** 2 * K)')

    @h_avg.setter
    def h_avg(self, v: Quantity) -> None:
        """Sets the average heat transfer coefficient between the finned surface
        and the surrounding fluid.
        """
        self._h_avg = v.to('W / (m ** 2 * K)').m
        self.fin.h_avg = v

    @property
    def T_b(self) -> Quantity:
        """Gets the base temperature of the finned surface."""
        return Q_(self._T_b, 'K')

    @T_b.setter
    def T_b(self, v: Quantity) -> None:
        """Sets the base temperature of the finned surface."""
        self._T_b = v.to('K').m

    @property
    def T_fluid(self) -> Quantity:
        """Gets the temperature of the fluid (at some distance from the finned
        surface).
        """
        return Q_(self._T_fluid, 'K')

    @T_fluid.setter
    def T_fluid(self, v: Quantity) -> None:
        """Sets the temperature of the fluid (at some distance from the finned
         surface).
         """
        self._T_fluid = v.to('K').m

    @property
    def resistance(self) -> Quantity:
        """Gets the thermal resistance between the finned surface and the
        surrounding fluid."""
        eta_fin = self.fin.efficiency
        U = self._h_avg * (eta_fin * self.A_f.m + self.A_p.m)
        R = 1 / U
        return Q_(R, 'K / W')

    @property
    def heat_transfer(self) -> Quantity:
        """Gets the heat transfer between the finned surface and the surrounding
        fluid."""
        R = self.resistance.m
        Q = (self._T_b - self._T_fluid) / R
        return Q_(Q, 'W')

    @property
    def efficiency(self) -> float:
        """Gets the overall surface efficiency of the finned surface, i.e. the
        ratio between the actual heat transfer and the heat transfer that would
        occur if the finned surface had a uniform temperature equal to the base
        temperature."""
        r = self.A_f.m / self.A_t.m
        eta_fin = self.fin.efficiency
        eta_o = 1 - r * (1 - eta_fin)
        return eta_o
