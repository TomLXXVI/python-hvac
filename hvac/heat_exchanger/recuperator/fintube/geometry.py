import math
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from hvac import Quantity


Q_ = Quantity
PI = math.pi


@dataclass
class TubeBank(ABC):
    """
    Models a general tube bank.

    Attributes
    ----------
    width: Quantity
        Frontal area width.
    height: Quantity
        Frontal area height.
    num_rows: int | None
        Number of rows.
    pitch_trv: Quantity
        Vertical spacing between tubes in one row (transversal pitch).
    pitch_lon: Quantity
        Horizontal spacing between the tubes in two adjacent rows (longitudinal
        pitch).
    d_o: Quantity
        Tube external diameter.
    d_i: Quantity
        Tube internal diameter.
    t_fin: Quantity
        Fin thickness.
    fin_density: Quantity
        Number of fins per unit length of tube (= inverse of fin spacing).
    k_fin: Quantity, default Q_(237, 'W / (m * K)')
        Thermal conductivity of fin material. Default value applies to aluminum.
    d_r: Quantity | None = None
        TODO: add description.

    internal: InternalGeometry
        Represents the internal side of the heat exchanging surface.
    external: ExternalGeometry | None = None.
        Represents the external side of the heat exchanging surface.
    """
    width: Quantity
    height: Quantity
    num_rows: int | None
    pitch_trv: Quantity
    pitch_lon: Quantity
    d_o: Quantity
    d_i: Quantity
    t_fin: Quantity
    fin_density: Quantity
    k_fin: Quantity = Q_(237, 'W / (m * K)')
    d_r: Quantity | None = None
    _length: Quantity = field(init=False, default=None)

    def __post_init__(self):
        if self.num_rows is not None:
            self._length = self.num_rows * self.pitch_lon
        if self.d_r is None:
            self.d_r = self.d_o
        self.internal = InternalGeometry(self)
        self.external: ExternalGeometry | None = None

    @property
    @abstractmethod
    def num_tubes(self) -> float:
        """Returns the number of straight tube pieces in the tube bank."""
        ...

    @property
    def num_tubes_1st_row(self) -> float:
        """
        Returns the number of straight tube pieces in the first row of the
        tube bank.
        """
        h = self.height.to('m').m
        p_t = self.pitch_trv.to('m').m
        n_1r = h / p_t
        return n_1r

    @property
    def volume(self) -> Quantity:
        """Returns the volume taken in by the tube bank."""
        w = self.width.to('m')
        h = self.height.to('m')
        l = self.length.to('m')
        V = w * h * l
        return V

    @property
    def length(self) -> Quantity:
        """
        Returns the length of the tube bank (number of rows x
        longitudinal pitch).
        """
        return self._length

    @length.setter
    def length(self, v: Quantity) -> None:
        """
        Sets the length of the tube bank and determines the number of rows
        based on the longitudinal pitch.
        """
        self._length = v
        self.num_rows = self._length.to('mm').m / self.pitch_lon.to('mm').m


@dataclass
class StaggeredTubeBank(TubeBank):

    @property
    def num_tubes(self) -> float:
        h = self.height.to('m').m
        p_t = self.pitch_trv.to('m').m
        l = self.length.to('m').m
        p_l = self.pitch_lon.to('m').m
        n_t = h / p_t * (l / p_l + 1) / 2
        n_t += (h / p_t - 1) * (l / p_l - 1) / 2
        return n_t


@dataclass
class InlineTubeBank(TubeBank):

    @property
    def num_tubes(self) -> float:
        l = self.length.to('m').m
        h = self.height.to('m').m
        p_t = self.pitch_trv.to('m').m
        p_l = self.pitch_lon.to('m').m
        n_t = l * h / (p_t * p_l)
        return n_t


@dataclass
class InternalGeometry:
    """
    Represents the internal side of the heat exchanging surface of the tube
    bank.
    """
    parent: TubeBank

    @property
    def A_tot(self) -> Quantity:
        """Returns the total area of the internal heat exchanging surface."""
        n_t = self.parent.num_tubes
        d_i = self.parent.d_i.to('m')
        w = self.parent.width.to('m')
        A = PI * d_i * w * n_t
        return A

    @property
    def A_min(self) -> Quantity:
        """
        Returns the total cross-sectional area of the tubes in the tube bank
        (cross-sectional area of one tube x number of tubes).
        """
        d_i = self.parent.d_i.to('m')
        A_tube = PI * d_i ** 2 / 4
        A = A_tube * self.parent.num_tubes
        return A

    @property
    def A_face(self) -> Quantity:
        """
        Returns the internal side face area of the tube bank (i.e. tube bank
        length x tube bank height -> perpendicular to the tubes).
        """
        l = self.parent.length.to('m')
        h = self.parent.height.to('m')
        A = l * h
        return A

    @property
    def d_h(self) -> Quantity:
        """
        Returns the hydraulic diameter = inside diameter of the tubes.
        """
        return self.parent.d_i.to('m')

    @property
    def sigma(self) -> float:
        """
        Ratio of the minimum free flow area to frontal area.
        """
        d_i = self.parent.d_i.to('m').m
        p_t = self.parent.pitch_trv.to('m').m
        p_l = self.parent.pitch_lon.to('m').m
        A_tube = PI * d_i ** 2 / 4
        A_front = p_t * p_l
        sigma = A_tube / A_front
        return sigma

    @property
    def alpha(self) -> Quantity:
        """
        Ratio of the total inside heat transfer area to the volume of the heat
        exchanger core.
        """
        return self.A_tot / self.parent.volume


@dataclass
class ExternalGeometry(ABC):
    """
    Represents the external side of the heat exchanging surface of the tube
    bank.
    """
    parent: TubeBank

    @property
    @abstractmethod
    def A_pri(self) -> Quantity:
        """Returns the primary surface area, i.e. without the fins."""
        ...

    @property
    @abstractmethod
    def A_fin(self) -> Quantity:
        """
        Returns the finned surface area, excluding the primary surface
        area.
        """
        ...

    @property
    @abstractmethod
    def A_tot(self) -> Quantity:
        """Returns the total surface area (primary + finned surface area)."""
        ...

    @property
    @abstractmethod
    def A_min(self) -> Quantity:
        """Returns the minimum free flow area."""
        ...

    @property
    def A_face(self) -> Quantity:
        """Returns the face area (width of tube bank x height of tube bank)."""
        w = self.parent.width.to('m')
        h = self.parent.height.to('m')
        A = w * h
        return A

    @property
    def d_h(self) -> Quantity:
        """Returns the hydraulic diameter of the external surface."""
        d_h = 4 * self.sigma / self.alpha
        return d_h

    @property
    def sigma(self) -> float:
        """Ratio of minimum free flow area to frontal area."""
        A_min = self.A_min.to('m ** 2').m
        A_face = self.A_face.to('m ** 2').m
        return A_min / A_face

    @property
    def alpha(self) -> Quantity:
        """
        Ratio of the total heat transfer area to the volume of the heat
        exchanger core.
        """
        A = self.A_tot.to('m ** 2')
        V = self.parent.volume.to('m ** 3')
        return A / V
