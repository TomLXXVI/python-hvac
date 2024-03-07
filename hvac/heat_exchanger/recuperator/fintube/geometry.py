from math import pi as PI
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from hvac import Quantity


Q_ = Quantity


@dataclass
class TubeBank(ABC):
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
        ...

    @property
    def num_tubes_1st_row(self) -> float:
        h = self.height.to('m').m
        p_t = self.pitch_trv.to('m').m
        n_1r = h / p_t
        return n_1r

    @property
    def volume(self) -> Quantity:
        w = self.width.to('m')
        h = self.height.to('m')
        l = self.length.to('m')
        V = w * h * l
        return V

    @property
    def length(self) -> Quantity:
        return self._length

    @length.setter
    def length(self, v: Quantity) -> None:
        self._length = v
        self.num_rows = self._length / self.pitch_lon


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
    parent: TubeBank

    @property
    def A_tot(self) -> Quantity:
        n_t = self.parent.num_tubes
        d_i = self.parent.d_i.to('m')
        w = self.parent.width.to('m')
        A = PI * d_i * w * n_t
        return A

    @property
    def A_min(self) -> Quantity:
        d_i = self.parent.d_i.to('m')
        A_tube = PI * d_i ** 2 / 4
        A = A_tube * self.parent.num_tubes
        return A

    @property
    def A_face(self) -> Quantity:
        l = self.parent.length.to('m').m
        h = self.parent.height.to('m').m
        A = l * h
        return A

    @property
    def d_h(self) -> Quantity:
        return self.parent.d_i.to('m')

    @property
    def sigma(self) -> float:
        """Ratio of minimum free flow area to frontal area."""
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
        Ratio of the total heat transfer area to the volume of the heat
        exchanger core.
        """
        return self.A_tot / self.parent.volume


@dataclass
class ExternalGeometry(ABC):
    parent: TubeBank

    @property
    @abstractmethod
    def A_pri(self) -> Quantity:
        ...

    @property
    @abstractmethod
    def A_fin(self) -> Quantity:
        ...

    @property
    @abstractmethod
    def A_tot(self) -> Quantity:
        pass

    @property
    @abstractmethod
    def A_min(self) -> Quantity:
        pass

    @property
    def A_face(self) -> Quantity:
        w = self.parent.width.to('m')
        h = self.parent.height.to('m')
        A = w * h
        return A

    @property
    def d_h(self) -> Quantity:
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
