import math
from abc import abstractmethod
from dataclasses import dataclass
from hvac import Quantity
from ..geometry import ExternalGeometry, StaggeredTubeBank, InlineTubeBank


Q_ = Quantity
PI = math.pi


class ContinuousFinExternalGeometry(ExternalGeometry):

    @property
    def A_pri(self) -> Quantity:
        d_o = self.parent.d_o.to('m')
        d_r = self.parent.d_r.to('m')
        w = self.parent.width.to('m')
        h = self.parent.height.to('m')
        l = self.parent.length.to('m')
        t_f = self.parent.t_fin.to('m')
        n_f = self.parent.fin_density.to('1 / m')
        n_t = self.parent.num_tubes
        A = PI * d_r * w * (1 - t_f * n_f) * n_t
        A += 2 * (l * h - (PI * d_o ** 2 / 4) * n_t)
        return A

    @property
    def A_fin(self) -> Quantity:
        w = self.parent.width.to('m')
        l = self.parent.length.to('m')
        h = self.parent.height.to('m')
        d_r = self.parent.d_r.to('m')
        n_t = self.parent.num_tubes
        n_f = self.parent.fin_density.to('1 / m')
        t_f = self.parent.t_fin.to('m')
        a = 2 * (l * h - (PI * d_r ** 2 / 4) * n_t) * n_f * w
        b = 2 * h * t_f * n_f * w
        A = a + b
        return A

    @property
    def A_tot(self) -> Quantity:
        A_p = self.A_pri.to('m ** 2')
        A_f = self.A_fin.to('m ** 2')
        A = A_p + A_f
        return A

    @property
    def ratio_A_fin_to_tot(self) -> float:
        d_r = self.parent.d_r.to('m')
        t_f = self.parent.t_fin.to('m')
        n_f = self.parent.fin_density.to('1 / m')
        p_l = self.parent.pitch_lon.to('m')
        p_t = self.parent.pitch_trv.to('m')
        d_o = self.parent.d_o.to('m')
        # prime area of 1 tube with length 1 m:
        A_p = PI * d_r * Q_(1, 'm') * (1 - t_f * n_f)
        A_p += 2 * (p_l * p_t - (PI * d_o ** 2 / 4))
        # finned area of 1 tube with length 1 m:
        a = 2 * (p_l * p_t - (PI * d_r ** 2 / 4)) * n_f * Q_(1, 'm')
        b = 2 * p_t * t_f * n_f * Q_(1, 'm')
        A_f = a + b
        Af_to_A = (1 + A_p / A_f) ** -1
        return Af_to_A.m

    @property
    @abstractmethod
    def A_min(self) -> Quantity:
        ...


class ContinuousFinInlineExternalGeometry(ContinuousFinExternalGeometry):

    @property
    def A_min(self) -> Quantity:
        p_t = self.parent.pitch_trv.to('m')
        d_r = self.parent.d_r.to('m')
        w = self.parent.width.to('m')
        h = self.parent.height.to('m')
        t_f = self.parent.t_fin.to('m')
        n_f = self.parent.fin_density.to('1 / m')
        a = (p_t - d_r) * w
        b = (p_t - d_r) * t_f * n_f * w
        A = (a - b) * (h / p_t)
        return A


class ContinuousFinStaggeredExternalGeometry(ContinuousFinExternalGeometry):

    @property
    def A_min(self) -> Quantity:
        p_t = self.parent.pitch_trv.to('m')
        p_l = self.parent.pitch_lon.to('m')
        d_r = self.parent.d_r.to('m')
        t_f = self.parent.t_fin.to('m')
        n_f = self.parent.fin_density.to('1 / m')
        w = self.parent.width.to('m')
        h = self.parent.height.to('m')
        x = 0.5 * (p_t - d_r) * (1 - t_f * n_f)
        p_d = ((p_t / 2) ** 2 + p_l ** 2) ** 0.5
        y = p_d - d_r - (p_t - d_r) * t_f * n_f
        z = 2 * min(x, y)
        a = (h / p_t - 1) * z
        b = p_t - d_r
        c = (p_t - d_r) * t_f * n_f
        A = (a + b - c) * w
        return A


@dataclass
class ContinuousFinStaggeredTubeBank(StaggeredTubeBank):

    def __post_init__(self):
        super().__post_init__()
        self.external = ContinuousFinStaggeredExternalGeometry(self)


@dataclass
class ContinuousFinInlineTubeBank(InlineTubeBank):

    def __post_init__(self):
        super().__post_init__()
        self.external = ContinuousFinInlineExternalGeometry(self)
