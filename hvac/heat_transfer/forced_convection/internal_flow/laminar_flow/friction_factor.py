"""
Implementation of the correlations for fully developed and average friction
factor in case of laminar flow.

The correlations were taken from:
Nellis G. F., & Klein S. A. (2021)
INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.
"""
import math
from hvac import Quantity


class CircularTube:

    @staticmethod
    def fully_developed_local_friction_factor(Re: float) -> float:
        f = 64 / Re
        return f

    @staticmethod
    def average_friction_factor(L_tilde: float, Re: float) -> float:
        f_fd = CircularTube.fully_developed_local_friction_factor(Re)
        r = 0.215 * (L_tilde ** -0.5)
        n = 0.0196 * (L_tilde ** -1.0) + 1.0 - 0.215 * (L_tilde ** -0.5)
        d = 1 + 0.00021 * (L_tilde ** -2.0)
        r += n / d
        f_avg = r * f_fd
        return f_avg


class RectangularTube:

    @staticmethod
    def fully_developed_local_friction_factor(
        Re: float,
        a: Quantity,
        b: Quantity
    ) -> float:
        a = a.to('m').m
        b = b.to('m').m
        ar = min(a, b) / max(a, b)  # aspect ratio
        r = 96 * (
            1
            - 1.3553 * ar
            + 1.9467 * ar ** 2
            - 1.7012 * ar ** 3
            + 0.9564 * ar ** 4
            - 0.2537 * ar ** 5
        )
        f = r / Re
        return f

    @staticmethod
    def average_friction_factor(
        L_tilde: float,
        Re: float,
        a: Quantity,
        b: Quantity
    ) -> float:
        f_fd = RectangularTube.fully_developed_local_friction_factor(Re, a, b)
        r = 13.76 * (L_tilde ** -0.5)
        n = 1.25 * (L_tilde ** -1.0) + (f_fd * Re) - 13.76 * (L_tilde ** -0.5)
        d = 1 + 0.00021 * (L_tilde ** -2.0)
        r += n / d
        f_avg = r / Re
        return f_avg


class AnnularDuct:

    @staticmethod
    def fully_developed_local_friction_factor(
        Re: float,
        ri: Quantity,
        ro: Quantity
    ) -> float:
        ri = ri.to('m').m
        ro = ro.to('m').m
        rr = ri / ro
        n = 64.0 * (1 - rr) ** 2
        d = 1 + (rr ** 2) - ((1 - rr ** 2) / math.log(1 / rr))
        r = n / d
        f_fd = r / Re
        return f_fd

    @staticmethod
    def average_friction_factor(
        L_tilde: float,
        Re: float,
        ri: Quantity,
        ro: Quantity
    ) -> float:
        f_fd = AnnularDuct.fully_developed_local_friction_factor(Re, ri, ro)
        r = 13.76 * (L_tilde ** -0.5)
        n = 1.25 * (L_tilde ** -1.0) + (f_fd * Re) - 13.76 * (L_tilde ** -0.5)
        d = 1 + 0.00021 * (L_tilde ** -2.0)
        r += n / d
        f_avg = r / Re
        return f_avg
