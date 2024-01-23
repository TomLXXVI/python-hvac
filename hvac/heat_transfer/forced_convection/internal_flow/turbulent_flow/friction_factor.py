"""
Implementation of the correlations for fully developed and average friction
factor in case of turbulent flow.

The correlations were taken from:
Nellis G. F., & Klein S. A. (2021)
INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.
"""
import math
from hvac import Quantity


class SmoothTube:

    @staticmethod
    def fully_developed_local_friction_factor(Re: float) -> float:
        """Correlation of Li et al. (2011) valid for:
        - e / Dh = 0, and
        - 4000 < Re < 1.0e7
        """
        C1 = -0.0015702
        C2 = 0.3942031
        C3 = 2.5341533
        x = math.log(Re)
        f_fd = 4 * (C1 / x + C2 / (x ** 2) + C3 / (x ** 3))
        return f_fd

    @staticmethod
    def average_friction_factor(
        Re: float,
        Dh: Quantity,
        L: Quantity
    ) -> float:
        f_fd = SmoothTube.fully_developed_local_friction_factor(Re)
        k = 1 + (Dh / L) ** 0.7
        f_avg = k * f_fd
        return f_avg


class RoughTube:

    @staticmethod
    def fully_developed_local_friction_factor(Re: float, e_r: float) -> float:
        """Correlation of Offor and Alabi (2016) valid for:
        - 1.0e-6 < e / Dh < 5.0e-2, and
        - 4000 < Re < 1.0e8
        """
        C1 = 3.71
        C2 = 1.975
        C3 = 3.93
        C4 = 1.092
        C5 = 7.627
        C6 = 395.9
        C7 = -2.0
        x1 = e_r / C1
        x2 = C2 / Re
        x3 = e_r / C3
        x4 = C5 / (Re + C6)
        x = x1 - x2 * (math.log(x3 ** C4 + x4))
        x = -2.0 * math.log10(x)
        f_fd = x ** C7
        return f_fd

    @staticmethod
    def average_friction_factor(
        Re: float,
        Dh: Quantity,
        L: Quantity,
        e: Quantity
    ) -> float:
        Dh = Dh.to('m').m
        L = L.to('m').m
        e = e.to('m').m
        e_r = e / Dh
        f_fd = RoughTube.fully_developed_local_friction_factor(Re, e_r)
        k = 1 + (Dh / L) ** 0.7
        f_avg = k * f_fd
        return f_avg

    @staticmethod
    def fully_rough_local_friction_factor(e_r: float) -> float:
        d = (1.74 + 2.0 * math.log10(0.5 * 1 / e_r)) ** 2
        f_fd = 1 / d
        return f_fd
