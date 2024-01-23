"""
Implementation of the correlations for fully developed and average Nusselt number
in case of turbulent flow.

The correlations were taken from:
Nellis G. F., & Klein S. A. (2021)
INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.
"""
import math
from hvac import Quantity


def fully_developed_local_nusselt_number(
    f_fd: float,
    Re: float,
    Pr: float
) -> float:
    """Gnielinski correlation valid for:
    - 0.5 < Pr < 2000, and
    - 2300 < Re < 5.e6
    """
    n = (f_fd / 8) * (Re - 1000) * Pr
    d = 1 + 12.7 * ((Pr ** (2 / 3)) - 1) * math.sqrt(f_fd / 8)
    Nu_fd = n / d
    return Nu_fd


def average_nusselt_number(
    f_fd: float,
    Re: float,
    Pr: float,
    Dh: Quantity,
    L: Quantity
) -> float:
    Nu_fd = fully_developed_local_nusselt_number(f_fd, Re, Pr)
    Dh = Dh.to('m').m
    L = L.to('m').m
    C, m = 1.0, 0.7
    try:
        Nu_avg = Nu_fd * (1 + C * (L / Dh) ** -m)
    except Exception:
        pass
    return Nu_avg
