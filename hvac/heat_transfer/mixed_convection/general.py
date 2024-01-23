"""
Mixed convection.

From:
Nellis G. F., & Klein S. A. (2021)
INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.
"""
from hvac import Quantity


def avg_heat_trf_coeff(
    h_avg_fc: Quantity,
    h_avg_nc: Quantity,
    sign: str,
    m: float = 3.0
) -> Quantity:
    """Returns the average heat transfer coefficient for mixed convection.

    Mixed convection applies if the ratio of the Grashof number for natural
    convection to the square of the Reynolds number for forced convection is
    around 1. If much less than 1, only natural convection is present. If much
    greater than 1, natural convection is negligible.

    Parameters
    ----------
    h_avg_fc:
        Average heat transfer coefficient in case of only forced convection.
    h_avg_nc:
        Average heat transfer coefficient in case of only natural or free
        convection.
    sign: str, {'+', '-'}
        The plus sign is used if the natural and forced convection flows
        augment one another (if the velocities have components in the same
        direction, i.e. the angle between the velocity vectors is between 0°
        and 90° included).
    m: float, default 3.0
        Exponent, usually taken to be 3.

    Returns
    -------
    h_avg: Quantity
        Average heat transfer coefficient for mixed convection.
    """
    a = max(h_avg_fc, h_avg_nc) ** m
    b = min(h_avg_fc, h_avg_nc) ** m
    if sign == '+':
        h_avg = (a + b) ** (1 / m)
    else:
        h_avg = (a - b) ** (1 / m)
    return h_avg
