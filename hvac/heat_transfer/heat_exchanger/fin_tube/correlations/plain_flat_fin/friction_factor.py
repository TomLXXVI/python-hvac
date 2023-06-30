"""Friction factor correlations for plain flat finned, staggered tube banks.

References
----------
[1]     Shah, R. K., & Sekulic, D. P. (2003). Fundamentals of Heat Exchanger
        Design. John Wiley & Sons. (§ 7.5.4.2)
"""
import math


def f_Wang_and_Chi(
    Re_D_r: float,
    D_r: float,
    S_t: float,
    S_l: float,
    N_f: float,
    N_r: float
) -> float:
    """Friction factor correlation by Wang and Chi (2000) for plain flat
    fins on staggered tube banks.

    Parameters
    ----------
    Re_D_r:
        Reynolds number based on collar diameter.
    D_r:
        Collar diameter [m].
    S_t:
        Lateral or transverse pitch, i.e. spacing between tubes of the
        same row [m].
    S_l:
        Longitudinal pitch, i.e. spacing between tubes of two adjacent tube
        rows [m].
    N_f:
        The number of fins per unit length [1/m]
    N_r:
        Number of rows.

    Notes
    -----
    Valid ranges of parameters:
    Re_D_r: 300 - 20_000
    D_r: 6.9 - 13.6 mm
    D_h: 1.30 - 9.37 mm
    S_t: 20.4 - 31.8 mm
    S_l: 12.7 - 32 mm
    p_f: 1.0 - 8.7 mm
    N_r: 1 - 6

    Returns
    -------
    Fanning friction factor.
    """
    N_r = min(N_r, 6)
    p_f = 1 / N_f  # fin pitch
    ln_Re = math.log(Re_D_r)
    c1 = -0.764 + 0.739 * (S_t / S_l) + 0.177 * (p_f / D_r) - 0.00758 / N_r
    c2 = -15.689 + 64.021 / ln_Re
    c3 = 1.696 - 15.695 / ln_Re
    k1 = (S_t / S_l) ** c2
    k2 = (p_f / D_r) ** c3
    f = 0.0267 * (Re_D_r ** c1) * k1 * k2
    return f
