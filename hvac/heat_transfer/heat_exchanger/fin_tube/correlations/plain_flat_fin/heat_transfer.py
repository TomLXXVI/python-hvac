"""Heat transfer correlations for plain flat finned, staggered tube banks.

References
----------
[1]     Shah, R. K., & Sekulic, D. P. (2003). Fundamentals of Heat Exchanger
        Design. John Wiley & Sons. (§ 7.5.4.2)
"""
import math


def j_Wang_and_Chi(
    Re_D_r: float,
    D_r: float,
    D_h: float,
    S_t: float,
    S_l: float,
    N_f: float,
    N_r: float
) -> float:
    """Colburn j-factor correlation by Wang and Chi (2000) for plain flat
    fins on staggered tube banks.

    Parameters
    ----------
    Re_D_r:
        Reynolds number based on collar diameter.
    D_r:
        Collar diameter [m].
    D_h:
        Hydraulic diameter [m].
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
    Colburn j-factor.
    """
    N_r = min(N_r, 6)
    p_f = 1 / N_f  # fin pitch
    ln_Re = math.log(Re_D_r)
    if N_r < 2:
        c1 = 1.9 - 0.23 * ln_Re
        c2 = -0.236 + 0.126 * ln_Re
        k1 = (S_t / S_l) ** c1
        k2 = (p_f / D_r) ** -1.084
        k3 = (p_f / D_h) ** -0.786
        k4 = (p_f / S_t) ** c2
        j = 0.108 * (Re_D_r ** -0.29) * k1 * k2 * k3 * k4
    else:
        c3 = (
            -0.361 - (0.042 * N_r / ln_Re)
            + 0.158 * math.log(N_r * (p_f / D_r) ** 0.41)
        )
        c4 = -1.224 - (0.076 * (S_l / D_h) ** 1.42) / ln_Re
        c5 = -0.083 + 0.058 * N_r / ln_Re
        c6 = -5.735 + 1.21 * math.log(Re_D_r / N_r)
        k5 = N_r ** c4
        k6 = (p_f / D_r) ** c5
        k7 = (p_f / D_h) ** c6
        k8 = (p_f / S_t) ** -0.93
        j = 0.086 * (Re_D_r ** c3) * k5 * k6 * k7 * k8
    return j
