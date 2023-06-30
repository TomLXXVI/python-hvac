"""Heat transfer correlations for corrugated flat finned, staggered tube banks.

References
----------
[1]     Shah, R. K., & Sekulic, D. P. (2003). Fundamentals of Heat Exchanger
        Design. John Wiley & Sons. (§ 7.5.4.3)
"""
import math


def j_Wang_large(
    Re_D_r: float,
    D_r: float,
    S_l: float,
    t_f: float,
    p_d: float,
    x_f: float,
    N_f: float,
    N_r: float
) -> float:
    """Colburn j-factor correlation by Wang (2000) for corrugated (herringbone
    or sharp wave) flat fins on staggered tube banks with large-diameter tubes
    (12.7 - 15.88 mm).

    Parameters
    ----------
    Re_D_r:
        Reynolds number based on collar diameter.
    D_r:
        Collar diameter [m].
    S_l:
        Longitudinal pitch, i.e. spacing between tubes of two adjacent tube
        rows [m].
    t_f:
        Thickness of fin [m].
    p_d:
        Fin pattern depth (peak-to-valley distance, excluding fin thickness) [m].
    x_f:
        Projected fin pattern length for one-half wavelength [m]
    N_f:
        The number of fins per unit length [1 / m]
    N_r:
        Number of rows.

    Notes
    -----
    Valid ranges of parameters:
    Re_D_r: 500 - 10_000
    D_r: 13.6 - 16.85 mm
    D_h: 3.63 - 7.23 mm
    S_t: 31.75 - 38.1 mm
    S_l: 27.5 - 33 mm
    p_f: 2.98 - 6.34 mm
    N_r: 1 - 6
    theta: 12.3 - 14.7 deg
    x_f: 6.87 - 8.25 mm
    p_d: 1.8 mm

    Returns
    -------
    Colburn j-factor.
    """
    s = (1 - N_f * t_f) / N_f  # length of bare tube between 2 fins
    p_f = s + t_f  # fin pitch
    k1 = (S_l / t_f) ** -0.493
    k2 = (p_f / D_r) ** -0.886
    k3 = (p_d / x_f) ** -0.0296
    c1 = -0.1707 - 1.374 * k1 * k2 * (N_r ** -0.134) * k3
    k8 = (S_l / t_f) ** -0.456
    k9 = N_r ** -0.27
    k10 = (p_f / D_r) ** -1.343
    k11 = (p_d / x_f) ** 0.317
    j = 1.7910 * (Re_D_r ** c1) * k8 * k9 * k10 * k11
    return j


def j_Wang_small(
    Re_D_r: float,
    D_r: float,
    D_h: float,
    S_t: float,
    S_l: float,
    t_f: float,
    N_f: float,
    N_r: float,
    theta: float
) -> float:
    """Colburn j-factor correlation by Wang (2000) for corrugated (herringbone
    or sharp wave) flat fins on staggered tube banks with smaller-diameter tubes
    (7.94 - 9.53 mm).

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
    t_f:
        Thickness of fin [m].
    N_f:
        The number of fins per unit length [1 / m]
    N_r:
        Number of rows.
    theta:
        Slope angle of waves [deg]

    Notes
    -----
    Valid ranges of parameters:
    Re_D_r: 300 - 8000
    D_r: 8.58 - 10.38 mm
    D_h: 1.53 - 4.52 mm
    S_t: 25.4 mm
    S_l: 19.05 - 25.04 mm
    p_f: 1.21 - 3.66 mm
    N_r: 1 - 6
    theta: 14.5 - 18.5 deg
    x_f: 4.76 - 6.35 mm
    p_d: 1.18 - 1.68 mm

    Returns
    -------
    Colburn j-factor.
    """
    s = (1 - N_f * t_f) / N_f  # length of bare tube between 2 fins
    p_f = s + t_f  # fin pitch
    theta = math.radians(theta)
    k1 = (p_f / D_r) ** 0.6
    k2 = (S_l / D_h) ** 0.54
    k3 = N_r ** -0.284
    k4 = math.log(0.5 * math.tan(theta))
    c1 = -0.229 + 0.115 * k1 * k2 * k3 * k4
    c2 = -0.251 + 0.232 * (N_r ** 1.37) / (math.log(Re_D_r) - 2.303)
    k5 = (p_f / D_h) ** 0.09
    k6 = (S_l / S_t) ** -1.75
    k7 = N_r ** -0.93
    c3 = -0.439 * k5 * k6 * k7
    c4 = 0.502 * (math.log(Re_D_r) - 2.54)
    k12 = (p_f / S_l) ** c2
    k13 = math.tan(theta) ** c3
    k14 = (S_l / S_t) ** c4
    k15 = N_r ** 0.428
    j = 0.324 * (Re_D_r ** c1) * k12 * k13 * k14 * k15
    return j
