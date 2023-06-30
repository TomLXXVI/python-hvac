"""Friction factor correlations for corrugated flat finned, staggered tube banks.

References
----------
[1]     Shah, R. K., & Sekulic, D. P. (2003). Fundamentals of Heat Exchanger
        Design. John Wiley & Sons. (§ 7.5.4.3)
"""
import math


def f_Wang_large(
    Re_D_r: float,
    D_r: float,
    D_h: float,
    S_t: float,
    S_l: float,
    t_f: float,
    p_d: float,
    x_f: float,
    N_f: float,
    N_r: float,
    A: float,
    A_pt: float
) -> float:
    """Friction correlation by Wang (2000) for corrugated (herringbone
    or sharp wave) flat fins on staggered tube banks with large-diameter tubes
    (12.7 - 15.88 mm).

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
    p_d:
        Fin pattern depth (peak-to-valley distance, excluding fin thickness) [m].
    x_f:
        Projected fin pattern length for one-half wavelength [m]
    N_f:
        The number of fins per unit length [1 / m]
    N_r:
        Number of rows.
    A:
        External heat transfer area [m ** 2]
    A_pt:
        External prime area [m ** 2], i.e. the tube outside surface area when
        there are no fins.

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
    Fanning friction factor.
    """
    s = (1 - N_f * t_f) / N_f  # length of bare tube between 2 fins
    p_f = s + t_f  # fin pitch
    k4 = (p_f / S_l) ** 0.25
    k5 = math.log(A / A_pt)
    k6 = (p_d / x_f) ** -0.2
    k7 = (p_f / S_t) ** 0.3
    c2 = 0.1714 - 0.07372 * k4 * k5 * k6
    c3 = 0.426 * k7 * k5
    c4 = -10.2192 / math.log(Re_D_r)
    k12 = (p_d / x_f) ** c3
    k13 = (p_f / S_t) ** c4
    k14 = math.log(A / A_pt) ** -2.726
    k15 = (D_h / D_r) ** 0.1325
    k16 = N_r ** 0.02305
    f = 0.05273 * (Re_D_r ** c2) * k12 * k13 * k14 * k15 * k16
    return f


def f_Wang_small(
    Re_D_r: float,
    D_r: float,
    D_h: float,
    S_t: float,
    S_l: float,
    t_f: float,
    N_f: float,
    N_r: float,
    theta: float,
    A: float,
    A_pt: float
) -> float:
    """Friction factor correlation by Wang (2000) for corrugated (herringbone
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
    A:
        External heat transfer area [m ** 2]
    A_pt:
        External prime area [m ** 2], i.e. the tube outside surface area when
        there are no fins.

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
    Fanning friction factor.
    """
    s = (1 - N_f * t_f) / N_f  # length of bare tube between 2 fins
    p_f = s + t_f  # fin pitch
    theta = math.radians(theta)
    k8 = (p_f / S_l) ** 0.58
    k9 = math.log(A / A_pt)
    k10 = math.tan(theta) ** -1.5
    c5 = 0.4604 - 0.01336 * k8 * k9 * k10
    k11 = (p_f / S_t) ** 1.4
    c6 = 3.247 * k11 * k9
    c7 = -20.113 / math.log(Re_D_r)
    k16 = math.tan(theta) ** c6
    k17 = (p_f / S_l) ** c7
    k18 = math.log(A / A_pt) ** -5.35
    k19 = (D_h / D_r) ** 1.3796
    k20 = N_r ** -0.0916
    f = 0.01915 * (Re_D_r ** c5) * k16 * k17 * k18 * k19 * k20
    return f
