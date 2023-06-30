"""Friction factor correlations for cross-flow over circular finned, staggered
tube banks.

References
----------
[1]     Thulukkanam, K. (2013). Heat Exchanger Design Handbook (2nd ed.).
        CRC Press. (§ 4.5.1)
"""


def convert_friction_factor(
    f_tb: float,
    S_l: float,
    D_h: float
) -> float:
    """Convert tube bank friction factor to Fanning friction factor.

    Parameters
    ----------
    f_tb:
        Tube bank friction factor [-]
    S_l:
        Longitudinal pitch, i.e. the spacing between tubes of two adjacent tube
        rows [m].
    D_h:
        Hydraulic diameter of flow path [m]

    Returns
    -------
    Fanning friction factor
    """
    f = (D_h / S_l) * f_tb
    return f


def f_Robinson_and_Briggs_equi_triangular(
    Re_D_o: float,
    S_t: float,
    D_o: float
) -> float:
    """Robinson and Briggs correlation for the tube bank friction factor
    applicable to high-finned tubes with equilateral triangular pitch.

    Parameters
    ----------
    Re_D_o:
        Reynolds number based on tube outer diameter [-]
    S_t:
        Lateral or transverse pitch, i.e. spacing between tubes of the
        same row [m]
    D_o:
        Outside diameter of the tubes [m]

    Notes
    -----
    The correlation is applicable for:
    2000 < Re < 50_000
    1.562 < D_f < 2.750 mm
    1.687 < S_t < 4.50 mm

    Returns
    -------
    Tube bank friction factor
    """
    k1 = (S_t / D_o) ** -0.927
    f_tb = 9.465 * (Re_D_o ** -0.316) * k1
    return f_tb


def f_Robinson_and_Briggs_isosceles_triangular(
    Re_D_o: float,
    S_t: float,
    S_l: float,
    D_r: float
) -> float:
    """Robinson and Briggs correlation for the tube bank friction factor
    applicable to high-finned tubes with isosceles triangular pitch.

    Parameters
    ----------
    Re_D_o:
        Reynolds number based on tube outer diameter [-]
    S_t:
        Lateral or transverse pitch, i.e. spacing between tubes of the
        same row [m]
    S_l:
        Longitudinal pitch, i.e. spacing between tubes of two adjacent tube
        rows [m]
    D_r:
        Effective diameter of the root of a circular fin. Depending upon the
        manufacturing techniques, it may be the tube outside diameter or tube
        outside diameter plus two times the fin collar thickness [m].

    Notes
    -----

    Returns
    -------
    Tube bank friction factor
    """
    S_d = (S_t ** 2 + S_l ** 2) ** 0.5
    k1 = (S_t / D_r) ** -0.927
    k2 = (S_t / S_d) ** 0.515
    f_tb = 9.465 * (Re_D_o ** -0.316) * k1 * k2
    return f_tb


def f_Rabas_Eckels_Sabatino(
    Re_D_o: float,
    S_t: float,
    S_l: float,
    D_r: float,
    D_f: float,
    t_f: float,
    N_f: float
) -> float:
    """Rabas, Eckels and Sabatino correlation for the tube bank friction factor
    applicable to low-finned tubes with equilateral triangular pitch.

    Parameters
    ----------
    Re_D_o:
        Reynolds number based on tube outer diameter [-]
    S_t:
        Lateral or transverse pitch, i.e. spacing between tubes of the
        same row [m]
    S_l:
        Longitudinal pitch, i.e. spacing between tubes of two adjacent tube
        rows [m]
    D_r:
        Effective diameter of the root of a circular fin. Depending upon the
        manufacturing techniques, it may be the tube outside diameter or tube
        outside diameter plus two times the fin collar thickness [m].
    D_f:
        Fin tip diameter [m].
    t_f:
        Thickness of a circular fin [m].
    N_f:
        Number of fins per unit length [1/m]

    Notes
    -----

    Returns
    -------
    Tube bank friction factor
    """
    s = (1 - N_f * t_f) / N_f  # length of bare tube between 2 fins
    h_f = (D_f - D_r) / 2
    k1 = (s / D_f) ** 0.251
    k2 = (h_f / s) ** 0.759
    k3 = (D_r / D_f) ** 0.729
    k4 = (D_r / S_t) ** 0.709
    k5 = (S_t / S_l) ** 0.379
    f_tb = 3.805 * (Re_D_o ** -0.234) * k1 * k2 * k3 * k4 * k5
    return f_tb


def f_Chai(
    Re_D_o: float,
    S_t: float,
    S_l: float,
    D_r: float,
    D_f: float,
    t_f: float,
    N_f: float
) -> float:
    """Chai correlation for the tube bank friction factor.

    Parameters
    ----------
    Re_D_o:
        Reynolds number based on tube outer diameter [-]
    S_t:
        Lateral or transverse pitch, i.e. spacing between tubes of the
        same row [m]
    S_l:
        Longitudinal pitch, i.e. spacing between tubes of two adjacent tube
        rows [m]
    D_r:
        Effective diameter of the root of a circular fin. Depending upon the
        manufacturing techniques, it may be the tube outside diameter or tube
        outside diameter plus two times the fin collar thickness [m].
    D_f:
        Fin tip diameter [m].
    t_f:
        Thickness of a circular fin [m].
    N_f:
        Number of fins per unit length [1/m]

    Notes
    -----

    Returns
    -------
    Tube bank friction factor
    """
    s = (1 - N_f * t_f) / N_f  # length of bare tube between 2 fins
    h_f = (D_f - D_r) / 2
    k1 = (h_f / s) ** 0.552
    k2 = (D_r / S_t) ** 0.599
    k3 = (D_r / S_l) ** 0.1738
    f = 1.748 * (Re_D_o ** -0.233) * k1 * k2 * k3
    return f
