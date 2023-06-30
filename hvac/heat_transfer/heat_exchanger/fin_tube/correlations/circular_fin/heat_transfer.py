"""Heat transfer correlations for cross-flow over circular finned, staggered
tube banks.

References
----------
[1]     Thulukkanam, K. (2013). Heat Exchanger Design Handbook (2nd ed.).
        CRC Press. (§ 4.5.1)
[2]     Shah, R. K., & Sekulic, D. P. (2003). Fundamentals of Heat Exchanger
        Design. John Wiley & Sons. (§ 7.5.4)
"""


def Nu_Briggs_and_Young_low_fin(
    Re_D_o: float,
    Pr: float,
    N_f: float,
    t_f: float,
    h_f: float
) -> float:
    """Briggs and Young correlation for circular low-fin staggered tube bank.

    Parameters
    ----------
    Re_D_o:
        Reynolds number based on tube outer diameter [-]
    Pr:
        Prandtl number [-]
    N_f:
        Number of fins per unit length [1/m]
    t_f:
        Thickness of a circular fin [m]
    h_f:
        Fin height [m]

    Notes
    -----
    The correlation is applicable for a circular finned tube with low fin height
    and high fin density, 1000 < Re_D_o < 20_000, N_r >= 6 and equilateral
    triangular pitch.

    Returns
    -------
    Nusselt number
    """
    s = (1 - N_f * t_f) / N_f  # length of bare tube between 2 fins
    k1 = (s / h_f) ** 0.164
    k2 = (s / t_f) ** 0.075
    Nu = 0.1507 * Re_D_o ** 0.667 * Pr ** (1 / 3) * k1 * k2
    return Nu


def Nu_Briggs_and_Young_high_fin(
    Re_D_o: float,
    Pr: float,
    N_f: float,
    t_f: float,
    h_f: float
) -> float:
    """Briggs and Young correlation for circular high-fin staggered tube bank.

    Parameters
    ----------
    Re_D_o:
        Reynolds number based on tube outer diameter [-]
    Pr:
        Prandtl number [-]
    N_f:
        Number of fins per unit length [1/m]
    t_f:
        Thickness of a circular fin [m]
    h_f:
        Fin height [m]

    Returns
    -------
    Nusselt number
    """
    s = (1 - N_f * t_f) / N_f  # length of bare tube between 2 fins
    k = (s / h_f) ** 0.296
    Nu = 0.1378 * Re_D_o ** 0.718 * Pr ** (1 / 3) * k
    return Nu


def Nu_Briggs_and_Young(
    Re_D_o: float,
    Pr: float,
    N_f: float,
    t_f: float,
    h_f: float
) -> float:
    """Briggs and Young correlation for circular finned staggered tube
    bank.

    Parameters
    ----------
    Re_D_o:
        Reynolds number based on tube outer diameter [-]
    Pr:
        Prandtl number [-]
    N_f:
        Number of fins per unit length [1/m]
    t_f:
        Thickness of a circular fin [m]
    h_f:
        Fin height [m]

    Returns
    -------
    Nusselt number
    """
    s = (1 - N_f * t_f) / N_f  # length of bare tube between 2 fins
    k1 = (s / h_f) ** 0.200
    k2 = (s / t_f) ** 0.1134
    Nu = 0.134 * Re_D_o ** 0.681 * Pr ** (1 / 3) * k1 * k2
    return Nu


def j_Briggs_and_Young(
    Re_D_o: float,
    N_f: float,
    t_f: float,
    h_f: float
) -> float:
    """Briggs and Young correlation for circular finned staggered tube
    bank.

    Parameters
    ----------
    Re_D_o:
        Reynolds number based on tube outer diameter [-]
    N_f:
        Number of fins per unit length [1/m]
    t_f:
        Thickness of a circular fin [m]
    h_f:
        Fin height [m]

    Returns
    -------
    Colburn j factor
    """
    s = (1 - N_f * t_f) / N_f  # length of bare tube between 2 fins
    k1 = (s / h_f) ** 0.2
    k2 = (s / t_f) ** 0.11
    j = 0.134 * Re_D_o ** -0.319 * k1 * k2
    return j
