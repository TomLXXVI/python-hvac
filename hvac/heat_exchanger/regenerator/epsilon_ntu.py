"""
Implementation of equations used in the epsilon-NTUo method for the heat
transfer analysis of rotary regenerators.

References
----------
Shah, R. K., & Sekulic, D. P. (2003). Fundamentals of Heat Exchanger Design.
John Wiley & Sons. Chapter 5: Thermal design theory for regenerators, § 5.2.
"""
import numpy as np
from hvac import Quantity

Q_ = Quantity


def UA_o(
    h_h: Quantity,
    h_c: Quantity,
    A_h: Quantity,
    A_c: Quantity,
    t_w: Quantity | None = None,
    k_w: Quantity | None = None
) -> Quantity:
    """Returns the modified overall heat transfer coefficient of the
    regenerator.

    Parameters
    ----------
    h_h:
        Heat transfer coefficient between the matrix wall and the hot fluid.
    h_c:
        Heat transfer coefficient between the matrix wall and the cold fluid.
    A_h:
        Matrix wall surface area on the hot side.
    A_c:
        Matrix wall surface area on the cold side.
    t_w: optional
        Wall thickness.
    k_w: optional
        Heat conduction transfer coefficient of the matrix material.
    """
    R_h = 1 / (h_h * A_h)
    R_c = 1 / (h_c * A_c)
    R_w = Q_(0.0, 'K / W')
    if t_w is not None and k_w is not None:
        R_w = t_w / (6 * k_w) * (1 / A_h + 1 / A_c)
    UA_o = 1 / (R_h + R_c + R_w)
    return UA_o.to('W / K')


def NTU_o(UA_o: Quantity, C_dot_min: Quantity) -> Quantity:
    """Returns the modified number of transfer units of the regenerator."""
    C_dot_min = C_dot_min.to('W / K')
    UA_o = UA_o.to('W / K')
    NTU_o = UA_o / C_dot_min
    return NTU_o.to('frac')


def C_dot_star(C_dot_min: Quantity, C_dot_max: Quantity) -> Quantity:
    C_dot_min = C_dot_min.to('W / K')
    C_dot_max = C_dot_max.to('W / K')
    C_dot_star = C_dot_min / C_dot_max
    return C_dot_star.to('frac')


def C_r_dot_star(C_r_dot: Quantity, C_dot_min: Quantity) -> Quantity:
    C_r_dot = C_r_dot.to('W / K')
    C_dot_min = C_dot_min.to('W / K')
    C_r_dot_star = C_r_dot / C_dot_min
    return C_r_dot_star.to('frac')


def hA_star(
    h_h: Quantity,
    C_dot_h: Quantity,
    theta_h: Quantity,
    h_c: Quantity,
    C_dot_c: Quantity,
    theta_c: Quantity
) -> Quantity:
    C_dot_min = min(C_dot_c.to_base_units(), C_dot_h.to_base_units())
    hA_C_dot_min = (h_c * theta_c).to_base_units()
    hA_C_dot_max = (h_h * theta_h).to_base_units()
    hA_star = hA_C_dot_min / hA_C_dot_max
    if C_dot_min == C_dot_c:
        return hA_star.to('frac')
    else:
        return (1 / hA_star).to('frac')


def _eps_KaysLondon(
    NTU_o: Quantity,
    C_dot_star: Quantity,
    C_r_dot_star: Quantity,
    C_lamda: Quantity | None = None
) -> Quantity:
    """Returns the regenerator effectiveness using an empirical correlation
    by Kays and London (1998) valid for eps <= 90 %.
    """
    if C_dot_star.magnitude == 1.0:
        eps_cf = NTU_o / (1 + NTU_o)
        # i.e. the counterflow recuperator effectiveness
    else:
        n = 1.0 - np.exp(-NTU_o * (1.0 - C_dot_star))
        d = 1.0 - C_dot_star * np.exp(-NTU_o * (1.0 - C_dot_star))
        eps_cf = n / d
    eps = eps_cf * (1 - 1 / (9.0 * (C_r_dot_star ** 1.93)))
    if C_lamda is not None:
        # correct effectiveness for longitudinal heat conduction between
        # hot and cold side of the regenerator:
        eps *= (1 - (C_lamda / (2 - C_dot_star)))
    return eps.to('frac')


def _eps_Razelos(
    NTU_o: Quantity,
    C_dot_star: Quantity,
    C_r_dot_star: Quantity,
    C_lamda: Quantity | None = None
) -> Quantity:
    """Returns the regenerator effectiveness using a procedure proposed by
    Razelos (1980) for the case `C_dot_star` < 1.
    """
    NTU_om = 2 * NTU_o * C_dot_star / (1 + C_dot_star)
    C_rm_dot_star = 2 * C_r_dot_star * C_dot_star / (1 + C_dot_star)
    eps_r = NTU_om / (1 + NTU_om)
    eps_r *= (1 - 1 / (9.0 * (C_rm_dot_star ** 1.93)))
    if C_lamda is not None:
        # correct effectiveness for longitudinal heat conduction between
        # hot and cold side of the regenerator:
        eps_r *= (1 - C_lamda)
    k = eps_r * (C_dot_star ** 2 - 1) / (2 * C_dot_star * (1 - eps_r))
    eps = (1 - np.exp(k)) / (1 - C_dot_star * np.exp(k))
    return eps.to('frac')


def eps(
    NTU_o: Quantity,
    C_dot_star: Quantity,
    C_r_dot_star: Quantity,
    C_lamda: Quantity | None = None
) -> Quantity:
    """Returns the regenerator effectiveness."""
    if C_dot_star.m > 1.0:
        eps = _eps_KaysLondon(NTU_o, C_dot_star, C_r_dot_star, C_lamda)
    else:
        eps = _eps_Razelos(NTU_o, C_dot_star, C_r_dot_star, C_lamda)
    return eps


def Q_dot_max(
    C_dot_min: Quantity,
    T_h_i: Quantity,
    T_c_i: Quantity
) -> Quantity:
    """Returns the theoretical maximum heat transfer rate in the regenerator."""
    T_h_i = T_h_i.to('K')
    T_c_i = T_c_i.to('K')
    Q_dot_max = C_dot_min * (T_h_i - T_c_i)
    return Q_dot_max.to('W')


def Q_dot(eps: Quantity, Q_dot_max: Quantity) -> Quantity:
    """Returns the heat transfer rate from the hot to the cold fluid in the
    regenerator.
    """
    Q_dot = eps * Q_dot_max
    return Q_dot.to('W')


def lamda(
    k_w: Quantity,
    A_kt: Quantity,
    L: Quantity,
    C_dot_min: Quantity
) -> Quantity:
    """Returns parameter "lambda" used to include the effect of longitudinal
    heat conduction between the hot and cold side of the regenerator.

    Parameters
    ----------
    k_w:
        Heat conduction coefficient of the matrix wall.
    A_kt:
        Total surface area for longitudinal conduction.
    L:
        Rotor length.
    C_dot_min:
        Minimum of the hot and cold fluid heat capacity rate.
    """
    lamda = (k_w * A_kt) / (L * C_dot_min)
    return lamda.to('frac')


def A_kt(A_fr: Quantity, sigma: Quantity) -> Quantity:
    """Returns the total surface area for longitudinal heat conduction from the
    hot to the cold side of the regenerator.

    Parameters
    ----------
    A_fr:
        Frontal area of wheel.
    sigma:
        Matrix porosity.
    """
    A_kt = A_fr * (1 - sigma)
    return A_kt.to('m ** 2')


def A_k_star(A_k_C_dot_min: Quantity, A_k_C_dot_max: Quantity) -> Quantity:
    A_k_star = A_k_C_dot_min / A_k_C_dot_max
    return A_k_star.to('frac')


def Fi(lamda: Quantity, NTU_o: Quantity) -> Quantity:
    a = (lamda * NTU_o / (1 + lamda * NTU_o)) ** 0.5
    b = np.tanh(NTU_o / a)
    Fi = a * b
    return Fi


def C_lamda(NTU_o: Quantity, lamda: Quantity, Fi: Quantity) -> Quantity:
    a = NTU_o * (1 + lamda * Fi) / (1 + lamda * NTU_o)
    C_lamda = 1 / (1 + a) - 1 / (1 + NTU_o)
    return C_lamda

