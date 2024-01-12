# Source: "SOLAR ENGINEERING OF THERMAL PROCESSES, PHOTOVOLTAICS AND WIND"
# (5th. Ed.) by John A. Duffie and William A. Beckman, §1.6 (p. 12).
import numpy as np


# instead of warnings, we want an exception
np.seterr(all='raise')


def _B(n: float | np.ndarray) -> float | np.ndarray:
    return np.radians((n - 1) * (360 / 365))


def declination(n: float | np.ndarray) -> float | np.ndarray:
    """Returns the declination angle (delta) of the sun for day of the year `n`
    in radians.
    """
    # more accurate equation from Spencer (1971):
    B = _B(n)
    delta_ = (
        0.006918 - 0.399912 * np.cos(B) + 0.070257 * np.sin(B)
        - 0.006758 * np.cos(2 * B) + 0.000907 * np.sin(2 * B)
        - 0.002697 * np.cos(3 * B) + 0.00148 * np.sin(3 * B)
    )
    return delta_
    # approximate equation of Cooper (1969):
    # return np.radians(23.45) * np.sin(2 * np.pi * (284 + n) / 365)


def hour_angle(solar_time: float | np.ndarray) -> float | np.ndarray:
    """Returns the angular displacement of the sun east or west of the local
    meridian due to rotation of the earth on its axis at 15° per hour; morning
    negative, afternoon positive.

    Parameters
    ----------
    solar_time:
        solar time expressed in decimals hours

    Returns
    -------
    hour angle in radians
    """
    omega = np.radians(15 * (solar_time - 12))
    return omega


def _theta_01(
    delta: float | np.ndarray,
    fi: float | np.ndarray,
    beta: float | np.ndarray,
    gamma: float | np.ndarray,
    omega: float | np.ndarray
) -> float | np.ndarray:
    """Returns the angle of incidence of beam radiation on a surface.

    Parameters
    ----------
    delta: radians
        declination angle of the sun
    fi: radians
        latitude; north positive; -pi/2 <= fi <= pi/2
    beta: radians
        slope angle between the plane of the surface and the horizontal;
        0 <= beta <= pi
    gamma: radians
        surface azimuth angle, with zero due south, east negative, and west
        positive; -pi <= gamma <= pi.
    omega: radians
        hour angle; morning negative, afternoon positive

    Returns
    -------
    incidence angle in radians
    """
    theta_ = np.arccos(
        np.sin(delta) * np.sin(fi) * np.cos(beta)
        - np.sin(delta) * np.cos(fi) * np.sin(beta) * np.cos(gamma)
        + np.cos(delta) * np.cos(fi) * np.cos(beta) * np.cos(omega)
        + np.cos(delta) * np.sin(fi) * np.sin(beta) * np.cos(gamma) * np.cos(omega)
        + np.cos(delta) * np.sin(beta) * np.sin(gamma) * np.sin(omega)
    )
    return theta_


def _theta_02(
    theta_z: float | np.ndarray,
    beta: float | np.ndarray,
    gamma_s: float | np.ndarray,
    gamma: float | np.ndarray
) -> float | np.ndarray:
    """Returns the angle of incidence of beam radiation on a surface.

    Parameters
    ----------
    theta_z: radians
        zenith angle
    beta: radians
        slope angle between the plane of the surface and the horizontal;
        0 <= beta <= pi
    gamma_s: radians
        solar azimuth angle; east of south negative, west of south positive
    gamma: radians
        surface azimuth angle, with zero due south, east negative, and west
        positive; -pi <= gamma <= pi.

    Returns
    -------
    incidence angle in radians
    """
    theta_ = np.arccos(
        np.cos(theta_z) * np.cos(beta)
        + np.sin(theta_z) * np.sin(beta) * np.cos(gamma_s - gamma)
    )
    return theta_


def incidence_angle(
    beta: float | np.ndarray,
    gamma: float | np.ndarray,
    delta: float | np.ndarray | None = None,
    fi: float | np.ndarray | None = None,
    omega: float | np.ndarray | None = None,
    theta_z: float | np.ndarray | None = None,
    gamma_s: float | np.ndarray | None = None
) -> float | np.ndarray:
    """Returns the angle of incidence of beam radiation on a surface.

    Parameters
    ----------
    beta: radians
        slope angle between the plane of the surface and the horizontal;
        0 <= beta <= pi
    gamma: radians
        surface azimuth angle, with zero due south, east negative, and west
        positive; -pi <= gamma <= pi.
    delta: radians, default None
        declination angle of the sun
    fi: radians, default None
        latitude; north positive; -pi/2 <= fi <= pi/2
    omega: radians, default None
        hour angle; morning negative, afternoon positive
    theta_z: radians, default None
        zenith angle
    gamma_s: radians, default None
        solar azimuth angle; east of south negative, west of south positive

    Notes
    -----
    Parameters `beta` and `gamma` are always required. Depending on which other
    parameters are given, another formula will be used internally to calculate
    the incidence angle. Two sets of other parameters are possible: (1) either
    also `delta`, `fi`, and `omega` are specified, (2) or also `theta_z` and
    `gamma_s` are specified.

    If the incidence angle is greater than 90°, it means that the sun is behind
    the surface. Also, it is necessary to ensure that the earth is not blocking
    the sun, i.e. that the hour angle `omega` is between sunrise and sunset.

    Returns
    -------
    incidence angle in radians
    """
    set_1 = [
        True if elem is not None else False
        for elem in [delta, fi, beta, gamma, omega]
    ]
    set_2 = [
        True if elem is not None else False
        for elem in [theta_z, beta, gamma_s, gamma]
    ]
    if all(set_1):
        return _theta_01(delta, fi, beta, gamma, omega)
    elif all(set_2):
        return _theta_02(theta_z, beta, gamma_s, gamma)
    else:
        raise ValueError('incorrect function call')


def zenith_angle(
    delta: float | np.ndarray,
    fi: float | np.ndarray,
    omega: float | np.ndarray
) -> float | np.ndarray:
    """Returns the zenith angle of the sun.

    Parameters
    ----------
    delta: radians, default None
        declination angle of the sun
    fi: radians, default None
        latitude; north positive; -pi/2 <= fi <= pi/2
    omega: radians, default None
        hour angle; morning negative, afternoon positive

    Returns
    -------
    zenith angle in radians
    """
    return incidence_angle(0.0, 0.0, delta, fi, omega, None, None)


def solar_altitude_angle(
    delta: float | np.ndarray,
    fi: float | np.ndarray,
    omega: float | np.ndarray
) -> float | np.ndarray:
    """Returns the altitude angle of the sun.

    Parameters
    ----------
    delta: radians, default None
        declination angle of the sun
    fi: radians, default None
        latitude; north positive; -pi/2 <= fi <= pi/2
    omega: radians, default None
        hour angle; morning negative, afternoon positive

    Returns
    -------
    solar altitude angle in radians
    """
    theta_z_ = zenith_angle(delta, fi, omega)
    alpha_s_ = np.pi / 2 - theta_z_
    return alpha_s_


def solar_azimuth_angle(
    omega: float | np.ndarray,
    theta_z: float | np.ndarray,
    fi: float | np.ndarray,
    delta: float | np.ndarray
) -> float | np.ndarray:
    """Returns the solar azimuth angle.

    Parameters
    ----------
    omega: radians
        hour angle; morning negative, afternoon positive
    theta_z: radians
        zenith angle of sun
    fi: radians
        latitude; north positive; -pi/2 <= fi <= pi/2
    delta: radians
        declination angle of the sun

    Returns
    -------
    solar azimuth angle in radians
    """
    try:
        gamma_s_ = np.sign(omega) * np.abs(np.arccos(
            (np.cos(theta_z) * np.sin(fi) - np.sin(delta))
            / (np.sin(theta_z) * np.cos(fi))
        ))
    except FloatingPointError:
        return 0.0
    else:
        return gamma_s_


def sunset_hour_angle(
    fi: float | np.ndarray,
    delta: float | np.ndarray
) -> float | np.ndarray:
    """Returns the sunset hour angle.

    Parameters
    ----------
    fi: radians
        latitude; north positive; -pi/2 <= fi <= pi/2
    delta: radians
        declination angle of the sun

    Notes
    -----
    The sunrise hour angle is the negative of the sunset hour angle.

    Returns
    -------
    sunset hour angle in radians
    """
    omega_s_ = np.arccos(-np.tan(fi) * np.tan(delta))
    return omega_s_


def daylight_time(
    fi: float | np.ndarray,
    delta: float | np.ndarray
) -> float | np.ndarray:
    """Returns the number of daylight hours.

    Parameters
    ----------
    fi: radians
        latitude; north positive; -pi/2 <= fi <= pi/2
    delta: radians
        declination angle of the sun

    Returns
    -------
    number of daylight hours
    """
    omega_s_ = sunset_hour_angle(fi, delta)
    N = (2 / 15) * np.degrees(omega_s_)
    return N


def profile_angle(
    alpha_s: float | np.ndarray,
    gamma_s: float | np.ndarray,
    gamma: float | np.ndarray
) -> float | np.ndarray:
    """Returns the profile angle of beam radiation on a receiver plane.

    Parameters
    ----------
    alpha_s: radians
        the solar altitude angle, which is the complement of the zenith angle
    gamma_s: radians
        the solar azimuth angle
    gamma: radians
        surface azimuth angle, with zero due south, east negative, and west
        positive; -pi <= gamma <= pi.

    Returns
    -------
    the profile angle in radians
    """
    alpha_p_ = np.arctan(
        np.tan(alpha_s) / np.cos(gamma_s - gamma)
    )
    return alpha_p_
