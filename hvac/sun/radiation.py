import typing
from datetime import date as Date
import pandas as pd
import numpy as np
from hvac.sun.time import day_number
from hvac.sun.geometry import declination, daylight_time

G_sc = 1367.0  # W / m ** 2


class ClimateType:
    """
    Correction factors for atmospheric transmittance for beam radiation `tau_b`
    depending on climate type.

    References
    ----------
    Duffie, J. A., Beckman, W. A., & Blair, N. (2020).
    SOLAR ENGINEERING OF THERMAL PROCESSES, PHOTOVOLTAICS AND WIND.
    John Wiley & Sons., p. 70, table 2.8.1.
    """
    TROPICAL = 'tropical'
    MID_LATITUDE_SUMMER = 'mid-latitude summer'
    SUBARCTIC_SUMMER = 'subarctic summer'
    MID_LATITUDE_WINTER = 'mid-latitude winter'
    _corr_factors = pd.DataFrame(
        data=[
            [0.95, 0.98, 1.02],
            [0.97, 0.99, 1.02],
            [0.99, 0.99, 1.01],
            [1.03, 1.01, 1.00]
        ],
        index=[
            TROPICAL,
            MID_LATITUDE_SUMMER,
            SUBARCTIC_SUMMER,
            MID_LATITUDE_WINTER
        ],
        columns=['r0', 'r1', 'rk']
    )

    @classmethod
    def correction_factors(cls, type_: str) -> tuple[float, ...]:
        return typing.cast(tuple, (
            cls._corr_factors.loc[type_, 'r0'],
            cls._corr_factors.loc[type_, 'r1'],
            cls._corr_factors.loc[type_, 'rk']
        ))


class ReferenceDates:
    _dates: dict[str, Date] = {
        'jan': Date(2023, 1, 17),
        'feb': Date(2023, 2, 16),
        'mar': Date(2023, 3, 16),
        'apr': Date(2023, 4, 15),
        'may': Date(2023, 5, 15),
        'jun': Date(2023, 6, 11),
        'jul': Date(2023, 7, 17),
        'aug': Date(2023, 8, 16),
        'sep': Date(2023, 9, 15),
        'oct': Date(2023, 10, 15),
        'nov': Date(2023, 11, 14),
        'dec': Date(2023, 12, 10)
    }

    @classmethod
    def get_date_for(cls, m: str | int) -> Date:
        if isinstance(m, int) and 1 <= m <= 12:
            month = list(cls._dates.keys())[m - 1]
        elif isinstance(m, str):
            for month in cls._dates.keys():
                if month in m.lower():
                    month = month
                    break
            else:
                raise ValueError(f'month {m} does not exist')
        else:
            raise ValueError(f'month {m} does not exist')
        return cls._dates[month]

    @classmethod
    def months(cls):
        for m in cls._dates.keys():
            yield m


def air_mass(theta_z: float | np.ndarray) -> float | np.ndarray:
    """Close approximation for air Mass at sea level when zenith angle is
    between 0° to 70°.

    Air mass is the ratio of the mass of atmosphere trough which beam radiation
    passes to the mass it would pass through if the sun were directly overhead
    (i.e. at the zenith).
    """
    return 1 / np.cos(theta_z)


def ext_irradiance_normal(n: float | np.ndarray) -> float | np.ndarray:
    """Extraterrestrial irradiance incident on the plane normal to the
    radiation (beam radiation).
    """
    G_on = G_sc * (1.0 + 0.033 * np.cos(np.radians(360 * n / 365)))
    return G_on


def ext_irradiance_hor_surf(
    n: float | np.ndarray,
    theta_z: float | np.ndarray
) -> float | np.ndarray:
    """Extraterrestrial irradiance incident on a horizontal plane (outside
    the atmosphere).
    """
    G_o = ext_irradiance_normal(n) * np.cos(theta_z)
    return G_o


def clear_sky_beam_transmittance(
    A: float | np.ndarray,
    climate_type: str,
    theta_z: float | np.ndarray
) -> float | np.ndarray:
    """Atmospheric transmittance for clear-sky beam radiation.

    This is the ratio of beam radiation to the extraterrestrial beam
    radiation on a horizontal or tilted surface.
    """
    a0_star = 0.4237 - 0.00821 * (6 - A) ** 2
    a1_star = 0.5055 + 0.00595 * (6.5 - A) ** 2
    k_star = 0.2711 + 0.01858 * (2.5 - A) ** 2
    r0, r1, rk = ClimateType.correction_factors(climate_type)
    a0 = r0 * a0_star
    a1 = r1 * a1_star
    k = rk * k_star
    try:
        tau_b = a0 + a1 * np.exp(-k / np.cos(theta_z))
    except FloatingPointError:
        return a0
    else:
        return tau_b


def clear_sky_diffuse_transmittance(tau_b: float | np.ndarray) -> float | np.ndarray:
    """Atmospheric transmittance for clear-sky diffuse radiation.

    This is the ratio of diffuse radiation `G_cd` to the extraterrestrial
    beam radiation on the horizontal plane `G_o`.
    """
    return 0.271 - 0.294 * tau_b


def estimate_avg_K_T(
    fi: float,
    month: str | int,
    a: float,
    b: float,
    sunshine_hours: float
) -> float:
    """Get estimation of monthly average daily clearness index at the given
    latitude from hours of sunshine for the given month.

    The monthly average daily clearness index (`K_T_avg`) is the ratio of
    monthly average daily total radiation (`_H_avg`) to the extraterrestrial
    monthly average daily radiation (`H_o_avg`) on a horizontal surface:
    ```
    K_T_avg = _H_avg / H_o_avg
    ```

    Parameters
    ----------
    fi:
        latitude
    month:
        month of the year
    a, b:
        constants depending on location. See "SOLAR ENGINEERING OF THERMAL
        PROCESSES, PHOTOVOLTAICS AND WIND" (5th. Ed.), p. 69, table 2.7.2
    sunshine_hours:
        monthly average daily hours of bright sunshine.
        See "SOLAR ENGINEERING OF THERMAL PROCESSES, PHOTOVOLTAICS AND WIND"
        (5th. Ed.), p. 66, table 2.7.1.
    """
    date = ReferenceDates.get_date_for(month)
    delta = declination(day_number(date))
    N = daylight_time(fi, delta)
    K_T_avg = a + b * sunshine_hours / N
    return K_T_avg


def cumulative_K_T_histogram(
    K_T: float | np.ndarray,
    K_T_avg: float
) -> float | np.ndarray:
    """Given the monthly average clearness index `K_T_avg` returns the fraction
    of days in the given month that are less clear than the given daily
    clearness index `K_T`.

    See "SOLAR ENGINEERING OF THERMAL PROCESSES, PHOTOVOLTAICS AND WIND"
    (5th. Ed.), p. 76, Eqs. 2.9.4, 2.9.5 and 2.9.6.
    """
    K_T_min = 0.05
    K_T_max = 0.6313 + 0.267 * K_T_avg - 11.9 * (K_T_avg - 0.75) ** 8
    zeta = (K_T_max - K_T_min) / (K_T_max - K_T_avg)
    gamma = (
        -1.498 + (1.184 * zeta - 27.182 * np.exp(-1.5 * zeta))
        / (K_T_max - K_T_min)
    )
    f = (
        (np.exp(gamma * K_T_min) - np.exp(gamma * K_T))
        / (np.exp(gamma * K_T_min) - np.exp(gamma * K_T_max))
    )
    return f


def estimate_I_d_fraction(k_T: float) -> float:
    """Get estimation of the fraction of the hourly total radiation on a
    horizontal plane which is diffuse (`I_d / I`) according to correlation of
    Erbs et al. (1982).

    See "SOLAR ENGINEERING OF THERMAL PROCESSES, PHOTOVOLTAICS AND WIND"
    (5th. Ed.), p. 78, Eq. 2.10.1

    Parameters
    ----------
    k_T:
        hourly clearness index (ratio of hourly total radiation to hourly
        extraterrestrial radiation on a horizontal surface, `I / I_o`)
    """
    if k_T <= 0.22:
        frac_I_d = 1.0 - 0.09 * k_T
    elif 0.22 < k_T <= 0.8:
        frac_I_d = (
            0.9511 - 0.1604 * k_T + 4.388 * (k_T ** 2)
            - 16.638 * (k_T ** 3) + 12.336 * (k_T ** 4)
        )
    else:
        frac_I_d = 0.165
    return frac_I_d


def estimate_H_d_fraction(K_T: float, omega_ss: float) -> float:
    """Get estimation of the fraction of the daily total radiation on a
    horizontal plane which is diffuse (`H_d / H`) according to correlation of
    Erbs et al. (1982).

    See "SOLAR ENGINEERING OF THERMAL PROCESSES, PHOTOVOLTAICS AND WIND"
    (5th. Ed.), p. 80, Eq. 2.11.1

    Parameters
    ----------
    K_T:
        daily clearness index (`H / H_o`)
    omega_ss:
        sunset hour angle in radians
    """
    if omega_ss <= np.radians(81.4):
        if K_T < 0.715:
            frac_H_d = (
                1.0 - 0.2727 * K_T + 2.4495 * (K_T ** 2)
                - 11.9514 * (K_T ** 3) + 9.3879 * (K_T ** 4)
            )
        else:
            frac_H_d = 0.143
    else:
        if K_T < 0.722:
            frac_H_d = (
                1.0 + 0.2832 * K_T - 2.5557 * (K_T ** 2)
                + 0.8448 * (K_T ** 3)
            )
        else:
            frac_H_d = 0.175
    return min(frac_H_d, 1.0)


def estimate_avg_H_d_fraction(
    K_T_avg: float,
    omega_ss: float
) -> float:
    """Get estimation of the fraction of the monthly average daily total
    radiation on a horizontal plane (`H_d_avg / _H_avg`) which is diffuse
    according to correlation of Erbs et al. (1982).

    See "SOLAR ENGINEERING OF THERMAL PROCESSES, PHOTOVOLTAICS AND WIND"
    (5th. Ed.), p. 82, Eq. 2.12.1

    Parameters
    ----------
    K_T_avg:
        monthly average daily clearness index
    omega_ss
        sunset hour angle in radians
    """
    if 0.3 <= K_T_avg <= 0.8:
        if omega_ss <= np.radians(81.4):
            frac_H_d_avg = (
                1.391 - 3.560 * K_T_avg + 4.189 * (K_T_avg ** 2)
                - 2.137 * (K_T_avg ** 3)
            )
        else:
            frac_H_d_avg = (
                1.311 - 3.022 * K_T_avg + 3.427 * (K_T_avg ** 2)
                - 1.821 * (K_T_avg ** 3)
            )
        return frac_H_d_avg
    else:
        raise ValueError('parameter `K_T_avg` must be between 0.3 and 0.8')


def estimate_r_t(omega: float | np.ndarray, omega_ss: float) -> float | np.ndarray:
    """Get estimation of the ratio of hourly total radiation to daily total
    radiation on a horizontal surface (`I / H` ) as function of day length
    (sunset hour angle) and the hour (hour angle) in question (i.e. the fraction
    of daily total radiation received on horizontal surface within a given hour
    of the day).

    See "SOLAR ENGINEERING OF THERMAL PROCESSES, PHOTOVOLTAICS AND WIND"
    (5th. Ed.), p. 83, Eqs. 2.13.2

    Parameters
    ----------
    omega:
        hour angle in radians
    omega_ss:
        sunset hour angle in radians
    """
    a = 0.409 + 0.5016 * np.sin(omega_ss - np.radians(60))
    b = 0.6609 - 0.4767 * np.sin(omega_ss - np.radians(60))
    r_t = (
        np.pi / 24 * (a + b * np.cos(omega))
        * (np.cos(omega) - np.cos(omega_ss))
        / (np.sin(omega_ss) - omega_ss * np.cos(omega_ss))
    )
    return r_t


def estimate_r_d(omega: float | np.ndarray, omega_ss: float) -> float | np.ndarray:
    """Get estimation of the ratio of hourly diffuse radiation to daily diffuse
    radiation on a horizontal surface (`I_d / H_d` ) as function of day length
    (sunset hour angle) and the hour (hour angle) in question (i.e. the fraction
    of daily diffuse radiation received on horizontal surface within a given
    hour of the day).

    See "SOLAR ENGINEERING OF THERMAL PROCESSES, PHOTOVOLTAICS AND WIND"
    (5th. Ed.), p. 82, Eq. 2.13.4

    Parameters
    ----------
    omega:
        hour angle in radians
    omega_ss:
        sunset hour angle in radians
    """
    r_d = (
        np.pi / 24
        * (np.cos(omega) - np.cos(omega_ss))
        / (np.sin(omega_ss) - omega_ss * np.cos(omega_ss))
    )
    return r_d
