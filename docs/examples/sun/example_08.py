"""
EXAMPLE 8
---------
From: Duffie, J. A., Beckman, W. A., & Blair, N. (2020).
SOLAR ENGINEERING OF THERMAL PROCESSES, PHOTOVOLTAICS AND WIND.
John Wiley & Sons.

Example 2.11.1, p.80:
The day's total radiation on a horizontal surface for St. Louis, Missouri
(latitude 38.6°) on September 3 is 23 MJ/m². Estimate the fraction and amount
that is diffuse.
"""
from datetime import date
from hvac import Quantity
from hvac.sun.surface import Location
from hvac.sun.radiation import estimate_H_d_fraction

Q_ = Quantity

H = Q_(23.0, 'MJ / m**2')

location = Location(
    fi=Q_(38.6, 'deg'),
    date=date(2023, 9, 3),
)

K_T = H / location.sun.H_o
H_d_frac = estimate_H_d_fraction(
    K_T=K_T.to('frac').m,
    omega_ss=location.omega_ss.to('rad').m
)
H_d = H_d_frac * H
print(
    f"Amount of diffuse radiation: {H_d.to('MJ / m**2'):~P.1f}"
)
