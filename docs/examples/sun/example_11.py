"""
EXAMPLE 11
----------
From: Duffie, J. A., Beckman, W. A., & Blair, N. (2020).
SOLAR ENGINEERING OF THERMAL PROCESSES, PHOTOVOLTAICS AND WIND.
John Wiley & Sons.

Example 2.13.2, p. 85:
The total radiation for Madison on August 23 was 31.4 MJ/m². Estimate the
radiation received between 1 and 2 p.m.
"""
from datetime import time, date
from hvac import Quantity
from hvac.sun.radiation import estimate_r_t
from hvac.sun.surface import Location

Q_ = Quantity

location = Location(
    fi=Q_(43, 'deg'),
    date=date(2023, 8, 23),
    solar_time=time(13, 30, 0)
)
r_t = estimate_r_t(
    location.omega.to('rad').m,
    location.omega_ss.to('rad').m
)
H = Q_(31.4, 'MJ / m ** 2')
I = r_t * H
print(I)
