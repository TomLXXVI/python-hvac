"""
EXAMPLE 4
---------
From: Duffie, J. A., Beckman, W. A., & Blair, N. (2020).
SOLAR ENGINEERING OF THERMAL PROCESSES, PHOTOVOLTAICS AND WIND.
John Wiley & Sons.

Example 1.8.2, p. 28
Calculate Rb for a surface at latitude 40 °N at a tilt 30° toward the south for
the hour 9 to 10 solar time on February 16.

Example 1.8.3, p. 29
Calculate Rb for a latitude 40 °N at a tilt 50° toward the south for the hour
9 to 10 solar time on February 16.
"""
from datetime import date, time
from hvac import Quantity
from hvac.sun.surface import Location, Surface

Q_ = Quantity

location = Location(
    fi=Q_(40, 'deg'),
    date=date(2023, 2, 16)
)

surface1 = Surface(
    location=location,
    gamma=Q_(0, 'deg'),
    beta=Q_(30, 'deg')
)

surface2 = Surface(
    location=location,
    gamma=Q_(0, 'deg'),
    beta=Q_(50, 'deg')
)

for i, surface in enumerate([surface1, surface2]):
    print(f"Rb surface {i+1}: {surface.R_b_ave(time(9), time(10)):~P.2f}")
