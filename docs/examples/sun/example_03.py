"""
EXAMPLE 3
---------
From: Duffie, J. A., Beckman, W. A., & Blair, N. (2020).
SOLAR ENGINEERING OF THERMAL PROCESSES, PHOTOVOLTAICS AND WIND.
John Wiley & Sons.

Example 1.6.1, p. 15:
Calculate the angle of incidence of beam radiation on a surface located at
Madison, Wisconsin, at 10:30 (solar time) on February 13 if the surface is
tilted 45° from the horizontal and pointed 15° west of South.

Example 1.8.1, p. 24:
What is the ratio of beam radiation to that on a horizontal surface for the
surface and time specified in example 1.6.1?
"""
from datetime import date, time
from hvac import Quantity
from hvac.sun.surface import Location, Surface

Q_ = Quantity

location = Location(
    fi=Q_(43, 'deg'),
    date=date(2023, 2, 13),
    solar_time=time(10, 30)
)

surface = Surface(
    location=location,
    gamma=Q_(15, 'deg'),
    beta=Q_(45, 'deg')
)

print(
    f"angle of incidence: {surface.theta.to('deg'):~P.1f}",
    f"ratio of beam radiation: {surface.R_b:~P.2f}",
    sep='\n'
)
