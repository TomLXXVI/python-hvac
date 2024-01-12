"""
EXAMPLE 13
----------
From: Duffie, J. A., Beckman, W. A., & Blair, N. (2020).
SOLAR ENGINEERING OF THERMAL PROCESSES, PHOTOVOLTAICS AND WIND.
John Wiley & Sons.

Example 2.14.1, p. 90:
On March 4 at a latitude of 45° and a surface slope of 60°, determine Rb at
6:30 a.m. and Rb,ave for the hour 6 to 7 a.m.
"""
from datetime import date, time
from hvac import Quantity
from hvac.sun.surface import Location, Surface

Q_ = Quantity

location = Location(
    fi=Q_(45, 'deg'),
    date=date(2023, 3, 4),
    solar_time=time(6, 30, 0)
)
print(f"sunrise at location: {location.sunrise}")

surface = Surface(
    location,
    gamma=Q_(0, 'deg'),
    beta=Q_(60, 'deg')
)
print(f"sunrise at surface: {surface.sunrise}")

print(
    f"Rb: {surface.R_b.to('frac'):~P.2f}",
    f"Rb,ave: {surface.R_b_ave(location.sunrise, time(7)).to('frac'):~P.2f}"
)

# Around sunrise (and sunset) the Rb value is unreasonably high. A more
# reasonable value can be attained by determining the average Rb value between
# the time of sunrise and the next whole hour, or to neglect the hours that
# contain sunrise or sunset.
