"""
EXAMPLE 1
---------
From: Duffie, J. A., Beckman, W. A., & Blair, N. (2020).
SOLAR ENGINEERING OF THERMAL PROCESSES, PHOTOVOLTAICS AND WIND.
John Wiley & Sons.

Example 1.6.3, p. 19:
Calculate the time of sunrise, solar altitude, zenith, solar azimuth, and
profile angles for a 60° sloped surface facing 25° west of south at 4:00 p.m.
solar time on March 16 at a latitude of 43°. Also calculate the time of sunrise
and sunset on the surface.
"""
from datetime import date, time
from hvac import Quantity
from hvac.sun.surface import Location, Surface

Q_ = Quantity

location = Location(fi=Q_(43, 'deg'))  # fi: latitude
location.date = date(2023, 3, 16)
location.solar_time = time(16, 0, 0)

surface = Surface(
    location=location,
    gamma=Q_(25, 'deg'),  # surface azimuth angle (S is 0°, E is -90°, W is +90°)
    beta=Q_(60, 'deg')    # surface slope angle
)

print(
    f"solar time of sunrise at location: {location.sunrise}",
    f"solar time of sunset at location: {location.sunset}",
    f"solar time of sunrise on surface: {surface.sunrise}",
    f"solar time of sunset on surface: {surface.sunset}",
    f"solar altitude angle at given date and solar time: {location.sun.alpha.to('deg')}",
    f"zenith angle at given date and solar time: {location.sun.theta.to('deg')}",
    f"solar azimuth angle at given date and solar time: {location.sun.gamma.to('deg')}",
    f"profile angle of surface: {surface.alpha_p.to('deg')}",
    sep='\n'
)
