"""
EXAMPLE 2
---------
From: Duffie, J. A., Beckman, W. A., & Blair, N. (2020).
SOLAR ENGINEERING OF THERMAL PROCESSES, PHOTOVOLTAICS AND WIND.
John Wiley & Sons.

Example 1.7.1, p. 23:
Calculate the angle of incidence of beam radiation, the slope of the surface,
and the surface azimuth angle for a surface at latitude 40° on July 17
a.  at hour angle 30°
b.  at hour angle 100°
if it is continuously rotated about an east-west axis to minimize the incidence
angle.
"""
from datetime import date
from hvac import Quantity
from hvac.sun.surface import Location, TrackingSurfaceEWC
from hvac.sun.time import hour_angle_to_solar_time, decimal_hour_to_time


Q_ = Quantity

omega_lst = [Q_(30, 'deg'), Q_(100, 'deg')]

for omega in omega_lst:
    location = Location(
        fi=Q_(40, 'deg'),
        date=date(2023, 7, 17),
        solar_time=decimal_hour_to_time(hour_angle_to_solar_time(omega.to('rad').m))
    )
    tracking_surface = TrackingSurfaceEWC(location)
    print(
        f"incidence angle: {tracking_surface.theta.to('deg')}",
        f"slope angle: {tracking_surface.beta.to('deg')}",
        f"surface azimuth angle: {tracking_surface.gamma.to('deg')}",
        sep='\n', end='\n\n'
    )
