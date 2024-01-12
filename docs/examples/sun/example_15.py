"""
EXAMPLE 15
----------
From: Duffie, J. A., Beckman, W. A., & Blair, N. (2020).
SOLAR ENGINEERING OF THERMAL PROCESSES, PHOTOVOLTAICS AND WIND.
John Wiley & Sons.

Example 2.18.2, p. 103:
A cylindrical concentrating collector is to be oriented so that it rotates about
a horizontal east-west axis to constantly minimize the angle of incidence
and thus maximize the incident beam radiation. It is to be located at 35° N
latitude. On April 13, the day's total radiation on a horizontal surface is
22.8 MJ/m². Estimate the beam radiation on the aperture (the moving plane) of
this collector for the hour 1 to 2 p.m.
"""
from datetime import date, time
from hvac import Quantity
from hvac.sun.surface import Location, TrackingSurfaceEWC


Q_ = Quantity


location = Location(
    fi=Q_(35, 'deg'),
    date=date(2023, 4, 13),
    solar_time=time(13, 30, 0)
)

surface = TrackingSurfaceEWC(location)

I_bT = surface.estimate_I_Tb(
    H=Q_(22.8, 'MJ / m**2'),
    t1=time(13),
    t2=time(14)
)
print(I_bT.to('MJ / m**2'))
