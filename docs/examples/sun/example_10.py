"""
EXAMPLE 10
----------
From: Duffie, J. A., Beckman, W. A., & Blair, N. (2020).
SOLAR ENGINEERING OF THERMAL PROCESSES, PHOTOVOLTAICS AND WIND.
John Wiley & Sons.

Example 2.13.1, p. 84:
What is the fraction of the average January daily radiation that is received at
Melbourne, Australia, in the hour between 8:00 and 9:00?
"""
from datetime import time
from hvac import Quantity
from hvac.sun.radiation import ReferenceDates, estimate_r_t
from hvac.sun.surface import Location

Q_ = Quantity


day = ReferenceDates.get_date_for('Jan')
melbourne = Location(
    fi=Q_(-38, 'deg'),  # southern hemisphere
    date=day,
    solar_time=time(8, 30, 0)
)

r_t = estimate_r_t(
    melbourne.omega.to('rad').m,
    melbourne.omega_ss.to('rad').m
)
print(r_t)
