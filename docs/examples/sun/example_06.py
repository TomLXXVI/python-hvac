"""
EXAMPLE 6
---------
From: Duffie, J. A., Beckman, W. A., & Blair, N. (2020).
SOLAR ENGINEERING OF THERMAL PROCESSES, PHOTOVOLTAICS AND WIND.
John Wiley & Sons.

Example 1.10.1, p. 38:
What is `Ho`, the day's solar radiation on a horizontal surface in the absence
of the atmosphere, at latitude 43 °N on April 15?

Example 1.10.2, p. 41:
What is the solar radiation on a horizontal surface in the absence of the
atmosphere at latitude 43 °N on April 15 between the hours of 10 and 11?
"""
from datetime import date, time
from hvac import Quantity
from hvac.sun.surface import Location

Q_ = Quantity

location = Location(
    fi=Q_(43, 'deg'),
    date=date(2023, 4, 15)
)

print(location.sun.H_o.to('MJ / m**2'))
print(location.sun.I_o(time(10), time(11)).to('MJ / m**2'))
