"""
EXAMPLE 5
---------
From: Duffie, J. A., Beckman, W. A., & Blair, N. (2020).
SOLAR ENGINEERING OF THERMAL PROCESSES, PHOTOVOLTAICS AND WIND.
John Wiley & Sons.

Example 1.9.2, p. 34:
It is proposed to install a solar collector at a level 4.0 m above the ground.
A rectangular building 30 m high is located 45 m to the south, has its long
dimension on an east-west axis, and has dimensions 52 + 8 m x 18 m. The latitude
is 45°. Diagram this building on the solar position plot to show the times of
day and year when it would shade the proposed collector.
"""
from hvac import Quantity
import numpy as np
from hvac.sun.sun_chart import SunChart, ShadingPoint, ShadingProfile, Location

Q_ = Quantity

# point A:
h_A = Q_(30 - 4, 'm')
d_A = Q_(45, 'm')
alpha_A = np.arctan(h_A / d_A)
gamma_A = Q_(0, 'deg')

# point B:
h_B = h_A
d_B = np.sqrt(d_A**2 + Q_(52, 'm')**2)
alpha_B = np.arctan(h_B / d_B)
gamma_B = np.arctan(Q_(52, 'm') / d_A)  # angle to the west: +

# point C:
h_C = h_A
d_C = np.sqrt(d_A ** 2 + Q_(8, 'm') ** 2)
alpha_C = np.arctan(h_C / d_C)
gamma_C = -np.arctan(Q_(8, 'm') / d_A)  # angle to the east: -

sun_chart = SunChart(Location(fi=Q_(45, 'deg')))
sun_chart.add_shading_profile(ShadingProfile(
    ID='building',
    pts=[
        ShadingPoint(gamma_A, alpha_A, d_A),
        ShadingPoint(gamma_B, alpha_B, d_B),
        ShadingPoint(gamma_C, alpha_C, d_C)
    ]
))
chart = sun_chart.plot()
chart.show()


