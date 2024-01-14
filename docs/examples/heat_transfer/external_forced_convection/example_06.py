"""
EXAMPLE 6
---------
From: Nellis G. F., & Klein S. A. (2021)
INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.

Example 8.6
An uninsulated pipe is carrying water at an average temperature of T_w = 160 °F
outdoors on a day when the wind is blowing with a velocity u_inf = 5 mph and
temperature T_inf = -5 °F. The pipe has an outer diameter of D = 0.5 inch and
a length of L = 10 ft. Estimate the rate of heat loss from the pipe.
"""
from math import pi as PI
from hvac import Quantity
from hvac.fluids import Fluid
from hvac.heat_transfer.forced_convection.external_flow.cylinder import Cylinder

Q_ = Quantity


Air = Fluid('Air')
T_inf = Q_(-5, 'degF').to('K')
T_water = Q_(160, 'degF').to('K')
T_flm = (T_water + T_inf) / 2
air = Air(T=T_flm, P=Q_(101.325, 'kPa'))
cylinder = Cylinder(D=Q_(0.5, 'inch'), L=Q_(10, 'feet'), fluid=air)
cylinder.u_inf = Q_(5, 'miles / hr')
h_avg = cylinder.avg_heat_trf_coeff()
Q = h_avg * PI * cylinder.D * cylinder.L * (T_water - T_inf)
print(Q.to('W'))
