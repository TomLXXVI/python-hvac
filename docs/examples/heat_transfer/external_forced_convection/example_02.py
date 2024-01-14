"""
EXAMPLE 2
---------
From: Nellis G. F., & Klein S. A. (2021)
INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.

Example 8.2
A flat plate is maintained at a temperature of Ts = 380 K while being subjected
on one side to a flow of water at T_inf = 300 K and u_inf = 3 m/s. The dimension
of the plate in the direction that the water flows is L = 1 m and the dimension
perpendicular to the water flow is W = 0.5 m. Estimate the rate of heat transfer
from the plate.
"""
from hvac import Quantity
from hvac.fluids import Fluid
from hvac.heat_transfer.forced_convection.external_flow.plate import Plate

Q_ = Quantity

# plate dimensions
L = Q_(1, 'm')
W = Q_(0.5, 'm')
T_plate = Q_(380, 'K')

# fluid
T_water = Q_(300, 'K')
T_film = (T_plate + T_water) / 2
Water = Fluid('Water')
water = Water(T=T_film, P=Q_(101.325, 'kPa'))

plate = Plate(L_char=L, b=W, fluid=water)
plate.u_inf = Q_(3, 'm / s')
h_avg = plate.avg_heat_trf_coeff(therm_bound_cond='ust')

Q = h_avg * (W * L) * (T_plate - T_water)
print(Q.to('kW'))
