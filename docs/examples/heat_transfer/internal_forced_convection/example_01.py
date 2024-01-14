"""
EXAMPLE 1
---------
From: Nellis G. F., & Klein S. A. (2021)
INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.

Example 9.1
A circular pipe with inner diameter D = 2 inch carries V_dot = 1.5 gal/min of
water at atmospheric pressure and room temperature. Estimate the hydrodynamic
entry length associated with the flow.
"""
from hvac import Quantity
from hvac.fluids import Fluid
from hvac.heat_transfer.forced_convection.internal_flow import CircularTube

Q_ = Quantity

Water = Fluid('Water')
water = Water(T=Q_(20, 'degC'), P=Q_(101.3, 'kPa'))

tube = CircularTube(Di=Q_(2, 'inch'), L=Q_(1, 'm'), fluid=water)
tube.V_dot = Q_(1.5, 'gal / min')
x_fd_h_lam = tube.laminar_hydrodynamic_entry_length()
print(x_fd_h_lam.to('m'))
