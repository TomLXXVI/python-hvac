"""
EXAMPLE 1
---------
From: Nellis G. F., & Klein S. A. (2021)
INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.

Example 8.1
A car has a sign mounted on its roof that is L = 4 ft long and H = 1 ft high.
The sign is thin and mounted so that it is parallel to the direction of the
car's motion. Estimate the additional drag force experienced by the car due to
the sign when it is traveling at u_inf = 55 mph.
"""
from hvac import Quantity
from hvac.fluids import Fluid
from hvac.heat_transfer.forced_convection.external_flow.plate import Plate


Q_ = Quantity
Air = Fluid('Air')


air = Air(T=Q_(27, 'degC'), P=Q_(101.325, 'kPa'))

# plate dimensions
L = Q_(4, 'ft')
H = Q_(1, 'ft')

plate = Plate(L_char=L, b=H, fluid=air)
plate.u_inf = Q_(55, 'miles / hr')
tau_shear_avg = plate.avg_shear_stress()
print(tau_shear_avg)
