"""
EXAMPLE 3
---------
From: Nellis G. F., & Klein S. A. (2021)
INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.

Example 8.3
A plate has a heater mounted on its underside that provides q = 25 W. The plate
is cooled by a flow of air at u_inf = 2.8 m/s and T_inf = 20 °C. The plate is
L = 10 cm long in the direction that the air is flowing and W = 10 cm wide
perpendicular to the flow. Determine the surface temperature of the plate.
"""
from hvac import Quantity
from hvac.fluids import Fluid
from hvac.heat_transfer.forced_convection.external_flow.plate import Plate

Q_ = Quantity

# plate:
L = Q_(10, 'cm').to('m')
W = Q_(10, 'cm').to('m')
plate = Plate(L_char=L, b=W, fluid=None)

# fluid:
Air = Fluid('Air')
u_air = Q_(2.8, 'm / s')
T_air = Q_(20.0, 'degC').to('K')

Q = Q_(25, 'W')
plate.u_inf = u_air

# initial guess of plate surface temperature:
T_surf = Q_(50.0, 'degC').to('K')
i_max = 20
i = 0
while i <= i_max:
    T_flm = (T_surf + T_air) / 2  # air film temperature
    air = Air(T=T_flm, P=Q_(101.325, 'kPa'))
    plate.fluid = air
    h_avg = plate.avg_heat_trf_coeff()
    T_surf_new = T_air + Q / (h_avg * (L * W))
    if abs(T_surf_new - T_surf) < Q_(0.01, 'K'):
        break
    T_surf = T_surf_new
    i += 1
else:
    print(f'maximum number of iterations was reached')

print(T_surf.to('degC'))