"""
EXAMPLE 4
---------
From: Nellis G. F., & Klein S. A. (2021)
INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.

Example 8.4
The heat transfer fluid Therminol 66 flows over a flat plate with u_inf = 3.2
m/s and T_inf = 60 °C. The entire length of the plate is L = 20 cm. In the flow
direction, the initial Luh = 10 cm of the plate is insulated and the remainder
of the plate is maintained at Ts = 15 °C. The width of the plate (perpendicular
to the flow direction) is W = 10 cm. Determine the rate of heat transfer from
the fluid to the plate.
"""
from hvac import Quantity
from hvac.fluids import Fluid
from hvac.heat_transfer.forced_convection.external_flow.plate import Plate

Q_ = Quantity

# plate:
L = Q_(20, 'cm')
W = Q_(10, 'cm')
L_uh = Q_(10, 'cm')
T_surf = Q_(15, 'degC').to('K')

# fluid:
Therminol66 = Fluid('T66', backend='INCOMP')
u_inf = Q_(3.2, 'm / s')
T_inf = Q_(60.0, 'degC').to('K')
T_flm = (T_surf + T_inf) / 2
therminol66 = Therminol66(T=T_flm, P=Q_(101.325, 'kPa'))

print(
    f"rho: {therminol66.rho}\n"
    f"mu: {therminol66.mu}\n"
    f"k: {therminol66.k}"
)

plate = Plate(L_char=L, b=W, fluid=therminol66, L_uh=L_uh)
plate.u_inf = u_inf
h_avg = plate.avg_heat_trf_coeff(therm_bound_cond='uhs')
Q = h_avg * (L - L_uh) * W * (T_inf - T_surf)
print(Q.to('W'))
