"""
EXAMPLE 5
---------
From: Nellis G. F., & Klein S. A. (2021)
INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.

Example 8.5
A flat plate experiences a uniform heat flux on its surface of q_s = 5400 W/m².
The plate is L = 1.5 m long in the flow direction. The free stream temperature
is T_inf = 20 °C and the free stream velocity is u_inf = 5 m/s. The properties
of the fluid are rho = 850 kg/m³, mu = 0.005 Pa.s, k = 0.32 W/(m.K), and
Pr = 5.4. Determine the average plate temperature as well as the hottest local
temperature on the plate.
"""
from hvac import Quantity
from hvac.fluids import FluidState
from hvac.heat_transfer.forced_convection.external_flow.plate import Plate, LaminarFlow
from hvac.heat_transfer.forced_convection.external_flow.general import Re_crit, reynolds_number

Q_ = Quantity

# Fluid:
T_inf = Q_(20, 'degC').to('K')
u_inf = Q_(5, 'm / s')
rho = Q_(850, 'kg / m ** 3')
mu = Q_(0.005, 'Pa * s')
k = Q_(0.32, 'W / (m * K)')
Pr = 5.4
cp = k * Pr / mu
fluid = FluidState(None, state_dict={'rho': rho, 'mu': mu, 'k': k, 'cp': cp})

# Plate:
plate = Plate(
    L_char=Q_(1.5, 'm'),
    b=Q_(1, 'm'),
    fluid=fluid
)
plate.u_inf = u_inf

h_avg = plate.avg_heat_trf_coeff(therm_bound_cond='uhf')
print(h_avg)

T_surf_avg = T_inf + Q_(5400, 'W / m ** 2') / h_avg
print(T_surf_avg.to('degC'))

# Distance from leading edge where flow transitions from laminar to turbulent:
x_crit = Re_crit * mu / (rho * u_inf)
print(x_crit.to('m'))

# Laminar flow just before transition:
Re_x = reynolds_number(rho, mu, u_inf, x_crit)
Nu_loc = LaminarFlow.UniformHeatFlux.local_nusselt_number(Re_x, Pr)
print(Nu_loc)
