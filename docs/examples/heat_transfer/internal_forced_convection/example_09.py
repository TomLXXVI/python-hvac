"""
EXAMPLE 9
---------
From: Nellis G. F., & Klein S. A. (2021)
INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.

Example 9.10
A single fuel rod is being cooled by a molten salt. The fuel rod is a solid
cylinder of radius r_in = 0.25 inch and the molten salt flows through an annular
gap formed by a smooth outer tube that is concentric to the rod. The outer
radius of the gap formed between the fuel rod and the outer tube is r_out =
0.375 inch. The mass flow rate of molten salt is m_dot = 0.27 kg/s and the inlet
temperature is T_in = 500 °C. The length of the fuel rod is L = 6 ft. The heat
flux at the surface of the fuel rod varies with position according to the
equation:
```
    q_s(x) = q_s_max * sin(pi * x / L) where q_s_max = 395 kW/m²
```
Determine the pressure drop associated with the molten salt as it flows through
the tube. Plot the mean temperature and the fuel rod surface temperature as a
function of position.
"""
import warnings
import numpy as np
from hvac import Quantity
from hvac.fluids import Fluid, CoolPropWarning
from hvac.charts import LineChart
from hvac.heat_transfer.forced_convection.internal_flow import AnnularDuct
from hvac.heat_transfer.forced_convection.internal_flow.energy_balance import (
    get_fluid_temp_from_heat_flux,
    get_wall_temp_from_heat_flux
)

warnings.filterwarnings('ignore', category=CoolPropWarning)

Q_ = Quantity


# Function that returns the heat flux as a function of position:
def q_fun(x: Quantity) -> Quantity:
    q_max = Q_(395, 'kW / m ** 2')
    x = x.to('m').m
    L = duct.L.to('m').m
    q = q_max * np.sin(np.pi * x / L)
    return q


# Define the molten salt:
MoltenSalt = Fluid('NaK', backend='INCOMP')


# Inlet temperature of the molten salt:
T_in = Q_(500, 'degC').to('K')


# Create the duct and set the mass flow rate:
duct = AnnularDuct(
    ri=Q_(0.25, 'inch'),
    ro=Q_(0.375, 'inch'),
    L=Q_(6, 'ft'),
    fluid=None
)
duct.m_dot = Q_(0.27, 'kg / s')


# Initial guess of the exit temperature of the molten salt:
T_out = Q_(600, 'degC').to('K')


i_max = 20
i = 0
while i <= i_max:
    T_avg = (T_in + T_out) / 2
    duct.fluid = MoltenSalt(T=T_avg, P=Q_(101.325, 'kPa'))
    # Get a new value of the fluid temperature at the exit:
    T_out_new = get_fluid_temp_from_heat_flux(
        x=duct.L,
        q_fun=q_fun,
        cp=duct.fluid.cp,
        m_dot=duct.m_dot,
        P_h=2 * np.pi * duct.ri,
        T_fl_in=T_in
    )
    # Get the corresponding wall temperature at the exit of the tube:
    T_w = get_wall_temp_from_heat_flux(
        q=q_fun(duct.L),
        h=duct.avg_heat_transfer_coefficient(duct.L),
        T_fl=T_out_new,
    )
    # Repeat the loop until the solution converges:
    if abs(T_out_new - T_out) < Q_(0.0001, 'K'):
        print(f'result after {i} iterations')
        break
    T_out = T_out_new
    i += 1
else:
    print('result after maximum iterations')

# Fluid exit temperature:
print(f"fluid exit temperature = {T_out.to('degC')}")

# Pressure drop:
print(
    f"pressure drop = {duct.pressure_drop().to('Pa')}",
    f"flow regime = {duct.get_flow_condition()}",
    f"reynolds number = {duct.reynolds_number()}",
    sep='\n'
)


# Determine the mean fluid temperature and the wall temperature as function
# of position:
x_rng = Q_(np.linspace(0, duct.L.to('ft')), 'ft')
T_fl_arr, T_w_arr = [], []
for x in x_rng:
    T_fl = get_fluid_temp_from_heat_flux(
        x=x,
        q_fun=q_fun,
        cp=duct.fluid.cp,
        m_dot=duct.m_dot,
        P_h=2 * np.pi * duct.ri,
        T_fl_in=T_in
    )
    T_w = get_wall_temp_from_heat_flux(
        q=q_fun(x),
        h=duct.avg_heat_transfer_coefficient(x),
        T_fl=T_fl,
    )
    T_fl_arr.append(T_fl)
    T_w_arr.append(T_w)


# Plot the functions:
chart = LineChart()
chart.add_xy_data(
    label='wall temperature',
    x1_values=x_rng.m.tolist(),
    y1_values=[T_w.to('degC').m for T_w in T_w_arr]
)
chart.add_xy_data(
    label='fluid temperature',
    x1_values=x_rng.m.tolist(),
    y1_values=[T_fl.to('degC').m for T_fl in T_fl_arr]
)
chart.add_legend()
chart.show()
