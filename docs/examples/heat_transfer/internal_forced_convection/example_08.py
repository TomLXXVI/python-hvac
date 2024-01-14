"""
EXAMPLE 8
---------
From: Nellis G. F., & Klein S. A. (2021)
INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.

Example 9.9
A device is being designed to heat a flow of air from an inlet temperature of
T_in = 20 °C to a mean outlet temperature of T_out = 80 °C. The inlet pressure
is p_in = 100 psi(a) and the mass flow rate of air m_dot = 0.01 kg/s. The
preliminary design concept is a drawn copper tube wrapped with a heater that
provides a uniform heat flux to the air. The thickness of the tube wall is 0.035
inch, the length is L = 4 ft long, and the outer diameter is D_o = 0.25 inch.
Determine the pressure drop, the required heat flux on the outer surface of the
tube, and the maximum temperature of the tube material.
"""
import math
from hvac import Quantity
from hvac.fluids import Fluid
from hvac.heat_transfer.forced_convection.internal_flow import CircularTube


Q_ = Quantity

Air = Fluid('Air')

# Inlet conditions of the air:
T_in = Q_(20, 'degC').to('K')
P_in = Q_(100, 'psi').to('Pa')

# Outlet conditions of the air:
T_out = Q_(80, 'degC').to('K')

# Mass flow rate of the air:
m_dot = Q_(0.01, 'kg / s')

# Tube specs:
wt = Q_(0.035, 'inch').to('m')  # wall thickness
Do = Q_(0.25, 'inch').to('m')   # outer tube diameter
L = Q_(4, 'ft').to('m')

# Determine the average specific heat of the air:
air = Air(T=(T_in + T_out) / 2, P=P_in)
cp_air = air.cp

# Heat rate to the air:
Q = m_dot * cp_air * (T_out - T_in)
print(Q.to('W'))

# Outer surface area of the tube and the heat flux to the air:
A_tube = math.pi * Do * L
q = Q / A_tube
print(q.to('W / m ** 2'))

# Create the `CircularTube` object with the inner diameter and set the given
# mass flow rate through the tube:
tube = CircularTube(
    Di=Do - 2 * wt,
    L=L,
    fluid=air,
    e=Q_(0.0015, 'mm')
)
tube.m_dot = m_dot

# Get the pressure drop across the tube:
dp = tube.pressure_drop()
print(dp.to('Pa'))

# Determine the maximum wall temperature:
h = tube.fully_dev_heat_transfer_coefficient()
T_w_max = T_out + q / h
print(T_w_max.to('degC'))
