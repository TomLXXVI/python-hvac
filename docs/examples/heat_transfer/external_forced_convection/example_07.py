"""
EXAMPLE 7
---------
From: Nellis G. F., & Klein S. A. (2021)
INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.

Example 8.7
A hot-wire anemometer is used to measure the velocity of air flowing in the duct
of a heating system. The anemometer consists of an electrically heated tungsten
wire with diameter D = 10 µm and length L = 0.75 cm. The air temperature is
T_inf = 20 °C. The anemometer is operated in constant temperature mode, which
means that the current passing through the wire is controlled so that the wire
temperature remains always deltaT = 40 K above the ambient temperature. As the
air velocity increases, the heat transfer coefficient increases and therefore
more current is required to maintain this temperature difference. The voltage
measured across the wire is therefore related to the air velocity.
If the air velocity is u_inf = 5 m/s, determine the voltage across the wire.
Prepare a plot showing the voltage as a function of air velocity.
"""
from math import pi as PI
import numpy as np
from hvac import Quantity
from hvac.fluids import Fluid
from hvac.heat_transfer.forced_convection.external_flow.cylinder import Cylinder


Q_ = Quantity


def tungsten_resistivity(T: Quantity) -> Quantity:
    rho_ref = Q_(5.60e-8, 'ohm * m').m
    T_ref = Q_(20, 'degC').to('K').m
    alpha_ref = Q_(4.50e-3, '1 / K').m
    T = T.to('K').m
    rho = rho_ref * (1 + alpha_ref * (T - T_ref))
    return Q_(rho, 'ohm * m')


D = Q_(10, 'µm')
L = Q_(0.75, 'cm')

T_air = Q_(20, 'degC').to('K')
u_air = Q_(5, 'm / s')

dT = Q_(40, 'K')
T_wire = T_air + dT
T_flm = (T_air + T_wire) / 2

Air = Fluid('Air')
air = Air(T=T_flm, P=Q_(101.325, 'kPa'))

wire = Cylinder(D=D, L=L, fluid=air)
wire.u_inf = u_air
h_avg = wire.avg_heat_trf_coeff()
Q = h_avg * PI * wire.D * wire.L * (T_wire - T_air)

rho_tungsten = tungsten_resistivity(T_wire)
R_tungsten = rho_tungsten * wire.L / (PI * wire.D ** 2 / 4)
I_wire = np.sqrt(Q / R_tungsten)
V_wire = R_tungsten * I_wire
print(I_wire.to('A'), V_wire.to('V'))
