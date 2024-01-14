"""
EXAMPLE 5
---------
From: Nellis G. F., & Klein S. A. (2021)
INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.

Example 9.6
A duct with an annular cross-section is formed by sliding an inner tube inside
an outer tube. Engine oil (10 W) is forced to flow through this annular gap
in a heat exchanger. The outer radius of the annular space is r_out = 1.25 inch
and the inner radius is r_in = 1.0 inch. The total length of the passage is
L = 40 ft. The oil pump that provides the flow is able to produce a pressure
rise of deltaP = 5 psi. Determine the volumetric flow rate of engine oil that
will result.
"""
from hvac import Quantity
from hvac.fluids import FluidState
from hvac.heat_transfer.forced_convection.internal_flow import (
    AnnularDuct,
    find_volume_flow_rate
)

Q_ = Quantity

engine_oil = FluidState(None, state_dict={
    'rho': Q_(885, 'kg / m ** 3'),
    'mu': Q_(0.133, 'Pa * s')
})

annular_duct = AnnularDuct(
    ri=Q_(1.0, 'inch'),
    ro=Q_(1.25, 'inch'),
    L=Q_(40, 'ft'),
    fluid=engine_oil
)

dp_pump = Q_(5.0, 'psi')
V_dot = find_volume_flow_rate(
    tube=annular_duct,
    dp=dp_pump,
    V_dot_max=Q_(1.0, 'm ** 3 / s')
)
print(f"{V_dot.to('m ** 3 / s'):~P.3g}")
