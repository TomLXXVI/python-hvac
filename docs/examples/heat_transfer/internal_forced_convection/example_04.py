"""
EXAMPLE 4
---------
From: Nellis G. F., & Klein S. A. (2021)
INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.

Example 9.5
A drug delivery system must deliver steady flow at rate of V = 24 cm³/min
through a tube with inner diameter D = 500 µm and length L = 2.5 cm. The drug
has properties rho = 1030 kg/m³ and mu = 0.004 Pa.s. Determine the pressure
elevation at the inlet of the tube that is required to provide this flow.
"""
from hvac import Quantity
from hvac.fluids import FluidState
from hvac.heat_transfer.forced_convection.internal_flow import CircularTube


Q_ = Quantity


drug = FluidState(None, state_dict={
        'rho': Q_(1030, 'kg / m ** 3'),
        'mu': Q_(0.004, 'Pa * s')
})

tube = CircularTube(Di=Q_(500, 'µm'), L=Q_(2.5, 'cm'), fluid=drug)
tube.V_dot = Q_(24, 'cm ** 3 / min')
print(
    f"flow regime = {tube.get_flow_condition()}",
    f"mean flow velocity = {tube.mean_velocity()}",
    f"friction factor = {tube.friction_factor()}",
    f"pressure drop = {tube.pressure_drop()}",
    sep='\n'
)
