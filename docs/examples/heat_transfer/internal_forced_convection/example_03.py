"""
EXAMPLE 3
---------
From: Nellis G. F., & Klein S. A. (2021)
INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.

Example 9.4
The heat transfer fluid Dowtherm A flows through a tube with inner diameter
D = 1.5 cm and length L = 12 m. The mass flow rate of the fluid is m_dot = 0.04
kg/s. Determine the hydrodynamic and thermal entry lengths for this situation.
"""
from hvac import Quantity
from hvac.fluids import FluidState
from hvac.heat_transfer.forced_convection.internal_flow import CircularTube


Q_ = Quantity

dowtherm_A = FluidState(None, state_dict={
    'rho': Q_(1056, 'kg / m**3'),
    'k': Q_(0.1379, 'W '),
    'mu': Q_(4.316e-3, 'Pa * s'),
    'Pr': 49.66
})
dowtherm_A.cp = dowtherm_A.Pr * dowtherm_A.k / dowtherm_A.mu

tube = CircularTube(Di=Q_(1.5, 'cm'), L=Q_(12, 'm'), fluid=dowtherm_A)
tube.m_dot = Q_(0.04, 'kg / s')
Re = tube.reynolds_number()
x_fd_h_lam = tube.laminar_hydrodynamic_entry_length()
x_fd_t_lam = tube.laminar_thermal_entry_length()
print(
    f"Re = {Re}",
    f"hydrodynamic entry length = {x_fd_h_lam}",
    f"thermal entry length = {x_fd_t_lam}",
    sep='\n'
)
