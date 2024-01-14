"""
EXAMPLE 7
---------
From: Nellis G. F., & Klein S. A. (2021)
INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.

Example 9.8
A circular pipe carries an ethylene glycol/water solution (40 percent by volume).
The pipe is wrapped in electrical heat tape that provides approximately a
constant heat flux at the surface. The inner diameter if the pipe is D = 1 cm
and the length is L = 5 m. The mass flow rate of the ethylene glycol solution
is m_dot = 0.035 kg/s. Determine the average heat transfer coefficient between
the fluid and the inner surface of the pipe.
"""
import warnings
from hvac import Quantity
from hvac.fluids import Fluid, CoolPropWarning
from hvac.heat_transfer.forced_convection.internal_flow import CircularTube

warnings.filterwarnings('ignore', category=CoolPropWarning)
Q_ = Quantity


EthyleneGlycol = Fluid('AEG', backend='INCOMP', vol_fractions=[Q_(0.4, 'frac')])
glycol = EthyleneGlycol(T=Q_(20, 'degC'), P=Q_(101.325, 'kPa'))

tube = CircularTube(Di=Q_(1, 'cm'), L=Q_(5, 'm'), fluid=glycol)
tube.m_dot = Q_(0.035, 'kg / s')

h_avg = tube.avg_heat_transfer_coefficient(tube.L)
if isinstance(h_avg, tuple):
    print(
        f"h_avg with uniform heat flux = {h_avg[0]}",
        f"h_avg with uniform wall temperature = {h_avg[1]}",
        sep='\n'
    )
else:
    print(f"h_avg = {h_avg}")
