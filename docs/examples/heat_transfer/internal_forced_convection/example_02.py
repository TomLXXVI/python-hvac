"""
EXAMPLE 2
---------
From: Nellis G. F., & Klein S. A. (2021)
INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.

Example 9.2
A duct with a rectangular cross-section is carrying a flow of helium gas at
room temperature and atmospheric pressure. The cross-section of the duct is
W = 1 cm x H = 0.8 cm and the length of the duct is L = 10 m. The mass flow rate
of helium is equal to m_dot = 1 g/s and the pressure drop across the duct is
deltaP = 17 kPa.
Determine the Reynolds number and the average friction factor associated with
this situation.
"""
import warnings
from hvac import Quantity
from hvac.fluids import Fluid, CoolPropWarning
from hvac.heat_transfer.forced_convection.internal_flow import RectangularTube


warnings.filterwarnings('ignore', category=CoolPropWarning)

Q_ = Quantity


Helium = Fluid('Helium')
helium = Helium(T=Q_(20, 'degC'), P=Q_(101.325, 'kPa'))
duct = RectangularTube(
    a=Q_(1, 'cm'), b=Q_(0.8, 'cm'), L=Q_(10, 'm'),
    fluid=helium
)
duct.m_dot = Q_(1, 'g / s')
print(duct.reynolds_number())
print(duct.friction_factor())
