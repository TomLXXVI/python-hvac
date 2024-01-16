"""CREATING A PUMP CURVE MODEL
------------------------------
In this script a second order polynomial is derived from the pump
curve of a circulation pump GRUNDFOS MAGNA3 32-120F operating at 4830 rpm.
Points on the pump curve were retrieved from the website of GRUNDFOS and
downloaded in an Excel-file `pump-curve-data.xlsx` residing in the same folder
as this script.

The coefficients of the second order polynomial that models the pump curve will
be used to analyze the pipe network of which the schematic diagram is added to
the same folder as this script. This analysis is done in the script
`analysis_network.py`.
"""
from hvac import Quantity
from hvac.fluids import Fluid
from hvac.fluid_flow import PumpCurve


Q_ = Quantity

Water = Fluid('Water')
water = Water(T=Q_((7 + 12) / 2, 'degC'), P=Q_(1.5, 'bar'))

g = Q_(9.81, 'm / s**2')


# GRUNDFOS MAGNA3 32-120F - rpm 4830

# Measured volume flow rates:
V_dot_rng = Q_([
    0,
    3.6325,
    6.7269483509464,
    9.20794976890782,
    11.2623968040139,
    13.6107776932948,
    16.0456668909742,
    18.6850648005623,
    21.7034806868725
], 'm**3 / hr')

# Corresponding measured pressure heads:
H_rng = Q_([
    11.99,
    12,
    10.29,
    8.568,
    7.212,
    5.702,
    4.168,
    2.535,
    0.688
], 'm')

# Convert pressure head to pressure difference:
dP_rng = water.rho * g * H_rng

# Determine by curve fitting a second order polynomial
# dP = a2 * V_dot**2 + a1 * V_dot + a0 to model the pump curve:
pump_curve_01 = PumpCurve.curve_fitting(
    coordinates=list(zip(V_dot_rng, dP_rng)),
    name='rpm 4830'
)

# Plot the curve-fitted pump curve model and the given, measured points in a
# diagram:
diagram = pump_curve_01.diagram(
    V_ini=V_dot_rng[0],
    V_fin=V_dot_rng[-1],
    V_unit='m**3 / hr',
    dP_unit='kPa'
)
diagram.add_xy_data(
    label='measured',
    x1_values=V_dot_rng.magnitude,
    y1_values=dP_rng.to('kPa').magnitude,
    style_props={'marker': 'o', 'linestyle': 'none'}
)
diagram.show()

# Get the coefficients of the curve-fitted pump curve model:
a0, a1, a2 = pump_curve_01.coefficients
print(
    f"a0 = {a0.to('Pa'):.0f}",
    f"a1 = {a1.to('Pa / (m**3 / s)'):.0f}",
    f"a2 = {a2.to('Pa / (m**3 / s)**2'):.0f}",
    sep='\n'
)
