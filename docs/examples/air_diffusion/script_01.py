"""
ANALYSIS OF A COMPACT FREE AIR JET
- local Archimedes number of the jet
- jet trajectory
- jet centerline velocity
- jet centerline temperature
"""
import warnings
import numpy as np
from hvac import Quantity
from hvac.fluids import Fluid, CoolPropWarning
from hvac.charts import LineChart
from hvac.air_diffusion.air_jets.air_jets import CompactJet

warnings.filterwarnings('ignore', category=CoolPropWarning)

Q_ = Quantity
Air = Fluid('Air')
P_atm = Q_(101_325, 'Pa')

# Specify:
# - state of supply air
# - state of room air
# - effective area of supply opening
# - supply air velocity
supply_air = Air(T=Q_(13, 'degC'), P=P_atm)
room_air = Air(T=Q_(26, 'degC'), P=P_atm)
A_o = Q_(1.2215, 'm**2')
U_o = Q_(1.35, 'm / s')

# Create jet:
jet = CompactJet(
    A_o=A_o,
    U_o=U_o,
    supply_air=supply_air,
    room_air=room_air,
    K1=5.8
)

# The horizontal distance according to the equation of the zone 3 centerline
# velocity where the centerline velocity equals the given supply air velocity:
print(jet.x_o)

# Define a range of horizontal distances measured from the supply opening:
x = Q_(np.linspace(0.50, 11.0, 100), 'm')

# LOCAL ARCHIMEDES NUMBER OF AIR JET
Ar_x = abs(jet.local_archimedes_number_zone3(x))
# The jet can be considered as isothermal as long as |Ar_x| < 0.1.

chart1 = LineChart()
chart1.add_xy_data(
    label='local Ar',
    x1_values=x,
    y1_values=Ar_x
)
chart1.add_xy_data(
    label='limit local Ar',
    x1_values=[x[0].m, x[-1].m],
    y1_values=[-0.1, -0.1],
    style_props={'color': 'red'}
)
chart1.x1.add_title('horizontal distance, m')
chart1.y1.add_title('local Archimedes number')
chart1.show()

# TRAJECTORY OF THE JET
y = jet.trajectory(x, alpha=Q_(0, 'deg'))

chart2 = LineChart()
chart2.add_xy_data(
    label='trajectory',
    x1_values=x,
    y1_values=y
)
chart2.x1.add_title('horizontal distance, m')
chart2.y1.add_title('vertical distance, m')
chart2.y1.scale(-6.0, 1.0, 0.5)
chart2.show()

# CENTER-LINE TEMPERATURE OF THE JET
T_m = jet.centerline_temperature_zone3(x)

chart3 = LineChart()
chart3.add_xy_data(
    label='center-line temperature',
    x1_values=x,
    y1_values=T_m.to('degC')
)
chart3.x1.add_title('horizontal distance, m')
chart3.y1.add_title('center-line temperature, °C')
chart3.show()


# CENTER-LINE VELOCITY OF THE JET
U_x = jet.centerline_velocity_zone3(x)

chart4 = LineChart()
chart4.add_xy_data(
    label='center-line velocity',
    x1_values=x,
    y1_values=U_x.to('m / s')
)
chart4.x1.add_title('horizontal distance, m')
chart4.y1.add_title('center-line velocity, m / s')
chart4.show()
