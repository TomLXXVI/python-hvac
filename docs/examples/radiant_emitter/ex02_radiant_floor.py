"""
EXAMPLE 2: RADIANT FLOOR PANEL
------------------------------

This script demonstrates the usage of class `RadiantFloorPanel` in subpackage
`hvac.radiant_heat_emitter`. This class can be used to size radiant floor panel
circuits and to analyze the performance of a radiant floor panel.
"""
import warnings
import numpy as np
from hvac import Quantity
from hvac.radiant_emitter import RadiantFloorPanel
from hvac.charts import LineChart

warnings.filterwarnings('ignore')  # ignore all warnings

Q_ = Quantity

# Create a `RadiantFloorPanel` passing the necessary design data (see docstring
# for more explanation about the meaning of the design parameters):
rfp = RadiantFloorPanel(
    floor_cover='tiles',
    Q_dot_load_des=Q_(4125, 'W'),
    A_fl=Q_(60, 'm ** 2'),
    T_int_des=Q_(22, 'degC'),
    S=Q_(200, 'mm'),
    T_fl_max=Q_(32, 'degC'),
    L_cir_max=Q_(140, 'm'),
    L_ldr=Q_(1, 'm'),
    DT_w_des=Q_(8, 'K')
)

# Calling method `design` returns a dict with the design results:
d = rfp.design()

print(
    "design conditions\n"
    "-----------------\n"
    f"- water supply temperature: {d['T_w_sup'].to('degC'):~P.1f}\n"
    f"- volume flow rate per circuit: {d['V_dot_w'].to('L / min'):~P.1f}\n"
    f"- average floor temperature: {d['T_fl'].to('degC'):~P.1f}\n"
    f"- number of circuits: {d['N_cir']}\n"
    f"- circuit length: {d['L_cir'].to('m'):~P.2f}\n"
    f"- additional heat supply: {d['Q_dot_add'].to('W'):~P.0f}"
)

# Analyze the heat output of the radiant floor panel for different supply water
# temperatures and water volume rates:

# Set a range of water supply temperatures:
T_w_sup_rng = Q_(np.arange(25, 50, 5), 'degC')

# Set a range of water volume flow rates:
V_dot_w_frac_rng = np.arange(0.10, 1.00, 0.05)
V_dot_w_rng = V_dot_w_frac_rng * Q_(2.5, 'L / min')

# Calculate the heat output for each value of the supply water temperature
# together with each value in the range of water volume flow rates:
Q_dot_rng = {}
for T_w_sup in T_w_sup_rng:
    Q_lst = []
    for V_dot_w in V_dot_w_rng:
        Q = rfp.Q_dot(T_w_sup=T_w_sup, V_dot_w=V_dot_w, T_i=Q_(24, 'degC'))
        Q_lst.append(Q)
    Q_dot_rng[T_w_sup.to('degC').m] = Q_lst

# For each supply water temperature, draw the curve showing the heat output as
# function of the water volume flow rate:
chart = LineChart()
for k in Q_dot_rng.keys():
    chart.add_xy_data(
        label=f'{k} °C',
        x1_values=[V_dot_w_frac * 100 for V_dot_w_frac in V_dot_w_frac_rng],
        y1_values=[Q_dot.to('W').m for Q_dot in Q_dot_rng[k]]
    )
chart.x1.add_title('percent design flow rate')
chart.y1.add_title('heat output, W')
chart.add_legend(columns=5)
chart.show()
