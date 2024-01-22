from pathlib import Path
import numpy as np
from hvac import Quantity
from hvac.fluids import Fluid
from hvac.vapor_compression import VariableSpeedCompressor
from hvac.charts import LineChart


Q_ = Quantity
R134a = Fluid('R134a')


compressor = VariableSpeedCompressor(
    coeff_file=Path("compressor_data/VTZ054-G_R134a.csv"),
    refrigerant_type=R134a,
    units={'m_dot': 'kg / hr', 'speed': '1 / s'}
)

compressor.dT_sh = Q_(10, 'K')
compressor.speed = Q_(4343, '1 / min')

T_evp_rng = Q_(np.linspace(0, 10, endpoint=True), 'degC')
T_cnd_rng = Q_(np.arange(50, 80, 5), 'degC')

output = {}
for T_cnd in T_cnd_rng:
    compressor.Tc = T_cnd
    key = f"{T_cnd.to('degC'):~P.0f}"
    output[key] = {'cmp_m_dot': [], 'cmp_W_dot': []}
    for T_evp in T_evp_rng:
        compressor.Te = T_evp
        output[key]['cmp_m_dot'].append(compressor.m_dot)
        output[key]['cmp_W_dot'].append(compressor.Wc_dot)


chart1 = LineChart(size=(8, 6))
for key in output.keys():
    chart1.add_xy_data(
        label=key,
        x1_values=T_evp_rng.to('degC').m,
        y1_values=[m_dot.to('kg / hr').m for m_dot in output[key]['cmp_m_dot']]
    )
chart1.add_xy_data(
    label='WP1',
    x1_values=[4.972],
    y1_values=[169.313],
    style_props={'marker': 'o', 'linestyle': 'none', 'color': 'red'}
)
chart1.add_xy_data(
    label='WP2',
    x1_values=[6.225],
    y1_values=[144.303],
    style_props={'marker': 'o', 'linestyle': 'none', 'color': 'red'}
)
chart1.x1.add_title('evaporation temperature, °C')
chart1.y1.add_title('refrigerant mass flow rate, kg/h')
chart1.add_legend(columns=3)
chart1.show()

chart2 = LineChart(size=(8, 6))
for key in output.keys():
    chart2.add_xy_data(
        label=key,
        x1_values=T_evp_rng.to('degC').m,
        y1_values=[W_dot.to('kW').m for W_dot in output[key]['cmp_W_dot']]
    )
chart2.add_xy_data(
    label='WP1',
    x1_values=[4.972],
    y1_values=[2.553],
    style_props={'marker': 'o', 'linestyle': 'none', 'color': 'red'}
)
chart2.add_xy_data(
    label='WP2',
    x1_values=[6.225],
    y1_values=[3.163],
    style_props={'marker': 'o', 'linestyle': 'none', 'color': 'red'}
)
chart2.x1.add_title('evaporation temperature, °C')
chart2.y1.add_title('compressor power, kW')
chart2.add_legend(columns=3)
chart2.show()
