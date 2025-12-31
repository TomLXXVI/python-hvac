"""
STEADY-STATE ANALYSIS OF A SINGLE-STAGE VAPOR COMPRESSION MACHINE

This script continues from the previous script *ex07_multi_vcm_analysis_part_a.py*.
In the previous script we analyzed the performance of a modeled vapor
compression machine at different compressor speeds, while keeping the other
input variables constant, i.e., the state of air entering the evaporator,
the mass flow rate of air through the evaporator, the state of air entering
the condenser, and the mass flow rate of air through the condenser.
The outputs of the analysis were written to an Excel-file.
In this script we read the Excel-file back into a Pandas DataFrame object, and
we make a plot of the dry-bulb temperature of the air leaving the evaporator
as a function of compressor speed. Using Scipy's `curve_fit` function, we try
to find a function that relates the leaving air temperature at the evaporator
exit to the compressor speed.
The same procedure can also be used to look at other output variables of the
simulation.
"""
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from hvac.charts import LineChart

# Read the Excel-file into a Pandas DataFrame object:
path = Path(
    "C:/Users/Tom/PycharmProjects/python-hvac/docs/examples/"
    "vapor_compression/analysis_results/07_multi_analysis_1.ods"
)

df = pd.read_excel(path)
df = df.dropna()

# Define an estimation function for curve fitting:
def _estimation_fun(x, a, b, c):
    return a * x**2 + b * x + c


# Take the data points and run Scipy's `curve_fit` function:
x_data = df['n_cmp'].values
y_data = df['T_evp_air_out'].values

popt, *_ = curve_fit(_estimation_fun, x_data, y_data)
a, b, c = popt
print(f"a = {a}, b = {b}, c= {c}")

# Create a new set of data points using our estimation functions with the
# parameters we get back from Scipy's `curve_fit` function:
x_fit_data = np.linspace(x_data[0], x_data[-1])
y_fit_data = _estimation_fun(x_fit_data, *popt)

# Plot the calculated data points taken from the Excel-file and draw a line
# through curve-fitted data points:
chart = LineChart()
chart.add_title('Evaporator Leaving Air Temperature')
chart.add_xy_data(
    label='calculated points',
    x1_values=x_data,
    y1_values=y_data,
    style_props={'marker': 'o', 'linestyle': 'none'}
)
chart.add_xy_data(
    label='curve fit',
    x1_values=x_fit_data,
    y1_values=y_fit_data,
)

# Mention the other input variables that were held fixed in the analysis in the
# previous script *ex07_multi_vcm_analysis_part_a.py*.
text = (
    "evaporator air in:\n"
    f"{df.loc[0, 'T_evp_air_in']} °C DB, "
    f"{df.loc[0, 'RH_evp_air_in']} % RH\n"
    "mass flow rate = "
    f"{df.loc[0, 'evp_air_m_dot']} kg/h\n\n"
    "condenser air in:\n"
    f"{df.loc[0, 'T_cnd_air_in']} °C DB, "
    f"{df.loc[0, 'RH_cnd_air_in']} % RH\n"
    "mass flow rate = "
    f"{df.loc[0, 'cnd_air_m_dot']} kg/h"
)
chart.x1.add_title('compressor speed, rpm')
chart.y1.add_title('temperature, °C')
chart.add_note(
    y_pos=0.05,
    text=text,
    font_size=10,
    vert_align='bottom'
)
chart.show()
