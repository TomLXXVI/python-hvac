"""
EXAMPLE 16
----------
From: Duffie, J. A., Beckman, W. A., & Blair, N. (2020).
SOLAR ENGINEERING OF THERMAL PROCESSES, PHOTOVOLTAICS AND WIND.
John Wiley & Sons.

Example 2.20.1, p. 110:
A collector is to be installed in Madison, latitude 43°, at a slope of 60° to
the south. Average daily radiation data are shown in Appendix D. The ground
reflectance is 0.2 for all months except December and March (rho_g = 0.4) and
January and February (rho_g = 0.7). Using the isotropic diffuse assumption,
estimate the monthly average radiation incident on the collector.
"""
from hvac import Quantity
from hvac.charts import LineChart
from hvac.sun.surface import Location, Surface
from hvac.sun.transposition import KTMethod

Q_ = Quantity

location = Location(fi=Q_(43, 'deg'))
surface = Surface(location, gamma=Q_(0, 'deg'), beta=Q_(60, 'deg'))

# Monthly average daily radiation data from Appendix D:
H_avg_arr = Q_([
    6.44, 9.89, 12.86, 16.05, 21.36, 23.04,
    22.58, 20.33, 14.59, 10.48, 6.37, 5.74
], 'MJ / m**2')

H_T_avg_arr = []
for m in range(1, 13):
    match m:
        case 1 | 2:
            rho_g = Q_(0.7, 'frac')
        case 3 | 12:
            rho_g = Q_(0.4, 'frac')
        case _:
            rho_g = Q_(0.2, 'frac')
    KT_model = KTMethod(
        month=m,
        surface=surface,
        H_avg=H_avg_arr[m - 1],
        rho_g=rho_g,
    )
    H_T_avg_arr.append(KT_model.H_T_avg.to('MJ / m**2'))
    print(
        f"Average daily radiation on collector in month {m} = "
        f"{H_T_avg_arr[-1]:~P.3f}"
    )

chart = LineChart()
chart.add_xy_data(
    label='average daily radiation',
    x1_values=[m for m in range(1, 13)],
    y1_values=[H_T_avg.magnitude for H_T_avg in H_T_avg_arr]
)
chart.x1.add_title('month index')
chart.y1.add_title('average daily radiation, MJ/m²')
chart.show()
